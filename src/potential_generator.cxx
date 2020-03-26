#include <cstdio> // printf, std::sprintf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::floor
#include <vector> // std::vector

#include "potential_generator.hxx"

#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set, pow2
#include "inline_tools.hxx" // align<n>
#include "constants.hxx" // ::sqrtpi, ::pi
#include "solid_harmonics.hxx" // ::Y00
#include "real_space_grid.hxx" // ::grid_t, ::add_function
#include "radial_grid.hxx" // ::radial_grid_t
#include "chemical_symbol.h" // element_symbols
#include "sho_projection.hxx" // ::sho_add, ::sho_project
#include "exchange_correlation.hxx" // ::lda_PZ81_kernel
#include "boundary_condition.hxx" // ::periodic_images
#include "data_view.hxx" // view2D<T>
#include "fourier_poisson.hxx" // ::fourier_solve
#include "iterative_poisson.hxx" // ::solve
#include "finite_difference.hxx" // ::Laplacian
#include "geometry_analysis.hxx" // ::read_xyz_file
#include "simple_timer.hxx" // // SimpleTimer
#include "control.hxx" // control::get

#ifdef DEVEL
    #include "lossful_compression.hxx" // RDP_lossful_compression
    #include "debug_output.hxx" // dump_to_file
    #include "radial_r2grid.hxx" // radial_r2grid_t
    #include "radial_r2grid.hxx" // r2_axis
#endif
#include "sho_unitary.hxx" // ::Unitary_SHO_Transform<real_t>

// #define OLD_SINGLE_ATOM_UPDATE_INTERFACE
#ifdef  OLD_SINGLE_ATOM_UPDATE_INTERFACE
  #include "single_atom.hxx" // ::update
#else
  #include "single_atom.hxx" // ::atom_update
#endif

// ToDo: restructure: move this into a separate compilation unit
#include "atom_image.hxx"// ::sho_atom_t
#include "grid_operators.hxx" // ::grid_operator_t, ::list_of_atoms
#include "conjugate_gradients.hxx" // ::eigensolve
#include "multi_grid.hxx" // ::restrict3D, ::interpolate3D
#include "density_generator.hxx" // ::density

#include "data_list.hxx" // data_list<T> // ToDo: replace the std::vector<double*> with new constructions


namespace potential_generator {
  // this module makes a DFT calculation based on atoms
  // that live in a spherical potential which is found
  // by projecting the 3D potential.
  // their wave functions do not hybridize but they 
  // feel the effect of the density of neighboring atoms
  
  double constexpr Y00sq = pow2(solid_harmonics::Y00);
  
  double print_stats(double const values[], size_t const all, double const dV=1, char const prefix=' ') {
      double gmin{9e307}, gmax{-gmin}, gsum{0}, gsum2{0};
      for(size_t i = 0; i < all; ++i) {
          gmin = std::min(gmin, values[i]);
          gmax = std::max(gmax, values[i]);
          gsum  += values[i];
          gsum2 += pow2(values[i]);
      } // i
      printf("%c grid stats min %g max %g integral %g avg %g\n", prefix, gmin, gmax, gsum*dV, gsum/all);
      return gsum*dV;
  } // print_stats

  template <typename real_x_t, typename real_y_t>
  void print_compressed(real_x_t const x[], real_y_t const y[], int const n, float const thr=1e-6) { 
      auto const mask = RDP_lossful_compression(x, y, n, thr);
      for(int i = 0; i < n; ++i) {
          if (mask[i]) printf("%g %g\n", x[i], y[i]);
      } // i
      printf("\n");
  } // print_compressed
 
  
//   inline int n_grid_points(double const suggest) { return (int)align<1>((int)std::ceil(suggest)); }
  inline int even(int const any) { return (((any - 1) >> 1) + 1) << 1;}
  inline int n_grid_points(double const suggest) { return (int)even((int)std::ceil(suggest)); }
  
  inline double fold_back(double x, double const cell_extend) { 
      while(x > 0.5*cell_extend) x -= cell_extend; 
      while(x < -.5*cell_extend) x += cell_extend;
      return x;
  } // fold_back
  
  status_t init_geometry_and_grid(real_space_grid::grid_t<1> & g, double **coordinates_and_Z, 
                                  int & natoms, int const echo=0) {
      // SimpleTimer init_function_timer(__FILE__, __LINE__, __func__, echo);
      status_t stat{0};
      
      // compute the self-consistent solution of a single_atom, all states in the core
      // get the spherical core_density and bring it to the 3D grid
      // get the ell=0 compensator charge and add it to the 3D grid
      // envoke exchange_correlation and fourier_poisson
      // add XC and electrostatic potential and zero potential contributions
      // project the total effective potential to each center using bessel_transforms
      // feed back potential shifts into single_atom

      natoms = 0;
      double cell[3];
      auto const geo_file = control::get("geometry.file", "atoms.xyz");
      int bc[3]; // boundary conditions
      stat += geometry_analysis::read_xyz_file(coordinates_and_Z, &natoms, geo_file, cell, bc, echo);

      double const h = control::get("potential_generator.grid.spacing", 0.2378); // works for GeSbTe with alat=6.04
      g = real_space_grid::grid_t<1>(n_grid_points(cell[0]/h), n_grid_points(cell[1]/h), n_grid_points(cell[2]/h));
      if (echo > 1) printf("# use  %d x %d x %d  grid points\n", g[0], g[1], g[2]);
      g.set_boundary_conditions(bc[0], bc[1], bc[2]);
      g.set_grid_spacing(cell[0]/g[0], cell[1]/g[1], cell[2]/g[2]);
      if (echo > 1) printf("# use  %g %g %g  %s grid spacing\n", g.h[0]*Ang, g.h[1]*Ang, g.h[2]*Ang, _Ang);
      for(int d = 0; d < 3; ++d) {
          if (std::abs(g.h[d]*g[d] - cell[d]) >= 1e-6) {
              warn("# grid in %c-direction seems inconsistent, %d * %g differs from %g %s", 
                             'x'+d, g[d], g.h[d]*Ang, cell[d]*Ang, _Ang);
          }
          cell[d] = g.h[d]*g[d];
          assert(std::abs(g.h[d]*g.inv_h[d] - 1) < 4e-16);
      } // d
      if (echo > 1) printf("# cell is  %g %g %g  %s\n", cell[0]*Ang, cell[1]*Ang, cell[2]*Ang, _Ang);
      
      return stat;
  } // init_geometry_and_grid
  
  
  status_t add_smooth_quantities(double values[] // add to this function on a 3D grid
                , real_space_grid::grid_t<1> const & g 
                , int const na, int32_t const nr2[], float const ar2[]
                , view2D<double> const & center
                , int const n_periodic_images, view2D<double> const & periodic_images
                , double const *const *const atom_qnt
                , int const echo=0, int const echo_q=0
                , double const factor=1
                , char const *quantity="smooth quantity") {

          status_t stat(0);
          // add contributions from smooth core densities
          for(int ia = 0; ia < na; ++ia) {
#ifdef DEVEL
              if (echo > 11) {
                  printf("\n## r, %s of atom #%i\n", quantity, ia);
                  print_compressed(r2_axis(nr2[ia], ar2[ia]).data(), atom_qnt[ia], nr2[ia]);
              } // echo
#endif
              double q_added{0};
              for(int ii = 0; ii < n_periodic_images; ++ii) {
                  double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                  double q_added_image = 0;
                  stat += real_space_grid::add_function(values, g, &q_added_image, atom_qnt[ia], nr2[ia], ar2[ia], cnt, factor);
                  if (echo_q > 11) printf("# %g electrons %s of atom #%d added for image #%i\n", q_added_image, quantity, ia, ii);
                  q_added += q_added_image;
              } // periodic images
#ifdef DEVEL
              if (echo_q > 0) {
                  printf("# after adding %g electrons %s of atom #%d:", q_added, quantity, ia);
                  print_stats(values, g.all(), g.dV());
              } // echo
              if (echo_q > 3) printf("# added %s for atom #%d is  %g electrons\n", quantity, ia, q_added);
//            if (echo_q > 3) printf("#    00 compensator charge for atom #%d is %g electrons\n", ia, atom_qlm[ia][00]*Y00inv);
#endif
          } // ia
          return stat;
  } // add_smooth_quantities
  
  status_t init(float const ion=0.f, int const echo=0) {
    
      double constexpr Y00 = solid_harmonics::Y00;
      double constexpr Y00inv = solid_harmonics::Y00inv;

      status_t stat{0};
      
      int na{0};
      double *coordinates_and_Z{nullptr};
      real_space_grid::grid_t<1> g;
      stat += init_geometry_and_grid(g, &coordinates_and_Z, na, echo);

      double const cell[3] = {g[0]*g.h[0], g[1]*g.h[1], g[2]*g.h[2]};
     
      std::vector<float> ionization(na, 0.f);
      if ((ion != 0.0) && (na > 1)) {
          if (echo > 2) printf("# %s distribute ionization of %g electrons between first and last atom\n", __func__, ion);
          ionization[0] = ion; ionization[na - 1] = -ionization[0];
      } // ionized
      
      std::vector<float> Za(na, 0.f); // list of atomic numbers
      view2D<double> center(na, 4, 0.0); // get memory for a list of atomic centers
      { // scope: prepare atomic coordinates
          view2D<double const> const xyzZ(coordinates_and_Z, 4); // wrap as (na,4)
          if (echo > 1) printf("# %s List of Atoms: (coordinates in %s)\n", __func__,_Ang);
          for(int ia = 0; ia < na; ++ia) {
              int const iZ = (int)std::round(xyzZ[ia][3]);
              char const *El = &(element_symbols[2*iZ]); // warning, this is not a null-determined C-string
              if (echo > 4) printf("# %c%c   %15.9f %15.9f %15.9f\n", El[0],El[1],
                              xyzZ[ia][0]*Ang, xyzZ[ia][1]*Ang, xyzZ[ia][2]*Ang);
              Za[ia] = (float)xyzZ[ia][3]; // list of atomic numbers
              for(int d = 0; d < 3; ++d) {
                  center[ia][d] = fold_back(xyzZ[ia][d], cell[d]) + 0.5*(g[d] - 1)*g.h[d]; // w.r.t. to the center of grid point (0,0,0)
              }   center[ia][3] = 0; // 4th component is not used
              if (echo > 1) printf("# relative%12.3f%16.3f%16.3f\n", center[ia][0]*g.inv_h[0],
                                           center[ia][1]*g.inv_h[1], center[ia][2]*g.inv_h[2]);
          } // ia
      } // scope

      float const rcut = 32; // radial grids usually end at 9.45 Bohr
      
      double *periodic_images_ptr{nullptr};
      int const n_periodic_images = boundary_condition::periodic_images(&periodic_images_ptr, cell, g.boundary_conditions(), rcut, echo);
      if (echo > 1) printf("# %s consider %d periodic images\n", __FILE__, n_periodic_images);
      view2D<double> periodic_images(periodic_images_ptr, 4); // wrap
      
      std::vector<double>  sigma_cmp(na, 1.); //
      std::vector<int32_t> numax(na, 3);
      std::vector<int32_t> lmax_qlm(na, -1);
      std::vector<int32_t> lmax_vlm(na, -1);

      // for each atom get sigma, lmax
#ifdef  OLD_SINGLE_ATOM_UPDATE_INTERFACE
      stat += single_atom::update(na, Za.data(), ionization.data(), 0, numax.data(), 
                       sigma_cmp.data(), 0, 0, 0, lmax_vlm.data(), lmax_qlm.data());
#else
      std::vector<double> Za_double(na); set(Za_double.data(), na, Za.data());
      stat += single_atom::atom_update("initialize", na, Za_double.data(), numax.data(), ionization.data());
      stat += single_atom::atom_update("lmax qlm",   na, nullptr,    lmax_qlm.data());
      stat += single_atom::atom_update("lmax vlm",   na, (double*)1, lmax_vlm.data());
      stat += single_atom::atom_update("sigma cmp",  na, sigma_cmp.data());
#endif

      data_list<double> atom_qlm;
      data_list<double> atom_vlm;
      data_list<double> atom_rho;
      data_list<double> atom_mat;
      if (1) { // scope: set up data_list items
          view2D<int32_t> num(4, na, 0); // how many
          for(int ia = 0; ia < na; ++ia) {
              num(0,ia) = pow2(1 + lmax_qlm[ia]); // qlm
              num(1,ia) = pow2(1 + lmax_vlm[ia]); // vlm
              int const ncoeff = sho_tools::nSHO(numax[ia]);
              num(2,ia) =   pow2(ncoeff);         // aDm
              num(3,ia) = 2*pow2(ncoeff);         // aHm and aSm
          } // ia
          // get memory in the form of data_list containers
          atom_qlm = data_list<double>(na, num[0], 0.0);
          atom_vlm = data_list<double>(na, num[1], 0.0);
          atom_rho = data_list<double>(na, num[2], 0.0);
          atom_mat = data_list<double>(na, num[3], 0.0);
      } // scope

      // the r^2-grid is used to bring radial quantities to a Cartesian grid
      std::vector<int32_t> nr2(na, 1 << 12); // 4096
      std::vector<float>   ar2(na, 16.f); // with nr2 == 4096 rcut = 15.998 Bohr
      data_list<double> atom_vbar(na, nr2.data(), 0.0); // zero potentials
      data_list<double> atom_rhoc(na, nr2.data(), 0.0); // core_densities

      sho_unitary::Unitary_SHO_Transform<double> const unitary(9);

#ifdef DEVEL
      std::vector<double> Laplace_Ves(g.all());
#endif
      std::vector<double>  rho(g.all());
      std::vector<double>  Vxc(g.all());
      std::vector<double>  cmp(g.all());
      std::vector<double>  Ves(g.all(), 0.0);
      std::vector<double> Vtot(g.all());

      char const *es_solver_name = control::get("electrostatic.solver", "iterative"); // {"fourier", "iterative", "none"}
      char const es_solver_method = *es_solver_name | 32; // should be one of {'f', 'i', 'n'}

      std::vector<double> rho_valence(g.all(), 0.0);
      
      int const max_scf_iterations = control::get("potential_generator.max.scf", 1.);
      for(int scf_iteration = 0; scf_iteration < max_scf_iterations; ++scf_iteration) {
          // SimpleTimer scf_iteration_timer(__FILE__, __LINE__, "scf_iteration", echo);
          if (echo > 1) printf("\n\n#\n# %s  SCF-Iteration #%d:\n#\n\n", __FILE__, scf_iteration);

#ifdef  OLD_SINGLE_ATOM_UPDATE_INTERFACE
          stat += single_atom::update(na, 0, 0, 0, 0, 0, atom_rhoc.data(), atom_qlm.data());
#else
          stat += single_atom::atom_update("core densities", na, 0, nr2.data(), ar2.data(), atom_rhoc.data());
          stat += single_atom::atom_update("qlm charges", na, 0, 0, 0, atom_qlm.data());
#endif

          set(rho.data(), g.all(), rho_valence.data());
          
          // add contributions from smooth core densities
          stat += add_smooth_quantities(rho.data(), g, na, nr2.data(), ar2.data(), 
                                center, n_periodic_images, periodic_images, atom_rhoc.data(),
                                echo, echo, Y00sq, "smooth core density");
#ifdef DEVEL
          for(int ia = 0; ia < na; ++ia) {
              if (echo > 3) printf("#    00 compensator charge for atom #%d is %g electrons\n", ia, atom_qlm[ia][00]*Y00inv);
          } // ia
#endif

          { // scope: eval the XC potential and energy
              double Exc{0}, Edc{0};
              for(size_t i = 0; i < g.all(); ++i) {
                  auto const exc_i = exchange_correlation::lda_PZ81_kernel(rho[i], Vxc[i]);
                  Exc += rho[i]*exc_i;
                  Edc += rho[i]*Vxc[i]; // double counting correction
              } // i
              Exc *= g.dV(); Edc *= g.dV(); // scale with volume element
              if (echo > 2) printf("# exchange-correlation energy on grid %.12g %s, double counting %.12g %s\n", Exc*eV,_eV, Edc*eV,_eV);
          } // scope

          set(cmp.data(), g.all(), 0.0); // init compensation charge density
          { // scope: solve the Poisson equation
            
              // add compensation charges cmp
              for(int ia = 0; ia < na; ++ia) {
                  if (echo > 6) printf("# 00 compensator charge for atom #%d is %g\n", ia, atom_qlm[ia][00]*Y00inv);
                  double const sigma = sigma_cmp[ia];
                  int    const ellmax = lmax_qlm[ia];
                  std::vector<double> coeff(sho_tools::nSHO(ellmax), 0.0);
                  stat += sho_projection::denormalize_electrostatics(coeff.data(), atom_qlm[ia], ellmax, sigma, unitary, echo);
                  if (echo > 7) printf("# before SHO-adding compensators for atom #%d coeff[000] = %g\n", ia, coeff[0]);
                  for(int ii = 0; ii < n_periodic_images; ++ii) {
                      double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                      stat += sho_projection::sho_add(cmp.data(), g, coeff.data(), ellmax, cnt, sigma, 0);
                  } // periodic images
#ifdef DEVEL
                  if (echo > 1) { // report extremal values of the density on the grid
                      printf("# after adding %g electrons compensator density for atom #%d:", atom_qlm[ia][00]*Y00inv, ia);
                      print_stats(cmp.data(), g.all(), g.dV());
                  } // echo
#endif
              } // ia

              // add compensators cmp to rho
              add_product(rho.data(), g.all(), cmp.data(), 1.);
              if (echo > 1) {
                  printf("\n# augmented charge density grid stats:");
                  print_stats(rho.data(), g.all(), g.dV());
              } // echo

              if ('f' == es_solver_method) {
                  // solve the Poisson equation using a Fast Fourier Transform
                  int ng[3]; double reci[3][4]; 
                  for(int d = 0; d < 3; ++d) { 
                      ng[d] = g[d];
                      set(reci[d], 4, 0.0);
                      reci[d][d] = 2*constants::pi/(ng[d]*g.h[d]);
                  } // d
                  stat += fourier_poisson::fourier_solve(Ves.data(), rho.data(), ng, reci);

              } else if ('n' == es_solver_method) { // "none"
                  warn("electrostatic.solver = %s may lead to unphysical results!", es_solver_name); 

              } else { // default
                  if (echo > 2) printf("# electrostatic.solver = %s\n", es_solver_name);
                  stat += iterative_poisson::solve(Ves.data(), rho.data(), g, echo);

              } // es_solver_method

              // test the potential in real space, find ves_multipoles
              for(int ia = 0; ia < na; ++ia) {
                  double const sigma = sigma_cmp[ia];
                  int    const ellmax = lmax_vlm[ia];
                  int const nc = sho_tools::nSHO(ellmax);
                  std::vector<double> coeff(nc, 0.0);
                  for(int ii = 0; ii < n_periodic_images; ++ii) {
                      std::vector<double> coeff_image(nc, 0.0);
                      double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                      stat += sho_projection::sho_project(coeff_image.data(), ellmax, cnt, sigma, Ves.data(), g, 0);
                      add_product(coeff.data(), nc, coeff_image.data(), 1.0); // need phase factors? no since we only test the electrostatic
                  } // periodic images
                  // SHO-projectors are brought to the grid unnormalized, i.e. p_{000}(0) = 1.0 and p_{200}(0) = -.5

                  stat += sho_projection::renormalize_electrostatics(atom_vlm[ia], coeff.data(), ellmax, sigma, unitary, echo);
#ifdef DEVEL
                  if (echo > 3) {
                      printf("# potential projection for atom #%d v_00 = %.9f %s\n", ia, atom_vlm[ia][00]*Y00*eV,_eV);
                      int const ellmax_show = std::min(ellmax, 2);
                      for(int ell = 1; ell <= ellmax_show; ++ell) {
                          printf("# potential projection for atom #%d v_%im =", ia, ell);
                          double const unitfactor = Y00 * eV * std::pow(Ang, -ell);
                          for(int emm = -ell; emm <= ell; ++emm) {
                              int const lm = sho_tools::lm_index(ell, emm);
                              printf(" %.6f", atom_vlm[ia][lm]*unitfactor);
                          } // emm
                          printf(" %s %s^%i\n", _eV, _Ang, -ell);
                      } // ell
                  } // echo
#endif // DEVEL
              } // ia

              if (echo > 3) {
                  double const Ees = 0.5*dot_product(g.all(), rho.data(), Ves.data())*g.dV();
                  printf("# inner product between rho and Ves = %g %s\n", 2*Ees*eV,_eV);
              } // echo
          } // scope

          // communicate vlm to the atoms, get zero potential, atom-centered Hamiltonian and overlap
#ifdef  OLD_SINGLE_ATOM_UPDATE_INTERFACE
          stat += single_atom::update(na, 0, 0, 0, 0, 0, 0, 0, atom_vlm.data(), 0, 0, atom_vbar.data(), atom_mat.data());
#else
          float mixing_ratio[] = {.49, .51}; // {potential, density}
          stat += single_atom::atom_update("update", na, 0, 0, mixing_ratio, atom_vlm.data());
          stat += single_atom::atom_update("hamiltonian", na, 0, 0, 0, atom_mat.data());
          stat += single_atom::atom_update("zero potential", na, 0, nr2.data(), ar2.data(), atom_vbar.data());
#endif
          
          set(Vtot.data(), g.all(), Vxc.data()); add_product(Vtot.data(), g.all(), Ves.data(), 1.);

          if (echo > 1) {
              printf("\n# Total effective potential (before adding zero potentials)");
              print_stats(Vtot.data(), g.all(), g.dV());
          } // echo

          // now also add the zero potential to Vtot
#if 0          
          for(int ia = 0; ia < na; ++ia) {
#ifdef DEVEL
              if (echo > 11) {
                  printf("\n## r, zero_potential of atom #%i\n", ia);
                  print_compressed(r2_axis(nr2[ia], ar2[ia]).data(), atom_vbar[ia], nr2[ia]);
              } // echo
#endif          
              for(int ii = 0; ii < n_periodic_images; ++ii) {
                  double dummy = 0, cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                  stat += real_space_grid::add_function(Vtot.data(), g, &dummy, atom_vbar[ia], nr2[ia], ar2[ia], cnt, Y00);
              } // ii periodic images
          } // ia
#else
          stat += add_smooth_quantities(Vtot.data(), g, na, nr2.data(), ar2.data(), 
                                center, n_periodic_images, periodic_images, atom_vbar.data(),
                                echo, 0, Y00, "zero potential");
#endif

          if (echo > 1) {
              printf("\n# Total effective potential  (after adding zero potentials)");
              print_stats(Vtot.data(), g.all(), g.dV());
          } // echo

          std::vector<double> rhov_new(g.all(), 0.0); // new valence density

          { // scope: solve the Kohn-Sham equation with the given Hamiltonian
              for(int ia = 0; ia < na; ++ia) {
                  if (echo > 6) {
                      double const *const aHm = atom_mat[ia];
                      int const ncoeff = sho_tools::nSHO(numax[ia]);
                      printf("\n# atom-centered %dx%d Hamiltonian (%s) for atom index %i\n", 
                                            ncoeff, ncoeff, _eV, ia);
                      for(int i = 0; i < ncoeff; ++i) {
                          printf("#%3i  ", i);
                          for(int j = 0; j < ncoeff; ++j) {
                              printf(" %.9g", aHm[i*ncoeff + j]*eV);
                          }   printf("\n");
                      } // i
                  } // echo
              } // ia

              // create a coarse grid descriptor
              real_space_grid::grid_t<1> gc(g[0]/2, g[1]/2, g[2]/2);
              gc.set_grid_spacing(cell[0]/gc[0], cell[1]/gc[1], cell[2]/gc[2]);
              gc.set_boundary_conditions(g.boundary_conditions());

              // restrict the local effective potential to the coarse grid
              std::vector<double> Veff(gc.all());
              multi_grid::restrict3D(Veff.data(), gc, Vtot.data(), g);
              if (echo > 1) {
                  printf("\n# Total effective potential (restricted to coarse grid)");
                  print_stats(Veff.data(), gc.all(), gc.dV());
              } // echo

              std::vector<atom_image::sho_atom_t> a;
              view2D<double> xyzZinso(na, 8);
              for(int ia = 0; ia < na; ++ia) {
                  set(xyzZinso[ia], 4, &coordinates_and_Z[4*ia]); // copy
                  xyzZinso[ia][4] = ia;  // global_atom_id
                  xyzZinso[ia][5] = 3;   // numax
                  xyzZinso[ia][6] = .55; // sigma, ToDo: get sigma_prj from LiveAtom
                  xyzZinso[ia][7] = 0;   // __not_used__
              } // ia
              stat += grid_operators::list_of_atoms(a, xyzZinso.data(), na, 8, gc, echo, atom_mat.data());
              
              // construct grid-based Hamiltonian and overlap operator descriptor
              grid_operators::grid_operator_t<double> op(gc, a, Veff.data());

              std::vector<double> rho_valence_new(gc.all(), 0.0); // new valence density
              std::vector<double*> atom_rho(na, nullptr); // atomic density matrix
              for(int ia = 0; ia < na; ++ia) {
                  int const numax = op.get_numax(ia);
                  int const ncoeff = sho_tools::nSHO(numax);
                  atom_rho[ia] = new double[ncoeff*ncoeff];
                  set(atom_rho[ia], ncoeff*ncoeff, 0.0); // init clear
              } // ia

              int const nbands = 8;
              view2D<double> waves(nbands, gc.all());
              view2D<double> energies(1, nbands);

              auto const KS_solver_method = control::get("eigensolver", "cg");
              for(int ikpoint = 0; ikpoint < 1; ++ikpoint) { // ToDo: implement k-points

                  // solve the Kohn-Sham equation using various solvers
                  if ('c' == KS_solver_method[0]) { // "cg" or "conjugate_gradients"
                      // need start wave-functions, ToDo
                      stat += conjugate_gradients::eigensolve(waves.data(), nbands, op, echo, 1e-6, energies[ikpoint]);
                  } else
                  if ('d' == KS_solver_method[0]) { // "dav" or "davidson"
                      ++stat; error("eigensolver method \'davidson\' needs to be implemented");
                  } else {
                      ++stat; error("unknown eigensolver method \'%s\'", KS_solver_method);
                  }

                  // add to density
                  stat += density_generator::density(rho_valence_new.data(), atom_rho.data(), waves.data(), op, nbands, 1, echo);

              } // ikpoint

              // generate a new density from the eigenstates
              // with occupation numbers from eigenenergies

              stat += multi_grid::interpolate3D(rho_valence.data(), g, rho_valence_new.data(), gc);

              // ToDo: density mixing

          } // scope: Kohn-Sham
          
          
      } // scf_iteration

#ifdef DEVEL
//   return 1; // warning! no cleanup has been run
  
//   printf("\n\n# Early exit in %s line %d\n\n", __FILE__, __LINE__); exit(__LINE__);

      { // scope: compute the Laplacian using high-order finite-differences
          finite_difference::finite_difference_t<double> const fd(g.h, 8);
          {   SimpleTimer timer(__FILE__, __LINE__, "finite-difference", echo);
              stat += finite_difference::Laplacian(Laplace_Ves.data(), Ves.data(), g, fd, -.25/constants::pi);
              { // Laplace_Ves should match rho
                  double res_a{0}, res_2{0};
                  for(size_t i = 0; i < g.all(); ++i) {
                      res_a += std::abs(Laplace_Ves[i] - rho[i]);
                      res_2 +=     pow2(Laplace_Ves[i] - rho[i]);
                  } // i
                  res_a *= g.dV(); res_2 = std::sqrt(res_2*g.dV());
                  if (echo > 1) printf("# Laplace*Ves - rho: residuals abs %.2e rms %.2e\n", res_a, res_2);
              }
          } // timer


          if (0) {
              // Jacobi solver for the Poisson equation: try if Ves is already the solution
              // Problem A*x == f    Jacobi update:  x_{k+1} = D^{-1}(f - T*x_{k}) where D+T==A and D is diagonal

              // SimpleTimer timer(__FILE__, __LINE__, "Jacobi solver", echo);
              finite_difference::finite_difference_t<double> fdj(g.h, 8);
              fdj.scale_coefficients(-.25/constants::pi);
              double const diagonal = fdj.clear_diagonal_elements();
              double const inverse_diagonal = 1/diagonal;
              auto & Ves_copy = Ves;
              std::vector<double> jac_work(g.all());
              for(int k = 0; k < 3; ++k) {
                  double const *x_k = Ves_copy.data(); // read from x_k
                  double      *Tx_k = jac_work.data();
                  double     *x_kp1 = Ves_copy.data();
                  stat += finite_difference::Laplacian(Tx_k, x_k, g, fdj);
                  double res_a{0}, res_2{0}; // internal change between x_{k} and x_{k+1}
                  for(size_t i = 0; i < g.all(); ++i) {
                      double const old_value = x_k[i];
                      double const new_value = inverse_diagonal*(rho[i] - Tx_k[i]);
                      res_a += std::abs(new_value - old_value);
                      res_2 +=     pow2(new_value - old_value);
                      x_kp1[i] = new_value;
                  } // i
                  res_a *= g.dV(); res_2 = std::sqrt(res_2*g.dV());
                  if (echo > 1) {
                      printf("# Jacobi iteration %i: residuals abs %.2e rms %.2e\n", k, res_a, res_2);
                      if (echo > 9) { printf("# it%i Ves grid stats:", k); print_stats(x_kp1, g.all(), g.dV()); }
                  } // echo
              } // k

              {   SimpleTimer timer(__FILE__, __LINE__, "finite-difference onto Jacobi solution", echo);
                  stat += finite_difference::Laplacian(Laplace_Ves.data(), Ves.data(), g, fd, -.25/constants::pi);
                  { // Laplace_Ves should match rho
                      double res_a{0}, res_2{0};
                      for(size_t i = 0; i < g.all(); ++i) {
                          res_a += std::abs(Laplace_Ves[i] - rho[i]);
                          res_2 +=     pow2(Laplace_Ves[i] - rho[i]);
                      } // i
                      res_a *= g.dV(); res_2 = std::sqrt(res_2*g.dV());
                      if (echo > 1) printf("# Laplace*Ves - rho: residuals abs %.2e rms %.2e\n", res_a, res_2);
                  }
              } // timer
          } // Jacobi solver
          
          
      } // scope

//       return 1; // warning! no cleanup has been run


      std::vector<radial_grid_t*> rg(na, nullptr); // pointers to smooth radial grid descriptors
#ifdef  OLD_SINGLE_ATOM_UPDATE_INTERFACE
      stat += single_atom::update(na, 0, 0, rg.data()); // get pointers
#else
      stat += single_atom::atom_update("radial grids", na, 0, 0, 0, reinterpret_cast<double**>(rg.data()));
#endif

      double* const value_pointers[] = {Ves.data(), rho.data(), Laplace_Ves.data(), cmp.data(), Vxc.data(), Vtot.data()};
      //     Ves; // analyze the electrostatic potential
      //     rho; // analyze the augmented density
      //     Laplace_Ves; // analyze the augmented density computed as Laplacian*Ves
      //     cmp; // analyze only the compensator density
      //     Vxc; // analyze the xc potential
      //     Vtot; // analyze the total potential: Vxc + Ves

      for(int iptr = 0; iptr < 1; iptr += 1) { // only loop over the first 1 for electrostatics
          // SimpleTimer timer(__FILE__, __LINE__, "Bessel-projection-analysis", echo);
          auto const values = value_pointers[iptr];

          // report extremal values of what is stored on the grid
          if (echo > 1) { printf("\n# real-space grid stats:"); print_stats(values, g.all(), g.dV()); }

          for(int ia = 0; ia < na; ++ia) {
    //           int const nq = 200; float const dq = 1.f/16; // --> 199/16 = 12.4375 sqrt(Rydberg) =~= pi/(0.25 Bohr)
              float const dq = 1.f/16; int const nq = (int)(constants::pi/(g.smallest_grid_spacing()*dq));
              std::vector<double> qc(nq, 0.0);

//            printf("\n\n# start bessel_projection:\n"); // DEBUG
              {
                  std::vector<double> qc_image(nq, 0.0);
                  for(int ii = 0; ii < n_periodic_images; ++ii) {
                      double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                      stat += real_space_grid::bessel_projection(qc_image.data(), nq, dq, values, g, cnt);
                      add_product(qc.data(), nq, qc_image.data(), 1.0);
                  } // ii
              }
//            printf("\n# end bessel_projection.\n\n"); // DEBUG

              scale(qc.data(), nq, pow2(solid_harmonics::Y00));
              
              std::vector<double> qcq2(nq, 0.0);
              for(int iq = 1; iq < nq; ++iq) {
                  qcq2[iq] = 4*constants::pi*qc[iq]/pow2(iq*dq); // cheap Poisson solver in Bessel transform
              } // iq
    
//               if (echo > 6) {
//                   printf("\n## Bessel coeff for atom #%d:\n", ia);
//                   for(int iq = 0; iq < nq; ++iq) {
//                       printf("%g %g %g\n", iq*dq, qc[iq], qcq2[iq]);
//                   }   printf("\n\n");
//               } // echo

              if (echo > 3) {
                  std::vector<double> rs(rg[ia]->n);
                  bessel_transform::transform_s_function(rs.data(), qc.data(), *rg[ia], nq, dq, true); // transform back to real-space again
                  printf("\n## Real-space projection for atom #%d:\n", ia);
//                   {   auto const mask = RDP_lossful_compression(rg[ia]->r, rs.data(), rg[ia]->n);
//                       for(int ir = 0; ir < rg[ia]->n; ++ir) {
//                           if (mask[ir]) printf("%g %g\n", rg[ia]->r[ir], rs[ir]);
//                       }   printf("\n\n");
//                   } // scope: RDP-mask
                  print_compressed(rg[ia]->r, rs.data(), rg[ia]->n);

                  if ((values == rho.data()) || (values == Laplace_Ves.data())) {
                      bessel_transform::transform_s_function(rs.data(), qcq2.data(), *rg[ia], nq, dq, true); // transform electrostatic solution to real-space
                      printf("\n## Hartree potential computed by Bessel transform for atom #%d:\n", ia);
//                       for(int ir = 0; ir < rg[ia]->n; ++ir) {
//                           printf("%g %g\n", rg[ia]->r[ir], rs[ir]); 
//                       }   printf("\n\n");
                      print_compressed(rg[ia]->r, rs.data(), rg[ia]->n);
                  } // density
              } // echo
          } // ia

      } // iptr loop for different quantities represented on the grid.


      // generate a file which contains the full potential Vtot
      if (echo > 0) {
          char title[96]; std::sprintf(title, "%i x %i x %i", g[2], g[1], g[0]);
          dump_to_file("vtot.dat", Vtot.size(), Vtot.data(), nullptr, 1, 1, title, echo);
      } // unless all output is suppressed

#endif // DEVEL

      delete[] periodic_images_ptr;
      delete[] coordinates_and_Z;

#ifdef  OLD_SINGLE_ATOM_UPDATE_INTERFACE
      stat += single_atom::update(-na); // memory cleanup
#else
      stat += single_atom::atom_update("memory cleanup", na);
#endif
      return stat;
  } // init


#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_init(int const echo=3) {
      float const ion = control::get("potential_generator.test.ion", 0.0);
      return init(ion, echo); // ionization of Al-P dimer by -ion electrons
  } // test_init

  status_t all_tests(int const echo) {
    status_t status(0);
    status += test_init(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace potential_generator
