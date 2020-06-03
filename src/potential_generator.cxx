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
#include "real_space.hxx" // ::grid_t, ::add_function
#include "radial_grid.hxx" // ::radial_grid_t
#include "chemical_symbol.hxx" // ::get
#include "sho_projection.hxx" // ::sho_add, ::sho_project
#include "exchange_correlation.hxx" // ::lda_PZ81_kernel
#include "boundary_condition.hxx" // ::periodic_images
#include "data_view.hxx" // view2D<T>
#include "fourier_poisson.hxx" // ::fourier_solve
#include "iterative_poisson.hxx" // ::solve
#include "finite_difference.hxx" // ::stencil_t, ::derive
#include "geometry_analysis.hxx" // ::read_xyz_file
#include "simple_timer.hxx" // // SimpleTimer
#include "control.hxx" // control::get

#ifdef DEVEL
    #include "lossful_compression.hxx" // print_compressed
    #include "debug_output.hxx" // dump_to_file
    #include "debug_tools.hxx" // ::read_from_file
    #include "radial_r2grid.hxx" // radial_r2grid_t
    #include "radial_r2grid.hxx" // r2_axis
#endif
#include "sho_unitary.hxx" // ::Unitary_SHO_Transform<real_t>

#include "single_atom.hxx" // ::atom_update

// ToDo: restructure: move this into a separate compilation unit
#include "atom_image.hxx"// ::sho_atom_t
#include "grid_operators.hxx" // ::grid_operator_t, ::list_of_atoms
#include "conjugate_gradients.hxx" // ::eigensolve
#include "davidson_solver.hxx" // ::rotate, ::eigensolve
#include "multi_grid.hxx" // ::restrict3D, ::interpolate3D
#include "density_generator.hxx" // ::density

#include "data_list.hxx" // data_list<T> // ToDo: replace the std::vector<double*> with new constructions


namespace potential_generator {
  // this module makes a DFT calculation based on atoms
  // that live in a spherical potential which is found
  // by projecting the 3D potential.
  // their wave functions do not hybridize but they 
  // feel the effect of the density of neighboring atoms
  
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
  
//   inline int n_grid_points(double const suggest) { return (int)align<1>((int)std::ceil(suggest)); }
  inline int even(int const any) { return (((any - 1) >> 1) + 1) << 1;}
  inline int n_grid_points(double const suggest) { return (int)even((int)std::ceil(suggest)); }
  
  inline double fold_back(double x, double const cell_extend) { 
      while(x > 0.5*cell_extend) x -= cell_extend; 
      while(x < -.5*cell_extend) x += cell_extend;
      return x;
  } // fold_back
  
  status_t init_geometry_and_grid(real_space::grid_t & g, double **coordinates_and_Z, 
                                  int & natoms, int const echo=0) {
      // SimpleTimer init_function_timer(__FILE__, __LINE__, __func__, echo);
      status_t stat{0};
      
      natoms = 0;
      double cell[3];
      auto const geo_file = control::get("geometry.file", "atoms.xyz");
      int bc[3]; // boundary conditions
      stat += geometry_analysis::read_xyz_file(coordinates_and_Z, &natoms, geo_file, cell, bc, echo);

      double const h = control::get("potential_generator.grid.spacing", 0.2378); // works for GeSbTe with alat=6.04
      g = real_space::grid_t(n_grid_points(cell[0]/h), n_grid_points(cell[1]/h), n_grid_points(cell[2]/h));
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
                , real_space::grid_t const & g 
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
                  print_compressed(radial_r2grid::r_axis(nr2[ia], ar2[ia]).data(), atom_qnt[ia], nr2[ia]);
              } // echo
#endif
              double q_added{0};
              for(int ii = 0; ii < n_periodic_images; ++ii) {
                  double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                  double q_added_image = 0;
                  stat += real_space::add_function(values, g, &q_added_image, atom_qnt[ia], nr2[ia], ar2[ia], cnt, factor);
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
    
      // compute the self-consistent solution of a single_atom, all states in the core
      // get the spherical core_density and bring it to the 3D grid
      // get the ell=0 compensator charge and add it to the 3D grid
      // envoke exchange_correlation and fourier_poisson
      // add XC and electrostatic potential and zero potential contributions
      // project the total effective potential to each center using bessel_transforms
      // feed back potential shifts into single_atom

      double constexpr Y00 = solid_harmonics::Y00;
      double constexpr Y00inv = solid_harmonics::Y00inv;
      double constexpr Y00sq = pow2(solid_harmonics::Y00);

      status_t stat{0};
      
      int na{0};
      double *coordinates_and_Z{nullptr};
      real_space::grid_t g;
      stat += init_geometry_and_grid(g, &coordinates_and_Z, na, echo);

      double const cell[3] = {g[0]*g.h[0], g[1]*g.h[1], g[2]*g.h[2]};
     
      std::vector<float> ionization(na, 0.f);
      if ((ion != 0.0) && (na > 1)) {
          if (echo > 2) printf("# %s distribute ionization of %g electrons between first and last atom\n", __func__, ion);
          ionization[0] = ion; ionization[na - 1] = -ionization[0];
      } // ionized

      
      float const rcut = 32; // radial grids usually end at 9.45 Bohr
      double *periodic_images_ptr{nullptr};
      int const n_periodic_images = boundary_condition::periodic_images(&periodic_images_ptr, cell, g.boundary_conditions(), rcut, echo);
      if (echo > 1) printf("# %s consider %d periodic images\n", __FILE__, n_periodic_images);
      view2D<double> periodic_images(periodic_images_ptr, 4); // wrap

      
      std::vector<double> Za(na);       // list of atomic numbers
      view2D<double> center(na, 4, 0.0); // get memory for a list of atomic centers
      { // scope: prepare atomic coordinates
          view2D<double const> const xyzZ(coordinates_and_Z, 4); // wrap as (na,4)
          if (echo > 1) printf("# %s List of Atoms: (coordinates in %s)\n", __func__,_Ang);
          for(int ia = 0; ia < na; ++ia) {
              double const Z = xyzZ[ia][3];
              char Symbol[4]; chemical_symbol::get(Symbol, Z, ' ');
              if (echo > 4) printf("# %s  %15.9f %15.9f %15.9f\n", Symbol,
                              xyzZ[ia][0]*Ang, xyzZ[ia][1]*Ang, xyzZ[ia][2]*Ang);
              Za[ia] = Z;
              for(int d = 0; d < 3; ++d) {
                  center[ia][d] = fold_back(xyzZ[ia][d], cell[d]) + 0.5*(g[d] - 1)*g.h[d]; // w.r.t. to the center of grid point (0,0,0)
              }   center[ia][3] = 0; // 4th component is not used
              if (echo > 1) printf("# relative%12.3f%16.3f%16.3f\n", center[ia][0]*g.inv_h[0],
                                           center[ia][1]*g.inv_h[1], center[ia][2]*g.inv_h[2]);
          } // ia
      } // scope

      std::vector<double>  sigma_cmp(na, 1.);
      std::vector<int32_t> numax(na, 3);
      std::vector<int32_t> lmax_qlm(na, -1);
      std::vector<int32_t> lmax_vlm(na, -1);

      // for each atom get sigma, lmax
      stat += single_atom::atom_update("initialize", na, Za.data(), numax.data(), ionization.data(), (double**)1);
      stat += single_atom::atom_update("lmax qlm",   na, nullptr,    lmax_qlm.data());
      stat += single_atom::atom_update("lmax vlm",   na, (double*)1, lmax_vlm.data());
      stat += single_atom::atom_update("sigma cmp",  na, sigma_cmp.data());

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
      data_list<double> atom_vbar(nr2, 0.0); // zero potentials
      data_list<double> atom_rhoc(nr2, 0.0); // core_densities

      sho_unitary::Unitary_SHO_Transform<double> const unitary(9);

#ifdef DEVEL
      std::vector<double> Laplace_Ves(g.all(), 0.0);
#endif
      std::vector<double>  rho(g.all());
      std::vector<double>  Vxc(g.all());
      std::vector<double>  cmp(g.all());
      std::vector<double>  Ves(g.all(), 0.0);
      std::vector<double> Vtot(g.all());

      char const *es_solver_name = control::get("electrostatic.solver", "multi-grid"); // {"fourier", "multi-grid", "CG", "SD", "none"}
      char const es_solver_method = *es_solver_name; // should be one of {'f', 'i', 'n'}

      std::vector<double> q00_valence;
      std::vector<double> rho_valence(g.all(), 0.0);
      char const *init_val_rho_name = control::get("initial.valence.density", "atomic"); // {"atomic", "load", "none"}
      if (echo > 0) printf("\n# initial.valence.density = \'%s\'\n", init_val_rho_name);
      if ('a' == init_val_rho_name[0]) {
          auto & atom_rhov = atom_rhoc; // temporarily rename the existing allocation
          q00_valence.resize(na);
          stat += single_atom::atom_update("valence densities", na, q00_valence.data(), nr2.data(), ar2.data(), atom_rhov.data());
          // add smooth spherical atom-centered valence densities
          stat += add_smooth_quantities(rho_valence.data(), g, na, nr2.data(), ar2.data(), 
                                center, n_periodic_images, periodic_images, atom_rhov.data(),
                                echo, echo, Y00sq, "smooth atomic valence density");

      } else
      if ('n' == init_val_rho_name[0]) {
          warn("initial.valence.density = none may cause problems");
      } // switch

      int const max_scf_iterations = control::get("potential_generator.max.scf", 1.);
      for(int scf_iteration = 0; scf_iteration < max_scf_iterations; ++scf_iteration) {
          // SimpleTimer scf_iteration_timer(__FILE__, __LINE__, "scf_iteration", echo);
          if (echo > 1) printf("\n\n#\n# %s  SCF-Iteration #%d:\n#\n\n", __FILE__, scf_iteration);

          stat += single_atom::atom_update("core densities", na, 0, nr2.data(), ar2.data(), atom_rhoc.data());
          stat += single_atom::atom_update("qlm charges", na, 0, 0, 0, atom_qlm.data());

          if (echo > 4) { // report extremal values of the valence density on the grid
              printf("# valence density in iteration #%i:", scf_iteration);
              print_stats(rho_valence.data(), g.all(), g.dV());
          } // echo
          
          set(rho.data(), g.all(), rho_valence.data());
          
          // add contributions from smooth core densities
          stat += add_smooth_quantities(rho.data(), g, na, nr2.data(), ar2.data(), 
                                center, n_periodic_images, periodic_images, atom_rhoc.data(),
                                echo, echo, Y00sq, "smooth core density");

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
                  if (q00_valence.size() == na) { 
                      atom_qlm[ia][00] += Y00*q00_valence[ia]; // add the valence charge deficit from the initial rhov
                      if (echo > -1) printf("# add %.6f electrons valence charge deficit of atom #%i\n", q00_valence[ia], ia);
                  } // q00_valence
                  double const sigma = sigma_cmp[ia];
                  int    const ellmax = lmax_qlm[ia];
                  std::vector<double> coeff(sho_tools::nSHO(ellmax), 0.0);
                  stat += sho_projection::denormalize_electrostatics(coeff.data(), atom_qlm[ia], ellmax, sigma, unitary, echo);
//                if (echo > 7) printf("# before SHO-adding compensators for atom #%d coeff[000] = %g\n", ia, coeff[0]);
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
              q00_valence.clear(); // only in the 1st iteration

              // add compensators cmp to rho
              add_product(rho.data(), g.all(), cmp.data(), 1.);
              if (echo > 1) {
                  printf("\n# augmented charge density:");
                  print_stats(rho.data(), g.all(), g.dV());
              } // echo

              { SimpleTimer timer(__FILE__, __LINE__, "Poisson equation", echo);
              // solve the Poisson equation: Laplace Ves == -4 pi rho
              if ('f' == (es_solver_method | 32)) { // "fft", "fourier" 
                  // solve the Poisson equation using a Fast Fourier Transform
                  int ng[3]; double reci[3][4]; 
                  for(int d = 0; d < 3; ++d) { 
                      ng[d] = g[d];
                      set(reci[d], 4, 0.0);
                      reci[d][d] = 2*constants::pi/(ng[d]*g.h[d]);
                  } // d
                  stat += fourier_poisson::fourier_solve(Ves.data(), rho.data(), ng, reci);

              } else if ('M' == es_solver_method) { // "Multi-grid" (upper case!)
                
                  // create a denser grid descriptor
                  real_space::grid_t gd(g[0]*2, g[1]*2, g[2]*2);
                  if (echo > 2) printf("# electrostatic.solver = %s is a multi-grid solver"
                          " on a %d x %d x %d grid\n", es_solver_name, gd[0], gd[1], gd[2]);
                  gd.set_grid_spacing(g.h[0]/2, g.h[1]/2, g.h[2]/2);
                  gd.set_boundary_conditions(g.boundary_conditions());

                  std::vector<double> Ves_dense(gd.all()), rho_dense(gd.all());
                  multi_grid::interpolate3D(rho_dense.data(), gd, rho.data(), g, echo);
                  
                  iterative_poisson::solve(Ves_dense.data(), rho_dense.data(), gd, 'M', echo);
                  
                  // restrict the electrostatic potential to grid g
                  multi_grid::restrict3D(Ves.data(), g, Ves_dense.data(), gd, echo);

              } else if ('n' == (es_solver_method | 32)) { // "none"
                  warn("electrostatic.solver = %s may lead to unphysical results!", es_solver_name); 
                  
              } else if ('l' == (es_solver_method | 32)) { // "load"
#ifdef DEVEL                  
                  auto const Ves_in_filename = control::get("electrostatic.potential.from.file", "v_es.dat");
                  auto const nerrors = debug_tools::read_from_file(Ves.data(), Ves_in_filename, g.all(), 1, 1, "electrostatic potential", echo);
                  if (nerrors) warn("electrostatic.solver = %s from file %s had %d errors", es_solver_name, Ves_in_filename, nerrors); 
#else
                  warn("electrostatic.solver = %s failed!", es_solver_name); 
#endif                  
              } else { // default
                  if (echo > 2) printf("# electrostatic.solver = %s\n", es_solver_name);
                  stat += iterative_poisson::solve(Ves.data(), rho.data(), g, es_solver_method, echo);

              } // es_solver_method
              } // timer

#ifdef DEVEL
              { // scope: export electrostatic potential to ASCII file
                  auto const Ves_out_filename = control::get("electrostatic.potential.to.file", "");
                  if (*Ves_out_filename) stat += dump_to_file(Ves_out_filename, g.all(), Ves.data(), nullptr, 1, 1, "electrostatic potential", echo);
              } // scope
#endif

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

#ifdef DEVEL
//        exit(__LINE__);
#endif
          
          // communicate vlm to the atoms, get zero potential, atom-centered Hamiltonian and overlap
          float mixing_ratio[] = {.5, .5}; // {potential, density}
          stat += single_atom::atom_update("update", na, 0, 0, mixing_ratio, atom_vlm.data());
          stat += single_atom::atom_update("hamiltonian", na, 0, 0, 0, atom_mat.data());
          stat += single_atom::atom_update("zero potentials", na, 0, nr2.data(), ar2.data(), atom_vbar.data());

          set(Vtot.data(), g.all(), Vxc.data()); add_product(Vtot.data(), g.all(), Ves.data(), 1.);

          if (echo > 1) {
              printf("\n# Total effective potential (before adding zero potentials)");
              print_stats(Vtot.data(), g.all(), g.dV());
          } // echo

          // now also add the zero potential vbar to Vtot
          stat += add_smooth_quantities(Vtot.data(), g, na, nr2.data(), ar2.data(), 
                                center, n_periodic_images, periodic_images, atom_vbar.data(),
                                echo, 0, Y00, "zero potential");

          if (echo > 1) {
              printf("\n# Total effective potential  (after adding zero potentials)");
              print_stats(Vtot.data(), g.all(), g.dV());
          } // echo

         
          /**  
           *  Potential generation done
           */

          
          std::vector<double> rhov_new(g.all(), 0.0); // new valence density

          { // scope: solve the Kohn-Sham equation with the given Hamiltonian
#ifdef DEVEL
              for(int ia = 0; ia < na; ++ia) {
                  if (echo > 6) {
                      int const n = sho_tools::nSHO(numax[ia]);
                      view2D<double const> const aHm(atom_mat[ia], n);
                      printf("\n# atom-centered %dx%d Hamiltonian (%s) for atom index #%i\n", n, n, _eV, ia);
                      for(int i = 0; i < n; ++i) {
                          printf("#%3i  ", i);
                          for(int j = 0; j < n; ++j) {
                              printf(" %.3f", aHm(i,j)*eV);
                          }   printf("\n");
                      } // i
                  } // echo
              } // ia
#endif
              // create a coarse grid descriptor
              real_space::grid_t gc(g[0]/2, g[1]/2, g[2]/2); // divide the dense grid numbers by two
              gc.set_grid_spacing(cell[0]/gc[0], cell[1]/gc[1], cell[2]/gc[2]); // alternative: 2*g.h[]
              gc.set_boundary_conditions(g.boundary_conditions());

              // restrict the local effective potential to the coarse grid
              std::vector<double> Veff(gc.all());
              multi_grid::restrict3D(Veff.data(), gc, Vtot.data(), g, 0); // mute
              if (echo > 1) {
                  printf("\n# Total effective potential  (restricted to coarse grid)   ");
                  print_stats(Veff.data(), gc.all(), gc.dV());
              } // echo
              
#ifdef DEVEL
              if (0) { // scope: interpolate the effective potential to the dense grid again and compare it to the original version Vtot
                // in order to test the interpolation routine
                  std::vector<double> v_dcd(g.all(), 0.0);
                  multi_grid::interpolate3D(v_dcd.data(), g, Veff.data(), gc, 0); // mute
                  if (echo > 1) {
                      printf("\n# Total effective potential (interpolated to dense grid)   ");
                      print_stats(v_dcd.data(), g.all(), g.dV());
                  } // echo
              } // scope
#endif

              // create a list of atoms
              view2D<double> xyzZinso(na, 8);
              for(int ia = 0; ia < na; ++ia) {
                  set(xyzZinso[ia], 4, &coordinates_and_Z[4*ia]); // copy
                  xyzZinso[ia][4] = ia;  // global_atom_id
                  xyzZinso[ia][5] = 3;   // numax
                  xyzZinso[ia][6] = .55; // sigma, ToDo: get sigma_prj from LiveAtom
                  xyzZinso[ia][7] = 0;   // __not_used__
              } // ia
              std::vector<atom_image::sho_atom_t> a;
              stat += grid_operators::list_of_atoms(a, xyzZinso.data(), na, 8, gc, echo, atom_mat.data());
              
              // construct grid-based Hamiltonian and overlap operator descriptor
              grid_operators::grid_operator_t<double> const op(gc, a, Veff.data());

              std::vector<double> rho_valence_new(gc.all(), 0.0); // new valence density
              std::vector<int> ncoeff_squared(na), ncoeff_a(na);
              for(int ia = 0; ia < na; ++ia) {
                  int const numax = op.get_numax(ia);
                  int const ncoeff = sho_tools::nSHO(numax);
                  ncoeff_a[ia] = ncoeff;
                  ncoeff_squared[ia] = pow2(ncoeff);
              } // ia
              data_list<double> atom_rho(ncoeff_squared, 0.0); // atomic density matrices

              int const nbands = 8;
              view2D<double> waves(nbands, gc.all()); // Kohn-Sham states
              { // scope: generate start waves from atomic orbitals
                  uint8_t qn[20][4]; // quantum numbers [nx, ny, nz, nu] with nu==nx+ny+nz
                  sho_tools::construct_index_table<sho_tools::order_Ezyx>(qn, 3);
                  data_list<double> single_atomic_orbital(ncoeff_a, 0.0);
                  for(int iband = 0; iband < nbands; ++iband) {
                      int const ia = iband % na, io = iband / na; // which atom? which orbital?
                      if (echo > 7) printf("# %s initialize band #%i as atomic orbital %x%x%x of atom #%i\n", 
                                              __func__, iband, qn[io][0], qn[io][1], qn[io][2], ia);
                      int const isho = sho_tools::zyx_index(op.get_numax(ia), qn[io][0], qn[io][1], qn[io][2]);
                      single_atomic_orbital[ia][isho] = 1;
                      op.get_start_waves(waves[iband], single_atomic_orbital.data(), echo);
                      single_atomic_orbital[ia][isho] = 0;
//                    print_stats(waves[iband], gc.all(), gc.dV());
                  } // iband
              } // scope
              view2D<double> energies(1, nbands); // Kohn-Sham eigenenergies

              auto const eigensolver_method = control::get("eigensolver", "cg");
              for(int ikpoint = 0; ikpoint < 1; ++ikpoint) { // ToDo: implement k-points

                  // solve the Kohn-Sham equation using various solvers
                  if ('c' == eigensolver_method[0]) { // "cg" or "conjugate_gradients"
                      stat += davidson_solver::rotate(waves.data(), nbands, op, echo);
                      stat += conjugate_gradients::eigensolve(waves.data(), nbands, op, echo, 1e-6, energies[ikpoint]);
                  } else
                  if ('d' == eigensolver_method[0]) { // "davidson"
                      stat += davidson_solver::eigensolve(waves.data(), nbands, op, echo + 9);
                  } else {
                      ++stat; error("unknown eigensolver method \'%s\'", eigensolver_method);
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
//    return 1; // warning! no cleanup has been run
//    printf("\n\n# Early exit in %s line %d\n\n", __FILE__, __LINE__); exit(__LINE__);

      int const verify_Poisson = int(control::get("potential_generator.verify.poisson", 0.));
      if (verify_Poisson)
      { // scope: compute the Laplacian using high-order finite-differences

          for(int nfd = 1; nfd < verify_Poisson; ++nfd) {
              finite_difference::stencil_t<double> const fd(g.h, nfd, -.25/constants::pi);
              {   SimpleTimer timer(__FILE__, __LINE__, "finite-difference", echo);
                  stat += finite_difference::apply(Laplace_Ves.data(), Ves.data(), g, fd);
                  { // Laplace_Ves should match rho
                      double res_a{0}, res_2{0};
                      for(size_t i = 0; i < g.all(); ++i) {
                          res_a += std::abs(Laplace_Ves[i] - rho[i]);
                          res_2 +=     pow2(Laplace_Ves[i] - rho[i]);
                      } // i
                      res_a *= g.dV(); res_2 = std::sqrt(res_2*g.dV());
                      if (echo > 1) printf("# Laplace*Ves - rho: residuals abs %.2e rms %.2e (FD-order=%i)\n", res_a, res_2, nfd);
                  }
              } // timer
          } // nfd

      } // scope

      int const use_Bessel_projection = int(control::get("potential_generator.use.bessel.projection", 0.0));
      if (use_Bessel_projection) 
      { // scope: use a Bessel projection around each atom position to compare 3D and radial quantities
        

          std::vector<radial_grid_t*> rg(na, nullptr); // pointers to smooth radial grid descriptors
          stat += single_atom::atom_update("radial grids", na, 0, 0, 0, reinterpret_cast<double**>(rg.data()));

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
              if (echo > 1) { printf("\n# real-space:"); print_stats(values, g.all(), g.dV()); }

              for(int ia = 0; ia < na; ++ia) {
                  float const dq = 1.f/16; int const nq = int(constants::pi/(g.smallest_grid_spacing()*dq));
                  std::vector<double> qc(nq, 0.0);

                  { // scope: Bessel core
                      std::vector<double> qc_image(nq, 0.0);
                      for(int ii = 0; ii < n_periodic_images; ++ii) {
                          double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                          stat += real_space::bessel_projection(qc_image.data(), nq, dq, values, g, cnt);
                          add_product(qc.data(), nq, qc_image.data(), 1.0);
                      } // ii
                  } // scope

                  scale(qc.data(), nq, Y00sq);
                  
                  std::vector<double> qcq2(nq, 0.0);
                  for(int iq = 1; iq < nq; ++iq) {
                      qcq2[iq] = 4*constants::pi*qc[iq]/pow2(iq*dq); // cheap Poisson solver in Bessel transform
                  } // iq
        
                  if (echo > 9) {
                      printf("\n## Bessel coeff for atom #%d:\n", ia);
                      for(int iq = 0; iq < nq; ++iq) {
                          printf("%g %g %g\n", iq*dq, qc[iq], qcq2[iq]);
                      }   printf("\n\n");
                  } // echo

                  if (echo > 3) {
                      std::vector<double> rs(rg[ia]->n);
                      bessel_transform::transform_s_function(rs.data(), qc.data(), *rg[ia], nq, dq, true); // transform back to real-space again
                      printf("\n## Real-space projection for atom #%d:\n", ia);
                      print_compressed(rg[ia]->r, rs.data(), rg[ia]->n);

                      if ((values == rho.data()) || (values == Laplace_Ves.data())) {
                          bessel_transform::transform_s_function(rs.data(), qcq2.data(), *rg[ia], nq, dq, true); // transform electrostatic solution to real-space
                          printf("\n## Hartree potential computed by Bessel transform for atom #%d:\n", ia);
                          print_compressed(rg[ia]->r, rs.data(), rg[ia]->n);
                      } // density
                  } // echo
              } // ia

          } // iptr loop for different quantities represented on the grid.

      } // scope Bessel

      int const export_Vtot = int(control::get("potential_generator.export.vtot", 1.));
      if (export_Vtot) {
          // generate a file which contains the full potential Vtot
          if (echo > 0) {
              char title[96]; std::sprintf(title, "%i x %i x %i", g[2], g[1], g[0]);
              dump_to_file("vtot.dat", Vtot.size(), Vtot.data(), nullptr, 1, 1, title, echo);
          } // unless all output is suppressed
      } // export_Vtot

#endif // DEVEL

      delete[] periodic_images_ptr;
      delete[] coordinates_and_Z;

      stat += single_atom::atom_update("memory cleanup", na);
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
