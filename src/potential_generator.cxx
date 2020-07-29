#include <cstdio> // printf, std::sprintf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::floor
#include <vector> // std::vector
#include <complex> // std::complex

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

#include "boundary_condition.hxx" // Periodic_Boundary, Isolated_Boundary
#include "bessel_transform.hxx" // ::Bessel_j0

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
#include "sho_hamiltonian.hxx" // ::solve

#include "data_list.hxx" // data_list<T> // ToDo: replace the std::vector<double*> with new constructions

#define DEBUG
#ifdef  DEBUG
    #include "debug_output.hxx" // here
#else
    #define here
#endif

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
  
  inline double fold_back(double const position, double const cell_extend) { 
      double x{position};
      while(x >= 0.5*cell_extend) x -= cell_extend;
      while(x < -0.5*cell_extend) x += cell_extend;
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
  
  template <typename real_t>
  status_t write_array_to_file(char const *filename, real_t const array[], 
               int const nx, int const ny, int const nz, int const echo=0, char const *arrayname="") {
      char title[128]; std::sprintf(title, "%i x %i x %i  %s", nz, ny, nx, arrayname);
      auto const size = size_t(nz) * size_t(ny) * size_t(nx);
      return dump_to_file(filename, size, array, nullptr, 1, 1, title, echo);
  } // write_array_to_file

  status_t add_smooth_quantities(double values[] // add to this function on a 3D grid
                , real_space::grid_t const & g 
                , int const na, int32_t const nr2[], float const ar2[]
                , view2D<double> const & center // (natoms, 4)
                , int const n_periodic_images, view2D<double const> const & periodic_images
                , double const *const *const atom_qnt
                , int const echo=0, int const echo_q=0
                , double const factor=1
                , char const *quantity="???") {

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

  
  
  status_t Bessel_Poisson_solver(double *const Ves, real_space::grid_t const & g
                            , double const *const rho, double const center[3]
                            , int const echo=0) {
      if (g.number_of_boundary_conditions(Isolated_Boundary) < 3) {
          warn("Bessel Poisson solver requires 3 isolated boundary conditions");
          return -1;
      } // failed

      status_t stat(0);
        
      // report extremal values of what is stored on the grid
      if (echo > 1) {
          printf("# real-space stats of input density:   ");
          print_stats(rho, g.all(), g.dV());
      } // echo
      
      if (echo > 5) printf("# Bessel j0 projection around position %g %g %g %s\n",
                              center[0]*Ang, center[1]*Ang, center[2]*Ang, _Ang);
      float const dq = 1.f/16; int const nq = int(constants::pi/(g.smallest_grid_spacing()*dq));
      std::vector<double> qc(nq, 0.0);
      stat += real_space::bessel_projection(qc.data(), nq, dq, rho, g, center);
      qc[0] = 0; // stabilize charge neutrality, q=0-component must vanish

      double const by4pi = .25/constants::pi; // 1/(4*pi)

      if (echo > 19) {
          printf("\n## Bessel coeff of {density, Ves}:\n");
          for(int iq = 0; iq < nq; ++iq) {
              double const q = iq*dq;
              printf("%g %g %g\n", q, qc[iq]*by4pi, qc[iq]/(q*q + 1e-12));
          } // iq
          printf("\n\n");
      } // echo

      double const sqrt2pi = std::sqrt(2./constants::pi); // this makes the transform symmetric
      scale(qc.data(), nq, sqrt2pi*dq);

      // expand electrostatic potential onto grid
      double rho_abs{0}, rho_squ{0}, rho_max{0};
      for(int iz = 0; iz < g[2]; ++iz) {
          double const z = iz*g.h[2] - center[2], z2 = z*z;
          for(int iy = 0; iy < g[1]; ++iy) {
              double const y = iy*g.h[1] - center[1], y2 = y*y; 
              for(int ix = 0; ix < g[0]; ++ix) {
                  double const x = ix*g.h[0] - center[0], x2 = x*x;
                  double const r = std::sqrt(x2 + y2 + z2);
                  int const izyx = (iz*g[1] + iy)*g[0] + ix;
                  double ves_r{0}, rho_r{0};
                  for(int iq = 0; iq < nq; ++iq) {
                      double const q = iq*dq;
                      double const bessel_kernel = bessel_transform::Bessel_j0(q*r);
                      ves_r += qc[iq]*bessel_kernel; // cheap Poisson solver
                      rho_r += qc[iq]*bessel_kernel * by4pi * q*q; // reconstruct a spherical density
                  } // iq
                  Ves[izyx] = ves_r; // store
                  // check how much the density deviates from a spherical density
                  double const rho_dif = std::abs(rho[izyx] - rho_r);
                  rho_max = std::max(rho_max, rho_dif);
                  rho_abs += rho_dif;
                  rho_squ += pow2(rho_dif);
              } // ix
          } // iy
      } // iz
      rho_abs *= g.dV();
      rho_squ = std::sqrt(rho_squ)*g.dV();

      if (echo > 1) { 
          printf("# real-space stats of output potential:"); 
          print_stats(Ves, g.all(), g.dV());
          printf("\n");
      } // echo
      
      if (echo > 3) printf("# Bessel_Poisson deviation from a spherical density: "
          "max %.1e a.u., abs %.1e e, rms %.1e e\n", rho_max, rho_abs, rho_squ);

      if (echo > 17) {
          printf("\n## Bessel_Poisson: r, Ves, rho (in a.u.):\n");
          for(int ir = 0; ir < 99; ++ir) {
              double const r = .1*ir;
              double ves_r{0}, rho_r{0};
              for(int iq = 0; iq < nq; ++iq) {
                  double const q = iq*dq;
                  double const bessel_kernel = bessel_transform::Bessel_j0(q*r);
                  ves_r += qc[iq]*bessel_kernel; // cheap Poisson solver
                  rho_r += qc[iq]*bessel_kernel * by4pi * q*q; // reconstruct a spherical density
              } // iq
              printf("%g %g %g \n", r, ves_r, rho_r);
          } // ir
          printf("\n\n");
      } // echo
      
      return stat;
   } // Bessel_Poisson_solver
  
  
  
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
      
      char line[32]; set(line, 31, '='); line[31] = '\0'; // a line of 31x '='
      if (echo > 0) printf("\n\n# %s\n# Initialize\n# %s\n\n", line, line);      
      
      double *coordinates_and_Z{nullptr};
      real_space::grid_t g;
      int na_noconst{0};
      stat += init_geometry_and_grid(g, &coordinates_and_Z, na_noconst, echo);
      int const na{na_noconst};

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
      view2D<double const> periodic_images(periodic_images_ptr, 4); // wrap

      
      std::vector<double> Za(na);        // list of atomic numbers
      view2D<double const> const xyzZ(coordinates_and_Z, 4); // wrap as (na,4)
      view2D<double> center(na, 4, 0.0); // get memory for a list of atomic centers
      { // scope: prepare atomic coordinates
          if (echo > 1) printf("# %s List of Atoms: (coordinates in %s)\n", __func__,_Ang);
          for(int ia = 0; ia < na; ++ia) {
              double const Z = xyzZ(ia,3);
              char Symbol[4]; chemical_symbol::get(Symbol, Z, ' ');
              if (echo > 4) printf("# %s  %15.9f %15.9f %15.9f", Symbol,
                              xyzZ(ia,0)*Ang, xyzZ(ia,1)*Ang, xyzZ(ia,2)*Ang);
              Za[ia] = Z;
              for(int d = 0; d < 3; ++d) {
                  center(ia,d) = fold_back(xyzZ(ia,d), cell[d]) + 0.5*(g[d] - 1)*g.h[d]; // w.r.t. to the center of grid point (0,0,0)
              }   center(ia,3) = 0; // 4th component is not used
              if (echo > 1) printf("  relative%12.3f%16.3f%16.3f", center(ia,0)*g.inv_h[0],
                                          center(ia,1)*g.inv_h[1], center(ia,2)*g.inv_h[2]);
              if (echo > 4) printf("\n");
          } // ia
      } // scope


      std::vector<double> rho_valence(g.all(), 0.0);
     
      char const *initial_valence_density_method = control::get("initial.valence.density", "atomic"); // {"atomic", "load", "none"}
      if (echo > 0) printf("\n# initial.valence.density = \'%s\'\n", initial_valence_density_method);
      
      auto const spherical_valence_decay = control::get("atomic.valence.decay", 10.); // after SCF iteration # 10, take_atomic_valence_densities is zero., never if 0
      float take_atomic_valence_densities = {0};

      if ('a' == *initial_valence_density_method) { // atomic
          if (echo > 1) printf("# include spherical atomic valence densities in the smooth core densities\n");
          take_atomic_valence_densities = 1; // 100% of the smooth spherical atomic valence densities is included in the smooth core densities
          
      } else if ('l' == *initial_valence_density_method) { // load
          error("initial.valence.density = load has not been implemented yet");
      } else if ('n' == *initial_valence_density_method) { // none
          warn("initial.valence.density = none may cause problems");
      } else {
          warn("initial.valence.density = %s is unknown, valence density is zero", initial_valence_density_method);
      } // switch
      
      
      std::vector<double>  sigma_cmp(na, 1.); // spread of the Gaussian used in the compensation charges
      std::vector<int32_t> numax(na, -1); // init with 0 projectors
      std::vector<int32_t> lmax_qlm(na, -1);
      std::vector<int32_t> lmax_vlm(na, -1);

      // for each atom get sigma, lmax
      stat += single_atom::atom_update("initialize", na, Za.data(), numax.data(), ionization.data(), (double**)1);
      stat += single_atom::atom_update("lmax qlm",   na, nullptr,    lmax_qlm.data(), &take_atomic_valence_densities);
      stat += single_atom::atom_update("lmax vlm",   na, (double*)1, lmax_vlm.data());
      stat += single_atom::atom_update("sigma cmp",  na, sigma_cmp.data());

      data_list<double> atom_qlm, atom_vlm, atom_rho, atom_mat;
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

      // the r^2-grids are used to bring radial quantities to a Cartesian grid
      std::vector<int32_t> nr2(na, 1 << 12); // 4096
      std::vector<float>   ar2(na, 16.f); // with nr2 == 4096 rcut = 15.998 Bohr
      data_list<double> atom_vbar(nr2, 0.0); // zero potentials
      data_list<double> atom_rhoc(nr2, 0.0); // core_densities

      int constexpr numax_unitary = 9;
      sho_unitary::Unitary_SHO_Transform<double> const unitary(numax_unitary);

      std::vector<double>  rho(g.all()); // [augmented] charge density
      std::vector<double>  Vxc(g.all()); // exchange-correlation potential
      std::vector<double>  cmp(g.all()); // compensation charge densities
      std::vector<double>  Ves(g.all(), 0.0); // electrostatic potential
      std::vector<double> Vtot(g.all()); // total effective potential

      char const *es_solver_name = control::get("electrostatic.solver", "multi-grid"); // {"fourier", "multi-grid", "CG", "SD", "none"}
      char const es_solver_method = *es_solver_name; // should be one of {'f', 'i', 'n'}

      // prepare for solving the Kohn-Sham equation on the real-space grid
      auto const basis_method = control::get("basis", "grid");
      bool const psi_on_grid = ((basis_method[0] | 32) == 'g');
              // create a coarse grid descriptor
              real_space::grid_t gc(g[0]/2, g[1]/2, g[2]/2); // divide the dense grid numbers by two
              gc.set_grid_spacing(cell[0]/gc[0], cell[1]/gc[1], cell[2]/gc[2]); // alternative: 2*g.h[]
              gc.set_boundary_conditions(g.boundary_conditions());

              // create a list of atoms
              std::vector<atom_image::sho_atom_t> a(0);
              std::vector<double> sigma_a(na, .5);
              { // scope: collect information for projectors and construct a list of atoms
                  std::vector<int32_t> numax_a(na, 3);
                  stat += single_atom::atom_update("projectors", na, sigma_a.data(), numax_a.data());
                  if (psi_on_grid) {
                      view2D<double> xyzZinso(na, 8);
                      for(int ia = 0; ia < na; ++ia) {
                          set(xyzZinso[ia], 4, &coordinates_and_Z[4*ia]); // copy
                          xyzZinso(ia,4) = ia;  // global_atom_id
                          assert( numax_a[ia] == numax[ia] ); // check consistency between atom_update("i") and ("p")
                          xyzZinso(ia,5) = numax_a[ia];
                          xyzZinso(ia,6) = sigma_a[ia];
                          xyzZinso(ia,7) = 0;   // __not_used__
                      } // i
                      stat += grid_operators::list_of_atoms(a, xyzZinso.data(), na, xyzZinso.stride(), gc, echo);
                  } // psi_on_grid
              } // scope

              // construct grid-based Hamiltonian and overlap operator descriptor
              grid_operators::grid_operator_t<double> op(gc, a);
              // Mind that local potential and atom matrices are still unset!


              int const nkpoints = 1; // ToDo
              int const nbands_per_atom = int(control::get("bands.per.atom", 4.)); // s- and p-states
              int const nbands = nbands_per_atom*na;
              view3D<double> psi; // Kohn-Sham states in real-space grid representation
              if (psi_on_grid) { // scope: generate start waves from atomic orbitals
                  psi = view3D<double>(nkpoints, nbands, gc.all()); // get memory
                  float const scale_sigmas = control::get("start.waves.scale.sigma", 5.); // how much more spread in the start waves compared to sigma_prj
                  uint8_t qn[20][4]; // first 20 sets of quantum numbers [nx, ny, nz, nu] with nu==nx+ny+nz
                  sho_tools::construct_index_table<sho_tools::order_Ezyx>(qn, 3); // nu-ordered, take 1, 4, 10 or 20
                  std::vector<int32_t> ncoeff_a(na);
                  for(int ia = 0; ia < na; ++ia) {
                      ncoeff_a[ia] = sho_tools::nSHO(op.get_numax(ia));
                  } // ia
                  data_list<double> single_atomic_orbital(ncoeff_a, 0.0); // get memory and initialize
                  for(int iband = 0; iband < nbands; ++iband) {
                      int const ia = iband % na; // which atom?
                      int const io = iband / na; // which orbital?
                      if (echo > 7) printf("# %s initialize band #%i as atomic orbital %x%x%x of atom #%i\n", 
                                              __func__, iband, qn[io][0], qn[io][1], qn[io][2], ia);
                      int const isho = sho_tools::zyx_index(op.get_numax(ia), qn[io][0], qn[io][1], qn[io][2]);
                      single_atomic_orbital[ia][isho] = 1;
                      op.get_start_waves(psi(0,iband), single_atomic_orbital.data(), scale_sigmas, echo);
                      single_atomic_orbital[ia][isho] = 0;
//                    print_stats(psi(0,iband), gc.all(), gc.dV());
                  } // iband
              } // scope

              view2D<double> energies(nkpoints, nbands, 0.0); // Kohn-Sham eigenenergies

              auto const eigensolver_method = control::get("eigensolver", "cg");
              int const nrepeat = int(control::get("repeat.eigensolver", 1.)); // repetitions of the solver
      // KS solver prepared
      
      

      int const max_scf_iterations = control::get("potential_generator.max.scf", 1.);
      for(int scf_iteration = 0; scf_iteration < max_scf_iterations; ++scf_iteration) {
          SimpleTimer scf_iteration_timer(__FILE__, __LINE__, "scf_iteration", echo);
          if (echo > 1) printf("\n\n# %s\n# SCF-iteration step #%i:\n# %s\n\n", line, scf_iteration, line);

          stat += single_atom::atom_update("x densities", na, 0, nr2.data(), ar2.data(), atom_rhoc.data());
          stat += single_atom::atom_update("qlm charges", na, 0, 0, 0, atom_qlm.data());
          here;

          if (echo > 4) { // report extremal values of the valence density on the grid
              printf("# valence density in iteration #%i:", scf_iteration);
              print_stats(rho_valence.data(), g.all(), g.dV());
          } // echo
          
          set(rho.data(), g.all(), rho_valence.data(), 1. - take_atomic_valence_densities);

          // add contributions from smooth core densities, and optionally spherical valence densities
          stat += add_smooth_quantities(rho.data(), g, na, nr2.data(), ar2.data(), 
                                center, n_periodic_images, periodic_images, atom_rhoc.data(),
                                echo, echo, Y00sq, "smooth core density");
          here;

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
          here;

          set(cmp.data(), g.all(), 0.0); // init compensation charge density
          { // scope: solve the Poisson equation

              // add compensation charges cmp
              for(int ia = 0; ia < na; ++ia) {
                  double const sigma = sigma_cmp[ia];
                  int    const ellmax = lmax_qlm[ia];
                  if (echo > -1) printf("# use generalized Gaussians (sigma= %g %s, lmax=%d) as compensators for atom #%i\n", sigma*Ang, _Ang, ellmax, ia);
                  std::vector<double> coeff(sho_tools::nSHO(ellmax), 0.0);
                  stat += sho_projection::denormalize_electrostatics(coeff.data(), atom_qlm[ia], ellmax, sigma, unitary, echo);
//                if (echo > 7) printf("# before SHO-adding compensators for atom #%i coeff[000] = %g\n", ia, coeff[0]);
                  for(int ii = 0; ii < n_periodic_images; ++ii) {
                      double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                      stat += sho_projection::sho_add(cmp.data(), g, coeff.data(), ellmax, cnt, sigma, 0);
                  } // periodic images
#ifdef DEVEL
                  if (echo > 1) { // report extremal values of the density on the grid
                      printf("# after adding %g electrons compensator density for atom #%i:", atom_qlm[ia][00]*Y00inv, ia);
                      print_stats(cmp.data(), g.all(), g.dV());
                  } // echo
#endif
              } // ia

              // add compensators cmp to rho
              add_product(rho.data(), g.all(), cmp.data(), 1.);
              if (echo > 1) {
                  printf("\n# augmented charge density:");
                  print_stats(rho.data(), g.all(), g.dV());
              } // echo

              
              { // scope: solve the Poisson equation: Laplace Ves == -4 pi rho
                  SimpleTimer timer(__FILE__, __LINE__, "Poisson equation", echo);
#ifdef DEVEL
                  if (echo > 0) {
                      printf("\n\n# %s\n# Solve Poisson equation\n# %s\n\n", line, line);
                      fflush(stdout); // flush stdout so if the Poisson solver takes long, we can already see the output up to here
                  } // echo
#endif
              
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
#ifdef DEVEL
                      // create a 2x denser grid descriptor
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
#else
                      error("electrostatic.solver = %s only available in the development branch!", es_solver_name); 
#endif
                  } else if ('B' == es_solver_method) { // "Bessel0" (upper case!)
#ifdef DEVEL
                      if (echo > 0) printf("# use a spherical Bessel solver for the Poisson equation\n");
                      assert(1 == na); // must be exactly one atom
                      auto const st = Bessel_Poisson_solver(Ves.data(), g, rho.data(), center[0], echo);
                      if (st) warn("Bessel Poisson solver failed with status=%i", int(st));
                      stat += st;
#else
                      error("electrostatic.solver = %s only available in the development branch!", es_solver_name); 
#endif
                  } else if ('l' == (es_solver_method | 32)) { // "load"
#ifdef DEVEL
                      auto const Ves_in_filename = control::get("electrostatic.potential.from.file", "v_es.dat");
                      auto const nerrors = debug_tools::read_from_file(Ves.data(), Ves_in_filename, g.all(), 1, 1, "electrostatic potential", echo);
                      if (nerrors) warn("electrostatic.solver = %s from file %s had %d errors", es_solver_name, Ves_in_filename, nerrors); 
#else
                      error("electrostatic.solver = %s only available in the development branch!", es_solver_name); 
#endif
                  } else if ('n' == (es_solver_method | 32)) { // "none"
                      warn("electrostatic.solver = %s may lead to unphysical results!", es_solver_name); 

                  } else { // default
                      if (echo > 2) printf("# electrostatic.solver = %s\n", es_solver_name);
                      stat += iterative_poisson::solve(Ves.data(), rho.data(), g, es_solver_method, echo);
                  } // es_solver_method
                  
                  if (echo > 1) {
                      printf("\n# electrostatic potential");
                      print_stats(Ves.data(), g.all(), g.dV());
                  } // echo
                  
              } // scope
              here;

#ifdef DEVEL
              { // scope: export electrostatic potential to ASCII file
                  auto const Ves_out_filename = control::get("electrostatic.potential.to.file", "");
                  if (*Ves_out_filename) stat += write_array_to_file(Ves_out_filename, Ves.data(), g[0], g[1], g[2], echo, "electrostatic potential");
//                     dump_to_file(Ves_out_filename, g.all(), Ves.data(), nullptr, 1, 1, "electrostatic potential", echo);
              } // scope
#endif

              // probe the electrostatic potential in real space, find ves_multipoles
              for(int ia = 0; ia < na; ++ia) {
                  double const sigma = sigma_cmp[ia];
                  int    const ellmax = lmax_vlm[ia];
                  int const nc = sho_tools::nSHO(ellmax);
                  std::vector<double> coeff(nc, 0.0);
                  for(int ii = 0; ii < n_periodic_images; ++ii) {
                      std::vector<double> coeff_image(nc, 0.0);
                      double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                      stat += sho_projection::sho_project(coeff_image.data(), ellmax, cnt, sigma, Ves.data(), g, 0);
                      // alternative to projecting Ves we could use Vtot (constructed later) but we need consistency with the operations inside the spheres
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
                  printf("# inner product between rho_aug and Ves = %g %s\n", 2*Ees*eV,_eV);
              } // echo
          } // scope
          here;

#ifdef DEVEL
//        exit(__LINE__);
#endif
          
          // communicate vlm to the atoms, get zero potential, atom-centered Hamiltonian and overlap
          float mixing_ratios[] = {.5, .5, .5, .5}; // {potential, core_density, semicore_density, valence_density}
          stat += single_atom::atom_update("update", na, 0, 0, mixing_ratios, atom_vlm.data());
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
          here;

         
          /**  
           *  Potential generation done
           */

          
          std::vector<double> rhov_new(g.all(), 0.0); // new valence density

          { // scope: solve the Kohn-Sham equation with the given Hamiltonian
              SimpleTimer KS_timer(__FILE__, __LINE__, "Solve KS-equation", echo);
#ifdef DEVEL
              if (echo > 6) {
                  for(int ia = 0; ia < na; ++ia) {
                      int const n = sho_tools::nSHO(numax[ia]);
                      view2D<double const> const aHm(atom_mat[ia], n);
                      printf("\n# atom-centered %dx%d Hamiltonian (in %s) for atom index #%i\n", n, n, _eV, ia);
                      for(int i = 0; i < n; ++i) {
                          printf("#%3i  ", i);
                          for(int j = 0; j < n; ++j) {
                              printf(" %.3f", aHm(i,j)*eV);
                          }   printf("\n");
                      } // i
                  } // ia
              } // echo
#endif 

              if (psi_on_grid) {
                
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

                  // copy the local potential and non-local atom matrices into the grid operator descriptor
                  op.set_potential(Veff.data(), gc.all(), atom_mat.data(), echo);

                  std::vector<double> rho_valence_new(gc.all(), 0.0); // new valence density

                  for(int ikpoint = 0; ikpoint < nkpoints; ++ikpoint) { // ToDo: implement k-points
                      auto psi_k = psi[ikpoint]; // get a sub-view

                      // solve the Kohn-Sham equation using various solvers
                      if ('c' == *eigensolver_method) { // "cg" or "conjugate_gradients"
                          stat += davidson_solver::rotate(psi_k.data(), energies[ikpoint], nbands, op, echo);
                          for(int irepeat = 0; irepeat < nrepeat; ++irepeat) {
                              stat += conjugate_gradients::eigensolve(psi_k.data(), nbands, op, echo - 5, 1e-6, energies[ikpoint]);
                              stat += davidson_solver::rotate(psi_k.data(), energies[ikpoint], nbands, op, echo);
                          } // irepeat
                      } else
                      if ('d' == *eigensolver_method) { // "davidson"
                          for(int irepeat = 0; irepeat < nrepeat; ++irepeat) {
                              stat += davidson_solver::eigensolve(psi_k.data(), energies[ikpoint], nbands, op, echo + 9);
                          } // irepeat
                      } else
                      if ('n' == *eigensolver_method) { // "none"
                          if (take_atomic_valence_densities < 1) warn("eigensolver=\'none\' generates no new valence density");
                      } else {
                          ++stat; error("unknown eigensolver method \'%s\'", eigensolver_method);
                      } // eigensolver_method

                      // add to density
                      stat += density_generator::density(rho_valence_new.data(), atom_rho.data(), psi_k.data(), op, nbands, 1, echo);

                  } // ikpoint

                  // generate a new density from the eigenstates
                  // with occupation numbers from eigenenergies

                  stat += multi_grid::interpolate3D(rho_valence.data(), g, rho_valence_new.data(), gc);
                  
              } else { // psi_on_grid
                  here;
                
                  stat += sho_hamiltonian::solve<std::complex<double>>(na, xyzZ, g, Vtot.data(), na, sigma_a.data(), numax.data(), atom_mat.data(), echo);
                
                  here;
              } // psi_on_grid

              // ToDo: density mixing
              
              stat += single_atom::atom_update("atomic density matrices", na, 0, 0, 0, atom_rho.data());              

          } // scope: Kohn-Sham

          if (spherical_valence_decay > 0) { // update take_atomic_valence_densities
              auto const progress = scf_iteration/spherical_valence_decay;
              take_atomic_valence_densities = (progress >= 1) ? 0 : 2*pow3(progress) - 3*pow2(progress) + 1; // smooth transition function
              stat += single_atom::atom_update("lmax qlm", na, 0, lmax_qlm.data(), &take_atomic_valence_densities);
          } // spherical_valence_decay
          
          here;
      } // scf_iteration
      here;

#ifdef DEVEL

      std::vector<double> Laplace_Ves(g.all(), 0.0);

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
        

          std::vector<radial_grid_t const*> rg(na, nullptr); // pointers to smooth radial grid descriptors
          {   // break the interface to get the radial grid descriptors
              auto const dcpp = reinterpret_cast<double const *const *>(rg.data());
              auto const dpp  =       const_cast<double       *const *>(dcpp);
              stat += single_atom::atom_update("radial grids", na, 0, 0, 0, dpp);
          }

          double* const value_pointers[] = {Ves.data(), Vxc.data(), Vtot.data(), rho.data(), cmp.data(), Laplace_Ves.data()};
          char const *  value_names[] = {"Ves", "Vxc", "Vtot", "rho", "cmp", "LVes"};
          //     Ves; // analyze the electrostatic potential
          //     rho; // analyze the augmented density
          //     Laplace_Ves; // analyze the augmented density computed as Laplacian*Ves
          //     cmp; // analyze only the compensator density
          //     Vxc; // analyze the xc potential
          //     Vtot; // analyze the total potential: Vxc + Ves

          for(int iptr = 0; iptr < std::min(6, use_Bessel_projection); ++iptr) {
              // SimpleTimer timer(__FILE__, __LINE__, "Bessel-projection-analysis", echo);
              auto const values = value_pointers[iptr];

              // report extremal values of what is stored on the grid
              if (echo > 1) { printf("\n# real-space stats of %s:", value_names[iptr]); print_stats(values, g.all(), g.dV()); }

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
                  for(int iq = 1; iq < nq; ++iq) { // start from 1 to avoid the q=0 term
                      qcq2[iq] = 4*constants::pi*qc[iq]/pow2(iq*dq); // cheap Poisson solver in Bessel transform
                  } // iq

                  if (echo > 11) {
                      printf("\n# Bessel coeff of %s for atom #%d:\n", value_names[iptr], ia);
                      for(int iq = 0; iq < nq; ++iq) {
                          printf("# %g %g %g\n", iq*dq, qc[iq], qcq2[iq]);
                      }   printf("\n\n");
                  } // echo

                  if (echo > 3) {
                      std::vector<double> rs(rg[ia]->n);
                      bessel_transform::transform_s_function(rs.data(), qc.data(), *rg[ia], nq, dq, true); // transform back to real-space again
                      printf("\n## Real-space projection of %s for atom #%d:\n", value_names[iptr], ia);
                      float const compression_threshold = 1e-4;
                      print_compressed(rg[ia]->r, rs.data(), rg[ia]->n, compression_threshold);

                      if ((values == rho.data()) || (values == Laplace_Ves.data())) {
                          bessel_transform::transform_s_function(rs.data(), qcq2.data(), *rg[ia], nq, dq, true); // transform electrostatic solution to real-space
                          printf("\n## Electrostatics computed by Bessel transform of %s for atom #%d:\n", value_names[iptr], ia);
                          print_compressed(rg[ia]->r, rs.data(), rg[ia]->n, compression_threshold);
                      } // density
                  } // echo
              } // ia

          } // iptr loop for different quantities represented on the grid.

      } // scope Bessel

//       int const export_Vtot = int(control::get("potential_generator.export.vtot", 1.));
//       if (export_Vtot) {
//           generate a file which contains the full potential Vtot
//           if (echo > 0) {
//               char title[96]; std::sprintf(title, "%i x %i x %i", g[2], g[1], g[0]);
//               dump_to_file("vtot.dat", Vtot.size(), Vtot.data(), nullptr, 1, 1, title, echo);
//           } // unless all output is suppressed
//       } // export_Vtot

      { // scope: export total potential to ASCII file
          auto const Vtot_out_filename = control::get("total.potential.to.file", "vtot.dat");
          if (*Vtot_out_filename) stat += write_array_to_file(Vtot_out_filename, Vtot.data(), g[0], g[1], g[2], echo, "total effective potential");
      } // scope

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
      float const ion = control::get("potential_generator.test.ion", 0.);
      return init(ion, echo); // ionization of Al-P dimer by -ion electrons
  } // test_init

  status_t all_tests(int const echo) {
    status_t status(0);
    status += test_init(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace potential_generator
