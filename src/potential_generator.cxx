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
#include "sho_tools.hxx" // ::quantum_number_table
#include "exchange_correlation.hxx" // ::lda_PZ81_kernel
#include "boundary_condition.hxx" // ::periodic_images
#include "data_view.hxx" // view2D<T>

#include "finite_difference.hxx" // ::stencil_t, ::derive
#include "geometry_analysis.hxx" // ::read_xyz_file, ::fold_back
#include "simple_timer.hxx" // // SimpleTimer
#include "control.hxx" // ::get

#include "boundary_condition.hxx" // Periodic_Boundary, Isolated_Boundary
#include "bessel_transform.hxx" // ::Bessel_j0
#include "debug_tools.hxx" // ::read_from_file

#ifdef DEVEL
    #include "real_space.hxx" // ::Bessel_projection
    #include "lossful_compression.hxx" // print_compressed
    #include "debug_output.hxx" // dump_to_file
    #include "radial_r2grid.hxx" // radial_r2grid_t
    #include "radial_r2grid.hxx" // r2_axis
#endif // DEVEL
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
#include "pw_hamiltonian.hxx" // ::solve
#include "fermi_distribution.hxx" // ::FermiLevel_t

#include "poisson_solver.hxx" // ::solve, ::solver_method
// #include "fourier_poisson.hxx" // ::fourier_solve
// #include "iterative_poisson.hxx" // ::solve

#include "brillouin_zone.hxx" // ::get_kpoint_mesh

#include "data_list.hxx" // data_list<T> // ToDo: replace the std::vector<double*> with new constructions
#include "print_tools.hxx" // print_stats, printf_vector

#define DEBUG
#ifdef  DEBUG
    #include "debug_output.hxx" // here
#else  // DEBUG
    #define here
#endif // DEBUG

namespace potential_generator {
  // this module makes a DFT calculation based on atoms
  // that live in a spherical potential which is found
  // by projecting the 3D potential.
  // their wave functions do not hybridize but they 
  // feel the effect of the density of neighboring atoms


  inline int even(int const any) { return (((any - 1) >> 1) + 1) << 1;}
  inline int n_grid_points(double const suggest) { return (int)even((int)std::ceil(suggest)); }
  
  
  status_t init_geometry_and_grid(
        real_space::grid_t & g // output grid descriptor
      , view2D<double> & xyzZ // output atom coordinates and core charges Z
      , int & natoms // output number of atoms found
      , int const echo=0 // log-level
  ) {
      // SimpleTimer init_function_timer(__FILE__, __LINE__, __func__, echo);
      status_t stat(0);
      
      natoms = 0;
      double cell[3];
      auto const geo_file = control::get("geometry.file", "atoms.xyz");
      int bc[3]; // boundary conditions
      stat += geometry_analysis::read_xyz_file(xyzZ, natoms, geo_file, cell, bc, echo);

      double const h = control::get("potential_generator.grid.spacing", 0.23622);
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
  status_t write_array_to_file(
        char const *filename // file name to write to
      , real_t const array[]  // array pointer
      , int const nx, int const ny, int const nz // grid dimensions
      , int const echo=0 // log-level
      , char const *arrayname="" // some description to appear in the file
  ) {
      char title[128]; std::sprintf(title, "%i x %i x %i  %s", nz, ny, nx, arrayname);
      auto const size = size_t(nz) * size_t(ny) * size_t(nx);
      return dump_to_file(filename, size, array, nullptr, 1, 1, title, echo);
  } // write_array_to_file

  status_t add_smooth_quantities(
        double values[] // add to this function on a 3D grid
      , real_space::grid_t const & g // Cartesian real-space grid descriptor
      , int const na, int32_t const nr2[], float const ar2[] // r2-grid descriptor
      , view2D<double> const & center // (natoms, 4) center coordinates
      , int const n_periodic_images, view2D<double> const & periodic_images
      , double const *const *const atom_qnt // atom data on r2-grids
      , int const echo=0 // log-level
      , int const echo_q=0 // log-level for the charge
      , double const factor=1 // multipliyer
      , char const *quantity="???" // description for log-messages
  ) {
      // add contributions from smooth core densities

      status_t stat(0);
      for(int ia = 0; ia < na; ++ia) {
#ifdef DEVEL
          if (echo > 11) {
              printf("\n## r, %s of atom #%i\n", quantity, ia);
              print_compressed(radial_r2grid::r_axis(nr2[ia], ar2[ia]).data(), atom_qnt[ia], nr2[ia]);
          } // echo
#endif // DEVEL
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
//        if (echo_q > 3) printf("#    00 compensator charge for atom #%d is %g electrons\n", ia, atom_qlm[ia][00]*Y00inv);
#endif // DEVEL
      } // ia
      return stat;
  } // add_smooth_quantities

  
  template <typename real_t>
  void print_direct_projection(
        real_t const array[]
      , real_space::grid_t const &g
      , double const factor=1
      , double const *center=nullptr
  ) {
      // write all values of a grid array to stdout 
      // as function of their distance to a given center

      double cnt[3];
      if (nullptr != center) { 
          set(cnt, 3, center); // copy
      } else {
          for(int d = 0; d < 3; ++d) {
              cnt[d] = 0.5*(g[d] - 1)*g.h[d];
          } // d
      } // center given
      printf("# projection center (relative to grid point (0,0,0) is %g %g %g in units of grid spacings\n",
                cnt[0]*g.inv_h[0], cnt[1]*g.inv_h[1], cnt[2]*g.inv_h[2]);
      for(int iz = 0; iz < g[2]; ++iz) {
          double const z = iz*g.h[2] - cnt[2], z2 = z*z;
          for(int iy = 0; iy < g[1]; ++iy) {
              double const y = iy*g.h[1] - cnt[1], y2 = y*y; 
              for(int ix = 0; ix < g[0]; ++ix) {
                  double const x = ix*g.h[0] - cnt[0], x2 = x*x;
                  double const r = std::sqrt(x2 + y2 + z2);
                  int const izyx = (iz*g[1] + iy)*g[0] + ix;
                  printf("%g %g\n", r*Ang, array[izyx]*factor);
              } // ix
          } // iy
      } // iz
      printf("# radii in %s\n\n", _Ang);
  } // print_direct_projection


  status_t init(
        float const ion=0.f // ionization between first and last atom
      , int const echo=0 // log-level
  ) {
      // compute the self-consistent solution of a single_atom, all states in the core
      // get the spherical core_density and bring it to the 3D grid
      // get the ell=0 compensator charge and add it to the 3D grid
      // envoke exchange_correlation and fourier_poisson
      // add XC and electrostatic potential and zero potential contributions
      // project the total effective potential to each center using bessel_transforms
      // feed back potential shifts into single_atom

      auto const check = int(control::get("check", 0.0));
      auto const run = std::min(std::max(0, 1 - check), 1);

      double constexpr Y00 = solid_harmonics::Y00;
      double constexpr Y00inv = solid_harmonics::Y00inv;
      double constexpr Y00sq = pow2(solid_harmonics::Y00);

      status_t stat(0);
      
      char h_line[32]; set(h_line, 31, '='); h_line[31] = '\0'; // a horizontal line of 31x '='
      if (echo > 0) printf("\n\n# %s\n# Initialize\n# +check=%i run= %i\n# %s\n\n",
                            h_line, check, run, h_line);      

      view2D<double> xyzZ_noconst;
      real_space::grid_t g;
      int na_noconst{0};
      stat += init_geometry_and_grid(g, xyzZ_noconst, na_noconst, echo);
      int const na{na_noconst};
      view2D<double const> const xyzZ(xyzZ_noconst.data(), xyzZ_noconst.stride()); // wrap as (na,4)

      double const cell[3] = {g[0]*g.h[0], g[1]*g.h[1], g[2]*g.h[2]};
     
      std::vector<float> ionization(na, 0.f);
#ifdef DEVEL      
      if ((ion != 0.0) && (na > 1)) {
          if (echo > 2) printf("# %s distribute ionization of %g electrons between first and last atom\n", __func__, ion);
          ionization[0] = ion; ionization[na - 1] = -ionization[0];
          if (echo > 2) {
              printf("# %s ionizations:", __func__);
              printf_vector(" %g", ionization.data(), na);
          } // echo
      } // ionized
#endif // DEVEL

      float const rcut = 32; // radial grids usually end at 9.45 Bohr
      view2D<double> periodic_images;
      int const n_periodic_images = boundary_condition::periodic_images(periodic_images,
                                       cell, g.boundary_conditions(), rcut, echo - 4);
      if (echo > 1) printf("# %s consider %d periodic images\n", __FILE__, n_periodic_images);


      std::vector<double> Za(na);        // list of atomic numbers
      view2D<double> center(na, 4, 0.0); // get memory for a list of atomic centers
      { // scope: prepare atomic coordinates
          double const grid_offset[3] = {0.5*(g[0] - 1)*g.h[0],
                                         0.5*(g[1] - 1)*g.h[1],
                                         0.5*(g[2] - 1)*g.h[2]};
          if (echo > 1) printf("\n# %s List of Atoms: (coordinates in %s)\n", __func__, _Ang);
          for(int ia = 0; ia < na; ++ia) {
              double const Z = xyzZ(ia,3);
              char Symbol[4]; chemical_symbol::get(Symbol, Z, ' ');
              if (echo > 4) printf("# %s  %15.9f %15.9f %15.9f", Symbol,
                              xyzZ(ia,0)*Ang, xyzZ(ia,1)*Ang, xyzZ(ia,2)*Ang);
              Za[ia] = Z;
              for(int d = 0; d < 3; ++d) {
                  center(ia,d) = geometry_analysis::fold_back(xyzZ(ia,d), cell[d]) + grid_offset[d]; // w.r.t. to the center of grid point (0,0,0)
              } // d
              center(ia,3) = 0; // 4th component is not used
#ifdef DEVEL
              if (echo > 6) printf("  relative %g %g %g", center(ia,0)*g.inv_h[0],
                                 center(ia,1)*g.inv_h[1], center(ia,2)*g.inv_h[2]);
#endif // DEVEL
              if (echo > 4) printf("\n");
          } // ia
      } // scope


      std::vector<double> rho_valence(run*g.all(), 0.0);
     
      char const *initial_valence_density_method = control::get("initial.valence.density", "atomic"); // {"atomic", "load", "none"}
      if (echo > 0) printf("\n# initial.valence.density=%s\n", initial_valence_density_method);

      auto const spherical_valence_decay = control::get("atomic.valence.decay", 10.); // after SCF iteration # 10, take_atomic_valence_densities is zero., never if 0
      float take_atomic_valence_densities{0};

      if ('a' == *initial_valence_density_method) { // atomic
          if (echo > 1) printf("# include spherical atomic valence densities in the smooth core densities\n");
          take_atomic_valence_densities = 1; // 100% of the smooth spherical atomic valence densities is included in the smooth core densities
      } else if ('l' == *initial_valence_density_method) { // load
          error("initial.valence.density=load has not been implemented yet");
      } else if ('n' == *initial_valence_density_method) { // none
          warn("initial.valence.density=none may cause problems");
      } else {
          warn("initial.valence.density=%s is unknown, valence density is zero", initial_valence_density_method);
      } // initial_valence_density_method

      
      std::vector<double>  sigma_cmp(na, 1.); // spread of the Gaussian used in the compensation charges
      std::vector<int32_t> numax(na, -1); // init with 0 projectors
      std::vector<int32_t> lmax_qlm(na, -1);
      std::vector<int32_t> lmax_vlm(na, -1);

      // initialize and get sigma, lmax for each atom
      stat += single_atom::atom_update("initialize", na, Za.data(), numax.data(), ionization.data(), (double**)1);
      stat += single_atom::atom_update("lmax qlm",   na, nullptr,    lmax_qlm.data(), &take_atomic_valence_densities);
      stat += single_atom::atom_update("lmax vlm",   na, (double*)1, lmax_vlm.data());
      stat += single_atom::atom_update("sigma cmp",  na, sigma_cmp.data());

      double nve{0}; // determine the number of valence electrons
      if ('a' == (*control::get("valence.electrons", "auto") | 32)) {
          std::vector<double> n_electrons_a(na, 0.); // number of valence electrons brought in by each atom
          stat += single_atom::atom_update("#valence electrons", na, n_electrons_a.data());
          nve = std::accumulate(n_electrons_a.begin(), n_electrons_a.end(), 0.0);
          if (echo > 0) printf("\n# valence.electrons=auto --> %g valence electrons\n\n", nve);
      } else {
          nve = control::get("valence.electrons", 0.0);
          if (echo > 0) printf("\n# valence.electrons=%g\n\n", nve);
      } // auto
      double const n_valence_electrons = nve; // do not use nve beyond this point

      data_list<double> atom_qlm, atom_vlm, atom_rho, atom_mat;
      if (1) { // scope: set up data_list items
          view2D<int32_t> num(4, na, 0); // how many
          for(int ia = 0; ia < na; ++ia) {
              num(0,ia) = pow2(1 + lmax_qlm[ia]); // qlm
              num(1,ia) = pow2(1 + lmax_vlm[ia]); // vlm
              int const ncoeff = sho_tools::nSHO(numax[ia]);
              num(2,ia) =   pow2(ncoeff);         // aDm
              num(3,ia) = 2*pow2(ncoeff);         // aHm+aSm
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

      // allocations of grid quantities
      std::vector<double>  rho(run*g.all()); // [augmented] charge density
      std::vector<double>  Vxc(run*g.all()); // exchange-correlation potential
      std::vector<double>  cmp(run*g.all()); // compensation charge densities
      std::vector<double>  Ves(run*g.all(), 0.0); // electrostatic potential
      std::vector<double> Vtot(run*g.all()); // total effective potential

      char const *es_solver_name = control::get("electrostatic.solver", "multi-grid"); // {"fft", "multi-grid", "MG", "CG", "SD", "none"}
      auto const es_solver_method = poisson_solver::solver_method(es_solver_name);

      // create a FermiLevel object
      fermi_distribution::FermiLevel_t Fermi(n_valence_electrons, 2,
              control::get("electronic.temperature", 1e-3), echo);

      // prepare for solving the Kohn-Sham equation on the real-space grid
      auto const basis_method = control::get("basis", "grid");
      bool const plane_waves = ((*basis_method | 32) == 'p');
      bool const psi_on_grid = ((*basis_method | 32) == 'g');

      // ============================================================================================
      // == prepaprations for the KS Solver on the real-space grid ==================================
      // ============================================================================================
              // create a coarse grid descriptor
              real_space::grid_t gc(g[0]/2, g[1]/2, g[2]/2); // divide the dense grid numbers by two
              gc.set_grid_spacing(cell[0]/gc[0], cell[1]/gc[1], cell[2]/gc[2]); // alternative: 2*g.h[]
              gc.set_boundary_conditions(g.boundary_conditions());
              if (echo > 1) {
                  printf("# use  %d x %d x %d  coarse grid points\n", gc[0], gc[1], gc[2]);
                  printf("# use  %g %g %g  %s  coarse grid spacing\n", gc.h[0]*Ang, gc.h[1]*Ang, gc.h[2]*Ang, _Ang);
              } // echo

              // create a list of atoms
              auto list_of_atoms = grid_operators::empty_list_of_atoms();
              std::vector<double> sigma_a(na, .5);
              { // scope: collect information for projectors and construct a list of atoms
                  { // scope: get the numax for each atom
                      std::vector<int32_t> numax_a(na, 3);
                      stat += single_atom::atom_update("projectors", na, sigma_a.data(), numax_a.data());
                      for(int ia = 0; ia < na; ++ia) {
                          assert( numax_a[ia] == numax[ia] ); // check consistency between atom_update("i") and ("p")
                      } // ia
                  } // scope: numax_a

                  if (psi_on_grid) {
                      view2D<double> xyzZinso(na, 8);
                      for(int ia = 0; ia < na; ++ia) {
                          set(xyzZinso[ia], 4, xyzZ[ia]); // copy
                          xyzZinso(ia,4) = ia;  // global_atom_id
                          xyzZinso(ia,5) = numax[ia];
                          xyzZinso(ia,6) = sigma_a[ia];
                          xyzZinso(ia,7) = 0;   // __not_used__
                      } // ia
                      list_of_atoms = grid_operators::list_of_atoms(xyzZinso.data(), na, xyzZinso.stride(), gc, echo);
                  } // psi_on_grid
              } // scope

              // construct grid-based Hamiltonian and overlap operator descriptor
              using real_wave_function_t = float; // decide here if float or double precision
//            using real_wave_function_t = double;
//            using wave_function_t = std::complex<real_wave_function_t>; // decide here if real or complex
              using wave_function_t = real_wave_function_t;               // decide here if real or complex
              grid_operators::grid_operator_t<wave_function_t, real_wave_function_t> op(gc, list_of_atoms);
              // Mind that local potential and atom matrices of op are still unset!
              list_of_atoms.clear();


              view2D<double> kmesh; // kmesh(nkpoints, 4);
              int const nkpoints = brillouin_zone::get_kpoint_mesh(kmesh);
              if (echo > 0) printf("# k-point mesh has %d points\n", nkpoints);

              double const nbands_per_atom = control::get("bands.per.atom", 10.); // 1: s  4: s,p  10: s,p,ds*  20: s,p,ds*,fp*
              int const nbands = int(nbands_per_atom*na);
              view3D<wave_function_t> psi; // Kohn-Sham states in real-space grid representation
              if (psi_on_grid) { // scope: generate start waves from atomic orbitals
                  auto const start_wave_file = control::get("start.waves", "");
                  psi = view3D<wave_function_t>(run*nkpoints, nbands, gc.all(), wave_function_t(0)); // get memory (potentially large)
                  if (0 == *start_wave_file) {
                      if (echo > 1) printf("# initialize grid wave functions as %d atomic orbitals, %g orbitals per atom\n", nbands, nbands_per_atom);
                      float const scale_sigmas = control::get("start.waves.scale.sigma", 10.); // how much more spread in the start waves compared to sigma_prj
                      uint8_t qn[20][4]; // first 20 sets of quantum numbers [nx, ny, nz, nu] with nu==nx+ny+nz
                      sho_tools::quantum_number_table(qn[0], 3, sho_tools::order_Ezyx); // Ezyx-ordered, take 1, 4, 10 or 20
                      std::vector<int32_t> ncoeff_a(na, 20);
                      data_list<wave_function_t> single_atomic_orbital(ncoeff_a, 0.0); // get memory and initialize
                      for(int iband = 0; iband < nbands; ++iband) {
                          int const ia = iband % na; // which atom?
                          int const io = iband / na; // which orbital?
                          if (io >= 20) error("requested more than 20 start wave functions per atom! bands.per.atom=%g", nbands_per_atom);
                          auto const q = qn[io];
                          if (echo > 7) printf("# initialize band #%i as atomic orbital %x%x%x of atom #%i\n", iband, q[2], q[1], q[0], ia);
                          int const isho = sho_tools::zyx_index(3, q[0], q[1], q[2]); // isho in order_zyx w.r.t. numax=3
                          single_atomic_orbital[ia][isho] = 1./std::sqrt((q[3] > 0) ? ( (q[3] > 1) ? 53. : 26.5 ) : 106.); // set normalization depending on s,p,ds*
                          if (run) op.get_start_waves(psi(0,iband), single_atomic_orbital.data(), scale_sigmas, echo);
                          single_atomic_orbital[ia][isho] = 0; // reset
    //                    if (echo > 17) print_stats(psi(0,iband), gc.all(), gc.dV(), "# single band stats:");
                      } // iband
                  } else {
//                    if (is_complex<wave_function_t>()) error("not prepared for complex wave function reading!");
                      if (echo > 1) printf("# try to read start waves from file \'%s\'\n", start_wave_file);
                      if (run) {
                          auto const nerrors = debug_tools::read_from_file(psi.data(), start_wave_file, nbands, gc.all(), gc.all(), "wave functions", echo);
                          if (nerrors) {
                              error("failed to read start wave functions from file \'%s\'", start_wave_file);
                          } else {
                              if (echo > 1) printf("# read %d bands x %ld numbers from file \'%s\'\n", nbands, gc.all(), start_wave_file);
                          } 
                      } // run
                  } // start wave method
              } // scope

              view2D<double> energies(nkpoints, nbands, 0.0); // Kohn-Sham eigenenergies

              auto const grid_eigensolver_method = control::get("grid.eigensolver", "cg");
              auto const nrepeat = int(control::get("grid.eigensolver.repeat", 1.)); // repetitions of the solver
      // KS solver prepared
      // ============================================================================================
      
      
      if (run != 1) { 
          if (echo > 0) printf("# +check=%i --> done\n", check);
          stat += single_atom::atom_update("memory cleanup", na);
          return stat;
      } // run

      int const max_scf_iterations = control::get("potential_generator.max.scf", 1.);
      for(int scf_iteration = 0; scf_iteration < max_scf_iterations; ++scf_iteration) {
          SimpleTimer scf_iteration_timer(__FILE__, __LINE__, "scf_iteration", echo);
          if (echo > 1) printf("\n\n# %s\n# SCF-iteration step #%i:\n# %s\n\n", h_line, scf_iteration, h_line);

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
                                echo, 0, Y00sq, "smooth core density");
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
                  if (echo > 6) printf("# use generalized Gaussians (sigma= %g %s, lmax=%d) as compensators for atom #%i\n", sigma*Ang, _Ang, ellmax, ia);
                  std::vector<double> coeff(sho_tools::nSHO(ellmax), 0.0);
                  stat += sho_projection::denormalize_electrostatics(coeff.data(), atom_qlm[ia], ellmax, sigma, unitary, echo);
//                if (echo > 7) printf("# before SHO-adding compensators for atom #%i coeff[000] = %g\n", ia, coeff[0]);
                  for(int ii = 0; ii < n_periodic_images; ++ii) {
                      double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                      stat += sho_projection::sho_add(cmp.data(), g, coeff.data(), ellmax, cnt, sigma, 0);
                  } // periodic images
#ifdef DEVEL
                  if (echo > 6) { // report extremal values of the density on the grid
                      printf("# after adding %g electrons compensator density for atom #%i:", atom_qlm[ia][00]*Y00inv, ia);
                      print_stats(cmp.data(), g.all(), g.dV());
                  } // echo
#endif // DEVEL
              } // ia

              // add compensators cmp to rho
              add_product(rho.data(), g.all(), cmp.data(), 1.);
              if (echo > 1) print_stats(rho.data(), g.all(), g.dV(), "\n# augmented charge density:");

              { // scope: solve the Poisson equation: Laplace Ves == -4 pi rho
                  SimpleTimer timer(__FILE__, __LINE__, "Poisson equation", echo);
#ifdef DEVEL
                  if (echo > 0) {
                      printf("\n\n# %s\n# Solve Poisson equation\n# %s\n\n", h_line, h_line);
                      std::fflush(stdout); // if the Poisson solver takes long, we can already see the output up to here
                  } // echo
#endif // DEVEL     
                  stat += poisson_solver::solve(Ves.data(), rho.data(), g, es_solver_method, echo, (na > 0)? center[0] : nullptr);
              } // scope
              here;

#ifdef DEVEL
              { // scope: export electrostatic potential to ASCII file
                  auto const Ves_out_filename = control::get("electrostatic.potential.to.file", "");
                  if (*Ves_out_filename) stat += write_array_to_file(Ves_out_filename, Ves.data(), g[0], g[1], g[2], echo, "electrostatic potential");
//                     dump_to_file(Ves_out_filename, g.all(), Ves.data(), nullptr, 1, 1, "electrostatic potential", echo);
              } // scope
#endif // DEVEL

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
                          double const unitfactor = Y00 * eV * std::pow(Ang, -ell);
                          int const ilm0 = sho_tools::lm_index(ell, -ell);
                          printf("# potential projection for atom #%d v_%im =", ia, ell);
                          printf_vector(" %.6f", &atom_vlm[ia][ilm0], 2*ell + 1, nullptr, unitfactor);
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
#endif // DEVEL
          
          // communicate vlm to the atoms, get zero potential, atom-centered Hamiltonian and overlap
          float potential_mixing_ratio[] = {.5}; // {potential}
          stat += single_atom::atom_update("update", na, 0, 0, potential_mixing_ratio, atom_vlm.data());
          stat += single_atom::atom_update("hamiltonian", na, 0, 0, 0, atom_mat.data());
          stat += single_atom::atom_update("zero potentials", na, 0, nr2.data(), ar2.data(), atom_vbar.data());

          set(Vtot.data(), g.all(), Vxc.data()); add_product(Vtot.data(), g.all(), Ves.data(), 1.);

          if (echo > 1) print_stats(Vtot.data(), g.all(), g.dV(), "\n# Total effective potential (before adding zero potentials)", eV);

          // now also add the zero potential vbar to Vtot
          stat += add_smooth_quantities(Vtot.data(), g, na, nr2.data(), ar2.data(), 
                                center, n_periodic_images, periodic_images, atom_vbar.data(),
                                echo, 0, Y00, "zero potential");

          if (echo > 1) print_stats(Vtot.data(), g.all(), g.dV(), "\n# Total effective potential  (after adding zero potentials)", eV);
          here;

         
          /**  
           *  Potential generation done
           */

          
          std::vector<double> rhov_new(g.all(), 0.0); // new valence density

          { // scope: solve the Kohn-Sham equation with the given Hamiltonian
              SimpleTimer KS_timer(__FILE__, __LINE__, "solving KS-equation", echo);
#ifdef DEVEL
              if (echo > 8) {
                  for(int ia = 0; ia < na; ++ia) {
                      int const n = sho_tools::nSHO(numax[ia]);
                      view2D<double const> const aHm(atom_mat[ia], n);
                      printf("\n# atom-centered %d x %d Hamiltonian (in %s) for atom #%i\n", n, n, _eV, ia);
                      for(int i = 0; i < n; ++i) {
                          printf("#%3i  ", i);
                          printf_vector(" %.6f", aHm[i], n, "\n", eV);
                      } // i
                      view2D<double const> const aSm(atom_mat[ia] + n*n, n);
                      printf("\n# atom-centered %d x %d overlap matrix for atom #%i\n", n, n, ia);
                      for(int i = 0; i < n; ++i) {
                          printf("#%3i  ", i);
                          printf_vector(" %.6f", aSm[i], n);
                      } // i
                  } // ia
              } // echo
#endif // DEVEL

              std::vector<double> rho_valence_new(g.all(), 0.0); // new valence density
              if (psi_on_grid) {
                
                  // restrict the local effective potential to the coarse grid
                  std::vector<double> Veff(gc.all());
                  multi_grid::restrict3D(Veff.data(), gc, Vtot.data(), g, 0); // mute
                  if (echo > 1) print_stats(Veff.data(), gc.all(), gc.dV(), "\n# Total effective potential  (restricted to coarse grid)   ", eV);
#ifdef DEVEL
                  if (0) { // scope: interpolate the effective potential to the dense grid again and compare it to the original version Vtot
                    // in order to test the interpolation routine
                      std::vector<double> v_dcd(g.all(), 0.0);
                      multi_grid::interpolate3D(v_dcd.data(), g, Veff.data(), gc, 0); // mute
                      if (echo > 1) print_stats(v_dcd.data(), g.all(), g.dV(), "\n# Total effective potential (interpolated to dense grid)   ", eV);
                  } // scope

                  if (echo > 0) {
                      if (control::get("potential_generator.direct.projection", 0.) > 0) {
                          double cnt0[3]; if (na > 0) for(int d = 0; d < 3; ++d) cnt0[d] = center(0,d) + 0.5*(g.h[d] - gc.h[d]);
                          printf("\n## all values of Vtot in %s (on the coarse grid, unordered) as function of the distance to %s\n",
                                                            _eV, (na > 0) ? "atom #0" : "the cell center");
                          print_direct_projection(Veff.data(), gc, eV, (na > 0) ? cnt0 : nullptr);
                      } // control
                  } // echo
#endif // DEVEL

                  // copy the local potential and non-local atom matrices into the grid operator descriptor
                  op.set_potential(Veff.data(), gc.all(), atom_mat.data(), echo);

                  std::vector<double> rho_valence_gc(gc.all(), 0.0); // new valence density on the coarse grid

                  double charges[4] = {0, 0, 0, 0};
                  for(int ikpoint = 0; ikpoint < nkpoints; ++ikpoint) { // ToDo: implement k-points
                      auto psi_k = psi[ikpoint]; // get a sub-view

                      // solve the Kohn-Sham equation using various solvers
                      if ('c' == *grid_eigensolver_method) { // "cg" or "conjugate_gradients"
                          stat += davidson_solver::rotate(psi_k.data(), energies[ikpoint], nbands, op, echo);
                          for(int irepeat = 0; irepeat < nrepeat; ++irepeat) {
                              if (echo > 6) { printf("# SCF cycle #%i, CG repetition #%i\n", scf_iteration, irepeat); std::fflush(stdout); }
                              stat += conjugate_gradients::eigensolve(psi_k.data(), energies[ikpoint], nbands, op, echo - 5);
                              stat += davidson_solver::rotate(psi_k.data(), energies[ikpoint], nbands, op, echo);
                          } // irepeat
                      } else
                      if ('d' == *grid_eigensolver_method) { // "davidson"
                          for(int irepeat = 0; irepeat < nrepeat; ++irepeat) {
                              if (echo > 6) { printf("# SCF cycle #%i, DAV repetition #%i\n", scf_iteration, irepeat); std::fflush(stdout); }
                              stat += davidson_solver::eigensolve(psi_k.data(), energies[ikpoint], nbands, op, echo);
                          } // irepeat
                      } else
                      if ('n' == *grid_eigensolver_method) { // "none"
                          if (take_atomic_valence_densities < 1) warn("eigensolver=none generates no new valence density");
                      } else {
                          ++stat; error("unknown grid.eigensolver method \'%s\'", grid_eigensolver_method);
                      } // grid_eigensolver_method

                      // add to density
                      { // scope
                          std::vector<uint32_t> coeff_starts;
                          auto const atom_coeff = density_generator::atom_coefficients(coeff_starts,
                                                    psi_k.data(), gc.all(), na, op, nbands, 1, echo);
                          stat += density_generator::density(rho_valence_gc.data(), atom_rho.data(), Fermi,
                                                    energies[ikpoint], psi_k.data(), atom_coeff.data(), 
                                                    coeff_starts.data(), na, gc, nbands, 1, echo, nullptr, charges);
                      } // scope
                  } // ikpoint
                  if (echo > 2) printf("# %s: total charge %g electrons and derivative %g\n", __func__, 
                                  charges[1]/charges[0], charges[2]/charges[0]*Fermi.get_temperature());

                  stat += multi_grid::interpolate3D(rho_valence_new.data(), g, rho_valence_gc.data(), gc);

              } else if (plane_waves) {
                  std::vector<pw_hamiltonian::DensityIngredients> export_rho;
                  here;

                  stat += pw_hamiltonian::solve(na, xyzZ, g, Vtot.data(), sigma_a.data(), numax.data(), atom_mat.data(), echo, &export_rho);

                  double charges[4] = {0, 0, 0, 0};
                  for(auto & x : export_rho) {
                      if (echo > 1) { printf("\n# Generate valence density for %s\n", x.tag); std::fflush(stdout); }
                      stat += density_generator::density(rho_valence_new.data(), atom_rho.data(), Fermi,
                                                x.energies.data(), x.psi_r.data(), x.coeff.data(), 
                                                x.offset.data(), x.natoms, g, x.nbands, 1, echo, nullptr, charges);
                  } // ikpoint
                  if (echo > 2) printf("# %s: total charge %g electrons and derivative %g\n", __func__, 
                                  charges[1]/charges[0], charges[2]/charges[0]*Fermi.get_temperature());
                  if (charges[0] > 0 && std::abs(charges[0] - 1.0) > 1e-15) {
                      double const sf = 1./charges[0];
                      if (echo > 2) printf("# %s: renormalize density and density matrices by %.15f\n", __func__, sf);
                      scale(rho_valence_new.data(), g.all(), sf);
                      for(int ia = 0; ia < na; ++ia) {
                          scale(atom_rho[ia], pow2(sho_tools::nSHO(numax[ia])), sf);
                      } // ia
                      scale(charges, 3, sf);
                  } // normalize by sum of k-weights

                  here;
              } else { // SHO local orbitals
                  here;

                  stat += sho_hamiltonian::solve(na, xyzZ, g, Vtot.data(), na, sigma_a.data(), numax.data(), atom_mat.data(), echo);
                  warn("with basis=%s no new density is generated", basis_method); // ToDo: implement this

                  here;
              } // psi_on_grid
              
              if (echo > 1) print_stats(rho_valence_new.data(), g.all(), g.dV(), "\n# Total new valence density");

#if 1
              if (1) { // update take_atomic_valence_densities
                  float take_some_spherical_valence_density = 0.5f;
                  stat += single_atom::atom_update("lmax qlm", na, 0, lmax_qlm.data(), &take_some_spherical_valence_density);
              }
              float rho_mixing_ratios[] = {.5, .5, .5}; // {core_density, semicore_density, valence_density}             
              stat += single_atom::atom_update("atomic density matrices", na, 0, 0, rho_mixing_ratios, atom_rho.data());
#endif // 0

          } // scope: Kohn-Sham

          if (spherical_valence_decay > 0) { // update take_atomic_valence_densities
#if 0
              auto const progress = scf_iteration/spherical_valence_decay;
              take_atomic_valence_densities = (progress >= 1) ? 0 : 2*pow3(progress) - 3*pow2(progress) + 1; // smooth transition function
              stat += single_atom::atom_update("lmax qlm", na, 0, lmax_qlm.data(), &take_atomic_valence_densities);
#endif // 0
          } // spherical_valence_decay
          here;
          
          Fermi.correct_Fermi_level(nullptr, echo); // ToDo: pass in charge and its derivative w.r.t. the Fermi level here

      } // scf_iteration
      here;

      if (psi_on_grid) {
          auto const store_wave_file = control::get("store.waves", "");
          if ((nullptr != store_wave_file) && ('\0' != store_wave_file[0])) {
              auto const nerrors = dump_to_file(store_wave_file,
                  nbands, psi.data(), 0, gc.all(), gc.all(), "wave functions", echo);
              if (nerrors) warn("%d errors occured writing file \'%s\'", nerrors, store_wave_file); 
          } // filename not empty
      } // wave functions on Cartesian real-space grid

#ifdef DEVEL

      std::vector<double> Laplace_Ves(g.all(), 0.0);

      auto const verify_Poisson = int(control::get("potential_generator.verify.poisson", 0.));
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

      int const use_Bessel_projection = int(control::get("potential_generator.use.bessel.projection", 0.));
      if (use_Bessel_projection) 
      { // scope: use a Bessel projection around each atom position to compare 3D and radial quantities
        

          std::vector<radial_grid_t const*> rg(na, nullptr); // pointers to smooth radial grid descriptors
          {   // break the interface to get the radial grid descriptors
              auto const dcpp = reinterpret_cast<double const *const *>(rg.data());
              auto const dpp  =       const_cast<double       *const *>(dcpp);
              stat += single_atom::atom_update("radial grids", na, 0, 0, 0, dpp);
          }

          double const* const value_pointers[] = {Ves.data(), Vxc.data(), Vtot.data(), rho.data(), cmp.data(), Laplace_Ves.data()};
          char const *  array_names[] = {"Ves", "Vxc", "Vtot", "rho", "cmp", "LVes"};
          //     Ves; // analyze the electrostatic potential
          //     rho; // analyze the augmented density
          //     Laplace_Ves; // analyze the augmented density computed as Laplacian*Ves
          //     cmp; // analyze only the compensator density
          //     Vxc; // analyze the xc potential
          //     Vtot; // analyze the total potential: Vxc + Ves

          for(int iptr = 0; iptr < std::min(6, use_Bessel_projection); ++iptr) {
              // SimpleTimer timer(__FILE__, __LINE__, "Bessel-projection-analysis", echo);
              auto const values = value_pointers[iptr];
              auto const array_name = array_names[iptr];

              // report extremal values of what is stored on the grid
              if (echo > 1) { printf("\n# real-space stats of %s:", array_name); print_stats(values, g.all(), g.dV()); }

              for(int ia = 0; ia < na; ++ia) {
                  float const dq = 1.f/16;
                  int const nq = int(constants::pi/(g.smallest_grid_spacing()*dq));
                  std::vector<double> qc(nq, 0.0);

                  { // scope: Bessel core
                      std::vector<double> qc_image(nq, 0.0);
                      for(int ii = 0; ii < n_periodic_images; ++ii) {
                          double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                          stat += real_space::Bessel_projection(qc_image.data(), nq, dq, values, g, cnt);
                          add_product(qc.data(), nq, qc_image.data(), 1.0);
                      } // ii
                  } // scope

                  scale(qc.data(), nq, Y00sq);
                  
                  std::vector<double> qcq2(nq, 0.0);
                  for(int iq = 1; iq < nq; ++iq) { // start from 1 to avoid the q=0 term
                      qcq2[iq] = 4*constants::pi*qc[iq]/pow2(iq*dq); // cheap Poisson solver in Bessel transform
                  } // iq

                  if (echo > 11) {
                      printf("\n# Bessel coeff of %s for atom #%d:\n", array_name, ia);
                      for(int iq = 0; iq < nq; ++iq) {
                          printf("# %g %g %g\n", iq*dq, qc[iq], qcq2[iq]);
                      } // iq
                      printf("\n\n");
                  } // echo

                  if (echo > 3) {
                      std::vector<double> rs(rg[ia]->n);
                      bessel_transform::transform_s_function(rs.data(), qc.data(), *rg[ia], nq, dq, true); // transform back to real-space again
                      printf("\n## Real-space projection of %s for atom #%d:\n", array_names[iptr], ia);
                      float const compression_threshold = 1e-4;
                      print_compressed(rg[ia]->r, rs.data(), rg[ia]->n, compression_threshold);

                      if ((values == rho.data()) || (values == Laplace_Ves.data())) {
                          bessel_transform::transform_s_function(rs.data(), qcq2.data(), *rg[ia], nq, dq, true); // transform electrostatic solution to real-space
                          printf("\n## Electrostatics computed by Bessel transform of %s for atom #%d:\n", array_names[iptr], ia);
                          print_compressed(rg[ia]->r, rs.data(), rg[ia]->n, compression_threshold);
                      } // density
                  } // echo
              } // ia

          } // iptr loop for different quantities represented on the grid.

      } // scope Bessel

      if (echo > 1) {
          if (control::get("potential_generator.direct.projection", 0.) > 0) {
              printf("\n## all values of Vtot in %s (unordered) as function of the distance to %s\n",
                                                 _eV, (na > 0) ? "atom #0" : "the cell center");
              print_direct_projection(Vtot.data(), g, eV, (na > 0) ? center[0] : nullptr);          
          } // control
      } // echo

      { // scope: export total potential to ASCII file
          auto const Vtot_out_filename = control::get("total.potential.to.file", "vtot.dat");
          if (*Vtot_out_filename) stat += write_array_to_file(Vtot_out_filename, Vtot.data(), g[0], g[1], g[2], echo);
      } // scope

#endif // DEVEL

      stat += single_atom::atom_update("memory cleanup", na);

      return stat;
  } // init


#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_init(int const echo=3) {
      float const ion = control::get("potential_generator.test.ion", 0.);
      return init(ion, echo); // ionization of Al-P dimer by -ion electrons
  } // test_init

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_init(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace potential_generator
