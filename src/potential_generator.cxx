#include <cstdio> // std::printf, std::sprintf
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
#include "exchange_correlation.hxx" // ::LDA_kernel
#include "boundary_condition.hxx" // ::periodic_images
#include "data_view.hxx" // view2D<T>

#include "finite_difference.hxx" // ::stencil_t, ::derive
#include "geometry_analysis.hxx" // ::read_xyz_file, ::fold_back
#include "simple_timer.hxx" // // SimpleTimer
#include "control.hxx" // ::get

#include "boundary_condition.hxx" // Periodic_Boundary, Isolated_Boundary
#include "bessel_transform.hxx" // ::Bessel_j0
#include "debug_tools.hxx" // ::read_from_file, ::manage_stop_file
#include "debug_output.hxx" // dump_to_file

#ifdef DEVEL
    #include "real_space.hxx" // ::Bessel_projection
    #include "lossful_compression.hxx" // print_compressed
    #include "radial_r2grid.hxx" // radial_r2grid_t
    #include "radial_r2grid.hxx" // r2_axis
#endif // DEVEL
#include "sho_unitary.hxx" // ::Unitary_SHO_Transform<real_t>

#include "single_atom.hxx" // ::atom_update
#include "energy_contribution.hxx" // ::TOTAL, ::KINETIC, ::ELECTROSTATIC, ...

// ToDo: restructure: move this into a separate compilation unit
#include "atom_image.hxx"// ::sho_atom_t
#include "grid_operators.hxx" // ::grid_operator_t, ::list_of_atoms
#include "conjugate_gradients.hxx" // ::eigensolve
#include "davidson_solver.hxx" // ::rotate, ::eigensolve
#include "multi_grid.hxx" // ::restrict3D, ::interpolate3D
#include "density_generator.hxx" // ::density
#include "sho_hamiltonian.hxx" // ::solve
#include "plane_waves.hxx" // ::solve
#include "fermi_distribution.hxx" // ::FermiLevel_t
#include "unit_system.hxx" // ::length_unit

#include "poisson_solver.hxx" // ::solve, ::solver_method

#include "brillouin_zone.hxx" // ::get_kpoint_mesh

#include "data_list.hxx" // data_list<T> // ToDo: replace the std::vector<double*> with new constructions
#include "print_tools.hxx" // print_stats, printf_vector
#include "dense_solver.hxx" // ::display_spectrum
#include "dense_solver.hxx" // ::solve
#include "complex_tools.hxx" // complex_name

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
      status_t stat(0);

      int const check = control::get("check", 0.);
      int const run = std::min(std::max(0, 1 - check), 1);

      double constexpr Y00 = solid_harmonics::Y00;
      double constexpr Y00sq = pow2(Y00);
#ifdef DEVEL
      double constexpr Y00inv = solid_harmonics::Y00inv;
#endif // DEVEL

      char h_line[32]; set(h_line, 31, '='); h_line[31] = '\0'; // a horizontal line of 31x '='
      if (echo > 0) std::printf("\n\n# %s\n# Initialize\n# +check=%i run= %i\n# %s\n\n",
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
      if ((0 != ion) && (na > 1)) {
          if (echo > 2) std::printf("# %s distribute ionization of %g electrons between first and last atom\n", __func__, ion);
          ionization[0] = ion; ionization[na - 1] = -ionization[0];
          if (echo > 2) {
              std::printf("# %s ionizations:", __func__);
              printf_vector(" %g", ionization.data(), na);
          } // echo
      } // ionized
#endif // DEVEL

      float const rcut = 32; // radial grids usually end at 9.45 Bohr
      view2D<double> periodic_images;
      int const n_periodic_images = boundary_condition::periodic_images(periodic_images,
                                       cell, g.boundary_conditions(), rcut, echo - 4);
      if (echo > 1) std::printf("# %s consider %d periodic images\n", __FILE__, n_periodic_images);


      std::vector<double> Za(na);        // list of atomic numbers
      view2D<double> center(na, 4, 0.0); // get memory for a list of atomic centers
      { // scope: prepare atomic coordinates
          double const grid_offset[3] = {0.5*(g[0] - 1)*g.h[0],
                                         0.5*(g[1] - 1)*g.h[1],
                                         0.5*(g[2] - 1)*g.h[2]};
          if (echo > 1) std::printf("\n# %s List of Atoms: (coordinates in %s)\n", __func__, _Ang);
          for(int ia = 0; ia < na; ++ia) {
              double const Z = xyzZ(ia,3);
              char Symbol[4]; chemical_symbol::get(Symbol, Z, ' ');
              if (echo > 4) std::printf("# %s  %15.9f %15.9f %15.9f", Symbol,
                              xyzZ(ia,0)*Ang, xyzZ(ia,1)*Ang, xyzZ(ia,2)*Ang);
              Za[ia] = Z;
              for(int d = 0; d < 3; ++d) {
                  center(ia,d) = geometry_analysis::fold_back(xyzZ(ia,d), cell[d]) + grid_offset[d]; // w.r.t. to the center of grid point (0,0,0)
              } // d
              center(ia,3) = 0; // 4th component is not used
#ifdef DEVEL
              if (echo > 6) std::printf("  relative %g %g %g", center(ia,0)*g.inv_h[0],
                                 center(ia,1)*g.inv_h[1], center(ia,2)*g.inv_h[2]);
#endif // DEVEL
              if (echo > 4) std::printf("\n");
          } // ia
      } // scope


      std::vector<double> rho_valence(run*g.all(), 0.0);

      char const *initial_valence_density_method = control::get("initial.valence.density", "atomic"); // {"atomic", "load", "none"}
      if (echo > 0) std::printf("\n# initial.valence.density=%s\n", initial_valence_density_method);

      auto const spherical_valence_decay = control::get("atomic.valence.decay", 10.); // after SCF iteration # 10, take_atomic_valence_densities is zero, never if 0
      float take_atomic_valence_densities{0};

      if ('a' == *initial_valence_density_method) { // atomic
          if (echo > 1) std::printf("# include spherical atomic valence densities in the smooth core densities\n");
          take_atomic_valence_densities = 1; // 100% of the smooth spherical atomic valence densities is included in the smooth core densities
      } else if ('l' == *initial_valence_density_method) { // load
          error("initial.valence.density=load has not been implemented yet, found %s", initial_valence_density_method);
      } else if ('n' == *initial_valence_density_method) { // none
          warn("initial.valence.density=none may cause problems, found %s", initial_valence_density_method);
      } else {
          warn("initial.valence.density=%s is unknown, valence density is zero", initial_valence_density_method);
      } // initial_valence_density_method


      // how to solve the KS-equation
      auto const basis_method = control::get("basis", "grid");
      bool const plane_waves = ('p' == (*basis_method | 32));
      bool const psi_on_grid = ('g' == (*basis_method | 32));

      // ============================================================================================
      // prepare for solving the Kohn-Sham equation on the real-space grid
      // create a coarse grid descriptor (this could be done later but we want a better overview in the log file)
      real_space::grid_t gc(g[0]/2, g[1]/2, g[2]/2); // divide the dense grid numbers by two
      gc.set_grid_spacing(cell[0]/gc[0], cell[1]/gc[1], cell[2]/gc[2]); // alternative: 2*g.h[]
      gc.set_boundary_conditions(g.boundary_conditions());
      if (psi_on_grid && echo > 1) {
          std::printf("# use  %d x %d x %d  coarse grid points\n", gc[0], gc[1], gc[2]);
          double const max_grid_spacing = std::max(std::max(gc.h[0], gc.h[1]), gc.h[2]);
          std::printf("# use  %g %g %g  %s  coarse grid spacing, corresponds to %.2f Ry\n",
                    gc.h[0]*Ang, gc.h[1]*Ang, gc.h[2]*Ang, _Ang, pow2(constants::pi/max_grid_spacing));
      } // echo
      // ============================================================================================
      
      
      
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
          if (echo > 0) std::printf("\n# valence.electrons=auto --> %g valence electrons\n\n", nve);
      } else {
          nve = control::get("valence.electrons", 0.0);
          if (echo > 0) std::printf("\n# valence.electrons=%g\n\n", nve);
      } // auto
      double const n_valence_electrons = nve; // do not use nve beyond this point

      std::vector<int32_t> n_atom_rho(na, 0);
      data_list<double> atom_qlm, atom_vlm, atom_rho, atom_mat;
      if (1) { // scope: set up data_list items
          view2D<int32_t> num(4, na, 0); // how many
          for(int ia = 0; ia < na; ++ia) {
              num(0,ia) = pow2(1 + lmax_qlm[ia]); // qlm
              num(1,ia) = pow2(1 + lmax_vlm[ia]); // vlm
              int const ncoeff = sho_tools::nSHO(numax[ia]);
              num(2,ia) =   pow2(ncoeff);         // aDm
              num(3,ia) = 2*pow2(ncoeff);         // aHm+aSm
              n_atom_rho[ia] = num(2,ia); // number of coefficients for the atomic density matrices
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

      char const occupation_method = *control::get("fermi.level", "exact"); // {"exact", "linearized"}

      // create a FermiLevel object
      fermi_distribution::FermiLevel_t Fermi(n_valence_electrons, 2,
              control::get("electronic.temperature", 1e-3), echo);

      double density_mixing{1};

      // ============================================================================================
      // == preparations for the KS Solver on the real-space grid ===================================
      // ============================================================================================

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
                          set(xyzZinso[ia], 4, xyzZ[ia]); // copy x,y,z,Z
                          xyzZinso(ia,4) = ia;  // global_atom_id
                          xyzZinso(ia,5) = numax[ia];
                          xyzZinso(ia,6) = sigma_a[ia];
                          xyzZinso(ia,7) = 0;   // __not_used__
                      } // ia
                      list_of_atoms = grid_operators::list_of_atoms(xyzZinso.data(), na, xyzZinso.stride(), gc, echo);
                  } // psi_on_grid
              } // scope

              // construct grid-based Hamiltonian and overlap operator descriptor
//               using real_wave_function_t = float;  // decide here if single or double precision
              using real_wave_function_t = double; // decide here if single or double precision
              using wave_function_t = std::complex<real_wave_function_t>; // decide here if real or complex
//               using wave_function_t = real_wave_function_t;               // decide here if real or complex
              grid_operators::grid_operator_t<wave_function_t, real_wave_function_t> op(gc, list_of_atoms);
              // Mind that local potential and atom matrices of op are still unset!
              list_of_atoms.clear();
              if (psi_on_grid) {
                  if (echo > 0) std::printf("# real-space grid wave functions are of type %s\n", complex_name<wave_function_t>());
                  auto const floating_point_bits = int(control::get("hamiltonian.floating.point.bits", 64.)); // double by default
                  if ((floating_point_bits == 32) != (sizeof(real_wave_function_t) == 4)) {
                      warn("hamiltonian.floating.point.bits=%d but wave_function_t is %s", 
                          floating_point_bits, complex_name<wave_function_t>());
                  } // requested precision does not match chose precision
              } // psi_on_grid

              view2D<double> kmesh; // kmesh(nkpoints, 4);
              int const nkpoints = brillouin_zone::get_kpoint_mesh<is_complex<wave_function_t>()>(kmesh);
              if (echo > 0) std::printf("# k-point mesh has %d points\n", nkpoints);
              // ToDo: warn if boundary_condition is isolated but there is more than 1 kpoint

              double const nbands_per_atom = control::get("bands.per.atom", 10.); // 1: s  4: s,p  10: s,p,ds*  20: s,p,ds*,fp*
              int const nbands = int(nbands_per_atom*na);
              view3D<wave_function_t> psi; // Kohn-Sham states in real-space grid representation
              if (psi_on_grid) { // scope: generate start waves from atomic orbitals
                  auto const start_wave_file = control::get("start.waves", "");
                  psi = view3D<wave_function_t>(run*nkpoints, nbands, gc.all(), 0.0); // get memory (potentially large)
                  if ('\0' == *start_wave_file) {
                      if (echo > 1) std::printf("# initialize grid wave functions as %d atomic orbitals, %g orbitals per atom\n", nbands, nbands_per_atom);
                      float const scale_sigmas = control::get("start.waves.scale.sigma", 10.); // how much more spread in the start waves compared to sigma_prj
                      uint8_t qn[20][4]; // first 20 sets of quantum numbers [nx, ny, nz, nu] with nu==nx+ny+nz
                      sho_tools::quantum_number_table(qn[0], 3, sho_tools::order_Ezyx); // Ezyx-ordered, take 1, 4, 10 or 20
                      std::vector<int32_t> ncoeff_a(na, 20);
                      data_list<wave_function_t> single_atomic_orbital(ncoeff_a, 0.0); // get memory and initialize
                      for(int ikpoint = 0; ikpoint < nkpoints; ++ikpoint) {
                          op.set_kpoint(kmesh[ikpoint], echo);
                          for(int iband = 0; iband < nbands; ++iband) {
                              int const iatom    = iband % na; // which atom?
                              int const iorbital = iband / na; // which orbital?
                              if (iorbital >= 20) error("requested more than 20 start wave functions per atom! bands.per.atom=%g", nbands_per_atom);
                              auto const q = qn[iorbital];
                              if (echo > 7) std::printf("# initialize band #%i as atomic orbital %x%x%x of atom #%i\n", iband, q[2], q[1], q[0], iatom);
                              int const isho = sho_tools::zyx_index(3, q[0], q[1], q[2]); // isho in order_zyx w.r.t. numax=3
                              single_atomic_orbital[iatom][isho] = 1./std::sqrt((q[3] > 0) ? ( (q[3] > 1) ? 53. : 26.5 ) : 106.); // set normalization depending on s,p,ds*
                              if (run) op.get_start_waves(psi(ikpoint,iband), single_atomic_orbital.data(), scale_sigmas, echo);
                              single_atomic_orbital[iatom][isho] = 0; // reset
                          } // iband
                      } // ikpoint
                      op.set_kpoint(); // reset k-point
                  } else {
                      if (echo > 1) std::printf("# try to read start waves from file \'%s\'\n", start_wave_file);
                      if (run) {
                          auto const nerrors = debug_tools::read_from_file(psi(0,0), start_wave_file, nbands, psi.stride(), gc.all(), "wave functions", echo);
                          if (nerrors) {
                              error("failed to read start wave functions from file \'%s\'", start_wave_file);
                          } else {
                              if (echo > 1) std::printf("# read %d bands x %ld numbers from file \'%s\'\n", nbands, gc.all(), start_wave_file);
                          }
                      } // run
                      for(int ikpoint = 1; ikpoint < nkpoints; ++ikpoint) {
                          if (echo > 3) { std::printf("# copy %d bands for k-point #%i from k-point #0\n", nbands, ikpoint); std::fflush(stdout); }
                          if (run) set(psi(ikpoint,0), psi.dim1()*psi.stride(), psi(0,0)); // copy, ToDo: include Bloch phase factors
                      } // ikpoints
                  } // start wave method

              } // scope

              view2D<double> energies(nkpoints, nbands, 0.0); // Kohn-Sham eigenenergies

              auto const grid_eigensolver_method = control::get("grid.eigensolver", "cg");
              auto const nrepeat = int(control::get("grid.eigensolver.repeat", 1.)); // repetitions of the solver
      // KS solver on grid prepared
      // ============================================================================================
      

      // total energy contributions
      double grid_electrostatic_energy{0}, grid_kinetic_energy{0}, grid_xc_energy{0}, total_energy{0};
              
      if (run != 1) { 
          if (echo > 0) std::printf("# +check=%i --> done\n", check);
          stat += single_atom::atom_update("memory cleanup", na);
          return stat;
      } // run

      
      int const max_scf_iterations_input = control::get("potential_generator.max.scf", 1.);
      int max_scf_iterations{max_scf_iterations_input}; // non-const, may be modified by the stop file during the run
      { // scope: create a stop file with standard name
          auto const stat = debug_tools::manage_stop_file<'w'>(max_scf_iterations, echo);
          if (stat != 0) error("failed to write/create a stop file, status= %i", int(stat));
      } // scope

      for(int scf_iteration = 0; scf_iteration < max_scf_iterations; ++scf_iteration) {
          SimpleTimer scf_iteration_timer(__FILE__, __LINE__, "scf_iteration", echo);
          if (echo > 1) std::printf("\n\n# %s\n# SCF-iteration step #%i:\n# %s\n\n", h_line, scf_iteration, h_line);

          if (take_atomic_valence_densities > 0) {
              if (echo > 4) print_stats(rho_valence.data(), g.all(), g.dV(), "# previous valence density");

              if (echo > 4) std::printf("# compose valence density with %g %% of the atomic valence densities\n", take_atomic_valence_densities*100);
              scale(rho_valence.data(), g.all(), 1. - take_atomic_valence_densities); // mix old
              // add contributions from smooth core densities, and optionally spherical valence densities
              stat += single_atom::atom_update("valence densities", na, 0, nr2.data(), ar2.data(), atom_rhoc.data());
              stat += add_smooth_quantities(rho_valence.data(), g, na, nr2.data(), ar2.data(), 
                                    center, n_periodic_images, periodic_images, atom_rhoc.data(),
                                    echo, 0, Y00sq*take_atomic_valence_densities, "smooth valence density");
          } // take_atomic_valence_densities

          if (echo > 4) print_stats(rho_valence.data(), g.all(), g.dV(), "# valence density");

          set(rho.data(), g.all(), rho_valence.data());

          if (echo > 4) print_stats(rho.data(), g.all(), g.dV(), "# density before adding smooth core densities:");

          stat += single_atom::atom_update("core densities", na, 0, nr2.data(), ar2.data(), atom_rhoc.data());
          stat += single_atom::atom_update("qlm charges", na, 0, 0, 0, atom_qlm.data());
          // add contributions from smooth core densities
          stat += add_smooth_quantities(rho.data(), g, na, nr2.data(), ar2.data(), 
                                center, n_periodic_images, periodic_images, atom_rhoc.data(),
                                echo, 0, Y00sq, "smooth core density");

          if (echo > 4) print_stats(rho.data(), g.all(), g.dV(), "# density  after adding smooth core densities:");

          here;

          { // scope: eval the XC potential and energy
              double E_xc{0}, E_dc{0};
              for(size_t i = 0; i < g.all(); ++i) {
                  auto const exc_i = exchange_correlation::LDA_kernel(rho[i], Vxc[i]);
                  E_xc += rho[i]*exc_i;
                  E_dc += rho[i]*Vxc[i]; // double counting correction 
                  // E_dc is computed just for display so we can compare E_dc between grid and atomic[SMT] contributions in calculation with a single atom
              } // i
              E_xc *= g.dV(); E_dc *= g.dV(); // scale with volume element
              if (echo > 2) std::printf("# exchange-correlation energy on grid %.9f %s, double counting %.9f %s\n", E_xc*eV,_eV, E_dc*eV,_eV);
              grid_xc_energy = E_xc;
          } // scope

          here;

          set(cmp.data(), g.all(), 0.0); // init compensation charge density, contains a smooth proton density and charge deficit compensators
          { // scope: solve the Poisson equation

              // add compensation charges cmp
              for(int ia = 0; ia < na; ++ia) {
                  double const sigma = sigma_cmp[ia];
                  int    const ellmax = lmax_qlm[ia];
                  if (echo > 6) std::printf("# use generalized Gaussians (sigma= %g %s, lmax=%d) as compensators for atom #%i\n", sigma*Ang, _Ang, ellmax, ia);
                  std::vector<double> coeff(sho_tools::nSHO(ellmax), 0.0);
                  stat += sho_projection::denormalize_electrostatics(coeff.data(), atom_qlm[ia], ellmax, sigma, unitary, echo);
//                if (echo > 7) std::printf("# before SHO-adding compensators for atom #%i coeff[000] = %g\n", ia, coeff[0]);
                  for(int ii = 0; ii < n_periodic_images; ++ii) {
                      double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                      stat += sho_projection::sho_add(cmp.data(), g, coeff.data(), ellmax, cnt, sigma, 0);
                  } // periodic images
#ifdef DEVEL
                  if (echo > 6) { // report extremal values of the density on the grid
                      std::printf("# after adding %g electrons compensator density for atom #%i:", atom_qlm[ia][00]*Y00inv, ia);
                      print_stats(cmp.data(), g.all(), g.dV());
                  } // echo
#endif // DEVEL
              } // ia

              // add compensators cmp to rho, so now rho == rho_aug
              add_product(rho.data(), g.all(), cmp.data(), 1.);
              if (echo > 1) print_stats(rho.data(), g.all(), g.dV(), "\n# augmented charge density:");

              { // scope: solve the Poisson equation: Laplace Ves == -4 pi rho
                  SimpleTimer timer(__FILE__, __LINE__, "Poisson equation", echo);
#ifdef DEVEL
                  if (echo > 0) {
                      std::printf("\n\n# %s\n# Solve Poisson equation\n# %s\n\n", h_line, h_line);
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
                      std::printf("# potential projection for atom #%d v_00 = %.9f %s\n", ia, atom_vlm[ia][00]*Y00*eV,_eV);
                      int const ellmax_show = std::min(ellmax, 2);
                      for(int ell = 1; ell <= ellmax_show; ++ell) {
                          double const unitfactor = Y00 * eV * std::pow(Ang, -ell);
                          int const ilm0 = sho_tools::lm_index(ell, -ell);
                          std::printf("# potential projection for atom #%d v_%im =", ia, ell);
                          printf_vector(" %.6f", &atom_vlm[ia][ilm0], 2*ell + 1, nullptr, unitfactor);
                          std::printf(" %s %s^%i\n", _eV, _Ang, -ell);
                      } // ell
                  } // echo
#endif // DEVEL
              } // ia

              
              double const E_es = 0.5*dot_product(g.all(), rho.data(), Ves.data())*g.dV();
              grid_electrostatic_energy = E_es; // store
              if (echo > 3) std::printf("# smooth electrostatic grid energy %.9f %s\n", E_es*eV,_eV);

          } // scope
          here;

          
          // communicate vlm to the atoms, get zero potential, atom-centered Hamiltonian and overlap
          float potential_mixing_ratio[] = {.5}; // {potential}
          stat += single_atom::atom_update("update", na, 0, 0, potential_mixing_ratio, atom_vlm.data());
          stat += single_atom::atom_update("hamiltonian", na, 0, 0, 0, atom_mat.data());
          stat += single_atom::atom_update("zero potentials", na, 0, nr2.data(), ar2.data(), atom_vbar.data());

          set(Vtot.data(), g.all(), Vxc.data()); add_product(Vtot.data(), g.all(), Ves.data(), 1.);

          if (echo > 1) print_stats(Vtot.data(), g.all(), 0, "\n# Total effective potential (before adding zero potentials)", eV);

          // now also add the zero potential vbar to Vtot
          stat += add_smooth_quantities(Vtot.data(), g, na, nr2.data(), ar2.data(), 
                                center, n_periodic_images, periodic_images, atom_vbar.data(),
                                echo, 0, Y00, "zero potential");

          if (echo > 1) print_stats(Vtot.data(), g.all(), 0, "\n# Total effective potential  (after adding zero potentials)", eV);
          here;

         
          /**  
           *  Potential generation done
           */

          
          double double_counting_correction{0};

          { // scope: solve the Kohn-Sham equation with the given Hamiltonian
              SimpleTimer KS_timer(__FILE__, __LINE__, "solving KS-equation", echo);
              
#ifdef DEVEL
              if (echo > 0) {
                  std::printf("\n\n# %s\n# Solve Kohn-Sham equation\n# %s\n\n", h_line, h_line);
                  std::fflush(stdout);
              } // echo
#endif // DEVEL     

#ifdef DEVEL
              if (echo > 8) {
                  for(int ia = 0; ia < na; ++ia) {
                      int const n = sho_tools::nSHO(numax[ia]);
                      view2D<double const> const aHm(atom_mat[ia], n);
                      std::printf("\n# atom-centered %d x %d Hamiltonian (in %s) for atom #%i\n", n, n, _eV, ia);
                      for(int i = 0; i < n; ++i) {
                          std::printf("#%3i  ", i);
                          printf_vector(" %.6f", aHm[i], n, "\n", eV);
                      } // i
                      view2D<double const> const aSm(atom_mat[ia] + n*n, n);
                      std::printf("\n# atom-centered %d x %d overlap matrix for atom #%i\n", n, n, ia);
                      for(int i = 0; i < n; ++i) {
                          std::printf("#%3i  ", i);
                          printf_vector(" %.6f", aSm[i], n);
                      } // i
                  } // ia
              } // echo
#endif // DEVEL

              view2D<double> rho_valence_new(2, g.all(), 0.0); // new valence density and response
              data_list<double> atom_rho_new[2];
              atom_rho_new[0] = data_list<double>(n_atom_rho, 0.0); // new valence density matrices
              atom_rho_new[1] = data_list<double>(n_atom_rho, 0.0); // and valence response matrices
              double charges[4] = {0, 0, 0, 0}; // 0:kpoint_denominator, 1:charge, 2:d_charge, 3:unused

              if (psi_on_grid) {
                
                  // restrict the local effective potential to the coarse grid
                  std::vector<double> Veff(gc.all());
                  multi_grid::restrict3D(Veff.data(), gc, Vtot.data(), g, 0); // mute
                  if (echo > 1) print_stats(Veff.data(), gc.all(), 0, "\n# Total effective potential  (restricted to coarse grid)   ", eV);
#ifdef DEVEL
                  if (0) { // scope: interpolate the effective potential to the dense grid again and compare it to the original version Vtot
                    // in order to test the interpolation routine
                      std::vector<double> v_dcd(g.all(), 0.0);
                      multi_grid::interpolate3D(v_dcd.data(), g, Veff.data(), gc, 0); // mute
                      if (echo > 1) print_stats(v_dcd.data(), g.all(), 0, "\n# Total effective potential (interpolated to dense grid)   ", eV);
                  } // scope

                  if (echo > 0) {
                      if (control::get("potential_generator.use.direct.projection", 0.) > 0) {
                          double cnt0[3]; if (na > 0) for(int d = 0; d < 3; ++d) cnt0[d] = center(0,d) + 0.5*(g.h[d] - gc.h[d]);
                          std::printf("\n## all values of Vtot in %s (on the coarse grid, unordered) as function of the distance to %s\n",
                                                            _eV, (na > 0) ? "atom #0" : "the cell center");
                          print_direct_projection(Veff.data(), gc, eV, (na > 0) ? cnt0 : nullptr);
                      } // control
                  } // echo
#endif // DEVEL

                  // copy the local potential and non-local atom matrices into the grid operator descriptor
                  op.set_potential(Veff.data(), gc.all(), atom_mat.data(), echo*0); // muted

                  auto const export_Hamiltonian = control::get("export.hamiltonian", 0.0);
                  if (export_Hamiltonian) {
                      op.write_to_file(echo, control::get("export.hamiltonian.format", "xml"));
                      if (export_Hamiltonian < 0) abort("Hamiltonian exported, export.hamiltonian = %g < 0", export_Hamiltonian);
                  } // export_Hamiltonian

                  view2D<double> rho_valence_gc(2, gc.all(), 0.0); // new valence density on the coarse grid and response density
                  
                  for(int ikpoint = 0; ikpoint < nkpoints; ++ikpoint) {
                      op.set_kpoint(kmesh[ikpoint], echo);
                      char x_axis[96]; std::snprintf(x_axis, 95, "# %g %g %g spectrum ", kmesh(ikpoint,0),kmesh(ikpoint,1),kmesh(ikpoint,2));
                      auto psi_k = psi[ikpoint]; // get a sub-view
                      bool display_spectrum{true};

                      // solve the Kohn-Sham equation using various solvers
                      if ('c' == *grid_eigensolver_method) { // "cg" or "conjugate_gradients"
                          stat += davidson_solver::rotate(psi_k.data(), energies[ikpoint], nbands, op, echo);
                          for(int irepeat = 0; irepeat < nrepeat; ++irepeat) {
                              if (echo > 6) { std::printf("# SCF cycle #%i, CG repetition #%i\n", scf_iteration, irepeat); std::fflush(stdout); }
                              stat += conjugate_gradients::eigensolve(psi_k.data(), energies[ikpoint], nbands, op, echo - 5);
                              stat += davidson_solver::rotate(psi_k.data(), energies[ikpoint], nbands, op, echo);
                          } // irepeat
                      } else
                      if ('d' == *grid_eigensolver_method) { // "davidson"
                          for(int irepeat = 0; irepeat < nrepeat; ++irepeat) {
                              if (echo > 6) { std::printf("# SCF cycle #%i, DAV repetition #%i\n", scf_iteration, irepeat); std::fflush(stdout); }
                              stat += davidson_solver::eigensolve(psi_k.data(), energies[ikpoint], nbands, op, echo);
                          } // irepeat
                      } else
                      if ('e' == *grid_eigensolver_method) { // "explicit" dense matrix solver
                          view3D<wave_function_t> HSm(2, gc.all(), align<4>(gc.all()), 0.0); // get memory for the dense representation
                          op.construct_dense_operator(HSm(0,0), HSm(1,0), HSm.stride(), echo);
                          stat += dense_solver::solve(HSm, x_axis, echo, nbands, energies[ikpoint]);
                          display_spectrum = false; // the dense solver will display on its own
                          wave_function_t const factor = 1./std::sqrt(gc.dV()); // normalization factor? 
                          for(int iband = 0; iband < nbands; ++iband) {
                              set(psi(ikpoint,iband), gc.all(), HSm(0,iband), factor);
                          } // iband
                      } else
                      if ('n' == *grid_eigensolver_method) { // "none"
                          if (take_atomic_valence_densities < 1) warn("no new valence density for grid.eigensolver=%s", grid_eigensolver_method);
                      } else {
                          ++stat; error("unknown grid.eigensolver=%s", grid_eigensolver_method);
                      } // grid_eigensolver_method

                      if (display_spectrum) dense_solver::display_spectrum(energies[ikpoint], nbands, x_axis, eV, _eV);

//                       if we used the fermi.level=linearized option, 
//                       we could combine the KS equation solving for a k-point with the evaluation of the density

                  } // ikpoint
                  op.set_kpoint(); // reset to Gamma
                  

                  here;

                  if ('e' == (occupation_method | 32)) {
                      // determine the exact Fermi level as a function of all energies and kmesh_weights
                      view2D<double> kweights(nkpoints, nbands), occupations(nkpoints, nbands);
                      for(int ikpoint = 0; ikpoint < nkpoints; ++ikpoint) {
                          double const kpoint_weight = kmesh(ikpoint,brillouin_zone::WEIGHT);
                          set(kweights[ikpoint], nbands, kpoint_weight);
                      } // ikpoint
                      double const eF = fermi_distribution::Fermi_level(occupations.data(), 
                                      energies.data(), kweights.data(), nkpoints*nbands,
                                      Fermi.get_temperature(), Fermi.get_n_electrons(), Fermi.get_spinfactor(), echo);
                      Fermi.set_Fermi_level(eF, echo);
                  } // occupation_method "exact"

                  here;

                  // density generation
                  for(int ikpoint = 0; ikpoint < nkpoints; ++ikpoint) {
                      op.set_kpoint(kmesh[ikpoint], echo);
                      double const kpoint_weight = kmesh(ikpoint,brillouin_zone::WEIGHT);
                      std::vector<uint32_t> coeff_starts;
                      auto const atom_coeff = density_generator::atom_coefficients(coeff_starts,
                                                psi(ikpoint,0), op, nbands, echo, ikpoint);
                      stat += density_generator::density(rho_valence_gc[0], atom_rho_new[0].data(), Fermi,
                                                energies[ikpoint], psi(ikpoint,0), atom_coeff.data(),
                                                coeff_starts.data(), na, gc, nbands, kpoint_weight, echo - 4, ikpoint, 
                                                         rho_valence_gc[1], atom_rho_new[1].data(), charges);
                  } // ikpoint
                  op.set_kpoint(); // reset to Gamma

                  here;

                  auto const dcc_coarse = dot_product(gc.all(), rho_valence_gc[0], Veff.data()) * gc.dV();
                  if (echo > 4) std::printf("\n# double counting (coarse grid) %.9f %s\n", dcc_coarse*eV, _eV);
                  // beware: if fermi.level=linearized, dcc_coarse is computed from uncorrected densities

                  stat += multi_grid::interpolate3D(rho_valence_new[0], g, rho_valence_gc[0], gc);
                  stat += multi_grid::interpolate3D(rho_valence_new[1], g, rho_valence_gc[1], gc);

              } else 
              if (plane_waves) {
                  here;
                  std::vector<plane_waves::DensityIngredients> export_rho;
                  here;

                  stat += plane_waves::solve(na, xyzZ, g, Vtot.data(), sigma_a.data(), numax.data(), atom_mat.data(), echo, &export_rho);

                  here;
                  
                  if ('e' == (occupation_method | 32)) {
                      // determine the Fermi level exactly as a function of all export_rho.energies and .kpoint_weight
                      view2D<double> kweights(nkpoints, nbands, 0.0), occupations(nkpoints, nbands);
                      for(int ikpoint = 0; ikpoint < export_rho.size(); ++ikpoint) {
                          set(kweights[ikpoint], nbands, export_rho[ikpoint].kpoint_weight);
                          set(energies[ikpoint], nbands, export_rho[ikpoint].energies.data());
                      } // ikpoint
                      double const eF = fermi_distribution::Fermi_level(occupations.data(), 
                                      energies.data(), kweights.data(), nkpoints*nbands,
                                      Fermi.get_temperature(), Fermi.get_n_electrons(), Fermi.get_spinfactor(), echo);
                      Fermi.set_Fermi_level(eF, echo);
                  } // occupation_method == "exact"

                  here;

                  for(auto & x : export_rho) {
                      if (echo > 1) { std::printf("\n# Generate valence density for %s\n", x.tag); std::fflush(stdout); }
                      stat += density_generator::density(rho_valence_new[0], atom_rho_new[0].data(), Fermi,
                                                x.energies.data(), x.psi_r.data(), x.coeff.data(),
                                                x.offset.data(), x.natoms, g, x.nbands, x.kpoint_weight, echo - 4, x.kpoint_index, 
                                                         rho_valence_new[1], atom_rho_new[1].data(), charges);
                  } // x

                  here;
              } else { // SHO local orbitals
                  here;

                  stat += sho_hamiltonian::solve(na, xyzZ, g, Vtot.data(), na, sigma_a.data(), numax.data(), atom_mat.data(), echo);

                  warn("with basis=%s no new density is generated", basis_method); // ToDo: implement this

                  here;
              } // psi_on_grid

              
              
              double band_energy_sum = Fermi.get_band_sum();

              Fermi.correct_Fermi_level(nullptr, echo); // clear the accumulators
              
              // correct if k-point weights do not sum up to 1
              if (charges[0] > 0 && std::abs(charges[0] - 1.0) > 1e-15) {
                  double const renormalization_factor = 1./charges[0];
                  if (echo > 2) std::printf("# %s: renormalize density and density matrices by %.15f\n", __func__, renormalization_factor);
                  scale(rho_valence_new[0], g.all(), renormalization_factor);
                  scale(rho_valence_new[1], g.all(), renormalization_factor);
                  for(int ia = 0; ia < na; ++ia) {
                      scale(atom_rho_new[0][ia], n_atom_rho[ia], renormalization_factor);
                      scale(atom_rho_new[1][ia], n_atom_rho[ia], renormalization_factor);
                  } // ia
                  scale(charges, 3, renormalization_factor);
              } // normalize by sum of k-weights
              
              if (echo > 2) std::printf("# %s: total charge %g electrons and derivative %g\n", __func__, 
                              charges[1], charges[2]*Fermi.get_temperature());

              // now correct for using the old Fermi level during density generation
              if (std::abs(charges[1] - Fermi.get_n_electrons()) > 1e-15) {
                  if (std::abs(charges[2]) > 0) {
                      double const alpha = (Fermi.get_n_electrons() - charges[1]) / charges[2];
                      if (echo > 1) std::printf("# shift Fermi level by %g %s\n", alpha*eV, _eV);
                      // correct by energy difference alpha
                      add_product(rho_valence_new[0], g.all(), rho_valence_new[1], alpha);
                      for(int ia = 0; ia < na; ++ia) {
                          add_product(atom_rho_new[0][ia], n_atom_rho[ia], atom_rho_new[1][ia], alpha);
                      } // ia
                      charges[1] += charges[2]*alpha;
                      band_energy_sum += Fermi.get_band_sum(1)*alpha;
                      Fermi.set_Fermi_level(Fermi.get_Fermi_level() + alpha, echo);
                  } else {
                      warn("in SCF iteration #%i number of valence electrons %g deviates from %g but response is zero",
                               scf_iteration, charges[1], Fermi.get_n_electrons());
                  } // response is non-zero
              } // the number of valence electrons does not match the requested number

              if (echo > 2) std::printf("# %s: total charge %g electrons\n", __func__, charges[1]);

              if (echo > 1) print_stats(rho_valence_new[0], g.all(), g.dV(), "\n# Total new valence density");

              // display the sum of new valence eigenvalues
              if (echo > 1) std::printf("\n# sum of eigenvalues %.9f %s\n\n", band_energy_sum*eV, _eV);
              // in order to compute the kinetic energy of valence states, we need to subtract 
              // the expectation value of the potential which consists of two parts: the grid part 
              double_counting_correction = dot_product(g.all(), rho_valence_new[0], Vtot.data()) * g.dV();
              if (echo > 1) std::printf("\n# grid double counting %.9f %s\n\n", double_counting_correction*eV, _eV);
              // and the atom part: for each atom sum_ij D_ij v_ij where v_ij + t_ij = h_ij, 
              // with h_ij being the atomic non-local correction to the Hamiltonian
              // and D_ij the atomic valence density matrix.

              grid_kinetic_energy = band_energy_sum - double_counting_correction;
              if (echo > 1) std::printf("\n# grid kinetic energy %.9f %s (take %.1f %%)\n\n", 
                  grid_kinetic_energy*eV, _eV, (1. - take_atomic_valence_densities)*100);
              
              // valence density mixing
              double const mix_old = 1 - density_mixing;
              scale(rho_valence.data(), g.all(), mix_old); // mix old
              add_product(rho_valence.data(), g.all(), rho_valence_new[0], density_mixing); // mix new
              // valence density matrix mixing
              for(int ia = 0; ia < na; ++ia) {
                  scale(atom_rho[ia], n_atom_rho[ia], mix_old); // mix old
                  add_product(atom_rho[ia], n_atom_rho[ia], atom_rho_new[0][ia], density_mixing); // mix new
              } // ia


          } // scope: Kohn-Sham

          float rho_mixing_ratios[] = {.5, .5, .5}; // for spherical {core, semicore, valence} density             
          stat += single_atom::atom_update("atomic density matrices", na, 0, 0, rho_mixing_ratios, atom_rho.data());
          

          // compute the total energy
          std::vector<double> atomic_energy_diff(na, 0.0);

          bool const total_energy_details = true;
          if (total_energy_details) {
              int const nE = align<1>(energy_contribution::max_number);
              std::vector<int32_t> nEa(na, nE);
              data_list<double> atom_contrib(nEa, 0.0);
              stat += single_atom::atom_update("energies", na, atomic_energy_diff.data(), 0, 0, atom_contrib.data());
              std::vector<double> Ea(nE, 0.0);
              for(int ia = 0; ia < na; ++ia) {
                  add_product(Ea.data(), nE, atom_contrib[ia], 1.0); // all atomic weight factors are 1.0
              } // ia
              for(int i01 = 0; i01 < 2; ++i01) {
                  if (echo > 3 + 4*(1 - i01)) {
                      char const *const without_with = i01 ? "" : "out";
                      std::printf("\n# sum of atomic energy contributions with%s grid contributions (%s)\n", without_with, _eV);
                      std::printf("# kinetic       %32.9f\n", Ea[energy_contribution::KINETIC]*eV);
                      std::printf("# electrostatic %32.9f\n", Ea[energy_contribution::ELECTROSTATIC]*eV);
                      std::printf("# XC            %32.9f\n", Ea[energy_contribution::EXCHANGE_CORRELATION]*eV);
                      std::printf("# true total    %32.9f\n", Ea[energy_contribution::TOTAL]*eV);
                      std::printf("# reference     %32.9f\n", Ea[energy_contribution::REFERENCE]*eV);
                      double E_tot = Ea[energy_contribution::KINETIC] + Ea[energy_contribution::ELECTROSTATIC]
                                   + Ea[energy_contribution::EXCHANGE_CORRELATION] - Ea[energy_contribution::REFERENCE];
                      std::printf("# total         %32.9f\n", E_tot*eV);
                  } // echo
                  
                  // now add grid contributions
                  Ea[energy_contribution::KINETIC] += grid_kinetic_energy * (1. - take_atomic_valence_densities);
                  Ea[energy_contribution::EXCHANGE_CORRELATION] += grid_xc_energy;
                  Ea[energy_contribution::ELECTROSTATIC] += grid_electrostatic_energy;
                  // reconstruct total energy from its contributions KINETIC + ES + XC
                  Ea[energy_contribution::TOTAL] = Ea[energy_contribution::KINETIC] 
                                                 + Ea[energy_contribution::ELECTROSTATIC]
                                                 + Ea[energy_contribution::EXCHANGE_CORRELATION];
                                                 
              } // show without and with grid contributions
          } else {
              stat += single_atom::atom_update("energies", na, atomic_energy_diff.data());
          } // total_energy_details

          double atomic_energy_corrections{0};
          for(int ia = 0; ia < na; ++ia) {
              atomic_energy_corrections += atomic_energy_diff[ia];
          } // ia
          total_energy = grid_kinetic_energy * (1. - take_atomic_valence_densities)
                       + grid_xc_energy
                       + grid_electrostatic_energy
                       + atomic_energy_corrections;
          if (echo > 0) { std::printf("\n# total energy %.9f %s\n\n", total_energy*eV, _eV); std::fflush(stdout); }

          
          
          
#if 0
          if (1) { // update take_atomic_valence_densities
              take_atomic_valence_densities = 0.f; // 0.99f;
              if (echo > 0) std::printf("# set take_atomic_valence_densities = %g %%\n", take_atomic_valence_densities*100);
              stat += single_atom::atom_update("lmax qlm", na, 0, lmax_qlm.data(), &take_atomic_valence_densities);
          } // 1
#endif // 0

          if (spherical_valence_decay > 0) { // update take_atomic_valence_densities
#if 1
              auto const x = scf_iteration/spherical_valence_decay; // progress x
              take_atomic_valence_densities = (x >= 1) ? 0 : 2*pow3(x) - 3*pow2(x) + 1; // smooth transition function
              if (echo > 0) std::printf("# set take_atomic_valence_densities = %g %%\n", take_atomic_valence_densities*100);
              stat += single_atom::atom_update("lmax qlm", na, 0, lmax_qlm.data(), &take_atomic_valence_densities);
#endif // 1
          } // spherical_valence_decay
          here;
          
          
          
          { // scope: read the stop file with standard name, max_scf_iterations may be modified
              auto const stat = debug_tools::manage_stop_file<'r'>(max_scf_iterations, echo);
              if (stat) warn("failed to read the stop file, status= %i", int(stat));
              if (max_scf_iterations_input != max_scf_iterations) {
                  warn("the max. number of SCF iterations has been modified from %d to %d "
                       "by the stop file during SCF iteration #%i", 
                       max_scf_iterations_input, max_scf_iterations, scf_iteration);
              } // stop file has been used
          } // scope
          density_mixing = control::get("potential_generator.mix.density", 0.25); // for the next iteration

      } // scf_iteration

      here;

      {   auto const stat = debug_tools::manage_stop_file<'d'>(max_scf_iterations, echo);
          if (stat) warn("failed to delete stop file, status= %i", int(stat));
      }

      if (psi_on_grid) {
          auto const store_wave_file = control::get("store.waves", "");
          if ((nullptr != store_wave_file) && ('\0' != store_wave_file[0])) {
              auto const nerrors = dump_to_file(store_wave_file,
                  nbands, psi(0,0), nullptr, psi.stride(), gc.all(), "wave functions", echo);
              if (nerrors) warn("%d errors occured writing file \'%s\'", nerrors, store_wave_file); 
          } // filename not empty
      } // wave functions on Cartesian real-space grid

      
#ifdef DEVEL
      stat += potential_projections(g, cell, Ves.data(), Vxc.data(), Vtot.data(), rho.data(), cmp.data(), na, &center, rcut, echo);
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
      int n{0}; int const t = control::get("potential_generator.select.test", -1.); // -1:all
      if (t & (1 << n++)) stat += test_init(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace potential_generator
