#include <cstdio> // printf, std::snprintf
#include <cassert> // assert
// #include <algorithm> // std::copy
// #include <cmath> // std::
#include <vector> // std::vector
#include <complex> // std::complex

#include "structure_solver.hxx"

#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set, pow2
#include "inline_tools.hxx" // align<n>
// #include "constants.hxx" // ::sqrtpi, ::pi
#include "real_space.hxx" // ::grid_t
#include "sho_tools.hxx" // ::quantum_number_table
#include "data_view.hxx" // view2D<T>
#include "data_list.hxx" // data_list<T> // replaces std::vector<T*> constructions

// #include "simple_timer.hxx" // // SimpleTimer
#include "control.hxx" // ::get

#include "boundary_condition.hxx" // Periodic_Boundary, Isolated_Boundary
#include "debug_tools.hxx" // ::read_from_file

#include "atom_image.hxx"// ::sho_atom_t
#include "grid_operators.hxx" // ::grid_operator_t
#include "conjugate_gradients.hxx" // ::eigensolve
#include "davidson_solver.hxx" // ::rotate, ::eigensolve
#include "multi_grid.hxx" // ::restrict3D, ::interpolate3D
#include "density_generator.hxx" // ::density
#include "sho_hamiltonian.hxx" // ::solve
#include "plane_waves.hxx" // ::solve
#include "fermi_distribution.hxx" // ::FermiLevel_t, ::Fermi_level
// #include "unit_system.hxx" // ::length_unit

#include "brillouin_zone.hxx" // ::get_kpoint_mesh, ::needs_complex

#include "print_tools.hxx" // print_stats, printf_vector
#include "dense_solver.hxx" // ::display_spectrum, ::solve
#include "complex_tools.hxx" // complex_name

#define DEBUG
#ifdef  DEBUG
    #include "debug_output.hxx" // here
#else  // DEBUG
    #define here
#endif // DEBUG

namespace structure_solver {
  /** 
      Electronic structure solver 
      generates a new density given a DFT potential
   */

  template <typename wave_function_t>
  class KohnShamStates {
  public:
    using real_wave_function_t = decltype(std::real(wave_function_t(1)));
    
    KohnShamStates() {} // default constructor
    KohnShamStates(   // preferred constructor
          real_space::grid_t const & coarse_grid
        , std::vector<atom_image::sho_atom_t> & list_of_atoms
        , view2D<double> const & kpoint_mesh
        , int const number_of_kpoints=1
        , int const number_of_bands=1
        , int const run=1 // 1:run, 0:check
        , int const echo=0 // log-level
    )
      : gc(coarse_grid)
      , op(gc, list_of_atoms)
      , kmesh(kpoint_mesh)
      , grid_eigensolver_method(control::get("grid.eigensolver", "cg"))
      , nkpoints(number_of_kpoints)
      , nbands(number_of_bands)
      , nrepeat(1)
    {
        if (echo > 0) printf("# %s use  %d x %d x %d  coarse grid points\n", __func__, gc[0], gc[1], gc[2]);
        nrepeat = int(control::get("grid.eigensolver.repeat", 1.)); // repetitions of the solver

        // construct grid-based Hamiltonian and overlap operator descriptor
//         op = grid_operators::grid_operator_t<wave_function_t, real_wave_function_t>(gc, list_of_atoms);
        // Mind that local potential and atom matrices of op are still unset!

        if (echo > 0) printf("# real-space grid wave functions are of type %s\n", complex_name<wave_function_t>());
        psi = view3D<wave_function_t>(run*nkpoints, nbands, gc.all(), 0.0); // get memory (potentially large)

                  int const na = list_of_atoms.size();
                  auto const start_wave_file = control::get("start.waves", "");
                  if ('\0' == *start_wave_file) {
                      if (echo > 1) printf("# initialize grid wave functions as %d atomic orbitals on %d atoms\n", nbands, na);
                      float const scale_sigmas = control::get("start.waves.scale.sigma", 10.); // how much more spread in the start waves compared to sigma_prj
                      uint8_t qn[20][4]; // first 20 sets of quantum numbers [nx, ny, nz, nu] with nu==nx+ny+nz
                      sho_tools::quantum_number_table(qn[0], 3, sho_tools::order_Ezyx); // Ezyx-ordered, take 1, 4, 10 or 20
                      std::vector<int32_t> ncoeff_a(na, 20);
                      data_list<wave_function_t> single_atomic_orbital(ncoeff_a, 0.0); // get memory and initialize
                      for(int ikpoint = 0; ikpoint < nkpoints; ++ikpoint) {
                          op.set_kpoint(kmesh[ikpoint], echo);
                          for(int iband = 0; iband < nbands; ++iband) {
                              int const ia = iband % na; // which atom?
                              int const io = iband / na; // which orbital?
                              if (io >= 20) error("requested more than 20 start wave functions per atom! bands.per.atom=%g", nbands/double(na));
                              auto const q = qn[io];
                              if (echo > 7) printf("# initialize band #%i as atomic orbital %x%x%x of atom #%i\n", iband, q[2], q[1], q[0], ia);
                              int const isho = sho_tools::zyx_index(3, q[0], q[1], q[2]); // isho in order_zyx w.r.t. numax=3
                              single_atomic_orbital[ia][isho] = 1./std::sqrt((q[3] > 0) ? ( (q[3] > 1) ? 53. : 26.5 ) : 106.); // set normalization depending on s,p,ds*
                              if (run) op.get_start_waves(psi(ikpoint,iband), single_atomic_orbital.data(), scale_sigmas, echo);
                              single_atomic_orbital[ia][isho] = 0; // reset
                          } // iband
                      } // ikpoint
                      op.set_kpoint(); // reset k-point
                  } else {
                      if (echo > 1) printf("# try to read start waves from file \'%s\'\n", start_wave_file);
                      if (run) {
                          auto const nerrors = debug_tools::read_from_file(psi(0,0), start_wave_file, nbands, psi.stride(), gc.all(), "wave functions", echo);
                          if (nerrors) {
                              error("failed to read start wave functions from file \'%s\'", start_wave_file);
                          } else {
                              if (echo > 1) printf("# read %d bands x %ld numbers from file \'%s\'\n", nbands, gc.all(), start_wave_file);
                          }
                      } // run
                      for(int ikpoint = 1; ikpoint < nkpoints; ++ikpoint) {
                          if (echo > 3) { printf("# copy %d bands for k-point #%i from k-point #0\n", nbands, ikpoint); std::fflush(stdout); }
                          if (run) set(psi(ikpoint,0), psi.dim1()*psi.stride(), psi(0,0)); // copy, ToDo: include Bloch phase factors
                      } // ikpoints
                  } // start wave method

    } // constructor
    
    status_t solve(
          view2D<double> & rho_valence_gc
        , data_list<double> atom_rho_new[2]
        , view2D<double> & energies
        , double charges[]
        , fermi_distribution::FermiLevel_t & Fermi
        , view2D<double> const & Veff
        , data_list<double> const & atom_mat
        , char const occupation_method='e'
        , char const *grid_eigensolver_method="cg"
        , int const scf_iteration=-1
        , int const echo=0
    ) {
        status_t stat(0);
        int const na = atom_mat.nrows();

                  // copy the local potential and non-local atom matrices into the grid operator descriptor
                  op.set_potential(Veff.data(), gc.all(), atom_mat.data(), echo*0); // muted
                  
                  for(int ikpoint = 0; ikpoint < nkpoints; ++ikpoint) {
                      op.set_kpoint(kmesh[ikpoint], echo);
                      char x_axis[96]; std::snprintf(x_axis, 95, "# %g %g %g spectrum ", kmesh(ikpoint,0),kmesh(ikpoint,1),kmesh(ikpoint,2));
                      auto psi_k = psi[ikpoint]; // get a sub-view
                      bool display_spectrum{true};

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
//                           if (take_atomic_valence_densities < 1) warn("eigensolver=none generates no new valence density");
//                           ToDo: this warning in potential_generator.cxx
                      } else {
                          ++stat; error("unknown grid.eigensolver method \'%s\'", grid_eigensolver_method);
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
      
        return stat;
    } // solve

    ~KohnShamStates() {
        std::printf("# ~KohnShamStates<%s>\n", complex_name<wave_function_t>());
        psi = view3D<wave_function_t>(0,0,0, 0);
    } // destructor

  private: // members

    real_space::grid_t const & gc;
    grid_operators::grid_operator_t<wave_function_t, real_wave_function_t> op;
    view3D<wave_function_t> psi; // Kohn-Sham states in real-space grid representation
    view2D<double> const & kmesh; // kmesh(nkpoints, 4);
    char const *grid_eigensolver_method;
    int nkpoints;
    int nbands;
    int nrepeat;

  }; // class KohnShamStates<>

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  class RealSpaceKohnSham {
  public:
    
    RealSpaceKohnSham(
          real_space::grid_t const & g
        , std::vector<atom_image::sho_atom_t> & loa // list of atoms
        , int const run=1 // 1:run, 0:check
        , int const echo=0 // log-level
    )
      : gd(g)
    {
        // how to solve the KS-equation
        basis_method = control::get("basis", "grid");
        plane_waves = ((*basis_method | 32) == 'p');
        psi_on_grid = ((*basis_method | 32) == 'g');
        
        // get a default kmesh controlled by +hamiltonian.kmesh or +hamiltonian.kmesh.x, .y, .z
        nkpoints = brillouin_zone::get_kpoint_mesh<true>(kmesh);
        if (echo > 1) printf("# k-point mesh has %d points\n", nkpoints);
        // ToDo: warn if boundary_condition is isolated but there is more than 1 kpoint
        bool const needs_complex = brillouin_zone::needs_complex(kmesh, nkpoints);
        if (echo > 3) printf("# k-point %s wavefunctions\n", needs_complex?"needs complex":"allows real");

        occupation_method = *control::get("fermi.level", "exact"); // {"exact", "linearized"}

        int const na = loa.size(); // get the number of atoms from the list of atoms
        double const nbands_per_atom = control::get("bands.per.atom", 10.); // 1: s  4: s,p  10: s,p,ds*  20: s,p,ds*,fp*
        nbands = std::max(1, int(nbands_per_atom*na));

        if (psi_on_grid) {

            // ============================================================================================
            // prepare for solving the Kohn-Sham equation on the real-space grid
            // create a coarse grid descriptor (this could be done later but we want a better overview in the log file)
            gc = real_space::grid_t(g[0]/2, g[1]/2, g[2]/2); // divide the dense grid numbers by two
            gc.set_grid_spacing(2*g.h[0], 2*g.h[1], 2*g.h[2]); // alternative: cell[]/gc[]
            gc.set_boundary_conditions(g.boundary_conditions());

            if (echo > 1) {
                printf("# use  %d x %d x %d  coarse grid points\n", gc[0], gc[1], gc[2]);
                double const max_grid_spacing = std::max(std::max(gc.h[0], gc.h[1]), gc.h[2]);
                printf("# use  %g %g %g  %s  coarse grid spacing, corresponds to %.2f Ry\n",
                          gc.h[0]*Ang, gc.h[1]*Ang, gc.h[2]*Ang, _Ang, pow2(constants::pi/max_grid_spacing));
            } // echo

            int const floating_point_bits = control::get("hamiltonian.floating.point.bits", 64.); // double by default
            if (32 == floating_point_bits) {
                key = needs_complex ? 'c' : 's';
            } else if (64 == floating_point_bits) {
                key = needs_complex ? 'z' : 'd';
            } else {
                error("hamiltonian.floating.point.bits=%d must be 32 or 64 (default)", floating_point_bits);
            }

            if ('z' == key) z = new KohnShamStates<std::complex<double>>(gc, loa, kmesh, nkpoints, nbands, run, echo);
            if ('c' == key) c = new KohnShamStates<std::complex<float>> (gc, loa, kmesh, nkpoints, nbands, run, echo);
            if ('d' == key) d = new KohnShamStates<double>              (gc, loa, kmesh, nkpoints, nbands, run, echo);
            if ('s' == key) s = new KohnShamStates<float>               (gc, loa, kmesh, nkpoints, nbands, run, echo);

        } else {
            key = *basis_method;
        } // psi_on_grid

        energies = view2D<double>(nkpoints, nbands, 0.0); // Kohn-Sham eigenenergies

    } // constructor
    
    ~RealSpaceKohnSham() { // destructor
       if (psi_on_grid) {
            // call destructors
            if ('z' == key && nullptr != z) z->~KohnShamStates();
            if ('c' == key && nullptr != c) c->~KohnShamStates();
            if ('d' == key && nullptr != d) d->~KohnShamStates();
            if ('s' == key && nullptr != s) s->~KohnShamStates();
        } // psi_on_grid
    } // destructor

    status_t solve(
          fermi_distribution::FermiLevel_t & Fermi
        , real_space::grid_t const & g // dense grid
        , double const Vtot[] // local effective potential
        , std::vector<int32_t> const & n_atom_rho
        , data_list<double> const & atom_mat
        , char const occ='e' // occupation method
        , char const *solver="cg" // grid_eigensolver_method
        , int const scf=-1 // scf_iteration
        , int const echo=9 // log-level
    ) {
        status_t stat(0);
#ifdef DEVEL
              if (echo > 0) {
                  printf("\n\n#\n# Solve Kohn-Sham equation\n# \n\n");
                  std::fflush(stdout);
              } // echo
#endif // DEVEL

//           double double_counting_correction{0}; // ToDo

              view2D<double> rho_valence_new(2, g.all(), 0.0); // new valence density and response
              data_list<double> atom_rho_new[2];
              atom_rho_new[0] = data_list<double>(n_atom_rho, 0.0); // new valence density matrices
              atom_rho_new[1] = data_list<double>(n_atom_rho, 0.0); // and valence response matrices
              double charges[4] = {0, 0, 0, 0}; // 0:kpoint_denominator, 1:charge, 2:d_charge, 3:unused
      
        if (psi_on_grid) {

                  // restrict the local effective potential to the coarse grid
                  view2D<double> Veff(1, gc.all());
                  multi_grid::restrict3D(Veff[0], gc, Vtot, gd, 0); // mute
                  if (echo > 1) print_stats(Veff[0], gc.all(), 0, "\n# Total effective potential  (restricted to coarse grid)   ", eV);

                  view2D<double> rho_valence_gc(2, gc.all(), 0.0); // new valence density on the coarse grid and response density
                  
                  if ('z' == key) z->solve(rho_valence_gc, atom_rho_new, energies, charges, Fermi, Veff, atom_mat, occ, solver, scf, echo);
                  if ('c' == key) c->solve(rho_valence_gc, atom_rho_new, energies, charges, Fermi, Veff, atom_mat, occ, solver, scf, echo);
                  if ('d' == key) d->solve(rho_valence_gc, atom_rho_new, energies, charges, Fermi, Veff, atom_mat, occ, solver, scf, echo);
                  if ('s' == key) s->solve(rho_valence_gc, atom_rho_new, energies, charges, Fermi, Veff, atom_mat, occ, solver, scf, echo);

                  auto const dcc_coarse = dot_product(gc.all(), rho_valence_gc[0], Veff.data()) * gc.dV();
                  if (echo > 4) printf("\n# double counting (coarse grid) %.9f %s\n", dcc_coarse*eV, _eV);
                  // beware: if fermi.level=linearized, dcc_coarse is computed from uncorrected densities

                  stat += multi_grid::interpolate3D(rho_valence_new[0], g, rho_valence_gc[0], gc);
                  stat += multi_grid::interpolate3D(rho_valence_new[1], g, rho_valence_gc[1], gc);

        } // psi_on_grid
#if 0
        else {
            int const na = n_atom_rho.size();
            if (plane_waves) {

                  std::vector<plane_waves::DensityIngredients> export_rho;

                  stat += plane_waves::solve(na, xyzZ, g, Vtot, sigma_a.data(), numax.data(), atom_mat.data(), echo, &export_rho);

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
                  } // occupation_method "exact"

                  for(auto & x : export_rho) {
                      if (echo > 1) { printf("\n# Generate valence density for %s\n", x.tag); std::fflush(stdout); }
                      stat += density_generator::density(rho_valence_new[0], atom_rho_new[0].data(), Fermi,
                                                x.energies.data(), x.psi_r.data(), x.coeff.data(),
                                                x.offset.data(), x.natoms, g, x.nbands, x.kpoint_weight, echo - 4, x.kpoint_index, 
                                                         rho_valence_new[1], atom_rho_new[1].data(), charges);
                  } // ikpoint
            
            } else {

                  stat += sho_hamiltonian::solve(na, xyzZ, g, Vtot, na, sigma_a.data(), numax.data(), atom_mat.data(), echo);
                  warn("with basis=%s no new density is generated", basis_method); // ToDo: implement this

            }
        }
#endif
        here;
        return stat;
    } // solve
    
  private:
      char key;
      view2D<double> kmesh; // kmesh(nkpoints, 4);
      int nkpoints;
      int nbands;
      real_space::grid_t const & gd; // dense grid
      real_space::grid_t gc; // coarse grid descriptor

      bool psi_on_grid, plane_waves;

      // 2x2 versions for real space Kohn-Sham wave functions
      KohnShamStates<std::complex<double>> *z = nullptr;
      KohnShamStates<std::complex<float>>  *c = nullptr;
      KohnShamStates<double>               *d = nullptr;
      KohnShamStates<float>                *s = nullptr;

      char occupation_method;
      char const *basis_method;
  public:
      view2D<double> energies; // energies(nkpoints, nbands, 0.0); // Kohn-Sham eigenenergies

  }; // class RealSpaceKohnSham

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo) {
      status_t stat(0);
      real_space::grid_t g(7, 8, 9);
      std::vector<atom_image::sho_atom_t> list_of_atoms(0);
      RealSpaceKohnSham KS(g, list_of_atoms, 1, echo);
      return stat;
  } // test_create_and_destroy

  status_t test_free_electrons(int const echo) {
      status_t stat(0);
      real_space::grid_t g(8, 8, 8);
      std::vector<atom_image::sho_atom_t> list_of_atoms(0);
      RealSpaceKohnSham KS(g, list_of_atoms, 1, echo);
      std::vector<int32_t> n_atom_rho(0);
      data_list<double> atom_mat(n_atom_rho);
      fermi_distribution::FermiLevel_t Fermi;
      std::vector<double> Vtot(g.all(), 0.0);
      KS.solve(Fermi, g, Vtot.data(), n_atom_rho, atom_mat);
      return stat;
  } // test_free_electrons

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_create_and_destroy(echo);
      stat += test_free_electrons(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace structure_solver
