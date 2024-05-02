// This file is part of AngstromCube under MIT License

#include <cstdint> // int8_t
#include <cassert> // assert
#include <cstdio> // std::printf, ::snprintf
#include <algorithm> // std::copy, ::min, ::max
#include <vector> // std::vector

#include "self_consistency.hxx"

#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set, pow2, align<nBits>
#include "constants.hxx" // ::sqrtpi, ::pi
#include "solid_harmonics.hxx" // ::Y00, ::cleanup
#include "real_space.hxx" // ::grid_t, ::add_function
#include "chemical_symbol.hxx" // ::get
#include "sho_projection.hxx" // ::sho_add, ::sho_project
#include "sho_tools.hxx" // ::nSHO, ::lm_index
#include "exchange_correlation.hxx" // ::LDA_kernel
#include "boundary_condition.hxx" // ::periodic_images
#include "data_view.hxx" // view2D<T>
#include "data_list.hxx" // data_list<T>

#include "geometry_input.hxx" // ::read_xyz_file, ::get_temperature, ::init_geometry_and_grid
#include "geometry_analysis.hxx" // ::fold_back, length
#include "simple_timer.hxx" // SimpleTimer
#include "control.hxx" // ::get, ::set, ::echo_set_without_warning

#include "print_tools.hxx" // print_stats, printf_vector

#include "sho_unitary.hxx" // ::Unitary_SHO_Transform

#include "single_atom.hxx" // ::atom_update
#include "energy_contribution.hxx" // ::show, ::TOTAL, ::KINETIC, ::ELECTROSTATIC, ...

#include "structure_solver.hxx" // ::RealSpaceKohnSham
#include "potential_generator.hxx" // ::add_smooth_quantities
#ifdef    DEVEL
  #include "potential_generator.hxx" // ::potential_projections
#endif // DEVEL

#include "fermi_distribution.hxx" // ::FermiLevel_t
#include "unit_system.hxx" // ::length_unit

#include "poisson_solver.hxx" // ::solve, ::solver_method


#define   DEBUG
#ifdef    DEBUG
    #include "debug_output.hxx" // here
#else  // DEBUG
    #define here
#endif // DEBUG

namespace self_consistency {
  // This module makes an all-electron Density Functional Theory (DFT) calculation of atoms
  // described by the Projector Augmented Wave method (PAW) 

  status_t manage_stop_file(int & number, char const write_read_delete='w', int const echo=0, char const *filename="running.a43") {
      if ('w' == (write_read_delete | 32)) {
          auto *const f = std::fopen(filename, "w");
          if (nullptr == f) {
              if (echo > 1) std::printf("# %s<%c> unable to open file \"%s\"for writing!\n", __func__, write_read_delete, filename);
              return 1;
          } // failed to open
          std::fprintf(f, "%d   is the max. number of self-consistency iterations, user may modify", number);
          return std::fclose(f);
      } else
      if ('r' == (write_read_delete | 32)) {
          std::ifstream infile(filename, std::ifstream::in);
          if (infile.fail()) {
              if (echo > 0) std::printf("# %s<%c> unable to find file \"%s\" for reading\n", __func__, write_read_delete, filename);
              return 1; // error
          }
          infile >> number;
          if (echo > 1) std::printf("# %s<%c> found number=%d in file \"%s\" \n", __func__, write_read_delete, number, filename);
          return 0;
      } else
      if ('d' == (write_read_delete | 32)) {
          if (echo > 0) std::printf("# %s<%c> removes file \"%s\"\n", __func__, write_read_delete, filename);
          return std::remove(filename);
      } else {
          if (echo > 0) std::printf("# %s<%c> only 'w'/'W'/write, 'r'/'R'/read or 'd'/'D'/delete implemented!\n", __func__, write_read_delete);
          return -1; // error
      }
  } // manage_stop_file


  status_t SCF(int const echo) { // =0 // log-level
      // compute the self-consistent solution of a single_atom, all states in the core
      // get the spherical core_density and bring it to the 3D grid
      // get the ell=0 compensator charge and add it to the 3D grid
      // envoke exchange_correlation and fourier_poisson
      // add XC and electrostatic potential and zero potential contributions
      // feed back potential shifts into single_atom
      status_t stat(0);

      int const check = control::get("check", 0.);
      int const run = std::min(std::max(0, 1 - check), 1);

      double constexpr Y00 = solid_harmonics::Y00;
      double constexpr Y00sq = pow2(Y00);
#ifdef    DEVEL
      double constexpr Y00inv = solid_harmonics::Y00inv;
#endif // DEVEL

      char h_line[32]; set(h_line, 31, '='); h_line[31] = '\0'; // a horizontal line of 31x '='
      if (echo > 0) std::printf("\n\n# %s\n# Initialize\n# +check=%i  (run=%i)\n# %s\n\n",
                            h_line, check, run, h_line);

      view2D<double> xyzZ;
      real_space::grid_t g;
      int na_noconst{0};
      stat += geometry_input::init_geometry_and_grid(g, xyzZ, na_noconst, 2, echo);
      int const na{na_noconst}; // total number of atoms

      std::vector<float> ionization(na, 0.f);
      { // scope: ionization between first and last atom
          auto const ion = control::get("self_consistency.test.ion", 0.0);
#ifdef    DEVEL
          if (0 != ion && na > 0) {
              if (echo > 2) std::printf("# %s distribute ionization of %g electrons between first and last atom\n", __func__, ion);
              ionization[na - 1] = -ionization[0];
              ionization[0] = ion;
              if (echo > 2) {
                  std::printf("# %s ionizations:", __func__);
                  printf_vector(" %g", ionization.data(), na);
              } // echo
          } // na
#else  // DEVEL
        if (0 != ion) error("need to activate ionizations with -D DEVEL, requested ion= %g electrons", ion);
#endif // DEVEL
      } // ion

      float const rcut = 32; // radial grids usually end at 9.45 Bohr
      view2D<double> periodic_images;
      int const n_periodic_images = boundary_condition::periodic_images(periodic_images,
                                       g.cell, g.boundary_conditions(), rcut, echo - 4);
      if (echo > 1) std::printf("# %s consider %d periodic images\n", __func__, n_periodic_images);


      std::vector<double> Za(na); // list of atomic numbers
      view2D<double> center(na, 4, 0.0); // get memory for a list of atomic centers
      { // scope: prepare atomic coordinates
          double const grid_offset[3] = {0.5*(g[0] - 1)*g.h[0],
                                         0.5*(g[1] - 1)*g.h[1],
                                         0.5*(g[2] - 1)*g.h[2]};
          if (echo > 4) std::printf("\n# %s List of Atoms: (coordinates in %s)\n", __func__, _Ang);
          for (int ia = 0; ia < na; ++ia) {
              double const Z = xyzZ(ia,3);
              char Symbol[4]; chemical_symbol::get(Symbol, Z, ' ');
              if (echo > 4) std::printf("# %s  %15.9f %15.9f %15.9f", Symbol,
                              xyzZ(ia,0)*Ang, xyzZ(ia,1)*Ang, xyzZ(ia,2)*Ang);
              Za[ia] = Z;
              for (int d = 0; d < 3; ++d) {
                  center(ia,d) = geometry_analysis::fold_back(xyzZ(ia,d), g.cell[d][d]) + grid_offset[d]; // w.r.t. to the center of grid point (0,0,0)
              } // d
              center(ia,3) = 0; // 4th component is not used
              if (echo > 4) std::printf("\n");
          } // ia
          if (echo > 4) std::printf("\n");
      } // scope


      std::vector<double> rho_valence(run*g.all(), 0.0);

      char const *initial_valence_density_method = control::get("initial.valence.density", "atomic"); // {"atomic", "load", "none"}
      if (echo > 0) std::printf("\n# initial.valence.density=%s\n", initial_valence_density_method);

      auto const spherical_valence_decay = control::get("atomic.valence.decay", 10.); // after SCF iteration # 10, take_atomic_valence_densities is zero, never if 0
      float take_atomic_valence_densities{0};

      if ('a' == *initial_valence_density_method) { // "atomic"
          if (echo > 1) std::printf("# include spherical atomic valence densities in the smooth core densities\n");
          take_atomic_valence_densities = 1; // 100% of the smooth spherical atomic valence densities is included in the smooth core densities
      } else if ('l' == *initial_valence_density_method) { // "load"
          error("initial.valence.density=load has not been implemented yet, found %s", initial_valence_density_method);
      } else if ('n' == *initial_valence_density_method) { // "none"
          warn("initial.valence.density=none may cause problems, found %s", initial_valence_density_method);
      } else {
          warn("initial.valence.density=%s is unknown, valence density is zero", initial_valence_density_method);
      } // initial_valence_density_method


      std::vector<double> sigma_cmp(na, 1.); // spread of the Gaussian used in the compensation charges
      char const pawdata_from = (*control::get("pawdata.from", "auto")) | 32; // 'a': auto, 'f': pawxml_import
      std::vector<int32_t> numax(na, ('f' == pawdata_from)?-9:-1); // -1: LivePAW, -9: load from pawxml files
      std::vector<int32_t> lmax_qlm(na, -1);
      std::vector<int32_t> lmax_vlm(na, -1);

      // initialize and get sigma, lmax for each atom
      stat += single_atom::atom_update("initialize", na, Za.data(), numax.data(), ionization.data(), (double**)1);
      stat += single_atom::atom_update("lmax qlm",   na,    nullptr, lmax_qlm.data(), &take_atomic_valence_densities);
      stat += single_atom::atom_update("lmax vlm",   na, (double*)1, lmax_vlm.data());
      stat += single_atom::atom_update("sigma cmp",  na, sigma_cmp.data());

      double nve{0}; // non-const total number of valence electrons
      { // scope: determine the number of valence electrons
          auto const keyword_valence_electrons = "valence.electrons";
          if ('a' == (*control::get(keyword_valence_electrons, "auto") | 32)) {
              std::vector<double> n_electrons_a(na, 0.); // number of valence electrons added by each atom
              stat += single_atom::atom_update("#valence electrons", na, n_electrons_a.data());
              nve = std::accumulate(n_electrons_a.begin(), n_electrons_a.end(), 0.0);
              if (echo > 0) std::printf("\n# %s=auto --> %g valence electrons\n\n", keyword_valence_electrons, nve);
              control::set(keyword_valence_electrons, nve, control::echo_set_without_warning);
          } else {
              nve = control::get(keyword_valence_electrons, 0.0);
              if (echo > 0) std::printf("\n# %s=%g\n\n", keyword_valence_electrons, nve);
          } // auto
      } // scope
      double const n_valence_electrons = nve; // do not use nve beyond this point

      std::vector<int32_t> n_atom_rho(na, 0);
      data_list<double> atom_qlm, atom_vlm, atom_rho, atom_mat;
      if (1) { // scope: set up data_list items
          view2D<int32_t> num(4, na, 0); // how many
          for (int ia = 0; ia < na; ++ia) {
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

      sho_unitary::Unitary_SHO_Transform const unitary(9);

      // allocations of grid quantities
      std::vector<double>  rho(run*g.all()); // [augmented] charge density
      std::vector<double>  Vxc(run*g.all()); // exchange-correlation potential
      std::vector<double>  cmp(run*g.all()); // compensation charge densities
      std::vector<double>  Ves(run*g.all(), 0.0); // electrostatic potential
      std::vector<double> Vtot(run*g.all()); // total effective potential

      // configuration
      auto const *es_solver_name = control::get("electrostatic.solver", "fft"); // {"fft", "multi-grid", "MG", "CG", "SD", "none"}
      if (echo > 2) std::printf("# electrostatic.solver=%s from {fft, multi-grid, MultiGrid, CG, SD, none}\n", es_solver_name);
      auto const es_solver_method = poisson_solver::solver_method(es_solver_name);

      auto const compensator_method = *control::get("electrostatic.compensator", "factorizable") | 32; // {'f', 'g'}
      if (echo > 2) std::printf("# electrostatic.compensator=%c from {factorizable, generalizedGaussian}\n", compensator_method);

      auto const occupation_method = *control::get("fermi.level", "exact") | 32; // {'e':"exact", 'l':"linearized"}
      if (echo > 2) std::printf("# fermi.level=%c from {exact, linearized}\n", occupation_method);

      // create a FermiLevel object
      fermi_distribution::FermiLevel_t Fermi(n_valence_electrons, 2, geometry_input::get_temperature(echo), echo);

      double const density_mixing_fixed = control::get("self_consistency.mix.density", 0.25);
      double density_mixing{1.0}; // initialize with 100% since we have no previous density

      std::vector<double> sigma_a(na, .5);
      { // scope: collect information for projectors and construct a list of atoms
          std::vector<int32_t> numax_a(na, 3);
          stat += single_atom::atom_update("projectors", na, sigma_a.data(), numax_a.data());
          for (int ia = 0; ia < na; ++ia) {
              assert(numax_a[ia] == numax[ia] && "consistency between atom_update('i') and ('p')");
          } // ia
      } // scope

      view2D<double> xyzZinso(na, 8);
      for (int ia = 0; ia < na; ++ia) {
          set(xyzZinso[ia], 4, xyzZ[ia]); // copy x,y,z,Z
          xyzZinso(ia,4) = ia;            // global_atom_id
          xyzZinso(ia,5) = numax[ia];
          xyzZinso(ia,6) = sigma_a[ia];
          xyzZinso(ia,7) = 0;             // __not_used__
      } // ia
      // ============================================================================
      // prepare for solving the Kohn-Sham equation
      // ============================================================================

      structure_solver::RealSpaceKohnSham KS(g, xyzZinso, na, run, echo);

      if (check) {
          // do not start SCF iterations!
          if (echo > 0) std::printf("# +check=%i --> done\n", check);
          stat += single_atom::atom_update("memory cleanup", na);
          return stat;
      } // check only, no run

      // total energy contributions
      double grid_electrostatic_energy{0}, grid_kinetic_energy{0}, grid_xc_energy{0}, total_energy{0};

      int const max_scf_iterations_input = control::get("self_consistency.max.scf", 1.);
      int max_scf_iterations{max_scf_iterations_input}; // non-const, may be modified by the stop file during the run
      { // scope: create a stop file with standard name
          auto const stat = manage_stop_file(max_scf_iterations, 'w', echo);
          if (stat != 0) error("failed to write/create a stop file, status= %i", int(stat));
      } // scope

      for (int scf_iteration = 0; scf_iteration < max_scf_iterations; ++scf_iteration) {
          SimpleTimer scf_iteration_timer(__FILE__, __LINE__, "scf_iteration", echo);
          if (echo > 1) std::printf("\n\n# %s\n# SCF-iteration step #%i:\n# %s\n\n", h_line, scf_iteration, h_line);

          if (take_atomic_valence_densities > 0) {
              if (echo > 4) print_stats(rho_valence.data(), g.all(), g.dV(), "# previous valence density");

              if (echo > 4) std::printf("# compose valence density with %g %% of the atomic valence densities\n", take_atomic_valence_densities*100);
              scale(rho_valence.data(), g.all(), 1. - take_atomic_valence_densities); // mix old
              // add contributions from smooth core densities, and optionally spherical valence densities
              stat += single_atom::atom_update("valence densities", na, 0, nr2.data(), ar2.data(), atom_rhoc.data());
              stat += potential_generator::add_smooth_quantities(rho_valence.data(), g, na, nr2.data(), ar2.data(),
                                    center, n_periodic_images, periodic_images, atom_rhoc.data(),
                                    echo, 0, Y00sq*take_atomic_valence_densities, "smooth valence density");
          } // take_atomic_valence_densities

          if (echo > 4) print_stats(rho_valence.data(), g.all(), g.dV(), "# valence density");

          set(rho.data(), g.all(), rho_valence.data());

          if (echo > 4) print_stats(rho.data(), g.all(), g.dV(), "# density before adding smooth core densities:");

          stat += single_atom::atom_update("core densities", na, 0, nr2.data(), ar2.data(), atom_rhoc.data());
          stat += single_atom::atom_update("qlm charges", na, 0, 0, 0, atom_qlm.data());
          // add contributions from smooth core densities
          stat += potential_generator::add_smooth_quantities(rho.data(), g, na, nr2.data(), ar2.data(),
                                center, n_periodic_images, periodic_images, atom_rhoc.data(),
                                echo, 0, Y00sq, "smooth core density");

          if (echo > 4) print_stats(rho.data(), g.all(), g.dV(), "# density  after adding smooth core densities:");

          here;

          { // scope: eval the XC potential and energy
              double E_xc{0}, E_dc{0};
              for (size_t i = 0; i < g.all(); ++i) {
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
              for (int ia = 0; ia < na; ++ia) {
                  double const sigma = sigma_cmp[ia];
                  int    const ellmax = lmax_qlm[ia];
                  if (echo > 6) std::printf("# use generalized Gaussians (sigma= %g %s, lmax=%d) as compensators for atom #%i\n", sigma*Ang, _Ang, ellmax, ia);
                  std::vector<double> coeff;
                  if ('f' == compensator_method) { // "factorizable"
                      coeff.resize(sho_tools::nSHO(ellmax), 0.0);
                      stat += sho_projection::denormalize_electrostatics(coeff.data(), atom_qlm[ia], ellmax, sigma, unitary, echo);
#ifdef    DEVEL
//                    if (echo > 7) std::printf("# before SHO-adding compensators for atom #%i coeff[000] = %g\n", ia, coeff[0]);
#endif // DEVEL
                  } // factorizable
                  for (int ii = 0; ii < n_periodic_images; ++ii) {
                      double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                      if ('f' == compensator_method) { // "factorizable"
                          stat += sho_projection::sho_add(cmp.data(), g, coeff.data(), ellmax, cnt, sigma, 0);
                      } else {
                          stat += potential_generator::add_generalized_Gaussian(cmp.data(), g, atom_qlm[ia], ellmax, cnt, sigma, echo);
//                           double const one[] = {1};
//                           stat += potential_generator::add_generalized_Gaussian(cmp.data(), g, one, 0, cnt, sigma, echo);
                      }
                  } // periodic images
#ifdef    DEVEL
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
//                SimpleTimer timer(__FILE__, __LINE__, "Poisson equation", echo);
                  if (echo > 3) std::printf("\n\n# %s\n# Solve Poisson equation\n# %s\n\n", h_line, h_line);
                  stat += poisson_solver::solve(Ves.data(), rho.data(), g, es_solver_method, echo, (na > 0)?center[0]:nullptr);
              } // scope
              here;

              // probe the electrostatic potential in real space, find ves_multipoles
              for (int ia = 0; ia < na; ++ia) {
                  double const sigma = sigma_cmp[ia];
                  int    const ellmax = lmax_vlm[ia];
                  int const nc = sho_tools::nSHO(ellmax);
                  std::vector<double> coeff(nc, 0.0);
                  for (int ii = 0; ii < n_periodic_images; ++ii) {
                      std::vector<double> coeff_image(nc, 0.0);
                      double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                      stat += sho_projection::sho_project(coeff_image.data(), ellmax, cnt, sigma, Ves.data(), g, 0);
                      // alternative to projecting Ves we could use Vtot (constructed later) but we need consistency with the operations inside the spheres
                      add_product(coeff.data(), nc, coeff_image.data(), 1.0); // need phase factors? no since we only test the electrostatic
                  } // periodic images
                  // SHO-projectors are brought to the grid unnormalized, i.e. p_{000}(0) = 1.0 and p_{200}(0) = -.5

                  stat += sho_projection::renormalize_electrostatics(atom_vlm[ia], coeff.data(), ellmax, sigma, unitary, echo);
#ifdef    DEVEL
                  if (echo > 3) {
                      std::printf("# potential projection for atom #%d v_00 = %.9f %s\n", ia, atom_vlm[ia][00]*Y00*eV,_eV);
                      int const ellmax_show = std::min(ellmax, 2);
                      for (int ell = 1; ell <= ellmax_show; ++ell) {
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

          if (echo > 1) print_stats(Vtot.data(), g.all(), 0, "\n# Total effective potential (before adding zero potentials)", eV, _eV);

          // now also add the zero potential vbar to Vtot
          stat += potential_generator::add_smooth_quantities(Vtot.data(), g, na, nr2.data(), ar2.data(),
                                center, n_periodic_images, periodic_images, atom_vbar.data(),
                                echo, 0, Y00, "zero potential");

          if (echo > 1) print_stats(Vtot.data(), g.all(), 0,   "# Total effective potential  (after adding zero potentials)", eV, _eV);
          here;


          //
          //  Potential generation done
          //


          double double_counting_correction{0};

          { // scope: solve the Kohn-Sham equation with the given Hamiltonian
              SimpleTimer KS_timer(__FILE__, __LINE__, "solving KS-equation", echo);

#ifdef    DEVEL
              if (echo > 0) {
                  std::printf("\n\n# %s\n# Solve Kohn-Sham equation\n# %s\n\n", h_line, h_line);
                  std::fflush(stdout);
              } // echo
#endif // DEVEL

              view2D<double> rho_valence_new(2, g.all(), 0.0); // new valence density and response
              data_list<double> atom_rho_new[2];
              atom_rho_new[0] = data_list<double>(n_atom_rho, 0.0); // new valence density matrices
              atom_rho_new[1] = data_list<double>(n_atom_rho, 0.0); // and valence response matrices
              double charges[3] = {0, 0, 0}; // 0:kpoint_denominator, 1:charge, 2:response

              KS.solve(rho_valence_new, atom_rho_new, charges, Fermi,
                        g, Vtot.data(), na, atom_mat, occupation_method, scf_iteration, echo);

              double band_energy_sum = Fermi.get_band_sum(); // non-const since we might need to correct this

              Fermi.correct_Fermi_level(nullptr, echo); // clear the accumulators

              // correct if k-point weights do not sum up to 1
              if (charges[0] > 0 && std::abs(charges[0] - 1.0) > 1e-15) {
                  double const renormalization_factor = 1./charges[0];
                  if (echo > 2) std::printf("# %s: renormalize density and density matrices by %.15f\n", __func__, renormalization_factor);
                  scale(rho_valence_new[0], g.all(), renormalization_factor);
                  scale(rho_valence_new[1], g.all(), renormalization_factor);
                  for (int ia = 0; ia < na; ++ia) {
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
                      for (int ia = 0; ia < na; ++ia) {
                          add_product(atom_rho_new[0][ia], n_atom_rho[ia], atom_rho_new[1][ia], alpha);
                      } // ia
                      charges[1] += charges[2]*alpha;
                      band_energy_sum += Fermi.get_band_sum(1)*alpha;
                      Fermi.set_Fermi_level(Fermi.get_Fermi_level() + alpha, echo);
                  } else {
                      auto const n_e = Fermi.get_n_electrons();
                      warn("in SCF iteration #%i number of valence electrons %g deviates from %g by %g but response is zero",
                               scf_iteration, charges[1], n_e, charges[1] - n_e);
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
              for (int ia = 0; ia < na; ++ia) {
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
              for (int ia = 0; ia < na; ++ia) {
                  add_product(Ea.data(), nE, atom_contrib[ia], 1.0); // all atomic weight factors are 1.0
              } // ia
              if (echo > 7) std::printf("\n# sum of atomic energy contributions without grid contributions:\n");
              energy_contribution::show(Ea.data(), echo - 7, eV, _eV);
              // now add grid contributions
              Ea[energy_contribution::KINETIC] += grid_kinetic_energy * (1. - take_atomic_valence_densities);
              Ea[energy_contribution::EXCHANGE_CORRELATION] += grid_xc_energy;
              Ea[energy_contribution::ELECTROSTATIC] += grid_electrostatic_energy;
              // reconstruct total energy from its contributions KINETIC + ES + XC
              Ea[energy_contribution::TOTAL] = Ea[energy_contribution::KINETIC]
                                             + Ea[energy_contribution::ELECTROSTATIC]
                                             + Ea[energy_contribution::EXCHANGE_CORRELATION];
              if (echo > 7) std::printf("\n# energy contributions including grid contributions:\n");
              energy_contribution::show(Ea.data(), echo - 3, eV, _eV);
          } else {
              stat += single_atom::atom_update("energies", na, atomic_energy_diff.data());
          } // total_energy_details

          double atomic_energy_corrections{0};
          for (int ia = 0; ia < na; ++ia) {
              atomic_energy_corrections += atomic_energy_diff[ia];
          } // ia
          total_energy = grid_kinetic_energy * (1. - take_atomic_valence_densities)
                       + grid_xc_energy
                       + grid_electrostatic_energy
                       + atomic_energy_corrections;
          if (echo > 0) { std::printf("\n# total energy %.9f %s\n\n", total_energy*eV, _eV); std::fflush(stdout); }



          if (spherical_valence_decay > 0) { // update take_atomic_valence_densities
#if       1
              auto const x = scf_iteration/spherical_valence_decay; // progress x
              take_atomic_valence_densities = (x >= 1) ? 0 : 2*pow3(x) - 3*pow2(x) + 1; // smooth transition function
              if (echo > 0) std::printf("# set take_atomic_valence_densities = %g %%\n", take_atomic_valence_densities*100);
              stat += single_atom::atom_update("lmax qlm", na, 0, lmax_qlm.data(), &take_atomic_valence_densities);
#endif // 1
          } // spherical_valence_decay
          here;



          { // scope: read the stop file with standard name, max_scf_iterations may be modified
              auto const stat = manage_stop_file(max_scf_iterations, 'r', echo);
              if (stat) warn("failed to read the stop file, status= %i", int(stat));
              if (max_scf_iterations_input != max_scf_iterations) {
                  warn("the max. number of SCF iterations has been modified from %d to %d "
                       "by the stop file during SCF iteration #%i",
                       max_scf_iterations_input, max_scf_iterations, scf_iteration);
              } // stop file has been used
          } // scope
          density_mixing = density_mixing_fixed;

      } // scf_iteration

      here;

      { // scope: delete the stop file
          auto const stat = manage_stop_file(max_scf_iterations, 'd', echo);
          if (stat) warn("failed to delete stop file, status= %i", int(stat));
      } // scope

      KS.store(control::get("store.waves", ""), echo);

#ifdef    DEVEL
      stat += potential_generator::potential_projections(g, Ves.data(), Vxc.data(), Vtot.data(), rho.data(), cmp.data(),
                  na, &center, rcut, echo);
#endif // DEVEL

      stat += single_atom::atom_update("memory cleanup", na);

      solid_harmonics::cleanup();

      return stat;
  } // SCF

#ifdef    NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

  status_t test_scf(int const echo=3) {
      return SCF(echo);
  } // test_scf

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_scf(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace self_consistency
