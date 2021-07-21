#include <cstdio> // std::printf
#include <cassert> // assert
#include <cmath> // std::sqrt, std::pow, std::exp, std::abs, std::sqrt, std::round
#include <fstream> // std::ifstream, std::ofstream
#include <algorithm> // std::min, std::max
#include <iomanip> // std::setprecision
#include <vector> // std::vector<T>

#include "atom_core.hxx" // guess_energy, nl_index, ellchar 

#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // ::create_default_radial_grid, ::destroy_radial_grid
#include "display_units.h" // eV, _eV
#include "radial_potential.hxx" // ::Hartree_potential
#include "exchange_correlation.hxx" // ::LDA_kernel
#include "inline_math.hxx" // pow2, set, scale, product, add_product, dot_product
#include "radial_eigensolver.hxx" // ::shooting_method
#include "constants.hxx" // ::pi
#include "control.hxx" // ::get
#include "lossful_compression.hxx" // RDP_lossful_compression
#include "sigma_config.hxx" // ::get, element_t
#include "recorded_warnings.hxx" // warn
#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

// #define FULL_DEBUG
// #define DEBUG

#ifdef FULL_DEBUG
    #define DEBUG
    #define full_debug(print) print
#else
    #define full_debug(print)
#endif

// #ifdef  DEBUG
    #include "debug_output.hxx" // dump_to_file
// #endif

#ifdef DEBUG
    #define debug(print) print
#else
    #define debug(print)
#endif

namespace atom_core {

  enum energy_contributions {
      E_tot = 0, // total energy: E_eig - E_Htr - E_vxc + E_exc
      E_kin,     // kinetic energy: E_eig - 2 E_Htr - E_vxc - E_Cou
      E_exc,     // exchange correlation energy: 4\pi int r^2 dr E_xc(r) * rho(r)
      E_vxc,     // double counting correction:  4\pi int r^2 dr V_xc(r) * rho(r)
      E_eig,     // contribution from eigenvalues: sum_i occ_i * eig_i
      E_Htr,     // Hartree energy: 1/2 4\pi int r^2 dr V_Htr(r) * rho(r)
      E_Cou,     // Coulomb energy: -Z * 4\pi int r dr rho(r)
      E_est,     // electrostatic: E_Htr + E_Cou
      NumEnergyContributions };

  void rad_pot(
        double rV[] // result: r*V(r)
      , radial_grid_t const & g // radial grid descriptor
      , double const rho4pi[] // 4*\pi*rho(r)
      , double const Z // number of protons
      , double energies[] // =nullptr // energy contribution break down
  ) {
      double Eexc{0}, Evxc{0}, EHtr{0};
      double const ECou = -Z*radial_potential::Hartree_potential(rV, g, rho4pi); // set rV to the Hartree potential
      double const fpi = .25/constants::pi; // 1/(4*pi)
      for (int ir = 0; ir < g.n; ++ir) {
          EHtr += rho4pi[ir]*rV[ir]*g.rdr[ir];
          double Vxc{0};
          double const Exc = exchange_correlation::LDA_kernel(fpi*rho4pi[ir], Vxc);
          rV[ir] += g.r[ir]*Vxc; // add the exchange correlation potential
          Eexc += Exc*rho4pi[ir]*g.r2dr[ir]; // exchange correlation energy
          Evxc += Vxc*rho4pi[ir]*g.r2dr[ir]; // double counting correction
          rV[ir] -= Z; // add the Coulomb potential -Z/r to V(r) in r*V(r) representation
      } // ir
      EHtr *= 0.5;
      if (nullptr != energies) {
          energies[E_vxc] = Evxc;
          energies[E_exc] = Eexc;
          energies[E_Cou] = ECou;
          energies[E_Htr] = EHtr;
      } // export energy contributions
  } // rad_pot

  double initial_density( // returns total charge
        double r2rho[] // result: r^2*rho_initial(r)
      , radial_grid_t const & g // radial grid descriptor
      , double const Z // number of protons
      , double const charged // =0 // excess electrons
  ) {
      auto const alpha = 0.3058*std::pow(Z, 1./3.);
      auto const beta = std::sqrt(108./constants::pi)*std::max(0., Z - 2)*alpha;
      if (4 + 3.2*charged < 0) return -1; // error
      auto const gamma = std::sqrt(4 + 3.2*charged);
      auto const g2 = pow3(gamma) * ((Z < 2) ? 0.5*Z : 1);
      double q{0};
      for (auto ir = 0; ir < g.n; ++ir) {
          auto const r = g.r[ir];
          auto const x = alpha*r;
          auto const exa = (x < 50)? std::exp(-3*x) : 0;
          auto const xb = gamma*r;
          auto const exb = (xb < 150)? std::exp(-xb) : 0;
          r2rho[ir] = beta*std::sqrt(x)*exa + g2*exb*r*r;
  //      std::printf("%g %g\n", r, r2rho[ir]);
          q += r2rho[ir] * g.dr[ir]; // integrate total charge
      } // ir
      return q; // charge
  } // initial_density

  status_t read_Zeff_from_file(
        double Zeff[]
      , radial_grid_t const & g
      , double const Z
      , char const basename[]
      , double const factor
      , int const echo
      , char const prefix[]
  ) {
      status_t stat(0);
      char filename[96];
      get_Zeff_file_name(filename, basename, Z);
      if (echo > 3) std::printf("# %s %s Z=%g  try to read file %s\n", prefix, __func__, Z, filename);
      std::ifstream infile(filename);
      int ngr{0}, ngu{0}; // number of radial grid points read and used
      if (infile.is_open()) {
          double r_min{9e9}, r_max{-9e9};
          int ir{0};
          double r, Ze, r_prev{0}, Ze_prev = Z*factor;
          while (infile >> r >> Ze) {
              ++ngr;
              if (r >= 0) {
                  r_min = std::min(r, r_min);
                  r_max = std::max(r, r_max);
                  if (r <= g.rmax + 1e-6) {
                      full_debug(std::printf("# %s r=%g Zeff=%g\n",  __func__, r, Ze));
                      Ze *= factor;
                      while ((g.r[ir] < r) && (ir < g.n)) {
                        // interpolate linearly
                        Zeff[ir] = Ze_prev + (Ze - Ze_prev)*(g.r[ir] - r_prev)/std::max(r - r_prev, 1e-24);
                        ++ir;
                      } // while
                      ++ngu;
                  } // r <= rmax
              } // r >= 0
              r_prev = r; Ze_prev = Ze; // pass
              stat = (r_max < r_min);
          } // while
          if (echo > 3) std::printf("# %s %s Z=%g  use %d of %d values from file %s, interpolate to %d values\n",
                                  prefix, __func__, Z, ngu, ngr, filename, ir);
      } else {
          if (echo > 1) std::printf("# %s %s Z=%g  failed to open file %s\n", prefix, __func__, Z, filename);
          stat = -1; // failure
      } // is_open
      return stat;
  } // read_Zeff_from_file

  status_t store_Zeff_to_file(
        double const Zeff[]
      , double const r[]
      , int const nr
      , double const Z
      , char const basename[]
      , double const factor
      , int const echo
  ) {
      char filename[96];
      get_Zeff_file_name(filename, basename, Z);
      if (echo > 3) std::printf("# %s  Z=%g  try to write file %s\n",  __func__, Z, filename);
      std::ofstream outfile(filename);
      if (outfile.is_open()) {
          outfile << std::setprecision(15);
          for (int ir = 0; ir < nr; ++ir) {
              outfile << r[ir] << " " << Zeff[ir]*factor << "\n";
          } // write to file
          return 0; // success
      } else {
          if (echo > 1) std::printf("# %s  Z=%g  failed to open file %s for writing\n",  __func__, Z, filename);
          return 1; // failure
      } // is_open
  } // store_Zeff_to_file


  typedef struct {
      double E; // energy
      double occ; // occupation number
      enn_QN_t enn; // principal quantum number
      ell_QN_t ell; // angular momentum
  } orbital_t;


  status_t scf_atom(
        radial_grid_t const & g // radial grid descriptor
      , double const Z // atomic number
      , int const echo // log output level
      , double const occupations[][2] // =nullptr
  ) {
      status_t stat(0);
      if (echo > 7) std::printf("\n# %s:%d  %s \n\n", __FILE__, __LINE__, __func__);

      int constexpr sra = 1; // scalar-relativistic
      int constexpr MAXCYCLES = 200;
      int constexpr MINCYCLES = 3;
      double constexpr THRESHOLD = 1e-11;

      int imax{-1};
      assert(Z <= 120);
      orbital_t orb[20];
      {   // prepare orbitals
          int iZ{0}; // integer atomic number needed for occupations
          int i{0}; // init shell index i
          for (int m = 0; m < 8; ++m) { // auxiliary quantum number
              int enn{(m + 1)/2};
              for (int ell = m/2; ell >= 0; --ell) {
                  ++enn; // main quantum number
                  int const max_occ = 2*ell + 1;
                  orb[i].enn = enn;
                  orb[i].ell = ell;
                  int const inl = nl_index(enn, ell);
                  auto const auto_occ = std::min(std::max(0.0, Z - iZ), 2.*max_occ);
                  orb[i].occ = occupations ? std::min(2.*max_occ, std::abs(occupations[inl][0]) + std::abs(occupations[inl][1])) : auto_occ;
                  if (std::abs(orb[i].occ - auto_occ) > 2e-16) {
                      warn("For Z=%g the occupation of the %d%c-orbital differs: %g vs %g (auto)",
                                      Z, enn, ellchar(ell), orb[i].occ, auto_occ);
                  } // warning
                  orb[i].E = guess_energy(Z, enn);
                  if (orb[i].occ > 0) {
                      if (echo > 6) std::printf("# %s  Z=%g  i=%d %d%c f= %g  guess E= %g %s\n",
                           __func__, Z, i, enn, ellchar(ell), orb[i].occ, orb[i].E*eV, _eV);
                      imax = std::max(i, imax);
                  } // occupied
                  iZ += 2*max_occ; // iZ jumps between atomic numbers of atoms with full shells
                  ++i; // next shell
              } // ell
          } // m
      } // orb
      double previous_eigenvalues[20];

      std::vector<double> r2rho(g.n);
      std::vector<double> rho4pi(g.n);
      std::vector<double> r2rho4pi(g.n);
      std::vector<double> rV_new(g.n, -Z); // init as Hydrogen-like potential
      std::vector<double> rV_old(g.n, -Z); // init as Hydrogen-like potential
      
      double mix{0}; // init mixing coefficient for the potential
      double const alpha = 0.33; // limit case potential mixing coefficient

      double energies[NumEnergyContributions];
      double eigenvalue_sum{0};
      double previous_energy{0};

      enum { Task_Solve, Task_ChkRho, Task_GenPot, Task_MixPot, Task_Energy } next_task{Task_Solve};

      int icyc{0};
      { // start scope
          bool loading_failed = (0 != int(read_Zeff_from_file(rV_old.data(), g, Z, "pot/Zeff", -1.)));
          full_debug(dump_to_file("rV_loaded.dat", g.n, rV_old.data(), g.r));

          if (Z != std::round(Z))

          if (loading_failed) {
              if (Z != std::round(Z)) {
                  // maybe loading failed because there is no file for non-integer core charges
                  loading_failed = (0 != int(read_Zeff_from_file(rV_old.data(), g, std::round(Z), "pot/Zeff", -1)));
              }
              if (loading_failed) {
                // use guess density if loading failed
                auto const q = initial_density(r2rho4pi.data(), g, Z);
                if (echo > 2) std::printf("# %s  Z=%g  guess rho with %.6f electrons\n",  __func__, Z, q);
                next_task = Task_ChkRho; // different entry point into the loop, skip the 1st solver call
              } // loading_failed (still)
          } // loading_failed
      } // start scope

      double res{9e9}; // residual
      bool run{true};
      while (run) {
//           run = ((res > THRESHOLD) || (icyc <= MINCYCLES)) && (icyc < MAXCYCLES);

          switch (next_task) {
              ///////////////////////////////////////////////////////
              case Task_Solve: {
                  full_debug(std::printf("# Task_Solve\n"));

                  set(r2rho4pi.data(), g.n, 0.0); // init accumulator density
                  eigenvalue_sum = 0; // reset energy accumulator

                  for (int i = 0; i <= imax; ++i) {
                      if (orb[i].occ > 0) {
                          previous_eigenvalues[i] = orb[i].E; // copy
                          auto const stat = radial_eigensolver::shooting_method(sra, g, rV_old.data(), orb[i].enn, orb[i].ell, orb[i].E, nullptr, r2rho.data());
                          if (stat) {
                              std::printf("# %s  Z=%g  failed solving for %d%c, status = %d\n",  __func__, Z, orb[i].enn, ellchar(orb[i].ell), int(stat));
                              return stat;
                          } // failed
                          if (echo > 6) {
                              std::printf("# %s  Z=%g  %d%c E=%15.6f %s\n",  __func__, Z, orb[i].enn, ellchar(orb[i].ell), orb[i].E*eV, _eV);
                          } // echo
                          // add orbital density
                          double const q = dot_product(g.n, r2rho.data(), g.dr);
                          assert(q > 0);
                          double const f = orb[i].occ/q;
                          add_product(r2rho4pi.data(), g.n, r2rho.data(), f);
                          if (echo > 9) {
                              std::printf("# %s  %d%c f= %g q= %g dE= %g %s\n",
                                     __func__, orb[i].enn, ellchar(orb[i].ell), orb[i].occ, q,
                                     (orb[i].E - previous_eigenvalues[i])*eV,_eV);
                          } // echo
                          // add orbital energy
                          eigenvalue_sum += orb[i].occ * orb[i].E;
                      } // occupied
                  } // i, orbitals

                  next_task = Task_ChkRho;
              } break; // Task_Solve
              ///////////////////////////////////////////////////////
              case Task_ChkRho: {
                  full_debug(std::printf("# Task_Check_Rho\n"));

                  product(rho4pi.data(), g.n, r2rho4pi.data(), g.rinv, g.rinv);

                  double const q  = dot_product(g.n, rho4pi.data(), g.r2dr); // can be merged with the loop above
                  double const qq = dot_product(g.n, r2rho4pi.data(), g.dr); // can be merged with the loop above
                  if (std::abs(q - Z) > .1 && echo > 0) {
                      warn("for Z=%g rho has %g (or %g) electrons", Z, q, qq);
                  } // not charge neutral

                  next_task = Task_GenPot;
              } break; // Task_ChkRho
              ///////////////////////////////////////////////////////
              case Task_GenPot: {
                  full_debug(std::printf("# Task_Generate_Pot\n"));

                  rad_pot(rV_new.data(), g, rho4pi.data(), Z, energies);
                  debug({ char fn[99]; std::sstd::printf(fn, "rV_new-icyc%d.dat", icyc); dump_to_file(fn, g.n, rV_new, g.r); });

                  next_task = Task_Energy;
              } break; // Task_GenPot
              ///////////////////////////////////////////////////////
              case Task_Energy: {
                  full_debug(std::printf("# Task_Total_Energy\n"));

                  energies[E_eig] = eigenvalue_sum;
                  energies[E_kin] = energies[E_eig] - energies[E_Htr]*2
                                  - energies[E_Cou] - energies[E_vxc]; // total kinetic energy
                  energies[E_est] = energies[E_Htr] + energies[E_Cou]; // total electrostatic energy
                  energies[E_tot] = energies[E_eig] - energies[E_Htr]
                                  + energies[E_exc] - energies[E_vxc]; // total energy

                  ++icyc; // transit to the next SCF iteration
                  auto const which = E_est; // monitor the change on some energy contribution which depends on the density only
                  res = std::abs(energies[which] - previous_energy);
                  previous_energy = energies[which]; // store for the next iteration
                  if (echo > 3) {
                      int const display_every = 1 << std::max(0, 2*(6 - echo)); // echo=4:every 16th, 5:every 4th, 6:every cycle
                      if (0 == (icyc & (display_every - 1))) {
                          std::printf("# %s  Z=%g  icyc=%d  residual=%.1e  E_tot=%.9f %s\n",
                                  __func__, Z, icyc, res, energies[E_tot]*eV, _eV);
                      } // display?
                  } // echo
                  run = ((res > THRESHOLD) || (icyc <= MINCYCLES)) && (icyc < MAXCYCLES);

                  next_task = Task_MixPot;
              } break; // Task_Energy
              ///////////////////////////////////////////////////////
              case Task_MixPot: {
                  debug(std::printf("# Task_Mix_Potentials with %.3f %%\n", run?100*mix:0.));

                  if (run) {
                      scale(rV_old.data(), g.n, 1. - mix);
                      add_product(rV_old.data(), g.n, rV_new.data(), mix);
                  } else {
                      mix = 0; // do not mix in the last iteration
                      run = true; // run one more iteration to do this task: 'solve'
                  } // last iteration

                  mix = alpha * icyc/(icyc + .5*Z + 5.); // set the mixing coefficient for the next iteration

                  next_task = Task_Solve;
              } break; // Task_MixPot
              ///////////////////////////////////////////////////////
          } // next_task

      } // while run

      if (res > THRESHOLD) {
          if (echo > 0) std::printf("# %s  Z=%g  Warning! Did not converge in %d iterations, res=%.1e\n", __func__, Z, icyc, res);
          stat += 9;
      } else {
          if (echo > 1) {
              std::printf("# %s  Z=%g  converged in %d iterations to res=%.1e, E_tot= %.9f %s\n",
                      __func__, Z, icyc, res, energies[E_tot]*eV, _eV);
              std::printf("# %s  Z=%g  E_kin= %.9f E_xc= %.9f E_es= %.9f %s\n", __func__, Z, 
                      energies[E_kin]*eV, energies[E_exc]*eV, energies[E_est]*eV, _eV);
              std::printf("# %s  Z=%g  E_Coulomb= %.9f E_Hartree= %.9f %s\n", __func__, Z, 
                      energies[E_Cou]*eV, energies[E_Htr]*eV, _eV);
          } // echo
          if (echo > 2) {
              for (int i = 0; i <= imax; ++i) {
                  std::printf("# %s  Z=%g  %d%c  E=%15.6f %s  f=%g\n",  __func__, Z, orb[i].enn, ellchar(orb[i].ell), orb[i].E*eV, _eV, orb[i].occ);
              } // i
//            std::printf("\n"); // blank line after displaying spectrum
          } // echo
          stat += store_Zeff_to_file(rV_old.data(), g.r, g.n, Z, "pot/Zeff", -1);
      } // converged?

      // dump_to_file("rV_converged.dat", g.n, rV_old, g.r);
      // ToDo: export the converged potential rV_old and the position of energies
#if 0
      for (int i = 0; i < 20; ++i) {
          if (orb[i].occ <= 0) {
              // solve for the unoccupied states
              auto const stat = radial_eigensolver::shooting_method(sra, g, rV_old, orb[i].enn, orb[i].ell, orb[i].E);
              if (0 == stat) {
                  if (echo > 2)
                      std::printf("# %s  Z=%g  %d%c  E=%15.6f %s  f=%g\n",  __func__, Z, orb[i].enn, ellchar(orb[i].ell), orb[i].E*eV, _eV, orb[i].occ);
              } else { // worked
                  if (echo > 0)
                      std::printf("# %s  Z=%g  %d%c  E=%15.6f %s  f=%g  status=%d\n",  __func__, Z, orb[i].enn, ellchar(orb[i].ell), orb[i].E*eV, _eV, orb[i].occ, stat);
              }
          } // unoccupied
      } // i
#endif

      int const show_state_diagram = control::get("show.state.diagram", 0.); // 0:do not show, 1:show
      if (show_state_diagram) {
          float state_diagram[30][4]; // [inl][start,end,energy,spare]
          for (int inl = 0; inl < 30; ++inl) { 
              state_diagram[inl][0] = 0; // start
              state_diagram[inl][1] = 0; // end
              state_diagram[inl][2] = -1e5;
          } // inl

          for (int i = 0; i <= imax; ++i) {
              if (orb[i].occ > 0) {
                  radial_eigensolver::shooting_method(sra, g, rV_old.data(), orb[i].enn, orb[i].ell, orb[i].E, nullptr, r2rho.data());
                  double const norm = dot_product(g.n, r2rho.data(), g.dr);
                  assert(norm > 0);
                  scale(r2rho.data(), g.n, 1./norm);
                  int const inl = atom_core::nl_index(orb[i].enn, orb[i].ell);
                  // find the radius where the state starts (95%) and where it ends (95%) such that 90% are inside that range
                  double q{0}; int ir{1};
                  while (q < .05) { q += r2rho[ir]*g.dr[ir]; ++ir; }
                  state_diagram[inl][0] = g.r[ir];
                  q = 0; ir = g.n;
                  while (q < .05) { --ir; q += r2rho[ir]*g.dr[ir]; }
                  state_diagram[inl][1] = g.r[ir];
                  state_diagram[inl][2] = orb[i].E;
              } // occupied
          } // i, orbitals

          std::printf("## state diagram for Z=%g in %s (sqrt on radius, log10 on energies)\n", Z, _eV);
          for (int inl = 0; inl < 30; ++inl) {
              std::printf("%g %g\n%g %g\n\n", std::sqrt(state_diagram[inl][0]*Ang), -std::log10(-state_diagram[inl][2]*eV), 
                                         std::sqrt(state_diagram[inl][1]*Ang), -std::log10(-state_diagram[inl][2]*eV)); // energy level
          } // inl
          std::printf("## potential state diagram for Z=%g in %s\n", Z, _eV);
          for (int ir = 1; ir < g.n; ++ir) {
              double const pot = rV_old[ir]*g.rinv[ir];
              if (pot > -1e5 && pot < 0) std::printf("%g %g\n", std::sqrt(g.r[ir]*Ang), -std::log10(-pot*eV));
          } // ir
      } // show_state_diagram

      if (echo > 2) std::printf("\n"); // blank line after atoms is computed
      
      return stat;
  } // scf_atom



  status_t solve(
        double const Z
      , int const echo // =0
      , char const config // ='a'
      , radial_grid_t const *rg // =nullptr
  ) {
      bool const my_radial_grid = (nullptr == rg);
      radial_grid_t       *const m = my_radial_grid ? radial_grid::create_default_radial_grid(Z) : nullptr;
      radial_grid_t const *const g = my_radial_grid ? m : rg;
      if ('a' == config) { // "auto"
          return scf_atom(*g, Z, echo);
      } else {
          auto const element = sigma_config::get(Z, echo);
          return scf_atom(*g, Z, echo, element.occ);
      }
      if (my_radial_grid) radial_grid::destroy_radial_grid(m);
  } // solve

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_neutral_atom_total_energy(int const echo=0) {
      if (echo < 7) return 0;
      std::printf("\n\n## %s\n", __func__);
      for (int iZ = -1; iZ <= 120*2; ++iZ) {
          double const Z = iZ*0.5;
          std::printf("%g %.9f\n", Z, neutral_atom_total_energy(Z));
      } // Z
      std::printf("\n# %s\n\n", __func__);
      return 0;
  } // test_neutral_atom_total_energy

  status_t test_initial_density(int const echo=0) {
      if (echo > 3) std::printf("\n# %s:%d  %s \n\n", __FILE__, __LINE__, __func__);
      double maxdev{0};
      for (double Z = 0; Z < 128; Z += 1) {
          auto & g = *radial_grid::create_default_radial_grid(Z);
          std::vector<double> r2rho(g.n, 0.0);
          double const q = initial_density(r2rho.data(), g, Z);
          double const dev = Z - q;
          if (echo > 5) std::printf("# %s:%d Z = %g charge = %.3f electrons, diff = %g\n",
                                  __FILE__, __LINE__, Z, q, dev);
          maxdev = std::max(maxdev, std::abs(dev));
          radial_grid::destroy_radial_grid(&g);
      } // Z
      if (echo > 1) std::printf("# %s: max. deviation of %g electrons\n", __func__, maxdev);
      return (maxdev > 2e-5);
  } // test_initial_density

  status_t test_nl_index(int const echo=2) {
      int inl{0};
      for (int enn = 1; enn < 9; ++enn) { // principal quantum number of an atom
          for (int ell = 0; ell < enn; ++ell) { // angular momentum quantum number
              int const nl = nl_index(enn, ell);
              if ((inl != nl) || (echo > 8)) std::printf("# %s: %s n=%d l=%d (%d%c) nl_index %d %d\n", __FILE__, __func__, enn, ell, enn, ellchar(ell), inl, nl);
              assert(inl == nl);
              ++inl;
          } // ell
      } // enn
      return 0;
  } // test_nl_index

  status_t test_core_solver(int const echo=3) {
      double const Z_begin = control::get("atom_core.test.Z", 29.); // default copper
      double const Z_inc   = control::get("atom_core.test.Z.inc", 1.); // default: sample only integer values
      double const Z_end   = control::get("atom_core.test.Z.end", Z_begin + Z_inc); // default: only one core
      if (echo > 0) std::printf("\n# %s:%d  %s(echo=%d) from %g to %g in steps of %g\n\n",
                      __FILE__, __LINE__, __func__, echo, Z_begin, Z_end - Z_inc, Z_inc);
      status_t stat{0};
      char const custom_config = 32 | *control::get("atom_core.occupations", "custom");
      int const i_end = std::round((Z_end - Z_begin)/std::max(1e-9, Z_inc));
      // ToDo: OpenMP loop
      for (int i = 0; i < i_end; ++i) {
          double const Z = Z_begin + i*Z_inc;
          if (echo > 1) std::printf("\n# atom_core solver for Z= %g\n", Z);
          auto const stat_Z = solve(Z, echo, custom_config);
          if (stat_Z) warn("atom_core.test for Z=%g returned status=%i", Z, int(stat_Z));
          stat += stat_Z;
      } // i
      return stat;
  } // test_core_solver

  
  status_t simplify_Zeff_file(
        double const Z
      , float const epsilon=1e-6
      , int const echo=3
  )
    // Apply Ramer-Douglas-Peucker lossful compression
    // to the input files full_pot/Zeff.00Z
    // with  output files      pot/Zeff.00Z
  {
      status_t stat(0);
      auto & g = *radial_grid::create_default_radial_grid(Z);
      std::vector<double> y(g.n, 0.0);
      stat += read_Zeff_from_file(y.data(), g, Z, "full_pot/Zeff", 1., echo);
      auto const mask = RDP_lossful_compression(g.r, y.data(), g.n, epsilon);
      int const new_n = std::count(mask.begin(), mask.end(), true);
      if (echo > 4) std::printf("# Ramer-Douglas-Peucker for Z=%g reduced %d to %d points, ratio= %.3f\n", 
                                    Z, g.n, new_n, g.n/std::max(1., 1.*new_n));
      std::vector<double> new_y(new_n, 0.0), new_r(new_n, 0.0);
      { // scope: compress y(r)
          int i{0};
          for (int ir = 0; ir < g.n; ++ir) {
              if (mask[ir]) {
                  new_r[i] = g.r[ir];
                  new_y[i] = y[ir];
                  ++i;
              } // active
          } // ir
          assert(i == new_n); // after running RDP, there must be exactly new_n true entries left
      } // scope
      stat += std::abs(store_Zeff_to_file(new_y.data(), new_r.data(), new_n, Z, "pot/Zeff", 1., echo));
      radial_grid::destroy_radial_grid(&g);
      return stat;
  } // simplify_Zeff_file

  status_t test_Zeff_file_compression(int const echo=0) {
      float const threshold = control::get("atom_core.compression.threshold", -1.);
      if (threshold < 0) return 0; // do not do anything, silent return
      if (echo > 0) std::printf("\n# %s:%d  %s(echo=%d)\n\n", __FILE__, __LINE__, __func__, echo);
      status_t stat{0};
      for (int Z = 120; Z >= 0; --Z) { // test all atoms, backwards
          stat += simplify_Zeff_file(Z, threshold, echo); // apply Ramer-Douglas-Peucker reduction to Z_eff(r)
      } // Z
      return stat;
  } // test_Zeff_file_compression

  
  status_t all_tests(int const echo) {
      status_t stat(0);
      int n{0}; int const t = control::get("atom_core.select.test", -1.); // -1:all
      if (t & (1 << n++)) stat += test_neutral_atom_total_energy(echo);
      if (t & (1 << n++)) stat += test_initial_density(echo);
      if (t & (1 << n++)) stat += test_nl_index(echo);
      if (t & (1 << n++)) stat += test_core_solver(echo);
      if (t & (1 << n++)) stat += test_Zeff_file_compression(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace atom_core
