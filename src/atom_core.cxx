// #include <vector> // std::vector
#include <cstdio> // printf
#include <cassert> // assert
#include <cmath> // sqrt, pow, exp, fabs, sqrt
#include <fstream> // ifstream

#include "atom_core.hxx"

#include "min_and_max.h" // min, max
#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // create_exponential_radial_grid
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
#include "output_units.h" // eV, _eV
#include "radial_potential.hxx" // Hartree_potential, lda_PZ81_kernel
#include "radial_eigensolver.hxx" // shooting_method

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

typedef int status_t;


namespace atom_core {
  using namespace radial_grid;

  double constexpr C_PI = 3.14159265358979323846; // pi
  
  typedef struct {
      float occ; // occupation number
      enn_QN_t enn; // principal quantum number
      ell_QN_t ell; // angular momentum
      double E; // energy
  } orbital_t;

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

  void rad_pot(double rV[], radial_grid_t const &g, double const rho4pi[], double const Z=0, double *energies=nullptr) {
      double Eexc = 0, Evxc = 0, EHtr = 0;
      double const ECou = -Z*radial_potential::Hartree_potential(rV, g, rho4pi); // set rV to the Hartree potential
      double const fpi = .25/C_PI;
      for(int ir = 0; ir < g.n; ++ir) {
          EHtr += rho4pi[ir]*rV[ir]*g.rdr[ir];
          double Vxc;
          double const Exc = radial_potential::lda_PZ81_kernel(fpi*rho4pi[ir], Vxc);
          rV[ir] += g.r[ir]*Vxc; // add the exchange correlation potential
          Eexc += Exc*rho4pi[ir]*g.r2dr[ir]; // exchange correlation energy
          Evxc += Vxc*rho4pi[ir]*g.r2dr[ir]; // double counting correction
          rV[ir] -= Z; // add the Coulomb potential -Z/r
      } // ir
      EHtr *= 0.5;
      if (nullptr != energies) {
          energies[E_vxc] = Evxc;
          energies[E_exc] = Eexc;
          energies[E_Cou] = ECou;
          energies[E_Htr] = EHtr;
      } // export energy contributions
  } // rad_pot
  
  
  
  char ellchar(ell_QN_t const ell) {
      char constexpr special_ellchars[] = "spd";
      if (ell < 0) return '?';
      if (ell < 3) return special_ellchars[ell];
      return 99 + ell; // "fghijk ..."
  } // ellchar

  double dot_product(int const n, double const bra[], double const ket[]) {
      double d = 0;
      for(int i = 0; i < n; ++i) {
          d += bra[i]*ket[i];
      } // i
      return d;
  } // dot_product

    
  double initial_density(double r2rho[], radial_grid_t const &g, double const Z, double const charged=0) {
    auto const alpha = 0.3058*pow(Z, 1./3.);
    auto const beta = sqrt(108./C_PI)*max(0, Z - 2)*alpha;
    if (4 + 3.2*charged < 0) return 1; // error
    auto const gamma = sqrt(4 + 3.2*charged);
    auto g2 = gamma*gamma*gamma;
    if (Z < 2) g2 *= 0.5*Z;
    double q = 0;
    for(auto ir = 0; ir < g.n; ++ir) {
        auto const r = g.r[ir];
        auto const x = alpha*r;
        auto const exa = (x < 50)? exp(-3*x) : 0;
        auto const xb = gamma*r;
        auto const exb = (xb < 150)? exp(-xb) : 0;
        r2rho[ir] = beta*sqrt(x)*exa + g2*exb*r*r;
//      printf("%g %g\n", r, r2rho[ir]);
        q += r2rho[ir] * g.dr[ir]; // integrate total charge
    } // ir
    return q; // charge
  } // initial_density
  
  
  status_t scf_atom(radial_grid_t const &g, // radial grid descriptor
               float const Z, // atomic number
               int const echo) { // log output level
      debug(printf("\n# %s:%d  %s \n\n", __FILE__, __LINE__, __func__));

      int constexpr sra = 1; // scalar-relativistic
      int constexpr MAXCYCLES = 200; 
      int constexpr MINCYCLES = 3;
      double constexpr THRESHOLD = 1e-11;
      
      int imax = -1;
      assert(Z <= 120);
      orbital_t orb[20]; 
      {   // prepare orbitals
          int iZ = 0; // integer atomic number needed for occupations
          int i = 0; // init shell index i
          for(int m = 0; m < 8; ++m) {
              int enn = (m + 1)/2;
              for(int ell = m/2; ell >= 0; --ell) {
                  ++enn; // main quantum number
                  int const max_occ = 2*(ell + 1 + ell); // spin degenerate
                  orb[i].enn = enn;
                  orb[i].ell = ell;
                  orb[i].occ = min(max(0, Z - iZ), max_occ);
                  orb[i].E = -.5*(Z/enn)*(Z/enn) *  // Hydrogen-like energies in the Hartree unit system
                            (.783517 + 2.5791E-5*(Z/enn)*(Z/enn)) * // fit for the correct 1s energy
                            exp(-.01*(enn - 1)*Z); // guess energy
                  if (orb[i].occ > 0) {
                      if (echo > 4) {
                          printf("# %s  i=%d %d%c f= %g  guess E= %g %s\n", __func__, i, enn, ellchar(ell), orb[i].occ, orb[i].E*eV, _eV);
                      } // echo
                      imax = max(i, imax);
                  } // occupied
                  iZ += max_occ; // iZ jumps between atomic numbers of atoms with full shells
                  ++i; // next shell
              } // ell
          } // m
      } // orb
      double previous_eigenvalues[20];

      auto const r2rho = new double[g.n];
      auto const rho4pi = new double[g.n];
      auto const r2rho4pi = new double[g.n];
      
      auto rV_new = new double[g.n];
      auto rV_old = new double[g.n];
      for(int ir = 0; ir < g.n; ++ir) {
          rV_old[ir] = -Z; // Hydrogen-like potential
          rV_new[ir] = -Z; // Hydrogen-like potential
      } // ir
      double mix = 0; // current mixing coefficient for the potential 
      double const alpha = 0.33; // limit case potential mixing coefficient
      
      double energies[NumEnergyContributions];
      double eigenvalue_sum = 0;
      double previous_energy = 0;

      enum { Task_Solve, Task_ChkRho, Task_GenPot, Task_MixPot, Task_Energy } next_task = Task_Solve;

      int icyc = 0;
      { // start scope
          bool loading_failed = true;
          // ToDo: try to load rV_old from a file
          char const path[] = "pot/"; // ToDo: should be an external argument
          char filename[99]; sprintf(filename, "%s/rV.%03d", path, (int)Z);
          {
              if (echo > 3) printf("# %s  Z=%g  try to read from file %s\n",  __func__, Z, filename);
              std::ifstream infile(filename);
              double r_min = 9e9, r_max = - 9e9;
              int ir = 1;
              double r, rV, r_prev=0, rV_prev=-Z;
              while (infile >> r >> rV) {
                  if (r >= 0) { 
                      r_min = min(r, r_min);
                      r_max = max(r, r_max);
                      if (r <= g.rmax) {
                          full_debug(printf("# %s  r=%g rV=%g\n",  __func__, r, rV));
                          while ((g.r[ir] < r) && (ir < g.n - 2)) {
                            // interpolate
                            rV_old[ir] = rV_prev + (rV - rV_prev)*(g.r[ir] - r_prev)/max(r - r_prev, 1e-24);
                            ++ir;
                          } // while
                      } // r <= rmax
                  } // r >= 0
                  r_prev = r; rV_prev = rV;
              } // input
              loading_failed = (r_max < r_min);
          }
          for(int ir = 1; ir < g.n; ++ir) {
              rV_old[ir] = max(rV_old[ir], rV_old[ir - 1]);
          } // monotoneous
          full_debug(dump_to_file("rV_loaded.dat", g.n, rV_old, g.r));
  
          if (loading_failed) {
              // use guess density if loading failed
              auto const q = initial_density(r2rho4pi, g, Z);
              if (echo > 2) printf("# %s  Z=%g  guess rho with %.6f electrons\n",  __func__, Z, q);
              next_task = Task_ChkRho; // different entry point into the loop, skip the 1st solver call
          } // loading_failed
      } // start scope

      double res = 9e9; // residual
      bool run = true;
      while (run) {
//           run = ((res > THRESHOLD) || (icyc <= MINCYCLES)) && (icyc < MAXCYCLES);

          switch (next_task) {
              ///////////////////////////////////////////////////////
              case Task_Solve: {
                  full_debug(printf("# Task_Solve\n"));
                
                  for(int ir = 0; ir < g.n; ++ir) {
                      r2rho4pi[ir] = 0; // init accumulator density
                  } // ir
                  eigenvalue_sum = 0; // init energy accumulator

                  for(int i = 0; i <= imax; ++i) {
                      if (orb[i].occ > 0) {
                          previous_eigenvalues[i] = orb[i].E; // copy
                          auto const stat = radial_eigensolver::shooting_method(sra, g, rV_old, orb[i].enn, orb[i].ell, orb[i].E, nullptr, r2rho);
                          if (0 != stat) {
                              printf("# %s  Z=%g  failed solving for %d%c\n",  __func__, Z, orb[i].enn, ellchar(orb[i].ell));
                              return stat;
                          } // failed
                          if (echo > 6) {
                              printf("# %s  Z=%g  %d%c E=%15.6f %s\n",  __func__, Z, orb[i].enn, ellchar(orb[i].ell), orb[i].E*eV, _eV);
                          } // echo
                          // add orbital density
                          double const q = dot_product(g.n, r2rho, g.dr);
                          assert(q > 0);
                          double const f = orb[i].occ/q;
                          for(int ir = 0; ir < g.n; ++ir) {
                              r2rho4pi[ir] += f*r2rho[ir]; // DAXPY
                          } // ir
                          if (echo > 9) {
                              printf("# %s  %d%c f= %g q= %g dE= %g %s\n",  
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
                  full_debug(printf("# Task_Check_Rho\n"));
                  
                  for(int ir = 0; ir < g.n; ++ir) {
                      rho4pi[ir] = r2rho4pi[ir]*g.rinv[ir]*g.rinv[ir];
                  } // ir
                  double const q = dot_product(g.n, rho4pi, g.r2dr); // can be merged with the loop above
                  double const qq = dot_product(g.n, r2rho4pi, g.dr); // can be merged with the loop above
                  if (fabs(q - Z) > .1 && echo > 0) {
                      printf("# %s  Z=%g  icyc=%d  Warning: rho has %.6f (or %.6f) electrons\n",  __func__, Z, icyc, q, qq);
                  } // not charge neutral
                  
                  next_task = Task_GenPot;
              } break; // Task_ChkRho
              ///////////////////////////////////////////////////////
              case Task_GenPot: {
                  full_debug(printf("# Task_Generate_Pot\n"));
                  
//                for(int ir = 0; ir < g.n; ++ir) { rV_new[ir] = -Z; } // create a Hydrogen-like potential
                  rad_pot(rV_new, g, rho4pi, Z, energies);
                  debug({ char fn[99]; sprintf(fn, "rV_new-icyc%d.dat", icyc); dump_to_file(fn, g.n, rV_new, g.r); });

                  next_task = Task_Energy;
              } break; // Task_GenPot
              ///////////////////////////////////////////////////////
              case Task_Energy: {
                  full_debug(printf("# Task_Total_Energy\n"));

                  energies[E_eig] = eigenvalue_sum;
                  energies[E_kin] = energies[E_eig] - energies[E_Htr]*2
                                  - energies[E_Cou] - energies[E_vxc];
                  energies[E_est] = energies[E_Htr] + energies[E_Cou];
                  energies[E_tot] = energies[E_eig] - energies[E_Htr] 
                                  + energies[E_exc] - energies[E_vxc];

                  ++icyc;
                  auto const which = E_est; // monitor the change on some energy contribution which depends on the density only
                  res = fabs(energies[which] - previous_energy); 
                  previous_energy = energies[which]; // store for the next iteration
                  if (echo > 3) {
                      int const display_every = 1 << max(0, 2*(6 - echo)); // echo=4:every 16th, 5:every 4th, 6:every cycle
                      if (0 == (icyc & (display_every - 1))) {
                          printf("# %s  Z=%g  icyc=%d  residual=%.1e  E_tot=%.9f %s\n", 
                                  __func__, Z, icyc, res, energies[E_tot]*eV, _eV);
                      } // display?
                  } // echo
                  run = ((res > THRESHOLD) || (icyc <= MINCYCLES)) && (icyc < MAXCYCLES);

                  next_task = Task_MixPot;
              } break; // Task_Energy
              ///////////////////////////////////////////////////////
              case Task_MixPot: {
                  debug(printf("# Task_Mix_Potentials with %.3f %%\n", run?100*mix:0.));
                
                  if (run) {
                      for(int ir = 0; ir < g.n; ++ir) {
                          rV_old[ir] = (1 - mix) * rV_old[ir] + mix * rV_new[ir];
                      } // ir
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
      
      if (echo > 0) {
          if (res > THRESHOLD) {
              printf("# %s  Z=%g  Warning! Did not converge in %d iterations, res=%.1e\n", 
                        __func__, Z, icyc, res);
          } else {
              printf("# %s  Z=%g  converged in %d iterations to res=%.1e, E_tot= %.9f %s\n", 
                      __func__, Z, icyc, res, energies[E_tot]*eV, _eV);
              if (echo > 1) {
                  for(int i = 0; i <= imax; ++i) {
                      printf("# %s  Z=%g  %d%c  E=%15.6f %s  f=%g\n",  __func__, Z, orb[i].enn, ellchar(orb[i].ell), orb[i].E*eV, _eV, orb[i].occ);
                  } // i
                  printf("\n");
              }
          } // converged?
      } // echo
      
      // dump_to_file("rV_converged.dat", g.n, rV_old, g.r);
      // ToDo: export the converged potential rV_old and the position of energies
#if 0      
      for(int i = 0; i < 20; ++i) {
          if (orb[i].occ <= 0) {
              // solve for the unoccupied states
              auto const stat = radial_eigensolver::shooting_method(sra, g, rV_old, orb[i].enn, orb[i].ell, orb[i].E);
              if (0 == stat) {
                  if (echo > 2)
                      printf("# %s  Z=%g  %d%c  E=%15.6f %s  f=%g\n",  __func__, Z, orb[i].enn, ellchar(orb[i].ell), orb[i].E*eV, _eV, orb[i].occ);
              } else { // worked
                  if (echo > 0)
                      printf("# %s  Z=%g  %d%c  E=%15.6f %s  f=%g  status=%d\n",  __func__, Z, orb[i].enn, ellchar(orb[i].ell), orb[i].E*eV, _eV, orb[i].occ, stat);
              }
          } // unoccupied
      } // i
#endif     
     
      return (res > THRESHOLD);
  } // scf_atom
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_initial_density(radial_grid_t const &g, int const echo=0) {
    printf("\n# %s:%d  %s \n\n", __FILE__, __LINE__, __func__);
    auto rho = new double[g.n];
    double diff = 0;
    for(float zz = 0; zz < 128; zz += 1) {
        initial_density(rho, g, zz);
        double q = 0;
        for(auto ir = 0; ir < g.n; ++ir) {
            q += g.r2dr[ir]*rho[ir];
        } // ir
        diff = max(diff, fabs(zz - q));
        if (echo > 0)  printf("# %s:%d Z = %g charge = %.3f electrons, diff = %g\n", 
                                 __FILE__, __LINE__, zz, q, diff);
    } // zz
    return (diff > 1e-3);
  } // test_initial_density

  status_t test_core_solver(radial_grid_t const &g, float const Z) {
    int const echo = 3;
    printf("\n# %s:%d  %s(echo=%d)\n\n", __FILE__, __LINE__, __func__, echo);
    return scf_atom(g, Z, echo);
  } // test_core_solver

  status_t all_tests() {
    auto status = 0;
//  status += test_initial_density(*create_exponential_radial_grid(512));
//     for(int Z = 1; Z < 120; ++Z)
    for(int Z = 29; Z <= 29; ++Z)
        status += test_core_solver(*create_exponential_radial_grid(250*sqrt(Z + 9.)+.5), Z);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace atom_core
