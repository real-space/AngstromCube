// #include <vector> // std::vector
#include <cstdio> // printf
#include <cassert> // assert
#include <cmath> // sqrt, pow, exp, fabs, sqrt
#include <fstream> // ifstream, ofstream
#include <algorithm> // min, max
#include <iomanip> // setprecision

#include "atom_core.hxx"

#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // create_default_radial_grid, dot_product<real_t>
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
#include "display_units.h" // eV, _eV
#include "radial_potential.hxx" // Hartree_potential
#include "exchange_correlation.hxx" // lda_PZ81_kernel
#include "inline_math.hxx" // pow2, set, scale, product, add_product
#include "radial_eigensolver.hxx" // shooting_method
#include "constants.hxx" // constants::pi

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

  void rad_pot(double rV[], radial_grid_t const &g, double const rho4pi[], double const Z, double *energies) {
      double Eexc = 0, Evxc = 0, EHtr = 0;
      double const ECou = -Z*radial_potential::Hartree_potential(rV, g, rho4pi); // set rV to the Hartree potential
      double const fpi = .25/constants::pi;
      for(int ir = 0; ir < g.n; ++ir) {
          EHtr += rho4pi[ir]*rV[ir]*g.rdr[ir];
          double Vxc;
          double const Exc = exchange_correlation::lda_PZ81_kernel(fpi*rho4pi[ir], Vxc);
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

    
  double initial_density(double r2rho[], radial_grid_t const &g, double const Z, double const charged) {
    auto const alpha = 0.3058*std::pow(Z, 1./3.);
    auto const beta = std::sqrt(108./constants::pi)*std::max(0., Z - 2)*alpha;
    if (4 + 3.2*charged < 0) return 1; // error
    auto const gamma = std::sqrt(4 + 3.2*charged);
    auto g2 = gamma*gamma*gamma;
    if (Z < 2) g2 *= 0.5*Z;
    double q = 0;
    for(auto ir = 0; ir < g.n; ++ir) {
        auto const r = g.r[ir];
        auto const x = alpha*r;
        auto const exa = (x < 50)? std::exp(-3*x) : 0;
        auto const xb = gamma*r;
        auto const exb = (xb < 150)? std::exp(-xb) : 0;
        r2rho[ir] = beta*std::sqrt(x)*exa + g2*exb*r*r;
//      printf("%g %g\n", r, r2rho[ir]);
        q += r2rho[ir] * g.dr[ir]; // integrate total charge
    } // ir
    return q; // charge
  } // initial_density
  
  
  
  status_t read_Zeff_from_file(double Zeff[], radial_grid_t const &g, float const Z, 
          char const name[], float const factor, int const echo) {
      status_t stat = 0;
      char filename[99]; sprintf(filename, "%s.%03d", name, (int)std::ceil(Z));
      if (echo > 3) printf("# %s  Z=%g  try to read file %s\n",  __func__, Z, filename);
      std::ifstream infile(filename);
      int ngr = 0, ngu = 0; // number of radial grid points read and used
      if (infile.is_open()) {
          double r_min = 9e9, r_max = - 9e9;
          int ir = 0;
          double r, Ze, r_prev=0, Ze_prev = Z*factor;
          while (infile >> r >> Ze) {
              ++ngr;
              if (r >= 0) {
                  r_min = std::min(r, r_min);
                  r_max = std::max(r, r_max);
                  if (r <= g.rmax + 1e-6) {
                      full_debug(printf("# %s  r=%g Zeff=%g\n",  __func__, r, Ze));
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
          if (echo > 3) printf("# %s  Z=%g  use %d of %d values from file %s, interpolate to %d values\n",
                                  __func__, Z, ngu, ngr, filename, ir);
      } else {
          if (echo > 1) printf("# %s  Z=%g  failed to open file %s\n",  __func__, Z, filename);
          stat = -1; // failure
      } // is_open     
      return stat;
  } // read_Zeff_from_file

  status_t store_Zeff_to_file(double const Zeff[], radial_grid_t const &g, float const Z, 
      char const name[]="pot/Zeff", float const factor=1, int const echo=9) {
      char filename[99]; sprintf(filename, "%s.%03d", name, (int)std::ceil(Z));
      if (echo > 3) printf("# %s  Z=%g  try to write file %s\n",  __func__, Z, filename);
      std::ofstream outfile(filename);
      if (outfile.is_open()) {
          outfile << std::setprecision(15);
          for(int ir = 0; ir < g.n; ++ir) {
              outfile << g.r[ir] << " " << Zeff[ir]*factor << "\n";
          } // write to file
          return 0; // success
      } else {
          if (echo > 1) printf("# %s  Z=%g  failed to open file %s for writing\n",  __func__, Z, filename);
      } // is_open
      return -1; // failure
  } // store_Zeff_to_file
  


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
                  orb[i].occ = std::min(std::max(0.f, Z - iZ), max_occ*1.f);
                  orb[i].E = guess_energy(Z, enn);
                  if (orb[i].occ > 0) {
                      if (echo > 4) {
                          printf("# %s  i=%d %d%c f= %g  guess E= %g %s\n", __func__, i, enn, ellchar(ell), orb[i].occ, orb[i].E*eV, _eV);
                      } // echo
                      imax = std::max(i, imax);
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
      set(rV_old, g.n, -1.*Z); // Hydrogen-like potential
      set(rV_new, g.n, -1.*Z); // Hydrogen-like potential
      double mix = 0; // current mixing coefficient for the potential 
      double const alpha = 0.33; // limit case potential mixing coefficient
      
      double energies[NumEnergyContributions];
      double eigenvalue_sum = 0;
      double previous_energy = 0;

      enum { Task_Solve, Task_ChkRho, Task_GenPot, Task_MixPot, Task_Energy } next_task = Task_Solve;

      int icyc = 0;
      { // start scope
          bool const loading_failed = (0 != read_Zeff_from_file(rV_old, g, Z, "pot/Zeff", -1));
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
                
                  set(r2rho4pi, g.n, 0.0); // init accumulator density
                  eigenvalue_sum = 0; // init energy accumulator

                  for(int i = 0; i <= imax; ++i) {
                      if (orb[i].occ > 0) {
                          previous_eigenvalues[i] = orb[i].E; // copy
                          auto const stat = radial_eigensolver::shooting_method(sra, g, rV_old, orb[i].enn, orb[i].ell, orb[i].E, nullptr, r2rho);
                          if (stat) {
                              printf("# %s  Z=%g  failed solving for %d%c, status = %d\n",  __func__, Z, orb[i].enn, ellchar(orb[i].ell), stat);
                              return stat;
                          } // failed
                          if (echo > 6) {
                              printf("# %s  Z=%g  %d%c E=%15.6f %s\n",  __func__, Z, orb[i].enn, ellchar(orb[i].ell), orb[i].E*eV, _eV);
                          } // echo
                          // add orbital density
                          double const q = radial_grid::dot_product(g.n, r2rho, g.dr);
                          assert(q > 0);
                          double const f = orb[i].occ/q;
                          add_product(r2rho4pi, g.n, r2rho, f);
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
                  
                  product(rho4pi, g.n, r2rho4pi, g.rinv, g.rinv);

                  double const q  = radial_grid::dot_product(g.n, rho4pi, g.r2dr); // can be merged with the loop above
                  double const qq = radial_grid::dot_product(g.n, r2rho4pi, g.dr); // can be merged with the loop above
                  if (std::abs(q - Z) > .1 && echo > 0) {
                      printf("# %s  Z=%g  icyc=%d  Warning: rho has %.6f (or %.6f) electrons\n",  __func__, Z, icyc, q, qq);
                  } // not charge neutral
                  
                  next_task = Task_GenPot;
              } break; // Task_ChkRho
              ///////////////////////////////////////////////////////
              case Task_GenPot: {
                  full_debug(printf("# Task_Generate_Pot\n"));
                  
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
                  res = std::abs(energies[which] - previous_energy); 
                  previous_energy = energies[which]; // store for the next iteration
                  if (echo > 3) {
                      int const display_every = 1 << std::max(0, 2*(6 - echo)); // echo=4:every 16th, 5:every 4th, 6:every cycle
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
                      scale(rV_old, g.n, 1. - mix);
                      add_product(rV_old, g.n, rV_new, mix);
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
              
              store_Zeff_to_file(rV_old, g, Z, "pot/Zeff", -1);

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

  template<typename real_t>
  int RamerDouglasPeucker(char active[] // used as boolean, length on[end]
        , real_t const x_list[], real_t const y_list[], int const end
        , int const begin=0, float const epsilon=1e-6) {
    // Find the point with the maximum distance
    double d2max = 0;
    int index = -1;
    // construct a straight line through the 1st and last point
    int const p0 = begin, pl = end - 1;
    assert(active[p0]); assert(active[pl]); // both points still have to be active
    double const x0 = x_list[p0], y0 = y_list[p0];
    double const x_dist = x_list[pl] - x0;
    double const y_dist = y_list[pl] - y0;
    double const det = x_dist*x_dist + y_dist*y_dist;
    double const det_inv = 1./det;
    
    for(int i = p0 + 1; i < pl; ++i) { // loop over all points in between
        if (active[i]) {
            double const xi = x_list[i], yi = y_list[i];
            /*
            // solve
            //    p0x + s0*x_dist == pix + si*y_dist == 0
            //    p0y + s0*y_dist == piy - si*x_dist == 0
            // for s, t:
            //    / x_dist   -y_dist \   / s0 \     / pix - p0x \
            //    |                  | * |    |  == |           |
            //    \ y_dist    x_dist /   \ si /     \ piy - p0y /
            //
            // solution: adjoint matrix divided by determinant
            //    / s0 \     / x_dist    y_dist \   / pix - p0x \     1
            //    |    |  == |                  | * |           |  * ----
            //    \ si /     \ -y_dist   x_dist /   \ piy - p0y /    det
            */
//             double const s0 = (x_dist*(xi - x0) + y_dist*(yi - y0))*det_inv;
//             double const xf0 = x0 + s0*x_dist;
//             double const yf0 = y0 + s0*y_dist;
            
            double const si = (x_dist*(yi - y0) - y_dist*(xi - x0))*det_inv;
//             double const xf = xi + si*y_dist;
//             double const yf = yi - si*x_dist;
//             assert( pow2(xf - xf0) + pow2(yf - yf0) < 1e-12);
          
//             double const d2 = pow2(yf - yi) + pow2(xf - xi); // perpendicularDistance^2
            double const d2 = si*si*det;
    //         printf("# DouglasPeucker distance^2 = %g at point #%d\n", d2, i);
            if (d2 > d2max) {
                index = i;
                d2max = d2;
            } // find the maximum
        } // active?
    } // i
    
    // If max distance is greater than epsilon, recursively simplify
    if (d2max > epsilon*epsilon) {
        assert(index > p0); assert(index < pl);
//      printf("#      DouglasPeucker case N, distance^2 = %g at point #%d\n", d2max, index);
        // Recursive call
//      printf("# call DouglasPeucker on interval [%d, %d] and  [%d, %d] later \n", begin, index, index, end - 1);
        int const n0 = RamerDouglasPeucker(active, x_list, y_list, index + 1, begin, epsilon);
        int const n1 = RamerDouglasPeucker(active, x_list, y_list, end,       index, epsilon);
        return n0 + n1 - 1;
    } else {
        int n_off = 0;
        for(int i = p0 + 1; i < pl; ++i) {
            n_off += (0 != active[i]);
            active[i] = 0; // switch off every point in between since the curve is sufficiently smooth there
        } // i
//      if (n_off > 0) printf("# DouglasPeucker eleminate %d points, largest distance^2 is %g\n", n_off, d2max);
        return 2;
    } // if

  } // RamerDouglasPeucker
  
  status_t simplify_Zeff_file(int const iZ, float const epsilon=1e-6, int const echo=9) {
      float const Z = iZ;
      auto const g = radial_grid::create_default_radial_grid(Z);
      auto const y = new double[g->n];
      auto stat = read_Zeff_from_file(y, *g, Z, "full_pot/Zeff", 1., echo);
      
      auto const active = new char[g->n];
      set(active, g->n, (char)1); // init all entries as active
      int const new_n = RamerDouglasPeucker(active, g->r, y, g->n, 0, epsilon);
      if (echo > 0) printf("# RamerDouglasPeucker reduced %d points to %d points\n", g->n, new_n);
      auto const new_y = new double[new_n];
      radial_grid_t new_g; new_g.n = new_n; new_g.r = new double[new_n]; // this radial grid is incomplete and only provides r[]
      {
          int i = 0;
          for(int ir = 0; ir < g->n; ++ir) {
              if (active[ir]) {
                  new_g.r[i] = g->r[ir];
                  new_y[i] = y[ir];
                  ++i;
              } // active
          } // ir
          assert(i == new_n); // after running RamerDouglasPeucker, there must be exactly new_n active entries left
          delete[] active;
          delete[] y;
      }
      stat += store_Zeff_to_file(new_y, new_g, Z, "pot/Zeff", 1., echo);
      delete[] new_g.r;
      delete[] new_y;
      return stat;
  } // simplify_Zeff_file
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_initial_density(radial_grid_t const &g, int const echo=0) {
    printf("\n# %s:%d  %s \n\n", __FILE__, __LINE__, __func__);
    auto rho = new double[g.n];
    double diff = 0;
    for(float zz = 0; zz < 128; zz += 1) {
        initial_density(rho, g, zz);
        double const q = radial_grid::dot_product(g.n, g.r2dr, rho);
        diff = std::max(diff, std::abs(zz - q));
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

  status_t test_nl_index(int const echo=2) {
      int inl = 0;
      for(int enn = 1; enn < 9; ++enn) {
          for(int ell = 0; ell < enn; ++ell) {
              int const nl = nl_index(enn, ell);
              if ((inl - nl) || (echo > 8)) printf("# %s: %s n=%d l=%d nl_index %d %d\n", __FILE__, __func__, enn, ell, inl, nl);
              assert(inl == nl);
              ++inl;
          } // ell
      } // enn
      return 0;
  } // test_nl_index
  
  status_t all_tests() {
    auto status = 0;
//  status += test_initial_density(*radial_grid::create_exponential_radial_grid(512));
    status += test_nl_index();
//     for(int Z = 120; Z >= 0; --Z) { // test all atoms, backwards
//         printf("\n\n# Z = %d\n\n", Z);
//         status += simplify_Zeff_file(Z, 1e-8);
    
//     for(int Z = 119; Z > 0; --Z) { // test all atoms, backwards
    // for(int Z = 0; Z < 120; ++Z) { // test all atoms, backwards
        // printf("\n\n# Z = %d\n\n", Z);
        // status += test_core_solver(*radial_grid::create_default_radial_grid(Z), Z);
//     { int const Z = 1;  // 1:hydrogen
//     { int const Z = 5;  // 5:boron
//        { int const Z = 6;  // 6:carbon
    { int const Z = 29; // 29:copper
//     { int const Z = 70; // 70:ytterbium
//     { int const Z = 79; // 79:gold
//     { int const Z = 120; // very heavy non-existing element
        status += test_core_solver(*radial_grid::create_default_radial_grid(Z), Z);
//         status += simplify_Zeff_file(Z);
    } // Z
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace atom_core
