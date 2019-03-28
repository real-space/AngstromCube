#include <vector> // std::vector
#include <cstdio> // printf
#include <cstdlib> // abs
#include <cmath> // sqrt

#include "radial_eigensolver.hxx"

#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // create_exponential_radial_grid
#include "inline_math.hxx" // sgn, pow2
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
#include "output_units.h" // eV, _eV, Ang, _Ang
#include "radial_integrator.hxx" // shoot, integrate_inwards, integrate_outwards

// #define FULL_DEBUG
// #define DEBUG

#ifdef FULL_DEBUG
    #define full_debug(print) print 
    #define DEBUG
#else  
    #define full_debug(print)
#endif

#ifdef  DEBUG
    #include "debug_output.hxx" // dump_to_file
#endif


namespace radial_eigensolver {
  // solves the radial eigenvalue problem of the spherical potential with the shooting method
  using namespace radial_grid;
  using namespace radial_integrator;

  status_t shooting_method(
      int const sra, // 1:scalar relativistic approximation, 0:Schroedinger equation
      radial_grid_t const &g, // radial grid descriptor
      double const rV[], // radial potential r*V_Hxc(r) - e^2*Z
      enn_QN_t const enn, // energy quantum number
      ell_QN_t const ell, // angular momentum quantum number
      double &E, // energy (eigen-)value in Hartree
      double* rf, // [optional, out] radial wave function*r
      double* r2rho, // [optional, out] density of that wave function*r^2
      int const maxiter, // [optional] maximum number of iterations
      float const threshold) // [optional] threshold for eigenvalue convergence
  {
#ifdef DEBUG
      if (enn < 1) exit(__LINE__); // "rINT shooting_method: ENN < 1 unphysical."
      if (ell >= enn) exit(__LINE__); // "rINT shooting_method: ELL must be < ENN"
#endif

      auto const max_dE = 1.25e-4 * 0.5*(pow2(-rV[0]/enn) + .1); // for jumps, if the number of nodes is wrong, -rV[0]==Z
      int const nno = enn - 1 - ell; // number of nodes requested

      int nn = 0;
      double kink = shoot(sra, g, rV, ell, E, nn);

      int constexpr MaxIter_node_count = 999;
      int iit = 1;
      while ((nno != nn) && iit < MaxIter_node_count) { // while number of nodes incorrect
          ++iit;
          E += (nno - nn) * max_dE;
          kink = shoot(sra, g, rV, ell, E, nn);
#ifdef  FULL_DEBUG
          printf("# %s: find-correct-node for n=%d l=%d E= %.9f %s, %d nodes expected, %d nodes found\n", __func__, enn, ell, E*eV, _eV, nno, nn);
#endif
      } // while

#ifdef DEBUG
      printf("# %s: needed %d iterations to find the correct number of nodes for n=%d l=%d E= %.9f %s\n", __func__, iit, enn, ell, E*eV, _eV);
#endif

      
      double ene[2], knk[2], mdE[2] = {max_dE, max_dE};
      int nnn[2];
#ifdef  DEBUG
      int itn[2] = {0, 0};
#endif
      for(auto ib = 0; ib < 2; ++ib) {
          ene[ib] = E + (2*ib - 1) * max_dE; // initialize with an energy which produces the correct number of nodes
          knk[ib] = shoot(sra, g, rV, ell, ene[ib], nnn[ib]);
          int constexpr MaxIter_kink_sign = 999;
          iit = 0;
          while (((2*ib - 1)*knk[ib] > 0.) && (iit < MaxIter_kink_sign)) { // while sign of kink is incorrect
              ++iit;
              ene[ib] += (2*ib - 1) * mdE[ib];
              mdE[ib] *= 1.1; // growing exponentially
              knk[ib] = shoot(sra, g, rV, ell, ene[ib], nnn[ib]);
#ifdef  FULL_DEBUG
              printf("# %s: get-correct-kink-sign for (%s) n=%d l=%d E=%g %s, kink= %g, %d nodes\n", __func__, (ib)? "upper" : "lower", enn, ell, ene[ib]*eV, _eV, knk[ib], nnn[ib]);
#endif
#ifdef  FULL_DEBUG
              if (nno != nnn[ib])
                  printf("# %s: warning for n=%d l=%d E= %.9f %s, %d nodes expected, %d nodes found\n", __func__, enn, ell, E*eV, _eV, nno, nnn[ib]);
#endif
          } // while
#ifdef  DEBUG
          itn[ib] = iit;
#endif            
      } // ib
#ifdef  DEBUG
      printf("# %s: needed %d and %d iterations for the start values n=%d l=%d\n", __func__, itn[0], itn[1], enn, ell);
      printf("# %s: start interval [%g, %g] %s, kinks = [%.6f, %.6f] for n=%d l=%d\n", __func__, ene[0]*eV, ene[1]*eV, _eV, knk[0], knk[1], enn, ell);
#endif

// #ifdef DEBUG
//       printf("# %s: set start energy interval for n=%d l=%d to [%.9f, %.9f] %s\n", __func__, enn, ell, ene[0]*eV, ene[1]*eV, _eV);
// #endif

      double E_prev = 0; // previous energy
      iit = 0;
      bool run = true, converged=false;
      while (run) {
          ++iit;

          // new energy in the middle of ene[0] and ene[1]
          E = 0.5*(ene[0] + ene[1]);
          
#ifdef  FULL_DEBUG
          printf("# %s: converge-energy for n=%d l=%d, try E= %.9f %s\n", __func__, enn, ell, E*eV, _eV);
#endif
          kink = shoot(sra, g, rV, ell, E, nn);

          // bisection: reset either the left or the right boundary
          int const ib = (kink < 0)? 1 : 0;
          ene[ib] = E;
#ifdef  FULL_DEBUG
          printf("# %s: found kink = %g, adjust %s limit: [%.9f, %.9f] %s\n", __func__, kink, (ib)? "upper" : "lower", ene[0]*eV, ene[1]*eV, _eV);
#endif
          auto const res = fabs(E - E_prev); // calculate energy change
          E_prev = E;

          converged = (res < threshold);
          run = (!converged) && (iit < maxiter);
      } // while(run)

#ifdef  FULL_DEBUG
      printf("# %s: needed %d iterations to converge to E= %.9f %s\n", __func__, iit, E*eV, _eV);
#endif

      // call the shooting once more, potentially with export of rf and r2rho
      kink = shoot(sra, g, rV, ell, E, nn, rf, r2rho);

      if (converged) {
          return 0; // success
      } else if (iit >= maxiter) {
          return maxiter;
      } else if (nn != nno) {
          return -33;
      } else {
          return -66; // unknown error
      }
  } // shooting_method

  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_hydrogen_like_potential(
    radial_grid_t const g, // radial grid descriptor
    float const Z) { // number of protons in the nucleus
      
    status_t status = 0;
    auto const ellchar = "spdfghijkl";
    auto rV = std::vector<double>(g.n, -Z); // fill all potential values with r*V(r) == -Z
    auto const rf = new double[g.n];
    for(auto sra = 1; sra <= 1; ++sra) { // 0:non-relativistic, 1:scalar-relativistic, 2:scalar-rel-with-linearized-sqrt
      printf("\n\n# %s %s SRA approximation level = %d\n", __FILE__, __func__, sra);
      for(auto enn = 1; enn <= 9; ++enn) {
        for(auto ell = 0; ell < enn; ++ell) {
          double E = -.5*pow2(Z/enn); // guess energy for hydrogen like atoms
          status += abs(shooting_method(sra, g, rV.data(), enn, ell, E, rf));
          printf("%2d%c energy for Z = %.3f found at E = %.12f %s\n", enn, ellchar[ell], Z, E*eV, _eV);
#ifdef  DEBUG
          { char filename[32]; sprintf(filename, "Z%d%c_radial_wave_function.dat", enn, ellchar[ell]);
            dump_to_file(filename, g.n, rf, g.r);
          }
#endif
        } // ell
      } // enn
    } // sra
    return status;
  } // test_hydrogen_like_potential
  
  
  status_t all_tests() {
    auto status = 0;
    status += test_hydrogen_like_potential(*create_exponential_radial_grid(2610), 100);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace radial_eigensolver
