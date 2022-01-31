#include <vector> // std::vector
#include <cstdio> // std::printf, std::snprintf
#include <cstdlib> // std::abs

#include "radial_eigensolver.hxx"

#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // ::create_radial_grid, ::destroy_radial_grid
#include "inline_math.hxx" // sgn, pow2
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
#include "display_units.h" // eV, _eV, Ang, _Ang
#include "radial_integrator.hxx" // ::shoot
#include "recorded_warnings.hxx" // warn

#ifdef FULL_DEBUG
    #define DEBUG
#endif // FULL_DEBUG

#ifdef  DEBUG
    #include "debug_output.hxx" // dump_to_file
#endif // DEBUG


namespace radial_eigensolver {
  // solves the radial eigenvalue problem of the spherical potential with the shooting method
  
  auto const ellchar = "spdfghijkl";

  status_t shooting_method(
        int const sra // 1:scalar relativistic approximation, 0:Schroedinger equation
      , radial_grid_t const & g // radial grid descriptor
      , double const rV[] // radial potential r*V_Hxc(r) - e^2*Z
      , enn_QN_t const enn // energy quantum number
      , ell_QN_t const ell // angular momentum quantum number
      , double & E // energy (eigen-)value in Hartree
      , double* rf // [optional, out] radial wave function*r
      , double* r2rho // [optional, out] density of that wave function*r^2
      , int const maxiter // [optional] maximum number of iterations
      , float const threshold // [optional] threshold for eigenvalue convergence
  ) {
#ifdef DEBUG
      if (enn < 1) exit(__LINE__); // "rINT shooting_method: ENN < 1 unphysical."
      if (ell >= enn) exit(__LINE__); // "rINT shooting_method: ELL must be < ENN"
#endif // DEBUG

      double constexpr lower_energy_stop = -1e5;
      double max_dE = 1.25e-4 * 0.5*(pow2(-rV[0]/enn) + .1); // for jumps, if the number of nodes is wrong, -rV[0]==Z
      int const nno = enn - 1 - ell; // number of nodes requested

      int nn{0}; // number of nodes
      double kink = radial_integrator::shoot(sra, g, rV, ell, E, nn);
      
      // first make sure that we are on the right branch with the correct node count
      int constexpr MaxIter_node_count = 999;
      int iit{1}; // start from 1 since we already envoked shoot once above
      while ((nno != nn) && iit < MaxIter_node_count) { // while number of nodes incorrect
          ++iit;
          E += (nno - nn) * max_dE;
          max_dE *= 1.0625; // growing exponentially
          kink = radial_integrator::shoot(sra, g, rV, ell, E, nn);
#ifdef  FULL_DEBUG
          std::printf("# %s: find-correct-node for n=%d l=%d E= %.9f %s, %d nodes expected, %d nodes found\n", __func__, enn, ell, E*eV, _eV, nno, nn);
#endif // FULL_DEBUG
          if (E < lower_energy_stop) {
              warn("%d%c-energy E= %g Ha fell below %g Ha while searching the branch", enn,ellchar[ell], E, lower_energy_stop);
              return -99; // something went wrong
          }
      } // while

#ifdef DEBUG
      std::printf("\n# %s: needed %d iterations to find the correct number of nodes for n=%d l=%d E= %.9f %s\n", __func__, iit, enn, ell, E*eV, _eV);
#endif // DEBUG

      // next we have to ensure that we find two start energies for which the kinks have different signs
      // and which both lead to the correct node count
      double ene[2], knk[2], mdE[2] = {.01, .01};
      double const inc[2] = {1.125, 1.01};
      int nnn[2];
#ifdef  DEBUG
      int itn[2] = {0, 0};
#endif // DEBUG
      for (auto ib = 0; ib < 2; ++ib) { // ib in {0, 1} == {lower, upper}
          auto const sgn_ib = (2*ib - 1); // sgn in {-1, +1}
          ene[ib] = E; // initialize with an energy which produces the correct number of nodes ...
          knk[ib] = kink; // ... and the corresponding kink value from the node-count search
          int constexpr MaxIter_kink_sign = 999;
          iit = 0;
          while ((sgn_ib*knk[ib] > 0) && (iit < MaxIter_kink_sign)) { // while sign of kink is incorrect
              ++iit;
              ene[ib] += sgn_ib * mdE[ib];
              mdE[ib] *= inc[ib]; // growing exponentially for the next iteration
              knk[ib] = radial_integrator::shoot(sra, g, rV, ell, ene[ib], nnn[ib]);
#ifdef  FULL_DEBUG
              std::printf("# %s: get-correct-kink-sign for (%s) n=%d l=%d E=%g %s, kink= %g, %d nodes\n", __func__, ib?"upper":"lower", enn, ell, ene[ib]*eV, _eV, knk[ib], nnn[ib]);
#endif // FULL_DEBUG
#ifdef  DEBUG
              if (nno != nnn[ib]) std::printf("# %s: Warning for n=%d l=%d E= %.15f %s, %d nodes expected, %d nodes found\n", 
                                            __func__, enn, ell, ene[ib]*eV, _eV, nno, nnn[ib]);
#endif // DEBUG
              if (E < lower_energy_stop) {
                  warn("%d%c-energy E= %g Ha fell below %g Ha while searching the interval", enn,ellchar[ell], E, lower_energy_stop);
                  return -98; // something went wrong
              }
          } // while
#ifdef  DEBUG
          itn[ib] = iit;
#endif // DEBUG
      } // ib
#ifdef  DEBUG
      std::printf("# %s: needed %d and %d iterations for the start values n=%d l=%d\n", __func__, itn[0], itn[1], enn, ell);
      std::printf("# %s: start interval [%g, %g] %s, kinks = [%.6f, %.6f] for n=%d l=%d\n", __func__, ene[0]*eV, ene[1]*eV, _eV, knk[0], knk[1], enn, ell);
#endif // DEBUG

// #ifdef DEBUG
//       std::printf("# %s: set start energy interval for n=%d l=%d to [%.9f, %.9f] %s\n", __func__, enn, ell, ene[0]*eV, ene[1]*eV, _eV);
// #endif // DEBUG

      double E_prev{0}; // previous energy
      iit = 0;
      bool run{true}, converged{false};
      while (run) {
          ++iit;

          // new energy in the middle of ene[0] and ene[1]
          E = 0.5*(ene[0] + ene[1]);
          
#ifdef  FULL_DEBUG
          std::printf("# %s: converge-energy for n=%d l=%d, try E= %.9f %s\n", __func__, enn, ell, E*eV, _eV);
#endif // FULL_DEBUG
          kink = radial_integrator::shoot(sra, g, rV, ell, E, nn);

          // bisection: reset either the left or the right boundary
          int const ib = (kink < 0)? 1 : 0;
          ene[ib] = E;
#ifdef  FULL_DEBUG
          std::printf("# %s: found kink = %g, adjust %s limit: [%.9f, %.9f] %s\n", __func__, kink, (ib)? "upper" : "lower", ene[0]*eV, ene[1]*eV, _eV);
#endif // FULL_DEBUG
          auto const res = std::abs(E - E_prev); // calculate energy change
          E_prev = E;

          converged = (res < threshold);
          run = (!converged) && (iit < maxiter);
      } // while(run)

#ifdef  FULL_DEBUG
      std::printf("# %s: needed %d iterations to converge to E= %.9f %s\n", __func__, iit, E*eV, _eV);
#endif // FULL_DEBUG

      // call the shooting once more, potentially with export of rf and r2rho
      kink = radial_integrator::shoot(sra, g, rV, ell, E, nn, rf, r2rho);

      if (converged) {
          if (nn != nno) warn("%d%c-eigenstate at E= %g %s has %d nodes but expected %d", enn,ellchar[ell], E*eV, _eV, nn, nno);
          return nn - nno; // success if the number of nodes is correct
      } else if (iit >= maxiter) {
          warn("Number of iterations for %d%c-eigenstate exceeded max of %d", enn,ellchar[ell], maxiter);
          return maxiter;
      } else if (nn != nno) {
          warn("%d%c-eigenstate did not converge", enn,ellchar[ell]);
          return -33;
      } else {
          warn("Unknown error in %d%c-eigenstate E= %g %s and %d nodes after %d iterations", enn,ellchar[ell], E*eV, _eV, nn, iit);
          return -66; // unknown error
      }
  } // shooting_method

  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_hydrogen_like_potential(
        int const echo=5 // log-level
      , double const Z=100 // atomic number, number of protons in the nucleus
  ) {
      status_t stat(0);
      auto & g = *radial_grid::create_radial_grid(2610);
      std::vector<double> const rV(g.n, -Z); // fill all potential values with r*V(r) == -Z
      std::vector<double> rf(g.n); // radial wave function
      char const SRA_name[][21] = {"non-relativistic", "scalar-relativistic", "linearized-sqrt"};
      double maxdev{0}; // deviation of the non-relativistic final energies from the reference energy
      for (auto sra = 0; sra <= 2; ++sra) { // 0:, 1:, 2:
          if (echo > 0) std::printf("\n\n# %s %s (Z= %g) %s\n", __FILE__, __func__, Z, SRA_name[sra]);
          for (auto enn = 1; enn <= 9; ++enn) {
              auto const E_ref = -0.5*pow2(Z/double(enn)); // for sra == 0 we can compute a reference energy for a hydrogen-like potential
              if (echo > 2 && 0 == sra) std::printf("# %2d -energy %.12f %s reference\n", enn, E_ref*eV, _eV);
              for (auto ell = 0; ell < enn; ++ell) {
                  double E{E_ref}; // guess a start energy for hydrogen like atoms
                  stat += std::abs(int(shooting_method(sra, g, rV.data(), enn, ell, E, rf.data())));
                  if (echo > 1) {
                      if (0 == sra) {
                          std::printf("# %2d%c-energy %.12f  dev %.1e %s\n", enn, ellchar[ell], E*eV, (E - E_ref)*eV, _eV);
                          maxdev = std::max(maxdev, std::abs(E - E_ref));
                      } else {
                          std::printf("# %2d%c-energy %.12f %s\n", enn, ellchar[ell], E*eV, _eV);
                      }
                  } // echo
#ifdef  DEBUG
                  char filename[32]; std::snprintf(filename, 31, "Z%d%c_radial_wave_function.dat", enn, ellchar[ell]);
                  dump_to_file(filename, g.n, rf.data(), g.r);
#endif // DEBUG
              } // ell
          } // enn
      } // sra
      if (echo > 1) std::printf("# %s: largest energy deviation for %s solver is %.1e %s\n", __func__, SRA_name[0], maxdev*eV, _eV);
      radial_grid::destroy_radial_grid(&g);
      return stat + (maxdev > 1e-8*pow2(Z));
  } // test_hydrogen_like_potential


  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_hydrogen_like_potential(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace radial_eigensolver
