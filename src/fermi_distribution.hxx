#pragma once

#include <cstdio> // printf
#include <cassert> // assert
#include <cmath> // std::exp
#include <vector> // std::vector

#include "inline_math.hxx" // set
#include "simple_math.hxx" // ::random
#include "display_units.h" // eV, _eV
#include "recorded_warnings.hxx" // warn

#include "status.hxx" // status_t

namespace fermi_distribution {
  
  template <typename real_t>
  inline real_t FermiDirac(real_t const x, real_t *derivative=nullptr) { 
      if (std::abs(x) > 36) {
          if (derivative) *derivative = 0;
          return real_t(x < 0);
      } // |x| large
      real_t const d = 1 + std::exp(x); // denominator
      real_t const f = 1/d;
      if (derivative) *derivative = (1 - d)*f*f;
      return f;
  } // FermiDirac

  template <typename real_t>
  real_t FermiDirac_broadening(real_t const x) {
      real_t der;
      FermiDirac(x, &der);
      return -der;
  } // FermiDirac_broadening

  template <typename real_t>
  double count_electrons(int const n, real_t const energies[], 
      double const eF, double const kTinv, 
      double const weights[]=nullptr, double *derivative=nullptr, double occupations[]=nullptr) {
      double ne{0}, dnde{0};
      assert(nullptr == weights); // not implemented!
      for(int i = 0; i < n; ++i) {
          double dfde;
          double const occ = FermiDirac((energies[i] - eF)*kTinv, &dfde);
          double const w8 = weights ? weights[i] : 1;
          ne   += occ  *w8;
          dnde -= dfde *w8;
          if (occupations) occupations[i] = occ;
      } // i
      if (derivative) *derivative = dnde*kTinv;
      return ne;
  } // count_electrons

  template <typename real_t>
  status_t Fermi_level(
        double & eF
      , double occupations[]
      , real_t const energies[] // energies
      , int const nb // number of bands
      , double const kT
      , double const number_of_electrons
      , int const spin_factor=2 // 2:spin_paired, 1:spin_resolved
      , int const echo=9
  ) {
      double const n_electrons = number_of_electrons/spin_factor;
      // ToDo: count the states with weights, check if there are enough states to host n_electrons
      std::vector<real_t> ene(nb); // energies
      set(ene.data(), nb, energies); // copy
      std::sort(ene.begin(), ene.end()); // sort to ascending order
      auto const e_min = ene[0], e_max = ene[nb - 1];
      if (echo > 0) printf("# %s E_min= %g E_max= %g %s\n", __func__, e_min*eV, e_max*eV, _eV);
      
      { // scope: find T=0 Fermi level:
          int ieF{0};
          while (ieF < n_electrons) ++ieF; // ignores weights
          eF = (ieF < nb) ? ((ieF > 0) ? 0.5*(ene[ieF] + ene[ieF - 1]) : e_min) : e_max;
          if (echo > 0) printf("# %s T=0 Fermi level at %g %s\n", __func__, eF*eV, _eV);
      } // scope

      double const kTinv = 1./std::max(1e-9, kT);
      
      double e[2] = {eF, eF}, ne[2];
      double delta_e = std::max(kT, 1e-3);
      for(int i01 = 0; i01 <= 1; ++i01) {
          double const sgn = 2*i01 - 1; // {-1, 1}
          int change{0};
          // prepare start value for bisection
          do {
              e[i01] += sgn*change*delta_e; change = 1; delta_e *= 1.125;
              ne[i01] = count_electrons(nb, energies, e[i01], kTinv);
              if (echo > 5) { printf("# %s correct %ser start energy to %g %s --> %g electrons\n", 
                            __func__, i01?"upp":"low", e[i01]*eV, _eV, spin_factor*ne[i01]); fflush(stdout); }
          } while(sgn*ne[i01] < sgn*n_electrons);
      } // i01

      // start bisection
      double res{1};
      int const maxiter = 99;
      int iter{0}; // init with max. number of iterations, hard limit
      while(res > 1e-9 && iter < maxiter) {
          ++iter;
          auto const em = 0.5*(e[0] + e[1]);
          auto const nem = count_electrons(nb, energies, em, kTinv);
          if (echo > 7) { printf("# %s with energy %g %s --> %g electrons\n", __func__, em*eV, _eV, spin_factor*nem); fflush(stdout); }
          int const i01 = (nem > n_electrons);
          e[i01] = em;
          ne[i01] = nem;
          res = std::abs(nem - n_electrons);
      }
      eF = 0.5*(e[0] + e[1]);
      double DoS_at_eF;
      auto const nem = count_electrons(nb, energies, eF, kTinv, nullptr, &DoS_at_eF, occupations);
      if (echo > 6) printf("# %s with energy %g %s --> %g electrons\n", __func__, eF*eV, _eV, spin_factor*nem); 
      if (res > 1e-9) {
          warn("Fermi level converged only to +/- %.1e electrons in %d iterations", res*spin_factor, iter); 
      } else {
          if (echo > 3) printf("# %s at %g %s has a DOS of %g states\n", __func__, eF*eV, _eV, DoS_at_eF*spin_factor);
      }
      return 0;
  } // Fermi_level

  inline
  status_t density_of_states(int const n, double const energies[], double const kT, double const eF=0, double const de=3e-4, int const echo=9) {
      if (n < 1) return 0;
      double const kTinv = 1./std::max(1e-9, kT);
      double e_min{9e307}, e_max{-e_min}; //, e_stat[4] = {0,0,0,0};
      for(int i = 0; i < n; ++i) {
          auto const ei = energies[i] - eF;
          e_min = std::min(e_min, ei);
          e_max = std::max(e_max, ei);
//           e_stat[0] += 1;
//           e_stat[1] += ei;
//           e_stat[2] += pow2(ei);
//           e_stat[3] += pow3(ei);
      } // i
      if (echo > 4) printf("# %s E_min= %g E_max= %g %s\n", __func__, e_min*eV, e_max*eV, _eV);

      double const per_eV = 1./eV;
      if (echo > 8) printf("\n## energy(%s), DoS, dDoS/dE_Fermi\n", _eV);
      double integral{0};
      for(int ie = (e_min - 18*kT)/de; ie <= (e_max + 18*kT)/de; ++ie) {
          double const e = ie*de;
          double dos{0}, ddos{0};
          for(int i = 0; i < n; ++i) {
              auto const ei = energies[i] - eF;
              dos += FermiDirac_broadening((e - ei)*kTinv);
          } // i
          if (echo > 8) printf("%g %g %g\n", e*eV, dos*kTinv*per_eV, ddos);
          integral += dos*FermiDirac(e*kTinv);
      } // ie energy
      if (echo > 6) printf("# integrated DoS is %g\n", integral*kTinv*de);
      return 0;
  } // density_of_states

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_integration(int const echo=3, int const n=580) {
      if (echo > 6) printf("\n## %s %s (%d points)\n", __FILE__, __func__, n);
      double s0{0}, s1{0}, dfde;
      double const de = .0625;
      for(int ie = -n; ie < n; ++ie) { double const e = (ie + .5)*de;
          double const f = FermiDirac(e, &dfde);
          if (echo > 6) printf("%g %g %g\n", e, f, -dfde);
          s0 += f;
          s1 -= dfde;
      } // e
      if (echo > 2) printf("# Fermi-Dirac distribution deviation %.1e %.1e\n", s0/n - 1, s1*de - 1);
      return 0;
  } // test_integration

  inline status_t test_bisection(int const echo=3, int const n=99) {
      if (echo > 6) printf("\n## %s %s (%d states)\n", __FILE__, __func__, n);
      std::vector<double> ene(n);
      for(int i = 0; i < n; ++i) ene[i] = simple_math::random(-10., 10.);
      double eF{0};
      auto const stat = Fermi_level(eF, nullptr, ene.data(), n, 3e-2, n*.75, 1, echo);
      if (echo > 8) density_of_states(n, ene.data(), 3e-2, eF);
      return stat;
  } // test_bisection

  inline status_t all_tests(int const echo=0) {
      if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
      status_t status(0);
      status += test_integration(echo);
      status += test_bisection(echo);
      return status;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace fermi_distribution
