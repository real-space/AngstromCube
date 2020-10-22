#pragma once

#include <cstdio> // printf
#include <cassert> // assert
#include <cmath> // std::exp
#include <vector> // std::vector<T>

#include "inline_math.hxx" // set
#include "simple_math.hxx" // ::random
#include "display_units.h" // eV, _eV
#include "recorded_warnings.hxx" // warn

#include "status.hxx" // status_t

namespace fermi_distribution {
  
  template <typename real_t> inline
  real_t FermiDirac(real_t const x, real_t *derivative=nullptr) { 
      if (std::abs(x) > 36) {
          if (derivative) *derivative = 0;
          return real_t(x < 0);
      } // |x| large
      real_t const d = 1 + std::exp(x); // denominator
      real_t const f = 1/d;
      if (derivative) *derivative = (1 - d)*f*f;
      return f;
  } // FermiDirac

  template <typename real_t> inline
  real_t FermiDirac_broadening(real_t const x) {
      real_t der;
      FermiDirac(x, &der);
      return -der;
  } // FermiDirac_broadening

  template <typename real_t>
  double count_electrons(
        int const n
      , real_t const energies[]
      , double const eF
      , double const kTinv
      , double const weights[]=nullptr
      , double *derivative=nullptr
      , double occupations[]=nullptr
      , double ddeF_occupations[]=nullptr
  ) {
      double const ddeF_x = -kTinv;
      double ne{0}, ddx_ne{0};
      assert(nullptr == weights); // not implemented!
      for(int i = 0; i < n; ++i) {
          double ddx_f;
          double const x = (energies[i] - eF)*kTinv;
          double const f = FermiDirac(x, &ddx_f);
          double const w8 = weights ? weights[i] : 1;
          ne     += f     *w8;
          ddx_ne += ddx_f *w8;
          if (occupations) occupations[i] = f;
          if (ddeF_occupations) ddeF_occupations[i] = ddx_f * ddeF_x;
      } // i
      if (derivative) *derivative = ddx_ne * ddeF_x; // derivative w.r.t. eF
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
      , double response_occ[]=nullptr
      , double *DoS_at_eF=nullptr
  ) {
      double const n_electrons = number_of_electrons/spin_factor;
      // ToDo: count the states with weights, check if there are enough states to host n_electrons
      std::vector<real_t> ene(nb); // energies
      set(ene.data(), nb, energies); // copy
      std::sort(ene.begin(), ene.end()); // sort to ascending order
      auto const e_min = ene[0], e_max = ene[nb - 1];
      if (echo > 0) printf("# %s E_min= %g E_max= %g %s\n", __func__, e_min*eV, e_max*eV, _eV);

      { // scope: find T=0 Fermi level, works only for a sorted spectrum and ignores weights
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
      while(res > 2e-15 && iter < maxiter) {
          ++iter;
          auto const em = 0.5*(e[0] + e[1]);
          auto const nem = count_electrons(nb, energies, em, kTinv);
          if (echo > 7) {
//            printf("# %s with energy %g %s\t--> %g electrons\n", __func__, em*eV, _eV, spin_factor*nem);
              auto const ene = std::abs(e[1] - e[0]); // energy residual
              int const ndigits_energy = std::min(std::max(1, int(-std::log10(ene)) - 1), 15);
              int const ndigits_charge = std::min(std::max(1, int(-std::log10(res)) - 1), 15);
//            printf("# %s ndigits= %d and %d\n", __func__, ndigits_energy, ndigits_charge);
              char fmt[64]; std::snprintf(fmt, 63, "# %s with energy %%.%df %s --> %%.%df electrons\n",
                                                      __func__, ndigits_energy, _eV, ndigits_charge);
//            printf("# %s fmt= %s\n", __func__, fmt);
              printf(fmt, em*eV, spin_factor*nem);
              fflush(stdout);
          } // echo
          int const i01 = (nem > n_electrons);
          e[i01] = em;
          ne[i01] = nem;
          res = std::abs(nem - n_electrons); // charge residual
      }
      eF = 0.5*(e[0] + e[1]);
      double DoS; // density of states at the Fermi level
      auto const nem = count_electrons(nb, energies, eF, kTinv, nullptr, &DoS, occupations, response_occ);
      if (echo > 6) printf("# %s with energy %g %s --> %g electrons\n", __func__, eF*eV, _eV, spin_factor*nem); 
      if (res > 1e-9) {
          warn("Fermi level converged only to +/- %.1e electrons in %d iterations", res*spin_factor, iter); 
      } else {
          if (echo > 3) printf("# %s at %.9f %s has a DOS of %g states\n", __func__, eF*eV, _eV, DoS*spin_factor);
      }
      if (DoS_at_eF) *DoS_at_eF = DoS;
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

  
  class FermiLevel_t
  { 
      // group all the information needed for the determination
      // of the Fermi level during SCF cycles
    public:
      static double constexpr Fermi_level_not_initialized = -9e99;
      static double constexpr default_temperature = 9.765625e-4; // k*T in Hartree units

      FermiLevel_t( // constructor
            double const n_electrons
          , int const spin_factor=2
          , double const kT=default_temperature
          , int const echo=9 // log-level
      )
        : _ne(std::max(0.0, n_electrons))
        , _mu(Fermi_level_not_initialized)
        , _kT(kT)
        , _spin_factor(std::min(std::max(1, spin_factor), 2))
      {
          if (echo > 0) printf("\n# new %s(%g electrons, kT=%g %s)\n", __func__, _ne, kT*Kelvin, _Kelvin);
          set_temperature(kT, echo);
          set_accumulators();
      } // constructor

      bool set_temperature(double const kT=default_temperature, int const echo=0) {
          // change the temperature during SCF cycles
          _kT = std::max(minimum_temperature(), kT);
          if (_kT > 0) _kTinv = 1./_kT;
          if (echo > 0) printf("# FermiLevel %s(%g) --> kT= %g %s\n", __func__, kT*Kelvin, _kT*Kelvin, _Kelvin);
          return (kT >= 0);
      } // set_temperature

      inline double minimum_temperature() const { return 1e-9; }

      inline bool is_initialized() const { return (Fermi_level_not_initialized != _mu); }

      void set_accumulators(double const values[]=nullptr, int const echo=0) {
          if (nullptr == values) {
              set(_accu, 4, 0.0); // reset accumulators
              return;
          }
          set(_accu, 4, values);
          if (echo > 3) printf("# FermiLevel %s(%g, %g %s, %g, %g)\n",
                  __func__,  _accu[0], _accu[1]*eV,_eV, _accu[2], _accu[3]);
      } // set_accumulators
      
      void set_Fermi_level(double const E_Fermi, int const echo=0) {
          double const a[] = {1, E_Fermi, 0, 0};
          set_accumulators(a);
          if (echo > 0) printf("# FermiLevel set to %g %s\n", E_Fermi*eV,_eV);
      } // set_Fermi_level

      double average() const { return (_accu[0] > 0) ? _accu[1]/_accu[0] : 0; }
      
      double add(double const E_Fermi, double const weight=1, int const echo=0) {
          double const w8 = std::max(weight, 0.0);
          _accu[0] += w8;
          _accu[1] += w8 * E_Fermi;
          _accu[2] += w8 * pow2(E_Fermi);
          _accu[3] += w8 * pow3(E_Fermi);
          auto const avg = average();
          if (echo > 0) printf("# FermiLevel %s(%g, weight=%g), new average= %g %s\n",
                            __func__, E_Fermi*eV, w8, avg*eV,_eV);
          return avg;
      } // add
      
      template <typename real_t>
      double get_occupations( // returns the density of state contribution at the Fermi level
            double occupations[] // output occupation numbers[nbands]
          , real_t const energies[] // energies[nbands]
          , int const nbands // number of bands
          , int const echo=9
          , double response_occ[]=nullptr // optional output, derivative of occupation numbers[nbands]
      ) {
            if (echo > 0) printf("# FermiLevel %s for %d bands\n", __func__, nbands);
            double DoS{0}; // density of states at the Fermi energy
            if (Fermi_level_not_initialized == _mu) {
                if (echo > 0) printf("# FermiLevel has not been initialized\n", __func__, nbands);
                double eF;
                Fermi_level(eF, occupations, energies, nbands, _kT, _ne, _spin_factor, echo, response_occ, &DoS);
                _mu = eF;
            } else {
                if (echo > 0) printf("# FermiLevel at %g %s %s\n", _mu*eV,_eV, __func__);
                count_electrons(nbands, energies, _mu, _kTinv, nullptr, &DoS, occupations, response_occ);
            } // initialized?
            return DoS;
      } // get_occupations
      
      double get_Fermi_level() const { return _mu; }

    private:

      // member variables
      double _ne; // number of electrons
      double _mu; // the chemical potential
      double _kT; // temperature times Boltzmann factor
      double _kTinv;
      double _accu[4]; // Fermi level accumulator
      int _spin_factor;

  }; // class FermiLevel_t

  
  
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

  inline status_t test_FermiLevel_class(int const echo=0) {
      if (echo > 6) printf("\n# %s %s\n", __FILE__, __func__);
      double const n_electrons = 12;
      auto mu = FermiLevel_t(n_electrons);
      double const a[] = {1., -1., 2., 3.3};
      mu.set_temperature(.01, echo);
      mu.set_accumulators(a, echo);
      mu.set_accumulators(); // reset
      for(double eF = -1.5; eF < 2; eF += 0.25) {
          mu.add(eF, 1, echo - 5);
      } // eF
      mu.set_Fermi_level(mu.average(), echo);
      int const nb = n_electrons*1.5; // number of bands
      std::vector<double> ene(nb, 0.0), occ(nb), docc(nb);
      for(int ib = 0; ib < nb; ++ib) ene[ib] = ib*.01;
      mu.get_occupations(occ.data(), ene.data(), nb, echo, docc.data());
      double sum_occ{0};
      for(int ib = 0; ib < nb; ++ib) {
          if (echo > 9) printf("# %s: band #%d \tenergy= %.6f %s \toccupation= %.6f \td/dE occ= %.3f\n",
                __func__, ib, ene[ib]*eV,_eV, occ[ib], docc[ib]);
          sum_occ += occ[ib];
      } // ib
      if (echo > 6) printf("# %s %s eF= %g %s for %g electrons, sum(occupations)= %g\n",
                        __FILE__, __func__, mu.get_Fermi_level()*eV,_eV, n_electrons, sum_occ);
      mu.get_occupations(occ.data(), ene.data(), nb, echo, docc.data());
      return 0;
  } // test_FermiLevel_class

  inline status_t all_tests(int const echo=0) {
      if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
      status_t status(0);
      status += test_integration(echo);
      status += test_bisection(echo);
      status += test_FermiLevel_class(echo);
      return status;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace fermi_distribution
