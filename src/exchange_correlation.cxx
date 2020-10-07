#include <cmath> // std::sqrt, std::pow, std::log
#include <cassert> // assert
#include <cstdio> // printf

#ifndef  NO_UNIT_TESTS
  #include "display_units.h" // eV, _eV, Ang, _Ang
#endif

#include "exchange_correlation.hxx"

#include "constants.hxx" // ::pi
#include "inline_math.hxx" // pow2

namespace exchange_correlation {
  
  // ToDo: use libxc in the long run
  template<typename real_t>
  real_t lda_PZ81_kernel(real_t const rho, real_t &Vdn,
       real_t const mag, real_t *Vup) { // defaults: mag=0, *Vup=nullptr

      real_t constexpr THIRD = 1./3., TINYDEN = 1e-20;
      real_t constexpr pi = constants::pi;
      real_t const tpt5 = .6108870577108572; // (2.25/(pi*pi))**THIRD

      if (nullptr == Vup) { // LDA, no magnetiziation

          if (rho < TINYDEN) {
              Vdn = 0;
              return 0;
          } else { // negligible
              real_t Exc = 0;
              auto const rs = std::pow(3.0/(4.0*pi*rho), THIRD);
              if (rs > 1.) {
                  auto const srs = std::sqrt(rs);
                  Exc = -0.1423/(1. + 1.0529*srs + 0.3334*rs);
                  Vdn = Exc - rs/3.0*(0.1423*(0.3334 + 0.52645/srs)/pow2(1. + 1.0529*srs + 0.3334*rs));
              } else {
                  auto const lrs = std::log(rs);
                  Exc = -0.048 + 0.0311*lrs - 0.0116*rs + 0.002*rs*lrs;
                  Vdn = Exc - rs/3.0*(0.0311/rs - 0.0096 + 0.002*lrs);
              } //
              Exc -=  0.75*tpt5/rs;
              Vdn -= (0.75*tpt5/rs + rs/3.0*(0.75*tpt5/(rs*rs)));
              return Exc;
          } // negligible

      } else { // LDA, no magnetiziation
        
          assert(0 == __LINE__); // will fail!

          if (rho < TINYDEN) {
              *Vup = 0;
               Vdn = 0;
              return 0;
          } else { // negligible
              *Vup = 0;
               Vdn = 0;
              return 0;
          } // negligible

      } // LDA, no magnetiziation

  } // lda_PZ81_kernel

  template // explicit template instantiation for double
  double lda_PZ81_kernel<double>(double const rho, double &Vdn, double const mag, double *Vup);

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_PZ81_potential(int const echo=0) {
      if (echo > 0) printf("\n## density(%s^-3) energy(%s) potential(%s)\n", _Ang, _eV, _eV);
      for(double rho = .625e-20; rho < 1e9; rho *= 2) {
          double V; double const E = lda_PZ81_kernel(rho, V);
          if (echo > 0) printf("%g %g %g\n", rho/pow3(Ang), E*eV, V*eV);
      } // rho
      return 0;
  } // test_PZ81_potential

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_PZ81_potential(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace exchange_correlation
