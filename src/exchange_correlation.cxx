#include <cmath> // std::sqrt, std::pow, std::log
#include <cassert> // assert
#include <cstdio> // std::printf

#ifndef  NO_UNIT_TESTS
  #include "display_units.h" // eV, _eV, Ang, _Ang
#endif

#include "exchange_correlation.hxx"

#include "constants.hxx" // ::pi
#include "inline_math.hxx" // pow2

namespace exchange_correlation {

  float constexpr TINY_DENSITY = 1e-20;
  double constexpr THIRD = 1./3.;

  // ToDo: use libxc in the long run
  template <typename real_t>
  real_t lda_PZ81_kernel(
        real_t const rho
      , real_t & Vdn
      , real_t const mag  // =0
      , real_t *Vup       // =nullptr
  ) {
      real_t const tpt5 = .6108870577108572; // (2.25/(pi*pi))**THIRD

      if (nullptr == Vup) { // LDA, no magnetiziation

          if (rho < TINY_DENSITY) {
              Vdn = 0;
              return 0;
          } else { // negligible
              real_t Exc{0};
              auto const rs = std::pow(3.0/(4.0*constants::pi*rho), THIRD);
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

          *Vup = 0;
           Vdn = 0;
          return 0;

      } // LDA, no magnetiziation

  } // lda_PZ81_kernel



  template <typename real_t>
  inline real_t G(
        real_t const rtrs
      , real_t const A
      , real_t const alpha1
      , real_t const beta1
      , real_t const beta2
      , real_t const beta3
      , real_t const beta4
      , real_t & dGdrs
  ) {
      real_t const Q0 = -2.0 * A * (1.0 + alpha1 * rtrs * rtrs);
      real_t const Q1 =  2.0 * A *  rtrs * (beta1 + 
                                        rtrs * (beta2 + 
                                            rtrs * (beta3 + 
                                                rtrs * beta4)));
      real_t const G1 = Q0 * std::log(1.0 + 1.0 / Q1);
      real_t const dQ1drs = A * (beta1 / rtrs + 2.0 * beta2 +
                          rtrs * (3.0 * beta3 + 4.0 * beta4 * rtrs));
      dGdrs = -2.0 * A * alpha1 * G1 / Q0 - Q0 * dQ1drs / (Q1 * (Q1 + 1.0));
      return G1;
  } // G

  template <typename real_t>
  real_t pw91_exchange(real_t rs, real_t & dedrs) {
      real_t const e = -0.45816529328314287 / rs;
      dedrs = -e / rs;
      return e;
  } // pw91_exchange

  template <typename real_t>  
  real_t pw91_correlation(real_t rs, real_t & dedrs) {
      real_t const rtrs = std::sqrt(rs);
      return G(rtrs, 0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294, dedrs);
  } // pw91_correlation

  template <typename real_t>
  real_t lda_PW91_kernel(        
        real_t const rho
      , real_t & Vdn
      , real_t const mag  // =0
      , real_t *Vup       // =nullptr
  ) {

      if (nullptr == Vup) { // LDA, no magnetiziation

          if (rho < TINY_DENSITY) {
              Vdn = 0;
              return 0;
          } else { // negligible
              auto const rs = std::pow(3.0/(4.0*constants::pi*rho), THIRD);
              real_t dexdrs, decdrs;
              auto const ex = pw91_exchange(rs, dexdrs);
              auto const ec = pw91_correlation(rs, decdrs);
              Vdn = ex + ec - rs * (dexdrs + decdrs) / 3.0;
              return ex + ec;
          } // negligible

      } else { // LDA, no magnetiziation

          assert(0 == __LINE__); // will fail!

          *Vup = 0;
           Vdn = 0;
          return 0;

      } // LDA, no magnetiziation

  } // lda_PW91_kernel

  template // explicit template instantiation for double
  double lda_PW91_kernel<double>(double const rho, double &Vdn, double const mag, double *Vup);

  template // explicit template instantiation for double
  double lda_PZ81_kernel<double>(double const rho, double &Vdn, double const mag, double *Vup);


#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_LDA_potentials(int const echo=0) {
      if (echo > 2) std::printf("\n## density(%s^-3) energy(%s) potential(%s) for {PW91, PZ81}\n", _Ang, _eV, _eV);
      for (double rho = .625e-20; rho < 1e9; rho *= 2) {
          double vPZ; double const ePZ = lda_PZ81_kernel(rho, vPZ);
          double vPW; double const ePW = lda_PW91_kernel(rho, vPW);
          if (echo > 2) std::printf("%g %g %g %g %g\n", rho/pow3(Ang), ePW*eV, vPW*eV, ePZ*eV, vPZ*eV);
      } // rho
      if (echo > 2) std::printf("\n# use gnuplot or xmgrace -nxy to plot the columns above\n");
      return 0;
  } // test_LDA_potentials

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_LDA_potentials(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace exchange_correlation
