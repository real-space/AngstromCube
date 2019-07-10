#include <cstdio> // printf
#include <cmath> // sqrt, pow, log
#include <math.h> // C_PI
#include <cassert> // assert

#include "exchange_correlation.hxx"

#include "output_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // pow2

// #define FULL_DEBUG
// #define DEBUG

namespace exchange_correlation {
  
  // ToDo: use libxc in the long run
  template<typename real_t> 
  real_t lda_PZ81_kernel(real_t const rho, real_t &Vup, 
                         real_t const mag, real_t *Vdn) { // optional

    real_t constexpr C_PI = 3.14159265358979323846; // pi
    real_t constexpr THIRD  = 1./3., TINYDEN = 1e-20;
    real_t const tpt5 = .6108870577108572; // (2.25/(C_PI*C_PI))**THIRD

    if (nullptr == Vdn) { //  LDA

      if (rho < TINYDEN) {
        Vup = 0;
        return 0;
      } else { // negligible
        real_t Exc;
        auto const rs = pow(3.0/(4.0*C_PI*rho), THIRD);
        if (rs > 1.) {
          auto const srs = sqrt(rs);
          Exc = -0.1423/(1. + 1.0529*srs + 0.3334*rs);
          Vup = Exc - rs/3.0*(0.1423*(0.3334 + 0.52645/srs)/pow2(1. + 1.0529*srs + 0.3334*rs));
        } else {
          auto const lrs = log(rs);
          Exc = -0.048 + 0.0311*lrs - 0.0116*rs + 0.002*rs*lrs;
          Vup = Exc - rs/3.0*(0.0311/rs - 0.0096 + 0.002*lrs);
        } //
        Exc = Exc - 0.75*tpt5/rs;
        Vup = Vup - 0.75*tpt5/rs - rs/3.0*(0.75*tpt5/(rs*rs));
        return Exc;
      } // negligible

    } else { // spin resolved

      if (rho < TINYDEN) {
        *Vdn = 0;
         Vup = 0;
        return 0;
      } else { // negligible
        real_t Exc = 0;
         Vup = 0;
        *Vdn = 0;
        return Exc;
      } // negligible
    } // Magnetization

  } // lda_PZ81_kernel
  
  // explicit template instanciation for float and double
  template double lda_PZ81_kernel<double>(double const rho, double &Vup, double const mag, double *Vdn);
  template float  lda_PZ81_kernel<float> (float  const rho, float  &Vup, float  const mag, float  *Vdn);

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_PZ81_potential() {
    for(double rho = .625e-20; rho < 1e9; rho *= 2) {
        double V; double const E = lda_PZ81_kernel(rho, V);
        printf("%g %g %g\n", rho, E, V);
    } // rho
    return 0;
  } // test_PZ81_potential

  status_t all_tests() {
    auto status = 0;
    status += test_PZ81_potential();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace exchange_correlation
