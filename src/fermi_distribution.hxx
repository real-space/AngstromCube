#pragma once

#include <cstdio> // printf
#include <cassert> // assert
#include <cmath> // std::exp

typedef int status_t;

namespace fermi_distribution {
  
  template<typename real_t>
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
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
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

  inline status_t all_tests(int const echo=0) {
      if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
      auto status = 0;
      status += test_integration(echo);
      return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace fermi_distribution
