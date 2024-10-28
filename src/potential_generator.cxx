// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert
#include <vector> // std::vector

#include "potential_generator.hxx" // ::add_generalized_Gaussian

#include "status.hxx" // status_t
#include "control.hxx" // ::get
#include "display_units.h" // Ang, _Ang

namespace potential_generator {

#ifdef    NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

  status_t test_generalized_Gaussian(int const echo=0) {
      status_t stat(0);
      double const sigma = control::get("potential_generator.test.generalized.gaussian.sigma", 1.);
      if (echo > 1) std::printf(      "# potential_generator.test.generalized.gaussian.sigma=%g %s\n", sigma*Ang, _Ang);
      int const n        = control::get("potential_generator.test.generalized.gaussian.n", 24*sigma);
      if (echo > 1) std::printf(      "# potential_generator.test.generalized.gaussian.n=%d grid points\n", n);
      int const ellmax   = control::get("potential_generator.test.generalized.gaussian.ellmax", 3.);
      if (echo > 1) std::printf(      "# potential_generator.test.generalized.gaussian.ellmax=%d\n", ellmax);
      assert(ellmax >= 0);
      assert(sigma > 0);
      assert(n > 0);
      real_space::grid_t g(n, n, n); // grid spacing is unity
      double const center[] = {.5*(n - 1), .5*(n - 1), .5*(n - 1)};
      std::vector<double> coeff((1 + ellmax)*(1 + ellmax), 0.0);
      std::vector<double> values(g.all(), 0.0);
      int constexpr debug = 1;
      stat += add_generalized_Gaussian<double,debug>(values.data(), g, coeff.data(), ellmax, center, sigma, echo);
      return stat;
  } // test_generalized_Gaussian

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_generalized_Gaussian(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace potential_generator
