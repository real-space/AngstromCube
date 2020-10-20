#include <cstdio> // printf
#include <cassert> // assert
#include <cmath> // std::abs
#include <algorithm> // std::copy
#include <vector> // std::vector<T>

#include "finite_difference.hxx"

#include "real_space.hxx" // ::grid_t
#include "constants.hxx" // ::pi
#include "recorded_warnings.hxx" // ::clear_warnings
#include "boundary_condition.hxx" // Periodic_Boundary

namespace finite_difference {
  
#ifndef NO_UNIT_TESTS

  template<typename real_t>
  status_t test_coefficients(int const echo=2) {
      status_t stat = 0;
      int const mantissa_bits = (sizeof(real_t) > 4)? 52 : 23; // 23:float, 52:double
      double const precision = 4./(1ul << mantissa_bits);
      if (echo > 4) printf("# expected precision for real_%ld is %g\n", sizeof(real_t), precision);
      int const M = 16;
      double maxdev = 0; int maxdev_nn = -99;
      real_t c[M];
      for(int nn = M - 2; nn >= -1; --nn) {
          for(int i = 0; i < M; ++i) c[i] = 0; // clear
          set_Laplacian_coefficients(c, nn);
          double checksum = c[0]; for(int i = 1; i < M; ++i) checksum += 2*c[i];
          if (echo > 6) printf("# Laplacian with nn=%d has c0 = %g,   %.1e should be zero\n", nn, c[0], checksum);
          checksum = std::abs(checksum);
          stat += (checksum > precision);
          if (checksum > maxdev) { maxdev = checksum; maxdev_nn = nn; } // get the maximum and the nn where it occurred
      } // n
      if (stat && (echo > 0)) printf("# %s found %d errors!\n", __func__, stat);
      if (echo > 3) printf("# Laplacian with nn=%d has largest deviation in checksum: %.1e (should be zero)\n", maxdev_nn, maxdev);
      return stat;
  } // test_coefficients

  status_t test_create_and_destroy(int const echo=9) {
      auto const f = new stencil_t<float>();
      f->~stencil_t();
      stencil_t<double> d;
      return 0;
  } // test_create_and_destroy

  template<typename real_t>
  status_t test_Laplacian(int const echo=3) {
      status_t stat(0);
      double const h[3] = {1,1,1};
      for(int dir = 0; dir < 3; ++dir) {
          int nn[3] = {0,0,0}; nn[dir] = 12; // switch FD off for the two perpendicular directions
          stencil_t<real_t> Laplacian(h, nn);
          int dims[] = {1,1,1}; dims[dir] = 127 + dir;
          real_space::grid_t g(dims);
          g.set_boundary_conditions(Periodic_Boundary);
          double const k = (1 + dir)*2*constants::pi/g[dir]; // wave vector of a single plane wave
          std::vector<real_t> values(g.all()), result(g.all());
          for(size_t i = 0; i < g.all(); ++i) values[i] = std::cos(k*i); // fill with some non-zero values
          stat += finite_difference::apply(result.data(), values.data(), g, Laplacian);
          if (echo > 5) printf("\n# in, result, ref values:\n");
          double dev{0};
          for(size_t i = 0; i < g.all(); ++i) {
              double const ref = -k*k*values[i]; // analytic solution to the Laplacian operator applied to a plane wave
              if (echo > 5) printf("%ld %g %g %g\n", i, values[i], result[i], ref);
              // compare in the middle range result and ref values
              dev += std::abs(result[i] - ref);
          } // i
          if (echo > 2) printf("# %s %c-direction: dev = %g\n", __func__, 120+dir, dev);
      } // direction
      return stat;
  } // test_Laplacian

  status_t test_dispersion(int const echo=9) {
      if (echo < 7) return 0; // this function is only plotting
      for(int nn = 1; nn <= 13; ++nn) { // largest order implemented is 13
          stencil_t<double> const fd(1.0, nn);
          printf("\n## finite-difference dispersion for nn=%d\n", nn);
          for(int ik = 0; ik <= 100; ++ik) {
              double const k = 0.01 * ik * constants::pi;
              double E_k{-0.5*fd.c2nd[0][0]};
              for(int j = 1; j <= nn; ++j) {
                  E_k -= std::cos(k*j) * fd.c2nd[0][j];
              } // j
              printf("%g %g %g\n", k, E_k, 0.5*k*k); // in Hartree
          } // ik
      } // nn
      return 0;
  } // test_dispersion
  
  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_coefficients<double>(echo);
      stat += test_coefficients<float>(echo);
      stat += test_create_and_destroy(echo);
      stat += test_Laplacian<double>(echo);
      recorded_warnings::clear_warnings(); // clear
      stat += test_dispersion(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace finite_difference
