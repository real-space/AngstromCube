#include <cstdio> // printf
#include <cassert> // assert
#include <cmath> // std::abs
#include <algorithm> // std::copy

#include "finite_difference.hxx"

#include "real_space_grid.hxx" // grid_t
#include "constants.hxx" // pi

// #define FULL_DEBUG
// #define DEBUG

namespace finite_difference {
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template<typename real_t>
  status_t test_coefficients(int const echo=2) {
      status_t stat = 0;
      int const mantissa_bits = (sizeof(real_t) > 4)? 52 : 23; // 23:float, 52:double
      double const precision = 4./(1ul << mantissa_bits);
      if (echo > 4) printf("# expected precision for real_%ld is %g\n", sizeof(real_t), precision);
      int const M = 15;
      double maxdev = 0; int maxdev_nn = -99;
      real_t c[M];
      for(int nn = M - 1; nn >= -1; --nn) {
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
      auto const f = new finite_difference_t<float>();
      f->~finite_difference_t();
      finite_difference_t<double> d;
      return 0;
  } // test_create_and_destroy

  template<typename real_t>
  status_t test_Laplacian(int const echo=3) {
      status_t stat = 0;
      int const bc[3] = {1,1,1}; // periodic
      double const h[3] = {1,1,1};
      for(int dir = 0; dir < 3; ++dir) {
          int nn[3] = {0,0,0}; nn[dir] = 12; // switch FD off for the two perpendicular directions
          finite_difference_t<real_t> fd(h, bc, nn);
          int dims[] = {0,0,0}; dims[dir] = 127 + dir;
          real_space_grid::grid_t<real_t,1> g(dims);
          double const k = (1 + dir)*2*constants::pi/g.dim(dir);
          for(int i = 0; i < g.all(); ++i) g.values[i] = std::cos(k*i); // fill with some non-zero values
          real_t* out = new real_t[g.all()];
          stat += Laplacian(out, g, fd);
          if (echo > 5) printf("\n# in, out, ref values:\n");
          double dev = 0;
          for(int i = 0; i < g.all(); ++i) {
              double const ref = -k*k*g.values[i];
              if (echo > 5) printf("%d %g %g %g\n", i, g.values[i], out[i], ref);
              // compare in the middle range out and ref values
              dev += std::abs(out[i] - ref);
          } // i
          if (echo > 2) printf("# %c-direction: dev = %g\n", 120+dir, dev);
      } // direction
      return stat;
  } // test_Laplacian

  
  status_t all_tests() {
    auto status = 0;
    status += test_coefficients<double>();
    status += test_coefficients<float>();
    status += test_create_and_destroy();
    status += test_Laplacian<double>();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace finite_difference
