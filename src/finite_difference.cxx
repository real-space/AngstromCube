#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy

#include "finite_difference.hxx"

#include "real_space_grid.hxx" // grid_t
#include "constants.hxx" // pi

// #define FULL_DEBUG
// #define DEBUG

namespace finite_difference {
  
  template<typename real_t, int D0>
  status_t Laplacian(real_t *out, real_space_grid::grid_t<real_t,D0> const &g, 
                     finite_difference_t<real_t> const &fd) {
      int* list[3];
      for(int d = 0; d < 3; ++d) {
          int const n = g.dim(d);
          int const nh = 16 + n + 16; // number including largest halos
          list[d] = new int[nh]; // get memory
          for(int i = 0; i < n; ++i) list[d][16 + i] = i; // itself
          if (0 == fd.bc[d][0]) { // open lower BC
              for(int i = 0; i < 16; ++i) list[d][i] = -1; // grid points are zero / do not exist
          } else if (1 == fd.bc[d][0]) { // periodic lower BC
              for(int i = 0; i < 16; ++i) list[d][i] = n + i - 16; // wrap around
          } else if (-1 == fd.bc[d][0]) { // mirrored lower BC
              for(int i = 0; i < 16; ++i) list[d][i] = 16 - 1 - i; // wrap around and mirror
          }
          if (0 == fd.bc[d][1]) { // open upper BC
              for(int i = 0; i < 16; ++i) list[d][16 + n + i] = -1; // grid points are zero / do not exist
          } else if (1 == fd.bc[d][1]) { // periodic upper BC
              for(int i = 0; i < 16; ++i) list[d][16 + n + i] = i; // wrap around
          } else if (-1 == fd.bc[d][1]) { // mirrored upper BC
              for(int i = 0; i < 16; ++i) list[d][16 + n + i] = n - 1 - i; // wrap around and mirror
          }
      } // spatial direction d
      return 0; // success
  } // Laplacian
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template<typename real_t>
  status_t test_coefficients(int echo=2) {
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

  status_t test_create_and_destroy(int echo=9) {
      auto const f = new finite_difference_t<float>();
      f->~finite_difference_t();
      finite_difference_t<double> d;
      return 0;
  } // test_create_and_destroy
  
  status_t all_tests() {
    auto status = 0;
    status += test_coefficients<double>();
    status += test_coefficients<float>();
    status += test_create_and_destroy();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace finite_difference
