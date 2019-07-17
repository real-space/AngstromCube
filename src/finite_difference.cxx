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
  
  template<typename real_t, int D0>
  status_t Laplacian(real_t *out, real_space_grid::grid_t<real_t,D0> const &g, 
                     finite_difference_t<real_t> const &fd) {
      int const n16 = 16; // max number of finite difference neighbors
      int* list[3];
      for(int d = 0; d < 3; ++d) {
          int const n = g.dim(d);
          int const nf = fd.nn[d];
          int const nh = n16 + n + n16; // number including largest halos
          list[d] = new int[nh]; // get memory
          for(int i = 0; i < nh; ++i) list[d][i] = -1; // init as non-existing
          for(int i = 0; i < n; ++i) list[d][n16 + i] = i; // itself
          // lower boundary
          if (1 == fd.bc[d][0]) { // periodic BC
              for(int i = n16 - nf; i < n16; ++i) list[d][i] = n + i - n16; // wrap around
          } else if (-1 == fd.bc[d][0]) { // mirror BC
              for(int i = n16 - nf; i < n16; ++i) list[d][i] = n16 - 1 - i; // wrap around and mirror
          } // else open BC, list[:] = -1
          if (1 == fd.bc[d][1]) { // periodic BC
              for(int i = 0; i < nf; ++i) list[d][n16 + n + i] = i; // wrap around
          } else if (-1 == fd.bc[d][1]) { // mirror BC
              for(int i = 0; i < nf; ++i) list[d][n16 + n + i] = n - 1 - i; // wrap around and mirror
          } // else open BC, list[:] = -1
          printf("# indirection list for %c  ", 120+d);
          for(int i = n16 - nf; i < n16 + n + nf; ++i) {
              if ((n16 == i) || (n16 + n) == i) printf(" |");
              printf(" %d", list[d][i]);
          }   printf("\n");
      } // spatial direction d

      assert(1 == D0); // no vectorization active

      for(int z = 0; z < g.dim('z'); ++z) {
          for(int y = 0; y < g.dim('y'); ++y) {
              for(int x = 0; x < g.dim('x'); ++x) {
                  int const i_zyx = (z*g.dim('y') + y)*g.dim('x') + x;
                  out[i_zyx] = 0; // init result

                  for(int ddir = 0; ddir < 3; ++ddir) {
                      int zyx[3] = {x, y, z};
                      int const i_center = zyx[ddir];
                      for(int jmi = -fd.nn[ddir]; jmi <= fd.nn[ddir]; ++jmi) {
                          int const j = i_center + jmi;
                          int const index = list[ddir][n16 + j];
                          if (index > 0) {
                              zyx[ddir] = index;
                              int const zyx_prime = (zyx[2]*g.dim('y') + zyx[1])*g.dim('x') + zyx[0];
                              out[i_zyx] += g.values[zyx_prime] * fd.c2nd[ddir][std::abs(jmi)];
                          } // index exists
                      } // jmi
                  } // ddir direction of the derivative

              } // x
          } // y
      } // z

      return 0; // success
  } // Laplacian
  
  
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
  status_t test_Laplacian(int const echo=9) {
      finite_difference_t<real_t> fd;
      int const dims[] = {10, 11, 12};
      real_space_grid::grid_t<real_t,1> g(dims);
      for(int i = 0; i < g.all(); ++i) g.values[i] = std::cos(i/32.); // fill with some non-zero values
      printf("\n# in  values: \n"); for(int i = 0; i < g.all(); ++i) printf("%g\n", g.values[i]); 
      real_t* out = new real_t[g.all()];
      auto const stat = Laplacian(out, g, fd);
      printf("\n# out values: \n"); for(int i = 0; i < g.all(); ++i) printf("%g\n", out[i]); 
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
