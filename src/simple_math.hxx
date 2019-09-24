#pragma once

#include <random> // std::random_device, std::mt19937

  typedef int status_t;
  
namespace simple_math {

  // templated, performant, thread safe pseudo random number generator based on Mersenne twister
  template<typename real_t, typename Generator=std::mt19937>
  real_t random(real_t const from, real_t const to) {

      thread_local static Generator gen(std::random_device{}());

      using dist_type = typename std::conditional
      <   std::is_integral<real_t>::value
        , std::uniform_int_distribution<real_t>
        , std::uniform_real_distribution<real_t>
      >::type;

      thread_local static dist_type dist;

      return dist(gen, typename dist_type::param_type{from, to});
  } // random

  
  template<typename real_t>
  real_t determinant1x1(real_t const a[], int const as=1) { return a[0*as + 0]; } // determinant

  template<typename real_t>
  real_t invert1x1(real_t inv[], int const is, real_t const a[], int const as=1) {
    auto const det = determinant1x1(a, as);
    if (std::abs(det) < 1e-16) return 0;
    auto const inv_det = 1/det;
    inv[0*is + 0] = a[0*as + 0]*inv_det;
    return det;
  } // invert
  
  template<typename real_t>
  real_t determinant2x2(real_t const a[], int const as=2) {
    return a[0*as + 0]*a[1*as + 1] - a[1*as + 0]*a[0*as + 1]; 
  } // determinant

  template<typename real_t>
  real_t invert2x2(real_t inv[], int const is, real_t const a[], int const as=2) {
    auto const det = determinant2x2(a, as);
    if (std::abs(det) < 1e-16) return 0;
    auto const inv_det = 1/det;
    inv[0*is + 0] =   a[1*as + 1]*inv_det;
    inv[0*is + 1] = - a[0*as + 1]*inv_det;
    inv[1*is + 0] = - a[1*as + 0]*inv_det;
    inv[1*is + 1] =   a[0*as + 0]*inv_det;
    return det;
  } // invert

  template<typename real_t>
  real_t determinant3x3(real_t const a[], int const as=3) {
    return a[0*as + 0]*(a[1*as + 1]*a[2*as + 2] - a[1*as + 2]*a[2*as + 1])
         + a[0*as + 1]*(a[1*as + 2]*a[2*as + 0] - a[1*as + 0]*a[2*as + 2])
         + a[0*as + 2]*(a[1*as + 0]*a[2*as + 1] - a[1*as + 1]*a[2*as + 0]);
  } // determinant

  template<typename real_t>
  real_t invert3x3(real_t inv[], int const is, real_t const a[], int const as=3) {
    auto const det = determinant3x3(a, as);
    if (std::abs(det) < 1e-16) return 0;
    auto const inv_det = 1/det;
    for(int i = 0; i < 3; ++i) {      int const i1 = (i + 1)%3, i2 = (i + 2)%3;
        for(int j = 0; j < 3; ++j) {  int const j1 = (j + 1)%3, j2 = (j + 2)%3;
//             printf("# i=%d i1=%d i2=%d j=%d j1=%d j2=%d\n", i, i1, i2, j, j1, j2);
            inv[i*is + j] = (a[j1*as + i1]*a[j2*as + i2] - a[j2*as + i1]*a[j1*as + i2])*inv_det;
        } // j
    } // i
    return det;
  } // invert3x3
  
  template<typename real_t>
  inline void matrix_rotation(int const n, real_t c[], int const cs, real_t const a[], int const as, real_t const u[], int const us) {
      // compute C = U^transposed * A * U
      for(int i = 0; i < n; ++i) {
          for(int j = 0; j < n; ++j) {
              real_t c_ij = 0;
              for(int k = 0; k < n; ++k) {
                  real_t t_kj = 0;
                  for(int l = 0; l < n; ++l) {
                      t_kj += a[k*as + l] * u[j*us + l];
                  } // l
                  c_ij += u[i*us + k] * t_kj;
              } // k
              c[i*cs + j] = c_ij;
          } // j
      } // i
  } // matrix_rotation


#ifdef  NO_UNIT_TESTS
  inline status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template<typename real_t>
  inline double matmul(real_t c[], int const n, real_t const a[], real_t const b[]) {
      double dev = 0;
      for(int i = 0; i < n; ++i) {
          for(int j = 0; j < n; ++j) {
              real_t cc = 0;
              for(int k = 0; k < n; ++k) {
                  cc += a[i*n + k] * b[k*n + j];
              } // k
              c[i*n + j] = cc;
              dev += std::abs(cc - (i == j));
          } // j
      } // i
      return dev; // deviation
  } // matmul

  inline status_t test_inversion2(int const echo=3) {
    printf("\n# %s\n", __func__);
    int const n = 2;
    double a[n*n], inv[n*n], ainv[n*n], inva[n*n];
    for(int ij = 0; ij < n*n; ++ij) a[ij] = random<double>(-1, 1);
    auto const det = invert2x2(inv, n, a);
    if (0 == det) return -1;
    auto const dev_ia = matmul(inva, n, inv, a);
    auto const dev_ai = matmul(ainv, n, a, inv);
    printf("\n# %s deviations %g and %g\n", __func__, dev_ia, dev_ai);
    if (echo > 5)
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j) 
          printf("# i=%d j=%d a*inv_ij=%.6f inv*a_ij=%.6f a_ij=%g inv_ij=%g\n", i, j, ainv[i*n + j], inva[i*n + j], a[i*n + j], inv[i*n + j]);
    return (dev_ia + dev_ai > 1e-14);
  } // test inversion

  inline status_t test_inversion3(int const echo=3) {
    printf("\n# %s\n", __func__);
    int const n = 3;
    double a[n*n], inv[n*n], ainv[n*n], inva[n*n];
    for(int ij = 0; ij < n*n; ++ij) a[ij] = random<double>(-1, 1);
    auto const det = invert3x3(inv, n, a);
    if (0 == det) return -1;
    auto const dev_ia = matmul(inva, n, inv, a);
    auto const dev_ai = matmul(ainv, n, a, inv);
    printf("\n# %s deviations %g and %g\n", __func__, dev_ia, dev_ai);
    if (echo > 5)
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n; ++j) 
              printf("# i=%d j=%d a*inv_ij=%.6f inv*a_ij=%.6f a_ij=%g inv_ij=%g\n", i, j, ainv[i*n + j], inva[i*n + j], a[i*n + j], inv[i*n + j]);
    return (dev_ia + dev_ai > 1e-14);
  } // test inversion
  
  inline status_t all_tests(int const echo=3) {
    if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
    auto status = 0;
    status += std::abs(test_inversion2());
    status += std::abs(test_inversion3());
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace simple_math
