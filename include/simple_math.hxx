#pragma once


#include <cstdio> // std::printf
#include <random> // std::random_device, std::mt19937, std::is_integral, std::conditional, std::uniform_int_distribution, std::uniform_real_distribution
#include <complex> // std::complex<real_t>

#include "status.hxx" // status_t

namespace simple_math {

  // templated, performant, thread safe pseudo random number generator based on Mersenne twister
  template <typename real_t, typename Generator=std::mt19937>
  inline real_t random(real_t const from, real_t const to) {

      thread_local static Generator gen(std::random_device{}());

      using dist_type = typename std::conditional
      < std::is_integral<real_t>::value
        , std::uniform_int_distribution<real_t>
        , std::uniform_real_distribution<real_t>
      >::type;

      thread_local static dist_type dist;

      return dist(gen, typename dist_type::param_type{from, to});
  } // random


  // a00 a01
  // a10 a11

  template <typename T>
  inline T determinant(T a00, T a01, T a10, T a11) {
    return a00*a11 - a10*a01; 
  } // determinant
  
  // a00 a01 a02
  // a10 a11 a12
  // a20 a21 a22

  template <typename T>
  inline T determinant(T a00, T a01, T a02, T a10, T a11, T a12, T a20, T a21, T a22) {
    return a00*determinant(a11, a12, a21, a22) -
    	   a10*determinant(a01, a02, a21, a22) +
    	   a20*determinant(a01, a02, a11, a12);
  } // determinant

  // a00 a01 a02 a03
  // a10 a11 a12 a13
  // a20 a21 a22 a23
  // a30 a31 a32 a33

  template <typename T>
  inline T determinant(T a00, T a01, T a02, T a03, T a10, T a11, T a12, T a13, 
  	                   T a20, T a21, T a22, T a23, T a30, T a31, T a32, T a33) {
    return a00*determinant(a11, a12, a13, a21, a22, a23, a31, a32, a33) -
    	   a10*determinant(a01, a02, a03, a21, a22, a23, a31, a32, a33) +
    	   a20*determinant(a01, a02, a03, a11, a12, a13, a31, a32, a33) -
    	   a30*determinant(a01, a02, a03, a11, a12, a13, a21, a22, a23);
  } // determinant

  //       a11 a12 a13
  // a00*  a21 a22 a23
  //       a31 a32 a33

  //       a01 a02 a03
  // a10*  a21 a22 a23
  //       a31 a32 a33

  //       a01 a02 a03
  // a20*  a11 a12 a13
  //       a31 a32 a33

  //       a01 a02 a03
  // a30*  a11 a12 a13
  //       a21 a22 a23

  template <typename T>
  inline T determinant(int const n, T const a[], int const as) {

  	int constexpr i0 = 0, i1 = 1, i2 = 2, i3 = 3;
  	int constexpr j0 = 0, j1 = 1, j2 = 2, j3 = 3; 

  	if (1 == n) return  a[i0*as + j0];
  	if (2 == n) return determinant(
  						a[i0*as + j0], a[i0*as + j1], 
  			            a[i1*as + j0], a[i1*as + j1]);
  	if (3 == n) return determinant(
  						a[i0*as + j0], a[i0*as + j1], a[i0*as + j2],
                		a[i1*as + j0], a[i1*as + j1], a[i1*as + j2],
  			    		a[i2*as + j0], a[i2*as + j1], a[i2*as + j2]);
  	if (4 == n) return determinant(
  						a[i0*as + j0], a[i0*as + j1], a[i0*as + j2], a[i0*as + j3],
            			a[i1*as + j0], a[i1*as + j1], a[i1*as + j2], a[i1*as + j3],
  		    			a[i2*as + j0], a[i2*as + j1], a[i2*as + j2], a[i2*as + j3],
  						a[i3*as + j0], a[i3*as + j1], a[i3*as + j2], a[i3*as + j3]);
    return (n > 4) - (n < 1); // not implemented (n > 4) or impossible (n < 1)
  } // determinant

  float constexpr threshold = 1e-30;
  
  template <typename real_t>
  inline real_t invert1x1(real_t inv[], int const is, real_t const a[], int const as=1) {
    auto const det = determinant(1, a, as);
    if (std::abs(det) < threshold) return 0;
    auto const inv_det = real_t(1)/det;
    inv[0*is + 0] = inv_det;
    return det;
  } // invert

  template <typename real_t>
  inline real_t invert2x2(real_t inv[], int const is, real_t const a[], int const as=2) {
    auto const det = determinant(2, a, as);
    if (std::abs(det) < threshold) return 0;
    auto const inv_det = real_t(1)/det;
    inv[0*is + 0] =   a[1*as + 1]*inv_det;
    inv[0*is + 1] = - a[0*as + 1]*inv_det;
    inv[1*is + 0] = - a[1*as + 0]*inv_det;
    inv[1*is + 1] =   a[0*as + 0]*inv_det;
    return det;
  } // invert

  template <typename real_t>
  inline real_t invert3x3(real_t inv[], int const is, real_t const a[], int const as=3) {
    auto const det = determinant(3, a, as);
    if (std::abs(det) < threshold) return 0;
    auto const inv_det = real_t(1)/det;
    for (int i = 0; i < 3; ++i) {      int const i1 = (i + 1)%3, i2 = (i + 2)%3;
        for (int j = 0; j < 3; ++j) {  int const j1 = (j + 1)%3, j2 = (j + 2)%3;
//             std::printf("# i=%d i1=%d i2=%d j=%d j1=%d j2=%d\n", i, i1, i2, j, j1, j2);
            inv[i*is + j] = (a[j1*as + i1]*a[j2*as + i2] - a[j2*as + i1]*a[j1*as + i2])*inv_det;
        } // j
    } // i
    return det;
  } // invert

  template <typename real_t>
  inline real_t determinant4x4(real_t const a[], int const as=4) {
    return determinant(
    	 a[0*as+0], a[0*as+1], a[0*as+2], a[0*as+3],
         a[1*as+0], a[1*as+1], a[1*as+2], a[1*as+3],
    	 a[2*as+0], a[2*as+1], a[2*as+2], a[2*as+3],
    	 a[3*as+0], a[3*as+1], a[3*as+2], a[3*as+3]);
  } // determinant


  template <typename real_t>
  inline real_t invert4x4(real_t inv[], int const is, real_t const a[], int const as=4) {
    auto const det = determinant(4, a, as);
    if (std::abs(det) < threshold) return 0;
    auto const inv_det = real_t(1)/det;
    for (int j = 0; j < 4; ++j) {
	    int const j0 = 0 + (0 >= j), j1 = 1 + (1 >= j), j2 = 2 + (2 >= j);
	    for (int i = 0; i < 4; ++i) {
    		int const i0 = 0 + (0 >= i), i1 = 1 + (1 >= i), i2 = 2 + (2 >= i);
	    	real_t const sign = 1 - 2*( (i + j) & 1 );
            inv[j*is + i] = sign * inv_det * determinant(
		    	 a[i0*as + j0], a[i0*as + j1], a[i0*as + j2],
        		 a[i1*as + j0], a[i1*as + j1], a[i1*as + j2],
    	 		 a[i2*as + j0], a[i2*as + j1], a[i2*as + j2]);
        } // i
    } // j
    return det;
  } // invert


  template <typename T>
  inline T invert(int const n, T inv[], int const is, T const a[], int const as) {
  	  if (1 == n) return invert1x1(inv, is, a, as);
  	  if (2 == n) return invert2x2(inv, is, a, as);
  	  if (3 == n) return invert3x3(inv, is, a, as);
  	  if (4 == n) return invert4x4(inv, is, a, as);
 	  return 0; // returns zero
  } // invert

  
  template <typename real_t>
  inline void matrix_rotation(int const n, real_t c[], int const cs, real_t const a[], int const as, real_t const u[], int const us) {
      // compute C = U * A * U^transposed
      for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {
              real_t c_ij = 0;
              for (int k = 0; k < n; ++k) {
                  real_t t_kj = 0;
                  for (int l = 0; l < n; ++l) {
                      t_kj += a[k*as + l] * u[j*us + l];
                  } // l
                  c_ij += u[i*us + k] * t_kj;
              } // k
              c[i*cs + j] = c_ij;
          } // j
      } // i
  } // matrix_rotation


#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <typename real_t>
  inline double matmul(
        int const n // matrix dimension
      , real_t       c[], int const cs // result matrix, stride
      , real_t const a[], int const as // left  input matrix, stride
      , real_t const b[], int const bs // right input matrix, stride
  ) {
      double dev{0};
      for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {
              real_t cc(0);
              for (int k = 0; k < n; ++k) {
                  cc += a[i*as + k] * b[k*bs + j];
              } // k
              c[i*cs + j] = cc;
              dev += (i == j) ? std::abs(cc - real_t(1)) : std::abs(cc);
          } // j
      } // i
      return dev; // deviation
  } // matmul

  template <typename real_t>
  inline status_t test_invert(int const echo=3, double const threshold=1e-12) {
      if (echo > 5) std::printf("\n# %s\n", __func__);
      int constexpr N = 4;
      int constexpr M = 8;
      real_t a[N*M], inv[N*M], ainv[N*M], inva[N*M];
      double maxdev{0};
      for (int n = -1; n <= N; ++n) { // dimension
          for (int m = n; m <= M; ++m) { // stride
              // fill with random values
              for (int ij = 0; ij < n*m; ++ij) {
                  a[ij] = random<real_t>(-1, 1);
              } // ij
              auto const det = invert(n, inv, m, a, m);
              auto const dev_ia = matmul(n, inva, m, inv, m, a, m);
              auto const dev_ai = matmul(n, ainv, m, a, m, inv, m);
              maxdev = std::max(maxdev, std::max(dev_ia, dev_ai));
              if (echo > 5) std::printf("# %s n=%i stride=%i determinant %g deviations %g and %g\n", 
                                      __func__, n, m, det, dev_ia, dev_ai);
              if (echo > 8) {
                  for (int i = 0; i < n; ++i) {
                      for (int j = 0; j < n; ++j) {
                          std::printf("# i=%d j=%d \ta*inv_ij=%.6f \tinv*a_ij=%.6f \ta_ij=%g \tinv_ij=%g\n", 
                                            i,   j, ainv[i*m + j], inva[i*m + j], a[i*m + j], inv[i*m + j]);
                      } // j
                  } // i
              } // echo
          } // m
      } // n
      if (echo > 3) std::printf("# %s: largest deviation for inversions up to %ix%i is %g\n", __func__, N, N, maxdev);
      return (maxdev > threshold);
  } // test_invert

  inline status_t all_tests(int const echo=0) {
      if (echo > 0) std::printf("\n# %s %s\n", __FILE__, __func__);
      status_t stat(0);
      stat += test_invert<double>(echo);
      stat += test_invert<float>(echo, 1e-4);
//    stat += test_invert<std::complex<double>>(echo); // does not compile
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace simple_math
