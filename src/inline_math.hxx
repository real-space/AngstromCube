#pragma once


  #include "status.hxx" // status_t


  template<typename T> inline int constexpr sgn(T const val) {
      return (T(0) < val) - (val < T(0));
  } // sgn

  template<typename T> inline T constexpr pow2(T const x) { return x*x; } // pow2
  template<typename T> inline T constexpr pow3(T const x) { return x*pow2(x); } // pow3
  template<typename T> inline T constexpr pow4(T const x) { return pow2(pow2(x)); }
  template<typename T> inline T constexpr pow8(T const x) { return pow2(pow4(x)); }

  template<typename real_t> inline 
  real_t intpow(real_t const x, unsigned const nexp) {
      // power function using recursive doubling, only non-negative powers possible
      unsigned n = nexp;
      real_t xbin = x, xpow = (real_t)1;
      while (n) {
          if (n & 1) xpow *= xbin; // if n modulo 2 == 1
          n >>= 1; // divide n by 2
          xbin *= xbin; // square x
      } // while n is nonzero
      return xpow;
  } // intpow


  template<unsigned Step=1> inline
  double constexpr factorial(unsigned const n) {
      return (n > 1)? factorial<Step>(n - Step)*((double)n) : 1;
  } // factorial n! for Step=1 and double_factorial n!! for Step=2

  template<typename real_t> inline
  void set(real_t y[], size_t const n, real_t const a) {
      for(size_t i = 0; i < n; ++i) y[i] = a;
  } // set

  template<typename real_t, typename real_a_t, typename real_f_t=real_t> inline
  void set(real_t y[], size_t const n, real_a_t const a[], real_f_t const f=1) {
      for(size_t i = 0; i < n; ++i) y[i] = a[i]*f;
  } // set

  template<typename real_t, typename real_a_t> inline
  void scale(real_t y[], size_t const n, real_a_t const a[], real_t const f=1) {
      for(size_t i = 0; i < n; ++i) y[i] *= a[i]*f;
  } // scale

  template<typename real_t> inline
  void scale(real_t y[], size_t const n, real_t const f) {
      // no default value for f here, as scaling with 1.0 has no effect
      for(size_t i = 0; i < n; ++i) y[i] *= f;
  } // scale

  template<typename real_t, typename real_a_t, typename real_b_t> inline
  void product(real_t y[], size_t const n, real_a_t const a[], real_b_t const b[], real_t const f=1) {
      for(size_t i = 0; i < n; ++i) y[i] = a[i]*b[i]*f;
  } // product

  template<typename real_t, typename real_a_t, typename real_b_t, typename real_c_t> inline
  void product(real_t y[], size_t const n, real_a_t const a[], real_b_t const b[], real_c_t const c[], real_t const f=1) {
      for(size_t i = 0; i < n; ++i) y[i] = a[i]*b[i]*c[i]*f;
  } // product

  template<typename real_t, typename real_a_t> inline
  void add_product(real_t y[], size_t const n, real_a_t const a[], real_t const f) {
      // no default value for f here, otherwise the name add_product is missleading!
      for(size_t i = 0; i < n; ++i) y[i] += a[i]*f;
  } // add_product == axpy-type

  template<typename real_t, typename real_a_t, typename real_b_t> inline
  void add_product(real_t y[], size_t const n, real_a_t const a[], real_b_t const b[], real_t const f=1) {
      for(size_t i = 0; i < n; ++i) y[i] += a[i]*b[i]*f;
  } // add_product

  template<typename real_t, typename real_a_t> inline
  double dot_product(size_t const n, real_t const bra[], real_a_t const ket[]) {
      double dot = 0;
      for(size_t i = 0; i < n; ++i) {
          dot += bra[i]*ket[i];
      } // i
      return dot;
  } // dot_product

  template<typename real_t, typename real_a_t, typename real_b_t> inline
  double dot_product(size_t const n, real_t const bra[], real_a_t const ket[], real_b_t const metric[]) {
      double dot = 0;
      for(size_t i = 0; i < n; ++i) {
          dot += bra[i]*metric[i]*ket[i];
      } // i
      return dot;
  } // dot_product
  
  
    #include <cstdio> // printf
#ifndef NO_UNIT_TESTS
    #include <cassert> // assert
    #include <cmath> // std::round, std::pow (for reference)
    #include <algorithm> // std::max
   
   inline bool is_integer(double const f) { return (f == std::round(f)); }

#endif // NO_UNIT_TESTS

namespace inline_math {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template<typename real_t>
  inline status_t test_intpow(int const echo=4, double const threshold=9e-14) {
      if (echo > 2) printf("\n# %s %s \n", __FILE__, __func__);
      double max_all_dev = 0;
      for(int ipower = 0; ipower < 127; ++ipower) { // 127 for float, up to 1023 for double
          double max_dev = 0;
          for(real_t x = 0.5; x < 2.0; x *= 1.1) {
              // scan the relevant mantissa and exponent range of floating point numbers
              //   in 2*23596 (0.1% increment) or 2*2380 (1.0%) or 2*272 (10%) steps
              auto const ref = std::pow(x, ipower);
              auto const value = intpow(x, ipower);
              auto const rel_dev = std::abs(value - ref)/ref;
              if (echo > 8) printf("# %s: relative deviation for %.6f^%d is %.1e\n", __func__, x, ipower, rel_dev);
              max_dev = std::max(max_dev, rel_dev);
          } // x
          if (echo > 5) printf("# %s: max. relative deviation for x^%d is %.1e\n", __func__, ipower, max_dev);
          max_all_dev = std::max(max_all_dev, max_dev);
      } // ipower
      if (echo > 1) printf("\n# %s: max. relative deviation for <real_%ld>^<int> is %.1e\n", __func__, sizeof(real_t), max_all_dev);
      return (max_all_dev > threshold);
  } // test_intpow

  inline status_t test_factorials(int const echo=4, double const threshold=9e-14) {
      if (echo > 2) printf("\n# %s %s \n", __FILE__, __func__);
      status_t stat = 0;
      double fac = 1;
      for(int n = 0; n < 171; ++n) {
          assert( factorial(n) == fac );
//        printf("# %i! = %g\n", n, fac);
          if (!is_integer(fac)) { ++stat; if (echo > 0) printf("# %i! = %g is non-integer by %g\n", n, fac, fac - std::round(fac)); }
          fac *= (n + 1); // prepare next
      } // n
      double dfac[] = {1, 1}; // {even, odd}
      for(int n = 0; n < 301; ++n) {
          double & fac2 = dfac[n & 1];
          assert( factorial<2>(n) == fac2 );
//        printf("# %i!! = %g\n", n, fac2);
          if (!is_integer(fac2)) { ++stat; if (echo > 0) printf("# %i! = %g is non-integer by %g\n", n, fac2, fac2 - std::round(fac2)); }
          fac2 *= (n + 2); // prepare next
      } // n
      // found: exponent of double numbers exceeds limits before the mantissa produces non-integer representations    
      return stat;
  } // test_factorials

  inline status_t all_tests(int const echo=3) {
    if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
    status_t status(0);
    status += test_intpow<float>(echo, 6e-6);
    status += test_intpow<double>(echo);
    status += test_factorials(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace inline_math
