#pragma once

#ifndef NO_UNIT_TESTS
    #include <cstdio> // printf
    #include <cmath> // std::pow (for reference)
    #include <algorithm> // std::max
#endif // NO_UNIT_TESTS


  template<typename T> inline int constexpr sgn(T const val) {
      return (T(0) < val) - (val < T(0));
  } // sgn

  template<typename T> inline T constexpr pow2(T const x) { return x*x; } // pow2
  template<typename T> inline T constexpr pow3(T const x) { return x*pow2(x); } // pow3

  template<typename real_t> 
  inline real_t intpow(real_t const x, unsigned n) {
      // power function using recursive doubling, only non-negative powers possible
      real_t xbin = x, xpow = (real_t)1;
      while (n) {
          if (n & 1) xpow *= xbin; // if n modulo 2 == 1
          n >>= 1; // divide n by 2
          xbin *= xbin; // square x
      } // while n is nonzero
      return xpow;
  } // intpow


  template<typename real_t> 
  inline real_t constexpr factorial(unsigned const n) {
      return (n > 1)? factorial<real_t>(n - 1)*((real_t)n) : 1;
  } // factorial


  template<typename real_t>
  void set(real_t y[], int const n, real_t const a) {
      for(int i = 0; i < n; ++i) y[i] = a;
  } // set

  template<typename real_t>
  void set(real_t y[], int const n, real_t const a[], real_t const f=1) {
      for(int i = 0; i < n; ++i) y[i] = a[i]*f;
  } // set

  template<typename real_t>
  void scale(real_t y[], int const n, real_t const a[], real_t const f=1) {
      for(int i = 0; i < n; ++i) y[i] *= a[i]*f;
  } // scale

  template<typename real_t>
  void scale(real_t y[], int const n, real_t const f) {
      // no default value for f here, as scaling with 1.0 has no effect
      for(int i = 0; i < n; ++i) y[i] *= f;
  } // scale

  template<typename real_t>
  void product(real_t y[], int const n, real_t const a[], real_t const b[], real_t const f=1) {
      for(int i = 0; i < n; ++i) y[i] = a[i]*b[i]*f;
  } // product

  template<typename real_t>
  void product(real_t y[], int const n, real_t const a[], real_t const b[], real_t const c[], real_t const f=1) {
      for(int i = 0; i < n; ++i) y[i] = a[i]*b[i]*c[i]*f;
  } // product

  template<typename real_t>
  void add_product(real_t y[], int const n, real_t const a[], real_t const f) {
      // no default value for f here, otherwise the name add_product is missleading!
      for(int i = 0; i < n; ++i) y[i] += a[i]*f;
  } // add_product == axpy-type

  template<typename real_t>
  void add_product(real_t y[], int const n, real_t const a[], real_t const b[], real_t const f=1) {
      for(int i = 0; i < n; ++i) y[i] += a[i]*b[i]*f;
  } // add_product
  
namespace inline_math {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
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

  inline status_t all_tests(int const echo=3) {
    if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
    auto status = 0;
    status += test_intpow<float>(echo, 6e-6);
    status += test_intpow<double>(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace inline_math
