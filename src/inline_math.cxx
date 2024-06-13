// This file is part of AngstromCube under MIT License

#ifndef   NO_UNIT_TESTS

    #include <cstdio> // std::printf
    #include <cassert> // assert
    #include <algorithm> // std::max
    #include <cmath> // std::pow

    #include "inline_math.hxx"

#endif // NO_UNIT_TESTS

    #include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

namespace inline_math {

#ifdef    NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

  template <typename real_t>
  status_t test_intpow(int const echo=4, double const threshold=9e-14) {
      if (echo > 3) std::printf("\n# %s %s \n", __FILE__, __func__);
      double max_all_dev{0};
      for (int ipower = 0; ipower < 127; ++ipower) { // 127 for float, up to 1023 for double
          double max_dev{0};
          for (real_t x = 0.5; x < 2.0; x *= 1.1) {
              // scan the relevant mantissa and exponent range of floating point numbers
              //   in 2*23596 (0.1% increment) or 2*2380 (1.0%) or 2*272 (10%) steps
              auto const ref = std::pow(x, ipower);
              auto const value = intpow(x, ipower);
              auto const rel_dev = std::abs(value - ref)/ref;
              if (echo > 8) std::printf("# %s: relative deviation for %.6f^%d is %.1e\n",
                                            __func__,                x, ipower, rel_dev);
              max_dev = std::max(max_dev, rel_dev);
          } // x
          if (echo > 5) std::printf("# %s: max. relative deviation for x^%d is %.1e\n",
                                        __func__,                     ipower, max_dev);
          max_all_dev = std::max(max_all_dev, max_dev);
      } // ipower
      if (echo > 3) std::printf("\n# %s: max. relative deviation for <%s>^<int> is %.1e\n",
                            __func__, (4 == sizeof(real_t))?"float":"double", max_all_dev);
      return (max_all_dev > threshold);
  } // test_intpow

  status_t test_factorials(int const echo=4, double const threshold=9e-14) {
      if (echo > 3) std::printf("\n# %s %s \n", __FILE__, __func__);
      status_t stat(0);
      double fac{1};
      for (int n = 0; n < 29; ++n) { // this can go up to 171(-O0, -O1, -O2), 29(-Ofast)
          if (factorial(n) != fac) std::printf("# factorial(%d) deviates\n", n);
          assert( factorial(n) == fac );
//        std::printf("# %i! = %g\n", n, fac);
          if (!is_integer(fac)) {
              ++stat;
              if (echo > 0) std::printf("# %i! = %g is non-integer by %g\n",
                                       n, fac, fac - std::round(fac));
          } // is not integer
          fac *= (n + 1.); // prepare next
      } // n
      double dfac[] = {1, 1}; // {even, odd}
      for (int n = 0; n < 31; ++n) { // this can go up to 301(-O0), 31(-Ofast)
          double & fac2 = dfac[n & 0x1];
          if (factorial<2>(n) != fac2 && echo > 0) std::printf("# double factorial(%d) deviates\n", n);
          assert( factorial<2>(n) == fac2 );
//        std::printf("# %i!! = %g\n", n, fac2);
          if (!is_integer(fac2)) { 
              ++stat;
              if (echo > 0) std::printf("# %i! = %g is non-integer by %g\n", 
                                       n, fac2, fac2 - std::round(fac2));
          } // is not integer
          fac2 *= (n + 2); // prepare next
      } // n
      // found: exponent of double numbers exceeds limits before the mantissa produces non-integer representations
      if (stat + echo > 3) std::printf("# %s(threshold= %.1e): %d errors\n", __func__, threshold, int(stat)); 
      return stat;
  } // test_factorials

  status_t test_align(int const echo=1) {
      status_t stat(0);
      for (int i = 0; i < 99; ++i) {
          stat += (align<0>(i) != i); // error if align<0> does not reproduce itself
      } // i
      for (int i = 1; i < (1 << 30); i *= 2) {
          if (echo > 15) std::printf("# align<%d>(%d) = %ld\n", 1, i, align<1>(i));
          stat += (align<0>(i) != i);
          stat += (align<1>(i) != i && i > 1);
          stat += (align<2>(i) != i && i > 2);
          stat += (align<3>(i) != i && i > 4);
      } // i
      stat += (align<1>(1) != 2);
      stat += (align<1>(3) != 4);
      stat += (align<1>(7) != 8);
      stat += (align<2>(1) != 4);
      stat += (align<2>(3) != 4);
      stat += (align<2>(7) != 8);
      stat += (align<3>(1) != 8);
      stat += (align<3>(3) != 8);
      stat += (align<3>(7) != 8);
      if (stat + echo > 3) std::printf("# %s: %d errors\n", __func__, int(stat)); 
      return stat;
  } // test_align

  status_t test_sgn(int const echo=0) {
      status_t stat(0); // double             // float              // int
      stat += (sgn(-1.10) != -1.0); stat += (sgn(-1.1f) != -1.f); stat += (sgn(-1) != -1);
      stat += (sgn( 1.10) !=  1.0); stat += (sgn( 1.1f) !=  1.f); stat += (sgn( 1) !=  1);
      stat += (sgn( 0.00) !=  0.0); stat += (sgn( 0.f)  !=  0.f); stat += (sgn( 0) !=  0);
      if (stat + echo > 3) std::printf("# %s<double float int>: %d errors\n", __func__, int(stat)); 
      return stat;
  } // test_sgn

  status_t test_is_integer(int const echo=0) {
      status_t stat(0); // double             // float
      stat += (is_integer(-1.5) != false); stat += (is_integer(-1.5f) != false);
      stat += (is_integer( 1.0) != true ); stat += (is_integer( 1.0f) != true );
      if (stat + echo > 3) std::printf("# %s<double float>: %d errors\n", __func__, int(stat)); 
      return stat;
  } // test_is_integer

  template <typename real_t, size_t N>
  status_t test_scale(int const echo=0) {
      status_t stat(0);
      real_t y[N], a[N], yf[N], yaf[N], f{3};
      for (size_t i{0}; i < N; ++i) {
          y[i] = i - N/real_t(2);
          a[i] = (i + 1/real_t(4))*i;
          yf[i] = y[i];
          yaf[i] = y[i];
      } // i
      scale(yaf, N, a, f); // y[:] *= a[:]*f
      scale(yf, N, f);     // y[:] *= f
      for (size_t i{0}; i < N; ++i) {
          stat += (y[i]*a[i]*f != yaf[i]);
          stat += (y[i]*f      != yf[i]);
      } // i
      if (stat + echo > 3) std::printf("# %s(N= %ld): %d errors\n", __func__, N, int(stat)); 
      return stat;
  } // test_scale

  template <typename real_t, size_t N>
  status_t test_set(int const echo=0) {
      status_t stat(0);
      real_t a[N], ya[N], yaf[N], f{3};
      for (size_t i{0}; i < N; ++i) {
          a[i] = i - N/real_t(2);
      } // i
      set(ya, N, a); // y[:] = a[:]
      set(yaf, N, a, f); // y[:] = a[:]*f
      for (size_t i{0}; i < N; ++i) {
          stat += (a[i] != ya[i]);
          stat += (a[i]*f != yaf[i]);
      } // i
      if (stat + echo > 3) std::printf("# %s(N= %ld): %d errors\n", __func__, N, int(stat)); 
      return stat;
  } // test_set

  template <typename real_t, size_t N>
  status_t test_product(int const echo=0) {
      status_t stat(0);
      real_t a[N], b[N], c[N], f{3}, yabf[N], yabcf[N];
      for (size_t i{0}; i < N; ++i) {
          a[i] = i - N/real_t(2);
          b[i] = (i + 1/real_t(4))*i;
          c[i] = i*i + 9*i + 99;
      } // i
      product(yabcf, N, a, b, c, f); // y[:] = a[:]*b[:]*c[:]*f
      product(yabf, N, a, b, f);     // y[:] = a[:]*b[:]*f
      for (size_t i{0}; i < N; ++i) {
          stat += (a[i]*b[i]*c[i]*f != yabcf[i]);
          stat += (a[i]*b[i]*f      != yabf[i]);
      } // i
      if (stat + echo > 3) std::printf("# %s(N= %ld): %d errors\n", __func__, N, int(stat)); 
      return stat;
  } // test_product

  template <typename real_t, size_t N>
  status_t test_add_product(int const echo=0) {
      status_t stat(0);
      real_t a[N], b[N], y[N], f{3}, yaf[N], yabf[N];
      for (size_t i{0}; i < N; ++i) {
          a[i] = i - N/real_t(2);
          b[i] = (i + 1/real_t(4))*i;
          y[i] = i*i + 9*i + 99;
          yaf[i] = y[i]; yabf[i] = y[i];
      } // i
      add_product(yabf, N, a, b, f); // y[:] += a[:]*b[:]*f
      add_product(yaf, N, a, f);     // y[:] += a[:]*f
      for (size_t i{0}; i < N; ++i) {
          stat += (y[i] + a[i]*b[i]*f != yabf[i]);
          stat += (y[i] + a[i]*f      != yaf[i]);
      } // i
      if (stat + echo > 3) std::printf("# %s(N= %ld): %d errors\n", __func__, N, int(stat)); 
      return stat;
  } // test_add_product

  template <typename real_t, size_t N>
  status_t test_dot_product(int const echo=0) {
      status_t stat(0);
      real_t a[N], b[N], m[N];
      for (size_t i{0}; i < N; ++i) {
          a[i] = i + 1;
          b[i] = i + 2;
          m[i] = 1; // metric
      } // i
      auto const dot_ab  = dot_product(N, a, b);
      auto const dot_abm = dot_product(N, a, b, m);
      if (echo > 9) std::printf("# %s(N= %ld): results %g and %g\n", __func__, N, dot_ab*1., dot_abm*1.);
      auto const reference = 20*N*N - 327*N + 1653; 
      stat += (reference != dot_ab);
      stat += (reference != dot_abm);
      if (stat + echo > 3) std::printf("# %s(N= %ld): %d errors\n", __func__, N, int(stat)); 
      return stat;
  } // test_dot_product

  template <typename real_t>
  status_t test_pow234(size_t const n, int const echo=0) {
      status_t stat(0);
      for (size_t i{0}; i < n; ++i) {
          auto const x = i - n/real_t(2);
          auto const x2 = pow2(x), x3 = pow3(x), x4 = pow4(x);
          stat += (x*x != x2);
          stat += (x2*x != x3);
          stat += (x2*x2 != x4);
      } // i
      if (stat + echo > 3) std::printf("# %s(n= %ld): %d errors\n", __func__, n, int(stat)); 
      return stat;
  } // test_pow234

  status_t test_various(int const echo=0) {
      status_t stat(0);
      size_t constexpr N = 29;
      stat += test_pow234<float>(N, echo);
      stat += test_pow234<double>(N, echo);
      stat += test_pow234<int64_t>(N, echo);
      stat += test_scale<float,N>(echo);
      stat += test_scale<double,N>(echo);
      stat += test_scale<int64_t,N>(echo);
      stat += test_set<float,N>(echo);
      stat += test_set<double,N>(echo);
      stat += test_set<int64_t,N>(echo);
      stat += test_product<float,N>(echo);
      stat += test_product<double,N>(echo);
      stat += test_product<int64_t,N>(echo);
      stat += test_add_product<float,N>(echo);
      stat += test_add_product<double,N>(echo);
      stat += test_add_product<int64_t,N>(echo);
      stat += test_dot_product<float,N>(echo);
      stat += test_dot_product<double,N>(echo);
      stat += test_dot_product<int64_t,N>(echo);
      if (stat + echo > 3) std::printf("# %s<float double int64_t>: %d errors\n", __func__, int(stat)); 
      return stat;
  } // test_various

  status_t all_tests(int const echo) {
      if (echo > 3) std::printf("\n# %s %s\n", __FILE__, __func__);
      status_t stat(0);
      stat += test_intpow<float>(echo, 6e-6);
      stat += test_intpow<double>(echo);
      stat += test_align(echo);
      stat += test_factorials(echo);
      stat += test_sgn(echo);
      stat += test_is_integer(echo);
      stat += test_various(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace inline_math
