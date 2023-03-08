#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <complex> // std::complex<real_t>, ::conj
#include <cmath> // std::sqrt
#include <vector> // std::vector<T>

#include "constants.hxx" // ::pi
#include "inline_math.hxx" // pow2

namespace spherical_harmonics {

  template <typename real_t>
  inline void Ylm(
        std::complex<real_t> ylm[]
      , int const ellmax
      , double const v[3]
  ) {
// !************************************************************
// !     generate the spherical harmonics for the vector v
// !     using a stable upward recursion in l.  (see notes
// !     by m. weinert.)
// !          m.weinert   january 1982
// !     modified by R. Podloucky (added in ynorm); July 1989
// !     cleaned up    mw 1995
// !
// !     modified to make use of f90 constructs. note that
// !     the normalization is an internal subroutine and hence
// !     can only be called from here. also, no need to dimension
// !     arrays for ynorm, done dynamically.          mw 1999
// !************************************************************

      real_t constexpr small = 1e-12;


      // check whether  or not normalizations are needed
      static std::vector<real_t> ynorm;
      static int ellmaxd = -1; // -1:not_initalized

      if (ellmax > ellmaxd) {
#ifdef DEBUG
           std::printf("# %s: resize table of normalization constants from %d to %d\n",
              __func__, pow2(1 + ellmaxd), pow2(1 + ellmax));
#endif // DEBUG
          ynorm.resize(pow2(1 + ellmax));

// !********************************************************************
// !     normalization constants for ylm (internal subroutine has access
// !     to ellmax and ynorm from above)
// !********************************************************************
          { // scope to fill ynorm with values
              double const fpi = 4*constants::pi; // 4*pi
              for (int l = 0; l <= ellmax; ++l) {
                  int const lm0 = l*l + l;
                  double const a = std::sqrt((2*l + 1.)/fpi);
                  double cd{1}, sgn{-1};
                  ynorm[lm0] = a;
                  for (int m = 1; m <= l; ++m) {
                      cd /= ((l + 1. - m)*(l + m));
                      auto const yn = a*std::sqrt(cd);
                      ynorm[lm0 + m] = yn;
                      ynorm[lm0 - m] = yn * sgn;
                      sgn = -sgn; // prepare (-1)^m for the next iteration
                  } // m
              } // l
          } // scope
          ellmaxd = ellmax; // set static variable
      } else if (ellmax < 0) {
          ellmaxd = -1; // set static variable
          ynorm.resize(0); // cleanup
      }
      if (ellmax < 0) return;

      // calculate sin and cos of theta and phi
      auto const x = v[0], y = v[1], z = v[2];
      auto const xy2 = x*x + y*y;
      auto const r = std::sqrt(xy2 + z*z);
      auto const rxy = std::sqrt(xy2);

      real_t cth{1}, sth{0};
      if (r > small) {
         cth = z/r;
         sth = rxy/r;
      }
      real_t cph{1}, sph{0};
      if (rxy > small) {
         cph = x/rxy;
         sph = y/rxy;
      }

      int const S = (1 + ellmax); // stride for p, the array of associated legendre functions
      std::vector<real_t> p((1 + ellmax)*S);

      // generate associated legendre functions for m >= 0
      real_t fac{1};
      // loop over m values
      for (int m = 0; m < ellmax; ++m) {
          fac *= (1 - 2*m);
          p[m     + S*m] = fac;
          p[m + 1 + S*m] = (m + 1 + m)*cth*fac;
          // recurse upward in l
          for (int l = m + 2; l <= ellmax; ++l) {
              p[l + S*m] = ((2*l - 1)*cth*p[l - 1 + S*m] - (l + m - 1)*p[l - 2 + S*m])/real_t(l - m);
          } // l
          fac *= sth;
      } // m
      p[ellmax + S*ellmax] = (1 - 2*ellmax)*fac;

      std::vector<real_t> c(1 + ellmax, real_t(1)),
                          s(1 + ellmax, real_t(0));
      // determine sin and cos of phi
      if (ellmax > 0) {
          c[1] = cph; s[1] = sph;
          auto const cph2 = 2*cph;
          for (int m = 2; m <= ellmax; ++m) {
              s[m] = cph2*s[m - 1] - s[m - 2];
              c[m] = cph2*c[m - 1] - c[m - 2];
          } // m
      } // ellmax > 0

      // multiply the normalization factors
      for (int m = 0; m <= ellmax; ++m) {
          for (int l = m; l <= ellmax; ++l) {
              int const lm0 = l*l + l;
              auto const ylms = p[l + S*m]*std::complex<real_t>(c[m], s[m]);
              ylm[lm0 + m] = ynorm[lm0 + m]*ylms;
              ylm[lm0 - m] = ynorm[lm0 - m]*std::conj(ylms);
          } // l
      } // m

      return;
  } // Ylm

  template <typename real_t=double>
  void cleanup(int const echo=0) {
      if (echo > 5) std::printf("# %s %s<%s>: free internal memory\n",
          __FILE__, __func__, (sizeof(real_t) == 4)?"float":"double");
      std::complex<real_t> z{0};
      double v[3];
      Ylm(&z, -1, v);
  } // cleanup









#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <typename real_t>
  inline status_t test_memory_cleanup(int const echo=0, int const ellmax=9) {
      status_t stat(0);
      double const vec[] = {1, 2, 3};
      for (int ell = 0; ell <= ellmax; ++ell) {
          std::vector<std::complex<real_t>> ylm(pow2(1 + ell));
          Ylm(ylm.data(), ell, vec);
          cleanup<real_t>(echo);
      } // ell
      return stat;
  } // test_memory_cleanup

  inline status_t all_tests(int const echo=0) {
      if (echo > 0) std::printf("\n# %s %s\n", __FILE__, __func__);
      status_t stat(0);
      stat += test_memory_cleanup<float>(echo);
      stat += test_memory_cleanup<double>(echo);
      // ToDo: insufficient testing done here. Orthgonality of spherical_harmonics can be checked using angular_grid!
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace spherical_harmonics
