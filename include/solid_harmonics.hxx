#pragma once

#include <cstdio> // std::printf
#include <cmath> // std::sqrt, std::sin, std::cos
#include <cassert> // assert
#include <vector> // std::vector

#include "constants.hxx" // ::pi
#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

namespace solid_harmonics {

  double constexpr pi = constants::pi;
  double constexpr Y00inv = 3.5449077018110318, // == sqrt(4*pi)
                      Y00 = .28209479177387817; // == 1./Y00inv;

  template <typename real_t>
  void rlXlm_implementation(
        real_t xlm[]
      , int const ellmax
      , real_t const cth, real_t const sth
      , real_t const cph, real_t const sph
      , real_t const r2=1
      , bool const cos_trick=true
  ) {
// !************************************************************
// !     generate the spherical harmonics for the vector v
// !     using a stable upward recursion in l.  (see notes
// !     by m. weinert.)
// !          m.weinert   january 1982
// !     modified by R. Podloucky (added in xnorm); July 1989
// !     cleaned up    mw 1995
// !************************************************************

// check whether  or not normalizations are needed
      static std::vector<real_t> xnorm;
      static int ellmaxd = -1; // -1:not_initalized

      if (ellmax > ellmaxd) {
#ifdef DEBUG
          if (xnorm.size() > 0) std::printf("# %s resize table of normalization constants from %d to %d\n", __func__, (1 + ellmaxd)*(1 + ellmaxd), (1 + ellmax)*(1 + ellmax));
#endif
          xnorm.resize((1 + ellmax)*(1 + ellmax));

// !********************************************************************
// !     normalization constants for ylm (internal subroutine has access
// !     to ellmax and xnorm from above)
// !********************************************************************
          { // scope to fill xnorm with values
              double const fpi = 4.0*pi;
              for (int l = 0; l <= ellmax; ++l) {
                  int const lm0 = l*l + l;
                  double const a = std::sqrt((2*l + 1.)/fpi);
                  double cd{1};
                  xnorm[lm0] = a;
                  double sgn{-1};
                  for (int m = 1; m <= l; ++m) {
                      cd /= ((l + 1. - m)*(l + m));
                      auto const xn = a*std::sqrt(2*cd); // different from Ylm!
                      xnorm[lm0 + m] = xn;
                      xnorm[lm0 - m] = xn * sgn; // not used here
                      sgn = -sgn; // prepare (-1)^m for the next iteration
                  } // m
              } // l
          } // scope
          ellmaxd = ellmax; // set static variable
      } else if (ellmax < 0) {
          xnorm.resize(0); // cleanup
      }
      if (ellmax < 0) return;

      int const S = (1 + ellmax); // stride for p, the array of associated Legendre functions
      auto const p = new real_t[S*S]; // associated Legendre functions

      // generate associated Legendre functions for m >= 0
      real_t fac{1};
      // loop over m values
      for (int m = 0; m < ellmax; ++m) {
          fac *= (1 - 2*m);
          p[m     + S*m] = fac;
          p[m + 1 + S*m] = (m + 1 + m)*cth*fac;
          // recurse upward in l
          for (int l = m + 2; l <= ellmax; ++l) {
              p[l + S*m] = ((2*l - 1)*cth*p[l - 1 + S*m] - (l + m - 1)*r2*p[l - 2 + S*m])/real_t(l - m);
          } // l
          fac *= sth;
      } // m
      p[ellmax + S*ellmax] = (1 - 2*ellmax)*fac;

      auto const cs = new real_t[S];
      auto const sn = new real_t[S];
      // determine cos(m*phi) and sin(m*phi)
      cs[0] = 1;
      sn[0] = 0;
      if (cos_trick) {
          if (ellmax > 0) {
              cs[1] = cph; sn[1] = sph;
              auto const cph2 = 2*cph;
              for (int m = 2; m <= ellmax; ++m) {
                  // this two-step recursion formula is more accurate but does not work for r^l*X_lm
                  sn[m] = cph2*sn[m - 1] - sn[m - 2];
                  cs[m] = cph2*cs[m - 1] - cs[m - 2];
              } // m
          } // ellmax > 0
      } else {
          for (int m = 1; m <= ellmax; ++m) {
              // addition theorem for sine and cosine
              sn[m] = cph*sn[m - 1] + cs[m - 1]*sph;
              cs[m] = cph*cs[m - 1] - sn[m - 1]*sph;
          } // m
      } // cos_trick

      // multiply in the normalization factors
      for (int m = 0; m <= ellmax; ++m) {
          for (int l = m; l <= ellmax; ++l) {
              int const lm0 = l*l + l;
              real_t const Plm = p[l + S*m];
              xlm[lm0 + m] = Plm*xnorm[lm0 + m]*sn[m]; // imag part
              xlm[lm0 - m] = Plm*xnorm[lm0 + m]*cs[m]; // real part
          } // l
      } // m
      delete[] sn;
      delete[] cs;
      delete[] p;

  } // Xlm_implementation

  template <typename real_t>
  void Xlm(real_t xlm[], int const ellmax, double const theta, double const phi) {
      rlXlm_implementation(xlm, ellmax, std::cos(theta), std::sin(theta), std::cos(phi), std::sin(phi));
  } // Xlm

  template <typename real_t, typename vector_real_t>
  void Xlm(real_t xlm[], int const ellmax, vector_real_t const v[3]) {
      real_t constexpr small = 1e-12;
      auto const x = v[0], y = v[1], z = v[2];
      auto const xy2 = x*x + y*y;
      auto const r = std::sqrt(xy2 + z*z);
      auto const rxy = std::sqrt(xy2);

      // calculate sin and cos of theta and phi
      real_t cth, sth;
      if (r > small) {
         cth = z/r;
         sth = rxy/r;
      } else {
         cth = 1;
         sth = 0;
      }
      real_t cph, sph;
      if (rxy > small) {
         cph = x/rxy;
         sph = y/rxy;
      } else {
         cph = 1;
         sph = 0;
      }
      rlXlm_implementation(xlm, ellmax, cth, sth, cph, sph);
  } // Xlm

  template <typename real_t, typename vector_real_t>
  void rlXlm(real_t xlm[], int const ellmax, vector_real_t const v[3]) {
      real_t const x = v[0], y = v[1], z = v[2], r2 = x*x + y*y + z*z;
      rlXlm_implementation(xlm, ellmax, z, 1., x, y, r2, false); // ToDo: check: maybe something is still missing in Xlm_implementation
  } // rlXlm

  template <typename real_t=double>
  void cleanup() { real_t z{0}; rlXlm_implementation(&z, -1, z, z, z, z); } // free internal memory

  inline int find_ell(int const lm) { int lp1{0}; while (lp1*lp1 <= lm) ++lp1; return lp1 - 1; }
  inline int find_emm(int const lm, int const ell) { return lm - (ell*ell + ell); }
  inline int find_emm(int const lm) { return find_emm(lm, find_ell(lm)); }
  inline int lm_index(int const ell, int const emm) { return ell*ell + ell + emm; }

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_indices(int const echo=0) { // test interal consistency of find_-functions
      status_t stat(0);
      for (int lm = -3; lm < 64; ++lm) {
          int const ell = find_ell(lm),
                    emm = find_emm(lm, ell);
          if (echo > 4) std::printf("# %s    lm=%d -> ell=%d emm=%d\n", __FILE__, lm, ell, emm);
          stat += (lm_index(ell, emm) != lm);
      } // lm
      return stat;
  } // test_indices

  inline status_t test_Y00_inverse(int const echo=0) {
      return ( Y00 * Y00inv != 1.0 ); // should be exact
  } // test_Y00_inverse

  inline status_t all_tests(int const echo=0) {
      if (echo > 0) std::printf("\n# %s %s\n", __FILE__, __func__);
      status_t stat(0);
      stat += test_Y00_inverse(echo);
      stat += test_indices(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace solid_harmonics
