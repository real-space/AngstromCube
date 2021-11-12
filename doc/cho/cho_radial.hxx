#pragma once

#include <cmath> // std::sqrt
#include <cassert> // assert
#include <cstdint> // int8_t
#include <vector> // std::vector<T>
#ifndef NO_UNIT_TESTS
  #include <cstdio> // std::printf
  #include <algorithm> // std::max
#endif

namespace cho_radial {

  /*
   *  Circular Harmonic Oscillator == CHO --- 2-dimensional isotropic harmonics oscillator states
   */

  template <typename int_t>
  inline int_t constexpr nCHO_radial(int_t const numax) { return (numax*(numax + 4) + 4)/4; }

  template <typename real_t>
  void radial_eigenstates(
        real_t poly[] // coefficients of a polynomial in r^2
      , int const nrn // number of radial nodes
      , int const ell // angular momentum quantum number
      , real_t const factor=1 // if we know the normalization prefactor in advance, we can provide it here
  ) {

      // recursion relation of the radial CHO coefficients:
      // a_0 = 1, a_{k+1} = (k-nrn)/((k+1)*(ell+k+3/2)) a_k

      poly[0] = factor;
      for (int k = 0; k < nrn; ++k) {
          poly[k + 1] = (poly[k]*(k - nrn))/real_t((k + 1)*(k + 1 + ell));
      } // k

  } // radial_eigenstates

  double exponential_integral_k(int const k) {
      // I_k = int\limit_0^\infty dr r^(2k+1) exp(-r^2)
      if (0 == k) return 0.5;
      assert(k > 0);
      return k*exponential_integral_k(k - 1); // recursive invokation
  } // exponential_integral_k

  template <typename real_t>
  real_t inner_product(real_t const coeff0[], int const nrn0, real_t const coeff1[], int const nrn1, int const ell) {
      // radial function = sum_p coeff[p] r^(ell + 2p)
      std::vector<double> prod(nrn0 + nrn1 + 1, 0.0); // coefficients of a polynomial in r^2
      for (int p0 = 0; p0 <= nrn0; ++p0) {
          for (int p1 = 0; p1 <= nrn1; ++p1) {
              prod[p0 + p1] += double(coeff0[p0])*coeff1[p1]; // polynomial product with itself
          } // p1
      } // p0
      // product function = sum_p prod[p] r^(2ell + 2p)

      // Overlap integral for circular harmonics: int\limit_0^\infty dr r^(2k+1) exp(-r^2)
      double dot{0};
      double exp_int_k = exponential_integral_k(ell);
      for (int p = 0; p <= nrn0 + nrn1; ++p) { // loop must run serial forward
//        std::printf("# %s ell=%d nrn=%d p=%d\n", __func__, ell, nrn, p);
          assert(exp_int_k == exponential_integral_k(ell + p)); // DEBUG
          dot += prod[p]*exp_int_k;
          exp_int_k *= (ell + p + 1); // prepare exp_int_k for the next iteration
      } // pp
      return dot; // dot product
  } // inner_product

  template <typename real_t>
  real_t radial_normalization(real_t const coeff[], int const nrn, int const ell) {
      auto const norm = inner_product(coeff, nrn, coeff, nrn, ell);
      assert(norm > 0);
      return 1./std::sqrt(norm); // normalization prefactor
  } // radial_normalization

  template <typename real_t=double>
  real_t radial_normalization(int const nrn, int const ell) {
      std::vector<real_t> coeff(nrn + 1);  // coefficients of a polynomial in r^2
      radial_eigenstates(coeff.data(), nrn, ell); // get polynomial coefficients of the CHO eigenstates
      auto const result = radial_normalization(coeff.data(), nrn, ell);
      return result;
  } // radial_normalization

  template <typename real_t>
  real_t expand_poly(real_t const coeff[], int const ncoeff, double const x) {
      real_t value{0};
      double xpow{1};
      for (int p = 0; p < ncoeff; ++p) {
          value += coeff[p] * xpow;
          xpow *= x;
      } // p
      return value;
  } // expand_poly

#ifdef NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  char const *const ellchar = "spdfghijklmno";

  template <typename real_t>
  real_t numerical_inner_product(real_t const c0[], int const nrn0, 
                                 real_t const c1[], int const nrn1, int const ell) {
      double constexpr dr = 1./(1 << 12), rmax = 12.;
      int const nr = rmax / dr;
      real_t dot{0};
      for (int ir = 0; ir < nr; ++ir) {
          double const r = (ir - .5)*dr;
          double const r2 = r*r;
          double const Gauss = std::exp(-r2); // == exp(-0.5*r^2)*exp(-0.5*r^2)
          double const f0 = expand_poly(c0, 1 + nrn0, r2);
          double const f1 = expand_poly(c1, 1 + nrn1, r2);
          double const r2pow = r*std::pow(r2, ell); // == r * (r^ell)^2
          dot += Gauss * r2pow * f0 * f1;
      } // ir
      dot *= dr; // can be taken out of the loop since it is constant
      return dot;
  } // numerical_inner_product

  template <int numax=9>
  inline status_t test_orthonormality(int const echo=1) {
      int constexpr n = nCHO_radial(numax);
      if (echo > 1) std::printf("# %s  numax= %d has %d different radial CHO states\n", __func__, numax, n);

      double c[n][8]; // polynomial coefficients

      int8_t   ell_list[n];
      uint8_t  nrn_list[n];
      double   fac_list[n];

      int i = 0;
//       for (int ene = 0; ene <= numax; ++ene) { //  E_CHO = ene + 1
//           for (int nrn = ene/2; nrn >= 0; --nrn) { // start from low ell withing the block, i.e. many radial nodes
//               int const ell = ene - 2*nrn;
//       for (int ell = 0; ell <= numax; ++ell) { //
      for (int ell = numax; ell >= 0; --ell) { //
          for (int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {

              ell_list[i] = ell;
              nrn_list[i] = nrn;
              radial_eigenstates(c[i], nrn, ell);
              double const fac = radial_normalization(c[i], nrn, ell);
              assert(std::abs(fac - radial_normalization(nrn, ell)) < 1e-12); // check that the other interface produces the same value,with -O0 we can even check for ==
              fac_list[i] = fac;
              if (echo > 2) std::printf("# %s %3d state  nrn= %d  ell= %d  factor= %g\n", __func__, i, nrn, ell, fac);
              radial_eigenstates(c[i], nrn, ell, fac); // overwrite the coefficient series with the normalized onces

              ++i; // count the number of states
          } // nrn
      } // ene
      assert(n == i); // check that nCHO_radial(numax) agrees

      double deviation[][2] = {{0, 0}, {0, 0}};
      for (int method = 0; method < 2; ++method) { // 0:numerical, 1:analytical
          double *const dev = deviation[method];
          for (int i = 0; i < n; ++i) {
              int const ell = ell_list[i];
              for (int j = i; j >= 0; --j) { // triangular loop structure
    //        for (int j = 0; j < n; ++j) {  // block loop structure
                  if (ell == ell_list[j]) {
                      auto const delta_ij = method ? inner_product(c[i], nrn_list[i], c[j], nrn_list[j], ell):
                                           numerical_inner_product(c[i], nrn_list[i], c[j], nrn_list[j], ell);
                      int const i_equals_j = (i == j);
                      if (true) {
                          if (i != j) {
                              auto const delta_ji = method ? inner_product(c[j], nrn_list[j], c[i], nrn_list[i], ell):
                                                   numerical_inner_product(c[j], nrn_list[j], c[i], nrn_list[i], ell);
                              auto const asymmetry = delta_ji - delta_ij;
                              if (std::abs(asymmetry) > 1e-12) {
                                  if (echo > 0) std::printf("# asymmetry of radial CHO eigenfunctions by %.1e between %c%d and %c%d method=%d\n", 
                                                               asymmetry, ellchar[ell_list[i]], nrn_list[i], ellchar[ell_list[j]], nrn_list[j], method);
                                  assert(std::abs(asymmetry) < 1e-9);
                              } // asymmetry large
                          } // i != j
                      } // sanity check
                      if (echo > 4 - i_equals_j) std::printf("%g ", delta_ij - i_equals_j);
                      // measure the deviation from unity for diagonal elements
                      dev[i_equals_j] = std::max(dev[i_equals_j], std::abs(delta_ij - i_equals_j));
                  } // ell matches
              } // j
              if (echo > 3) std::printf("\n");
          } // i
          for (int diag = 0; diag < 2*(echo > 0); ++diag) {
              std::printf("# normalization of radial CHO eigenfunctions differs by %g from %s (method=%s)\n",
                                         dev[diag], diag?"unity":"zero ", method?"analytical":"numerical");
          } // diag
      } // method
      return (deviation[1][0] + deviation[1][1] > 5e-11);
  } // test_orthonormality

  
  template <int numax=7>
  inline status_t test_Gram_Schmidt(int const echo=1) {

      double dev{0};
      double c[1 + numax][1 + numax];
      for (int ell = 0; ell <= numax; ++ell) {
          if (echo > 3) std::printf("# %s ell=%d\n", __func__, ell);

          for (int i = 0; i < (1 + numax)*(1 + numax); ++i) {
              c[0][i] = 0; // clear polynomial coefficients
          } // i
          double fac[1 + numax]; // prefactors
          for (int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
              c[nrn][nrn] = (nrn % 2)?-1:1; // init as r^(ell + 2*nrn)
              for (int jrn = 0; jrn < nrn; ++jrn) {
                  // Gram-Schmidt-orthogonalize against previous
                  auto const a = inner_product(c[nrn], nrn, c[jrn], jrn, ell);
                  auto const d = inner_product(c[jrn], jrn, c[jrn], jrn, ell);
                  assert(d > 0);
                  auto const a_over_d = a/d;
                  for (int p = 0; p <= jrn; ++p) {
                      c[nrn][p] -= a_over_d*c[jrn][p]; // subtract state jrn from nrn
                  } // p
                  auto const check = inner_product(c[nrn], nrn, c[jrn], jrn, ell);
                  dev = std::max(dev, std::abs(check));
                  if (echo > 11) std::printf("# overlap between %c%d and %c%d is %g\n",
                                        ellchar[ell], nrn, ellchar[ell], jrn, check);
              } // jrn

              auto const f = radial_normalization(c[nrn], nrn, ell);
              fac[nrn] = f; // store prefactor for plotting

              if (echo > 0) {
                  std::printf("# %c%d-poly", ellchar[ell], nrn);
                  for (int p = nrn; p >= 0; --p) {
                      std::printf(" + %g*r^%d", f*c[nrn][p], ell+2*p);
                  } // p
                  std::printf("\n");
              } // echo
              
              if (1) { // check if radial_eigenstates delivers the same coefficients
                  double coeff[1 + numax];
                  radial_eigenstates(coeff, nrn, ell);
                  auto const fc = radial_normalization(coeff, nrn, ell);
                  std::printf("# %c%d-poly", ellchar[ell], nrn);
                  double diff{0};
                  for (int p = nrn; p >= 0; --p) {
                      std::printf(" + %g*r^%d", fc*coeff[p], ell+2*p);
                      diff = std::max(diff, std::abs(f*c[nrn][p] - fc*coeff[p]));
                  } // p
                  std::printf("  largest difference %.1e\n", diff);
              } // true

          } // nrn

          if (echo > 8) { // plot the functions
              bool const with_exp = true,
                         with_ell = true;
              std::printf("\n# Radial CHO functions for ell=%d\n", ell);
              for (int ir = 0; ir <= 1000; ++ir) {
                  auto const r = ir*0.01, r2 = r*r;
                  auto const Gauss = with_exp ? std::exp(-0.5*r2) : 1;
                  auto const r_ell = with_ell ? std::pow(r, ell)  : 1;
                  std::printf("%g", r);
                  for (int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
                      std::printf("\t%g", fac[nrn]*r_ell*expand_poly(c[nrn], 1 + nrn, r2)*Gauss);
                  } // nrn
                  std::printf("\n");
              } // ir
              std::printf("\n");
          } // echo

      } // ell

      if (echo > 1) std::printf("# %s largest deviation is %.1e\n", __func__, dev);
      return (dev > 1e-12);
  } // test_Gram_Schmidt
  
  
  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_orthonormality(0);
      stat += test_Gram_Schmidt(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS
  
} // namespace cho_radial
