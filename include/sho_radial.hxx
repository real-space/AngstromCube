#pragma once

#include <cmath> // std::sqrt
#include <cassert> // assert

#include "constants.hxx" // ::sqrtpi

#ifndef NO_UNIT_TESTS
  #include <cstdio> // std::printf
  #include <algorithm> // std::max
#endif

#include "status.hxx" // status_t
#include "quantum_numbers.h" // ell_QN_t

namespace sho_radial {

  template <typename real_t>
  void radial_eigenstates(
        real_t poly[] // coefficients of a polynomial in r^2
      , int const nrn // number of radial nodes
      , int const ell // angular momentum quantum number
      , real_t const factor=1 // if we know the normalization prefactor in advance, we can provide it here
  ) {

      // recursion relation of coefficients:
      // a_0 = 1, a_{k+1} = (k-nrn)/((k+1)*(ell+k+3/2)) a_k

      poly[0] = factor;
      for (int k = 0; k < nrn; ++k) {
          // from https://quantummechanics.ucsd.edu/ph130a/130_notes/node244.html
          poly[k + 1] = (poly[k]*(k - nrn)*2)/(real_t)((k + 1)*(2*ell + 2*k + 3));
      } // k

  } // radial_eigenstates

  template <typename real_t>
  real_t exponential_integral_k(int const k) {
      // I_k = int\limit_0^\infty dr r^k exp(-r^2)
      if (0 == k) return real_t(0.5*constants::sqrtpi);
      if (1 == k) return real_t(0.5); // just for completeness, odd cases are not relevant in this module
      assert(k > 0);
      return real_t(0.5)*(k - 1)*exponential_integral_k<real_t>(k - 2); // recursive invokation
  } // exponential_integral_k

  template <typename real_t>
  real_t radial_normalization(real_t const coeff[], int const nrn, int const ell) {

      auto const prod = new real_t[2*nrn + 1]; // coefficients of a polynomial in r^2
      for (int p = 0; p < 2*nrn + 1; ++p) {
          prod[p] = 0;
      } // p
      for (int k = 0; k <= nrn; ++k) {
          for (int p = 0; p <= nrn; ++p) {
              prod[k + p] += coeff[k]*coeff[p]; // polynomial product with itself
          } // p
      } // k

      real_t exp_int_k = exponential_integral_k<real_t>(2*ell + 2);
      real_t norm{0};
      for (int p = 0; p <= 2*nrn; ++p) { // loop must run serial forward
  //    assert(exp_int_k == exponential_integral_k<real_t>(2*p + 2*ell + 2)); // DEBUG
        norm += prod[p]*exp_int_k;
        exp_int_k *= (p + ell + 1.5); // prepare exp_int_k for the next iteration
      } // p
      delete[] prod;
      return 1./std::sqrt(norm); // normalization prefactor
  } // radial_normalization
  
  template <typename real_t>
  real_t radial_normalization(int const nrn, int const ell) {
      auto const coeff = new real_t[nrn + 1]; // coefficients of a polynomial in r^2
      radial_eigenstates(coeff, nrn, ell); // get polynomial coefficients of the SHO eigenstates
      auto const result = radial_normalization(coeff, nrn, ell);
      delete[] coeff;
      return result;
  } // radial_normalization

  template <typename real_t>
  real_t expand_poly(real_t const coeff[], int const ncoeff, double const x) {
      real_t value{0};
      double xpow{1};
      for (int i = 0; i < ncoeff; ++i) {
          value += coeff[i] * xpow;
          xpow *= x;
      } // i
      return value;
  } // expand_poly
  
#ifdef NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <typename real_t>
  real_t numerical_norm(real_t const c0[], int const nrn0, 
                        real_t const c1[], int const nrn1, int const ell) {
      double constexpr dr = 1/((real_t)(1 << 12)), rmax = 12.;
      int const nr = rmax / dr;
      real_t norm{0};
      for (int ir = 0; ir < nr; ++ir) {
          double const r = (ir - .5)*dr;
          double const r2 = r*r;
          double const Gauss = std::exp(-r2); // == exp(-0.5*r^2)*exp(-0.5*r^2)
          double const f0 = expand_poly(c0, 1 + nrn0, r2);
          double const f1 = expand_poly(c1, 1 + nrn1, r2);
          double const r2pow = std::pow(r2, 1 + ell); // == r^2 * (r^ell)^2
          norm += Gauss * r2pow * f0 * f1;
      } // ir
      norm *= dr; // can be taken out of the loop since it is constant
      return norm;
  } // numerical_norm

  inline status_t test_orthonormality(int const echo=1) {
      int const numax = 9;
      int const n = (numax*(numax + 4) + 4)/4; // sho_tools::nSHO_radial(numax);
      if (echo > 1) std::printf("# %s  numax= %d has %d different radial SHO states\n", __func__, numax, n);
      
      double c[n][8]; // polynomial coefficients

      ell_QN_t ell_list[n];
      enn_QN_t nrn_list[n];
      double   fac_list[n];
      
      int i = 0;
//       for (int ene = 0; ene <= numax; ++ene) { //  E_SHO = ene + 3/2
//           for (int nrn = ene/2; nrn >= 0; --nrn) { // start from low ell withing the block, i.e. many radial nodes
//               int const ell = ene - 2*nrn;
//       for (int ell = 0; ell <= numax; ++ell) { //
      for (int ell = numax; ell >= 0; --ell) { //
          for (int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
            
              ell_list[i] = ell;
              nrn_list[i] = nrn;
              radial_eigenstates(c[i], nrn, ell);
              double const fac = radial_normalization(c[i], nrn, ell);
              assert(std::abs(fac - radial_normalization<double>(nrn, ell)) < 1e-12); // check that the other interface produces the same value,with -O0 we can even check for ==
              fac_list[i] = fac;
              if (echo > 2) std::printf("# %s %3d state  nrn= %d  ell= %d  factor= %g\n", __func__, i, nrn, ell, fac);
              radial_eigenstates(c[i], nrn, ell, fac_list[i]); // overwrite the coefficent series with the normalized onces

              ++i; // count the number of states
          } // nrn
      } // ene
      assert(n == i); // check that nSHO_radial(numax) agrees

      double dev[] = {0, 0};
      for (int i = 0; i < n; ++i) {
          int const ell = ell_list[i];
          for (int j = i; j >= 0; --j) { // triangular loop structure
//        for (int j = 0; j < n; ++j) {  // block loop structure
              if (ell == ell_list[j]) {
                  auto const delta_ij = numerical_norm(c[i], nrn_list[i], c[j], nrn_list[j], ell);
                  int const i_equals_j = (i == j);
                  if (false) {
                      if (i != j) {
                          auto const delta_ji = numerical_norm(c[j], nrn_list[j], c[i], nrn_list[i], ell);
                          assert(std::abs(delta_ji - delta_ij) < 1e-15);
                      } // i != j
                  } // sanity check
                  if (echo > 4 - i_equals_j) std::printf("%g ", delta_ij - i_equals_j);
                  // measure the deviation from unity for diagonal elements
                  dev[i_equals_j] = std::max(dev[i_equals_j], std::abs(delta_ij - i_equals_j));
              } // ell matches
          } // j
          if (echo > 3) std::printf("\n");
      } // i
      if (echo > 0) std::printf("# normalization of radial SHO eigenfunctions differs by %g from unity\n", dev[1]); // summary
      if (echo > 0) std::printf("# orthogonality of radial SHO eigenfunctions differs by %g from zero\n",  dev[0]); // summary
      return (dev[0] + dev[1] > 5e-11);
  } // test_orthonormality

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_orthonormality(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS
  
} // namespace sho_radial
