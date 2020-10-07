#include <cstdio> // printf
#include <cassert> // assert
#include <cmath> // std::pow, std::exp, std::abs
#include <algorithm> // std::max

#include "sho_radial.hxx"

#include "constants.hxx" // ::sqrtpi
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t
#include "inline_tools.hxx" // align<nbits>

namespace sho_radial {

#ifndef NO_UNIT_TESTS

  template<typename real_t>
  real_t numerical_norm(real_t const c0[], int const nrn0, 
                        real_t const c1[], int const nrn1, int const ell) {
      double constexpr dr = 1/((real_t)(1 << 12)), rmax = 12.;
      int const nr = rmax / dr;
      real_t norm = 0;
      for(int ir = 0; ir < nr; ++ir) {
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

  int test_orthonormality(int const echo=1) {
      int const numax = 9;
      int const n = (numax*(numax + 4) + 4)/4; // sho_tools::nSHO_radial(numax);
      if (echo > 1) printf("# %s  numax= %d has %d different radial SHO states\n", __func__, numax, n);
      
      double c[n][8]; // polynomial coefficients

      ell_QN_t ell_list[n];
      enn_QN_t nrn_list[n];
      double   fac_list[n];
      
      int i = 0;
//       for(int ene = 0; ene <= numax; ++ene) { //  E_SHO = ene + 3/2
//           for(int nrn = ene/2; nrn >= 0; --nrn) { // start from low ell withing the block, i.e. many radial nodes
//               int const ell = ene - 2*nrn;
//       for(int ell = 0; ell <= numax; ++ell) { //
      for(int ell = numax; ell >= 0; --ell) { //
          for(int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
            
              ell_list[i] = ell;
              nrn_list[i] = nrn;
              radial_eigenstates(c[i], nrn, ell);
              double const fac = radial_normalization(c[i], nrn, ell);
              assert(std::abs(fac - radial_normalization<double>(nrn, ell)) < 1e-12); // check that the other interface produces the same value,with -O0 we can even check for ==
              fac_list[i] = fac;
              if (echo > 2) printf("# %s %3d state  nrn= %d  ell= %d  factor= %g\n", __func__, i, nrn, ell, fac);
              radial_eigenstates(c[i], nrn, ell, fac_list[i]); // overwrite the coefficent series with the normalized onces

              ++i; // count the number of states
          } // nrn
      } // ene
      assert(n == i); // check that nSHO_radial(numax) agrees

      double dev[] = {0, 0};
      for(int i = 0; i < n; ++i) {
          int const ell = ell_list[i];
          for(int j = i; j >= 0; --j) { // triangular loop structure
//        for(int j = 0; j < n; ++j) {  // block loop structure
              if (ell == ell_list[j]) {
                  auto const delta_ij = numerical_norm(c[i], nrn_list[i], c[j], nrn_list[j], ell);
                  int const i_equals_j = (i == j);
                  if (false) {
                      if (i != j) {
                          auto const delta_ji = numerical_norm(c[j], nrn_list[j], c[i], nrn_list[i], ell);
                          assert(std::abs(delta_ji - delta_ij) < 1e-15);
                      } // i != j
                  } // sanity check
                  if (echo > 4 - i_equals_j) printf("%g ", delta_ij - i_equals_j);
                  // measure the deviation from unity for diagonal elements
                  dev[i_equals_j] = std::max(dev[i_equals_j], std::abs(delta_ij - i_equals_j));
              } // ell matches
          } // j
          if (echo > 3) printf("\n");
      } // i
      if (echo > 0) printf("# normalization of radial SHO eigenfunctions differs by %g from unity\n", dev[1]); // summary
      if (echo > 0) printf("# orthogonality of radial SHO eigenfunctions differs by %g from zero\n",  dev[0]); // summary
      return (dev[0] + dev[1] > 5e-11);
  } // test_orthonormality

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_orthonormality(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  
  
} // namespace sho_radial
