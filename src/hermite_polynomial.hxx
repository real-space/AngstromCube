#pragma once

#include <cmath> // std::exp

  template<typename real_t>
  inline void hermite_polys(real_t H[], real_t const x, int const numax, real_t const rcut=9) {
      // Hermite-polynomials times Gauss function, not normalized!
      real_t const H0 = (x*x < rcut*rcut)? std::exp(-0.5*x*x) : 0; // Gaussian envelope function

      H[0] = H0; // store Gaussian H_0

      real_t Hnup1, Hnu = H0, Hnum1 = 0;
      for(int nu = 0; nu < numax; ++nu) {
          real_t const nuhalf = 0.5 * (real_t)nu; // nu/2.
          // use the two step recurrence relation to create the 1D Hermite polynomials
          // H[nu+1] = x * H[nu] - nu/2. * H[nu-1];
          Hnup1 = x * Hnu - nuhalf * Hnum1; // see snippets/3Hermite.F90 for the derivation of this
          H[nu + 1] = Hnup1; // store
          // rotate registers for the next iteration
          Hnum1 = Hnu; Hnu = Hnup1; // ordering is important here
      } // nu

  } // hermite_polys

#ifndef NO_UNIT_TESTS
    #include <cstdio> // printf
    #include <algorithm> // std::max
    #include <cmath> // std::abs, std::sqrt
    #include <cassert> // assert
    
    #include "constants.hxx" // constants::pi
#endif // NO_UNIT_TESTS

namespace hermite_polynomial {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template<typename real_t>
  inline status_t test_hermite_polynomials(int const echo=4, double const threshold=7e-15) {
      int constexpr n = 1 << 13;
      int constexpr numax = 7, M = (1 + numax);
      double const hg = 1.25/(1 << 9);
      if (echo > 1) printf("# %s: numax = %d\n", __func__, numax);
      real_t hh[M*M];
      for(int ij = 0; ij < M*M; ++ij) hh[ij] = 0; // clear

      real_t hp[n*M];
      for(int i = 0; i < n; ++i) {
          real_t const x = (i - .5*(n - 1))*hg;
          hermite_polys(&(hp[i*M]), x, numax);
          for(int ij = 0; ij < M*M; ++ij) {
              hh[ij] += hp[i*M + ij%M]*hp[i*M + ij/M]; // integrate overlap matrix
          } // ij
          if (echo > 6) {
              printf("%g  ", x);
              for(int nu = 0; nu <= numax; ++nu) {
                  printf(" %g", hp[i*M + nu]);
              }   printf("\n");
          } // echo
      } // i
      if (echo > 5) printf("\n\n");
      
      for(int ij = 0; ij < M*M; ++ij) hh[ij] *= hg; // grid spacing
      
      double const pi_factor = 1./std::sqrt(constants::pi);
      double max_dev[] = {0, 0}; // {diagonal, off-diagonal}
      for(int nu = 0; nu <= numax; ++nu) {
          double const diag = hh[nu*M + nu];
          if (echo > 2) printf("# nu = %d norm^2 = %g --> %g\n", nu, diag, diag*pi_factor*(1 << nu));
          assert(diag > 0);
          auto const scal = 1./std::sqrt(diag);
          for(int nup = 0; nup <= numax; ++nup) {
              hh[nu*M + nup] *= scal;
              hh[nup*M + nu] *= scal;
              int const d = (nu == nup);
              max_dev[d] = std::max(max_dev[d], std::abs((double)hh[nu*M + nup] - d));
          } // nup
          assert(std::abs(hh[nu*M + nu] - 1) < 1e-7);
      } // nu

      if (echo > 3) {
          printf("# overlap matrix %d x %d\n", M, M);
          for(int nu = 0; nu <= numax; ++nu) {
              printf("# nu = %d ", nu);
              for(int nup = 0; nup <= numax; ++nup) {
                  printf("%16.9f", hh[nu*M + nup]);
              }   printf("\n");
          }   printf("\n\n");
      } // echo

      if (echo > 1) printf("# %s largest deviation from unit matrix is %.1e and %.1e on the diagonal\n", 
                            __func__, max_dev[0], max_dev[1]);
      return (max_dev[0] + max_dev[1] > threshold);
  } // test_hermite_polynomials

  inline status_t all_tests(int const echo=3) {
    if (echo > 0) printf("\n# %s: %s\n\n", __FILE__, __func__);
    auto status = 0;
    status += test_hermite_polynomials<float>(echo, 3.5e-6);
    status += test_hermite_polynomials<double>(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace hermite_polynomial
