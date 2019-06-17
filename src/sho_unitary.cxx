#include <cstdio> // printf
#include <cassert> // assert
#include <cmath> // sqrt, pow, exp
#include <algorithm> // max

#include "inline_tools.hxx" // align
#include "sho_unitary.hxx"

namespace sho_unitary {
  
  void generate_unitary_transform(int const numax, int echo=9) {
//       int const mx = align<2>(1 + numax);

      int const lmax = numax;
      int const lx = (1 + lmax);
      int xory[lmax + 1 + lmax][lx]; // xory[lmax + m][k] --> for m <  0: x^(|m|- k)*y^k if k even
                                     //                   --> for m >= 0: y^(|m|- k)*x^k if k odd
      for(int k = 0; k < (2*lmax + 1)*lx; ++k) xory[0][k] = 0; // clear
      
      int Pl[lx][lx]; // generate (r^2 - u^2)^l
      for(int k = 0; k < lx*lx; ++k) Pl[0][k] = 0; // clear
      
      int rxy2pl[lx][lx]; // generate (x^2 + y^2)^l
      for(int k = 0; k < lx*lx; ++k) rxy2pl[0][k] = 0; // clear
      
      { // bnc-scope
          int bnc[lx]; bnc[0] = 1; for(int l = 1; l <= lmax; ++l) bnc[l] = 0; // init binomial coefficients
          for(int m = 0; m <= lmax; ++m) {
              
              if (echo > 2) printf("# Pl l=%d  ", m);
              int j4 = 1;
              for(int k = 0; k <= m; ++k) {
                  Pl[m][k] = bnc[k]*(1 - 2*(k & 1));
                  rxy2pl[m][k] = bnc[k];
                  if (echo > 2) printf(" %d", Pl[m][k]);
                  int const j4mod2 = j4 % 2;
                  int const msgn = 2*j4mod2 - 1;
                  int const sgn = 1 + j4mod2 - j4;
                  xory[lmax - m*msgn][k] = bnc[k]*sgn*msgn;
                  j4 = (j4 + 1) % 4; // mimique the behaviour of complex i^k
              } // k
              if (echo > 2) printf("\n");
              
              // prepare bnc
              for(int k = m + 1; k > 0; --k) {
                  bnc[k] += bnc[k - 1];
    //               if (echo > 8) printf("%6d", bnc[k]);
              } // k
    //           if (echo > 8) printf("%6d\n", bnc[0]);
          } // m
      } // bnc-scope

      if (echo > 2) printf("\n# (x^2 + y^2)^%d has highest coefficient %d\n", lmax, rxy2pl[lmax][lmax/2]);
      
      if (echo > 1) { 
          printf("\n\n");
          for(int m = -lmax; m <= lmax; ++m) {
              if (echo > 2) printf("# m=%3d   ", m);
              for(int k = 0; k <= std::abs(m); ++k) {
                  auto const xy = xory[lmax + m][k];
                  if (echo > 2) printf("%6d", xy);
              } // k
              if (echo > 2) printf("\n");
          } // m
          if (echo > 2) printf("\n\n");
      } // echo
      
      int64_t Plm[lx][lx][2*lx]; // data-type needs to capture factorial(2*lmax)
      for(int k = 0; k < lx*lx*2*lx; ++k) Plm[0][0][k] = 0; // clear
      
      for(int l = 0; l <= lmax; ++l) {
          int64_t poly[2*l + 1]; // polynomial in u
          for(int k = 0; k <= 2*l; ++k) poly[k] = 0; // clear
          for(int k = 0; k <= l;   ++k) poly[2*k] = Pl[l][k]; // copy underived

          for(int m = -l; m <= l; ++m) { // start at -l to derive l times for pure Legendre polynomials
              if (m >= 0) {          // before valid m range starts for associated Legendre Polynomials
                  if (echo > 2) printf("# Plm l=%d m=%d  ", l, m);
                  for(int k = 0; k <= l - m; ++k) {
                      Plm[l][m][k] = poly[k]; // store associated Legendre polynomial
                      if (echo > 2) printf(" %ldx^%d ", Plm[l][m][k], k);
                  } // k
                  if (echo > 2) printf("\n");
              }
              // derive (r^2 - u^2)^l w.r.t. u one time, i.e. for l+m times
              for(int k = 1; k <= 2*l; ++k) poly[k - 1] = k*poly[k]; poly[2*l] = 0; // derive in-place
          } // m
          for(int k = 0; k <= 2*l; ++k) assert(0 == poly[k]); // all coefficients of poly must vanish since we derive u^{2l} for 2l times
      } // l
      
      
      
      
  } // generate_unitary_transform
  

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  int test(int echo=9) {
      generate_unitary_transform(9, echo);
      return 0;
  } // test

  status_t all_tests() {
    auto status = 0;
    status += test();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  
  
} // namespace sho_unitary
