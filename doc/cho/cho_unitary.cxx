/*
 *   Unitary transformation for circular harmonics
 */


#include <cstdio> // std::printf
#include <cstdlib> // std::abs
#include <cassert> // assert

typedef int status_t;

#include "cho_radial.hxx" // ::all_tests

  template <typename T> int sgn(T val) {
      return (T(0) < val) - (val < T(0));
  } // sgn from https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c/10133700

  template <int lmax=9>
  status_t generate_unitary_transform(int const echo) {
      int constexpr Lcut = 1 + lmax;
      // ToDo: modify this for CHO (2D isotropic harmonic oscillator)
    
      // (x + iy)^ell --> decompose into real and imaginary part
      // 
      //               
      // (x + iy)^ell = sum_k bnc(ell over k) x^{ell - k} (iy)^{k}
      //
      int32_t xpiy[Lcut][Lcut][2]; for (int i = 0; i < Lcut*Lcut*2; ++i) xpiy[0][0][i] = 0;
      { // bnc-scope
          std::vector<int> bnc(Lcut + 1, 0); bnc[0] = 1; // init binomial coefficients for ell=0
          if (echo > 2) std::printf("\n# (x + iy)^ell\n");
          for (int ell = 0; ell <= lmax; ++ell) {

              if (echo > 2) std::printf("# (x + iy)^%d = ", ell);
              for (int k = 0; k <= ell; ++k) {
                  int const j4 = k & 0x3; // modulo 4, mimique the behaviour of (complex i)^k, i^0=1:0b00, i^1=i:0b01, i^2=-1:0b10, i^3=-i:0b11
                  int const reim = j4 & 0x1; // least significant bit, 0:real, 1:imaginary
                  int const bit1 = (j4 >> 1) & 0x1; // second bit
                  if (echo > 2) std::printf(" %c%d x^%d %sy^%d   ", bit1?'-':' ', bnc[k], ell - k, reim?"i":"", k);
                  int const sgn = 1 - 2*bit1;
                  xpiy[ell][k][reim] = sgn*bnc[k];
              } // k
              if (echo > 2) { std::printf("\n"); std::fflush(stdout); }

              // prepare bnc for ell+1
              for (int k = ell + 1; k > 0; --k) {
                  bnc[k] += bnc[k - 1];
              } // k
          } // ell
      } // bnc-scope


      int constexpr nr2cut = (Lcut + 1)/2;
      // construct (x + iy)^ell * (x^2 + y^2)^nrn
      int32_t xpiy_r2[nr2cut][Lcut][Lcut][2]; for (int i = 0; i < nr2cut*Lcut*Lcut*2; ++i) xpiy_r2[0][0][0][i] = 0;
      { // bnc-scope
          std::vector<int> bnc(1 + nr2cut, 0); bnc[0] = 1; // init binomial coefficients for nrn=0
          if (echo > 2) std::printf("\n# (x + iy)^ell (x^2 + y^2)^nrn\n");
          for (int nrn = 0; nrn < nr2cut; ++nrn) {
              // (x^2 + y^2)^nrn = sum_{irn=0}^nrn bnc(nrn over irn) x^{2(nrn - irn)} y^{2*irn}

              for (int ell = 0; ell <= lmax - 2*nrn; ++ell) {
                  for (int irn = 0; irn <= nrn; ++irn) { // y^{2*irn}
                      for (int k = 0; k <= ell; ++k) { // y^k
                          assert(k + 2*irn < Lcut);
                          for (int reim = 0; reim < 2; ++reim) {
                              xpiy_r2[nrn][ell][k + 2*irn][reim] += bnc[irn] * xpiy[ell][k][reim];
                          } // reim
                      } // k
                  } // irn

                  if (echo > 2) {
                      for (int reim = 0; reim < 2; ++reim) {
                          std::printf("# %s (x + iy)^%d * r^%d = ", reim?"Im":"Re", ell, 2*nrn);
                          int const kmax = ell + 2*nrn;
                          for (int k = 0; k <= kmax; ++k) {
                              auto const c = xpiy_r2[nrn][ell][k][reim];
                              if (0 != c) std::printf("  %d x^%d y^%d", c, kmax - k, k);
                          } // k
                          std::printf("\n");
                      } // reim
                  } // echo
              } // ell

              // prepare bnc for nrn+1
              for (int k = nrn + 1; k > 0; --k) {
                  bnc[k] += bnc[k - 1];
              } // k
          } // nrn
      } // bnc-scope
      

      // prepare the metric:
      //
      //  /infty
      //  | dx   exp(-x^2) x^(i + j)
      //  /-infty
      // 
      // (i + j) odd --> zero due to symmetry
      // (i + j) even, i + j = 2*n, then
      //  
      //  \sqrt{\pi} (2n - 1)!! / 2^n
// #define SIMPLIFY_PI
#ifdef  SIMPLIFY_PI
      double constexpr pi = 1, sqrtpi = 1; // pi, sqrt(pi) will be eleminated from the formulas, so we can set it to unity
#else  // SIMPLIFY_PI
      double constexpr pi = 3.14159265358979323846; // pi 
      double constexpr sqrtpi = 1.77245385090551602729816748334115; // sqrt(pi)
#endif // SIMPLIFY_PI

      double metric[Lcut]; // integral x^(2*ell) exp(-x^2) dx from -infty to infty
      metric[0] = sqrtpi;
      if (echo > 9) std::printf("\n# metric[0]\t%16.9f\n", metric[0]);
      for (int n = 1; n < Lcut; ++n) {
          metric[n] = metric[n - 1]*(n - 0.5);
          if (echo > 9) std::printf("# metric[%d]\t%16.9f\n", 2*n, metric[n]);
      } // n

      
      double H[Lcut][Lcut]; // Hermite polynomials
      for (int i = 0; i < Lcut; ++i) {
          for (int j = 0; j < Lcut; ++j) H[i][j] = 0; // clear
          if (i < 2) {
              H[i][i] = 1;
          } else { // i < 1
              // H[nu+1] = x * H[nu] - nu/2. * H[nu-1];
              int const nu = i - 1;
              for (int j = 0; j < i; ++j) {
                  H[nu + 1][j + 1] = H[nu][j]; // *x
              } // j
              for (int j = 0; j < nu; ++j) {
                  H[nu + 1][j] -= 0.5*nu * H[nu - 1][j]; // *sigma^2
              } // j
          } // i < 1
          
          std::printf("# Hermite polynomial #%i\t", i);
          for (int j = 0; j <= i; ++j) {
              std::printf("%10.4f", H[i][j]);
          } // j
          std::printf("\n");
      } // i


      double H_norm2[Lcut]; // normalization constant for Hermite polynomials
      bool constexpr check_overlap = true;  // check the overlap w.r.t. the metric
      for (int i = 0; i < Lcut; ++i) {
          if (echo > 7) std::printf("# Hermite overlap #%i\t", i);
          for (int j = 0; j < Lcut; ++j) {
              if (check_overlap || i == j) {
                  std::vector<double> product(2*Lcut, 0.0);
                  for (int ip = 0; ip <= i; ++ip) {
                      for (int jp = 0; jp <= j; ++jp) {
                          product[ip + jp] += H[i][ip] * H[j][jp];
                      } // jp
                  } // ip

                  double ovl{0};
                  for (int p = 0; p < Lcut; ++p) {
                      ovl += product[2*p]*metric[p];
                  } // p
                  if (echo > 7) std::printf("%9.3f", ovl);
                  if (i == j) H_norm2[i] = ovl; // store diagonal elements for normalization
              } // i == j
          } // j
          if (echo > 7) std::printf("\n");
      } // i

      // We need < (x + iy)^ell | (x + iy)^ell > with the metric exp(-r^2) for normalization
      //    /infty
      // 2pi|dr   r^{2ell + 1} exp(-r^2) --> see exponential_integral_k
      //    /0
      double xpiy_norm2[Lcut]; // normalization constant for (x+iy)^ell
//    double f{pi};
      for (int ell = 0; ell < Lcut; ++ell) {
          xpiy_norm2[ell] = 2*pi*cho_radial::exponential_integral_k(ell);
          if (echo > 7) std::printf("# <(x+iy)^%d|(x+iy)^%d> = %g\n", ell, ell, xpiy_norm2[ell]);
//        f *= (ell + 1);
      } // ell

      // What if we consider that we have only the real or only the imaginary part of it for ell > 0?
      for (int ell = 1; ell < Lcut; ++ell) {
          xpiy_norm2[ell] *= 0.5;
      } // ell


      // compute the inner products
      // < H_nx(x)*H_ny(y) | (x + iy)^ell (x^2 + y^2)^nrn > under the metric exp(-r^2)
      // and under the constraint that nx + ny == nu == ell + 2*nrn

      bool constexpr check_unitarity = true;
      double maxdev{0};
      
      
      int nmatrix{0};
//    double r2mat[Lcut][Lcut][Lcut]; for (int i = 0; i < Lcut*Lcut*Lcut; ++i) r2mat[0][0][i] = 0;
      double mat[Lcut][Lcut];
      for (int nu = 0; nu <= lmax; ++nu) {
          for (int i = 0; i < Lcut*Lcut; ++i) mat[0][i] = 0; // clear
          for (int nx = 0; nx <= nu; ++nx) {  // Cartesian quantum number for x
              int const ny = nu - nx;         // Cartesian quantum number for y

              int j{0};
              for (int ell = (nu & 0x1); ell <= nu; ell += 2) { // start from nu%2 to get the same parity

                  int const nrn = (nu - ell)/2; // number of radial nodes
                  assert(nrn < nr2cut);

                  double const one = 1;
                  auto const radial_norm0 = cho_radial::radial_normalization(&one, 0, ell);
                  std::vector<double> radial_coeff(1 + nrn, 0.0);
                  cho_radial::radial_eigenstates(radial_coeff.data(), nrn, ell);
                  auto const radial_norm = cho_radial::radial_normalization(radial_coeff.data(), nrn, ell) / radial_norm0;

                  if (echo > 22) {
                      std::printf("\n# radial eigenstate for nrn=%d is r^%d*(", nrn, ell);
                      for (int irn = 0; irn <= nrn; ++irn) {
                          std::printf("+ %g r^%d", radial_norm*radial_coeff[irn], 2*irn);
                      } // irn
                      std::printf(")\n");
                  } // echo

                  for (int reim = 0; reim <= (ell > 0); ++reim) { // 0:(m >= 0), 1:(m < 0)
                      int const m = ell*(1 - 2*reim);

                      double matrix{0};
                      for (int irn = 0; irn <= nrn; ++irn) { // loop over coefficients of (x^2 + y^2)^irn

                          double inner{0};
                          int const kmax = ell + 2*irn;
                          for (int k = 0; k <= kmax; ++k) {
                              double ovl[2] = {0, 0}; // inner products of Hermite polynomials with a pure power of x or y
                              for (int d = 0; d < 2; ++d) { // direction 0:x, 1:y
                                  int const n         = d ? ny : nx; // Cartesian quantum number
                                  int const addpower  = d ? k  : kmax - k;
                                  int const parity    = n & 0x1; // 0:even, 1:odd
                                  if (parity == (addpower & 0x1)) { // only terms of the same parity are nonzero
                                      for (int p = parity; p <= n; p += 2) {
                                          ovl[d] += H[n][p] * metric[(p + addpower)/2];
                                      } // p
                                  } // same parity
                              } // d
                              auto const add2inner = ovl[0] * ovl[1] * xpiy_r2[irn][ell][k][reim];
                              if (0 != add2inner && echo > 22)
                                  std::printf("# for nu=%d nx=%d ny=%d nr=%d m=%d ir=%d  ovl_x * ovl_y * xpiy_r2 = %g * %g * %d = %g\n", nu, nx, ny, nrn, m, irn, ovl[0], ovl[1], xpiy_r2[irn][ell][k][reim], add2inner);
                              inner += add2inner;
                          } // k
                          matrix += radial_coeff[irn] * inner;

                      } // irn

                      // normalize properly
                      matrix *= radial_norm/std::sqrt(H_norm2[nx] * H_norm2[ny] * xpiy_norm2[ell]);
                      if (std::abs(matrix) > 1e-14) {
                          if (echo > 22) std::printf("# for nu=%d nx=%d ny=%d nr=%d m=%d  normalize by sqrt of %g * %g * %g = %g\n",
                              nu, nx, ny, nrn, m, H_norm2[nx], H_norm2[ny], xpiy_norm2[ell], H_norm2[nx] * H_norm2[ny] * xpiy_norm2[ell]);
                          // show matrix elements as: with 9 digits, squared with 9 digits, in %g format, squared with sign in %g format
                          std::printf("# overlap of nu=%d nx=%d ny=%d nr=%d m=%d\t%16.9f\t%16.9f\t%g\t%g\n",
                                      nu, nx, ny, nrn, m, matrix, matrix*matrix, matrix, sgn(matrix)*matrix*matrix);
                          ++nmatrix;
                          mat[nx][j] = matrix;
                      } // nonzero

                      ++j;
                  } // reim
              } // ell
              assert(nu + 1 == j);

          } // nx

          if (check_unitarity) {
              double dev[][2] = {{0, 0}, {0, 0}};
              for (int i = 0; i <= nu; ++i) {
                  for (int j = 0; j <= nu; ++j) {
                      int const diag = (i == j);
                      double mmT{0}, mTm{0};
                      for (int k = 0; k <= nu; ++k) {
                          mmT += mat[i][k] * mat[j][k];
                          mTm += mat[k][i] * mat[k][j];
                      } // k
                      if (echo > 33) std::printf("# unitarity of nu=%d nx=%d j=%d %g %g\n", nu, i, j, mmT - diag, mTm - diag);
                      dev[0][diag] = std::max(dev[0][diag], std::abs(mmT - diag));
                      dev[1][diag] = std::max(dev[1][diag], std::abs(mTm - diag));
                  } // j
              } // i
              std::printf("# unitarity of nu=%d max deviations are %.1e, %.1e and %.1e, %.1e\n", nu, dev[0][0], dev[0][1], dev[1][0], dev[1][1]);
              maxdev = std::max(dev[0][0], std::max(dev[0][1], std::max(dev[1][0], std::max(dev[1][1], maxdev))));
          } // check_unitarity

      } // nu = nx + ny
      std::printf("# overlap up to lmax=%d has %d nonzero matrix elements\n", lmax, nmatrix);
      if (check_unitarity) std::printf("# unitarity of nu up to %d has max deviations of %.2e\n", lmax, maxdev);

      return (maxdev > 1e-13);
  } // generate_unitary_transform

int main(int argc, char *argv[]) {
    int const echo = (argc > 1) ? std::atoi(argv[1]) : 3;
    cho_radial::all_tests(echo);
    return generate_unitary_transform(echo);
} // main
