/*
 *   Unitary transformation for circular harmonics
 */


#include <cstdio> // std::printf
#include <cstdlib> // std::abs
#include <cassert> // assert
// #include <complex> // std::complex<T>

typedef int status_t;

#include "cho_radial.hxx" // ::all_tests


//   class Monom {
//   public:
//       double  c;
//       int16_t x,y,z,s;
//   public:
//       Monom(double const coeff=0, char const xyz='?', unsigned const power=0)
//         : c(coeff)
//         , x(('x' == xyz)*power)
//         , y(('y' == xyz)*power)
//         , z(('z' == xyz)*power)
//         , s(('s' == xyz)*power)
//       {} // constructor
//       
//   }; // Monom
// 
//   Monom operator* (Monom const & left, Monom const & right) {
//       Monom m(left.c * right.c);
//       m.x = left.x + right.x;
//       m.y = left.y + right.y;
//       m.z = left.z + right.z;
//       m.s = left.s + right.s;
//       return m;
//   } // operator*

  template <int lmax=9>
  status_t generate_unitary_transform(int const echo) {
      int constexpr Lcut = 1 + lmax;
      // ToDo: modify this for CHO (2D isotropic harmonic oscillator)
    
      // generate sin(m*theta) and cos(m*theta):
//       Monom c[1 + lmax], s[1 + lmax];
//       c[0] = 1;
//       s[0] = 0;
//       if (lmax > 0) {
//           c[1] = Monom(1,'x'); s[1] = Monom(1,'y');
//           Monom const cph2(2,'x');
//           for (int m = 2; m <= lmax; ++m) {
//               s[m] = cph2*s[m - 1] - s[m - 2];
//               c[m] = cph2*c[m - 1] - c[m - 2];
//           } // m
//       } // lmax > 0

//       int16_t cs[Lcut][2][Lcut][Lcut]; // cos and sin
//       for (int i = 0; i < Lcut*2*Lcut*Lcut; ++i) cs[0][0][0][i] = 0; // clear
//       cs[0][0][0][0] = 1; // c[0] = 1
//       cs[1][0][0][1] = 1; // c[1] = x
//       cs[1][1][1][0] = 1; // s[1] = y
//       for (int m = 2; m < Lcut; ++m) {
//           for (int i = 0; i < m; ++i) {
//               for (int j = 0; j < m; ++j) {
//                   cs[m][1][i][j] = 2*cs[m - 1][1][i][j] - cs[m - 2][1][i][j];
//                   cs[m][0][i][j] = 2*cs[m - 1][0][i][j] - cs[m - 2][0][i][j];
//               } // j
//           } // i
//       } // m
      
      // (x + iy)^m --> decompose into real and imaginary part
//       std::complex<int16_t> cs[Lcut][Lcut];
//       for (int i = 0; i < Lcut*Lcut; ++i) cs[0][i] = 0; // clear
//       cs[0][0] = 1;
//       cs[1][0] = 1;
//       cs[1][1] = std::complex<int16_t>(0,1);
//       for (int m = 2; m < Lcut; ++m) {
//           for (int j = 0; j < m; ++j) {
//               cs[m][j]     += cs[1][0]*cs[m - 1][j];
//               cs[m][j + 1] += cs[1][1]*cs[m - 1][j];
//           } // j
//       } // m

      // (x + iy)^m --> decompose into real and imaginary part
      // 
      //               
      // (x + iy)^m = sum_k bnc(m over k) x^{m - k} (iy)^{k}
      //
      int32_t xpiy[Lcut][Lcut][2]; for (int i = 0; i < Lcut*Lcut*2; ++i) xpiy[0][0][i] = 0;
      { // bnc-scope
          std::vector<int> bnc(Lcut, 0); bnc[0] = 1; // init binomial coefficients
          for (int m = 0; m <= lmax; ++m) {

              if (echo > 2) std::printf("# binomcoeff l=%d  ", m);
              for (int k = 0; k <= m; ++k) {
                  int const j4 = k & 0x3; // modulo 4, mimique the behaviour of (complex i)^k, i^0=1:0b00, i^1=i:0b01, i^2=-1:0b10, i^3=-i:0b11
                  int const reim = j4 & 0x1; // least significant bit, 0:real, 1:imaginary
                  int const bit1 = (j4 >> 1) & 0x1; // second bit
                  if (echo > 2) std::printf(" %c%d x^%d %sy^%d   ", bit1?'-':' ', bnc[k], m - k, reim?"i":"", k);
                  int const sgn = 2*bit1 - 1;
                  xpiy[m][k][reim] = sgn*bnc[k];  
              } // k
              if (echo > 2) std::printf("\n");

              // prepare bnc for m+1
              for (int k = m + 1; k > 0; --k) {
                  bnc[k] += bnc[k - 1];
              } // k
          } // m
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
      double constexpr pi = 3.14159265358979323846; // pi 
      double constexpr sqrtpi = 1.77245385090551602729816748334115; // sqrt(pi)

      double metric[Lcut];
      metric[0] = sqrtpi;
      for (int n = 1; n < Lcut; ++n) {
          metric[n] = metric[n - 1]*(n - 0.5);
          std::printf("# metric n=%d  %g\n", n, metric[n]);
      } // n

      
      // prepare Hermite polynomials
      double H[Lcut][Lcut]; 
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
                  H[nu + 1][j] -= 0.5*nu * H[nu - 1][j];
              } // j
          } // i < 1
          
          std::printf("# Hermite polynomial #%i\t", i);
          for (int j = 0; j <= i; ++j) {
              std::printf("%10.4f", H[i][j]);
          } // j
          std::printf("\n");
      } // i


      double H_diagonal[Lcut];
      double xpiy_norm[Lcut];
      bool constexpr check_overlap = true;  // check the overlap w.r.t. the metric
      for (int i = 0; i < Lcut; ++i) {
          if (echo > 11) std::printf("# Hermite overlap #%i\t", i);
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
                  if (echo > 11) std::printf("%9.3f", ovl);
                  if (i == j) H_diagonal[i] = ovl; // store diagonal elements for normalization
              } // i == j
          } // j
          if (echo > 11) std::printf("\n");
          xpiy_norm[i] = 2*pi*cho_radial::exponential_integral_k(i);
      } // i

      // compute the inner products
      // < H(x)*H(y) | (x + iy)^m > under the metric exp(-r^2)
      // furthermore, we need < (x + iy)^m | (x + iy)^m > for normalization
      //    /infty
      // 2pi|dr   r^{2m + 1} exp(-r^2) --> see exponential_integral_k
      //    /0
      
      for (int nu = 0; nu < Lcut; ++nu) {
          for (int nx = 0; nx <= nu; ++nx) {  // Cartesian quantum number for x
              int const ny = nu - nx;         // Cartesian quantum number for y

              for (int m = -nu; m <= nu; ++m) {
                  int const reim = (m < 0);
                  int const mabs = std::abs(m);

                  double inner{0};
                  for (int k = 0; k <= mabs; ++k) {
                      double ovl[2] = {0, 0};
                      for (int d = 0; d < 2; ++d) { // direction 0:x, 1:y
                          int const n        = d ? ny : nx;
                          int const addpower = d ? k  : mabs - k;
                          int const parity   = n & 0x1; // 0:even, 1:odd
                          if (false) {
                              std::vector<double> product(2*Lcut, 0.0);
                              for (int p = 0; p <= n; ++p) {
                                  product[p + addpower] += H[n][p];
                              } // p
                              for (int p = 0; p < Lcut; ++p) {
                                  ovl[d] += product[2*p] * metric[p];
                              } // p
                          } else
                          if (parity == (addpower & 0x1)) { // only terms of the same parity are nonzero
//                               for (int p = 0; p <= n; ++p) {
//                                   int const pp = p + addpower;
//                                   if (0 == (pp & 0x1)) { // only even powers
//                                       ovl[d] += H[n][p] * metric[pp >> 1];
//                                   } // even
//                               } // p
                              for (int p = parity; p <= n; p += 2) {
                                  ovl[d] += H[n][p] * metric[(p + addpower) >> 1];
                              } // p
                          } // same parity
                      } // d
                      inner += ovl[0] * ovl[1] * xpiy[mabs][k][reim];
                  } // k

                  // normalize properly
                  inner /= std::sqrt(xpiy_norm[mabs] * H_diagonal[nx] * H_diagonal[ny]);
                  if (std::abs(inner) > 1e-14)
                  std::printf("# overlap of nu=%d nx=%d ny=%d m=%d \t%16.9f  \t%16.9f  \t%g\n", nu, nx, ny, m, inner, inner*inner, inner);
              } // m

          } // nx
      } // nu = nx + ny
      
      

      return 0;
  } // generate_unitary_transform

int main(int argc, char *argv[]) {
    cho_radial::all_tests(9);
    return generate_unitary_transform(9);
} // main
