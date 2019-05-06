#include <cstdio> // printf
#include <cstdlib> // abs
#include <cmath> // sqrt
#include <algorithm> // max

#include "overlap.hxx"

// #include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
// #include "output_units.h" // eV, _eV, Ang, _Ang

// #define FULL_DEBUG
#define DEBUG

#ifdef  DEBUG
    #include "debug_output.hxx" // dump_to_file
#endif

#ifdef FULL_DEBUG
    #define full_debug(print) print 
#else
    #define full_debug(print)
#endif

#ifdef DEBUG
    #define debug(print) print 
#else
    #define debug(print)
#endif


namespace overlap {
  // computes the overlap between Gaussian-localized 3D-factorizable polynomials

  double constexpr C_PI = 3.14159265358979323846; // pi
  double constexpr sqrtpi = 1.77245385091;
  
  
  template<typename real_t>
  int multiply(real_t pxp[], int const n, // result
               real_t const p0[], int const n0,
               real_t const p1[], int const n1) {
    for(int d = 0; d < n; ++d) {
        pxp[d] = 0; // clear
    } // d
    int nloss = 0;
    for(int d0 = 0; d0 < n0; ++d0) {
        for(int d1 = 0; d1 < n1; ++d1) {
            int const d = d0 + d1;
            if (d < n) {
                pxp[d] += p0[d0] * p1[d1];
            } else {
                nloss += (0 != p0[d0] * p1[d1]);
            }
        } // d1
    } // d0
    return nloss; // return a positive number if potentially non-zero coefficients have been lost because n < n0 + n1 - 1
  } // multiply

  
  template<typename real_t>
  real_t integrate(real_t const p[], int const m, double const sigma=1) {
      real_t value = 0;
      real_t kern = sqrtpi * sigma;
      for(int d = 0; 2*d < m; ++d) {
          value += p[2*d] * kern;
          kern *= (d + 0.5) * sigma*sigma;
      } // d
      return value;
  } // integrate
  
  
  template<typename real_t>
  void prepare_centered_Hermite_polynomials(real_t H[], int const ncut,
                    double const siginv=1, double const normalize=1) {
      
      int const S = ncut; // stride
      for(int i = 0; i < S*S; ++i) H[i] = 0; // clear
      
      H[S*0 + 0] = 1; // H_0 is the Gaussian envelope function exp(-.5*(x/sigma)^2) which is implicit here
      for(int n = 1; n < ncut; ++n) {
          for(int d = 0; d < n; ++d) {
              H[S*n + (d + 1)] = H[S*(n - 1) + d]*siginv; // times (x/sigma)
          } // d
          for(int d = 0; d < n - 1; ++d) {
              H[S*n + d] -= (0.5*(n - 1)) * H[S*(n - 2) + d];
          } // d
      } // n

      if (0 != normalize) {
          double nfactorial = 1;
          for(int n = 0; n < ncut; ++n) {
              double nrmf = normalize*sqrt(siginv/(sqrtpi*nfactorial));
              for(int d = 0; d <= n; ++d) {
                  H[S*n + d] *= nrmf;
              } // d
              nfactorial *= (n + 1.)*0.5; // update nfactorial
          } // n
      } // normalize

  } // prepare
  
  template<typename real_t>
  void plot_poly(real_t const poly[], int const m, char const *name) {
      printf("Poly %s : ", name);
      for(int d = 0; d < m; ++d) {
          printf("%.6f  ", poly[d]);
      } // d
      printf("\n");
  } // plot
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_Hermite_polynomials(int const echo=1) {
    // see if the first ncut Hermite polynomials are orthogonal and normalized
    int constexpr ncut = 16;
    double const sigma = 1.4567; // any positive real number
    double H[ncut*ncut];
    prepare_centered_Hermite_polynomials(H, ncut, 1./sigma);
    double hh[2*ncut];
    int ndev = 0; double mdev = 0;
    for(int n = 0; n < ncut; ++n) {
        if (echo > 3) plot_poly(&H[ncut*n], 1+n, "H");
        if (echo > 1) printf("%s   %d   ortho", __func__, n);
        for(int m = 0; m < ncut; ++m) {
            multiply(hh, 2*ncut, &H[ncut*n], 1+n, &H[ncut*m], 1+m);
            if (echo > 3) plot_poly(hh, 2*n - 1, "H^2");
            double const norm = integrate(hh, 2*ncut, sigma);
            mdev = std::max(mdev, fabs(norm - (m == n)));
            if (echo > 1) printf(" %.1e", norm - (m == n));
            ndev += (fabs(norm - (m == n)) > 1e-10); 
        } // m
        if (echo > 1) printf("\n");
    } // n
    if (echo) printf("%s: up to %d the largest deviation from Kroecker is %.1e \n", __func__, ncut - 1, mdev);
    return ndev;
  } // test
  
  status_t all_tests() {
    auto status = 0;
    status += test_Hermite_polynomials();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace overlap
