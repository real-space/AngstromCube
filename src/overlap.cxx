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
  // computes the overlap between Gaussian-localized 1D polynomials

//   double constexpr C_PI = 3.14159265358979323846; // pi
  double constexpr sqrtpi = 1.77245385091;
  
  template<typename real_t>
  int multiply(real_t pxp[], int const n, // result
               real_t const p0[], int const n0,
               real_t const p1[], int const n1) {
    // multiplication of two polynomials
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
      //            / infty
      // kern_{n} = |   exp(-x^2/sigma^2) x^n dx
      //            /-infty
      real_t kern = sqrtpi * sigma; // init recursive computation
      for(int d = 0; 2*d < m; ++d) {
          value += p[2*d] * kern;
          kern *= (d + 0.5) * (sigma*sigma);
      } // d
      return value;
  } // integrate

  
  template<typename real_t>
  void prepare_centered_Hermite_polynomials(real_t H[], int const ncut,
                    double const siginv=1, double const normalize=1) {
      
      int const S = ncut; // access stride
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
  void derive_Hermite_Gauss_polynomials(real_t dH[], real_t const H[], int const ncut,
                    double const siginv=1) {
      // the Gaussian envelope function exp(-.5*(x/sigma)^2) is implicit here
      // but needs to be considered when deriving:
      // 
      // d/dx ( p(x)*exp(-x^2/2) ) = (d/dx p(x) - x*p(x)/sigma^2) * exp(-.5*(x/sigma)^2)
    
      // derive the polynomial first
      dH[ncut - 1] = 0;
      for(int d = 1; d < ncut; ++d) {
          dH[d - 1] = d*H[d];
      } // d
      
      // now add the terms coming from the inner derivative of exp(-.5*(x/sigma)^2)
      for(int d = 0; d < ncut - 1; ++d) {
          dH[d + 1] -= H[d]*siginv*siginv; // times (x/sigma^2)
      } // d

  } // derive
  
  template<typename real_t>
  void shift_polynomial_centers(real_t c_shifted[], // result: shifted polynomial
                                real_t const c[], // assume p(x) = sum_k=0...nmax-1 c[k] * x^k
                                int const nmax,
                                real_t const x_shift) {
    
      real_t c_old[nmax];
      for(int k = 0; k < nmax; ++k) {
          c_old[k] = c[k]; // get a work copy
      } // k

      double kfactorial = 1; // init kfactorial with 0! == 1
      for(int k = 0; k < nmax; ++k) { // loop MUST run forward from 0

          // evaluate the value of d^k p(x) / d x^k at x=x_shift
          real_t val = 0;
          {   real_t xsp = 1; // x_shift^p
              for(int p = 0; p < nmax - k; ++p) { // we only need to run up to nmax-k as the degree of the input poly is decreased with every k
                  val += xsp * c_old[p];
                  xsp *= x_shift; // update x_shift^p for the next p-iteration
              } // p
          } // ToDo: Horner-scheme could be used

          c_shifted[k] = val / kfactorial;

          // now derive the original polynomial, in-place, for the next k-iteration
          for(int p = 1; p < nmax - k; ++p) { // loop MUST run forward from 1
              c_old[p - 1] = p * c_old[p]; // d/dx x^p = p*x^{p-1}
          } // p
          c_old[nmax - k - 1] = 0;

          kfactorial *= (k + 1); // update kfactorial for the next k-iteration
      } // k

  } // shift_polynomial_centers
  
  
  template<typename real_t>
  real_t overlap_of_two_Hermite_Gauss_functions(
      real_t const H0[], int const n0, double const s0,
      real_t const H1[], int const n1, double const s1, 
      double const distance) {
      auto const k0 = 1/(s0*s0), k1 = 1/(s1*s1);
      auto const sigma = 1/sqrt(.5*(k0 + k1));

      auto const sh0 = -distance*k0/(k0 + k1);
      real_t H0s[n0]; // H0 shifted by sh0
      shift_polynomial_centers(H0s, H0, n0, sh0);

      auto const sh1 =  distance*k1/(k0 + k1);
      real_t H1s[n1]; // H1 shifted by sh1
      shift_polynomial_centers(H1s, H1, n1, sh1);

      int const m = n0 + n1;
      real_t h0xh1[m]; // product of H0s and H1s
      multiply(h0xh1, m, H0s, n0, H1s, n1);
      return integrate(h0xh1, m, sigma) * exp(-0.5*k0*sh0*sh0 -0.5*k1*sh1*sh1);
  } // overlap_of_two_Hermite_Gauss_functions

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
    int constexpr ncut = 8;
    double const sigma = 1.4567; // any positive real number
    double H[ncut*ncut];
    prepare_centered_Hermite_polynomials(H, ncut, 1./sigma);
    double hh[2*ncut];
    int ndev = 0; double mdev = 0;
    for(int n = 0; n < ncut; ++n) {
        if (echo > 3) plot_poly(&H[ncut*n], 1+n, "H");
        if (echo > 1) printf("# %s   %d   ortho", __func__, n);
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
    if (echo) printf("# %s: up to %d the largest deviation from Kroecker is %.1e \n", __func__, ncut - 1, mdev);
    return ndev;
  } // test

  status_t test_Hermite_Gauss_overlap(int const echo=1) {
    // show the overlap of the lowest 1D Hermite-Gauss functions as a function of distance
    int constexpr ncut = 4;
    double const sigma0 = 1.4567, sigma1 = sigma0 + .876; // any two positive real numbers
//     double const sigma0 = 1, sigma1 = sigma0;
    double H0[ncut*ncut], H1[ncut*ncut];
    prepare_centered_Hermite_polynomials(H0, ncut, 1/sigma0);
    prepare_centered_Hermite_polynomials(H1, ncut, 1/sigma1);
    for(auto dist = 0.0; dist < 11; dist += .1) {
        if (echo > 1) printf("# %s  distance=%.3f    ", __func__, dist);
        for(int n = 0; n < ncut; ++n) {
            for(int m = 0; m < ncut; ++m) {
                double const ovl = overlap_of_two_Hermite_Gauss_functions(
                  &H0[ncut*n], 1+n, sigma0,
                  &H1[ncut*m], 1+m, sigma1, dist);
                if (echo > 1) printf(" %.6f", ovl);
            } // m
        } // n
        if (echo > 1) printf("\n");
    } // dist
    return 0;
  } // test

  status_t test_kinetic_overlap(int const echo=3) {
    // show the kinetic energy of the lowest 1D Hermite-Gauss functions as a function of distance
    // test if the derivation operator can be cast to any side
    int constexpr ncut = 6;
    int constexpr mcut = ncut - 2;
    double const sigma0 = 1, sigma1 = sigma0 + .01; // fails for different sigmas
    double H0[ncut*ncut], H1[ncut*ncut];
    prepare_centered_Hermite_polynomials(H0, ncut, 1/sigma0);
    prepare_centered_Hermite_polynomials(H1, ncut, 1/sigma1);

    double dH0[ncut*mcut], dH1[ncut*mcut], d2H0[ncut*mcut], d2H1[ncut*mcut];
    for(int n = 0; n < mcut; ++n) {
        // first derivatives
        derive_Hermite_Gauss_polynomials(&dH0[ncut*n], &H0[ncut*n], ncut, 1/sigma0);
        derive_Hermite_Gauss_polynomials(&dH1[ncut*n], &H1[ncut*n], ncut, 1/sigma1);
        // second derivatives
        derive_Hermite_Gauss_polynomials(&d2H0[ncut*n], &dH0[ncut*n], ncut, 1/sigma0);
        derive_Hermite_Gauss_polynomials(&d2H1[ncut*n], &dH1[ncut*n], ncut, 1/sigma1);
    } // n
    double maxdev1 = 0, maxdev2 = 0;
    for(auto dist = 0.0; dist < 11; dist += .01) {
        if (echo > 1) printf("# %s  distance=%.3f    ", __func__, dist);
        for(int n = 0; n < mcut; ++n) {
            for(int m = 0; m < mcut; ++m) {
                auto const d2d0 = overlap_of_two_Hermite_Gauss_functions(&d2H0[ncut*n], ncut, sigma0, &H1[ncut*m], ncut, sigma1, dist);
                auto const d0d2 = overlap_of_two_Hermite_Gauss_functions(&H0[ncut*n], ncut, sigma0, &d2H1[ncut*m], ncut, sigma1, dist);
                auto const d1d1 = overlap_of_two_Hermite_Gauss_functions(&dH0[ncut*n], ncut, sigma0, &dH1[ncut*m], ncut, sigma1, dist);
                if (echo > 1) printf("  %.9f %.9f %.9f", d2d0, d0d2, -d1d1); // show 3 values
//              if (echo > 1) printf("  %.1e %.1e %.1e", d2d0 + d1d1, d0d2 + d1d1, d2d0 - d0d2); // show deviations
//              if (echo > 1) printf(" %.9f", -d1d1); // show 1 value
                maxdev2 = std::max(maxdev2, std::abs(d2d0 - d0d2));
                maxdev1 = std::max(maxdev1, std::abs(d2d0 + d1d1));
                maxdev1 = std::max(maxdev1, std::abs(d0d2 + d1d1));
            } // m
        } // n
        if (echo > 1) printf("\n");
    } // dist
    if (echo > 0) printf("# %s deviations %g and %g\n", __func__, maxdev1, maxdev2);
    return (maxdev1 > 2e-14) + (maxdev2 > 2e-14);
  } // test

  status_t all_tests() {
    auto status = 0;
    status += test_Hermite_polynomials();
    status += test_Hermite_Gauss_overlap();
    status += test_kinetic_overlap();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace overlap
