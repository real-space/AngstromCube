// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdlib> // std::abs
#include <cmath> // std::sqrt, ::exp, ::pow, ::abs
#include <algorithm> // std::max
#include <complex> // std::complex<real_t>
#include <vector> // std::vector<T>
#include <array> // std::array<T,n>
#include <cassert> // assert

#include "sho_overlap.hxx"

#include "vector_math.hxx" // ::vec<N,T>
#include "constants.hxx" // ::pi, ::sqrtpi
#include "control.hxx" // ::get
#include "inline_math.hxx" // pow2, set
#include "data_view.hxx" // view2D<T>
#include "sho_tools.hxx" // ::nSHO
#include "linear_algebra.hxx" // ::inverse
#include "recorded_warnings.hxx" // warn
#include "display_units.h" // Ang, _Ang, eV, _eV
#include "print_tools.hxx" // printf_vector

#ifndef   NO_UNIT_TESTS
    #include "simple_math.hxx" // ::random<real_or_int_t>
#endif // NO_UNIT_TESTS undefined

// #define FULL_DEBUG
// #define DEBUG

#ifdef    FULL_DEBUG
    #define full_debug(print) print
#else  // FULL_DEBUG
    #define full_debug(print)
#endif // FULL_DEBUG

#ifdef    DEBUG
    #include "debug_output.hxx" // dump_to_file
    #define debug(print) print
#else  // DEBUG
    #define debug(print)
#endif // DEBUG


namespace sho_overlap {
  // computes the overlap between Gaussian-localized 1D polynomials

  template <typename real_t>
  int multiply( // returns the number of non-zero coefficients that did not fit into p0xp1
        real_t    p0xp1[], int const n  // result
      , real_t const p0[], int const n0 // left input
      , real_t const p1[], int const n1 // right input
  ) {
      // multiplication of two polynomials

      set(p0xp1, n, real_t(0)); // clear result polynomial coefficients
      int nloss{0};
      for (int d0{0}; d0 < n0; ++d0) {
          for (int d1{0}; d1 < n1; ++d1) {
              int const d = d0 + d1;
              if (d < n) {
                  p0xp1[d] += p0[d0] * p1[d1];
              } else {
                  nloss += (0 != (p0[d0] * p1[d1]));
              } // d < n
          } // d1
      } // d0
      return nloss; // return a positive number if potentially non-zero coefficients have been lost because n < n0 + n1 - 1
  } // multiply


  template <typename real_t>
  double integrate( // returns the analytically integrated value
        real_t const p[] // input polynomial coefficients
      , int const m // number of coefficients in p[]
      , double const sigma=1.0 // Gaussian spread
      , int const moment=0 // shifting of powers
  ) {
      // contract the polynomial p[k] x^{k + moment} with kern_{k + moment}
      //
      //            / infty
      // kern_{n} = |   exp(-x^2/sigma^2) x^n dx  only non-zero for even n
      //            /-infty
      //
      assert( 0 <= moment );
      double value{0};
      double kern = constants::sqrtpi * sigma; // init recursive computation
      for (int d{0}; (2*d - moment) < m; ++d) {
          int const ip = 2*d - moment;
          if (ip >= 0) value += p[ip] * kern;
          kern *= (d + 0.5) * (sigma*sigma);
      } // d
      return value;
  } // integrate


  template <typename real_t>
  void prepare_centered_Hermite_polynomials(
        real_t H[] // result: Hermite polynomials H[i*ncut + j] up to order n < ncut
      , int const ncut // number of coefficients in H[]
      , double const sigmainv=1.0 // inverse of Gaussian spread sigma
      , double const normalize=1.0 // normalize to this number, except if 0.0 exactly
  ) {
      // prepare the polynomial coefficient matrix of all ncut-1 centered Hermite polynomials

      int const S = ncut; // access stride
      for (int i{0}; i < S*S; ++i) { H[i] = 0; } // clear

      H[0*S + 0] = 1; // H_0 is the Gaussian envelope function exp(-.5*(x/sigma)^2) which is implicit here
      for (int n{1}; n < ncut; ++n) { // recursion in n
          for (int d{0}; d < n; ++d) {
              H[n*S + (d + 1)] = H[(n - 1)*S + d] * sigmainv; // times (x/sigma)
          } // d
          for (int d{0}; d < n - 1; ++d) {
              H[n*S + d] -= (0.5*(n - 1)) * H[(n - 2)*S + d];
          } // d
      } // n

      if (0.0 != normalize) {
          double nfactorial{1};
          for (int n{0}; n < ncut; ++n) {
              double const nrmf = normalize*std::sqrt(sigmainv/(constants::sqrtpi*nfactorial));
              for (int d{0}; d <= n; ++d) {
                  H[n*S + d] *= nrmf; // scale
              } // d
              nfactorial *= (n + 1)*0.5; // update nfactorial to be n!/2^n
          } // n
      } // normalize

  } // prepare_centered_Hermite_polynomials


  template <typename real_t>
  int derive_Hermite_Gauss_polynomial( // returns 1 if the highest term could not be stored
        real_t dH[] // output derived Hermite polynomial coefficients
      , real_t const H[] // input Hermite polynomial coefficients H[n]
      , int const n // number of coefficients for both, H[] and dH[]
      , double const sigmainv=1.0 // inverse of the Gaussian spread sigma
  ) {
      // the Gaussian envelope function exp(-.5*(x/sigma)^2) is implicit here
      // but needs to be considered when deriving:
      //
      // d/dx ( p(x)*exp(-x^2/2) ) = (d/dx p(x) - p(x)*x/sigma^2) * exp(-.5*(x/sigma)^2)

      // derive the polynomial part first
      for (int d{1}; d < n; ++d) {
          dH[d - 1] = d*H[d];
      } // d
      dH[n - 1] = 0;

      // now add the terms from the inner derivative of exp(-.5*(x/sigma)^2)
      for (int d{0}; d < n - 1; ++d) {
          dH[d + 1] -= H[d]*pow2(sigmainv); // times -(x/sigma^2)
      } // d

      return (H[n - 1] != 0); // returns 1 if the highest term could not be stored
  } // derive

  template <typename real_t>
  void shift_polynomial_centers(
        real_t c_shifted[] // result: shifted polynomial
      , real_t const c[] // assume p(x) = sum_k=0...nmax-1 c[k] * x^k
      , int const nmax // number of coefficients for both, c[] and c_shifted[]
      , real_t const x_shift // how much is the center shifted
  ) {
      // shift the center of the polynomial by Taylor expansion

      std::vector<real_t> c_old(nmax);
      for (int k{0}; k < nmax; ++k) {
          c_old[k] = c[k]; // get a work copy
      } // k

      double kfactorial{1}; // init kfactorial with 0! == 1
      for (int k{0}; k < nmax; ++k) { // loop MUST run forward from 0

          // evaluate the value of d^k p(x) / d x^k at x=x_shift
          real_t val{0};
          {
              real_t xsp{1}; // x_shift^p
              for (int p{0}; p < nmax - k; ++p) { // we only need to run up to nmax-k as the degree of the input poly is decreased with every k
                  val += xsp * c_old[p];
                  xsp *= x_shift; // update x_shift^p for the next p-iteration
              } // p
          } // ToDo: Horner-scheme could be used

          c_shifted[k] = val / kfactorial;

          // now derive the original polynomial, in-place, for the next k-iteration
          for (int p{1}; p < nmax - k; ++p) { // loop MUST run forward from 1
              c_old[p - 1] = p * c_old[p]; // d/dx x^p = p*x^{p-1}
          } // p
          c_old[nmax - k - 1] = 0;

          kfactorial *= (k + 1); // update kfactorial for the next k-iteration
      } // k
  } // shift_polynomial_centers


  template <typename real_t>
  real_t overlap_of_poly_times_Gauss_with_pure_powers( // return the integral
        real_t const p[] // polynomial coefficients p[i] * x^i, i=0...n0-1
      , int const n0 // number of polynomial coefficients
      , double const sigma0 // Gaussian decay exp(-x^2/(2*s0^2))
      , int const moment // additional power x^moment
  ) {
      assert( 0 <= moment );
      assert( 0 < sigma0 );
      auto const sigma = std::sqrt(2.)*sigma0; // the integrate function assumes exp(-x^2/sigma^2)
      return integrate(p, n0, sigma, moment);
  } // overlap_of_poly_times_Gauss_with_pure_powers


  template <typename real_t>
  real_t overlap_of_two_Hermite_Gauss_functions(
      real_t const H0[], int const n0, double const sigma0, // polynomial H0[j < n0] times Gaussian exp(-x^2/(2*sigma_0^2))
      real_t const H1[], int const n1, double const sigma1, // polynomial H1[i < n1] times Gaussian exp(-x^2/(2*sigma_1^2))
      double const distance, // separating distance of the two centers
      int const moment=0 // multiply a moment x^m to the polynomial before integration
  ) {

      auto const k0 = .5/pow2(sigma0), k1 = .5/pow2(sigma1);
      auto const denom = 1./(k0 + k1);
      auto const sigma = std::sqrt(denom);

      int const n = n0 + n1; // ToDo: check if we can use n0 + n1 - 1 without changing the results
      std::vector<real_t> h0xh1(n); // product of H0s and H1s

      if (0 == distance) {
          multiply(h0xh1.data(), n, H0, n0, H1, n1);
          real_t const value = integrate(h0xh1.data(), n, sigma, moment);
          return value;
      } // distance zero

      auto const sh0 =  distance*k1*denom;
      auto const sh1 = -distance*k0*denom;
      std::vector<real_t> h0s(n0), h1s(n1); // H0 shifted by sh0 and H1 shifted by sh1
      shift_polynomial_centers(h0s.data(), H0, n0, sh0);
      shift_polynomial_centers(h1s.data(), H1, n1, sh1);
      multiply(h0xh1.data(), n, h0s.data(), n0, h1s.data(), n1);
      real_t const value = integrate(h0xh1.data(), n, sigma, moment) * std::exp(-k0*sh0*sh0 -k1*sh1*sh1);
      return value;
  } // overlap_of_two_Hermite_Gauss_functions

  template <int ncut, typename real_t>
  status_t generate_density_tensor(
        real_t tensor[] // data layout [2*ncut-1][ncut][ncut]
      , int const echo=9
      , float const sigma_over_sigmap_squared=2
  ) {
      // this structure can be used to describe the density generation
      // for the density we assume that it is sufficient to
      // represent the density in a SHO basis
      // with sigma_\rho = sigma/sqrt(2) and nu_max_\rho = 2\nu_max

      status_t stat(0);
      if (echo > 1) std::printf("\n\n\n# %s ncut=%d\n", __func__, ncut);
      double const sigma = 1; // typically == 1
      double const sigma_inv = 1./sigma; // typically == 1
      double const sigmapinv2 = sigma_over_sigmap_squared*sigma_inv*sigma_inv; // typically == 2
      double const sigmapinv = std::sqrt(sigmapinv2); // typically == 1.414
      double const alpha = 2*sigma_inv*sigma_inv + sigmapinv2; // == 4
      double const sqrt_alpha_inv = 1./std::sqrt(alpha); // typically == 0.5
      view2D<double> H(ncut, ncut, 0.0), Hp(2*ncut, 2*ncut, 0.0);
      prepare_centered_Hermite_polynomials(H.data(), ncut, sigma_inv); // unit spread sigma=1, L2-normalized
      prepare_centered_Hermite_polynomials(Hp.data(), 2*ncut, sigmapinv); // spread sigma_p = sigma/sqrt(2), L2-normalized
      view3D<real_t> t3(tensor, ncut, ncut); // wrapper
      for (int p{0}; p < 2*ncut - 1; ++p) {
          if (echo > 1) std::printf("\n# p = %d\n", p);
          for (int n{0}; n < ncut; ++n) {
              std::vector<double> HHp(3*ncut, 0.0);
              stat += multiply(HHp.data(), 3*ncut, H[n], ncut, Hp[p], 2*ncut);
              for (int m{0}; m < ncut; ++m) {
                  real_t tensor_value = 0;
                  if (0 == (p + n + m) % 2) { // odd contributions are zero by symmetry
                      std::vector<double> HHpH(4*ncut, 0.0);
                      stat += multiply(HHpH.data(), 4*ncut, H[m], ncut, HHp.data(), 3*ncut);
                      auto const P_pnm = integrate(HHpH.data(), 4*ncut, sqrt_alpha_inv);
//                    if (echo > 1) std::printf(" %d%d%d %.9f\n", p,n,m, P_pnm); // show tensor values as list
                      if (echo > 1) std::printf(" %.9f", P_pnm); // show tensor values
                      // tensor has shape P_pnm[2*ncut-1][ncut][ncut] with each second entry zero
                      tensor_value = P_pnm;
                  } // even?
                  if (tensor) t3(p,n,m) = tensor_value; // store only if output array pointer is non-zero
              } // m
              if (echo > 1) std::printf("\n");
          } // n
          if (echo > 1) std::printf("\n");
      } // p
      return stat; // non-zero if some polynomial coefficients got lost during multiplication
  } // generate_density_tensor


  template <int ncut, typename real_t>
  status_t generate_density_or_potential_tensor(
        real_t tensor[] // data layout [2*ncut-1][ncut][ncut]
      , int const echo=9
      , float const sigma_over_sigmap_squared=2 // 2:typical for density tensor
   ) {
      status_t stat(0);
      // this structure can be used to describe the density generation
      // for the density we assume that it is sufficient to
      // represent the density in a SHO basis
      // with sigma_\rho = sigma/sqrt(2) and nu_max_\rho = 2\nu_max
      double const sigma = 1; // typically == 1
      double const sigma_inv = 1./sigma; // typically == 1
      double const sigmapinv2 = sigma_over_sigmap_squared*pow2(sigma_inv); // typically == 2
      double const sigmapinv = std::sqrt(sigmapinv2); // typically == 1.414
      double const alpha = 2*sigma_inv*sigma_inv + sigmapinv2; // == 4
      double const sqrt_alpha_inv = 1./std::sqrt(alpha); // typically == 0.5
      view2D<double> H(ncut, ncut), Hp(2*ncut, 2*ncut);
      prepare_centered_Hermite_polynomials(H.data(), ncut, sigma_inv); // unit spread sigma=1, L2-normalized
      prepare_centered_Hermite_polynomials(Hp.data(), 2*ncut, sigmapinv); // spread sigma_p = sigma/sqrt(2), L2-normalized
      view3D<real_t> t3(tensor, ncut, ncut); // wrapper
      for (int n{0}; n < ncut; ++n) {
          for (int m{0}; m < ncut; ++m) {
              std::vector<double> HH(2*ncut, 0.0);
              stat += multiply(HH.data(), 2*ncut, H[n], ncut, H[m], ncut);
              for (int p{0}; p < 2*ncut - 1; ++p) {
                  real_t tensor_value = 0;
                  if (0 == (p + n + m) % 2) { // odd contributions are zero by symmetry
                      std::vector<double> HHHp(4*ncut, 0.0);
                      stat += multiply(HHHp.data(), 4*ncut, HH.data(), 2*ncut, Hp[p], 2*ncut);
                      auto const P_pnm = integrate(HHHp.data(), 4*ncut, sqrt_alpha_inv);
                      // tensor has shape P_pnm[2*ncut-1][ncut][ncut] with each second entry zero
                      tensor_value = P_pnm;
                  } // even?
                  if (tensor) t3(p,n,m) = tensor_value; // store only if output array pointer is non-zero
              } // p
          } // m
      } // n
      if (nullptr == tensor) return -1; // function had no effect
      if (echo > 1) {
          std::printf("\n\n\n# %s ncut=%d\n", __func__, ncut);
          for (int p{0}; p < 2*ncut - 1; ++p) {
              std::printf("\n# p = %d\n", p);
              for (int n{0}; n < ncut; ++n) {
                  for (int m{0}; m < ncut; ++m) {
                      if (0 == (p + n + m) % 2) { // odd contributions are zero by symmetry
                          std::printf(" %.9f", t3(p,n,m)); // show tensor values
                      } // even?
                  } // m
                  std::printf("\n");
              } // n
            std::printf("\n");
          } // p
      } // echo
      return stat; // non-zero if some polynomial coefficients got lost during multiplication
  } // generate_density_or_potential_tensor

  template <int ncut, typename real_t>
  status_t product_tensor_plain(
        real_t tensor[] // data layout [2*ncut-1][ncut][ncut]
      , double const sigma=2 // 2:typical for density tensor
      , double const sigma0=1
      , double const sigma1=1
  ) {

    status_t stat(0);
    double const sigma0inv = 1./sigma0;
    double const sigma1inv = 1./sigma1;
    double const sigmapinv = 1./sigma;
    double const alpha = pow2(sigma0inv) + pow2(sigma1inv) + pow2(sigmapinv);
    double const sqrt_alpha_inv = 1./std::sqrt(alpha); // typically == 0.5
    view2D<double> Hp(2*ncut, 2*ncut), H0(ncut, ncut), H1(ncut, ncut);
    prepare_centered_Hermite_polynomials(H0.data(), ncut, sigma0inv); // L2-normalized
    prepare_centered_Hermite_polynomials(H1.data(), ncut, sigma1inv); // L2-normalized
    prepare_centered_Hermite_polynomials(Hp.data(), 2*ncut, sigmapinv); // L2-normalized
    view3D<real_t> t3(tensor, ncut, ncut); // wrapper
    for (int n{0}; n < ncut; ++n) {
        for (int m{0}; m < ncut; ++m) {
            double HH[2*ncut];
            stat += multiply(HH, 2*ncut, H0[n], ncut, H1[m], ncut);
            for (int p{0}; p < 2*ncut - 1; ++p) {
                real_t tensor_value = 0;
                if (0 == (p + n + m) % 2) { // odd contributions are zero by symmetry
                    double HHHp[4*ncut];
                    stat += multiply(HHHp, 4*ncut, HH, 2*ncut, Hp[p], 2*ncut);
                    auto const P_pnm = integrate(HHHp, 4*ncut, sqrt_alpha_inv);
                    // tensor has shape P_pnm[2*ncut-1][ncut][ncut] with each second entry zero
                    tensor_value = P_pnm;
                } // even?
                if (tensor) t3(p,n,m) = tensor_value; // store only if output array pointer is non-zero
            } // p
        } // m
    } // n
    return stat; // non-zero if some polynomial coefficients got lost during multiplication
  } // product_tensor_plain


  template <typename real_t>
  status_t product_tensor(
        real_t tensor[]
      , int const ncut // data layout [2*ncut-1][ncut][ncut]
      , double const sigma  // =2: typical for density tensor
      , double const sigma1 // =1
      , double const sigma0 // =1
  ) {

      status_t stat(0);
      double const sigma0inv = 1./sigma0;
      double const sigma1inv = 1./sigma1;
      double const sigmapinv = 1./sigma;
      double const alpha = pow2(sigma0inv) + pow2(sigma1inv) + pow2(sigmapinv);
      double const sqrt_alpha_inv = 1./std::sqrt(alpha); // typically == 0.5
      view2D<double> H0(ncut, ncut);
      view2D<double> H1(ncut, ncut);
      view2D<double> Hp(2*ncut, 2*ncut);
      prepare_centered_Hermite_polynomials(H0.data(), ncut, sigma0inv); // L2-normalized
      prepare_centered_Hermite_polynomials(H1.data(), ncut, sigma1inv); // L2-normalized
      prepare_centered_Hermite_polynomials(Hp.data(), 2*ncut, sigmapinv); // L2-normalized
      view3D<real_t> t3(tensor, ncut, ncut); // wrapper
      for (int i{0}; i < ncut; ++i) {
          for (int j{0}; j < ncut; ++j) {
              std::vector<double> HH(2*ncut);
              stat += multiply(HH.data(), 2*ncut, H1[i], ncut, H0[j], ncut);
              for (int p{0}; p < 2*ncut - 1; ++p) {
                  real_t tensor_value = 0;
                  if (0 == (p + i + j) % 2) { // odd contributions are zero by symmetry
                      std::vector<double> HHHp(4*ncut);
                      stat += multiply(HHHp.data(), 4*ncut, HH.data(), 2*ncut, Hp[p], 2*ncut);
                      auto const P_pnm = integrate(HHHp.data(), 4*ncut, sqrt_alpha_inv);
                      // tensor has shape P_pnm[2*ncut-1][ncut][ncut] with each second entry zero
                      tensor_value = P_pnm;
                  } // even?
                  if (tensor) t3(p,i,j) = tensor_value; // store only if output array pointer is non-zero
              } // p
          } // j
      } // i
      return stat; // non-zero if some polynomial coefficients got lost during multiplication
  } // product_tensor

  template <typename real_t>
  status_t moment_tensor(
        real_t tensor[] // data layout [1 + max(0, maxmoment)][n1][n0]
      , double const distance
      , int const n1
      , int const n0
      , double const sigma1 // =1
      , double const sigma0 // =1
      , int const maxmoment // default=0 (overlap matrix)
  ) {

      double const sigma0inv = 1./sigma0, sigma1inv = 1./sigma1;
      view2D<double> H0(n0 + 1, n0 + 1); // polynomial coefficients
      view2D<double> H1(n1 + 1, n1 + 1); // polynomial coefficients
      double constexpr normalize = 1; // 1: L2-normalize Hermite polynomials with Gauss metric
      prepare_centered_Hermite_polynomials(H0.data(), n0 + 1, sigma0inv, normalize);
      prepare_centered_Hermite_polynomials(H1.data(), n1 + 1, sigma1inv, normalize);
      view3D<real_t> t3(tensor, n1, n0); // wrapper
      // ToDo: this function will shift each of the polynomials H0 (1 + maxmoment)*n1 times
      //                            and each of the polynomials H1 (1 + maxmoment)*n0 times
      //                            so there could be a lot of saving when the shifted polynomials
      //                            are constructed in advance and the moment loop becomes innermost

      // the underived Hermite polynomial H_n has non-zero coefficients up to 1+n
      if (-2 == maxmoment) {
          // derive Hermite polynomials to form the kinetic energy matrix
          // the Hermite of the derived Hermite Gauss function has one non-zero coefficient more
          std::vector<double> dH0(n0 + 1), dH1(n1 + 1);
          for (int j{0}; j < n0; ++j) {
              derive_Hermite_Gauss_polynomial(dH0.data(), H0[j], n0 + 1, sigma0inv);
              set(H0[j], n0 + 1, dH0.data()); // overwrite with the derived form
          } // j
          for (int i{0}; i < n1; ++i) {
              derive_Hermite_Gauss_polynomial(dH1.data(), H1[i], n1 + 1, sigma1inv);
              set(H1[i], n1 + 1, dH1.data()); // overwrite with the derived form
          } // i
      } // prepare derivatives

      for (int moment{0}; moment <= std::max(0, maxmoment); ++moment) {
          for (int i{0}; i < n1; ++i) {
              for (int j{0}; j < n0; ++j) {
                  t3(moment,i,j) = overlap_of_two_Hermite_Gauss_functions(
                      H0[j], n0 + 1, sigma0, H1[i], n1 + 1, sigma1, distance, moment);
              } // j
          } // i
      } // moment

      // modified to stay within allocated memory bounds, ToDo: check for correctness of the math!

      return 0; // success
  } // moment_tensor


  status_t moment_normalization(
        double matrix[] // output matrix, assumed layout [M][M]
      , int const M // dimensions
      , double const sigma // Gaussian spread
      , int const echo // log-level
  ) {

      // uses LAPACK for inversion although the matrix is already in upper triangular shape
      if (M < 1) return 0;
      view2D<double> mat1D(matrix, M); // wrapper for result, assume data layout [M][M]
      int constexpr check = 1;
      view2D<double> mcopy(check*M, M, 0.0); // get memory if check > 0

      view2D<double> H(M, M, 0.0); // polynomial coefficients
      double constexpr normalize = 0; // 1: L2-normalize Hermite polynomials with Gauss metric
      prepare_centered_Hermite_polynomials(H.data(), H.stride(), 1./sigma, normalize);
      if (echo > 4) std::printf("\n# %s numax=%i sigma=%g %s\n", __func__, M - 1, sigma*Ang, _Ang);
      for (int n{0}; n < M; ++n) {
          if (echo > 4) std::printf("# %s n=%i ", __func__, n);
          for (int moment{0}; moment < M; ++moment) { // square loop
              mat1D(n,moment) = overlap_of_poly_times_Gauss_with_pure_powers(H[n], M, sigma, moment);
              if (check > 0) mcopy(n,moment) = mat1D(n,moment); // store a copy
  //             if ((n & 1) == (moment & 1)) // pairity condition, ToDo: check if analytical parity improves it
              if (echo > 4) std::printf(" %g", mat1D(n,moment));
          } // moment
          if (echo > 4) std::printf("\n");
      } // n = 1D HO quantum number

      return 0;
  } // moment_normalization





  template // explicit template instantiation for double
  status_t moment_tensor(double*, double, int, int, double, double, int);

  template // explicit template instantiation for double
  status_t product_tensor(double*, int, double, double, double);





#ifdef    NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

#if 0
  template <typename real_t>
  void plot_poly(real_t const poly[], int const m, char const *name="") {
      std::printf("# Poly %s :", name);
      printf_vector("  %.6f", poly, m);
  } // plot
#endif // 0

  template <typename real_t>
  real_t eval_poly(real_t const poly[], int const m, double x) {
      double xpow{1}, val{0};
      for (int d{0}; d < m; ++d) {
          val += poly[d] * xpow;
          xpow *= x;
      } // d
      return val;
  } // eval

  status_t test_Hermite_polynomials(int const echo=1, int const ncut=8, double const sigma=1) {
      // see if the first ncut Hermite polynomials are orthogonal and normalized
      view2D<double> H(ncut, ncut, 0.0);
      prepare_centered_Hermite_polynomials(H.data(), ncut, 1./sigma);
      std::vector<double> hh(2*ncut, 0.0); // product polynomial
      int ndev{0}; double mdev{0};
      for (int n{0}; n < ncut; ++n) {
          if (echo > 4) std::printf("# %s n=%d check ortho-normality\n", __func__, n);
          if (echo > 3) {
//            plot_poly(H[n], 1+n, "H_n");
              std::printf("# H_%i: ", n); // no new line
              printf_vector(" %.6f", H[n], 1+n);
          } // echo
          for (int m{0}; m < ncut; ++m) {
              multiply(hh.data(), 2*ncut, H[n], 1+n, H[m], 1+m);
              if (echo > 5) {
//                plot_poly(hh.data(), 1+n+m, "H_n*H_m");
                  std::printf("# H_%i * H_%i: ", n, m); // no new line
                  printf_vector(" %.6f", hh.data(), 1+n+m);
              } // echo
              double const norm = integrate(hh.data(), 1+n+m, sigma);
              mdev = std::max(mdev, std::abs(norm - (m == n)));
              if (echo > 9) std::printf("%9.1e", norm - (m == n));
              ndev += (std::abs(norm - (m == n)) > 1e-10);
          } // m
          if (echo > 1) std::printf("\n");
      } // n
      if (echo > 0) std::printf("# %s: up to %d the largest deviation from Kroecker is %.1e \n", __func__, ncut - 1, mdev);
      return ndev;
  } // test_Hermite_polynomials

  status_t test_Hermite_Gauss_overlap(int const echo=1, int const numerical=999) {
      // show the overlap of the lowest 1D Hermite-Gauss functions as a function of distance
      int constexpr ncut = 4;
      double const sigma0 = 1.3, sigma1 = 0.75; // any two positive real numbers
  //     double const sigma0 = 1, sigma1 = sigma0; // both spreads are the same

      view2D<double> H0(ncut, ncut, 0.0), H1(ncut, ncut, 0.0);
      prepare_centered_Hermite_polynomials(H0.data(), ncut, 1/sigma0);
      prepare_centered_Hermite_polynomials(H1.data(), ncut, 1/sigma1);
      double maxdevall{0};
      if (echo > 3) std::printf("\n# %s sigma0=%g sigma1=%g %s\n", __func__, sigma0*Ang, sigma1*Ang, _Ang);
      for (auto dist{0.0}; dist < 11; dist += .1) {
          if (echo > 4) std::printf("# %s distance=%.3f  ", __func__, dist*Ang);
          double maxdev{0};
          for (int n{0}; n < ncut; ++n) {
              for (int m{0}; m < ncut; ++m) {
                  double const ovl = overlap_of_two_Hermite_Gauss_functions(H0[n], 1+n, sigma0,
                                                                            H1[m], 1+m, sigma1, dist);
                  if (echo > 4) std::printf(" %.6f", ovl);
                  if (numerical > 0) {
                      double const dx = 7.0/numerical;
                      double ovl_numerical = 0;
                      for (int ix{-numerical}; ix <= numerical + dist/dx; ++ix) {
                          double const x0 = ix*dx - 0; // position at zero
                          double const x1 = ix*dx - dist; // position at (plus) dist
                          ovl_numerical += eval_poly(H0[n], 1+n, x0) * std::exp(-0.5*pow2(x0/sigma0))
                                         * eval_poly(H1[m], 1+m, x1) * std::exp(-0.5*pow2(x1/sigma1));
                      } // ix
                      ovl_numerical *= dx;
                      if (echo > 5) std::printf(" %.6f", ovl_numerical);
                      maxdev = std::max(maxdev, std::abs(ovl - ovl_numerical));
                  } // numerical
              } // m
          } // n
          if (echo > 4) std::printf("\n");
          if (numerical > 0 && echo > 3) std::printf("# %s max deviation for distance=%.3f %s is %.1e\n", __func__, dist*Ang, _Ang, maxdev);
          maxdevall = std::max(maxdevall, maxdev);
      } // dist
      if (numerical > 0 && echo > 2) std::printf("# %s max deviation is %.1e\n", __func__, maxdevall);
      return (maxdevall > 1e-12); // only a valid automatized check if numerical > 0
  } // test_Hermite_Gauss_overlap

  status_t test_kinetic_overlap(int const echo=2) {
      // show the kinetic energy of the lowest 1D Hermite-Gauss functions as a function of distance
      // test if the derivation operator can be cast to any side
      // --> yes if sigma1 == sigma0, otherwise it breaks
      // --> we should use the first derivative applied to left and right for the kinetic energy
      int constexpr ncut = 6;
      int constexpr mcut = ncut - 2;
//    double const sigma0 = 1, sigma1 = sigma0; // same sigma
      double const sigma0 = 0.9, sigma1 = 1.11; // seems ok for different sigmas
      view2D<double> H0(ncut, ncut, 0.0), H1(ncut, ncut, 0.0);
      prepare_centered_Hermite_polynomials(H0.data(), ncut, 1./sigma0);
      prepare_centered_Hermite_polynomials(H1.data(), ncut, 1./sigma1);

      view2D<double> dH0(ncut, ncut, 0.0), dH1(ncut, ncut, 0.0), d2H0(ncut, ncut, 0.0), d2H1(ncut, ncut, 0.0);
      for (int n{0}; n < mcut; ++n) {
          // first derivatives
          derive_Hermite_Gauss_polynomial(dH0[n], H0[n], ncut, 1./sigma0);
          derive_Hermite_Gauss_polynomial(dH1[n], H1[n], ncut, 1./sigma1);
          // second derivatives
          derive_Hermite_Gauss_polynomial(d2H0[n], dH0[n], ncut, 1./sigma0);
          derive_Hermite_Gauss_polynomial(d2H1[n], dH1[n], ncut, 1./sigma1);
      } // n
      double maxdev1{0}, maxdev2{0}, maxdev3{0};
      if (echo > 4) std::printf("# %s  distance overlaps\n", __func__);
      for (auto dist{0.0}; dist < 11; dist += .01) {
          if (echo > 4) std::printf("%.3f", dist);
          for (int n{0}; n < mcut; ++n) {
              for (int m{0}; m < mcut; ++m) {
                  auto const d2d0 = overlap_of_two_Hermite_Gauss_functions(d2H0[n], ncut, sigma0, H1[m], ncut, sigma1, dist);
                  auto const d0d2 = overlap_of_two_Hermite_Gauss_functions(H0[n], ncut, sigma0, d2H1[m], ncut, sigma1, dist);
                  auto const d1d1 = overlap_of_two_Hermite_Gauss_functions(dH0[n], ncut, sigma0, dH1[m], ncut, sigma1, dist);
  //                 auto const ovl  = overlap_of_two_Hermite_Gauss_functions(H0[n], ncut, sigma0, H1[m], ncut, sigma1, dist);
  //                 if (echo > 1) std::printf(" %.9f", ovl); // show overlap
  //              if (echo > 1) std::printf("  %.9f %.9f %.9f", d2d0, d0d2, -d1d1); // show 3 values
  //              if (echo > 1) std::printf("  %.1e %.1e %.1e", d2d0 + d1d1, d0d2 + d1d1, d2d0 - d0d2); // show deviations
                  if (echo > 6) std::printf(" %.9f", -d1d1); // show 1 value
                  auto const d2avg = .5*d2d0 + .5*d0d2;
                  if (echo > 8) std::printf("  %.9f %.9f", d2avg, -d1d1); // show 2 values
                  maxdev3 = std::max(maxdev3, std::abs(d2avg + d1d1)); // one order better than dev1 and dev2
                  maxdev2 = std::max(maxdev2, std::abs(d2d0 - d0d2));
                  maxdev1 = std::max(maxdev1, std::abs(d2d0 + d1d1));
                  maxdev1 = std::max(maxdev1, std::abs(d0d2 + d1d1));
              } // m
          } // n
          if (echo > 4) std::printf("\n");
      } // dist
      if (echo > 0) std::printf("\n# %s deviations %g, %g and %g\n", __func__, maxdev1, maxdev2, maxdev3);
      return (maxdev3 > 2e-14);
  } // test_kinetic_overlap

  status_t test_density_or_potential_tensor(int const echo=2) {
      int constexpr ncut = 8, n = (2*ncut - 1)*ncut*ncut;
      std::vector<double> t(n), tp(n), tt(n);
      double df_max{0};
      float const ssp2_min = 1.f, ssp2_max = 3.f, ssp2_inc = 1.01f;
      for (float ssp2{ssp2_min}; ssp2 < ssp2_max; ssp2 *= ssp2_inc) {
          generate_density_tensor<ncut>(t.data(), 0, ssp2); // reference implementation
//           product_tensor_plain<ncut>(tt.data(), 1./std::sqrt(ssp2)); // old implementation
          product_tensor(tt.data(), ncut, 1./std::sqrt(ssp2)); // new implementation
          auto const & ts = tt;
          generate_density_or_potential_tensor<ncut>(tp.data(), 0, ssp2);
//           auto const & ts = tp;
          double df{0};
          for (int i{0}; i < n; ++i) {
              auto const ab = std::abs(t[i] - ts[i]);
              if ((ab > 1e-14) && (echo > 7)) std::printf("# %s deviations in element [%d] %g\n", __func__, i, ab);
              df = std::max(df, ab);
          } // i
          if (echo > 3) std::printf("# %s (%g) deviations %g\n", __func__, ssp2, df);
          df_max = std::max(df, df_max);
      } // ssp2
      if (echo > 0) std::printf("\n# %s (%.2f ... %.1f %% ... %.2f) largest deviation %g\n",
                    __func__, ssp2_min, (ssp2_inc - 1)*100, ssp2_max, df_max);
      if (echo > 2) generate_density_or_potential_tensor<ncut>(tp.data(), echo); // default ssp2=2
      return (df_max > 8e-9); // return error of the deviations are too strong
  } // test


  typedef std::complex<double> complex_t;
  extern "C" {
      // complex<double> hermitian generalized eigenvalue problem
      void zhegv_(int const*, char const*, char const*, int const*,
                  complex_t*, int const*, complex_t*, int const*,
                  double*, complex_t*, int const*, double*, int*);
      // complex<double> hermitian eigenvalue problem
      void zheev_(char const*, char const*, int const*, complex_t*,
                  int const*, double*, complex_t*, int const*, double*, int*);
  } // LAPACK

  status_t test_simple_crystal(int const echo=3) {
      auto const a0 = control::get("sho_overlap.lattice.constant", 8.0);
      if (echo > 0) std::printf("\n# %s\n", __func__);
      typedef vector_math::vec<3,double> vec3;
      typedef vector_math::vec<3,int>    vec3i;
      auto const numax = int(control::get("sho_overlap.crystal.numax", 4.));
      int const ncut = numax + 2;

      vec3 cv[3], bv[3]; // vectors of the cell and the Bravais matrix
      { // scope: lattice structure
          auto const structure_abbrev = control::get("sho_overlap.crystal.structure", "fcc"); // choice {sc, bcc, fcc}
          char const structure = structure_abbrev[0]; // only the 1st char counts
          if (echo > 1) std::printf("# %s %s (%c)\n", __func__, structure_abbrev, structure);
          assert('s' == structure || 'b' == structure || 'f' == structure);
          double const recip = (2*constants::pi)/a0;
          for (int d{0}; d < 3; ++d) {
              if ('s' == structure) { // simple-cubic
                  cv[d] = 0; cv[d][d] = a0;
                  bv[d] = 0; bv[d][d] = recip; // sc
              } else {
                  bool const bcc = ('b' == structure); // body-centered cubic
                  cv[d] = .5*a0; cv[d][d] *= bcc ? 0 : -1;
                  bv[d] = recip; bv[d][d] *= bcc ? -1 : 0;
              }
              if (echo > 4) std::printf("# cell %8.3f%8.3f%8.3f %s \treci %8.3f%8.3f%8.3f a.u.\n",
                  cv[d][0]*Ang, cv[d][1]*Ang, cv[d][2]*Ang,_Ang,  bv[d][0], bv[d][1], bv[d][2]);
          } // d
          if (echo > 7) {
              std::printf("\n# check that cell times reciprocal basis are a unit matrix\n");
              for (int i{0}; i < 3; ++i) {
                  std::printf("# cell * recip  ");
                  for (int j{0}; j < 3; ++j) std::printf(" %g", dot(cv[i], bv[j])/(2*constants::pi));
                  std::printf("  recip * cell\t");
                  for (int j{0}; j < 3; ++j) std::printf(" %g", dot(bv[i], cv[j])/(2*constants::pi));
                  std::printf("\n");
              } // i
              std::printf("\n");
          } // echo
      } // scope: lattice structure

      double shortest_bond2{9e99};
      for (int i3{-1}; i3 <= 1; ++i3) {
          for (int i2{-1}; i2 <= 1; ++i2) {
              for (int i1{-1}; i1 <= 1; ++i1) {
                  vec3 const pos = cv[0]*i1 + cv[1]*i2 + cv[2]*i3; // assume one atom per unit cell
                  double const d2 = norm(pos);
                  if (d2 > 0) shortest_bond2 = std::min(shortest_bond2, d2);
              } // i1
          } // i2
      } // i3
      double const shortest_bond = std::sqrt(shortest_bond2);

      if (echo > 0) std::printf("# shortest bond is %g Bohr\n", shortest_bond);

      bool const overlap_eigvals = false;
      // choose the return radius as a fraction of shortest_bond length
      double const sigma = .75*shortest_bond/std::sqrt(2.*numax + 3.),
                  sigma0 = sigma, sigma1 = sigma;
      if (echo > 0) std::printf("# SHO up to numax=%d, spread sigma = %.9f Bohr\n", numax, sigma);

      // return radius of the classical harmonic oscillator
      double const return_radius = sigma*std::sqrt(2.*numax + 3.);
      if (echo > 0) std::printf("# classical return radius at %g Bohr\n", return_radius);

      double const dmax = 12*sigma; // 12 sigma is converged for fcc

      double const normalize = 0; // 0:do not normalize, we have to deal with an overlap matrix anyway
      view2D<double> H0(ncut, ncut, 0.0), H1(ncut, ncut, 0.0);
      prepare_centered_Hermite_polynomials(H0.data(), ncut, 1./sigma0, normalize);
      prepare_centered_Hermite_polynomials(H1.data(), ncut, 1./sigma1, normalize);

      view2D<double> dH0(ncut, ncut, 0.0), dH1(ncut, ncut, 0.0);
      for (int n{0}; n < ncut; ++n) {
          // show the Hermite polynomial coefficients for H0
          if (echo > 3) {
              std::printf("# H[%x]: ", n);
              for (int m{0}; m <= n; ++m) {
                  std::printf("%8.4f", H0[n][m]);
              }   std::printf("\n");
          } // echo

          // construct first derivatives
          derive_Hermite_Gauss_polynomial(dH0[n], H0[n], ncut, 1./sigma0);
          derive_Hermite_Gauss_polynomial(dH1[n], H1[n], ncut, 1./sigma1);
      } // n

      int const n3D = sho_tools::nSHO(numax);
      if (echo > 5) {
          std::printf("# %d SHO functions up to numax=%d\n", n3D, numax);
          {   std::printf("# list %d SHO functions: ", n3D);
              for (int n0{0}; n0 <= numax; ++n0) {
                  for (int n1{0}; n1 <= numax - n0; ++n1) {
                      for (int n2{0}; n2 <= numax - n0 - n1; ++n2) {
                          std::printf("%x%x%x ", n0,n1,n2);
                      } // 2n2
                  } // n1
              } // n0
              std::printf("\n");
          } // scope
      } // echo

      bool const DoS = (0 < control::get("sho_overlap.test.DoS", 0.)); // 1: density of states, 0: bandstructure
      if (DoS && echo > 3) std::printf("# compute density of states (DoS)\n");
      bool const Ref = (0 < control::get("sho_overlap.test.Ref", 0.));
      if (Ref && echo > 3) std::printf("# show free electron parabolas as reference\n");
      int imx_ref = 9; if(Ref) imx_ref = int(control::get("sho_overlap.test.Ref.imx", 9.));

      vec3i const imax = std::ceil(dmax/a0);
      int const max_npi = 16*imax[2]*imax[1]*imax[0];
      if (echo > 2) std::printf("# assume at most %d periodic images up to %.3f Bohr\n", max_npi, dmax);
      view4D<double> mat(max_npi, 2, n3D, n3D);
      view2D<int> vpi(max_npi, 4); // periodic image shift vectors
      int npi{0};
      for (int i3{-imax[2]}; i3 <= imax[2]; ++i3) {
      for (int i2{-imax[1]}; i2 <= imax[1]; ++i2) {
      for (int i1{-imax[0]}; i1 <= imax[0]; ++i1) {
          vec3 const pos = cv[0]*i1 + cv[1]*i2 + cv[2]*i3;

          if (!Ref && norm(pos) < dmax*dmax) {
              if (echo > 9) std::printf("%f %f %f\n", pos[0],pos[1],pos[2]);
              int in{0};
              for (int n2{0}; n2 <= numax; ++n2) {
              for (int n1{0}; n1 <= numax - n2; ++n1) {
              for (int n0{0}; n0 <= numax - n2 - n1; ++n0) {
                  int const nv[] = {n0, n1, n2};
                  int im{0};
                  for (int m2{0}; m2 <= numax; ++m2) {
                  for (int m1{0}; m1 <= numax - m2; ++m1) {
                  for (int m0{0}; m0 <= numax - m2 - m1; ++m0) {
                      int const mv[] = {m0, m1, m2};
                      double ovl[3], lap[3];
                      // ToDo: overlap_of_two_Hermite_Gauss_functions
                      //       is called many more times than necessary
                      //       and the max. length of non-zero polynomial coefficients
                      //       can be shorter than ncut in many cases
                      for (int dir{0}; dir < 3; ++dir) {
                          ovl[dir] = overlap_of_two_Hermite_Gauss_functions(
                                          H0[nv[dir]], ncut, sigma0,
                                          H1[mv[dir]], ncut, sigma1, pos[dir]);
                          lap[dir] = overlap_of_two_Hermite_Gauss_functions(
                                        dH0[nv[dir]], ncut, sigma0,
                                        dH1[mv[dir]], ncut, sigma1, pos[dir]);
                      } // dir
                      double const o3D = ovl[0]*ovl[1]*ovl[2];
                      double const l3D = lap[0]*ovl[1]*ovl[2]
                                      + ovl[0]*lap[1]*ovl[2]
                                      + ovl[0]*ovl[1]*lap[2];
                      assert(n3D >= in);
                      assert(n3D >= im);
                      mat(npi,0,in,im) = o3D;
                      mat(npi,1,in,im) = l3D*0.5; // Hartree energy units
                      ++im;
                  }}} // m
                  ++in;
              }}} // n
              vpi(npi,0) = i1; vpi(npi,1) = i2; vpi(npi,2) = i3;
              ++npi; // count periodic images
              assert(max_npi >= npi);
          } // pos inside sphere
      }}} // i1 i2 i3
      if (echo > 2) std::printf("# account for %d periodic images up to %.3f Bohr\n", npi, dmax);
      int const num_periodic_images = npi;

      double smallest_eigval{9e99}, largest_eigval{-9e99};
      vec3 kv_smallest{-9};

      int const lwork = n3D*n3D;
      view2D<complex_t> ovl_mat(n3D, n3D), kin_mat(n3D, n3D);
      std::vector<complex_t> work(lwork);
      std::vector<double> rwork(lwork), eigvals(n3D);
      auto const jobz = 'n', uplo = 'u', jobv = 'v';

      std::vector<std::array<double,4>> kps;
      int diagonalization_failed{0};

      int const num_bins = 1 << 19;
  //     double const inv_bin_width = 27211.4; // 10 mili-electron-Volt
  //     double const bin_width = 1./inv_bin_width;
      double const bin_width = eV*std::min(5e-4, 0.01/eV);
      double const inv_bin_width = 1./bin_width;
      double const energy_offset = -1.5*bin_width;
      int ibin_out_of_range{0};
      std::vector<double> dos;

      if (DoS) {
          dos.assign(num_bins, 0); // allocate and clear
          // create a k-point set with weights
          auto const nkp_sampling = int(control::get("sho_overlap.kmesh.sampling", 2.)); // this is N, use a 2N x 2N x 2N k-point set
          double const inv_kp_sampling = 0.5/nkp_sampling;
          double w8sum{0};
          double const weight = 1;
          vec3 kvec;
          for (int iz{0}; iz < nkp_sampling; ++iz) {   kvec[2] = inv_kp_sampling*(iz + 0.25);
              for (int iy{0}; iy <= iz;      ++iy) {   kvec[1] = inv_kp_sampling*(iy + 0.25);
                  for (int ix{0}; ix <= iy;  ++ix) {   kvec[0] = inv_kp_sampling*(ix + 0.25);
                      kps.push_back({kvec[0], kvec[1], kvec[2], weight});
                      w8sum += weight;
                      if (echo > 8) std::printf("# new k-point %g %g %g weight %g\n", kvec[0], kvec[1], kvec[2], weight);
                  } // ix
              } // iy
          } // iz
          if (echo > 1) std::printf("# %ld k-points in the irriducible Brillouin zone, weight sum = %g\n", kps.size(), w8sum);
          double const w8scale = 1./w8sum; for (size_t ikp{0}; ikp < kps.size(); ++ikp) kps[ikp][3] *= w8scale; // rescale
      } else {
          int const nedges = 6;
          auto const sampling_density = control::get("sho_overlap.kpath.sampling", 1./32);
          double const kpath[nedges][3] = {{.0,.0,.0}, {.5,.0,.0}, {.5,.5,.0}, {.0,.0,.0}, {.5,.5,.5}, {.5,.5,.0}};
          double path_progress{0};
          for (int edge{0}; edge < nedges; ++edge) {
              int const e0 = edge % nedges, e1 = (edge + 1) % nedges;
              vec3 const v0 = kpath[e0], v1 = kpath[e1];
              vec3 const true_kdiff = bv[0]*(v1[0] - v0[0])
                                    + bv[1]*(v1[1] - v0[1])
                                    + bv[2]*(v1[2] - v0[2]);
              double const edge_length = std::sqrt(norm(true_kdiff));

              int const sampling = std::ceil(edge_length/sampling_density);
              double const frac = 1./sampling;
              if (echo > 1) std::printf("# k-point %.6f %.6f %.6f\n", v0[0],v0[1],v0[2]);
              for (int step{0}; step < sampling + (edge == (nedges - 1)); ++step) {
                  double const path_progress_edge = path_progress + (step*frac)*edge_length;
                  vec3 const kvec = v0 + (v1 - v0)*(step*frac);
                  kps.push_back({kvec[0], kvec[1], kvec[2], path_progress_edge});
              } // step
              if (echo > 2) std::printf("# k-point %.6f %.6f %.6f\n", v1[0],v1[1],v1[2]);
              path_progress += edge_length;
          } // edge
      } // DoS

      float progress_percent{.02}; // show first at 2%
      for (unsigned ik{0}; ik < kps.size(); ++ik) {
          vec3 kvec = &(kps[ik][0]); // copy first 3 doubles at this pointer
          vec3 const true_kv = bv[0]*kvec[0] + bv[1]*kvec[1] + bv[2]*kvec[2];

          int info{0};
          if (Ref) {

              int const imx = imx_ref;
              std::vector<double> free_E; free_E.reserve(9*pow3(imx));
              for (int iz{-imx}; iz <= imx; ++iz) {
                  for (int iy{-imx}; iy <= imx; ++iy) {
                      for (int ix{-imx}; ix <= imx; ++ix) {
                          auto const true_kv = bv[0]*(kvec[0] + ix)
                                            + bv[1]*(kvec[1] + iy)
                                            + bv[2]*(kvec[2] + iz);
                          free_E.push_back(0.5*norm(true_kv)); // energy parabolasin Hatree units
                      } // ix
                  } // iy
              } // iz
              std::sort(free_E.begin(), free_E.end()); // sort in-place, ascending
              set(eigvals.data(), n3D, free_E.data()); // copy lowest eigenvals

          } else { // Ref

              // clear matrices
              for (int in{0}; in < n3D; ++in) {
                  for (int im{0}; im < n3D; ++im) {
                      ovl_mat(in,im) = 0;
                      kin_mat(in,im) = 0;
                  } // im
              } // in

              for (int ipi{0}; ipi < num_periodic_images; ++ipi) {
                  vec3 const ipos = vpi[ipi];
                  complex_t const bloch_factor = std::polar(1.0, 2*constants::pi * dot(kvec, ipos));
                  if (echo > 9) std::printf("# periodic image%4d%4d%4d  Bloch-phase = %f + i %f\n",
                      vpi(ipi,0), vpi(ipi,1), vpi(ipi,2), bloch_factor.real(), bloch_factor.imag());
                  // add to matrixes
                  for (int in{0}; in < n3D; ++in) {
                      for (int im{0}; im < n3D; ++im) {
                          ovl_mat(in,im) += bloch_factor*mat(ipi,0,in,im);
                          kin_mat(in,im) += bloch_factor*mat(ipi,1,in,im);
                      } // im
                  } // in
              } // ipi

#if 0
              // check if matrices are hermitian
              auto const threshold = 1e-5;
              for (auto m{ovl_mat}; m == ovl_mat || m == kin_mat; m += (kin_mat - ovl_mat)) {
                  for (int in{0}; in < n3D; ++in) {
                      for (int im{0}; im < in; ++im) {
                          assert(std::abs(m(in,im).real() - m(im,in).real()) < threshold);
                          assert(std::abs(m(in,im).imag() + m(im,in).imag()) < threshold);
                      } // im
                      assert(std::abs(m(in,in).imag()) < threshold);
                  } // in
              } // m
#endif // 0
              // LAPACK call (Fortran77 interface);
              if (overlap_eigvals) {

                  // get the eigenvalues of the overlap operator only
                  zheev_(&jobv, &uplo, &n3D, ovl_mat.data(), &n3D,
                        eigvals.data(), work.data(), &lwork, rwork.data(), &info);
  //                 info = linear_algebra::eigenvalues(n3D,
  //                               kin_mat.data(), kin_mat.stride(), eigvals.data());
#if 0
                  // DEBUG
                  if (0 == info && eigvals[0] < .00315) {
                    std::printf("# lowest eigenvector ");
                    for (int i3D{0}; i3D < n3D; ++i3D) {
                        auto const c = ovl_mat(0,i3D);
                        std::printf("(%.9f,%.9f) ", c.real(), c.imag());
                    }  std::printf("\n");
                  } // DEBUG
#endif // 0
              } else { // overlap_eigvals

                  // solve generalized eigenvalue problem kin_mat*X == diag*ovl_mat*X
  //                 info = linear_algebra::generalized_eigenvalues(n3D,
  //                               kin_mat.data(), kin_mat.stride(),
  //                               ovl_mat.data(), ovl_mat.stride(), eigvals.data());
                  int const itype = 1;
                  zhegv_(&itype, &jobz, &uplo, &n3D, kin_mat.data(), &n3D, ovl_mat.data(), &n3D,
                        eigvals.data(), work.data(), &lwork, rwork.data(), &info);

              } // overlap_eigvals
          } // Ref

          for (int i3D{0}; i3D < n3D; ++i3D) {
              if (eigvals[i3D] < smallest_eigval) {
                  kv_smallest = kvec;
              } // store where the smallest eigenvalue was found
              smallest_eigval = std::min(smallest_eigval, eigvals[i3D]);
              largest_eigval  = std::max(largest_eigval,  eigvals[i3D]);
          } // i3D

          if (0 == info) {
              if (DoS) {
                  double const w8 = kps[ik][3];
                  // accumulate the density of states
                  for (int i3D{0}; i3D < n3D; ++i3D) {
                      double const E = eigvals[i3D];
                      double fbin = (E - energy_offset)*inv_bin_width;
                      int const ibin = (int)std::floor(fbin);
                      if ((ibin < num_bins - 1) && (ibin >= 0)) {
                          // linear interpolation
                          fbin -= ibin;
                          double const w1 = fbin*w8, w0 = (1 - fbin)*w8;
                          dos[ibin + 0] += w0;
                          dos[ibin + 1] += w1;
                      } else { // ibin in range
                          ++ibin_out_of_range;
                      } // ibin in range
                  } // i3D

              } else if(echo > 1) {
                  // show the bandstructure
                  double const path_progress_edge = kps[ik][3];
                  std::printf("%.6f ", path_progress_edge); // abscissa
                  for (int i3D{0}; i3D < n3D; ++i3D) {
                      std::printf("%g ", eigvals[i3D]);
                  } // i3D
                  std::printf("\n");
              } // echo
          } else {
              ++diagonalization_failed;
              if (echo > 2) std::printf("# ik=%i diagonalization failed, info = %d\n", ik, info);
          } // info

          if (progress_percent*kps.size() < ik) {
              if (echo > 3) {
                  std::printf("# progress = %.1f %%\n", ik/(.01*kps.size()));
                  fflush(stdout); // if we do not flush the output, it will be buffered an the progress report makes no sense
              }
              progress_percent = .1*std::ceil(ik/(.1*kps.size()));
          } // show percentage
      } // ik

      if (DoS) {
          if (echo > 2) {
              double const per_eV = inv_bin_width/eV;
              std::printf("\n## energy (%s), integrated DoS, density of states (%s^-1)\n", _eV, _eV);
              double dos_sum{0};
              std::printf("%.6f %g %g\n", 0., 0., 0.); // first entry
              for (int ibin{0}; ibin < num_bins; ++ibin) { // loop-carried dependency on dos_sum
                  double const Ebin = ibin*bin_width + energy_offset;
                  if (dos[ibin] > 0) {
                      dos_sum += dos[ibin];
                      std::printf("%.6f %g %g\n", Ebin*eV, dos_sum, dos[ibin]*per_eV);
                  } else {
                      dos[ibin] = 0;
                  }
              } // ibin
              std::printf("\n# Integrated density of states is %g\n\n", dos_sum);
          } // echo
          if (ibin_out_of_range > 0) warn("%d bin entries were out of range!", ibin_out_of_range);
      } // DoS

      if (echo > 1) std::printf("# diagonalized %d x %d Hamiltonian for %ld k-points\n", n3D, n3D, kps.size());

      if (diagonalization_failed > 0) {
          if (echo > 0) warn("%d diagonalizations failed!", diagonalization_failed);
      } else {
          if (echo > 1) std::printf("\n# smallest and largest eigenvalue%s are %g and %g\n",
              overlap_eigvals?" of the overlap operator":"", smallest_eigval, largest_eigval);
          if (echo > 1) std::printf("# smallest eigenvalue at kvec  %.6f %.6f %.6f\n", kv_smallest[0],kv_smallest[1],kv_smallest[2]);
      }
      return diagonalization_failed;
  } // test_simple_crystal

  status_t test_shifted_polynomial(int const echo=5) {
      status_t stat(0);
      int constexpr M = 8;
      double original[M], shifted[M];
      for (int it{10}; it-- > 0;) {
          for (int d{0}; d < M; ++d) {
              original[d] = simple_math::random(-1., 1.);
          } // d
          double const shift = simple_math::random(-3., 3.);
          shift_polynomial_centers(shifted, original, M, shift);
          double const Ps0 = eval_poly(shifted, M, 0);
          double const PoS = eval_poly(original, M, shift);
          if (echo > 7) std::printf("\n# %s P_shifted(0)=%g\n# %s P_ori(shift)=%g\n", __func__, Ps0, __func__, PoS);
          double const PsS = eval_poly(shifted, M, -shift);
          double const Po0 = eval_poly(original, M, 0);
          if (echo > 7) std::printf("# %s P_shifted(-s)=%g\n# %s P_original(0)=%g\n", __func__, PsS, __func__, Po0);
          stat += (std::abs(Ps0 - PoS) > 1e-9) + (std::abs(PsS - Po0) > 1e-9);
      } // it
      return stat;
  } // test_shifted_polynomial

  status_t test_pure_power_overlap(int const echo=1, int const numerical=999) {
      // show the overlap of the lowest 1D Hermite-Gauss functions with pure powers x^n
      status_t stat(0);
      int constexpr ncut = 8;
      view2D<double> Hs(ncut, ncut);
  //  view2D<double> H1(ncut, ncut);

      double maxreldevall{0};
      for (int isigma{-3}; isigma <= 3; ++isigma) {
          double maxreldev{0};
          double const sigma = std::pow(1.1, isigma);
          if (echo > 3) std::printf("# %s sigma = %g\n", __func__, sigma);
          double constexpr normalize = 0;
          prepare_centered_Hermite_polynomials(Hs.data(), Hs.stride(), 1./sigma, normalize);
  //      prepare_centered_Hermite_polynomials(H1.data(), H1.stride(), 1., normalize);
          for (int n{0}; n < ncut; ++n) {
              for (int moment{0}; moment < ncut; ++moment) { // pure powers x^m
                  double const ovl  = overlap_of_poly_times_Gauss_with_pure_powers(Hs[n], 1+n, sigma, moment);
  //              double const ovl1 = overlap_of_poly_times_Gauss_with_pure_powers(H1[n], 1+n, 1.0, moment) * intpow(sigma, moment);
                  if ((n & 1) ^ (moment & 1)) {
                      assert( 0 == ovl ); // results must be zero if n and moment have different parity!
                  } else {
                      if (echo > 4) std::printf(" %.6f", ovl);
  //                  if (echo > 4) std::printf(" %.6f", ovl1); // ToDo: ovl1 seems to differ from ovl
                      if (numerical > 0) {
                          double const dx = 9.0*sigma/numerical;
                          double ovl_numerical{0};
                          for (int ix{-numerical}; ix <= numerical; ++ix) {
                              double const x = ix*dx; // centered at zero
                              ovl_numerical += eval_poly(Hs[n], 1+n, x) * std::exp(-0.5*pow2(x/sigma)) * std::pow(x, moment);
                          } // ix
                          ovl_numerical *= dx;
                          if (echo > 5) std::printf(" %.6f", ovl_numerical);
                          auto const dev = ovl - ovl_numerical;
                          auto const absdev = std::abs(dev);
                          auto const ref = std::max(std::abs(ovl), std::abs(ovl_numerical));
                          maxreldev = std::max(maxreldev, absdev/ref);
                      } // numerical
                  } // even odd
              } // moment
              if (echo > 1) std::printf("\n");
          } // n
          maxreldevall = std::max(maxreldevall, maxreldev);
          if (numerical > 0 && echo > 2) std::printf("# %s max relative deviation for sigma= %g is %.1e\n", __func__, sigma, maxreldev);
          if (echo > 5) std::printf("\n");
          stat += (maxreldev > 2e-10);
      } // isigma
      if (numerical > 0 && echo > 1) std::printf("# %s max relative deviation between analytical and numerical is %.1e\n", __func__, maxreldevall);
      return stat;
  } // test_pure_power_overlap

  status_t test_moment_normalization(int const echo=1, int const m=8, int const numerical=999) {
      status_t stat(0);
      view2D<double> imat(m, m, 0.0);
      stat += moment_normalization(imat.data(), imat.stride(), 1.5, echo);
      if (echo < 4) return stat;
      for (int i{0}; i < m; ++i) {
          std::printf("%s# %s %i ", (0 == i)?"\n":"", __func__, i);
//        for (int j{i & 1}; j < m; j += 2) // only even/odd terms
          for (int j{0}; j < m; ++j) { // only even/odd terms
              std::printf(" %g", imat(i,j));
          } // j
          std::printf("\n");
      } // i
      // ToDo: can we derive a recursion relation for the imat coefficients?
      return stat; // will always return 0 as we do not invert any longer
  } // test_moment_normalization

  status_t all_tests(int const echo) {
      int n{0}; auto const t = int(control::get("sho_overlap.select.test", -1.)); // -1:all
      status_t stat(0);
      if (t & (1 << n++)) stat += test_moment_normalization(echo);
      if (t & (1 << n++)) stat += test_pure_power_overlap(echo);
      if (t & (1 << n++)) stat += test_shifted_polynomial(echo);
      if (t & (1 << n++)) stat += test_Hermite_polynomials(echo);
      if (t & (1 << n++)) stat += test_Hermite_Gauss_overlap(echo);
      if (t & (1 << n++)) stat += test_kinetic_overlap(echo);
      if (t & (1 << n++)) stat += test_density_or_potential_tensor(echo);
      if (t & (1 << n++)) stat += test_simple_crystal(echo); // expensive
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace sho_overlap
