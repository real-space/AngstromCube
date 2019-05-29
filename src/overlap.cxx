#include <cstdio> // printf
#include <cstdlib> // abs
#include <cmath> // sqrt
#include <algorithm> // max
#include <complex> // std::complex<real_t>
#include <complex>
#include <cmath>
#include <cassert>
 
#include "vec.hxx" // vector_math from exafmm

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

  double constexpr C_PI = 3.14159265358979323846; // pi
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
          value += p[2*d] * kern; // access only the even terms of p[]
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

  status_t test_kinetic_overlap(int const echo=1) {
    // show the kinetic energy of the lowest 1D Hermite-Gauss functions as a function of distance
    // test if the derivation operator can be cast to any side
    // --> yes if sigma1 == sigma0, otherwise it breaks
    // --> we should use the first derivative applied to left and right for the kinetic energy
    int constexpr ncut = 6;
    int constexpr mcut = ncut - 2;
    double const sigma0 = 1, sigma1 = sigma0 + .1; // fails for different sigmas
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
    double maxdev1 = 0, maxdev2 = 0, maxdev3 = 0;
    for(auto dist = 0.0; dist < 11; dist += .01) {
        if (echo > 1) printf("# %s  distance=%.3f    ", __func__, dist);
        for(int n = 0; n < mcut; ++n) {
            for(int m = 0; m < mcut; ++m) {
                auto const d2d0 = overlap_of_two_Hermite_Gauss_functions(&d2H0[ncut*n], ncut, sigma0, &H1[ncut*m], ncut, sigma1, dist);
                auto const d0d2 = overlap_of_two_Hermite_Gauss_functions(&H0[ncut*n], ncut, sigma0, &d2H1[ncut*m], ncut, sigma1, dist);
                auto const d1d1 = overlap_of_two_Hermite_Gauss_functions(&dH0[ncut*n], ncut, sigma0, &dH1[ncut*m], ncut, sigma1, dist);
//                 auto const ovl  = overlap_of_two_Hermite_Gauss_functions(&H0[ncut*n], ncut, sigma0, &H1[ncut*m], ncut, sigma1, dist);
//                 if (echo > 1) printf(" %.9f", ovl); // show overlap
//              if (echo > 1) printf("  %.9f %.9f %.9f", d2d0, d0d2, -d1d1); // show 3 values
//              if (echo > 1) printf("  %.1e %.1e %.1e", d2d0 + d1d1, d0d2 + d1d1, d2d0 - d0d2); // show deviations
                if (echo > 1) printf(" %.9f", -d1d1); // show 1 value
                auto const d2avg = .5*d2d0 + .5*d0d2;
//              if (echo > 1) printf("  %.9f %.9f", d2avg, -d1d1); // show 2 values
                maxdev3 = std::max(maxdev3, std::abs(d2avg + d1d1)); // one order better than dev1 and dev2
                maxdev2 = std::max(maxdev2, std::abs(d2d0 - d0d2));
                maxdev1 = std::max(maxdev1, std::abs(d2d0 + d1d1));
                maxdev1 = std::max(maxdev1, std::abs(d0d2 + d1d1));
            } // m
        } // n
        if (echo > 1) printf("\n");
    } // dist
    if (echo > 0) printf("# %s deviations %g, %g and %g\n", __func__, maxdev1, maxdev2, maxdev3);
    return (maxdev3 > 2e-14);
  } // test

  status_t test_density_tensor(int const echo=0) {
    // this structure can be used to describe the density generation
    // for the density we assume that it is sufficient to
    // represent the density in a SHO basis 
    // with sigma_\rho = sigma/sqrt(2) and nu_max_\rho = 2\nu_max
    int constexpr ncut = 8;
    double constexpr sqrt2 = sqrt(2.);
    double H[ncut*ncut], Hp[2*ncut*2*ncut];
    prepare_centered_Hermite_polynomials(H, ncut); // unit spread sigma=1, L2-normalized
    prepare_centered_Hermite_polynomials(Hp, 2*ncut, sqrt2); // spread sigma_p = sigma/sqrt(2), L2-normalized
    double HHp[3*ncut], HHpH[4*ncut];
    for(int p = 0; p < 2*ncut - 1; ++p) {
        if (echo > 1) printf("\n# p = %d\n", p);
        for(int n = 0; n < ncut; ++n) {
            multiply(HHp, 3*ncut, &H[ncut*n], ncut, &Hp[2*ncut*p], 2*ncut);
            for(int m = 0; m < ncut; ++m) {
                if (0 == (p + n + m) % 2) { // odd contributions are zero by symmetry
                    multiply(HHpH, 4*ncut, &H[ncut*m], ncut, HHp, 3*ncut);
                    auto const P_pnm = integrate(HHpH, 4*ncut, 0.5);
//                  if (echo > 1) printf(" %d%d%d %.9f\n", p,n,m, P_pnm); // show tensor values as list
                    if (echo > 1) printf(" %.9f", P_pnm); // show tensor values
                } // even?
            } // m
            if (echo > 1) printf("\n");
        } // n
    } // p
    return 0;
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
  
  template<typename real_t> 
  real_t pow2(real_t const base) { return base*base; }  
  
  status_t test_fcc(int const echo=5, float const a0=8) {
    typedef vector_math::vec<3,double> vec3;
    typedef vector_math::vec<3,int>    vec3i;
    int constexpr nmax = 6, ncut = nmax + 2;
    
    vec3 cv[3], bv[3]; // vectors of the cell and the Bravais matrix
    {
        double const recip = (2*C_PI)/a0;
        for(int dir = 0; dir < 3; ++dir) {
            if (0) { // simple-cubic
                cv[dir] = 0; cv[dir][dir] = a0;
                bv[dir] = 0; bv[dir][dir] = recip; // sc
            } else {
                int8_t const cell_bcc[3][3] = { {0,1,1},  {1,0,1},  {1,1,0}}; // bcc
                int8_t const cell_fcc[3][3] = {{-1,1,1}, {1,-1,1}, {1,1,-1}}; // fcc
                if (0) { // body-centered-cubic
                    cv[dir] = cell_bcc[dir];
                    bv[dir] = cell_fcc[dir];
                } else { // face-centered-cubic
                    cv[dir] = cell_fcc[dir];
                    bv[dir] = cell_bcc[dir];
                }
                cv[dir] *= (a0*.5);
                bv[dir] *= (recip);
            }
        } // dir
    }

    double shortest_bond2 = 9e99;
    for(int i3 = -1; i3 <= 1; ++i3) {
        for(int i2 = -1; i2 <= 1; ++i2) {
            for(int i1 = -1; i1 <= 1; ++i1) {
                vec3 const pos = cv[0]*i1 + cv[1]*i2 + cv[2]*i3; // assume one atom per unit cell
                double const d2 = norm(pos);
                if (d2 > 0) shortest_bond2 = std::min(shortest_bond2, d2);
            } // i1
        } // i2
    } // i3
    double const shortest_bond = std::sqrt(shortest_bond2);

    printf("# shortest bond is %g Bohr\n", shortest_bond);

    bool const overlap_eigvals = false;
    // choose the return radius as a fraction of shortest_bond length
    double const sigma = .75*shortest_bond/std::sqrt(2.*nmax + 3.), 
                 sigma0 = sigma, sigma1 = sigma;
    printf("# SHO up to numax=%d, spread sigma = %.9f Bohr\n", nmax, sigma);

    // return radius of the classical harmonic oscillator
    double const return_radius = sigma*std::sqrt(2.*nmax + 3.);
    printf("# classical return radius at %g Bohr\n", return_radius);
    
    double const dmax = 12*sigma; // 12 sigma is converged for fcc
    printf("# account for periodic images up to %.3f Bohr\n", dmax);

    double H0[ncut*ncut], H1[ncut*ncut], normalize=0; // do not normalize
    prepare_centered_Hermite_polynomials(H0, ncut, 1./sigma0, normalize);
    prepare_centered_Hermite_polynomials(H1, ncut, 1./sigma1, normalize);

    double dH0[ncut*ncut], dH1[ncut*ncut];
    for(int n = 0; n < ncut; ++n) {
        // show the Hermite polynomial coefficients for H0
        printf("# H[%x]: ", n);
        for(int m = 0; m <= n; ++m) {
            printf("%8.4f", H0[n*ncut + m]);
        }   printf("\n");
        
        // construct first derivatives
        derive_Hermite_Gauss_polynomials(&dH0[n*ncut], &H0[n*ncut], ncut, 1./sigma0);
        derive_Hermite_Gauss_polynomials(&dH1[n*ncut], &H1[n*ncut], ncut, 1./sigma1);
    } // n

    int const n3D = ((nmax + 1)*(nmax + 2)*(nmax + 3))/6;
    printf("# %d SHO functions up to numax=%d\n", n3D, nmax);
    {   printf("# list %d SHO functions: ", n3D);
        for(int n0 = 0; n0 <= nmax; ++n0) {
            for(int n1 = 0; n1 <= nmax - n0; ++n1) {
                for(int n2 = 0; n2 <= nmax - n0 - n1; ++n2) {
                    printf("%x%x%x ", n0,n1,n2);
                } // 2n2
            } // n1
        } // n0
        printf("\n");
    } // scope
    
    vec3i const imax = std::ceil(dmax/a0);
    int const max_npi = 8*imax[2]*imax[1]*imax[0];
    auto mat = new double[max_npi][2][n3D*n3D];
    auto vpi = new    int[max_npi][3]; // periodic image shift vectors
    int npi = 0;
    for(int i3 = -imax[2]; i3 <= imax[2]; ++i3) {
        for(int i2 = -imax[1]; i2 <= imax[1]; ++i2) {
            for(int i1 = -imax[0]; i1 <= imax[0]; ++i1) {
                vec3 const pos = cv[0]*i1 + cv[1]*i2 + cv[2]*i3;

                if (norm(pos) < dmax*dmax) {
                    if (echo > 9) printf("%f %f %f\n", pos[0],pos[1],pos[2]);
                    int in = 0;
                    for(int n0 = 0; n0 <= nmax; ++n0) {
                    for(int n1 = 0; n1 <= nmax - n0; ++n1) {
                    for(int n2 = 0; n2 <= nmax - n0 - n1; ++n2) {
                        int const nv[] = {n0, n1, n2};
                        int im = 0;
                        for(int m0 = 0; m0 <= nmax; ++m0) {
                        for(int m1 = 0; m1 <= nmax - m0; ++m1) {
                        for(int m2 = 0; m2 <= nmax - m0 - m1; ++m2) {
                            int const mv[] = {m0, m1, m2};
                            double ovl[] = {0, 0, 0};
                            double lap[] = {0, 0, 0};
                            // ToDo: overlap_of_two_Hermite_Gauss_functions 
                            //       is called many more times than necessary
                            //       and the max. length of non-zero polynomial coefficients 
                            //       can be shorter than ncut in many cases
                            for(int dir = 0; dir < 3; ++dir) {
                                ovl[dir] = overlap_of_two_Hermite_Gauss_functions(
                                               &H0[nv[dir]*ncut], ncut, sigma0,
                                               &H1[mv[dir]*ncut], ncut, sigma1, pos[dir]);
                                lap[dir] = overlap_of_two_Hermite_Gauss_functions(
                                              &dH0[nv[dir]*ncut], ncut, sigma0,
                                              &dH1[mv[dir]*ncut], ncut, sigma1, pos[dir]);
                            } // dir
                            double const o3D = ovl[0]*ovl[1]*ovl[2];
                            double const l3D = lap[0]*ovl[1]*ovl[2] + ovl[0]*lap[1]*ovl[2] + ovl[0]*ovl[1]*lap[2];
                            mat[npi][0][in*n3D + im] = o3D;
                            mat[npi][1][in*n3D + im] = l3D;
                            ++im;
                        }}} // m
                        ++in;
                    }}} // n
                    vpi[npi][0] = i1; vpi[npi][1] = i2; vpi[npi][2] = i3;
                    ++npi; // count periodic images
                } // pos inside sphere
            } // i1
        } // i2
    } // i3
    int const num_periodic_images = npi;


    double smallest_eigval = 9e99, largest_eigval = - 9e99;
    vec3 kv_smallest = -9;

    int const lwork = n3D*n3D;
    complex_t ovl_mat[n3D*n3D], lap_mat[n3D*n3D], work[lwork];
    double rwork[lwork], eigvals[n3D];
    auto const jobz = 'n', uplo = 'u', jobv = 'v';
    
    int diagonalization_failed = 0;
    int const nedges = 4;
    float const sampling_density = .03125/8;
    double const kpath[nedges][3] = {{.0,.0,.0}, {.5,.0,.0}, {.5,.5,.0}, {.5,.5,.5}};
    float path_progress = 0;
    for(int edge = 0; edge < nedges; ++edge) {
        int const e0 = edge % nedges, e1 = (edge + 1) % nedges;
        vec3 const v0 = kpath[e0], v1 = kpath[e1];
        vec3 const true_kdiff = bv[0]*(v1[0] - v0[0]) + bv[1]*(v1[1] - v0[1]) + bv[2]*(v1[2] - v0[2]);
        double const edge_length = std::sqrt(norm(true_kdiff));

        int const sampling = std::ceil(edge_length/sampling_density);
        double const frac = 1./sampling;
        if (echo > 1) printf("# k-point %.6f %.6f %.6f\n", v0[0],v0[1],v0[2]);
        for(int step = 0; step < sampling + (edge == (nedges - 1)); ++step) {
            float const path_progress_edge = path_progress + (step*frac)*edge_length;

            vec3 const kvec = v0 + (v1 - v0)*(step*frac);
            vec3 const true_kv = bv[0]*kvec[0] + bv[1]*kvec[1] + bv[2]*kvec[2];
            double const free_electron = norm(true_kv); // in Rydberg atomic units

            // clear matrixes
            for(int in = 0; in < n3D; ++in) {
                for(int im = 0; im < n3D; ++im) {
                    ovl_mat[in*n3D + im] = 0;
                    lap_mat[in*n3D + im] = 0;
                } // im
            } // in
            
            for(int ipi = 0; ipi < num_periodic_images; ++ipi) {
                vec3 ipos = vpi[ipi];
                complex_t const bloch_factor = std::polar(1.0, 2*C_PI * dot(kvec, ipos));
                if (echo > 9) printf("# periodic image%4d%4d%4d  Bloch-phase = %f + i %f\n", 
                    vpi[ipi][0], vpi[ipi][1], vpi[ipi][2], bloch_factor.real(), bloch_factor.imag());
                // add to matrixes
                for(int in = 0; in < n3D; ++in) {
                    for(int im = 0; im < n3D; ++im) {
                        ovl_mat[in*n3D + im] += bloch_factor*mat[ipi][0][in*n3D + im];
                        lap_mat[in*n3D + im] += bloch_factor*mat[ipi][1][in*n3D + im];
                    } // im
                } // in
            } // ipi

#if 0            
            // check if matrices are hermitian
            auto const threshold = 1e-5;
            for(auto m = ovl_mat; m == ovl_mat || m == lap_mat; m += (lap_mat - ovl_mat)) {
                for(int in = 0; in < n3D; ++in) {
                    for(int im = 0; im < in; ++im) {
                        assert(std::abs(m[in*n3D + im].real() - m[im*n3D + in].real()) < threshold);
                        assert(std::abs(m[in*n3D + im].imag() + m[im*n3D + in].imag()) < threshold);
                    } // im
                    assert(std::abs(m[in*n3D + in].imag()) < threshold);
                } // in
            } // m
#endif
            
            // LAPACK call (Fortran77 interface);
            int info = 0, itype = 1;
            if (overlap_eigvals) {
                // get the eigenvalues of the overlap operator only
                zheev_(&jobv, &uplo, &n3D, ovl_mat, &n3D, 
                       eigvals, work, &lwork, rwork, &info);
#if 0
                // DEBUG
                if (0 == info && eigvals[0] < .00315) {
                   printf("# lowest eigenvector "); 
                   for(int i3D = 0; i3D < n3D; ++i3D) {
                      auto const c = ovl_mat[0*n3D + i3D];
                      printf("(%.9f,%.9f) ", c.real(), c.imag());
                   }  printf("\n");
                } // DEBUG
#endif
            } else {
                // solve generalized eigenvalue problem lap_mat*X == diag*ovl_mat*X
                zhegv_(&itype, &jobz, &uplo, &n3D, lap_mat, &n3D, ovl_mat, &n3D, 
                       eigvals, work, &lwork, rwork, &info);
            }

            for(int i3D = 0; i3D < n3D; ++i3D) {
                if (eigvals[i3D] < smallest_eigval) { 
                    kv_smallest = kvec;
                } // store where the smallest eigenvalue was found
                smallest_eigval = std::min(smallest_eigval, eigvals[i3D]);
                largest_eigval  = std::max(largest_eigval,  eigvals[i3D]);
            } // i3D
            
            if (0 == info) {
                if(echo > 2) {
                    printf("%.6f ", path_progress_edge);
                    if (!overlap_eigvals) printf("%g ", free_electron);
                    for(int i3D = 0; i3D < n3D; ++i3D) {
                        printf("%g ", eigvals[i3D]);
                    } // i3D
                    printf("\n");
                } // echo
            } else {
                ++diagonalization_failed;
                if (echo > 2) printf("# %.6f diagonalization failed, info = %d\n", path_progress_edge, info);
            } // info

        } // step
        if (echo > 2) printf("# k-point %.6f %.6f %.6f\n", v1[0],v1[1],v1[2]);
        path_progress += edge_length;
    } // edge

    if (diagonalization_failed > 0) {
        if (echo > 0) printf("# Warning: %d diagonalizations failed!\n", diagonalization_failed);
    } else {
        if (echo > 1) printf("\n# smallest and largest eigenvalue%s are %g and %g\n", 
            overlap_eigvals?" of the overlap operator":"", smallest_eigval, largest_eigval);
        if (echo > 1) printf("# smallest eigenvalue at kvec  %.6f %.6f %.6f\n", kv_smallest[0],kv_smallest[1],kv_smallest[2]);
    }
    return diagonalization_failed;
  } // test_fcc

  status_t all_tests() {
    auto status = 0;
    status += test_Hermite_polynomials();
    status += test_Hermite_Gauss_overlap();
    status += test_kinetic_overlap();
    status += test_density_tensor();
    status += test_fcc();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace overlap
