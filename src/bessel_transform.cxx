#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::sin

#include "bessel_transform.hxx"

#include "display_units.h" // Ang, _Ang
#include "radial_grid.hxx" // create_exponential_radial_grid
#include "inline_math.hxx" // pow2

// #define FULL_DEBUG
// #define DEBUG

namespace bessel_transform {
  // a radial bessel transformation of s-functions (ell=0)
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_Gaussian(int const echo=9) {
      printf("# %s: %s\n", __FILE__, __func__);
      auto const rg = radial_grid::create_exponential_radial_grid(1 << 9);
      auto const in = new double[rg->n];
      for(int ir = 0; ir < rg->n; ++ir) {
          in[ir] = std::exp(-.5*pow2(rg->r[ir])); // this function is its own Bessel-transform
          printf("%g %g\n", rg->r[ir], in[ir]); // show the input function
      }   printf("\n\n");
      int const nq = 80; double const dq = 0.125;
      auto const bt = new double[nq];
      transform_s_function(bt, in, *rg, nq, dq); // transform to Bessel-space
      for(int iq = 0; iq < nq; ++iq) {
          printf("%g %g\n", iq*dq, bt[iq]); // show the Bessel-transformed function
      }   printf("\n\n");
      auto const rs = new double[rg->n];
      transform_s_function(rs, bt, *rg, nq, dq, true); // transform back to real-space again
      for(int ir = 0; ir < rg->n; ++ir) {
          printf("%g %g\n", rg->r[ir], rs[ir]);
      }   printf("\n\n");
      return 0;
  } // test_Gaussian

  status_t test_r2grid(int const echo=9) {
      if (echo > 0) printf("# %s: %s\n", __FILE__, __func__);
      auto const rg = radial_grid::create_exponential_radial_grid(1 << 9);
      auto const in = new double[rg->n];
      for(int ir = 0; ir < rg->n; ++ir) {
          double const r = rg->r[ir];
          in[ir] = std::exp(-.5*pow2(r)); // this function is its own Bessel-transform
          in[ir] *= std::cos(r*r); // modify it somehow
      } // ir
      
      if (echo > 4) {
          printf("\n# %s input (as function of r^2):\n", __func__);
          for(int ir = 0; ir < rg->n; ++ir) {
              printf("%g %g\n", pow2(rg->r[ir]), in[ir]); // show the input function with r^2 abscissa
          }   printf("\n\n");
      } // echo

      float const ar2 = (1 << 3); // ir2 = ar2*r^2
      int   const nr2 = std::ceil(ar2*pow2(rg->rmax));
      auto  const out = new double[nr2];

      auto const stat = transform_to_r2_grid(out, ar2, nr2, in, *rg);
      
      if (echo > 4) {
          printf("\n# %s output (as function of r^2):\n", __func__);
          for(int ir2 = 0; ir2 < nr2; ++ir2) {
              printf("%g %g\n", ir2/ar2, out[ir2]); // show the output function with r^2 abscissa
          }   printf("\n\n");
      } // echo
      
      delete[] out;
      return stat;
  } // test_r2grid

  status_t all_tests(int const echo) {
    auto status = 0;
//  status += test_Gaussian(echo);
    status += test_r2grid(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace bessel_transform
