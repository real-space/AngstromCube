#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::sin
#include <vector> // std::vector<T>

#include "bessel_transform.hxx"

#include "display_units.h" // Ang, _Ang
#include "radial_grid.hxx" // create_exponential_radial_grid
#include "inline_math.hxx" // pow2

// #define FULL_DEBUG
// #define DEBUG

namespace bessel_transform {
  // a radial bessel transformation of s-functions (ell=0)
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_Gaussian(int const echo=4) {
      if (echo > 3) printf("# %s: %s\n", __FILE__, __func__);
      auto const & rg = *radial_grid::create_exponential_radial_grid(1 << 9);
      int const nq = 80; double const dq = 0.125; auto const qcut = (nq - 1)*dq;
      std::vector<double> in(rg.n), bt(nq), out(rg.n); // get memory
      for(int ir = 0; ir < rg.n; ++ir) {
          in[ir] = std::exp(-.5*pow2(rg.r[ir])); // this function is its own Bessel-transform
      } // ir
      transform_s_function(bt.data(), in.data(), rg, nq, dq); // transform to Bessel-space
      if (echo > 6) {
          printf("## %s Bessel transformed function up to %g sqRyd\n", __func__, qcut);
          for(int iq = 0; iq < nq; ++iq) {
              printf("%g %g\n", iq*dq, bt[iq]); // show the Bessel-transformed function
          }   printf("\n\n");
      } // echo
      transform_s_function(out.data(), bt.data(), rg, nq, dq, true); // transform back to real-space again
      double dev[] = {0, 0, 0};
      if (echo > 5) printf("## %s real-space functions (in and out)\n", __func__);
      for(int ir = 0; ir < rg.n; ++ir) {
          dev[0] += rg.r2dr[ir];
          dev[1] += rg.r2dr[ir] * std::abs(out[ir] - in[ir]);
          dev[2] += rg.r2dr[ir] *     pow2(out[ir] - in[ir]);
          if (echo > 5) printf("%g %g %g\n", rg.r[ir], out[ir], in[ir]); // show the output and input vs r
      }   if (echo > 5) printf("\n\n");
      if (echo > 2) printf("# %s after filtering with cutoff %g sqRyd deviation is %.1e (abs) or %.1e (L2)\n",
                              __func__, qcut, dev[1]/dev[0], std::sqrt(dev[2]/dev[0]));
      return (dev[1]/dev[0] > 5e-15);
  } // test_Gaussian

  status_t test_r2grid(int const echo=9) {
      if (echo > 0) printf("# %s: %s\n", __FILE__, __func__);
      auto const & rg = *radial_grid::create_exponential_radial_grid(1 << 9);
      float const ar2 = (1 << 3); // ir2 = ar2*r^2
      int   const nr2 = std::ceil(ar2*pow2(rg.rmax));
      std::vector<double> in(rg.n), out(nr2);
      if (echo > 4) printf("\n## %s input (as function of r^2):\n", __func__);
      for(int ir = 0; ir < rg.n; ++ir) {
          double const r = rg.r[ir];
          in[ir] = std::exp(-.5*pow2(r)); // this function is its own Bessel-transform
          in[ir] *= std::cos(r*r); // modify it somehow
          if (echo > 4) printf("%g %g\n", r*r, in[ir]); // show the input function with r^2 abscissa
      } // ir
      if (echo > 4) printf("\n\n");
      auto const stat = transform_to_r2_grid(out.data(), ar2, nr2, in.data(), rg);
      if (echo > 4) {
          printf("\n## %s output (as function of r^2):\n", __func__);
          for(int ir2 = 0; ir2 < nr2; ++ir2) {
              printf("%g %g\n", ir2/ar2, out[ir2]); // show the output function with r^2 abscissa
          }   printf("\n\n");
      } // echo
      return stat;
  } // test_r2grid

  status_t all_tests(int const echo) {
    status_t stat{0};
    stat += test_Gaussian(echo);
    stat += test_r2grid(echo);
    return stat;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace bessel_transform
