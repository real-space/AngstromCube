// This file is part of AngstromCube under MIT License

#include <cmath> // std::cos, std::exp

#include "bessel_transform.hxx"

#include "inline_math.hxx" // pow2
#include "radial_grid.hxx" // ::create_radial_grid, ::destroy_radial_grid
#include "constants.hxx" // ::pi
#include "status.hxx" // STATUS_TEST_NOT_INCLUDED

namespace bessel_transform {

#ifdef    NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

  status_t test_Gaussian(int const echo=4) {
      if (echo > 3) std::printf("# %s: %s\n", __FILE__, __func__);
      auto & g = *radial_grid::create_radial_grid(1 << 9);
      int const nq = 80; double const dq = 0.125; auto const qcut = (nq - 1)*dq;
      std::vector<double> in(g.n), bt(nq), out(g.n); // get memory
      for (int ir = 0; ir < g.n; ++ir) {
          in[ir] = std::exp(-.5*pow2(g.r[ir])); // this function is its own Bessel-transform
      } // ir
      transform_s_function(bt.data(), in.data(), g, nq, dq); // transform to Bessel-space
      if (echo > 6) {
          std::printf("## %s Bessel transformed function up to %g sqRyd\n", __func__, qcut);
          for (int iq = 0; iq < nq; ++iq) {
              std::printf("%g %g\n", iq*dq, bt[iq]); // show the Bessel-transformed function
          } // iq
          std::printf("\n\n");
      } // echo
      transform_s_function(out.data(), bt.data(), g, nq, dq, true); // transform back to real-space again
      double dev[] = {0, 0, 0};
      if (echo > 5) std::printf("## %s real-space functions (in and out)\n", __func__);
      for (int ir = 0; ir < g.n; ++ir) {
          dev[0] += g.r2dr[ir];
          dev[1] += g.r2dr[ir] * std::abs(out[ir] - in[ir]);
          dev[2] += g.r2dr[ir] *     pow2(out[ir] - in[ir]);
          if (echo > 5) std::printf("%g %g %g\n", g.r[ir], out[ir], in[ir]); // show the output and input vs r
      } // ir
      if (echo > 5) std::printf("\n\n");
      if (echo > 2) std::printf("# %s after filtering with cutoff %g sqRyd deviation is %.1e (abs) or %.1e (L2)\n",
                              __func__, qcut, dev[1]/dev[0], std::sqrt(dev[2]/dev[0]));
      radial_grid::destroy_radial_grid(&g);
      return (dev[1]/dev[0] > 5e-15);
  } // test_Gaussian

  status_t test_r2grid(int const echo=9) {
      if (echo > 0) std::printf("# %s: %s\n", __FILE__, __func__);
      auto & g = *radial_grid::create_radial_grid(1 << 9);
      float const ar2 = (1 << 3); // ir2 = ar2*r^2
      int   const nr2 = std::ceil(ar2*pow2(g.rmax));
      std::vector<double> in(g.n), out(nr2);
      if (echo > 4) std::printf("\n## %s input (as function of r^2):\n", __func__);
      for (int ir = 0; ir < g.n; ++ir) {
          double const r = g.r[ir];
          in[ir] = std::exp(-.5*pow2(r)); // this function is its own Bessel-transform
          in[ir] *= std::cos(r*r); // modify it somehow
          if (echo > 4) std::printf("%g %g\n", r*r, in[ir]); // show the input function with r^2 abscissa
      } // ir
      if (echo > 4) std::printf("\n\n");

      auto const stat = transform_to_r2grid(out.data(), ar2, nr2, in.data(), g);

      if (echo > 4) {
          std::printf("\n## %s output (as function of r^2):\n", __func__);
          double const ar2inv = 1./ar2;
          for (int ir2 = 0; ir2 < nr2; ++ir2) {
              std::printf("%g %g\n", ir2*ar2inv, out[ir2]); // show the output function with r^2 abscissa
          } // ir2
          std::printf("\n\n");
      } // echo
      radial_grid::destroy_radial_grid(&g);
      return stat; // this test needs a human to check that inout and output are similar
  } // test_r2grid

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_Gaussian(echo);
      stat += test_r2grid(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace bessel_transform
