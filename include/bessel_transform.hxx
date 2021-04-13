#pragma once

#include <cstdio> // std::printf
#include <cmath> // std::sin
#include <vector> // std::vector<T>

#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // ::create_exponential_radial_grid
#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

#ifndef NO_UNIT_TESTS
    #include "inline_math.hxx" // pow2
//  #include "constants.hxx" // ::pi
//  #include <cmath> // std::cos, std::exp
#endif

namespace bessel_transform {
  
  inline double Bessel_j0(double const x) { return (x*x < 1e-16) ? 1.0 - x*x/6. : std::sin(x)/x; }

  template <typename real_t>
  inline status_t transform_s_function(
        real_t out[] // result
      , real_t const in[]
      , radial_grid_t const & g
      , int const nq
      , double const dq=.125
      , bool const back=false
      , int const echo=3 // log-level
  ) {
    
      if (echo > 8) std::printf("# %s(out=%p, in=%p, g=%p, nq=%d, dq=%.3f, back=%d, echo=%d);\n",
                           __func__, (void*)out, (void*)in, (void*)&g, nq, dq, back, echo);

      std::vector<double> q_lin(nq), dq_lin(nq);
      for(int iq = 0; iq < nq; ++iq) {
          double const q = iq*dq;
          q_lin[iq] = q;
          dq_lin[iq] = q*q*dq;
      } // iq

      int n_out, n_in;
      double *x_out, *x_in, *dx_in;
      if (back) {
          // from reciprocal space q to real-space r
          n_in  = nq;
          x_in  = q_lin.data();
          dx_in = dq_lin.data();
          n_out = g.n;
          x_out = g.r;
      } else {
          // from real-space r to reciprocal space q
          n_in  = g.n;
          x_in  = g.r;
          dx_in = g.r2dr;
          n_out = nq;
          x_out = q_lin.data();
      } // back-transform ?
      if (echo > 8) std::printf("# %s    n_in=%d x_in=%p dx_in=%p n_out=%d x_out=%p\n",
                      __func__, n_in, (void*)x_in, (void*)dx_in, n_out, (void*)x_out);

      double const sqrt2overpi = .7978845608028654; // this makes the transform symmetric
      //  assert(std::sqrt(2./constants::pi) == sqrt2overpi);
      for(int io = 0; io < n_out; ++io) {
          double tmp{0};
          for(int ii = 0; ii < n_in; ++ii) {
              double const qr = x_in[ii]*x_out[io];
              tmp += in[ii] * Bessel_j0(qr) * dx_in[ii];
          } // ii
          out[io] = tmp*sqrt2overpi;
      } // io

      return 0;
  } // transform_s_function

  template <typename real_t>
  inline status_t transform_to_r2_grid(
        real_t out[]
      , float const ar2
      , int const nr2
      , double const in[]
      , radial_grid_t const & g
      , int const echo=2
  ) {
      if (echo > 8) {
          std::printf("\n# %s input:\n", __func__);
          for(int ir = 0; ir < g.n; ++ir) {
              std::printf("%g %g\n", g.r[ir], in[ir]);
          } // ir
          std::printf("\n\n");
      } // echo

      int const nq = 256; double const dq = 0.125; // choose a q-grid: 255*0.125 = 32 = pi/0.1
//    int const nq = 128; double const dq = 0.125; // choose a q-grid: 127*0.125 = 16 = pi/0.2
//    int const nq = 64; double const dq = 0.125; // choose a q-grid: 64*0.125 = 8 = pi/0.39
      std::vector<double> bt(nq);
      auto stat = transform_s_function(bt.data(), in, g, nq, dq); // transform to Bessel-space

      // construct a r2-grid descriptor: ir2 = ar2*r^2
      double const ar2inv = 1./ar2;
      std::vector<double> r(nr2);
      for(int ir2 = 0; ir2 < nr2; ++ir2) {
          r[ir2] = std::sqrt(ir2*ar2inv);
      } // ir2
      radial_grid_t r2g; r2g.n = nr2; r2g.r = r.data();

      stat += transform_s_function(out, bt.data(), r2g, nq, dq, true); // transform back to real-space

      if (echo > 8) {
          std::printf("\n# %s output:\n", __func__);
          for(int ir2 = 0; ir2 < r2g.n; ++ir2) {
              std::printf("%g %g\n", r2g.r[ir2], out[ir2]);
          } // ir2
          std::printf("\n\n");
      } // echo

      return stat;
  } // transform_to_r2_grid

#ifdef NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_Gaussian(int const echo=4) {
      if (echo > 3) std::printf("# %s: %s\n", __FILE__, __func__);
      auto const & g = *radial_grid::create_exponential_radial_grid(1 << 9);
      int const nq = 80; double const dq = 0.125; auto const qcut = (nq - 1)*dq;
      std::vector<double> in(g.n), bt(nq), out(g.n); // get memory
      for(int ir = 0; ir < g.n; ++ir) {
          in[ir] = std::exp(-.5*pow2(g.r[ir])); // this function is its own Bessel-transform
      } // ir
      transform_s_function(bt.data(), in.data(), g, nq, dq); // transform to Bessel-space
      if (echo > 6) {
          std::printf("## %s Bessel transformed function up to %g sqRyd\n", __func__, qcut);
          for(int iq = 0; iq < nq; ++iq) {
              std::printf("%g %g\n", iq*dq, bt[iq]); // show the Bessel-transformed function
          } // iq
          std::printf("\n\n");
      } // echo
      transform_s_function(out.data(), bt.data(), g, nq, dq, true); // transform back to real-space again
      double dev[] = {0, 0, 0};
      if (echo > 5) std::printf("## %s real-space functions (in and out)\n", __func__);
      for(int ir = 0; ir < g.n; ++ir) {
          dev[0] += g.r2dr[ir];
          dev[1] += g.r2dr[ir] * std::abs(out[ir] - in[ir]);
          dev[2] += g.r2dr[ir] *     pow2(out[ir] - in[ir]);
          if (echo > 5) std::printf("%g %g %g\n", g.r[ir], out[ir], in[ir]); // show the output and input vs r
      } // ir
      if (echo > 5) std::printf("\n\n");
      if (echo > 2) std::printf("# %s after filtering with cutoff %g sqRyd deviation is %.1e (abs) or %.1e (L2)\n",
                              __func__, qcut, dev[1]/dev[0], std::sqrt(dev[2]/dev[0]));
      return (dev[1]/dev[0] > 5e-15);
  } // test_Gaussian

  inline status_t test_r2grid(int const echo=9) {
      if (echo > 0) std::printf("# %s: %s\n", __FILE__, __func__);
      auto const & g = *radial_grid::create_exponential_radial_grid(1 << 9);
      float const ar2 = (1 << 3); // ir2 = ar2*r^2
      int   const nr2 = std::ceil(ar2*pow2(g.rmax));
      std::vector<double> in(g.n), out(nr2);
      if (echo > 4) std::printf("\n## %s input (as function of r^2):\n", __func__);
      for(int ir = 0; ir < g.n; ++ir) {
          double const r = g.r[ir];
          in[ir] = std::exp(-.5*pow2(r)); // this function is its own Bessel-transform
          in[ir] *= std::cos(r*r); // modify it somehow
          if (echo > 4) std::printf("%g %g\n", r*r, in[ir]); // show the input function with r^2 abscissa
      } // ir
      if (echo > 4) std::printf("\n\n");

      auto const stat = transform_to_r2_grid(out.data(), ar2, nr2, in.data(), g);

      if (echo > 4) {
          std::printf("\n## %s output (as function of r^2):\n", __func__);
          double const ar2inv = 1./ar2;
          for(int ir2 = 0; ir2 < nr2; ++ir2) {
              std::printf("%g %g\n", ir2*ar2inv, out[ir2]); // show the output function with r^2 abscissa
          } // ir2
          std::printf("\n\n");
      } // echo
      return stat; // this test needs a human to check that inout and output are similar
  } // test_r2grid

  inline status_t all_tests(int const echo=0) {
      status_t stat{0};
      stat += test_Gaussian(echo);
      stat += test_r2grid(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace bessel_transform
