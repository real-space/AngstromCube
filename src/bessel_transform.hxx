#pragma once

#include <cmath> // std::sin
#include <vector> // std::vector<T>

#include "inline_math.hxx" // product
#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // ::create_exponential_radial_grid
#include "constants.hxx" // ::pi

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

namespace bessel_transform {
  
  inline double Bessel_j0(double const x) { return (x*x < 1e-16) ? 1.0 - x*x/6. : std::sin(x)/x; }

  template<typename real_t>
  inline status_t transform_s_function(real_t out[], // result
                  real_t const in[], radial_grid_t const &g, int const nq, 
                  double const dq=.125, bool const back=false, int const echo=3) {
      if (echo > 8) printf("# %s(out=%p, in=%p, g=%p, nq=%d, dq=%.3f, back=%d, echo=%d);\n",
                           __func__, (void*)out, (void*)in, (void*)&g, nq, dq, back, echo);
      
      auto const q_lin = new double[nq];
      for(int iq = 0; iq < nq; ++iq) {
          q_lin[iq] = iq*dq;
      } // iq

      int n_out, n_in;
      double *x_out, *x_in, *dx_in;
      if (back) {
          // from reciprocal space q to real-space r
          n_in  = nq;
          x_in  = q_lin;
          dx_in = new double[nq];
          product(dx_in, nq, q_lin, q_lin, dq);
          n_out = g.n;
          x_out = g.r;
      } else {
          // from real-space r to reciprocal space q
          n_in  = g.n;
          x_in  = g.r;
          dx_in = g.r2dr;
          n_out = nq;
          x_out = q_lin;
      } // back-transform ?
      if (echo > 8) printf("# %s    n_in=%d x_in=%p dx_in=%p n_out=%d x_out=%p\n",
                   __func__, n_in, (void*)x_in, (void*)dx_in, n_out, (void*)x_out);

      double const sqrt2pi = std::sqrt(2./constants::pi); // this makes the transform symmetric
      for(int io = 0; io < n_out; ++io) {
          double tmp{0};
          for(int ii = 0; ii < n_in; ++ii) {
              double const qr = x_in[ii]*x_out[io];
              tmp += in[ii] * Bessel_j0(qr) * dx_in[ii];
          } // ii
          out[io] = tmp*sqrt2pi;
      } // io
      
      if (back) delete[] dx_in;
      delete[] q_lin;
      
      return 0;
  } // transform_s_function

  template<typename real_t>
  inline status_t transform_to_r2_grid(real_t out[], float const ar2, int const nr2,
                                double const in[], radial_grid_t const &g, int const echo=2) {
      if (echo > 8) {
          printf("\n# %s input:\n", __func__);
          for(int ir = 0; ir < g.n; ++ir) {
              printf("%g %g\n", g.r[ir], in[ir]);
          }   printf("\n\n");
      } // echo
    
      int const nq = 256; double const dq = 0.125; // choose a q-grid: 255*0.125 = 32 = pi/0.1
//    int const nq = 128; double const dq = 0.125; // choose a q-grid: 127*0.125 = 16 = pi/0.2
//    int const nq = 64; double const dq = 0.125; // choose a q-grid: 64*0.125 = 8 = pi/0.39
      auto const bt = new double[nq];
      auto stat = transform_s_function(bt, in, g, nq, dq); // transform to Bessel-space

      // construct a r2-grid descriptor: ir2 = ar2*r^2
      double const ar2inv = 1./ar2;
      radial_grid_t r2g; r2g.n = nr2; r2g.r = new double[nr2];
      for(int ir2 = 0; ir2 < nr2; ++ir2) {
          r2g.r[ir2] = std::sqrt(ir2*ar2inv);
      } // ir2

      stat += transform_s_function(out, bt, r2g, nq, dq, true); // transform back to real-space

      if (echo > 8) {
          printf("\n# %s output:\n", __func__);
          for(int ir2 = 0; ir2 < r2g.n; ++ir2) {
              printf("%g %g\n", r2g.r[ir2], out[ir2]);
          }   printf("\n\n");
      } // echo

      delete[] r2g.r;
      delete[] bt;
      return stat;
  } // transform_to_r2_grid

#ifdef NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_Gaussian(int const echo=4) {
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

  inline status_t test_r2grid(int const echo=9) {
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

  inline status_t all_tests(int const echo=0) {
      status_t stat{0};
      stat += test_Gaussian(echo);
      stat += test_r2grid(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace bessel_transform
