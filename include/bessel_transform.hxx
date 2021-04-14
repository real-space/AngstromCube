#pragma once

#include <cstdio> // std::printf
#include <cmath> // std::sin
#include <vector> // std::vector<T>

#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // ::create_exponential_radial_grid
#include "status.hxx" // status_t

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
      for (int iq = 0; iq < nq; ++iq) {
          double const q = iq*dq;
          q_lin[iq] = q;
          dq_lin[iq] = q*q*dq;
      } // iq

      int n_out, n_in;
      double *x_out;
      double const *x_in, *dx_in;
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
      for (int io = 0; io < n_out; ++io) {
          double tmp{0};
          for (int ii = 0; ii < n_in; ++ii) {
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
          for (int ir = 0; ir < g.n; ++ir) {
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
      for (int ir2 = 0; ir2 < nr2; ++ir2) {
          r[ir2] = std::sqrt(ir2*ar2inv);
      } // ir2
      radial_grid_t r2g; r2g.n = nr2; r2g.r = r.data();

      stat += transform_s_function(out, bt.data(), r2g, nq, dq, true); // transform back to real-space

      if (echo > 8) {
          std::printf("\n# %s output:\n", __func__);
          for (int ir2 = 0; ir2 < r2g.n; ++ir2) {
              std::printf("%g %g\n", r2g.r[ir2], out[ir2]);
          } // ir2
          std::printf("\n\n");
      } // echo

      return stat;
  } // transform_to_r2_grid

  status_t all_tests(int const echo=0); // declaration only

} // namespace bessel_transform
