#pragma once

#include "inline_math.hxx" // product
#include "radial_grid.h" // radial_grid_t
#include "constants.hxx" // pi

typedef int status_t;

namespace bessel_transform {
  
  inline double j0(double const x) { return (x*x < 1e-16) ? 1.0 - x*x/6. : std::sin(x)/x; }
  
  template<typename real_t>
  status_t transform_s_function(real_t out[], // result
                  real_t const in[], radial_grid_t const &g, int const n, 
                  double const dq=.125, bool const back=false, int const echo=9) {

      auto q_lin = new double[n];
      for(int iq = 0; iq < n; ++iq) {
          q_lin[iq] = iq*dq;
      } // iq

      int n_out, n_in;
      double *x_out, *x_in, *dx_in;
      if (back) {
          // from reciprocal space q to real-space r
          n_in  = n;
          x_in  = q_lin;
          dx_in = new double[n];
          product(dx_in, n, q_lin, q_lin, dq);
          n_out = g.n;
          x_out = g.r;
      } else {
          // from real-space r to reciprocal space q
          n_in  = g.n;
          x_in  = g.r;
          dx_in = g.r2dr;
          n_out = n;
          x_out = q_lin;
      } // back-transform ?

      double const sqrt2pi = std::sqrt(2./constants::pi); // this makes the transform symmetric
      for(int io = 0; io < n_out; ++io) {
          double tmp = 0;
          for(int ii = 0; ii < n_in; ++ii) {
              double const qr = x_in[ii]*x_out[io];
              tmp += in[ii] * j0(qr) * dx_in[ii];
          } // ii
          out[io] = tmp*sqrt2pi;
      } // io
      
      if (back) delete[] dx_in;
      delete[] q_lin;
      
      return 0;
  } // transform_s_function
  
  status_t all_tests();

} // namespace bessel_transform
