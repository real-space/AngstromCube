#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::sin

#include "bessel_transform.hxx"

#include "display_units.h" // Ang, _Ang
#include "inline_math.hxx" // set
#include "constants.hxx" // pi
#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // create_exponential_radial_grid

// #define FULL_DEBUG
// #define DEBUG

namespace bessel_transform {

  inline double j0(double const x) { return (x*x < 1e-16) ? 1.0 - x*x/6. : std::sin(x)/x; }
  
  template<typename real_t>
  status_t transform_s_function(real_t out[], // result
                  real_t const in[], radial_grid_t const &g, int const n, 
                  float const dq=.125f, bool const back=false, int const echo=9) {
      double const sqrt2pi = std::sqrt(2./constants::pi);
      
//       int const n_out = g.n;
//       auto x_out = g.r;
// //    auto dx_out = g.r2dr;

      double const dx = dq;
      auto x_lin = new double[n];
      auto dx_lin = new double[n];
      for(int i = 0; i < n; ++i) {
          double const x = i*dx;
          x_lin[i] = x;
          dx_lin[i] = x*x*dx;
      } // i
      
      int n_out, n_in;
      double *x_out, *x_in, *dx_in;
      if (back) {
          // from reciprocal space q to real-space r
          n_out = g.n;
          x_out = g.r;
          n_in  = n;
          x_in  = x_lin;
          dx_in = dx_lin;
      } else {
          // from real-space r to reciprocal space q
          n_out = n;
          x_out = x_lin;
          n_in  = g.n;
          x_in  = g.r;
          dx_in = g.r2dr;
      } // back-transform ?
      
      for(int i_out = 0; i_out < n_out; ++i_out) {
          auto const xo = x_out[i_out];
          double tmp = 0;
          for(int ii = 0; ii < n_in; ++ii) {
              double const qr = x_in[ii]*xo;
              tmp += in[ii] * j0(qr) * dx_in[ii];
//               printf("%g %g\n", qr, j0(qr));
          } // ii
          out[i_out] = tmp*sqrt2pi;
      } // i_out
      
      return 0;
  } // transform_s_function
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo=9) {
      printf("# %s: %s\n", __FILE__, __func__);
      auto const rg = radial_grid::create_exponential_radial_grid(1 << 9);
      auto const in = new double[rg->n];
      for(int ir = 0; ir < rg->n; ++ir) {
          in[ir] = std::exp(-.5*pow2(rg->r[ir])); // this function is its own Bessel-transform
          printf("%g %g\n", rg->r[ir], in[ir]); // show the input function
      }   printf("\n\n");
      int const nq = 80; float const dq = 0.125f;
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
  } // test_create_and_destroy

  status_t all_tests() {
    auto status = 0;
    status += test_create_and_destroy();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace bessel_transform
