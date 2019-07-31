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

  template<typename real_t>
  status_t transform_to_r2_grid(real_t out[], float const ar2, int const nr2,
                                double const in[], radial_grid_t const &g, int const echo=9) {
      if (echo > 8) {
          printf("\n# %s input:\n", __func__);
          for(int ir = 0; ir < g.n; ++ir) {
              printf("%g %g\n", g.r[ir], in[ir]);
          }   printf("\n\n");
      } // echo
    
      int const nq = 64; double const dq = 0.125; // choose a q-grid
      auto const bt = new double[nq];
      transform_s_function(bt, in, g, nq, dq); // transform to Bessel-space
      
      // construct a r2-grid descriptor: ir2 = ar2*r^2
      double const ar2inv = 1./ar2;
      radial_grid_t r2g; r2g.n = nr2; r2g.r = new double[nr2];
      for(int ir2 = 0; ir2 < nr2; ++ir2) {
          r2g.r[ir2] = std::sqrt(ir2*ar2inv);
      } // ir2

      auto const stat = transform_s_function(out, bt, r2g, nq, dq, true); // transform back to real-space

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
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
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
      printf("# %s: %s\n", __FILE__, __func__);
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
  
  status_t all_tests() {
    auto status = 0;
//  status += test_Gaussian();
    status += test_r2grid();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace bessel_transform