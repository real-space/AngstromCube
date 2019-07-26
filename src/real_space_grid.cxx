#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::floor

#include "real_space_grid.hxx"

#include "display_units.h" // eV, _eV, Ang, _Ang
#include "constants.hxx" // pi

// #define FULL_DEBUG
// #define DEBUG

namespace real_space_grid {

  template<typename real_t, int D0> // inner dimension
  status_t add_function(grid_t<real_t,D0> &g, // grid values are modified
                        double const r2coeff[][D0], int const ncoeff, float const hcoeff,
                        double const center[3]=nullptr, float const rcut=-1, double const scale=1) {
  // Add a spherically symmetric regular function to the grid.
  // The function is tabulated as r2coeff[0 <= hcoeff*r^2 < ncoeff] 
      assert(D0 == g.dim('w'));
      double c[3] = {0,0,0}; if (center) std::copy(center, 3+center, c);
      bool const cutoff = (rcut >= 0); // negative values mean that we do not need a cutoff
      double const r2cut = cutoff ? rcut*rcut : 1.5e308; // infinity
      int imn[3], imx[3];
      size_t nwindow = 1;
      for(int i3 = 0; i3 < 3; ++i3) {
          int const M = g.dim(i3) - 1; // highest index
          if (cutoff) {
              imn[i3] = std::max(0, (int)std::floor((c[i3] - rcut)*g.inv_h[i3]));
              imx[i3] = std::min(M, (int)std::ceil ((c[i3] + rcut)*g.inv_h[i3]));
          } else {
              imn[i3] = 0;
              imx[i3] = M;
          } // cutoff
#ifdef  DEBUG
          printf("# %s window %c = %d elements from %d to %d\n", __func__, 'x'+i3, imx[i3] + 1 - imn[i3], imn[i3], imx[i3]);
#endif
          nwindow *= std::max(0, imx[i3] + 1 - imn[i3]);
      } // i3
      assert(hcoeff > 0);
      size_t modified = 0, out_of_range = 0;
      for(            int iz = imn[2]; iz <= imx[2]; ++iz) {  double const vz = iz*g.h[2] - c[2], vz2 = vz*vz;
          for(        int iy = imn[1]; iy <= imx[1]; ++iy) {  double const vy = iy*g.h[1] - c[1], vy2 = vy*vy;
              if (vz2 + vy2 < r2cut) {
                  for(int ix = imn[0]; ix <= imx[0]; ++ix) {  double const vx = ix*g.h[0] - c[0], vx2 = vx*vx;
                      double const r2 = vz2 + vy2 + vx2;
                      if (r2 < r2cut) {
                          int const ixyz = (iz*g.dim('y') + iy)*g.dim('x') + ix;
                          int const ir2 = (int)(hcoeff*r2);
                          if (ir2 < ncoeff) {
                              double const w8 = hcoeff*r2 - ir2; // linear interpolation weight
                              int const ir2p1 = ir2 + 1;
                              for(int i0 = 0; i0 < D0; ++i0) { // vectorize
                                  g.values[ixyz*D0 + i0] += scale*(r2coeff[ir2][i0]*(1 - w8) 
                                           + ((ir2p1 < ncoeff) ? r2coeff[ir2p1][i0] : 0)*w8);
                              } // i0
                              ++modified;
                          } else ++out_of_range;
                      } // inside rcut
                  } // ix
              } // rcut for (y,z)
          } // iy
      } // iz
#ifdef  DEBUG
      printf("# %s modified %.3f k inside a window of %.3f k on a grid of %.3f k grid values.\n", 
              __func__, modified*1e-3, nwindow*1e-3, g.dim('x')*g.dim('y')*g.dim('z')*1e-3); // show stats
#endif
      if (out_of_range)
          printf("# Warning! %s modified %ld points and found %ld entries out of range of the radial function!\n", 
              __func__, modified, out_of_range); // show stats
      return 0; // success
  } // add_function

#if 0
  template<typename real_t, int D0> // inner dimension
  status_t Fourier_transform(real_t *out, grid_t<real_t,D0> &in) {
      assert(-1 == D0); // no FFT implemented, break here

      for(int w = 0; w < in.dim('w'); ++w) { // outer dimension
          for(int z = 0; z < in.dim('z'); ++z) {
              for(int y = 0; y < in.dim('y'); ++y) {
                  for(int x = 0; x < in.dim('x'); ++x) {
                      int const wzyx = ((w*in.dim('z') + z)*in.dim('y') + y)*in.dim('x') + x;
                      for(int d0 = 0; d0 < D0; ++d0) { // vectorize
                          out[wzyx*D0 + d0] = in.values[wzyx*D0 + d0]; // plain copy
                      } // d0
                  } // x
              } // y
          } // z
      } // w

      return 0;
  } // Fourier transform
#endif

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int echo=9) {
      int const dims[] = {10, 20, 30};
      auto gp = new grid_t<float,1>(dims); // allocate and construct
      gp->~grid_t(); // explicity call the destructor
      int const d1ms[] = {15, 20, 25};
      grid_t<double,4> g(d1ms); // construct, destruction automatic
      return 0;
  } // test_create_and_destroy

  status_t test_add_function(int echo=9) {
      if (echo > 0) printf("\n# %s\n", __func__);
      int const dims[] = {32, 31, 30};
      grid_t<double,1> g(dims);
      std::fill(g.values, g.all() + g.values, 0.0);
      g.set_grid_spacing(0.333);
      double const cnt[] = {g.dim('x')*.42*g.h[0], 
                            g.dim('y')*.51*g.h[1], 
                            g.dim('z')*.60*g.h[2]}; // center is slightly shifted from exact grid point positions
      int const nr2 = 1 << 11;
      double r2c[nr2][1], rad_integral = 0;
      float const rcut = 4;
      float const inv_hr2 = nr2/(rcut*rcut);
      double const hr2 = 1./inv_hr2;
      if (echo > 4) printf("\n# values on the radial grid\n");
      for(int ir2 = 0; ir2 < nr2; ++ir2) { // sample r^2
          double const r2 = ir2*hr2, r = std::sqrt(r2);
          r2c[ir2][0] = std::exp(-r2); // function evaluation here
          if (echo > 4) printf("%g %g\n", r, r2c[ir2][0]); // plot function value vs r
          rad_integral += r2c[ir2][0] * r;
      } // ir2
      rad_integral *= 2*constants::pi/inv_hr2;
      if (echo > 2) printf("\n# add_function()\n\n");
      add_function(g, r2c, nr2, inv_hr2, cnt, rcut);
      if (echo > 6) printf("\n# non-zero values on the Cartesian grid\n");
      double xyz_integral = 0;
      int ixyz = 0;
      for(        int iz = 0; iz < g.dim('z'); ++iz) {  double const vz = iz*g.h[2] - cnt[2];
          for(    int iy = 0; iy < g.dim('y'); ++iy) {  double const vy = iy*g.h[1] - cnt[1];
              for(int ix = 0; ix < g.dim('x'); ++ix) {  double const vx = ix*g.h[0] - cnt[0];
                  auto const val = g.values[ixyz];
                  if (0 != val) {
                      if (echo > 6) printf("%g %g\n", std::sqrt(vz*vz + vy*vy + vx*vx), val); // plot function value vs r
                      xyz_integral += val;
                  } // non-zero
                  ++ixyz;
              } // ix
          } // iy
      } // iz
      xyz_integral *= g.dV(); // volume element
      auto const diff = xyz_integral - rad_integral;
      if (echo > 1) printf("# grid integral = %g  radial integral = %g  difference = %.1e (%.2f %%)\n", 
                                  xyz_integral, rad_integral, diff, 100*diff/rad_integral);
      return std::abs(diff/rad_integral) > 4e-4;
  } // test_add_function

  status_t all_tests() {
    auto status = 0;
    status += test_create_and_destroy();
    status += test_add_function(2);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace real_space_grid
