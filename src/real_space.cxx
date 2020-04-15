#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::floor

#include "real_space.hxx"

#include "display_units.h" // eV, _eV, Ang, _Ang
#include "constants.hxx" // pi

// #define FULL_DEBUG
// #define DEBUG

namespace real_space {

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo=9) {
      int const dims[] = {10, 20, 30};
      auto gp = new grid_t<1>(dims); // construct
      gp->~grid_t(); // explicity call the destructor
      int const d1ms[] = {15, 20, 25};
      grid_t<4> g(d1ms); // construct vectorized, destruction automatically
      return 0;
  } // test_create_and_destroy

  status_t test_add_function(int const echo=9) {
      if (echo > 0) printf("\n# %s\n", __func__);
      int const dims[] = {32, 31, 30};
      grid_t<1> g(dims);
      g.set_grid_spacing(0.333);
      double const cnt[] = {g[0]*.42*g.h[0], 
                            g[1]*.51*g.h[1], 
                            g[2]*.60*g.h[2]}; // center is slightly shifted from exact grid point positions
      int const nr2 = 1 << 11;
      double r2c[nr2], rad_integral = 0;
      float const rcut = 4;
      float const inv_hr2 = nr2/(rcut*rcut);
      double const hr2 = 1./inv_hr2;
      if (echo > 4) printf("\n# values on the radial grid\n");
      for(int ir2 = 0; ir2 < nr2; ++ir2) { // sample r^2
          double const r2 = ir2*hr2, r = std::sqrt(r2);
          r2c[ir2] = std::exp(-r2); // function evaluation here
          if (echo > 4) printf("%g %g\n", r, r2c[ir2]); // plot function value vs r
          rad_integral += r2c[ir2] * r;
      } // ir2
      rad_integral *= 2*constants::pi/inv_hr2;
      if (echo > 2) printf("\n# add_function()\n\n");
      double added;
      auto values = new double[g.all()];
      set(values, g.all(), 0.0);
      add_function(values, g, &added, r2c, nr2, inv_hr2, cnt);
      if (echo > 6) printf("\n# non-zero values on the Cartesian grid (sum = %g)\n", added);
      double xyz_integral = 0;
      for(        int iz = 0; iz < g('z'); ++iz) {  double const vz = iz*g.h[2] - cnt[2];
          for(    int iy = 0; iy < g('y'); ++iy) {  double const vy = iy*g.h[1] - cnt[1];
              for(int ix = 0; ix < g('x'); ++ix) {  double const vx = ix*g.h[0] - cnt[0];
                  auto const ixyz = (iz*g('y') + iy)*g('x') + ix;
                  auto const val = values[ixyz];
                  if (0 != val) {
                      if (echo > 6) printf("%g %g\n", std::sqrt(vz*vz + vy*vy + vx*vx), val); // plot function value vs r
                      xyz_integral += val;
                  } // non-zero
              } // ix
          } // iy
      } // iz
      xyz_integral *= g.dV(); // volume element
      auto const diff = xyz_integral - rad_integral;
      if (echo > 1) printf("# grid integral = %g  radial integral = %g  difference = %.1e (%.2f %%)\n", 
                                  xyz_integral, rad_integral, diff, 100*diff/rad_integral);
      return std::abs(diff/rad_integral) > 4e-4;
  } // test_add_function

  status_t all_tests(int const echo) {
    status_t status(0);
    status += test_create_and_destroy(echo);
    status += test_add_function(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace real_space