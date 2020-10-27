#pragma once

#include "status.hxx" // status_t

#include "real_space.hxx" // ::grid_t
#include "display_units.h" // Ang, _Ang

namespace poisson_solver {

  status_t solve(
        double Ves[] // result electrostatic potential
      , double const rho[] // charge density, typically augmented density, should be charge neutral
      , real_space::grid_t const & g // grid descriptor
      , int const echo=0 // log-level
      , double const Bessel_center[]=nullptr
  ); // declaration only

  inline
  char solver_method(char const *method) { return *method; }
  
  template <typename real_t>
  void print_direct_projection(
        real_t const array[]
      , real_space::grid_t const &g
      , double const factor=1
      , double const *center=nullptr
  ) {
      // write all values of a grid array to stdout 
      // as function of their distance to a given center

      double cnt[3];
      if (nullptr != center) { 
          set(cnt, 3, center); // copy
      } else {
          for(int d = 0; d < 3; ++d) {
              cnt[d] = 0.5*(g[d] - 1)*g.h[d];
          } // d
      } // center given
      printf("# projection center (relative to grid point (0,0,0) is %g %g %g in units of grid spacings\n",
                cnt[0]*g.inv_h[0], cnt[1]*g.inv_h[1], cnt[2]*g.inv_h[2]);
      for(int iz = 0; iz < g[2]; ++iz) {
          double const z = iz*g.h[2] - cnt[2], z2 = z*z;
          for(int iy = 0; iy < g[1]; ++iy) {
              double const y = iy*g.h[1] - cnt[1], y2 = y*y; 
              for(int ix = 0; ix < g[0]; ++ix) {
                  double const x = ix*g.h[0] - cnt[0], x2 = x*x;
                  double const r = std::sqrt(x2 + y2 + z2);
                  int const izyx = (iz*g[1] + iy)*g[0] + ix;
                  printf("%g %g\n", r*Ang, array[izyx]*factor);
              } // ix
          } // iy
      } // iz
      printf("# radii in %s\n\n", _Ang);
  } // print_direct_projection
  
  
  status_t all_tests(int const echo=0);

} // namespace poisson_solver
