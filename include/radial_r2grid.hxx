#pragma once

#include <cmath> // std::sqrt
#include <vector> // std::vector<T>

#include "status.hxx" // status_t

namespace radial_r2grid {

  /*
   *  The radial_r2grid is a radial grid descriptor
   *  with equidistant spacing in r^2, i.e.
   *        r_i = \sqrt(i/a)
   *  or
   *        ar^2 = i
   *  where a is restricted to float to reduce noise.
   *
   *  The r2grid is particularly useful when we want to
   *  bring spherical radial functions to a 3D grid
   *  since we can do the radial interpolation without
   *  the necessity to compute a sqrt.
   */

  template <typename real_t=double>
  std::vector<real_t> r_axis(int const nr2, float const ar2=1) {
      double const ar2inv = 1.0/ar2;
      std::vector<real_t> r(nr2);
      for (int ir2 = 0; ir2 < nr2; ++ir2) {
          r[ir2] = std::sqrt(ir2*ar2inv); // r^2-grid formula
          // r^2-grid inversion:      i0 = int(r2*ar2)
          //    --> interpolate between i0 and i0 + 1
      } // ir2
      return r;
  } // r_axis

// inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }

} // namespace radial_r2grid
