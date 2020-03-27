#pragma once

#include <cstdint> // uint32_t
#include <cmath> // std::sqrt
#include <vector> // std::vector<T>

// #include "status.hxx" // status_t

namespace radial_r2grid {

    template <typename real_t=float>
    std::vector<real_t> r_axis(int const nr2, float const ar2=1) {
        float const ar2inv = 1.0/ar2;
        std::vector<real_t> r(nr2);
        for(int ir2 = 0; ir2 < nr2; ++ir2) {
            r[ir2] = std::sqrt(ir2*ar2inv); // r^2-grid formula
            // r^2-grid inversion:      i0 = int(r2*ar2)
            //    --> interpolate between i0 and i0 + 1
        } // ir2
        return r;
    } // r_axis
 
//   status_t all_tests(int const echo=0);
  
} // namespace radial_r2grid
