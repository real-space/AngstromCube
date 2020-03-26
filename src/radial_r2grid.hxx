#pragma once

#include <cstdint> // uint32_t
#include <cmath> // std::sqrt
#include <vector> // std::vector<T>

// #include "status.hxx" // status_t


    template <typename real_t=float>
    std::vector<real_t> r2_axis(int const nr2, float const ar2=1) {
        real_t const ar2inv = 1.0/ar2;
        std::vector<real_t> r(nr2);
        for(int ir2 = 0; ir2 < nr2; ++ir2) {
            r[ir2] = std::sqrt(ir2*ar2inv); // r^2-grid
        } // ir2
        return r;
    } // r2_axis

#if 0
class radial_r2grid_t
{
public:
    radial_r2grid_t(int const n=0, float const a=1) : _ar2(a), _nr2(n) {} // constructor
    
    float    a() const { return _ar2; }
    uint32_t n() const { return _nr2; }
    float    a_inverse() const { return 1./_ar2; }
    
    template <typename real_t=float>
    std::vector<real_t> r_axis() const { return r2_axis(_nr2, _ar2); }

private:
    float    _ar2;
    uint32_t _nr2;
}; // class radial_r2grid_t
#endif

// namespace radial_r2grid {
 
//   status_t all_tests(int const echo=0);
  
// } // namespace radial_r2grid
