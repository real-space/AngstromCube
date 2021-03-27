#pragma once

#include <algorithm> // std::min, std::max
#include <cmath> // std::sqrt

namespace simple_stats {

  template <typename real_t=double>
  class Stats {
    public:

    Stats() { clear(); } // default constructor
    
    void clear() { mini = 1e38; maxi = -1e38; for(int p = 0; p < 4; ++p) v[p] = 0; times = 0; }

    void add(real_t const x, real_t const weight=1) {
        auto const w8 = std::abs(weight);
        v[0] += w8;
        v[1] += w8*x;
        v[2] += w8*x*x;
        v[3] += w8*x*x*x; // not used so far, beware of overflows
        maxi = std::max(maxi, x);
        mini = std::min(mini, x);
        ++times;
    } // add

    real_t min() const { return mini; }
    real_t max() const { return maxi; }
    real_t num() const { return v[0]; }
    real_t sum() const { return v[1]; }
    real_t avg() const { return (v[0] > 0) ? v[1]/v[0] : 0; } // aka mean
    real_t var() const { 
      auto const mean = avg();
      return (times > 0 && v[0] > 0) ?
          std::sqrt(std::max(real_t(0), v[2]/v[0] - mean*mean)) : 0;
    } // sqrt(variance)

    private:
      real_t v[4];
      real_t mini, maxi;
      size_t times;
  }; // class Stats<T>

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t all_tests(int const echo=0) {
    if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
    status_t status(0);
//  status += ToDo
    return status;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace simple_stats