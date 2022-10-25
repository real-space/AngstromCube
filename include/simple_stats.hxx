#pragma once

#include <algorithm> // std::min, ::max
#include <cmath> // std::sqrt

#ifndef NO_UNIT_TESTS
    #include <cstdio> // std::printf
#endif // NO_UNIT_TESTS

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
// #include "mpi_parallel.hxx" // ::max, ::sum

namespace simple_stats {

  template <typename real_t=double>
  class Stats {
    public:

    Stats(int const value=0) { set(); } // default constructor

//     void clear() {
//         times = 0;
//         for (int p = 0; p < 3; ++p) {
//             v[p] = 0;
//         } // p
//         mini =  1.7e38;
//         maxi = -1.7e38;
//     } // clear

    void add(real_t const x, real_t const weight=1) {
        auto const w8 = std::abs(weight);
        v[0] += w8;
        v[1] += w8*x;
        v[2] += w8*x*x;
//      v[3] += w8*x*x*x; // not used so far, beware of overflows
        maxi = std::max(maxi, x);
        mini = std::min(mini, x);
        ++times;
    } // add

//     int allreduce(MPI_Comm const comm=MPI_COMM_WORLD) {
// #ifdef    HAS_NO_MPI
//         return 0;
// #else  // HAS_NO_MPI
//         real_t minmax[2] = {-mini, maxi};
//         auto const status_max = mpi_parallel::max(minmax, 2, comm);
//         auto const status_sum = mpi_parallel::sum(     v, 3, comm);
//         auto const status_tim = mpi_parallel::sum(&times, 1, comm);
//         mini = -minmax[0];
//         maxi =  minmax[1];
//         return status_max + status_sum + status_tim;
// #endif // HAS_NO_MPI
//     } // allreduce

    void get(double values[8]) const { // export values for MPI_Allreduce
        values[0] = v[0];
        values[1] = v[1];
        values[2] = v[2];
        values[3] = 0; // reserved for extension to 3rd cumulant
        values[4] = times;
        values[5] = 0; // not in use
        values[6] = -mini; // these need MPI_MAX
        values[7] =  maxi; // these need MPI_MAX
    } // get

    void set(double const values[8]=nullptr) { // import values after MPI_Allreduce
        if (nullptr != values) {
            times = values[4];
            mini = -values[6];
            maxi =  values[7];
            v[0] =  values[0];
            v[1] =  values[1];
            v[2] =  values[2];
        } else { // nullptr
            times = 0;
            float constexpr LIM = 1.7e38;
            mini =  LIM;
            maxi = -LIM;
            for (int p = 0; p < 3; ++p) v[p] = 0;
        } // nullptr
    } // set

    real_t min() const { return mini; }
    real_t max() const { return maxi; }
    real_t num() const { return v[0]; }
    real_t sum() const { return v[1]; }
    size_t tim() const { return times; }
    double mean() const { return (v[0] > 0) ? v[1]/double(v[0]) : 0.0; }
    double variance() const {
        auto const mu = mean();
        return (times > 0 && v[0] > 0) ?
            std::max(0.0, v[2]/v[0] - mu*mu) : 0.0;
    } // variance
    double dev() const { return std::sqrt(variance()); } // standard deviation

    private:
      size_t times;
      real_t mini, maxi;
      real_t v[3];
  }; // class Stats<real_t>





#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <typename real_t>
  inline status_t test_basic(int const echo=0) {
      Stats<real_t> s;
      int const begin=0, end=100; double const ref[] = {49.5, 833.25}; // {mean, variance}
      for (int i = begin; i < end; ++i) {
          s.add(i);
      } // i
      if (echo > 3) std::printf("# %s: from %d to %d: %g +/- %g\n",
                      __func__, begin, end - 1, s.mean(), s.dev());
      if (echo > 8) std::printf("# %s: %ld Byte\n", __FILE__, sizeof(s));
      return (ref[0] != s.mean()) + (ref[1] != s.variance());
  } // test_basic

  inline status_t all_tests(int const echo=0) {
      if (echo > 0) std::printf("\n# %s %s\n", __FILE__, __func__);
      status_t stat(0);
      stat += test_basic<float>(echo);
      stat += test_basic<double>(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace simple_stats
