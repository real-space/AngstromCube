#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <algorithm> // std::min, ::max
#include <cmath> // std::sqrt
#include <string> // std::string

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

namespace simple_stats {

  template <typename real_t=double>
  class Stats {
    public:

    Stats(int const value=0) { set(); } // default constructor

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

//  int allreduce(MPI_Comm const comm=MPI_COMM_WORLD) --> moved to mpi_parallel.hxx

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

    real_t min() const { return mini*(times > 0); }
    real_t max() const { return maxi*(times > 0); }
    real_t num() const { return v[0]; }
    real_t sum() const { return v[1]; }
    size_t tim() const { return times; }
    double mean() const { return (v[0] > 0) ? v[1]/double(v[0]) : 0.0; }
    double variance() const {
        auto const mu = mean();
        return (times > 1 && v[0] > 0) ? std::max(0.0, v[2]/v[0] - mu*mu) : 0.;
    } // variance
    double dev() const { return std::sqrt(variance()); } // standard deviation

    std::string interval(double const factor=1) const {
        char buffer[96];
        std::snprintf(buffer, 96, "[%g, %g +/- %g, %g]",
            min()*factor, mean()*factor, dev()*factor, max()*factor);
        return std::string(buffer);
    } // interval

    private:
      size_t times;
      real_t mini, maxi;
      real_t v[3];
  }; // class Stats<real_t>





#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <typename real_t>
  inline status_t test_basic(int const echo=0, int offset=0, double const threshold=1e-6) {
      Stats<real_t> s;
      int const begin=offset, end=offset + 100;
      double const ref[] = {49.5 + offset, 833.25}; // {mean, variance} of integers in range [0, 99]
      for (int i = begin; i < end; ++i) {
          s.add(i);
      } // i
      auto const mean = s.mean();
      if (echo > 3) std::printf("\n# %s: from %d to %d: %g +/- %g\n", __func__, begin, end - 1, mean, s.dev());
      auto const dev_mean = std::abs(ref[0] - mean),
             dev_variance = std::abs(ref[1] - s.variance());
      if (echo > 7) std::printf("# %s: dev_mean= %g, dev_variance= %g\n", __func__, dev_mean, dev_variance);
      if (echo > 9) std::printf("# Stats<%s>: %ld Byte\n", (4 == sizeof(real_t))?"float":"double", sizeof(s));
      return (dev_mean > threshold*mean) + (dev_variance > 2*threshold*mean*mean);
  } // test_basic

  inline status_t all_tests(int const echo=0) {
      if (echo > 0) std::printf("\n# %s %s\n", __FILE__, __func__);
      status_t stat(0);
      stat += test_basic<float >(echo,         0, 4e-7);
      stat += test_basic<double>(echo,         0, 2e-16);
      stat += test_basic<float >(echo, 1000*1000, 4e-7);
      stat += test_basic<double>(echo, 1000*1000, 2e-16);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace simple_stats
