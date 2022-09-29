#pragma once

#include <chrono> // std::chrono::high_resolution_clock
#include <cstring> // std::strcpy
#include <cstdint> // int32_t
#include <cstdio> // std::printf
#include <algorithm> // std::max

#include "status.hxx" // status_t

  class ProgressReport {
    // This timer object prints the time elapsed between construction and the last iteration
    // if the time is larger than a given threshold
    // furthermore, it prints a report every time that the threshold is reached

    private:
      char  file[56];
      int32_t line;
      int32_t echo;
      std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
      std::chrono::time_point<std::chrono::high_resolution_clock> last_time;
      double delta;
    public:

      ProgressReport(char const *sourcefile, int const sourceline=0, double const interval=60, int const echo=1) 
      : line(sourceline), echo(echo), delta(interval) {
          std::strcpy(file, sourcefile);
          start_time = std::chrono::high_resolution_clock::now(); // start
          last_time = start_time;
      } // constructor

      double report(size_t const iteration, size_t const n_iterations) {
          auto const now = std::chrono::high_resolution_clock::now();
          auto const total = std::chrono::duration_cast<std::chrono::microseconds>(now - start_time).count()*1e-6;
          auto const since = std::chrono::duration_cast<std::chrono::microseconds>(now -  last_time).count()*1e-6;
          auto const average = total/(iteration + 1.);
          if (n_iterations - 1 == iteration) { // last iteration done
              auto const performance = (average > 0) ? 1./average : 0; 
              if (total > delta && echo > 0) {
                  std::printf("# timer started at %s:%d took %g sec for %ld iterations, %g sec/iteration, %g iterations/sec\n",
                                                  file, line, total, n_iterations,      average,          performance);
              }
          } else
          if (since > delta) {
              // estimated time of arrival (ETA)
              auto const done = (iteration + 1.)/std::max(1., 1.*n_iterations), left = 1 - done;
              auto const eta = total*left/done;
              if (echo > 0) {
                  std::printf("# timer started at %s:%d took %g sec for %ld of %ld iterations (%.2f %%), expect %3g sec more\n",
                                                  file, line, total, iteration + 1, n_iterations, done*100, eta);
              }
              last_time = now;
          }
          return total;
      } // report

  }; // class ProgressReport







namespace progress_report {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline int64_t fibonacci(int64_t const n) {
      if (n < 3) return 1;
      return fibonacci(n - 1) + fibonacci(n - 2);
  } // fibonacci

  inline status_t test_basic_usage(int const echo=3) {
      int64_t result;
      int64_t const inp = 40, reference = 102334155;
      { // scope: create a timer, do some iterations, destroy the timer
          double const every = 1.0; // report to stdout every second
          ProgressReport timer(__FILE__, __LINE__, every, echo);
          int const nits = 10;
          for (int it = 0; it < nits; ++it) {
              result = fibonacci(inp);
              timer.report(it, nits);
          } // it
      } // scope
      if (echo > 0) std::printf("# fibonacci(%lld) = %lld\n", inp, result);
      return (reference != result);
  } // test_basic_usage

  inline status_t all_tests(int const echo=0) {
      if (echo > 2) std::printf("\n# %s %s\n\n", __FILE__, __func__);
      status_t stat(0);
      stat += test_basic_usage(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace progress_report
