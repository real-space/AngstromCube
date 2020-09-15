#pragma once

#include <chrono>
#include <cstring> // strcpy

  #include "status.hxx" // status_t

  class SimpleTimer {
    // This timer object prints the time elapsed between construction and destruction 
    // into stdout immediately when destructed. See test_usage() below.
    private:
      std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
      char file[60];
      int  line, echo;
      char func[32];
    public:

      SimpleTimer(char const *sourcefile, int const sourceline=0, char const *function=nullptr, int const echo=1) 
      : line(sourceline), echo(echo) {
          strcpy(file, sourcefile);
          if (nullptr != function) strcpy(func, function); else func[0] = 0;
          start_time = std::chrono::high_resolution_clock::now(); // start
      } // constructor

      double stop(int const stop_echo=0) { 
          auto const stop_time = std::chrono::high_resolution_clock::now(); // stop
          auto const musec = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time).count();
          if (stop_echo > 0) printf("# timer started at %s:%d %s took %.5f sec\n", file, line, func, 1e-6*musec);
          return 1e-6*musec;
      } // stop

      ~SimpleTimer() {
          stop(echo);
      } // destructor

  }; // class SimpleTimer
  
namespace simple_timer {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  inline int64_t fibonacci(int64_t const n) { if (n < 3) return 1; return fibonacci(n - 1) + fibonacci(n - 2); }

  inline status_t test_usage(int const echo=3) {
    int64_t result, inp = 42;
    {   SimpleTimer timer(__FILE__, __LINE__, "comment=fibonacci", echo);
        result = fibonacci(inp);
    } // timer destructor is envoked at the end of this scope
    if (echo > 0) printf("# fibonacci(%ld) = %ld\n", inp, result);
    return 0;
  } // test_usage

  inline status_t all_tests(int const echo=0) {
    if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
    status_t status(0);
    status += test_usage(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace simple_timer
