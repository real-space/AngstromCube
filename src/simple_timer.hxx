#pragma once

#include <chrono>
#include <cstring> // strcpy

  typedef int status_t;

  class SimpleTimer {
    private:
      std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
      char file[60];
      int  line;
      char func[32];
    public:

      SimpleTimer(char const *sourcefile, int const sourceline=0, char const *function=nullptr) {
          start_time = std::chrono::high_resolution_clock::now(); // start
          strcpy(file, sourcefile);
          line = sourceline;
          if (nullptr != function) strcpy(func, function); else func[0] = 0;
      } // constructor

      ~SimpleTimer() {
          auto const stop_time = std::chrono::high_resolution_clock::now(); // stop
          auto const musec = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time).count();
          printf("# timer started at %s:%d %s took %.5f sec\n", file, line, func, 1e-6*musec);
      } // destructor

  }; // class SimpleTimer
  
namespace simple_timer {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  inline int64_t fibonacci(int64_t const n) { if (n < 3) return 1; return fibonacci(n - 1) + fibonacci(n - 2); }

  inline status_t test_usage() {
    int64_t const inp = 42;
    int64_t result;
    {   SimpleTimer timer(__FILE__, __LINE__);
        result = fibonacci(inp);
    } // timer
    printf("# fibonacci(%ld) = %ld\n", inp, result);
    return 0;
  } // test_usage

  inline status_t all_tests(int const echo=3) {
    if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
    auto status = 0;
    status += test_usage();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace simple_timer
