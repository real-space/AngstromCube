#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert

#ifdef    HAS_NO_OMP
  inline int  omp_get_max_threads() { return 1; }
  inline int  omp_get_num_threads() { return 1; }
  inline int  omp_get_thread_num()  { return 0; }
  inline int  omp_get_num_procs()   { return 1; }
  inline bool omp_in_parallel()     { return 0; }

  bool constexpr omp_replacement = true;
#else  // HAS_NO_OMP
  bool constexpr omp_replacement = false;
  #include <omp.h> // ...
#endif // HAS_NO_OMP

#include "status.hxx" // status_t

namespace omp_parallel {
  // Test the usage of OpenMP

  inline status_t all_tests(int const echo=0) {

      auto const max_threads = omp_get_max_threads(); // can be controlled by shell environment variable OMP_NUM_THREADS
      if (echo > 2) std::printf("# %sOpenMP with max %d threads\n", omp_replacement?"fake ":"", max_threads);
      assert(1 == omp_get_num_threads()); // outside a parallel region, only 1 thread should be running
      int sum{0};
      #pragma omp parallel for reduction(+:sum)
      for (int i = 0; i < max_threads; ++i) {
          auto const num_threads = omp_get_num_threads();
          auto const thread_id   = omp_get_thread_num();
          if (echo > 4) std::printf("# OpenMP thread#%i of %d works on item#%i\n", thread_id, num_threads, i);
          ++sum;
      } // i
      assert(max_threads == sum);

      #pragma omp parallel
      {
          auto const thread_id = omp_get_thread_num();
          auto const parallel  = omp_in_parallel();
          if (0 == thread_id && echo > 3) std::printf("# OpenMP %sactive\n", parallel?"":"in");
      } // parallel region

      return 0;
  } // all_tests

} // namespace omp_parallel
