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
#else  // HAS_NO_OMP
  #include <omp.h> // ...
#endif // HAS_NO_OMP

#include "status.hxx" // status_t

namespace omp_parallel {
  // Test the usage of OpenMP

  inline status_t all_tests(int const echo=0) {

      auto const nthreads = omp_get_num_threads(), mthreads = omp_get_max_threads();
      if (echo > 2) std::printf("# OpenMP with %d of max %d threads\n", nthreads, mthreads);
      int sum{0};
      #pragma omp parallel for reduction(+,sum)
      for (int i = 0; i < nthreads; ++i) {
          auto const thread_id = omp_get_thread_num();
          if (echo > 4) std::printf("# OpenMP thread#%i works on item#%i\n", thread_id, i);
          ++sum;
      } // i
      assert(nthreads == sum);

      #pragma omp parallel
      {
          auto const thread_id = omp_get_thread_num();
          auto const parallel  = omp_in_parallel();
          if (0 == thread_id && echo > 3) std::printf("# OpenMP %sactive\n", parallel?"":"in");
      } // parallel region

      return 0;
  } // all_tests

} // namespace omp_parallel
