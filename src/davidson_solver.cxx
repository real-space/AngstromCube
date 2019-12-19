#include <cstdio> // printf
#include <cassert> // assert
 
#include "davidson_solver.hxx"

// #define FULL_DEBUG
#define DEBUG

#ifdef  DEBUG
    #include "debug_output.hxx" // dump_to_file
#endif

#ifdef FULL_DEBUG
    #define full_debug(print) print 
#else
    #define full_debug(print)
#endif

#ifdef DEBUG
    #define debug(print) print 
#else
    #define debug(print)
#endif


namespace davidson_solver {
  // solve iteratively for the lowest eigenstates of an implicitly given Hamiltonian using the Davidson method
  
  template<typename real_t, int D0>
  status_t solve(real_t psi[] // on entry start wave functions, on exit improved eigenfunctions
    , size_t const nall // number of real_t in psi
    , int const echo=9 // log output level
  ) {
      return -1;
  } // solve

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t all_tests() {
    auto status = 0;
    printf("\n\n# %s: ToDo: implement tests\n", __FILE__);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace davidson_solver
