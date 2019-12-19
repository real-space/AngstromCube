#include <cstdio> // printf
// #include <cstdlib> // abs
// #include <cmath> // sqrt, exp
// #include <algorithm> // max
// #include <complex> // std::complex<real_t>
// #include <complex>
// #include <utility> // std::pair<T1,T2>, make_pair
// #include <vector> // std::vector<T>
// #include <array> // std::array<T,n>
// #include <cmath>
#include <cassert> // assert
 
#include "davidson_solver.hxx"

// #include "vector_math.hxx" // vector_math from exafmm
// #include "constants.hxx" // pi, sqrtpi


// #include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
// #include "display_units.h" // eV, _eV, Ang, _Ang

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
  // setup of the DFT Hamiltonian in SHO-basis
  
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
