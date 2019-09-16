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
 
#include "grid_operators.hxx"

// #include "vector_math.hxx" // vector_math from exafmm
// #include "constants.hxx" // pi, sqrtpi
#include "real_space_grid.hxx" // grid_t
#include "finite_difference.hxx" // finite_difference_t


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


namespace grid_operators {
  // setup of the real-space grid-based Hamiltonian and overlap operator

  template<typename real_t, int D0>
  status_t Hamiltonian(real_t Hpsi[], real_space_grid::grid_t<D0> const &g, real_t const psi[], 
                       int const na, int const strides[], double const **atom_part, double const kvec[3],
                       finite_difference::finite_difference_t<real_t> const &fd, double const *potential=nullptr) {
      return -1;
  } // Hamiltonian

  template<typename real_t, int D0>
  status_t Overlap(real_t Spsi[], real_space_grid::grid_t<D0> const &g, real_t const psi[],
                       int const na, int const strides[], double const **atom_part, double const kvec[3]) {
      return -1;
  } // Overlap operator

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t all_tests() {
    auto status = 0;
    printf("\n\n# %s: ToDo: implement tests\n", __FILE__);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace grid_operators
