#include <cstdio> // printf
#include <cassert> // assert
#include <numeric> // std::iota
#include <vector> // std::vector<T>

#include "density_generator.hxx"

#include "real_space_grid.hxx" // ::grid_t
#include "grid_operators.hxx" // ::grid_operator_t

namespace density_generator {
 
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_init(int const echo=3) {
      real_space_grid::grid_t<1> const g(4, 5, 6);
      grid_operators::grid_operator_t<float,double,1> const op(g);
      std::vector<float> wave(g.all());
      std::vector<double> rho(g.all());
      std::iota(wave.begin(), wave.end(), 0);
      return density(rho.data(), nullptr, wave.data(), op, 1, 1, echo);
  } // test_init

  status_t all_tests(int const echo) {
    status_t status(0);
    status += test_init(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace density_generator