#include <cstdio> // printf
#include <cassert> // assert
#include <numeric> // std::iota
#include <vector> // std::vector<T>

#include "density_generator.hxx"

#include "real_space.hxx" // ::grid_t
#include "grid_operators.hxx" // ::grid_operator_t

namespace density_generator {
 
#ifndef NO_UNIT_TESTS

  status_t test_init(int const echo=3) {
      real_space::grid_t const g(4, 5, 6);
      grid_operators::grid_operator_t<float,double> const op(g);
      std::vector<float> wave(g.all());
      std::vector<double> rho(g.all());
      std::iota(wave.begin(), wave.end(), 0);
      return density(rho.data(), nullptr, wave.data(), op, 1, 1, echo);
  } // test_init

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_init(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace density_generator
