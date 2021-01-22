#include <cstdio> // printf
#include <cassert> // assert
#include <numeric> // std::iota
#include <vector> // std::vector<T>

#include "density_generator.hxx"

#include "real_space.hxx" // ::grid_t
#include "grid_operators.hxx" // ::grid_operator_t, empty_list_of_atoms
#include "fermi_distribution.hxx" // ::FermiLevel_t

namespace density_generator {
 
#ifndef NO_UNIT_TESTS

  status_t test_init(int const echo=3) {
      real_space::grid_t const g(4, 5, 6);
      grid_operators::grid_operator_t<float,double> const op(g, grid_operators::empty_list_of_atoms());
      std::vector<float> wave(g.all());
      std::vector<double> rho(g.all());
      std::iota(wave.begin(), wave.end(), 0);
      fermi_distribution::FermiLevel_t eF(1);
      double spectrum[] = {0};
      return density(rho.data(), nullptr, eF, spectrum, wave.data(),
                     wave.data(), // this is only passed to match the interface
                     nullptr, op.get_natoms(), g, 1, 1, echo, -1, nullptr, nullptr);
  } // test_init

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_init(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace density_generator
