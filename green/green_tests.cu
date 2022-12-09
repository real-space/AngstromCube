/*
 *    Green function code, unit tests
 *    for the real-space grid based single-particle Hamiltonians
 *    as produced by density functional theory (DFT)
 */

#ifndef NO_UNIT_TESTS
  #include "simple_timer.hxx" // SimpleTimer
  // the following header files contain CUDA code
  #include "green_experiments.hxx" // ::all_tests
  #include "green_potential.hxx" // ::all_tests
  #include "green_function.hxx" // ::all_tests
  #include "green_kinetic.hxx" // ::all_tests
  #include "green_sparse.hxx" // ::all_tests
  #include "green_dyadic.hxx" // ::all_tests
  #include "green_action.hxx" // ::all_tests
#endif // not NO_UNIT_TESTS

#include "green_tests.hxx"

#include <cstdio> // std::printf
#include <vector> // std::vector
#include <string> // std::string
#include <tuple> // std::tuple<...>, ::make_tuple, ::get

namespace green_tests {

  void add_tests(
        std::vector<std::tuple<char const*, double, status_t>> & results
      , std::string const & input_name
      , bool const show
      , bool const all
      , int const echo
  ) {

#ifndef  NO_UNIT_TESTS

      { // testing scope

#include "add_module_test.h" // macro definition of add_module_test(MODULE_NAME), needs SimpleTimer, std::printf

          add_module_test(green_sparse);
          add_module_test(green_kinetic);
          add_module_test(green_potential);
          add_module_test(green_dyadic);
          add_module_test(green_action);
          add_module_test(green_function);
          add_module_test(green_experiments);

#undef    add_module_test

      } // testing scope

#endif // NO_UNIT_TESTS

  } // add_tests

} // namespace green_tests
