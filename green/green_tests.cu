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
  #include "green_input.hxx" // ::all_tests
#endif // not NO_UNIT_TESTS

#include <cstdio> // std::printf
#include <cassert> // assert
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

#define   add_module_test(MODULE_NAME) {                                            \
              auto const module_name = #MODULE_NAME;                                \
              if (all || (input_name == module_name)) {                             \
                  SimpleTimer timer(module_name, 0, "", 0);                         \
                  if (echo > 3) std::printf("\n\n\n# ============= Module test"     \
                     " for %s ==================\n\n", module_name);                \
                  auto const stat = show ? 0 : MODULE_NAME::all_tests(echo);        \
                  double const time = timer.stop();                                 \
                  results.push_back(std::make_tuple(module_name, time, stat));      \
              }                                                                     \
          } // add_module_test

          add_module_test(green_input);
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
