/*
 *    Green function code, unit tests
 *    for the real-space grid based single-particle Hamiltonians
 *    as produced by density functional theory (DFT)
 */

#include <vector> // std::vector
#include <string> // std::string
#include <tuple> // std::tuple<...>

namespace green_tests {

  void add_tests(
        std::vector<std::tuple<char const*, double, status_t>> & results // results of module tests: <module_name, time_needed, status_returned>
      , std::string const & input_name // name to be matched against module_name
      , bool const show //
      , bool const all // all modules need to be tested or if (show==true) to display their names
      , int const echo // verbosity
  ); // declaration only

} // namespace green_tests
