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
        std::vector<std::tuple<char const*, double, status_t>> & results
      , std::string const & input_name
      , bool const show
      , bool const all
      , int const echo
  ); // declaration only

} // namespace green_tests
