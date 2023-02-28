#include <cstdio> // std::printf, ::fprintf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <cstdint> // size_t, uint32_t
#include <algorithm> // std::max, ::min, ::swap

#include "export_upf.hxx"

#include "inline_math.hxx" // set, pow2
#include "recorded_warnings.hxx" // warn
#include "constants.hxx" // ::pi
#include "print_tools.hxx" // printf_vector
#include "control.hxx" // ::get

namespace export_upf {

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test1(int const echo=0) {
      status_t stat(0);

      return 0;
  } // test1

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test1(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace export_upf
