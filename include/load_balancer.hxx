#pragma once

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

namespace load_balancer {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS
  status_t all_tests(int const echo=0); // declaration only
#endif // NO_UNIT_TESTS

} // namespace load_balancer
