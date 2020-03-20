#pragma once

#include "status.hxx" // status_t

namespace conjugate_gradients {

  template<typename real_t> inline double tiny();
  template<> inline double tiny<double>() { return 2.25e-308; }
  template<> inline double tiny<float> () { return 1.18e-38; }

  status_t all_tests(int const echo=0);

} // namespace conjugate_gradients
