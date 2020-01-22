#pragma once

typedef int status_t;

#include <cstdio> // std::sprintf

#define warn(...) std::sprintf(recorded_warnings::_new_warning(__FILE__, __LINE__), __VA_ARGS__); 

namespace recorded_warnings {

  char* _new_warning(char const *file, int const line); // hidden, please use the macro above

  status_t show_warnings(int const echo=1);

  status_t clear_warnings(int const echo=1);

  status_t all_tests();

} // namespace recorded_warnings
