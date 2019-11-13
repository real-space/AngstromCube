#pragma once

#include <cstring> // strchr

typedef int status_t;

namespace control {

  int constexpr default_echo_level = 2;

  char const* set(char const *name, char const * value, int const echo=default_echo_level);
  char const* set(char const *name, double const value, int const echo=default_echo_level);

  char const* get(char const *name, char const * default_value, int const echo=default_echo_level);
  double      get(char const *name, double const default_value, int const echo=default_echo_level);

  inline char* find_equal_sign(char const * string) { return (char*)strchr(string, '='); }

  status_t all_tests();

} // namespace control
