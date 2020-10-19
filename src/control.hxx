#pragma once

#include "status.hxx" // status_t

namespace control {

  int constexpr default_echo_level = 2;

  status_t command_line_interface(char const *statement, int const echo=default_echo_level);

  char const* set(char const *name, char const * value, int const echo=default_echo_level);
  char const* set(char const *name, double const value, int const echo=default_echo_level);

  char const* get(char const *name, char const * default_value);//, int const echo=default_echo_level);
  double      get(char const *name, double const default_value);//, int const echo=default_echo_level);

  status_t read_control_file(char const *filename, int const echo=9);

  status_t all_tests(int const echo=0);

} // namespace control
