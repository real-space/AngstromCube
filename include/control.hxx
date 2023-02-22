#pragma once

#include "status.hxx" // status_t

namespace control {

  int constexpr default_echo_level = 2;

  status_t command_line_interface(char const *statement, int const iarg, int const echo=default_echo_level);

  void set(char const *name, char const *value, int const echo=default_echo_level);

  char const* get(char const *name, char const * default_value);
  double      get(char const *name, double const default_value);
  double      get(double vec[], char const *name, char const *xyz="xyz", double const default_value=0);

  status_t read_control_file(char const *filename, int const echo=0);

  status_t show_variables(int const echo=default_echo_level);

  status_t all_tests(int const echo=0); // declaration only

} // namespace control
