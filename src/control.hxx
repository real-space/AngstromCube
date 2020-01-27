#pragma once

typedef int status_t;

namespace control {

  int constexpr default_echo_level = 2;

  status_t cli(char const *statement, int const echo=default_echo_level); // command line interface

  char const* set(char const *name, char const * value, int const echo=default_echo_level);
  char const* set(char const *name, double const value, int const echo=default_echo_level);

  char const* get(char const *name, char const * default_value, int const echo=default_echo_level);
  double      get(char const *name, double const default_value, int const echo=default_echo_level);

  status_t all_tests(int const echo=0);

} // namespace control
