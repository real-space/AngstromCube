#pragma once

#include "status.hxx" // status_t

namespace control {

  int constexpr default_echo_level = 2;

  status_t command_line_interface(char const *statement, int const echo=default_echo_level);

  char const* set(char const *name, char const * value, int const echo=default_echo_level);
  char const* set(char const *name, double const value, int const echo=default_echo_level);

  char const* get(char const *name, char const * default_value);
  double      get(char const *name, double const default_value);

  status_t read_control_file(char const *filename, int const echo=0);

  status_t show_variables(int const echo=1);

  status_t all_tests(int const echo=0); // declaration only

#if 0
  template <typename T>
  void control_get_xyz(
        T vec[]
      , int const n
      , char const *const keyword
      , double const default
      , char const x='x'
  )
    // specify entire vectors, e.g. spacing.x, spacing.y, ... or coeff.a, coeff.b, ...
  {
      double const default_value = control::get(keyword, default);
      char keyword_x[96];
      for (int d = 0; d < n; ++d) {
          std::snprintf(keyword_x, 95, "%s.%c", keyword, x + d);
          vec[d] = control::get(keyword_x, default_value);
      } // d
  } // control_get_xyz
#endif

} // namespace control
