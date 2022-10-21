#pragma once

#include "status.hxx" // status_t

namespace control {

  int constexpr default_echo_level = 2;

  status_t command_line_interface(char const *statement, int const iarg, int const echo=default_echo_level);

  void set(char const *name, char const *value, int const echo=default_echo_level);

  char const* get(char const *name, char const * default_value);
  double      get(char const *name, double const default_value);

  status_t read_control_file(char const *filename, int const echo=0);

  status_t show_variables(int const echo=default_echo_level);

  status_t all_tests(int const echo=0); // declaration only

#if 0
  inline void control_get_xyz(
        double vec[] // result array
      , char const *const keyword // root keyword
      , int const n=3 // number of elements
      , double const default=0 // default value
      , char const x='x' // character to start with
  )
    // specify entire vectors, e.g. spacing.x, spacing.y, ... or coeff.a, coeff.b, ...
  {
      double const default_value = control::get(keyword, default);
      char keyword_x[96];
      for (int d = 0; d < n; ++d) {
          std::snprintf(keyword_x, 96, "%s.%c", keyword, x + d);
          vec[d] = control::get(keyword_x, default_value);
      } // d
  } // control_get_xyz
#endif // 0

} // namespace control
