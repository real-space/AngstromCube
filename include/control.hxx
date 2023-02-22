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

#if 1
  inline double get(
        double vec[] // result array
      , char const *const name // base keyword
      , char const *xyz="xyz" // characters to extract
      , double const default_value=0 // default value
  )
    // specify entire vectors, e.g. spacing.x, spacing.y, ... or coeff.a, coeff.b, ...
  {
      double const default_val = control::get(name, default_value);
      char name_x[96];
      for (int d = 0; '\0' != xyz[d]; ++d) {
          std::snprintf(name_x, 96, "%s.%c", name, xyz[d]);
          vec[d] = control::get(name_x, default_val);
      } // d
      return default_val;
  } // get
#endif

} // namespace control
