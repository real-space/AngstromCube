#pragma once
// This file is part of AngstromCube under MIT License

#include "status.hxx" // status_t

namespace control {
  /*
   *    This is an internal variable environment that facilitates code development.
   *    Variables can be defined in an input file, the command line or is the code.
   *    Variables can be read anywhere in the code.
   *    Internally, variables are stored as pairs of two strings, name and value.
   */

  int constexpr default_echo_level = 2;

  status_t command_line_interface(char const *statement, int const iarg // define a (name, value) pair
                                  , int const echo=default_echo_level); // by "name=value" syntax

  char const* set(char const *name, char const  *value, int const echo=default_echo_level); // define a (name, value) pair
  double      set(char const *name, double const value, int const echo=default_echo_level);
  int constexpr echo_set_without_warning = -1;

  char const* get(char const *name, char const * default_value); // read a string  value from the variable environment
  double      get(char const *name, double const default_value); // read a numeric value from the variable environment
  double      get(double vec[], char const *name
      , char const *xyz="xyz", double const default_value=0); // read several values for name.x, name.y, ...

  status_t read_control_file(char const *filename, int const echo=0); // read definitions given in an input file

  status_t show_variables(int const echo=default_echo_level); // show which variables were defined

  status_t all_tests(int const echo=0); // declaration only

} // namespace control
