#pragma once

#include <cstdint> // int8_t

#include "status.hxx" // status_t

namespace chemical_symbol {

  int8_t decode(char const S, char const y);
  status_t get(char* Sy, double const Z, char const blank=0); // blank=' ' makes all symbols have a length of 2 chars

  status_t all_tests(int const echo=0);

} // namespace chemical_symbol
