#pragma once

#include <cstdint> // int8_t

#include "status.hxx" // status_t

namespace chemical_symbol {

  int8_t decode(char const S, char const y); // declaration only
  int8_t get(char Sy[4], double const Z, char const blank='\0'); // declaration only, blank=' ' makes all symbols have a length of 3 non-zero chars

  status_t all_tests(int const echo=0); // declaration only

} // namespace chemical_symbol
