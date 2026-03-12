#pragma once
// This file is part of AngstromCube under MIT License

// #doc
// The chemical_symbol module offers functionality to translate
// a number of protons (which may be non-integer) into a chemical species
// and offers chemical two-letter symbols as well as decoding of those.
//

#include <cstdint> // int8_t
#include <cmath> // std::round

#include "status.hxx" // status_t

namespace chemical_symbol {

    // translate a number of protons into a species number
    inline int8_t get(double const Z) { return int(std::round(Z)); }

    // species number for unidentifiable symbols
    int8_t constexpr SYMBOL_UNKNOWN = -128;

    // translate two letters into a species number in [0, 127] or SYMBOL_UNKNOWN
    int8_t decode(char const S, char const y); // declaration only

    // retrieve the symbol, blank=' ' makes all symbols have a length of 3 non-zero chars
    int8_t get(char Sy[4], double const Z, char const blank='\0'); // declaration only

    // self-tests
    status_t all_tests(int const echo=0); // declaration only

} // namespace chemical_symbol
