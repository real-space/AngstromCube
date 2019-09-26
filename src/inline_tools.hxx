#pragma once

#include <cstdint> // uint64_t

template<int nbits> size_t inline align(int64_t const in) 
         { return ((((in - 1) >> nbits) + 1) << nbits); }
// template<int nbits> size_t inline align(int64_t const in) { return in; } // alignment switched off

int inline required_bits(uint64_t const in) {
    int nbits = 0;
    auto n = in - 1;
    while (n) { 
        ++nbits;
        n >>= 1;
    } // while
    return nbits;
} // required_bits
