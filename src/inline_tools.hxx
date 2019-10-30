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


// template<typename T> 
// T** inline rectangular_array(size_t const n1, size_t const n0) {
//     auto const mem = new T[n1*n0];
//     auto const ptr = new (T*)[n1];
//     for(int i1 = 0; i1 < n1; ++i1) ptr[i1] = &mem[i1*n0]; 
//     return ptr;
// } // rectangular_array --- difficult to clean up
