#pragma once

#include <cstdint> // uint64_t, int64_t
#include <cstddef> // size_t

#include "status.hxx" // status_t

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

template <typename real_t> inline char const * real_t_name(); // provide no implementation for the general case
template <> inline char const * real_t_name<double>() { return "double"; }
template <> inline char const * real_t_name<float> () { return "float"; }

namespace inline_tools {

	uint64_t constexpr inline eightchars(char const string[8]) { return *((uint64_t*)string); }
	void inline eightchars(char string[8], uint64_t const c8) { *((uint64_t*)string) = c8; }

#ifdef  NO_UNIT_TESTS
	status_t inline all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

 	status_t inline test_eightchars(int const echo=1) {
	    return (eightchars("01234567") != 0x3736353433323130);
  	} // test_eightchars

 	status_t inline all_tests(int const echo=2) {
    	status_t stat{0};
    	stat += test_eightchars(echo);
	    return stat;
  	} // all_tests

#endif // NO_UNIT_TESTS

} // namespace inline_tools

