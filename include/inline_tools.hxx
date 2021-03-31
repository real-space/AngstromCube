#pragma once

#include <cstdint> // uint64_t, int64_t
#include <cstddef> // size_t

#include "status.hxx" // status_t

  template <int nBits> inline 
  size_t align(int64_t const in) { return ((((in - 1) >> nBits) + 1) << nBits); }
// template<int nBits> size_t inline align(int64_t const in) { return in; } // alignment switched off

/** currently not needed
  inline int required_bits(uint64_t const in) {
      int nbits{0};
      auto n = in - 1;
      while (n) { 
          ++nbits;
          n >>= 1;
      } // while
      return nbits;
  } // required_bits
*/
  
  template <typename real_t> inline char const * real_t_name(); // provide no implementation for the general case
  template <> inline char const * real_t_name<double>() { return "double"; }
  template <> inline char const * real_t_name<float> () { return "float"; }

namespace inline_tools {

  inline uint64_t constexpr eightchars(char const string[8]) { return *((uint64_t*)string); }
  inline void eightchars(char string[8], uint64_t const c8) { *((uint64_t*)string) = c8; }

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_align(int const echo=1) {
      status_t stat(0);
      for (int i = 0; i < 99; ++i) {
          stat += (align<0>(i) != i); // error if align<0> does not reproduce itself
      } // i
      for (int i = 1; i < (1 << 30); i *= 2) {
          if (echo > 15) std::printf("# align<%d>(%d) = %ld\n", 1, i, align<1>(i));
          stat += (align<0>(i) != i);
          stat += (align<1>(i) != i && i > 1);
          stat += (align<2>(i) != i && i > 2);
          stat += (align<3>(i) != i && i > 4);
      } // i
      stat += (align<1>(1) != 2);
      stat += (align<1>(3) != 4);
      stat += (align<1>(7) != 8);
      stat += (align<2>(1) != 4);
      stat += (align<2>(3) != 4);
      stat += (align<2>(7) != 8);
      stat += (align<3>(1) != 8);
      stat += (align<3>(3) != 8);
      stat += (align<3>(7) != 8);
      return stat;
  } // test_align

  inline status_t test_eightchars(int const echo=1) {
      return (eightchars("01234567") != 0x3736353433323130);
  } // test_eightchars

  inline status_t all_tests(int const echo=2) {
      status_t stat(0);
      stat += test_eightchars(echo);
      stat += test_align(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace inline_tools

