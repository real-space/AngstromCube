#pragma once

#include <cstdint> // int64_t, int32_t, uint32_t

#ifndef NO_UNIT_TESTS
  #include <cstdio> // std::printf
  #include "simple_math.hxx" // ::random
#endif // NO_UNIT_TESTS

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

namespace global_coordinates {
  // global coordinates make a unique identifyer from each tuple of 3D integer coordinates.
  // Each coordinate is in the range [0, 2^21)
  // Negative values are allowed as input but will be folded back into the positive range.

  inline int64_t get(int32_t x, int32_t y, int32_t z) {
      // interleave bit pattern of the lowest 21 bits to 63 bits:
      // result: sign,z20,y20,x20,z19,y19,x19, ... ,z1,y1,x1,z0,y0,x0
      // with sign == 0 for any valid global coordinate
      int64_t i63{0};
      for (int b21 = 0; b21 < 21; ++b21) {
          int64_t const i3 = (x & 0x1) + 2*(y & 0x1) + 4*(z & 0x1);
          int64_t const i3_shifted = i3 << (b21*3);
          i63 |= i3_shifted;
          x >>= 1; y >>= 1; z >>= 1; // delete the least significant bit of the 3 coordinates
      } // b21
      return i63; // when printing i63, use format %22.22lo
  } // get

  template <typename int_t>
  inline int64_t get(int_t const xyz[3]) { return get(xyz[0], xyz[1], xyz[2]); }

  inline status_t get(uint32_t xyz[3], int64_t i63) {
      // retrieve global_coordinates from i63 identifyer
      // coordinates will be cast into the half-open range [0, 2^21)
      uint32_t x{0}, y{0}, z{0};
      for (int b21 = 0; b21 < 21; ++b21) {
          x |= (i63 & 0x1) << b21;   i63 >>= 1;
          y |= (i63 & 0x1) << b21;   i63 >>= 1;
          z |= (i63 & 0x1) << b21;   i63 >>= 1;
      } // b21
      xyz[0] = x; xyz[1] = y; xyz[2] = z;
      return status_t(i63 & 0x1); // this is 0 if i63 >= 0, i.e. if i63 was valid
  } // get

  inline int32_t get(uint32_t const xyz) {
      // retrieve signed global_coordinates in the half-open range [-2^20, 2^20)
      uint32_t constexpr b20 = 1 << 20;
      int32_t  constexpr b21 = 1 << 21;
      return int32_t(xyz) - (xyz >= b20)*b21; // subtract 2^21 if larger equal 2^20
  } // get

  inline status_t get(int32_t xyz[3], int64_t i63) {
      // retrieve global_coordinates from i63 identifyer
      // coordinates will be cast into the half-open range [-2^20, 2^20)
      uint32_t uxyz[3];
      auto const stat = get(uxyz, i63);
      for (int d = 0; d < 3; ++d) {
          xyz[d] = get(uxyz[d]);
      } // d
      return stat;
  } // get








#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_global_coordinates(int const echo=0) {
      status_t stat(0);
      int const n_tested = (1 << 12);
      for (int i = 0; i < n_tested; ++i) {
          int32_t ixyz[3]; // input coordinates
          for (int d = 0; d < 3; ++d) {
              ixyz[d] = simple_math::random(0, 0x1fffff); // [0, 2^21 - 1]
          } // d
          int64_t const i63 = get(ixyz);
          uint32_t oxyz[3]; // output coordinates
          status_t is = get(oxyz, i63);
          // the output format for i63 indices is octal with 22 digits.
          // exactly 21 octal digits represent the bit-interleaved 3D coordinates
          // and one leading 0 indices the octal notation for non-negative i63
          if (echo > 11) std::printf("# global_coordinates(%i, %i, %i)\t--> %22.22llo --> (%i, %i, %i)\n",
                                  ixyz[0], ixyz[1], ixyz[2], i63, oxyz[0], oxyz[1], oxyz[2]);
          for (int d = 0; d < 3; ++d) {
              is += (ixyz[d] != oxyz[d]);
          } // d
          if (is != 0) {
              if (echo > 1) std::printf("# global_coordinates(%i, %i, %i)\t--> %22.22llo --> (%i, %i, %i) failed!\n",
                                  ixyz[0], ixyz[1], ixyz[2], i63, oxyz[0], oxyz[1], oxyz[2]);
              ++stat;
          } // is
      } // i
      if (echo > 3) std::printf("# %s tested %.3f k random coordinate tuples, found %d errors\n", 
                              __func__, n_tested*.001, int(stat));
      if (echo > 9) {
          int64_t const i63 = -1; // show how invalid i63 indices are displayed
          std::printf("# global_coordinates(impossible coordinates)\t--> %22.22llo == %lld\n", i63, i63);
      } // echo
      return stat;
  } // test_global_coordinates

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_global_coordinates(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace global_coordinates
