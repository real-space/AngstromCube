#pragma once

#include <cstdio> // printf

  template<typename T>
  int printf_vector( // returns the total number of chars written
        char const *const format // printf-format, should contain '%' only once
      , T const *const vec // pointer to vector start
      , int const n // number of elements to write
      , char const *const final="\n"
      , T const scale=T(1) // scale the vector components
      , T const add=T(0)   // add to the scaled vector components
  ) {
      // write a vector with increment +1 to stdout

      int n_chars_written{0};
      for(int i{0}; i < n; ++i) {
          n_chars_written += printf(format, vec[i]*scale + add);
      } // i
      if (final) n_chars_written += printf(final);
      return n_chars_written;
  } // printf_vector

  template <typename real_t>
  double print_stats(
        real_t const values[] // input values
      , size_t const all // how many
      , double const dV=1 // volume element for the integration
      , char const *prefix="" // leading printf messages
      , real_t const unit=1 // unit conversion factor
  ) {
      real_t gmin{9e307}, gmax{-gmin}; double gsum{0}, gsum2{0};
      for(size_t i = 0; i < all; ++i) {
          gmin = std::min(gmin, values[i]);
          gmax = std::max(gmax, values[i]);
          gsum  += values[i];
          gsum2 += pow2(values[i]);
      } // i
      printf("%s grid stats min %g max %g integral %g avg %g\n", prefix, gmin*unit, gmax*unit, gsum*dV*unit, gsum/all*unit);
      return gsum*dV;
  } // print_stats
