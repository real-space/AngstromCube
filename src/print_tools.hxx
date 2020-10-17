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
