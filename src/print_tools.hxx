#pragma once

#include <cstdio> // printf

  template<typename T>
  int printf_vector(char const *const fmt, T const *const vec, int const n
                , T const scale=T(1), T const add=T(0), char const *const final="\n") {
      int n_chars_written{0};
      for(int i{0}; i < n; ++i) {
          n_chars_written += printf(fmt, vec[i]*scale + add);
      } // i
      if (final) n_chars_written += printf(final);
      return n_chars_written;
  } // printf_vector
