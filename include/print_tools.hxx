#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <vector> // std::vector<T>
#include <algorithm> // std::min, ::max

#include "inline_math.hxx" // pow2

  template <typename T>
  int printf_vector( // returns the total number of chars written
        char const *const format // printf-format, should contain '%' only once
      , T const vec[] // pointer to vector start
      , int const n // number of elements to write
      , char const *const final="\n"
      , T const scale=T(1) // scale the vector components
      , T const add=T(0)   // add to the scaled vector components
  ) {
      // write a vector with increment +1 to stdout
      int n_chars_written{0};
      for (int i{0}; i < n; ++i) {
          n_chars_written += std::printf(format, vec[i]*scale + add);
      } // i
      if (final) n_chars_written += std::printf("%s", final);
      return n_chars_written;
  } // printf_vector

  template <typename T>
  int printf_vector( // returns the total number of chars written
        char const *const format // printf-format, should contain '%' only once
      , std::vector<T> const & vec // pointer to vector start
      , char const *const final="\n"
      , T const scale=T(1) // scale the vector components
      , T const add=T(0)   // add to the scaled vector components
  ) {
      return printf_vector(format, vec.data(), vec.size(), final, scale, add);
  } // printf_vector

  template <typename real_t>
  double print_stats(
        real_t const values[] // input values
      , size_t const all // how many
      , double const dV=1 // volume element for the integration
      , char const *prefix="" // leading printf messages
      , double const unit=1 // unit conversion factor
      , char const *_unit="" // unit indicator
  ) {
      double gmin{9e307}, gmax{-gmin}, gsum{0}, gsum2{0};
      for (size_t i = 0; i < all; ++i) {
          gmin = std::min(gmin, double(values[i]));
          gmax = std::max(gmax, double(values[i]));
          gsum  += values[i];
          gsum2 += pow2(values[i]);
      } // i
      std::printf("%s grid stats min %g max %g avg %g", prefix, gmin*unit, gmax*unit, gsum/all*unit);
      if (dV > 0) std::printf(" integral %g", gsum*dV*unit);
      std::printf(" %s\n", _unit);
      return gsum*dV;
  } // print_stats

  template <typename int_t=int32_t>
  class SparsifyPlot {
    // This helper class allows to loop over a histogram
    // and plotting non-zero entries and, in addition, those entries
    // before and after a non-zero entry.
    // Usage:
    //    SparsifyPlot<int> spp;
    //    for (int_t i = 0; i < hist.size(); ++i) {
    //        auto const j = spp.previous(i, hist[i] > 0);
    //        if (j >= 0) printf("%d %g\n", j, hist[j]);
    //    } // i
    // Drawback: the last entry (nhist-1) is never plotted!
  public:

      SparsifyPlot(bool const show0=false) {
          assert(int_t(-1) < 0 && "int_t must be a signed integer type");
          // show the first data point always, -1: do not show
          last[0] = show0 ? 0 : -1; last[1] = -1; last[2] = -1; last[3] = -1;
      } // constructor

      int_t previous(int_t const ih, bool const nonzero=true) {
          last[(ih - 2) & 0x3] = -1; // clear after usage
          if (nonzero) {
              last[(ih - 1) & 0x3] = ih - 1;
              last[(ih    ) & 0x3] = ih; // also plot indices left and right of ih
              last[(ih + 1) & 0x3] = ih + 1;
          } // nonzero
          return last[(ih - 1) & 0x3]; // returns previous index
      } // previous

  private:
      int_t last[4]; // state
  }; // class SparsifyPlot
