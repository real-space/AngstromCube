#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, ::fprintf, ::fopen, ::fclose, ::fflush, stdout

#include "complex_tools.hxx" // is_complex
#include "display_units.h" // Ang, _Ang

#ifdef    DEBUG
    #define here \
        if (echo > 5) { \
            std::printf("\n# here: %s %s:%i\n\n", __func__, __FILE__, __LINE__); \
            std::fflush(stdout); \
        }
#else  // DEBUG
    #define here ;
#endif // DEBUG

template <typename real_t>
inline int dump_to_file(
      char const *filename
    , size_t const N // number of rows
    , real_t const y_data[] // 2D-data array [N*Stride]
    , double *x_axis=nullptr // x-axis[N], replaced by 0,1,2,... if 0x0
    , size_t const Stride=1 // Stride, default=1 for 1D functions
    , size_t const M=1 // number of columns to display
    , char const *title=nullptr // title line in the file
    , int const echo=0 // report writing to stdout, 0: suppress report
) {
    auto *const f = std::fopen(filename, "w");
    if (nullptr == f) {
        if (echo > 1) std::printf("# %s Error opening file %s!\n", __func__, filename);
        return 1;
    } // failed to open

    bool constexpr is_a_complex = is_complex<real_t>();
    std::fprintf(f, "#%s %s\n", is_a_complex?"complex":"", title); // print this line also if title==nullptr

    for (int i = 0; i < N; i++) {
        if (nullptr != x_axis) {
            std::fprintf(f, "%g ", x_axis[i]);
        } else {
            std::fprintf(f, "%d ", i);
        } // x_axis given
        for (int j = 0; j < M; ++j) {
            auto const y = y_data[i*Stride + j];
            if (is_a_complex) {
                std::fprintf(f, "  %g %g", double(std::real(y)), double(std::imag(y)));
            } else {
                std::fprintf(f, " %g", double(std::real(y)));
            } // is_a_complex
        } // j
        std::fprintf(f, "\n"); // new line
    } // i
    std::fprintf(f, "\n"); // new line at the end of the file

    std::fclose(f);
    if (echo > 3) std::printf("# file %s written with %lu x %lu (of %lu) data entries.\n", filename, N, M, Stride);
    return 0;
} // dump_to_file

namespace debug_output {

  template <typename real_t>
  inline status_t write_array_to_file(
        char const *filename // file name to write to
      , real_t const array[]  // array pointer
      , int const nx, int const ny, int const nz // grid dimensions
      , int const echo=0 // log-level
      , char const *arrayname="" // some description to appear in the file
  ) {
      char title[128]; std::snprintf(title, 128, "%i x %i x %i  %s", nz, ny, nx, arrayname);
      auto const size = size_t(nz) * size_t(ny) * size_t(nx);
      return dump_to_file(filename, size, array, nullptr, 1, 1, title, echo);
  } // write_array_to_file

  inline status_t all_tests(int echo=0) { return 0; }

} // namespace debug_output
