#pragma once

#include <cstdio> // printf, std::fprintf, std::FILE, std::fopen, std::fclose, std::fflush, stdout
#include "display_units.h" // Ang, _Ang
#include "complex_tools.hxx" // is_complex

#define here if (echo > 5) { printf("\n# here: %s %s:%i\n\n", __func__, __FILE__, __LINE__); std::fflush(stdout); }

template <typename real_t>
int dump_to_file(
      char const *filename
    , size_t const N // number of rows
    , real_t const y_data[] // 2D-data array [N*Stride]
    , double *x_axis=nullptr // x-axis[N], replaced by 0,1,2,... if 0x0
    , size_t const Stride=1 // Stride, default=1 for 1D functions
    , size_t const M=1 // number of columns to display
    , char const *title=nullptr // title line in the file
    , int const echo=0 // report writing to stdout, 0: suppress report
) {
    std::FILE *f = std::fopen(filename, "w");
    if (nullptr == f) {
        if (echo > 1) printf("# %s Error opening file %s!\n", __func__, filename);
        return 1;
    } // failed to open

    std::fprintf(f, "#%s %s\n", (is_complex<real_t>())?"complex":"", title); // print this line also if title==nullptr

    for(int i = 0; i < N; i++) {
        if (nullptr != x_axis) {
            std::fprintf(f, "%g ", x_axis[i]);
        } else {
            std::fprintf(f, "%d ", i);
        } // x_axis given
        for(int j = 0; j < M; ++j) {
            auto const y = y_data[i*Stride + j];
            if (is_complex<real_t>()) {
                std::fprintf(f, "  %g %g", std::real(y), std::imag(y));
            } else {
                std::fprintf(f, " %g", y);
            } // is_complex
        } // j
        std::fprintf(f, "\n"); // new line
    } // i
    std::fprintf(f, "\n"); // new line at the end of the file

    std::fclose(f);
    if (echo > 3) printf("# file %s written with %d x %d (of %d) data entries.\n", filename, N, M, Stride);
    return 0;
} // dump_to_file

namespace debug_output {

  template <typename real_t>
  status_t write_array_to_file(
        char const *filename // file name to write to
      , real_t const array[]  // array pointer
      , int const nx, int const ny, int const nz // grid dimensions
      , int const echo=0 // log-level
      , char const *arrayname="" // some description to appear in the file
  ) {
      char title[128]; std::sprintf(title, "%i x %i x %i  %s", nz, ny, nx, arrayname);
      auto const size = size_t(nz) * size_t(ny) * size_t(nx);
      return dump_to_file(filename, size, array, nullptr, 1, 1, title, echo);
  } // write_array_to_file
  
} // namespace debug_output
