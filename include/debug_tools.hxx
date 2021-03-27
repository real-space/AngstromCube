#pragma once

#include <cstdio> // printf, std::remove, std::FILE
#include <cassert> // assert
#include <fstream> // std::ifstream
#include <sstream> // std::istringstream
#include <string> // std::string, std::getline

#include <complex> // std::complex<T>
#include "complex_tools.hxx" // is_complex, to_complex_t
#include "recorded_warnings.hxx" // warn
#include "status.hxx" // status_t

namespace debug_tools {

  template <typename real_t>
  int read_from_file(
        real_t y_data[] // 2D-data array [N*Stride]
      , char const *filename
      , size_t const N // expected number of rows
      , size_t const Stride=1 // Stride, default=1 for 1D functions
      , size_t const M=1 // number of columns to fill
      , char const *title=nullptr // title (only for display to log)
      , int const echo=0 // report writing to stdout, 0: suppress report
  ) {
      assert(Stride >= M);
      std::ifstream infile(filename, std::ifstream::in);
      if (infile.fail()) {
          if (echo > 1) printf("# %s Error opening file %s!\n", __func__, filename);
          return 1;
      } // failed to open
      
      char const rc_name[][8] = {"real", "complex"};
      bool const read_complex = is_complex<real_t>();

      std::string line; // read ASCII file line-by-line
      unsigned linenumber = 1;
      int ix = 0;
      while (std::getline(infile, line)) {
          char const *const line_ptr = line.c_str();
          char const c0 = line_ptr ? line_ptr[0] : '\0';
          switch (c0) {
              case '#': {
                  if (echo > 3) printf("# %s:%d reads \'%s\'.\n", filename, linenumber, line.c_str());
                  char const c1 = line_ptr[1]; // this char should be 'c' if the file has complex numbers
                  bool const contains_complex = ('c' == c1);
                  if (contains_complex != read_complex) {
                      warn("file %s contains %s numbers but trying to read %s numbers",
                           filename, rc_name[contains_complex], rc_name[read_complex]);
                  } // mismatch!
              } break;
              case ' ': case '\n': case '\t': case '\0': {
                  if (echo > 9) printf("# %s:%d reads \'%s\'.\n", filename, linenumber, line.c_str());
              } break;
              default: {
                  std::istringstream iss(line);
                  double x{0};
                  iss >> x;
                  for (int iy = 0; iy < M; ++iy) {
                      double y{0};
                      iss >> y;
                      if (ix < N) {
                          if (read_complex) {
                              double y_imag{0};
                              iss >> y_imag; // read also an imaginary part
                              using real_real_t = decltype(std::real(real_t(1)));
                              auto const yc = std::complex<real_real_t>(y, y_imag);
                              y_data[ix*Stride + iy] = to_complex_t<real_t, real_real_t>(yc);
                          } else {
                              y_data[ix*Stride + iy] = y;
                          } // read_complex
                      } // ix < N
                  } // iy
                  ++ix;
              } // default
          } // switch
          ++linenumber;
      } // while
      if (echo > 3) printf("# %d (of %ld) x %ld (of %ld) data entries%s%s read from file %s\n", 
                              ix, N, M, Stride, title?" for ":"", title, filename);
      return ix - N; // 0 if number of lines + 1 matched the expected number of rows
  } // read_from_file


  template <char const write_read_delete='w'>
  status_t manage_stop_file(int & number, int const echo=0, char const *filename="running.a43") {
      if ('w' == (write_read_delete | 32)) {
          std::FILE *f = std::fopen(filename, "w");
          if (nullptr == f) {
              if (echo > 1) printf("# %s<%c> unable to open file \"%s\"for writing!\n", __func__, write_read_delete, filename);
              return 1;
          } // failed to open
          std::fprintf(f, "%d   is the max. number of self-consistency iterations, user may modify", number);
          return std::fclose(f);
      } else
      if ('r' == (write_read_delete | 32)) {
          std::ifstream infile(filename, std::ifstream::in);
          if (infile.fail()) {
              if (echo > 0) printf("# %s<%c> unable to find file \"%s\" for reading\n", __func__, write_read_delete, filename);
              return 1; // error
          }
          infile >> number;
          if (echo > 1) printf("# %s<%c> found number=%d in file \"%s\" \n", __func__, write_read_delete, number, filename);
          return 0;
      } else
      if ('d' == (write_read_delete | 32)) {
          if (echo > 0) printf("# %s<%c> removes file \"%s\"\n", __func__, write_read_delete, filename);
          return std::remove(filename);
      } else {
          if (echo > 0) printf("# %s<%c> only 'w'/'W'/write, 'r'/'R'/read or 'd'/'D'/delete implemented!\n", __func__, write_read_delete);
          return -1; // error
      }
  } // manage_stop_file

} // namespace debug_tools
