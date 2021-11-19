#pragma once

#include <cstdio> // std::printf
#include <cassert> // assert
#include <cmath> // std::abs, ::pow, ::exp, ::sqrt
#include <algorithm> // std::max, ::fill
#include <fstream> // std::ifstream
#include <sstream> // std::istringstream
#include <numeric> // std::iota
#include <vector> // std::vector

// #include "status.hxx" // status_t
typedef int status_t;

// #include "sho_tools.hxx" // ::n2HO, ::SHO_order_t, ::order_*
// #include "recorded_warnings.hxx" // warn
#define warn std::printf

namespace cho_unitary {

  template <typename real_t>
  real_t signed_sqrt(real_t const x) { return (x < 0)? -std::sqrt(-x) : std::sqrt(x); }

  template <typename real_t>
  status_t read_unitary_matrix_from_file(
        real_t **u
      , int const numax
      , int &nu_high
      , char const filename[]="cho_unitary.dat"
      , int const echo=7
  ) {
      //
      // Expected File format:
      //    Each line has the following 6 entries:
      //      nx ny emm nrn nom den
      //    where
      //      nx, ny          are the three Cartesian CHO indices, >= 0
      //      emm             is the circular harmonic index
      //      nrn             is the number or radial nodes, nrn >= 0
      //      nom den         encode the value of the matrix entry of u:
      //                      u = sgn(nom)*sqrt(abs(nom)/den)
      //                      den > 0, nom may be negative but since
      //                      only non-zero entries are given, nom != 0
      //
      std::ifstream infile(filename, std::ios_base::in);
      if (infile.fail()) {
          warn("file \'%s\' for Unitary_CHO_Transform cannot be opened!", filename);
          return -1; // file not existing
      } // failed
      bool const file_is_nu_ordered = true;

      int n_ignored{0}; // number of ignored entries
      std::string line;
      while (std::getline(infile, line)) {
          char const c0 = line[0];
          if ('#' != c0 && ' ' != c0 && '\n' != c0 && 0 != c0) {
              std::istringstream iss(line);
              int nx{-1}, ny{-1}, emm{0}, nrn{-1};
              int64_t nom{0}, den{1};
              if (!(iss >> nx >> ny >> emm >> nrn >> nom >> den)) {
                  std::printf("# Failed to read integer number from \"%s\"!\n", line.c_str());
                  break;
              } // error
              assert(nx >= 0);
              assert(ny >= 0);
              assert(nrn >= 0);
              assert(den > 0);
              assert(std::abs(nom) <= den);

              real_t const u_entry = signed_sqrt(nom/double(den));
              if (echo > 8) std::printf("%d %d    %2d %d  %.15f\n", nx, ny,   emm, nrn, u_entry);
              int const ell = std::abs(emm);
              int const nu_xy = nx + ny;
              int const nu_rad = ell + 2*nrn;
              if (nu_xy != nu_rad) {
                  std::printf("# off-block entry found in file <%s>: nx=%d ny=%d (nu=%d)  emm=%d nrn=%d (nu=%d)\n",
                                  filename, nx, ny, nu_xy, emm, nrn, nu_rad);
                  return 1; // error
              } // nu matches
              int const nu = nu_xy;
              nu_high = std::max(nu_high, nu);
              if (nu > numax) {
                  ++n_ignored; // ignore the entry
                  if (file_is_nu_ordered) return 0; // we can return already since the entries are nu-ordered,
                          // .. so we do not need to parse past the first nu which exceeds the searched range.
              } else {
                  int const nb = nu + 1; // dimension of block
                  int const nrad = (nu + emm)/2; // radial index inside the block
                  assert(nrad >= 0);
                  assert(nrad < nb);
                  assert(nx < nb);
                  u[nu][nx*nb + nrad] = u_entry; // set the entry
              } // nu in range
          } // if ...
          // process pair (a,b)
      } // while
      if (n_ignored && (echo > 2)) std::printf("# ignored %d lines in file <%s> reading up to nu=%d\n", n_ignored, filename, numax);
      return 0;
  } // read_unitary_matrix_from_file


  template <typename real_t=double>
  class Unitary_CHO_Transform {
  public:

      Unitary_CHO_Transform(int const lmax=7, int const echo=8)
        : numax_(lmax)
      {
          u_ = new real_t*[1 + numax_]; // allocate pointers to blocks
          for (int nu = 0; nu <= numax_; ++nu) { // run serial forward
              int const nb = nu + 1; // dimension of block
              u_[nu] = new real_t[nb*nb]; // allocate square blocks
              // ToDo: fill with more than pseudo-values
              std::fill(u_[nu], u_[nu] + nb*nb, 0); // clear
          } // nu

          int highest_nu{-1};
          auto const stat = read_unitary_matrix_from_file(u_, numax_, highest_nu);
          if (stat) { // an error has occured while reading it from file
              for (int nu = 0; nu <= numax_; ++nu) { // run serial forward
                  int const nb = nu + 1; // dimension of block
                  std::fill(u_[nu], u_[nu] + nb*nb, 0); // clear
                  for (int ib = 0; ib < nb; ++ib) {
                      u_[nu][ib*nb + ib] = 1; // diagonal
                  } // ib
              } // nu
              warn("I/O failed with status=%i, Unitary_SHO_Transform was initialized as unit operator!", int(stat));
          } // stat
          if (highest_nu < numax_) {
              warn("file for Unitary_SHO_Transform provided elements only up to numax=%d, requested %d", highest_nu, numax_);
          } // warn
      } // constructor

      ~Unitary_CHO_Transform() {
          for (int nu = 0; nu <= numax_; ++nu) {
              delete[] u_[nu];
          } // nu
          delete[] u_;
      } // destructor

      double test_unitarity(int const echo=7) const {
          double maxdevall{0};
          for (int nu = 0; nu <= numax_; ++nu) {
              // as the transform is block-diagonal, we can test each block for unitarity
              int const nb = nu + 1; // dimension of block
              double mxd[2][2] = {{0,0},{0,0}}; // mxd[{0:uuT 1:uTu}][{0:off 1:diag}]
              for (int ib = 0; ib < nb; ++ib) { // cartesian block index or radial block index
                  for (int jb = 0; jb < nb; ++jb) { //
                      double uuT{0}, uTu{0};
                      for (int kb = 0; kb < nb; ++kb) {
                          uuT += u_[nu][ib*nb + kb] * u_[nu][jb*nb + kb]; // contraction over radial index
                          uTu += u_[nu][kb*nb + ib] * u_[nu][kb*nb + jb]; // contraction over Cartesian index
                      } // contraction index for matrix matrix multiplication
                      if (echo > 8) std::printf("# nu=%d ib=%d jb=%d uuT=%g uTu=%g\n", nu, ib, jb, uuT, uTu);
                      int const diag = (ib == jb); // 0:offdiagonal, 1:diagonal
                      mxd[0][diag] = std::max(std::abs(uuT - diag), mxd[0][diag]);
                      mxd[1][diag] = std::max(std::abs(uTu - diag), mxd[1][diag]);
                      maxdevall = std::max(maxdevall, std::max(mxd[0][diag], mxd[1][diag]));
                  } // jb
              } // ib
              if (echo > 3) std::printf("# U<nu=%d> deviations: Radial uuT %.1e (%.1e on diagonal), Cartesian uTu %.1e (%.1e on diagonal)\n",
                                                nu, mxd[0][0], mxd[0][1], mxd[1][0], mxd[1][1]);
          } // nu
          return maxdevall;
      } // test_unitarity

      inline int numax() const { return numax_; }

  private: // members
      real_t **u_; // block diagonal matrix entries
      int numax_;  // largest ell

  }; // class Unitary_CHO_Transform

  status_t all_tests(int const echo=0); // declaration only

} // namespace cho_unitary
