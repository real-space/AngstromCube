#pragma once

#include <cstdio> // std::printf
#include <cassert> // assert
#include <cmath> // std::abs, std::pow, std::exp, std::sqrt
#include <algorithm> // std::max, std::fill
#include <fstream> // std::ifstream
#include <sstream> // std::istringstream
#include <numeric> // std::iota
#include <vector> // std::vector

#include "status.hxx" // status_t
#include "sho_tools.hxx" // ::n2HO, ::SHO_order_t, ::order_*
#include "recorded_warnings.hxx" // warn

namespace sho_unitary {

  template <typename real_t>
  real_t signed_sqrt(real_t const x) { return (x < 0)? -std::sqrt(-x) : std::sqrt(x); }

  template <typename real_t>
  status_t read_unitary_matrix_from_file(
        real_t **u
      , int const numax
      , int &nu_high
      , char const filename[]="sho_unitary.dat"
      , int const echo=7
  ) {
      //
      // Expected File format:
      //    Each line has the following 8 entries:
      //      nx ny nz ell emm nrn nom den
      //    where
      //      nx, ny, nz      are the three Cartesian SHO indices, >= 0
      //      ell, emm        are the spherical harmonic indices
      //                      ell >= 0, -ell <= emm <= ell
      //      nrn             is the number or radial nodes, nrn >= 0
      //      nom den         encode the value of the matrix entry of u:
      //                      u = sgn(nom)*sqrt(abs(nom)/den)
      //                      den > 0, nom may be negative but since
      //                      only non-zero entries are given, nom != 0
      //
      //  if we want to squeeze it into a data-type, this would do
      //  struct { uint8_t nx, ny, nz, _spare, ell, nrn; int16_t emm;
      //           int64_t nom; uint64_t den; }; // maybe reorder members
      //  this would support the index ranges up to numax=255
      std::ifstream infile(filename, std::ios_base::in);
      bool const file_is_nu_ordered = true;

      int n_ignored{0}; // number of ignored entries
      std::string line;
      while (std::getline(infile, line)) {
          char const c0 = line[0];
          if ('#' != c0 && ' ' != c0 && '\n' != c0 && 0 != c0) {
              std::istringstream iss(line);
              int nx{-1}, ny{-1}, nz{-1}, ell{-1}, emm{0}, nrn{-1};
              int64_t nom{0}, den{1};
              if (!(iss >> nx >> ny >> nz >> ell >> emm >> nrn >> nom >> den)) {
                  std::printf("# Failed to read integer number from \"%s\"!\n", line.c_str());
                  break;
              } // error
              assert(nx >= 0);
              assert(ny >= 0);
              assert(nz >= 0);
              assert(ell >= 0);
              assert(emm*emm <= ell*ell);
              assert(nrn >= 0);
              assert(den > 0);
              assert(std::abs(nom) <= den);

              real_t const u_entry = signed_sqrt(nom/(real_t(den)));
              if (echo > 8) std::printf("%d %d %d    %d %2d %d  %.15f\n", nx, ny, nz, ell, emm, nrn, u_entry);
              int const nzyx = sho_tools::Ezyx_index(nx, ny, nz);
              int const nlnm = sho_tools::Elnm_index(ell, nrn, emm);
              int const nu_xyz = sho_tools::get_nu(nx, ny, nz);
              int const nu_rad = sho_tools::get_nu(ell, nrn);
              if (nu_xyz != nu_rad) {
                  std::printf("# off-block entry found in file <%s>: nx=%d ny=%d nz=%d (nu=%d)  ell=%d emm=%d nrn=%d (nu=%d)\n",
                                  filename, nx, ny, nz, nu_xyz, ell, emm, nrn, nu_rad);
                  return 1; // error
              } // nu matches
              int const nu = nu_xyz;
              nu_high = std::max(nu_high, nu);
              if (nu > numax) {
                  ++n_ignored; // ignore the entry
                  if (file_is_nu_ordered) return 0; // we can return already since the entries are nu-ordered,
                          // .. so we do not need to parse past the first nu which exceeds the searched range.
              } else {
                  int const nb = sho_tools::n2HO(nu); // dimension of block
                  int const ioff = sho_tools::nSHO(nu - 1); // offset from previous blocks
                  u[nu][(nzyx - ioff)*nb + (nlnm - ioff)] = u_entry; // set the entry
              } // nu in range
          } // if ...
          // process pair (a,b)
      } // while
      if (n_ignored && (echo > 2)) std::printf("# ignored %d lines in file <%s> reading up to nu=%d\n", n_ignored, filename, numax);
      return 0;
  } // read_unitary_matrix_from_file


  template <typename real_t=double>
  class Unitary_SHO_Transform {
  public:

      Unitary_SHO_Transform(int const lmax=7, int const echo=8)
        : numax_(lmax)
      {
          u_ = new real_t*[1 + numax_]; // allocate pointers to blocks
          for (int nu = 0; nu <= numax_; ++nu) { // run serial forward
              int const nb = sho_tools::n2HO(nu); // dimension of block
              u_[nu] = new real_t[nb*nb]; // allocate square blocks
              // ToDo: fill with more than pseudo-values
              std::fill(u_[nu], u_[nu] + nb*nb, 0); // clear
          } // nu

          int highest_nu{-1};
          auto const stat = read_unitary_matrix_from_file(u_, numax_, highest_nu);
          if (stat) { // an error has occured while reading it from file
              for (int nu = 0; nu <= numax_; ++nu) { // run serial forward
                  int const nb = sho_tools::n2HO(nu); // dimension of block
                  std::fill(u_[nu], u_[nu] + nb*nb, 0); // clear
                  for (int ib = 0; ib < nb; ++ib) {
                      u_[nu][ib*nb + ib] = 1; // diagonal
                  } // ib
              } // nu
              warn("I/O failed with status=%i, Unitary_SHO_Transform was initialized as unit operator!\n", int(stat));
          } // stat
          if (highest_nu < numax_) {
              warn("file for Unitary_SHO_Transform provided elements only up to numax=%d, requested %d\n", highest_nu, numax_);
          } // warn
      } // constructor

      ~Unitary_SHO_Transform() {
          for (int nu = 0; nu <= numax_; ++nu) {
              delete[] u_[nu];
          } // nu
          delete[] u_;
      } // destructor

      real_t inline get_entry(int const nzyx, int const nlnm) const { // input must both be energy ordered indices
          int const nu = sho_tools::get_nu(nzyx);
          if (nu != sho_tools::get_nu(nlnm)) return 0;
//               if (nu > numax_) return 0; // warn and return, ToDo: warning
          if (nu > numax_) std::printf("# Assumption nu <= numax failed: <numax=%d> nzyx=%d nlnm=%d nu=%d\n", numax_, nzyx, nlnm, nu);
          assert(nu <= numax_);
          assert(nu >= 0);
          int const nb = sho_tools::n2HO(nu); // dimension of block
          int const ioff = sho_tools::nSHO(nu - 1); // offset from previous blocks
          assert((nzyx - ioff) < nb);
          assert((nlnm - ioff) < nb);
          assert((nzyx - ioff) > -1);
          assert((nlnm - ioff) > -1);
          return u_[nu][(nzyx - ioff)*nb + (nlnm - ioff)];
      } // get_entry

      template <typename real_out_t>
      status_t construct_dense_matrix(
            real_out_t matrix[]
          , int const nu_max
          , int const matrix_stride=-1
          , sho_tools::SHO_order_t const row_order=sho_tools::order_Ezyx // energy-ordered Cartesian
          , sho_tools::SHO_order_t const col_order=sho_tools::order_Elnm // energy-ordered radial
      ) const {
          status_t stat(0);
          int const nSHO = sho_tools::nSHO(nu_max);
          int const stride = (matrix_stride > 0)? matrix_stride : nSHO;

          std::vector<int16_t> row_index(nSHO), col_index(nSHO);
          stat += sho_tools::construct_index_table(row_index.data(), nu_max, row_order);
          stat += sho_tools::construct_index_table(col_index.data(), nu_max, col_order);

          for (int i = 0; i < nSHO; ++i) {
              for (int j = 0; j < nSHO; ++j) {
                  matrix[i*stride + j] = get_entry(row_index[i], col_index[j]);
              } // j col index
              for (int j = nSHO; j < stride; ++j) {
                  matrix[i*stride + j] = 0; // fill stride-gaps
              } // j col index
          } // i row index

          return stat; // success
      } // construct_dense_matrix

      template <typename real_vec_t>
      status_t transform_vector(real_vec_t out[], sho_tools::SHO_order_t const out_order
                        , real_vec_t const inp[], sho_tools::SHO_order_t const inp_order
                        , int const nu_max, int const echo=0) const {
          status_t stat(0);
          int const nSHO = sho_tools::nSHO(nu_max);

          bool const c2r = sho_tools::is_Cartesian(inp_order); // if input is Cartesian, tranform to radial and vice versa
          assert(    c2r ^ sho_tools::is_Cartesian(out_order)); // either input or output must be Cartesian
          if (sho_tools::is_emm_degenerate(inp_order) || sho_tools::is_emm_degenerate(out_order)) {
              if (echo > 0) std::printf("# %s cannot operate on emm-degenerate orders, inp:order_%s out:order_%s\n", __func__,
                  sho_tools::SHO_order2string(inp_order).c_str(), sho_tools::SHO_order2string(out_order).c_str());
              return -1;
          } // emm-degeneracy found

          if (echo > 4) std::printf("# %s: from order_%s to order_%s with numax=%i\n", __func__,
              sho_tools::SHO_order2string(inp_order).c_str(), sho_tools::SHO_order2string(out_order).c_str(), nu_max);

          std::vector<char> inp_label, e_inp_label, e_out_label, out_label;
          if (echo > 7) {
              inp_label   = std::vector<char>(nSHO*8);
              e_inp_label = std::vector<char>(nSHO*8);
              e_out_label = std::vector<char>(nSHO*8);
              out_label   = std::vector<char>(nSHO*8);
              sho_tools::construct_label_table(inp_label.data(), nu_max, inp_order);
              sho_tools::construct_label_table(out_label.data(), nu_max, out_order);
              sho_tools::construct_label_table(e_out_label.data(), nu_max, c2r ? sho_tools::order_Elnm : sho_tools::order_Ezyx);
              sho_tools::construct_label_table(e_inp_label.data(), nu_max, c2r ? sho_tools::order_Ezyx : sho_tools::order_Elnm);
          } // echo

          std::vector<real_vec_t> ti, to;
          auto e_inp = inp; // e_inp is an input vector with energy ordered coefficients
          if (!sho_tools::is_energy_ordered(inp_order)) {
              if (echo > 5) std::printf("# %s: energy-order input vector from order_%s\n",
                                    __func__, sho_tools::SHO_order2string(inp_order).c_str());
              ti.resize(nSHO);
              std::vector<int16_t> inp_index(nSHO);
              stat += sho_tools::construct_index_table(inp_index.data(), nu_max, inp_order);
              for (int i = 0; i < nSHO; ++i) {
                  ti[inp_index[i]] = inp[i]; // reorder to be energy ordered
                  if (echo > 8) std::printf("# %s: inp[%s] = %g\n", __func__, &inp_label[i*8], inp[i]);
              } // i
              e_inp = ti.data();
          } // input is not energy ordered

          auto e_out = out; // e_out is an output vector with energy ordered coefficients
          if (!sho_tools::is_energy_ordered(out_order)) {
              to.resize(nSHO);
              e_out = to.data();
          } // output is not energy ordered

          // transform
          stat += (nu_max > numax_); // error, not all vector elements could be transformed
          for (int nu = 0; nu <= std::min(nu_max, numax_); ++nu) {
              int const offset = sho_tools::nSHO(nu - 1);
              int const nb = sho_tools::n2HO(nu);
              int const ib = c2r ? 1 : nb;
              int const jb = c2r ? nb : 1;
              for (int i = 0; i < nb; ++i) {
                  real_vec_t tmp = 0;
                  for (int j = 0; j < nb; ++j) {
                      int const ij = i*ib + j*jb;
                      assert( 0 <= ij );
                      assert( ij < nb*nb );
//                    if (echo > 9) std::printf("# %s: u[nu=%i][%i*%i + %i*%i] = %g\n", __func__, nu, i, ib, j, jb, u_[nu][ij]);
                      tmp += u_[nu][ij] * e_inp[offset + j]; // matrix-vector multiplication, block-diagonal in nu
                  } // j
                  e_out[offset + i] = tmp;
                  if (echo > 7) std::printf("# %s: nu=%i e_inp[%s] = %g \t e_out[%s] = %g\n", __func__,
                                  nu, &e_inp_label[(offset + i)*8], e_inp[offset + i],
                                      &e_out_label[(offset + i)*8], e_out[offset + i]);
              } // i
          } // nu

          if (!sho_tools::is_energy_ordered(out_order)) {
              if (echo > 5) std::printf("# %s: restore output vector order_%s\n",
                                __func__, sho_tools::SHO_order2string(out_order).c_str());
              std::vector<int16_t> out_index(nSHO);
              stat += sho_tools::construct_index_table(out_index.data(), nu_max, out_order);
              for (int i = 0; i < nSHO; ++i) {
                  out[i] = e_out[out_index[i]];
                  if (echo > 8) std::printf("# %s: out[%s] = %g\n", __func__, &out_label[i*8], out[i]);
              } // i
          } // output is not energy ordered

          return stat;
      } // transform_vector

      double test_unitarity(int const echo=7) const {
          double maxdevall{0};
          for (int nu = 0; nu <= numax_; ++nu) {
              // as the transform is block-diagonal, we can test each block for unitarity
              int const nb = sho_tools::n2HO(nu); // dimension of block
              double mxd[2][2] = {{0,0},{0,0}}; // mxd[{0:uuT 1:uTu}][{0:off 1:diag}]
              for (int ib = 0; ib < nb; ++ib) { // cartesian block index or radial block index
                  for (int jb = 0; jb < nb; ++jb) { //
                      double uuT = 0, uTu = 0;
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
              if (echo > 3) std::printf("# U<nu=%d> deviations: Radial (uuT) %g (%g on diagonal), Cartesian (uTu) %g (%g)\n",
                                                nu, mxd[0][0], mxd[0][1], mxd[1][0], mxd[1][1]);
          } // nu
          return maxdevall;
      } // test_unitarity

      inline int numax() const { return numax_; }

  private: // members
      real_t **u_; // block diagonal matrix entries
      int numax_;  // largest ell

  }; // class Unitary_SHO_Transform

  status_t all_tests(int const echo=0); // declaration only

} // namespace sho_unitary
