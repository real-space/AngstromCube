#pragma once

#include <cmath> // std::abs
#include <vector> // std::vector<T>

#include "sho_tools.hxx" // ::SHO_order_t, ...

typedef int status_t;

namespace sho_unitary {

  template<typename real_t>
  status_t read_unitary_matrix_from_file(real_t **u, int const numax, int &nu_high,
                  char const filename[]="sho_unitary.dat", int const echo=7);

  template<typename real_t> // typically real_t=double
  class Unitary_SHO_Transform {
      public:

          Unitary_SHO_Transform(int const lmax=7, int const echo=8)
              : numax_(lmax) {
              u_ = new real_t*[1 + numax_]; // allocate pointers to blocks
              for(int nu = 0; nu <= numax_; ++nu) { // run serial forward
                  int const nb = sho_tools::n2HO(nu); // dimension of block
                  u_[nu] = new real_t[nb*nb]; // allocate square blocks
                  // ToDo: fill with more than pseudo-values
                  std::fill(u_[nu], u_[nu] + nb*nb, 0); // clear
              } // nu

              int highest_nu = -1;
              auto const stat = read_unitary_matrix_from_file(u_, numax_, highest_nu);
              if (stat) { // an error has occured while reading it from file
                  for(int nu = 0; nu <= numax_; ++nu) { // run serial forward
                      int const nb = sho_tools::n2HO(nu); // dimension of block
                      std::fill(u_[nu], u_[nu] + nb*nb, 0); // clear
                      for(int ib = 0; ib < nb; ++ib) {
                          u_[nu][ib*nb + ib] = 1; // diagonal
                      } // ib
                  } // nu
                  printf("# Warning: I/O failed, Unitary_SHO_Transform was initialized as unit operator!\n");
              } // stat
              if (highest_nu < numax_) printf("# Warning: file for Unitary_SHO_Transform did not provide enough elements!\n");
          } // constructor

          ~Unitary_SHO_Transform() {
              for(int nu = 0; nu <= numax_; ++nu) {
                  delete[] u_[nu];
              } // nu
              delete[] u_;
          } // destructor

          real_t inline get_entry(int const nzyx, int const nlnm) const { // input must both be energy ordered indices
              int const nu = sho_tools::get_nu(nzyx);
              if (nu != sho_tools::get_nu(nlnm)) return 0;
//               if (nu > numax_) return 0; // warn and return, ToDo: warning
              if (nu > numax_) printf("# Assumption nu <= numax failed: <numax=%d> nzyx=%d nlnm=%d nu=%d\n", numax_, nzyx, nlnm, nu);
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

          template<typename real_out_t>
          status_t construct_dense_matrix(real_out_t matrix[], int const nu_max, int const matrix_stride=-1
                            , sho_tools::SHO_order_t const row_order=sho_tools::order_Ezyx // energy-ordered Cartesian
                            , sho_tools::SHO_order_t const col_order=sho_tools::order_Elnm // energy-ordered radial
                    ) const {
              status_t stat = 0;
              int const nSHO = sho_tools::nSHO(nu_max);
              int const stride = (matrix_stride > 0)? matrix_stride : nSHO;

              std::vector<int16_t> row_index(nSHO), col_index(nSHO);
              stat += sho_tools::construct_index_table(row_index.data(), nu_max, row_order);
              stat += sho_tools::construct_index_table(col_index.data(), nu_max, col_order);

              for(int i = 0; i < nSHO; ++i) {
                  for(int j = 0; j < nSHO; ++j) {
                      matrix[i*stride + j] = get_entry(row_index[i], col_index[j]);
                  } // j col index
                  for(int j = nSHO; j < stride; ++j) {
                      matrix[i*stride + j] = 0; // fill stride-gaps
                  } // j col index
              } // i row index

              return stat; // success
          } // construct_dense_matrix

          template<typename real_vec_t>
          status_t transform_vector(real_vec_t out[], sho_tools::SHO_order_t const out_order
                            , real_vec_t const inp[], sho_tools::SHO_order_t const inp_order
                            , int const nu_max, int const echo=0) const {
              status_t stat = 0;
              int const nSHO = sho_tools::nSHO(nu_max);

              bool const c2r = sho_tools::is_Cartesian(inp_order); // input is Cartesian, tranform to radial
              assert(    c2r ^ sho_tools::is_Cartesian(out_order) ); // either input or output may be Cartesian
              if ( sho_tools::is_emm_degenerate(inp_order) || sho_tools::is_emm_degenerate(out_order) ) {
                  if (echo > 0) printf("# %s cannot operate on emm-degenerate orders, inp:order_%s out:order_%s\n", __func__,
                      sho_tools::SHO_order2string(inp_order).c_str(), sho_tools::SHO_order2string(out_order).c_str());
                  return -1;
              } // emm-degeneracy found

              if (echo > 1) printf("# %s: from order_%s to order_%s with numax=%i\n", __func__,
                  sho_tools::SHO_order2string(inp_order).c_str(), sho_tools::SHO_order2string(out_order).c_str(), nu_max);

              std::vector<char[8]> inp_label, e_inp_label, e_out_label, out_label;
              if (echo > 3) {
                  inp_label   = std::vector<char[8]>(nSHO);
                  e_inp_label = std::vector<char[8]>(nSHO);
                  e_out_label = std::vector<char[8]>(nSHO);
                  out_label   = std::vector<char[8]>(nSHO);
                  sho_tools::construct_label_table(inp_label.data(), nu_max, inp_order);
                  sho_tools::construct_label_table(out_label.data(), nu_max, out_order);
                  sho_tools::construct_label_table(e_out_label.data(), nu_max, c2r ? sho_tools::order_Elnm : sho_tools::order_Ezyx);
                  sho_tools::construct_label_table(e_inp_label.data(), nu_max, c2r ? sho_tools::order_Ezyx : sho_tools::order_Elnm);
              } // echo

              std::vector<real_vec_t> ti, to;
              auto e_inp = inp; // e_inp is an input vector with energy ordered coefficients
              if (!sho_tools::is_energy_ordered(inp_order)) {
                  if (echo > 3) printf("# %s: energy-order input vector from order_%s\n",
                                        __func__, sho_tools::SHO_order2string(inp_order).c_str());
                  ti.resize(nSHO);
                  std::vector<int16_t> inp_index(nSHO);
                  stat += sho_tools::construct_index_table(inp_index.data(), nu_max, inp_order);
                  for(int i = 0; i < nSHO; ++i) {
                      ti[inp_index[i]] = inp[i]; // reorder to be energy ordered
                      if (echo > 8) printf("# %s: inp[%s] = %g\n", __func__, inp_label[i], inp[i]);
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
              for(int nu = 0; nu <= std::min(nu_max, numax_); ++nu) {
                  int const offset = sho_tools::nSHO(nu - 1);
                  int const nb = sho_tools::n2HO(nu);
                  int const ib = c2r ? 1 : nb;
                  int const jb = c2r ? nb : 1;
                  for(int i = 0; i < nb; ++i) {
                      real_vec_t tmp = 0;
                      for(int j = 0; j < nb; ++j) {
                          int const ij = i*ib + j*jb;
                          assert( 0 <= ij );
                          assert( ij < nb*nb );
//                        if (echo > 9) printf("# %s: u[nu=%i][%i*%i + %i*%i] = %g\n", __func__, nu, i, ib, j, jb, u_[nu][ij]);
                          tmp += u_[nu][ij] * e_inp[offset + j]; // matrix-vector multiplication, block-diagonal in nu
                      } // j
                      e_out[offset + i] = tmp;
                      if (echo > 7) printf("# %s: nu=%i e_inp[%s] = %g \t e_out[%s] = %g\n", __func__,
                                      nu, e_inp_label[offset + i], e_inp[offset + i],
                                          e_out_label[offset + i], e_out[offset + i]);
                  } // i
              } // nu

              if (!sho_tools::is_energy_ordered(out_order)) {
                  if (echo > 3) printf("# %s: restore output vector order_%s\n",
                                    __func__, sho_tools::SHO_order2string(out_order).c_str());
                  std::vector<int16_t> out_index(nSHO);
                  stat += sho_tools::construct_index_table(out_index.data(), nu_max, out_order);
                  for(int i = 0; i < nSHO; ++i) {
                      out[i] = e_out[out_index[i]];
                      if (echo > 8) printf("# %s: out[%s] = %g\n", __func__, out_label[i], out[i]);
                  } // i
              } // output is not energy ordered

              return stat;
          } // transform_vector

          double test_unitarity(int const echo=7) const {
              double maxdevall = 0;
              for(int nu = 0; nu <= numax_; ++nu) {
                  // as the transform is block-diagonal, we can test each block for unitarity
                  int const nb = sho_tools::n2HO(nu); // dimension of block
                  double mxd[2][2] = {{0,0},{0,0}}; // mxd[{0:uuT 1:uTu}][{0:off 1:diag}]
                  for(int ib = 0; ib < nb; ++ib) { // cartesian block index or radial block index
                      for(int jb = 0; jb < nb; ++jb) { //
                          double uuT = 0, uTu = 0;
                          for(int kb = 0; kb < nb; ++kb) {
                              uuT += u_[nu][ib*nb + kb] * u_[nu][jb*nb + kb]; // contraction over radial index
                              uTu += u_[nu][kb*nb + ib] * u_[nu][kb*nb + jb]; // contraction over cartesian index
                          } // contraction index for matrix matrix multiplication
                          if (echo > 8) printf("# nu=%d ib=%d jb=%d uuT=%g uTu=%g\n", nu, ib, jb, uuT, uTu);
                          int const diag = (ib == jb); // 0:offdiagonal, 1:diagonal
                          mxd[0][diag] = std::max(std::abs(uuT - diag), mxd[0][diag]);
                          mxd[1][diag] = std::max(std::abs(uTu - diag), mxd[1][diag]);
                          maxdevall = std::max(maxdevall, std::max(mxd[0][diag], mxd[1][diag]));
                      } // jb
                  } // ib
                  if (echo > 1) printf("# U<nu=%d> radial deviations (uuT) %g (%g on diagonal), Cartesian (uTu) %g (%g)\n",
                                        nu, mxd[0][0], mxd[0][1], mxd[1][0], mxd[1][1]);
              } // nu
              return maxdevall;
          } // test_unitarity

      inline int numax() const { return numax_; }

      private:
          real_t **u_; // block diagonal matrix entries
          int numax_; // largest ell
  }; // class Unitary_SHO_Transform

  status_t all_tests(int const echo=0);

} // namespace sho_unitary
