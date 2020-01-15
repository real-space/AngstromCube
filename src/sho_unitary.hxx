#pragma once

typedef int status_t;

#include <cmath> // std::abs
#include <vector> // std::vector<T>

#include "sho_tools.hxx" // sho_tools::, SHO_order_t

namespace sho_unitary {

  template<typename real_t> 
  status_t read_unitary_matrix_from_file(real_t **u, int const numax, int &nu_high, 
                  char const filename[]="sho_unitary.dat", int const echo=7);

  status_t all_tests();
  
  template<typename real_t> // typically real_t=double
  class Unitary_SHO_Transform {
      public:
        
          Unitary_SHO_Transform(int const lmax=7, int const echo=8) {
              numax = lmax;
              u = new real_t*[1 + numax]; // allocate pointers to blocks
              for(int nu = 0; nu <= numax; ++nu) { // run serial forward
                  int const nb = sho_tools::n2HO(nu); // dimension of block
                  u[nu] = new real_t[nb*nb]; // allocate square blocks
                  // ToDo: fill with more than pseudo-values
                  std::fill(u[nu], u[nu] + nb*nb, 0); // clear
              } // nu
              
              int highest_nu = -1;
              auto const stat = read_unitary_matrix_from_file(u, numax, highest_nu);
              if (stat) { // an error has occured while reading it from file
                  for(int nu = 0; nu <= numax; ++nu) { // run serial forward
                      int const nb = sho_tools::n2HO(nu); // dimension of block
                      std::fill(u[nu], u[nu] + nb*nb, 0); // clear
                      for(int ib = 0; ib < nb; ++ib) {
                          u[nu][ib*nb + ib] = 1; // diagonal
                      } // ib
                  } // nu
                  printf("# Warning: I/O failed, Unitary_SHO_Transform was initialized as unit operator!\n");
              } // stat
              if (highest_nu < numax) printf("# Warning: file for Unitary_SHO_Transform did not provide enough elements!\n");
          } // constructor

          ~Unitary_SHO_Transform() {
              for(int nu = 0; nu <= numax; ++nu) {
                  delete[] u[nu];
              } // nu
              delete[] u;
          } // destructor
          
          real_t inline get_entry(int const nzyx, int const nlnm) const { // input must both be energy ordered indices
              int const nu = sho_tools::get_nu(nzyx);
              if (nu != sho_tools::get_nu(nlnm)) return 0;
//               if (nu > numax) return 0; // warn and return, ToDo: warning
              if (nu > numax) printf("# Assumption nu <= numax failed: <numax=%d> nzyx=%d nlnm=%d nu=%d\n", numax, nzyx, nlnm, nu);
              assert(nu <= numax);
              assert(nu >= 0);
              int const nb = sho_tools::n2HO(nu); // dimension of block
              int const ioff = sho_tools::nSHO(nu - 1); // offset from previous blocks
              assert((nzyx - ioff) < nb);
              assert((nlnm - ioff) < nb);
              assert((nzyx - ioff) > -1);
              assert((nlnm - ioff) > -1);
              return u[nu][(nzyx - ioff)*nb + (nlnm - ioff)];
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

          double test_unitarity(int const echo=9) {
              double maxdevall = 0;
              for(int nu = 0; nu <= numax; ++nu) {
                  // as the transform is block-diagonal, we can test each block for unitarity
                  int const nb = sho_tools::n2HO(nu); // dimension of block
                  double mxd[2][2] = {{0,0},{0,0}}; // mxd[{0:uuT 1:uTu}][{0:off 1:diag}]
                  for(int ib = 0; ib < nb; ++ib) { // cartesian block index or radial block index
                      for(int jb = 0; jb < nb; ++jb) { // 
                          double uuT = 0, uTu = 0;
                          for(int kb = 0; kb < nb; ++kb) {
                              uuT += u[nu][ib*nb + kb] * u[nu][jb*nb + kb]; // contraction over radial index
                              uTu += u[nu][kb*nb + ib] * u[nu][kb*nb + jb]; // contraction over cartesian index
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

      private:
          real_t **u; // block diagonal matrix entries
          int numax; // largest ell
  }; // class Unitary_SHO_Transform
  
} // namespace sho_unitary
