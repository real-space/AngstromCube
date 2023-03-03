#include <cstdio> // std::printf
#include <cassert> // assert
#include <vector> // std::vector
#include <numeric> // std::iota
#include <cassert> // assert
#include <cmath> // std::abs, ::pow, ::exp, ::sqrt
#include <algorithm> // std::max, ::fill
#include <fstream> // std::ifstream
#include <sstream> // std::istringstream
#include <numeric> // std::iota
#include <vector> // std::vector

#include "sho_unitary.hxx" // ::Unitary_SHO_Transform

#include "sho_tools.hxx" // ...
#include "recorded_warnings.hxx" // warn

namespace sho_unitary {

  inline double signed_sqrt(double const x) { return (x < 0)? -std::sqrt(-x) : std::sqrt(x); }

  status_t read_unitary_matrix_from_file(
        double *const *const u
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
      if (infile.fail()) {
          warn("file \'%s\' for Unitary_SHO_Transform cannot be opened!", filename);
          return -1; // file not existing
      } // failed
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

              double const u_entry = signed_sqrt(nom/double(den));
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

      Unitary_SHO_Transform::Unitary_SHO_Transform(int const lmax, int const echo) // constructor
        : numax_(lmax)
      {
          u_ = new double*[1 + numax_]; // allocate pointers to blocks
          for (int nu = 0; nu <= numax_; ++nu) { // run serial forward
              int const nb = sho_tools::n2HO(nu); // dimension of block
              u_[nu] = new double[nb*nb]; // allocate square blocks
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
              warn("I/O failed with status=%i, Unitary_SHO_Transform was initialized as unit operator!", int(stat));
          } // stat
          if (highest_nu < numax_) {
              warn("file for Unitary_SHO_Transform provided elements only up to numax=%d, requested %d", highest_nu, numax_);
          } // warn
      } // constructor
    
      double Unitary_SHO_Transform::get_entry(int const nzyx, int const nlnm) const
      { // input must both be energy ordered indices
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


      double Unitary_SHO_Transform::test_unitarity(int const echo) const
      {
          double maxdevall{0};
          for (int nu = 0; nu <= numax_; ++nu) {
              // as the transform is block-diagonal, we can test each block for unitarity
              int const nb = sho_tools::n2HO(nu); // dimension of block
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

      status_t Unitary_SHO_Transform::transform_vector(
                          double       out[], sho_tools::SHO_order_t const out_order
                        , double const inp[], sho_tools::SHO_order_t const inp_order
                        , int const nu_max, int const echo) const {
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

          std::vector<double> ti, to;
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
                  double tmp = 0;
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
  
      status_t Unitary_SHO_Transform::construct_dense_matrix(
            double matrix[]
          , int const nu_max
          , int const matrix_stride // =-1
          , sho_tools::SHO_order_t const row_order // =sho_tools::order_Ezyx // energy-ordered Cartesian
          , sho_tools::SHO_order_t const col_order // =sho_tools::order_Elnm // energy-ordered radial
      ) const
      {
          status_t stat(0);
          int const nSHO = sho_tools::nSHO(nu_max);
          int const stride = (matrix_stride > 0) ? matrix_stride : nSHO;

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


#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

#if 0
  template <int numax=9>
  status_t generate_unitary_transform(int const echo) {
//       int const mx = align<2>(1 + numax);

      int constexpr lmax = numax;
      int constexpr lx = (1 + lmax);
      int xory[lmax + 1 + lmax][lx]; // xory[lmax + m][k] --> for m <  0: x^(|m|- k)*y^k if k even
                                     //                   --> for m >= 0: y^(|m|- k)*x^k if k odd
      for (int k = 0; k < (2*lmax + 1)*lx; ++k) xory[0][k] = 0; // clear

      int Pl[lx][lx]; // generate (r^2 - u^2)^l
      for (int k = 0; k < lx*lx; ++k) Pl[0][k] = 0; // clear

      int rxy2pl[lx][lx]; // generate (x^2 + y^2)^l
      for (int k = 0; k < lx*lx; ++k) rxy2pl[0][k] = 0; // clear

      { // bnc-scope
          int bnc[lx]; bnc[0] = 1; for (int l = 1; l <= lmax; ++l) bnc[l] = 0; // init binomial coefficients
          for (int m = 0; m <= lmax; ++m) {

              if (echo > 2) std::printf("# Pl l=%d  ", m);
              int j4 = 1;
              for (int k = 0; k <= m; ++k) {
                  Pl[m][k] = bnc[k]*(1 - 2*(k & 1));
                  rxy2pl[m][k] = bnc[k];
                  if (echo > 2) std::printf(" %d", Pl[m][k]);
                  int const j4mod2 = j4 % 2;
                  int const msgn = 2*j4mod2 - 1;
                  int const sgn = 1 + j4mod2 - j4;
                  xory[lmax - m*msgn][k] = bnc[k]*sgn*msgn;
                  j4 = (j4 + 1) % 4; // mimique the behaviour of complex i^k
              } // k
              if (echo > 2) std::printf("\n");

              // prepare bnc
              for (int k = m + 1; k > 0; --k) {
                  bnc[k] += bnc[k - 1];
    //               if (echo > 8) std::printf("%6d", bnc[k]);
              } // k
    //           if (echo > 8) std::printf("%6d\n", bnc[0]);
          } // m
      } // bnc-scope

      if (echo > 2) std::printf("\n# (x^2 + y^2)^%d has highest coefficient %d\n", lmax, rxy2pl[lmax][lmax/2]);

      if (echo > 1) {
          std::printf("\n\n");
          for (int m = -lmax; m <= lmax; ++m) {
              if (echo > 2) std::printf("# m=%3d   ", m);
              for (int k = 0; k <= std::abs(m); ++k) {
                  auto const xy = xory[lmax + m][k];
                  if (echo > 2) std::printf("%6d", xy);
              } // k
              if (echo > 2) std::printf("\n");
          } // m
          if (echo > 2) std::printf("\n\n");
      } // echo

      int64_t Plm[lx][lx][2*lx]; // data-type needs to capture factorial(2*lmax)
      for (int k = 0; k < lx*lx*2*lx; ++k) Plm[0][0][k] = 0; // clear

      for (int l = 0; l <= lmax; ++l) {
          int64_t poly[2*lx + 1]; // polynomial in u
          for (int k = 0; k <= 2*l; ++k) poly[k] = 0; // clear
          for (int k = 0; k <= l;   ++k) poly[2*k] = Pl[l][k]; // copy underived

          for (int m = -l; m <= l; ++m) { // start at -l to derive l times for pure Legendre polynomials
              if (m >= 0) {          // before valid m range starts for associated Legendre Polynomials
                  if (echo > 2) std::printf("# Plm l=%d m=%d  ", l, m);
                  for (int k = 0; k <= l - m; ++k) {
                      Plm[l][m][k] = poly[k]; // store associated Legendre polynomial
                      if (echo > 2) std::printf(" %lldx^%d ", Plm[l][m][k], k);
                  } // k
                  if (echo > 2) std::printf("\n");
              }
              // derive (r^2 - u^2)^l w.r.t. u one time, i.e. for l+m times
              for (int k = 1; k <= 2*l; ++k) poly[k - 1] = k*poly[k]; poly[2*l] = 0; // derive in-place
          } // m
          for (int k = 0; k <= 2*l; ++k) assert(0 == poly[k]); // all coefficients of poly must vanish since we derive u^{2l} for 2l times
      } // l

      return 0;
  } // generate_unitary_transform

  status_t test_generation(int const echo) {
      return generate_unitary_transform(echo);
  } // test_generation
#else  // 1
  status_t test_generation(int const echo) {
      return 0; // not included
  } // test_generation
#endif // 1

  status_t test_loading(int const echo=1, int const numax=9) {
      sho_unitary::Unitary_SHO_Transform U(numax);
      auto const dev = U.test_unitarity(echo);
      if (echo > 2) std::printf("# Unitary_SHO_Transform.test_unitarity = %.1e\n", dev);
      return (dev > 2e-7); // error if deviations are too large
  } // test_loading

  status_t test_vector_transform(int const echo=9, int const numax=3) {
      sho_unitary::Unitary_SHO_Transform U(numax);
      if (echo > 3) std::printf("\n# %s %s(numax=%i, echo=%i)\n", __FILE__, __func__, numax, echo);
      int const nc = sho_tools::nSHO(numax);
      std::vector<double> vi(nc), vo(nc, 0);
      std::iota(vi.begin(), vi.end(), 0); // ascending sequence
      return U.transform_vector(vo.data(), sho_tools::order_nlm, vi.data(), sho_tools::order_zyx, numax, echo);
  } // test_vector_transform

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_generation(echo);
      stat += test_loading(echo);
      stat += test_vector_transform(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace sho_unitary

