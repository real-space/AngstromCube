#pragma once

#include "grid_operators.hxx" // ::grid_operator_t

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>, view3D<T>
#include "linear_algebra.hxx" // ::eigenvalues
#include "inline_math.hxx" // set, pow2
#include "complex_tools.hxx" // conjugate, is_complex, to_double_complex_t
#include "display_units.h" // eV, _eV
#include "print_tools.hxx" // printf_vector

namespace davidson_solver {


  template <typename doublecomplex_t, typename complex_t>
  void inner_products(
        doublecomplex_t s[] // result <bra|ket> [nstates][mstates]
      , int const stride // stride for the result matrix
      , size_t const ndof // also assumed as stride for bra and ket
      , complex_t const bra[] // assumed shape [nstates][ndof]
      , int const nstates  // number of bra states
      , complex_t const ket[] // assumed shape [mstates][ndof]
      , int const mstates  // number of ket states
      , double const factor=1
  ) {
      assert(stride >= mstates);
      for(int ibra = 0; ibra < nstates; ++ibra) {
          auto const bra_ptr = &bra[ibra*ndof];
          for(int jket = 0; jket < mstates; ++jket) {
              auto const ket_ptr = &ket[jket*ndof];
              doublecomplex_t tmp(0);
              for(size_t dof = 0; dof < ndof; ++dof) {
                  tmp += doublecomplex_t(conjugate(bra_ptr[dof])) * doublecomplex_t(ket_ptr[dof]);
              } // dof
              s[ibra*stride + jket] = tmp*factor; // init
          } // jket
#ifdef NEVER
          printf("\n# davidson_solver: inner_products: coeffs (%i,:) ", ibra);
          printf_vector("%g ", &s[ibra*stride], mstates);
#endif // NEVER
      } // ibra
  } // inner_products


  template <typename complex_t>
  void vector_norm2s(
        double s[] // result <ket|ket> [mstates]
      , size_t const ndof
      , complex_t const ket[] // assumed shape [nstates][ndof]
      , int const mstates  // number of ket states
      , complex_t const *bra=nullptr
      , double const factor=1
  ) {
      for(int jket = 0; jket < mstates; ++jket) {
          auto const ket_ptr = &ket[jket*ndof];
          auto const bra_ptr = bra ? &bra[jket*ndof] : ket_ptr;
          double tmp{0};
          for(size_t dof = 0; dof < ndof; ++dof) {
              tmp += std::real(conjugate(bra_ptr[dof]) * ket_ptr[dof]);
          } // dof
          s[jket] = tmp*factor; // init
      } // jket
  } // vector_norm2s
  
  
  template <typename real_t>
  void show_matrix(
        real_t const mat[]
      , int const stride
      , int const n
      , int const m
      , char const *name=nullptr
      , double const unit=1
      , char const *_unit="1"
  ) {
      if (n < 1) return;
      if (is_complex<real_t>()) return;
      if (1 == n) {
          printf("# Vector=%s (%s)", name, _unit);
      } else {
          printf("\n# %dx%d Matrix=%s (%s)\n", n, m, name, _unit);
      } // n == 1
      for(int i = 0; i < n; ++i) {
          if (n > 1) printf("#%4i ", i);
          for(int j = 0; j < m; ++j) {
              printf((1 == n)?" %.3f":" %7.3f", std::real(mat[i*stride + j])*unit);
          } // j
          printf("\n");
      } // i
      printf("\n");
  } // show_matrix

  template <class operator_t>
  status_t eigensolve(
      typename operator_t::complex_t waves[] // on entry start wave functions, on exit improved eigenfunctions
    , double energies[] // export eigenenergies
    , int const nbands // number of bands
    , operator_t const & op
    , int const echo=0 // log output level
    , float const mbasis=2
    , unsigned const niterations=2
    , float const threshold=1e-4
  ) {
      using complex_t = typename operator_t::complex_t; // abbreviate
      using doublecomplex_t = decltype(to_double_complex_t(complex_t(1))); // double or complex<double>
      using real_t = decltype(std::real(complex_t(1)));
      status_t stat(0);
      if (nbands < 1) return stat;

      double const dV = op.get_volume_element();
      size_t const ndof = op.get_degrees_of_freedom(); // bad naming since e.g. complex numbers bear 2 DoF in it

      int const max_space = int(std::ceil(mbasis*nbands));
      int sub_space{nbands}; // init with the waves from the input
      if (echo > 0) printf("# start Davidson with %d bands, subspace size up to %d bands\n", sub_space, max_space);

      double const threshold2 = pow2(threshold); // do not add residual vectors with a norm < 1e-4

      auto const op_echo = echo - 16; // lower log level for operator calls

      view3D<doublecomplex_t> matrices(2, max_space, max_space);
      auto Hmt = matrices[0], Ovl = matrices[1]; // named sub-views
      std::vector<double> eigval(max_space);
      std::vector<double> residual_norm2s(max_space);

      complex_t const zero(0);
      view2D<complex_t>  psi(max_space, ndof, zero); //    |psi>
      view2D<complex_t> hpsi(max_space, ndof, zero); // ^H*|psi>
      view2D<complex_t> spsi(max_space, ndof, zero); // ^S*|psi>
      view2D<complex_t> epsi(max_space, ndof, zero); // new eigenfunctions

      set(psi.data(), nbands*ndof, waves); // copy nbands initial wave functions into psi

      unsigned niter{niterations};
      for(unsigned iteration = 0; iteration < niter; ++iteration) {
          if (echo > 9) printf("# Davidson iteration %i\n", iteration);

          int n_drop{0};
          do {
              // apply Hamiltonian and Overlap operator
              for(int istate = 0; istate < sub_space; ++istate) {
                  stat += op.Hamiltonian(hpsi[istate], psi[istate], op_echo);
                  stat += op.Overlapping(spsi[istate], psi[istate], op_echo);
              } // istate

              // compute matrix representation in the sub_space
              inner_products(Ovl.data(), Ovl.stride(), ndof, psi.data(), sub_space, spsi.data(), sub_space, dV);
              inner_products(Hmt.data(), Hmt.stride(), ndof, psi.data(), sub_space, hpsi.data(), sub_space, dV);

              if (echo > 9) show_matrix(Ovl.data(), Ovl.stride(), sub_space, sub_space, "Overlap");
              if (echo > 8) show_matrix(Hmt.data(), Hmt.stride(), sub_space, sub_space, "Hamiltonian", eV, _eV);

              if (0) { // inspect eigenvalues of the overlap matrix
                  view2D<doublecomplex_t> Ovl_copy(sub_space, sub_space, doublecomplex_t(0));
                  for(int i = 0; i < sub_space; ++i) {
                      set(Ovl_copy[i], sub_space, Ovl[i]); // deep copy
                  } // i

                  if (1) { // check if the overlap matrix is symmetric/Hermitian
                      double dev2{0};
                      for(int i = 0; i < sub_space; ++i) {
                          for(int j = 0; j < i; ++j) {
                              dev2 = std::max(dev2, std::norm(Ovl_copy(i,j) - conjugate(Ovl_copy(j,i))));
                          } // j
                      } // i
                      auto const dev = std::sqrt(std::max(0.0, dev2)); // since std::norm(z) generates |z|^2
                      if (echo > 9 && dev > 1e-14) printf("# Davidson: the %d x %d overlap matrix deviates from %s by %.1e\n",
                          sub_space, sub_space, is_complex<doublecomplex_t>() ? "Hermitian" : "symmetric", dev);
                      if (dev > 1e-12) warn("the overlap matrix deviates by %.1e from symmetric/Hermitian", dev);
                  } // check if the overlap matrix is symmetric/Hermitian

                  auto const info = linear_algebra::eigenvalues(eigval.data(), sub_space, Ovl_copy.data(), Ovl_copy.stride());
                  if (1) {
                      printf("# Davidson: lowest eigenvalues of the %d x %d overlap matrix ", sub_space, sub_space);
                      for(int i = 0; i < std::min(9, sub_space) - 1; ++i) {
                          printf(" %.3g", eigval[i]);
                      } // i
                      if (sub_space > 9) printf(" ...");
                      printf(" %g", eigval[sub_space - 1]);
                      if (0 != info) printf(", info= %i\n", int(info));
                      printf("\n");
                  } // 1
                  if (eigval[0] <= 0.0) {
                      warn("overlap matrix is not positive definite, lowest eigenvalue is %g", eigval[0]);
                  } // one or more non-positive eigenvalues

                  int drop_bands{0}; while (eigval[drop_bands] < 1e-4) ++drop_bands;
                  n_drop = drop_bands;
                  if (n_drop > 0) {
                      if (echo > 0) printf("# Davidson: drop %d bands to stabilize the overlap\n", n_drop);
                      for(int i = 0; i < sub_space - n_drop; ++i) {
                          int const ii = i + n_drop;
                          set(epsi[i], ndof, zero);
                          for(int j = 0; j < sub_space; ++j) {
                              add_product(epsi[i], ndof, psi[j], complex_t(Ovl_copy(ii,j)));
                          } // j
                      } // i
                      std::swap(psi, epsi); // pointer swap instead of deep copy
                      sub_space -= n_drop;
                  } // drop n bands

              } // inspect eigenvalues of the overlap matrix
          } while(n_drop > 0); //while

          auto const info = linear_algebra::eigenvalues(eigval.data(), sub_space, Hmt.data(), Hmt.stride(), Ovl.data(), Ovl.stride());
          if (info) {
              warn("generalized eigenvalue problem returned INFO=%i", info);
              stat += info;
          } else {
              auto const & eigvec = Hmt; // on success eigenvectors are stored in the memory location of the Hamiltonian
              if (echo > 8) show_matrix(eigval.data(), 0, 1, sub_space, "Eigenvalues", eV, _eV);
              // if (echo > 8) show_matrix(eigvec.data(), eigvec.stride(), sub_space, sub_space, "Eigenvectors");

              // now rotate the basis into the eigenspace, ToDo: we should use DGEMM-style operations
              for(int i = 0; i < sub_space; ++i) {
                  set(epsi[i], ndof, zero);
                  for(int j = 0; j < sub_space; ++j) {
                      add_product(epsi[i], ndof, psi[j], complex_t(eigvec(i,j)));
                  } // j
              } // i
              std::swap(psi, epsi); // pointer swap instead of deep copy

              if (sub_space < max_space) {

                  // apply Hamiltonian and Overlap operator again
                  bool const with_overlap = true;
                  for(int i = 0; i < sub_space; ++i) {
                      stat += op.Hamiltonian(hpsi[i], psi[i], op_echo);
                      set(epsi[i], ndof, hpsi[i]);
                      stat += op.Overlapping(spsi[i], psi[i], op_echo);
                      add_product(epsi[i], ndof, spsi[i], complex_t(-eigval[i]));

                      if (with_overlap) stat += op.Overlapping(spsi[i], epsi[i], op_echo);
                  } // i
                  vector_norm2s(residual_norm2s.data(), ndof, epsi.data(), sub_space, 
                                               with_overlap ? spsi.data() : nullptr, dV);
                  if (1) {
                      printf("# Davidson: unsorted residual norms^2 ");
                      printf_vector(" %.1e", residual_norm2s.data(), sub_space);
                  } // 1

                  std::vector<double> res_norm2_sort(sub_space);              
                  set(res_norm2_sort.data(), sub_space, residual_norm2s.data()); // copy
                  std::sort(res_norm2_sort.rbegin(), res_norm2_sort.rend()); // sort in place, reversed
                  // reverse as we seek for the largest residuals
                  if (1) {
                      printf("# Davidson: largest residual norms^2 ");
                      printf_vector(" %.1e", res_norm2_sort.data(), sub_space);
                  } // 1

                  // now select the threshold epsi with the largest residual and add them to the basis
                  int const max_bands = max_space - sub_space; // not more than this many
                  int new_bands{0};
                  double thres2{9.999e99};
                  for(int i = 0; i < sub_space; ++i) {
                      auto const rn2 = res_norm2_sort[i];
                      if (rn2 > threshold2) {
                          if (new_bands < max_bands) {
                              ++new_bands;
                              thres2 = std::min(thres2, rn2);
                          } // new_bands < max_bands
                      } // rn2
                  } // i
                  int const add_bands = new_bands;
                  if (echo > 0) {
                      printf("# Davidson: in iteration #%i add %d residual vectors", iteration, add_bands);
                      if (add_bands > 0) printf(" with norm2 above %.3e", thres2);
                      printf("\n");
                  } // echo

                  std::vector<short> indices(add_bands);
                  int new_band{0};
                  for(int i = 0; i < sub_space; ++i) {
                      auto const rn2 = residual_norm2s[i];
                      if (rn2 >= thres2) {
                          if (new_band < max_bands) {
                              indices[new_band] = i;
                              ++new_band;
                          } // new_band < max_bands
                      } // rn2
                  } // i

                  if (echo > 9 && add_bands > 0) {
                      printf("# Davidson: add %d residual vectors: ", add_bands);
                      printf_vector(" %i", indices.data(), new_band);
                  } // echo
                  if (new_band != add_bands) error("new_bands=%d != %d=add_bands", new_band, add_bands);

                  for(int i = 0; i < add_bands; ++i) {
                      int const j = indices[i];
                      int const ii = sub_space + i;
                      real_t const f = 1./std::sqrt(residual_norm2s[j]);
                      set(psi[ii], ndof, epsi[j], f); // and normalize
                  } // i

                  sub_space += add_bands; // increase the subspace size by the number of non-vanishing residual vectors

                  // ToDo: make sure that we run at least one more iteration at the increased function space size
                  niter = iteration + 2 * (add_bands > 0);
              } else { // if the space can grow
                  niter = iteration; // stop
              } // sub_space < max_space

          } // info

      } // iteration

      set(waves, nbands*ndof, psi.data());  // copy result wave functions back
      set(energies, nbands, eigval.data()); // export last set of nbands eigenvalues

      if (echo > 4) show_matrix(eigval.data(), 0, 1, sub_space, "Eigenvalues", eV, _eV);

      return stat;
  } // eigensolve

  template<class operator_t>
  status_t rotate(
      typename operator_t::complex_t waves[] // on entry start wave functions, on exit improved eigenfunctions
    , double energies[] // export eigenenergies
    , int const nbands // number of bands
    , operator_t const & op
    , int const echo=0
  ) {
      return eigensolve(waves, energies, nbands, op, echo, 1, 1);
  } // rotate

#ifdef NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS
  status_t all_tests(int const echo=0); // declaration only
#endif // NO_UNIT_TESTS

} // namespace davidson_solver
