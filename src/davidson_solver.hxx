#pragma once

#include "grid_operators.hxx" // ::grid_operator_t

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>, view3D<T>
#include "linear_algebra.hxx" // ::eigenvalues
#include "inline_math.hxx" // set, pow2
#include "complex_tools.hxx" // conjugate, is_complex, to_double_complex_t
#include "display_units.h" // eV, _eV


namespace davidson_solver {


  template<typename doublecomplex_t, typename complex_t>
  void inner_products(doublecomplex_t s[] // result <bra|ket> [nstates][mstates]
                   , int const stride // same stride assumed for both, bra and ket
                   , size_t const ndof
                   , complex_t const bra[] // assumed shape [nstates][ndof]
                   , int const nstates  // number of bra states
                   , complex_t const ket[] // assumed shape [mstates][ndof]
                   , int const mstates  // number of ket states
                   , double const factor=1) {
      assert(stride >= mstates);
      for(int istate = 0; istate < nstates; ++istate) {
          auto const bra_ptr = &bra[istate*ndof];
          for(int jstate = 0; jstate < mstates; ++jstate) {
              auto const ket_ptr = &ket[jstate*ndof];
              doublecomplex_t tmp(0);
              for(size_t dof = 0; dof < ndof; ++dof) {
                  tmp += doublecomplex_t(conjugate(bra_ptr[dof])) * doublecomplex_t(ket_ptr[dof]);
              } // dof
              s[istate*stride + jstate] = tmp*factor; // init
          } // jstate
      } // istate
  } // inner_products


  template<typename complex_t>
  void vector_norm2s(double s[] // result <ket|ket> [mstates]
                   , size_t const ndof
                   , complex_t const ket[] // assumed shape [nstates][ndof]
                   , int const mstates  // number of ket states
                   , complex_t const *bra=nullptr
                   , double const factor=1) {
      for(int jstate = 0; jstate < mstates; ++jstate) {
          auto const ket_ptr = &ket[jstate*ndof];
          auto const bra_ptr = bra ? &bra[jstate*ndof] : ket_ptr;
          double tmp{0};
          for(size_t dof = 0; dof < ndof; ++dof) {
              tmp += std::real(conjugate(bra_ptr[dof]) * ket_ptr[dof]);
          } // dof
          s[jstate] = tmp*factor; // init
      } // jstate
  } // vector_norm2s
  
  
  template<typename real_t>
  void show_matrix(real_t const mat[], int const stride, int const n, int const m, char const *name=nullptr, double const unit=1, char const *_unit="1") {
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
          }   printf("\n");
      }   printf("\n");
  } // show_matrix

  template<class operator_t>
  status_t eigensolve(typename operator_t::complex_t waves[] // on entry start wave functions, on exit improved eigenfunctions
    , double *const energies // export eigenenergies
    , int const nbands // number of bands
    , operator_t const &op
    , int const echo=0 // log output level
    , float const mbasis=2
    , unsigned const niterations=2
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

      double const threshold2 = pow2(1e-4); // do not add residual vectors with a norm < 1e-4

      int const op_echo = echo - 16; // lower log level in operator calls

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

          auto const info = linear_algebra::eigenvalues(eigval.data(), sub_space, Hmt.data(), Hmt.stride(), Ovl.data(), Ovl.stride());
          if (info) {
              warn("generalized eigenvalue problem returned INFO=%i", info);
              stat += info;
          } else {
              auto const & eigvec = Hmt;
              if (echo > 8) show_matrix(eigval.data(), 1, 1, sub_space, "Eigenvalues", eV, _eV);
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

                  std::vector<double> res_norm2_sort(sub_space);              
                  set(res_norm2_sort.data(), sub_space, residual_norm2s.data()); // copy
                  std::sort(res_norm2_sort.begin(), res_norm2_sort.end()); // sort in place

                  // now select the threshold epsi with the largest residual and add them to the basis
                  int const max_bands = max_space - sub_space; // not more than this many
                  int new_bands{0};
                  double thres2{9e99};
                  for(int i = 0; i < sub_space; ++i) {
                      auto const rn2 = res_norm2_sort[sub_space - 1 - i]; // reverse as we seek for the largest residuals
                      if (rn2 > threshold2) {
                          if (new_bands < max_bands) {
                              ++new_bands;
                              thres2 = std::min(thres2, rn2);
                          }
                      } // rn2
                  } // i
                  int const add_bands = new_bands;
                  if (echo > 0) {
                      printf("# Davidson: in iteration #%i add %d residual vectors with norm2 above %.3e\n",
                                          iteration, add_bands, thres2);
                  } // echo

                  std::vector<short> indices(add_bands);
                  int new_band{0};
                  if (echo > 9) printf("# Davidson: add residual vectors: ");
                  for(int i = 0; i < sub_space; ++i) {
                      auto const rn2 = residual_norm2s[i];
                      if (rn2 >= thres2) {
                          if (new_band < max_bands) {
                              indices[new_band] = i;
                              if (echo > 9) printf(" %i", i);
                              ++new_band;
                          }
                      } // rn2
                  } // i
                  if (echo > 9) printf("\n");
                  if (new_band != add_bands) error("new_bands=%d != %d=add_bands", new_band, add_bands);

                  for(int i = 0; i < add_bands; ++i) {
                      int const j = indices[i];
                      int const ii = sub_space + i;
                      real_t const f = 1./std::sqrt(residual_norm2s[j]);
                      set(psi[ii], ndof, epsi[j], f); // and normalize
                  } // i

                  sub_space += add_bands; // increase the subspace size by the number of non-vanishing residual vectors

                  // ToDo: make sure that we run at least one more iteration at the increased function space size
                  niter = iteration + 2;
              } else { // if the space can grow
                  niter = iteration; // stop
              }

          } // info

      } // iteration

      set(waves, nbands*ndof, psi.data());  // copy result wave functions back
      set(energies, nbands, eigval.data()); // export last set of nbands eigenvalues

      if (echo > 4) show_matrix(eigval.data(), 1, 1, sub_space, "Eigenvalues", eV, _eV);

      return stat;
  } // eigensolve


  template<class operator_t> // , typename complex_t, typename doublecomplex_t>
  status_t rotate(typename operator_t::complex_t waves[] // on entry start wave functions, on exit improved eigenfunctions
    , double *const energies // export eigenenergies
    , int const nbands // number of bands
    , operator_t const &op
    , int const echo=0)   { return eigensolve(waves, energies, nbands, op, echo, 1, 1); }

  status_t all_tests(int const echo=0);

} // namespace davidson_solver
