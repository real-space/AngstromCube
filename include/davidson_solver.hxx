#pragma once

#include <cstdio> // std::printf
#include <complex> // std::real, std::norm
#include <vector> // std::vector<T>
#include <algorithm> // std::swap, std::max, std::min, std::sort
#include <cmath> // std::ceil, std::sqrt

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>, view3D<T>
#include "linear_algebra.hxx" // ::eigenvalues
#include "inline_math.hxx" // set, pow2
#include "complex_tools.hxx" // conjugate, is_complex, to_double_complex_t
#include "display_units.h" // eV, _eV
#include "print_tools.hxx" // printf_vector
#include "recorded_warnings.hxx" // warn, error

#ifndef NO_UNIT_TESTS
//  #include <cmath> // std::cos
//  #include <complex> // std::complex
//  #include "complex_tools.hxx" // complex_name
    #include "simple_math.hxx" // ::random
    #include "grid_operators.hxx" // ::grid_operator_t, ::empty_list_of_atoms
#endif

namespace davidson_solver {
  // An iterative eigensolver using the Davidson subspace method

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
          std::printf("\n# davidson_solver: inner_products: coeffs (%i,:) ", ibra);
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
          std::printf("# Vector=%s (%s)", name, _unit);
      } else {
          std::printf("\n# %dx%d Matrix=%s (%s)\n", n, m, name, _unit);
      } // n == 1
      for(int i = 0; i < n; ++i) {
          if (n > 1) std::printf("#%4i ", i);
          for(int j = 0; j < m; ++j) {
              std::printf((1 == n)?" %.3f":" %7.3f", std::real(mat[i*stride + j])*unit);
          } // j
          std::printf("\n");
      } // i
      std::printf("\n");
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

      int const max_space = std::ceil(mbasis*nbands);
      int sub_space{nbands}; // init with the waves from the input
      if (echo > 0) std::printf("# start Davidson with %d bands, subspace size up to %d bands\n", sub_space, max_space);

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
          if (echo > 9) std::printf("# Davidson iteration %i\n", iteration);

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
                      if (echo > 9 && dev > 1e-14) std::printf("# Davidson: the %d x %d overlap matrix deviates from %s by %.1e\n",
                          sub_space, sub_space, is_complex<doublecomplex_t>() ? "Hermitian" : "symmetric", dev);
                      if (dev > 1e-12) warn("the overlap matrix deviates by %.1e from symmetric/Hermitian", dev);
                  } // check if the overlap matrix is symmetric/Hermitian

                  auto const info = linear_algebra::eigenvalues(eigval.data(), sub_space, Ovl_copy.data(), Ovl_copy.stride());
                  if (1) {
                      std::printf("# Davidson: lowest eigenvalues of the %d x %d overlap matrix ", sub_space, sub_space);
                      for(int i = 0; i < std::min(9, sub_space) - 1; ++i) {
                          std::printf(" %.3g", eigval[i]);
                      } // i
                      if (sub_space > 9) std::printf(" ...");
                      std::printf(" %g", eigval[sub_space - 1]);
                      if (0 != info) std::printf(", info= %i\n", int(info));
                      std::printf("\n");
                  } // 1
                  if (eigval[0] <= 0.0) {
                      warn("overlap matrix is not positive definite, lowest eigenvalue is %g", eigval[0]);
                  } // one or more non-positive eigenvalues

                  int drop_bands{0}; while (eigval[drop_bands] < 1e-4) ++drop_bands;
                  n_drop = drop_bands;
                  if (n_drop > 0) {
                      if (echo > 0) std::printf("# Davidson: drop %d bands to stabilize the overlap\n", n_drop);
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
#ifdef DEBUG
                  std::printf("# Davidson: unsorted residual norms^2 ");
                  printf_vector(" %.1e", residual_norm2s.data(), sub_space);
#endif // DEBUG

                  std::vector<double> res_norm2_sort(sub_space);              
                  set(res_norm2_sort.data(), sub_space, residual_norm2s.data()); // copy
                  std::sort(res_norm2_sort.rbegin(), res_norm2_sort.rend()); // sort in place, reversed
                  // reverse as we seek for the largest residuals
#ifdef DEBUG
                  std::printf("# Davidson: largest residual norms^2 ");
                  printf_vector(" %.1e", res_norm2_sort.data(), sub_space);
#endif // DEBUG

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
                      std::printf("# Davidson: in iteration #%i add %d residual vectors", iteration, add_bands);
                      if (add_bands > 0) std::printf(" with norm2 above %.3e", thres2);
                      std::printf("\n");
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
                      std::printf("# Davidson: add %d residual vectors: ", add_bands);
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

  template <class operator_t>
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

  // ToDo: write this test (particle_in_box) such that 
  //       conjugate_gradients and
  //       davidson_solver can be tested, e.g. with a functor
  template <typename complex_t>
  inline status_t test_solver(int const echo=9) {
      int const nbands = std::min(8, int(control::get("davidson_solver.test.num.bands", 4)));
      if (echo > 3) std::printf("\n# %s %s<%s> with %d bands\n", __FILE__, __func__, complex_name<complex_t>(), nbands);
      status_t stat{0};
      // particle in a box: lowest mode: sin(xyz*pi/L)^3 --> k_x=k_y=k_z=pi/L
      // --> ground state energy = 3*(pi/9)**2 Rydberg = 0.182 Hartree
      //                           3*(pi/8.78)**2      = 0.192 (found)
      //                           3*(pi/8)**2         = 0.231
      // first excitation energies should be 2*(pi/9)**2 + (2*pi/9)**2 = 0.384 Hartree (3-fold degenerate)
      real_space::grid_t const g(8, 8, 8); // boundary conditions are isolated by default
      view2D<complex_t> psi(nbands, g.all(), complex_t(0)); // get wave functions
      std::vector<double> energies(nbands, 0.0); // vector for eigenenergies

      int const swm = control::get("davidson_solver.test.start.waves", 0.);
      if (0 == swm) { // scope: create good start wave functions
          double const k = constants::pi/8.78; // ground state wave vector
          double wxyz[8] = {1, 0,0,0, 0,0,0, 0};
          for(int iz = 0; iz < g('z'); ++iz) { wxyz[3] = iz - 3.5; double const cos_z = std::cos(k*wxyz[3]);
          for(int iy = 0; iy < g('y'); ++iy) { wxyz[2] = iy - 3.5; double const cos_y = std::cos(k*wxyz[2]);
          for(int ix = 0; ix < g('x'); ++ix) { wxyz[1] = ix - 3.5; double const cos_x = std::cos(k*wxyz[1]);
              if (nbands > 4) {
                  wxyz[4] = wxyz[1]*wxyz[2]; // x*y (ell=2)
                  wxyz[5] = wxyz[2]*wxyz[3]; // y*z (ell=2)
                  wxyz[6] = wxyz[3]*wxyz[1]; // z*x (ell=2)
                  wxyz[7] = wxyz[1]*wxyz[2]*wxyz[3]; // x*y*z (ell=3)
              } // nbands > 4
              int const izyx = (iz*g('y') + iy)*g('x') + ix;
              for(int iband = 0; iband < nbands; ++iband) {
                  psi(iband,izyx) = wxyz[iband]*cos_x*cos_y*cos_z; // good start wave functions
              } // iband
          }}} // ix iy iz
          if (echo > 2) std::printf("\n# %s: use cosine functions as start vectors\n", __func__);
      } else if (1 == swm) {
          for(int iband = 0; iband < nbands; ++iband) {
              for(int izyx = 0; izyx < g.all(); ++izyx) {
                  psi(iband,izyx) = simple_math::random(-1.f, 1.f); // random wave functions (most probably not very good)
              } // izyx
          } // iband
          if (echo > 2) std::printf("\n# %s: use random values as start vectors\n", __func__);
      } else {
          for(int iband = 0; iband < nbands; ++iband) {
              psi(iband,iband) = 1; // bad start wave functions
          } // iband
          if (echo > 2) std::printf("\n# %s: use as start vectors some delta functions at the boundary\n", __func__);
      } // swm (start wave method)

      // create a real-space grid Hamiltonian without atoms and with a flat local potential (at zero)
      using real_t = decltype(std::real(complex_t(1)));
      auto const loa = grid_operators::empty_list_of_atoms();
      grid_operators::grid_operator_t<complex_t,real_t> const op(g, loa); 

      int const nit = control::get("davidson_solver.test.max.iterations", 1.);
      for(int it = 0; it < nit; ++it) {
          stat += eigensolve(psi.data(), energies.data(), nbands, op, echo);
      } // it
      return stat;
  } // test_solver

  inline status_t all_tests(int const echo=0) {
      status_t stat{0};
      stat += test_solver<std::complex<double>>(echo);
      stat += test_solver<std::complex<float>> (echo);
      stat += test_solver<double>(echo);
      stat += test_solver<float> (echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace davidson_solver
