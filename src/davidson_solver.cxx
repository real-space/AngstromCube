#include <cstdio> // printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <algorithm> // std::swap<T>, std::sort
#include <cmath> // std::sin

#include "davidson_solver.hxx"

#include "real_space.hxx" // ::grid_t
#include "grid_operators.hxx" // ::grid_operator_t
#include "data_view.hxx" // view2D<T>, view3D<T>
#include "linear_algebra.hxx" // ::generalized_eigval
#include "control.hxx" // ::get
#include "inline_math.hxx" // set
#include "simple_math.hxx" // ::random<T>
#include "complex_tools.hxx" // conjugate

// #define FULL_DEBUG
#define DEBUG

#ifdef  DEBUG
    #include "debug_output.hxx" // dump_to_file
#endif

#ifdef FULL_DEBUG
    #define full_debug(print) print
#else
    #define full_debug(print)
#endif

#ifdef DEBUG
    #define debug(print) print
#else
    #define debug(print)
#endif


namespace davidson_solver {
  // solve iteratively for the lowest eigenstates of an implicitly given Hamiltonian using the Davidson method

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
                  tmp += conjugate(bra_ptr[dof]) * ket_ptr[dof];
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
  void show_matrix(real_t const mat[], int const stride, int const n, int const m, char const *name=nullptr, real_t const factor=1) {
      if (n < 1) return;
      if (1 == n) {
          printf("# Vector=%s ", name);
      } else {
          printf("\n# %dx%d Matrix=%s\n", n, m, name);
      }
      for(int i = 0; i < n; ++i) {
          if (n > 1) printf("#%4i ", i);
          for(int j = 0; j < m; ++j) {
              printf((1 == n)?" %.3f":" %7.3f", mat[i*stride + j]*factor);
          }   printf("\n");
      }   printf("\n");
  } // show_matrix


  template<typename complex_t, typename real_fd_t>
  status_t eigensolve(complex_t waves[] // on entry start wave functions, on exit improved eigenfunctions
    , double *const energies // export eigenenergies
    , int const nbands // number of bands
    , grid_operators::grid_operator_t<complex_t,real_fd_t> const &op
    , int const echo // =0 log output level
    , float const mbasis // =2
    , unsigned const niterations // =2
  ) {
//    float const mbasis = 2; // control::get("davidson.basis.multiplier", 2.0);
      status_t stat(0);
      if (nbands < 1) return stat;
      
      double const dV = op.get_volume_element();
      size_t const ndof = op.get_degrees_of_freedom(); // bad naming since e.g. complex numbers bear 2 DoF in it

      int const max_space = int(std::ceil(mbasis*nbands));
      int       sub_space = nbands; // init with the waves from the input
      if (echo > 0) printf("# start Davidson with %d bands, subspace size up to %d bands\n", sub_space, max_space);

      double const threshold2 = 1e-8; // do not add residual vectors with a norm2 smaller than this

      int const op_echo = echo - 16;
      

      view3D<double> matrices(2, max_space, max_space);
      auto Hmt = matrices[0], Ovl = matrices[1];
      std::vector<double> eigval(max_space);
      std::vector<double> residual_norm2s(max_space);
      
      view2D<complex_t>  psi(max_space, ndof, 0.0); //    |psi>
      view2D<complex_t> hpsi(max_space, ndof, 0.0); // ^H*|psi>
      view2D<complex_t> spsi(max_space, ndof, 0.0); // ^S*|psi>
      view2D<complex_t> epsi(max_space, ndof, 0.0); // new eigenfunctions
      
      set(psi.data(), nbands*ndof, waves); // copy nbands initial wave functions into psi

      unsigned niter = niterations;
      for(unsigned iteration = 0; iteration < niter; ++iteration) {
          if (echo > 0) printf("# Davidson iteration %i\n", iteration);

          // apply Hamiltonian and Overlap operator
          for(int istate = 0; istate < sub_space; ++istate) {
              stat += op.Hamiltonian(hpsi[istate], psi[istate], op_echo);
              stat += op.Overlapping(spsi[istate], psi[istate], op_echo);
          } // istate

          // compute matrix representation in the sub_space
          inner_products(Ovl.data(), Ovl.stride(), ndof, psi.data(), sub_space, spsi.data(), sub_space, dV);
          inner_products(Hmt.data(), Hmt.stride(), ndof, psi.data(), sub_space, hpsi.data(), sub_space, dV);

          if (echo > 9) show_matrix(Ovl.data(), Ovl.stride(), sub_space, sub_space, "Overlap");
          if (echo > 8) show_matrix(Hmt.data(), Hmt.stride(), sub_space, sub_space, "Hamiltonian");

          auto const info = linear_algebra::eigenvalues(eigval.data(), sub_space, Hmt.data(), Hmt.stride(), Ovl.data(), Ovl.stride());
          if (info) {
              warn("generalized eigenvalue problem returned INFO=%i", info);
              stat += info;
          } else {
              auto const & eigvec = Hmt;
              if (echo > 8) show_matrix(eigval.data(), 1, 1, sub_space, "Eigenvalues");
              // if (echo > 8) show_matrix(eigvec.data(), eigvec.stride(), sub_space, sub_space, "Eigenvectors");

              // now rotate the basis into the eigenspace, ToDo: we should use DGEMM-style operations
              for(int i = 0; i < sub_space; ++i) {
                  set(epsi[i], ndof, complex_t(0));
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
                      } //
                  } // i
                  int const add_bands = new_bands;
                  if (echo > 0) printf("# Davidson: add %d residual vectors with norm2 above %.3e\n", add_bands, thres2);

                  std::vector<short> indices(add_bands);
                  int new_band{0};
                  if (echo > 0) printf("# Davidson: add residual vectors: ");
                  for(int i = 0; i < sub_space; ++i) {
                      auto const rn2 = residual_norm2s[i];
                      if (rn2 >= thres2) {
                          if (new_band < max_bands) {
                              indices[new_band] = i;
                              if (echo > 0) printf(" %i", i);
                              ++new_band;
                          }
                      } //
                  } // i
                  if (echo > 0) printf("\n");
                  if (new_band != add_bands) error("new_bands=%d != %d=add_bands", new_band, add_bands);

                  for(int i = 0; i < add_bands; ++i) {
                      int const j = indices[i];
                      int const ii = sub_space + i;
                      double const f = 1./std::sqrt(residual_norm2s[j]);
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

      if (echo > 4) show_matrix(eigval.data(), 1, 1, sub_space, "Eigenvalues", eV);

      return stat;
  } // eigensolve

#ifdef  NO_UNIT_TESTS
  // maybe explicit template instantiations needed here, ToDo
  
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template <typename complex_t>
  status_t test_solver(int const echo=9) {
      status_t stat{0};
      int const nbands = std::min(8, int(control::get("davidson_solver.num.bands", 4)));
      // particle in a box: lowest mode: sin(xyz*pi/L)^3 --> k_x=k_y=k_z=pi/L
      // --> ground state energy = 3*(pi/9)**2 Rydberg = 0.182 Hartree
      //                           3*(pi/8.78)**2      = 0.192 (found)
      //                           3*(pi/8)**2         = 0.231
      // first excitation energies should be 2*(pi/9)**2 + (2*pi/9)**2 = 0.384 Hartree (3-fold degenerate)
      real_space::grid_t g(8, 8, 8);
      std::vector<complex_t> psi(nbands*g.all(), 0.0);
      std::vector<double> energies(nbands, 0.0);

      int const swm = control::get("davidson_solver.start.waves", 0.);
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
              int const ixyz = (iz*g('y') + iy)*g('x') + ix;
              for(int iband = 0; iband < nbands; ++iband) {
                  psi[iband*g.all() + ixyz] = wxyz[iband]*cos_x*cos_y*cos_z; // good start wave functions
              } // iband
          }}} // ix iy iz
          if (echo > 2) printf("# %s: use cosine solutions as start vectors\n", __func__);
      } else if (1 == swm) {
          for(size_t i = 0; i < nbands*g.all(); ++i) {
              psi[i] = simple_math::random(-1.f, 1.f); // random wave functions (most probably not very good)
          } // i
          if (echo > 2) printf("# %s: use random values as start vectors\n", __func__);
      } else {
          for(int iband = 0; iband < nbands; ++iband) {
              psi[iband*g.all() + iband] = 1; // bad start wave functions
          } // iband
          if (echo > 2) printf("# %s: use as start vectors some delta functions at the boundary\n", __func__);
      }

      grid_operators::grid_operator_t<complex_t> op(g);

      int const nit = control::get("davidson_solver.max.iterations", 1.);
      for(int it = 0; it < nit; ++it) {
          stat += eigensolve(psi.data(), energies.data(), nbands, op, echo);
      } // it
      return stat;
  } // test_solver

  status_t all_tests(int const echo) {
    status_t stat{0};
    stat += test_solver<double>(echo);
    stat += test_solver<float>(echo);
    return stat;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace davidson_solver
