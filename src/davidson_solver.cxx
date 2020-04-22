#include <cstdio> // printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <algorithm> // std::swap<T>
#include <cmath> // std::sin

#include "davidson_solver.hxx"

#include "real_space.hxx" // ::grid_t
#include "grid_operators.hxx" // ::grid_operator_t
#include "data_view.hxx" // view2D<T>, view3D<T>
#include "linear_algebra.hxx" // ::generalized_eigval
#include "control.hxx" // ::get
#include "inline_math.hxx" // set
#include "simple_math.hxx" // ::random<T>

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

  template<typename real_t>
  void inner_products(real_t s[] // result <bra|ket> [nstates][nstates]
                   , int const stride // same stride assumed for both, bra and ket
                   , size_t const ndof
                   , real_t const bra[] // assumed shape [nstates][ndof]
                   , int const nstates  // number of bra states
                   , real_t const ket[] // assumed shape [nstates][ndof]
                   , int const mstates  // number of ket states
                   , real_t const factor=1) {
      assert(stride >= mstates);
      for(int istate = 0; istate < nstates; ++istate) {
          auto const bra_ptr = &bra[istate*ndof];
          for(int jstate = 0; jstate < mstates; ++jstate) {
              auto const ket_ptr = &ket[jstate*ndof];
              real_t tmp = 0;
              for(size_t dof = 0; dof < ndof; ++dof) {
                  tmp += bra_ptr[dof] * ket_ptr[dof];
              } // dof
              s[istate*stride + jstate] = tmp*factor; // init
          } // jstate
      } // istate
  } // inner_products


  template<typename real_t>
  void vector_norm2s(real_t s[] // result <ket|ket> [mstates]
                   , size_t const ndof
                   , real_t const ket[] // assumed shape [nstates][ndof]
                   , int const mstates  // number of ket states
                   , real_t const factor=1) {
      for(int jstate = 0; jstate < mstates; ++jstate) {
          auto const ket_ptr = &ket[jstate*ndof];
          real_t tmp = 0;
          for(size_t dof = 0; dof < ndof; ++dof) {
              tmp += pow2(ket_ptr[dof]);
          } // dof
          s[jstate] = tmp*factor; // init
      } // jstate
  } // vector_norm2s


  template<typename real_t>
  void show_matrix(real_t const mat, int const stride, int const n, int const m, char const *name=nullptr) {
      printf("\n# %s=%s%c", (n > 1)?"Matrix":"Vector", name, (n > 1)?'\n':' ');
      for(int i = 0; i < n; ++i) {
          if (n > 1) printf("#%4i ", i);
          for(int j = 0; j < m; ++j) {
              printf("%9.3f", mat[i*stride + j]);
          }   printf("\n");
      }   printf("\n");
  } // show_matrix

  template<typename real_t, typename real_fd_t>
  status_t eigensolve(real_t waves[] // on entry start wave functions, on exit improved eigenfunctions
    , int const nbands // number of bands
    , grid_operators::grid_operator_t<real_t,real_fd_t> const &op
    , int const echo // =0 log output level
    , int const mbasis // =2
    , unsigned const niterations // =2
  ) {
//    int const mbasis = 2; // control::get("davidson.basis.multiplier", 2);
      status_t stat = 0;
      
      auto const g = op.get_grid();
      
      size_t const ndof = size_t(g[0]) * size_t(g[1]) * size_t(g[2]);
      if (echo > 0) printf("# start Davidson (subspace size up to %i x %i bands)\n", mbasis, nbands);

      double const threshold2 = 1e-8;

      int const op_echo = echo - 6;
      
      int const max_space = mbasis*nbands;
      int       sub_space = nbands; // init with the waves from the input

      view3D<real_t> matrices(2, max_space, max_space);
      auto Hmt = matrices[0], Ovl = matrices[1];
      
      view2D<real_t>  psi(max_space, ndof, 0.0);
      set(psi.data(), nbands*ndof, waves); // copy nbands initial wave functions
      view2D<real_t> hpsi(max_space, ndof, 0.0);
      view2D<real_t> spsi(max_space, ndof, 0.0);
      view2D<real_t> epsi(max_space, ndof, 0.0); // new eigenfunctions
      std::vector<real_t> eigval(max_space);
      std::vector<real_t> residual_norm2s(nbands);
      std::vector<int> band_index(nbands);

      for(unsigned iteration = 0; iteration < niterations; ++iteration) {
          if (echo > 0) printf("# Davidson iteration %i\n", iteration);

          // apply Hamiltonian and Overlap operator
          for(int istate = 0; istate < sub_space; ++istate) {
              stat += op.Hamiltonian(hpsi[istate], psi[istate], op_echo);
              stat += op.Overlapping(spsi[istate], psi[istate], op_echo);
          } // istate

          // compute matrix representation in the sub_space
          inner_products<real_t>(Ovl.data(), Ovl.stride(), ndof, psi.data(), sub_space, spsi.data(), sub_space, g.dV());
          inner_products<real_t>(Hmt.data(), Hmt.stride(), ndof, psi.data(), sub_space, hpsi.data(), sub_space, g.dV());

          if (echo > 8) show_matrix(Ovl.data(), Ovl.stride(), sub_space, sub_space, "Overlap");
          if (echo > 8) show_matrix(Hmt.data(), Hmt.stride(), sub_space, sub_space, "Hamiltonian");

          stat += linear_algebra::generalized_eigenvalues(sub_space, Hmt.data(), Hmt.stride(), Ovl.data(), Ovl.stride(), eigval.data());
          auto const & eigvec = Hmt;
          if (echo > 8) show_matrix(eigval.data(), sub_space, 1, sub_space, "Eigenvalues");
          // if (echo > 8) show_matrix(eigvec.data(), eigvec.stride(), sub_space, sub_space, "Eigenvectors");

          // now rotate the basis into the eigenspace, ToDo: we should use DGEMM-style operations
          for(int i = 0; i < sub_space; ++i) {
              set(epsi[i], ndof, real_t(0));
              for(int j = 0; j < sub_space; ++j) {
                  add_product(epsi[i], ndof, psi[j], eigvec[i][j]);
              } // j
          } // i


//        if (iteration < 1) 
          if (0) { // only in the first iteration

              // apply Hamiltonian and Overlap operator again
              for(int istate = 0; istate < sub_space; ++istate) {
                  stat += op.Hamiltonian(hpsi[istate], epsi[istate], op_echo);
                  stat += op.Overlapping(spsi[istate], epsi[istate], op_echo);
              } // istate

              // now compute the residual vectors and add them to the basis [ToDo: if their norm is not too small]
              for(int i = 0; i < nbands; ++i) {
                  int const ii = i + nbands;
                  set(epsi[ii], ndof, hpsi[i]); add_product(epsi[ii], ndof, spsi[i], -eigval[i]);
              } // i

              int new_bands{0};
              vector_norm2s<real_t>(residual_norm2s.data(), ndof, epsi[nbands], nbands);
              {
                  for(int ires = 0; ires < nbands; ++ires) { // serial, loop-carried dependency in new_bands
                      if (residual_norm2s[ires] > threshold2) {
                          band_index[new_bands] = ires;
                          ++new_bands;
                      }
                  } // ires
              }

              // re-pack the new residuals
              for(int ires = 0; ires < new_bands; ++ires) { // serial due to RAW dependencies
                  int const itarget = nbands + ires;
                  int const isource = nbands + band_index[ires];
                  if (itarget < isource) {
                      set(epsi[itarget], ndof, epsi[isource]);
                  } else stat += (itarget > isource); // errors
              } // ires

              sub_space += new_bands; // increase the subspace size by the number of non-vanishing residual vectors

              if (new_bands < nbands) { printf("# add %i residuals\n", new_bands); exit(__LINE__); }

          } // not in the second iteration

          std::swap(psi, epsi); // pointer swap instead of deep copy

      } // iteration

      set(waves, nbands*ndof, psi.data()); // copy result wave functions back
      // ToDo: export last set of nbands eigenvalues

      if (echo > 4) show_matrix(eigval.data(), sub_space, 1, sub_space, "Eigenvalues");

      return stat;
  } // eigensolve

#ifdef  NO_UNIT_TESTS
  // maybe explicit template instanciations needed here, ToDo
  
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template <typename real_t>
  status_t test_solver(int const echo=9) {
      status_t stat{0};
      int const nbands = std::min(8, int(control::get("davidson_solver.num.bands", 4)));
      // particle in a box: lowest mode: sin(xyz*pi/L)^3 --> k_x=k_y=k_z=pi/L
      // --> ground state energy = 3*(pi/9)**2 Rydberg = 0.182 Hartree
      //                           3*(pi/8.78)**2      = 0.192 (found)
      //                           3*(pi/8)**2         = 0.231
      // first excitation energies should be 2*(pi/9)**2 + (2*pi/9)**2 = 0.384 Hartree (3-fold degenerate)
      real_space::grid_t g(8, 8, 8);
      std::vector<real_t> psi(nbands*g.all(), 0.0);

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

      grid_operators::grid_operator_t<real_t> op(g);

      int const nit = control::get("davidson_solver.max.iterations", 1.);
      for(int it = 0; it < nit; ++it) {
          stat += eigensolve(psi.data(), nbands, op, echo);
      } // it
      return stat;
  } // test_solver

  status_t all_tests(int const echo) {
    status_t stat{0};
    stat += test_solver<double>(echo);
    return stat;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace davidson_solver
