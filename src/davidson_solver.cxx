#include <cstdio> // printf
#include <cassert> // assert
#include <vector> // std::vector<T>

#include "davidson_solver.hxx"

#include "real_space_grid.hxx" // ::grid_t
#include "grid_operators.hxx" // ::grid_Hamiltonian, ::grid_Overlapping
#include "data_view.hxx" // view2D<T>, view3D<T>
#include "linear_algebra.hxx" // ::generalized_eigval
#include "control.hxx" // ::get
#include "inline_math.hxx" // set

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
  
  template<typename real_t, int D0> // D0: vectorization
  void inner_products(real_t s[] // result <bra|ket> [nstates]
                   , int const stride
                   , size_t const ndof
                   , real_t const bra[] // assumed shape [nstates/D0][ndof][D0]
                   , int const nstates  // number of bra states
                   , real_t const ket[] // assumed shape [nstates/D0][ndof][D0]
                   , int const mstates  // number of ket states
                   , real_t const factor=1) {
      for(int istate = 0; istate < nstates; ++istate) {     int const ivec = istate/D0, i0 = istate % D0;
          for(int jstate = 0; jstate < mstates; ++jstate) { int const jvec = jstate/D0, j0 = jstate % D0;
              real_t tmp = 0;
              for(size_t dof = 0; dof < ndof; ++dof) {
                  tmp += bra[(ivec*ndof + dof)*D0 + i0]
                       * ket[(jvec*ndof + dof)*D0 + j0];
              } // dof
              s[istate*stride + jstate] = tmp*factor; // init
          } // jstate
      } // istate
  } // inner_products

  
  template<typename real_t>
  void show_matrix(real_t const mat, int const stride, int const n, int const m, char const *name=nullptr) {
      printf("\n# %s=%s\n", (n > 1)?"Matrix":"Vector", name);
      for(int i = 0; i < n; ++i) {
          if (n > 1) printf("#%4i ", i); else printf("#  ");
          for(int j = 0; j < m; ++j) {
              printf("%10.3f", mat[i*stride + j]);
          }   printf("\n");
      }   printf("\n");
  } // show_matrix

  template<typename real_t, int D0> // D0: vectorization
  status_t solve(real_t waves[] // on entry start wave functions, on exit improved eigenfunctions
    , int const nbands // number of bands
    , real_space_grid::grid_t<D0> const & g // grid descriptor
    , int const echo=9 // log output level
  ) {
      status_t stat = 0;
      size_t const ndof = (size_t)g.dim(0) * (size_t)g.dim(1) * (size_t)g.dim(2);
      int const mbasis = control::get("davidson.basis.multiplier", 3);
      if (echo > 0) printf("# start Davidson (subspace size up to %i x %i bands)\n", mbasis, nbands);
      
      int const max_space = mbasis*nbands;
      int       sub_space = nbands; // init with the waves from the input
      
      view3D<real_t> matrices(2, max_space, max_space);
      auto Hmt = matrices[0], Ovl = matrices[1];

      // prepare the Hamiltonian and Overlapping
      std::vector<double> potential(g.dim(2)*g.dim(1)*g.dim(0), 0.0);
      std::vector<atom_image::sho_atom_t> a(1);
      std::vector<atom_image::atom_image_t> ai(1);
      a[0]  = atom_image::sho_atom_t(3, 0.5, 999); // numax=3, sigma=0.5, atom_id=999
      ai[0] = atom_image::atom_image_t(g.dim(0)*g.h[0]/2, g.dim(1)*g.h[1]/2, g.dim(2)*g.h[2]/2, 999, 0);
      int const bc[] = {0, 0, 0}, nn[] = {8, 8, 8};
      finite_difference::finite_difference_t<double> kinetic(g.h, bc, nn);

      std::vector<real_t>  psi(max_space*g.all(), 0.0);
      set(psi.data(), nbands*g.all(), waves); // copy nbands initial wave functions
      std::vector<real_t> hpsi(max_space*g.all(), 0.0);
      std::vector<real_t> spsi(max_space*g.all(), 0.0);
      std::vector<real_t> epsi(max_space*g.all(), 0.0); // new eigenfunctions

      int iteration = 0;
      while(sub_space <= max_space) {
          if (echo > 0) printf("# Davidson iteration %i\n", iteration);

          // apply Hamiltonian and Overlap operator
          for(int istate = 0; istate < sub_space; ++istate) {
              assert( 1 == D0 );
              stat += grid_operators::grid_Hamiltonian(hpsi.data() + istate*g.all(), psi.data() + istate*g.all(), g, a, ai, kinetic, potential.data());
              stat += grid_operators::grid_Overlapping(spsi.data() + istate*g.all(), psi.data() + istate*g.all(), g, a, ai);
          } // istate

          // compute matrix representation in the sub_space
          inner_products<real_t,D0>(Ovl.data(), Ovl.stride(), ndof, psi.data(), sub_space, spsi.data(), sub_space, g.dV());
          inner_products<real_t,D0>(Hmt.data(), Hmt.stride(), ndof, psi.data(), sub_space, hpsi.data(), sub_space, g.dV());
          
          if (echo > 8) show_matrix(Ovl.data(), Ovl.stride(), sub_space, sub_space, "Overlap");
          if (echo > 8) show_matrix(Hmt.data(), Hmt.stride(), sub_space, sub_space, "Hamiltonian");

          std::vector<real_t> eigval(sub_space);
          stat += linear_algebra::generalized_eigenvalues(sub_space, Hmt.data(), Hmt.stride(), Ovl.data(), Ovl.stride(), eigval.data());
          auto const eigvec = Hmt;
          if (echo > 8) show_matrix(eigval.data(), sub_space, 1, sub_space, "Eigenvalues");
          if (echo > 8) show_matrix(eigvec.data(), eigvec.stride(), sub_space, sub_space, "Eigenvectors");

          // now rotate the basis into the eigenspace, ToDo: can we use DGEMM-style operations?
          int const stride = g.all();
          for(int i = 0; i < sub_space; ++i) {
              for(int j = 0; j < sub_space; ++j) {
                  for(int k = 0; k < stride; ++k) {
                      epsi[i*stride + k] += eigvec[i][j] * psi[j*stride + k];
                  } // k
              } // j
          } // i


          if (sub_space < max_space) {
            
              // apply Hamiltonian and Overlap operator again
              for(int istate = 0; istate < sub_space; ++istate) {
                  assert( 1 == D0 );
                  stat += grid_operators::grid_Hamiltonian(hpsi.data() + istate*g.all(), epsi.data() + istate*g.all(), g, a, ai, kinetic, potential.data());
                  stat += grid_operators::grid_Overlapping(spsi.data() + istate*g.all(), epsi.data() + istate*g.all(), g, a, ai);
              } // istate

              // now compute the residual vectors and add them to the basis if their norm is not too small
              for(int i = 0; i < nbands; ++i) {
                  for(int k = 0; k < stride; ++k) {
                      epsi[(i + sub_space)*stride + k] = hpsi[i*stride + k] - eigval[i] * spsi[i*stride + k];
                  } // k
              } // i
              
          } // not in the last iteration
              
          set(psi.data(), sub_space*g.all(), epsi.data()); // copy sub_space wave functions

          sub_space += nbands;
          
      } // iteration
     
      set(waves, nbands*g.all(), psi.data()); // copy result wave functions back
      // ToDo: export last set of nbands eigenvalues
      
      return stat;
  } // solve

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_solver(int const echo=9) {
      int constexpr D0 = 1; // 1: no vectorization
      int const nbands = 4;
      int const dims[] = {8, 8, 8};
      real_space_grid::grid_t<D0> g(dims);
      std::vector<double> psi(nbands*g.all(), 0.0);
      for(int iband = 0; iband < nbands; ++iband) psi[iband*g.all() + iband] = 1; // very bad start wave functions
      return solve(psi.data(), nbands, g, echo);
  } // test_solver

  status_t all_tests() {
    auto status = 0;
    status += test_solver();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace davidson_solver
