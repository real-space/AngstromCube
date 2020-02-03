#include <cstdio> // printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <algorithm> // std::swap<T>
#include <cmath> // std::sin

#include "conjugate_gradients.hxx"

#include "real_space_grid.hxx" // ::grid_t
#include "grid_operators.hxx" // ::grid_Hamiltonian, ::grid_Overlapping
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


namespace conjugate_gradients {
  // solve iteratively for the lowest eigenstates of an implicitly given Hamiltonian using the conjugate gradients method

  template<typename real_t, int D0=1> // D0: vectorization
  void inner_products(real_t s[] // result <bra|ket> [nstates][nstates]
                   , int const stride // same stride assumed for both, bra and ket
                   , size_t const ndof
                   , real_t const bra[] // assumed shape [nstates/D0][ndof][D0]
                   , int const nstates  // number of bra states
                   , real_t const ket[] // assumed shape [nstates/D0][ndof][D0]
                   , int const mstates  // number of ket states
                   , real_t const factor=1) {
      assert(stride >= mstates);
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

  template<typename real_t, int D0=1> // D0: vectorization
  double inner_product(size_t const ndof
                   , real_t const bra[] // assumed shape [nstates/D0][ndof][D0]
                   , real_t const ket[] // assumed shape [nstates/D0][ndof][D0]
                   , real_t const factor=1) {
              double tmp = 0;
              for(size_t dof = 0; dof < ndof; ++dof) {
                  tmp += bra[dof*D0] * ket[dof*D0];
              } // dof
              return tmp*factor; // init
  } // inner_product

  template<typename real_t, int D0=1> // D0: vectorization
  void vector_norm2s(real_t s[] // result <ket|ket> [mstates]
                   , size_t const ndof
                   , real_t const ket[] // assumed shape [nstates/D0][ndof][D0]
                   , int const mstates  // number of ket states
                   , real_t const factor=1) {
      for(int jstate = 0; jstate < mstates; ++jstate) {     int const jvec = jstate/D0, j0 = jstate % D0;
          real_t tmp = 0;
          for(size_t dof = 0; dof < ndof; ++dof) {
              tmp += ket[(jvec*ndof + dof)*D0 + j0]
                   * ket[(jvec*ndof + dof)*D0 + j0];
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

  template<typename complex_t>
  double submatrix2x2(double const sAs, complex_t sAp, double const pAp, complex_t & alpha, double & beta) {
      // lowest eigenvalue of the 2 x 2 matrix ( sAs , pAs )
      //                                       ( sAp , pAp )
      double const sAppAs = pow2(std::abs(sAp)); // == |sAp|^2
      double const gamma = 0.5*(sAs + pAp) - std::sqrt(0.25*pow2(sAs - pAp) + sAppAs);

      if (std::abs((sAs - gamma)/gamma) > 2e-16) {
          double const tm2 = sAppAs + pow2(sAs - gamma);
          // ToDo: check if tm2 > 0;
          double const tmp = 1.0/std::sqrt(tm2);
          alpha = sAp*tmp;
          beta  = (gamma - sAs)*tmp;
      } else {
          alpha = 1;
          beta  = 0;
      }
      return gamma;
  } // submatrix2x2
  
  template<typename real_t, int D0> // D0: vectorization
  status_t solve(real_t waves[] // on entry start wave functions, on exit improved eigenfunctions
    , int const nbands // number of bands
    , real_space_grid::grid_t<D0> const & g // grid descriptor
    , int const echo=9 // log output level
  ) {
      status_t stat = 0;
      int const maxiter = control::get("conjugate_gradients.max.iter", 4);
      typedef double complex_t;
      
      if (echo > 0) printf("# start CG (onto %i bands)\n", nbands);

      double const threshold2 = 1e-8;

      // prepare the Hamiltonian and Overlapping
      std::vector<double> potential(g.dim(2)*g.dim(1)*g.dim(0), 0.0); // flat effective local potential
      std::vector<atom_image::sho_atom_t> a(1); // sho_atoms
      std::vector<atom_image::atom_image_t> ai(1); // atom_images
      a[0]  = atom_image::sho_atom_t(3, 0.5, 999); // numax=3, sigma=0.5, atom_id=999
      ai[0] = atom_image::atom_image_t(g.dim(0)*g.h[0]/2, g.dim(1)*g.h[1]/2, g.dim(2)*g.h[2]/2, 999, 0); // image position at the center, index=0 maps into list of sho_atoms
      int const bc[] = {0, 0, 0}, nn[] = {8, 8, 8}, nn_precond[] = {1, 1, 1};
      finite_difference::finite_difference_t<double> kinetic(g.h, bc, nn), precond(g.h, bc, nn_precond);
      for(int d = 0; d < 3; ++d) { precond.c2nd[d][1] = 1/12.; precond.c2nd[d][0] = 1/6.; } // preconditioner is a diffusion stencil
      bool constexpr use_precond = false;
      bool constexpr use_cg = true; // or steepest descent otherwise

      int const nv = g.all();
      view2D<real_t> s(waves,   nv); // set wrapper
      view2D<real_t> Ss(nbands, nv, 0.0); // get memory

      int const northo[] = {1, 0, 1};
      
      // get memory for single vectors:
      view2D<real_t> mem(9, nv, 0.0);
      auto const Hs=mem[0], grd=mem[1], SPgrd=mem[2], dir=mem[3], Sdir=mem[4], con=mem[5], Scon=mem[6], Hcon=mem[7], 
          Pgrd=use_precond ? mem[8] : grd;

      std::vector<double> energy(nbands, 0.0), residual(nbands, 9e9);
                  
      for(int iteration = 0; iteration < 1; ++iteration) {
          if (echo > 0) printf("# CG outer iteration %i\n", iteration);
          
          for(int ib_ = 0; ib_ < nbands; ++ib_) {      int const ib = ib_; // write protection to ib
              if (echo > 0) printf("# start CG for band #%i\n", ib);
              
              bool converged{false};
            
              assert( 1 == D0 );
              // apply Hamiltonian and Overlap operator
              stat += grid_operators::grid_Hamiltonian(Hs, s[ib], g, a, ai, kinetic, potential.data());
              stat += grid_operators::grid_Overlapping(Ss[ib], s[ib], g, a, ai);

              for(int iortho = 0; iortho < northo[0]; ++iortho) {
                  for(int jb = 0; jb < ib; ++jb) {
                      auto const a = inner_product(nv, s[ib], Ss[jb], g.dV());
                      if (echo > 0) printf("# orthogonalize band #%i against band #%i\n", ib, jb);
                      add_product(Ss[ib], nv, Ss[jb], -a);
                      add_product( s[ib], nv,  s[jb], -a);
                  } // jb

                  auto const snorm = inner_product(nv, s[ib], Ss[ib], g.dV());
                  if (snorm < 1e-12) {
                      printf("\n# CG failed for band #%i\n\n", ib); // ToDo: convert to a warning!
                      return 1;
                  }
                  auto const snrm = 1./std::sqrt(snorm);
                  scale(Ss[ib], nv, snrm);
                  scale( s[ib], nv, snrm);
              } // orthogonalize s[ib] against s[0...ib-1]
              
              // apply Hamiltonian and Overlap operator
              stat += grid_operators::grid_Hamiltonian(Hs, s[ib], g, a, ai, kinetic, potential.data());
              stat += grid_operators::grid_Overlapping(Ss[ib], s[ib], g, a, ai);
              
              double res_new{1}, res_old;
              set(Sdir, nv, real_t(0));
              set( dir, nv, real_t(0));
              
              int iiter{0};
              bool run = (iiter < maxiter);
              if (!run) {
                  energy[ib] = inner_product(nv, s[ib], Hs, g.dV());
              } // evaluate only the energy so the display is correct
              
              while (run) {
                  
                  res_old = res_new; // pass
                  
                  energy[ib] = inner_product(nv, s[ib], Hs, g.dV());
                  
                  // compute gradient = - ( H\psi - E*S\psi )
                  set(grd, nv, Hs, real_t(-1));
                  add_product(grd, nv, Ss[ib], energy[ib]);
                  
                  if (use_precond) {
                      assert(0); // ToDo preconditioner
                  } else assert(Pgrd == grd);
                  
                  for(int iortho = 0; iortho < northo[1]; ++iortho) {
                      for(int jb = 0; jb <= ib; ++jb) {
                          auto const a = inner_product(nv, Pgrd, Ss[jb], g.dV());
                          if (echo > 0) printf("# orthogonalize gradient against band #%i\n", jb);
                          add_product(Pgrd, nv, s[jb], -a);
                      } // jb
                  } // orthogonalize Pgrd against s[0...ib]
                      
                  // apply Overlap operator
                  stat += grid_operators::grid_Overlapping(SPgrd, Pgrd, g, a, ai);
                  
                  res_new = inner_product(nv, Pgrd, SPgrd, g.dV());
                  residual[ib] = std::abs(res_new); // store
                  
                  if (use_cg) {
                      // conjugate gradients method
                      double const f = (std::abs(res_new) > tiny<real_t>()) ? res_old/res_new : 0;
                      
                      scale(dir, nv, f);
                      add_product(dir, nv, Pgrd, real_t(1));
                      
                      scale(Sdir, nv, f);
                      add_product(dir, nv, SPgrd, real_t(1));

                      set(con,  nv,  dir);
                      set(Scon, nv, Sdir);
                  } else {
                      // steepest descent method
                      set(con,  nv,  Pgrd);
                      set(Scon, nv, SPgrd);
                  } // CG or SD
                  
                      
                  for(int iortho = 0; iortho < northo[2]; ++iortho) {
                      for(int jb = 0; jb <= ib; ++jb) {
                          auto const a = inner_product(nv, con, Ss[jb], g.dV());
                          if (echo > 0) printf("# orthogonalize conjugate direction against band #%i\n", jb);
                          add_product(con,  nv,  s[jb], -a);
                          add_product(Scon, nv, Ss[jb], -a);
                      } // jb
                  } // orthogonalize con against s[0...ib]

                  double const cnorm = inner_product(nv, con, Scon, g.dV());
                  if (cnorm < 1e-12) {
                      converged = true;
                  } else { 
                      real_t const cf = 1./std::sqrt(cnorm);
                      scale(con,  nv, cf);
                      // apply Hamiltonian
                      stat += grid_operators::grid_Hamiltonian(Hcon, con, g, a, ai, kinetic, potential.data());
                      
                      scale(Scon, nv, cf); // instead of recomputing it from S*con
                      
                      double const sHs = energy[ib];
                      complex_t const sHc = inner_product(nv, s[ib], Hcon, g.dV()); // must be complex_t
                      double    const cHc = inner_product(nv, con, Scon, g.dV());
                      complex_t alpha{1};
                      double    beta{0};
                      submatrix2x2(sHs, sHc, cHc, alpha, beta);
                      
                      // update the band
                      scale(s[ib],  nv, alpha);
                      add_product(s[ib],  nv,  con, beta);
                      scale(s[ib],  nv, alpha);
                      add_product(Ss[ib], nv, Scon, beta);
                      
                      converged = (res_new > threshold2);
                  }

                  ++iiter;
                  run = (iiter < maxiter) && (!converged);
              } // while
          
          } // ib
          
          if (echo > 1) printf("# %s outer iteration #%i\n", __func__, iteration);
          if (echo > 2) show_matrix(energy.data(), nbands, 1, nbands, "Eigenvalues");
          
      } // outer iterations

      return stat;
  } // solve

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS
  
  template <typename real_t>
  status_t test_solver(int const echo=9) {
      status_t stat{0};
      int constexpr D0 = 1; // 1: no vectorization
      int const nbands = std::min(8, (int)control::get("conjugate_gradients.num.bands", 4));
      int const dims[] = {8, 8, 8};
      // particle in a box: lowest mode: sin(xyz*pi/L)^3 --> k_x=k_y=k_z=pi/L
      // --> ground state energy = 3*(pi/9)**2 Rydberg = 0.182 Hartree
      //                           3*(pi/8.78)**2      = 0.192 (found)
      //                           3*(pi/8)**2         = 0.231
      // first excitation energies should be 2*(pi/9)**2 + (2*pi/9)**2 = 0.384 Hartree (3-fold degenerate)
      real_space_grid::grid_t<D0> g(dims);
      std::vector<real_t> psi(nbands*g.all(), 0.0);

      int const swm = control::get("conjugate_gradients.start.waves", 0.);
      if (0 == swm) { // scope: create good start wave functions
          double const k = constants::pi/8.78; // ground state wave vector
          double wxyz[8] = {1, 0,0,0, 0,0,0, 0};
          for(int iz = 0; iz < dims[2]; ++iz) { wxyz[3] = iz - 3.5; double const cos_z = std::cos(k*wxyz[3]);
          for(int iy = 0; iy < dims[1]; ++iy) { wxyz[2] = iy - 3.5; double const cos_y = std::cos(k*wxyz[2]);
          for(int ix = 0; ix < dims[0]; ++ix) { wxyz[1] = ix - 3.5; double const cos_x = std::cos(k*wxyz[1]);
              if (nbands > 4) {
                  wxyz[4] = wxyz[1]*wxyz[2]; // x*y (ell=2)
                  wxyz[5] = wxyz[2]*wxyz[3]; // y*z (ell=2)
                  wxyz[6] = wxyz[3]*wxyz[1]; // z*x (ell=2)
                  wxyz[7] = wxyz[1]*wxyz[2]*wxyz[3]; // x*y*z (ell=3)
              } // nbands > 4
              int const ixyz = (iz*dims[1] + iy)*dims[0] + ix;
              for(int iband = 0; iband < nbands; ++iband) {
                  psi[iband*g.all() + ixyz] = wxyz[iband]*cos_x*cos_y*cos_z; // good start wave functions
              } // iband
          }}} // ix iy iz
          if (echo > 2) printf("# %s: use cosine solutions as start vectors\n", __func__);
      } else if (1 == swm) {
          for(size_t i = 0; i < nbands*g.all(); ++i) {
              psi[i] = simple_math::random(-1., 1.); // random wave functions (most probably not very good)
          } // i
          if (echo > 2) printf("# %s: use random values as start vectors\n", __func__);
      } else {
          for(int iband = 0; iband < nbands; ++iband) {
              psi[iband*g.all() + iband] = 1; // bad start wave functions
          } // iband
          if (echo > 2) printf("# %s: use as start vectors some delta functions at the boundary\n", __func__);
      }

      int const nit = control::get("conjugate_gradients.max.iterations", 1.);
      for(int it = 0; it < nit; ++it) {
          stat += solve(psi.data(), nbands, g, echo);
      } // it
      return stat;
  } // test_solver

  status_t all_tests(int const echo) {
    status_t stat{0};
    stat += test_solver<double>(echo);
//     stat += test_solver<float>(echo);
    return stat;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace conjugate_gradients
