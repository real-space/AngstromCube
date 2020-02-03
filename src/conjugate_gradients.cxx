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
#include "inline_tools.hxx" // real_t_name<real_t>
#include "simple_math.hxx" // ::random<T>
#include "display_units.h" // eV, _eV

  #include <cstdarg> // 
  #include <utility> // std::forward

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

  template <class... Args>
  char const * print(char const *const fmt, Args &&... args) {
      static char line[256]; // warning, this static field makes this function NOT threadsafe!
      std::sprintf(line, fmt, std::forward<Args>(args)... );
      return line;
  } // print

namespace conjugate_gradients {
  // solve iteratively for the lowest eigenstates of an implicitly given Hamiltonian using the conjugate gradients method

  template<typename real_t, int D0=1> // D0: vectorization
  double inner_product(size_t const ndof
                   , real_t const bra[] // assumed shape [nstates/D0][ndof][D0]
                   , real_t const ket[] // assumed shape [nstates/D0][ndof][D0]
                   , double const factor=1) {
              double tmp = 0;
              for(size_t dof = 0; dof < ndof; ++dof) {
                  tmp += bra[dof*D0] * ket[dof*D0];
              } // dof
              return tmp*factor; // init
  } // inner_product

  template<typename real_t>
  void show_matrix(real_t const mat[], int const stride, int const n, int const m, char const *name=nullptr
      , char const initchar='#', char const final_newline='\n') {
      printf("%c# %s=%s%c", initchar, (n > 1)?"Matrix":"Vector", name, (n > 1)?'\n':' ');
      for(int i = 0; i < n; ++i) {
          if (n > 1) printf("#%4i ", i);
          for(int j = 0; j < m; ++j) {
              printf("%9.3f", mat[i*stride + j]);
          }   printf("\n");
      }   
      if (final_newline) printf("%c", final_newline);
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
  status_t eigensolve(real_t eigenstates[] // on entry start wave functions, on exit improved eigenfunctions
    , int const nbands // number of bands
    , real_space_grid::grid_t<D0> const & g // grid descriptor
    , int const echo=9 // log output level
    , float const threshold=1e-8f
    , double *eigenvalues=nullptr // export results
  ) {
      status_t stat = 0;
      int const maxiter = control::get("conjugate_gradients.max.iter", 32);
      int const restart_every_n_iterations = std::max(1, 4);
      typedef real_t complex_t;
      
      if (echo > 0) printf("# start CG (onto %d bands)\n", nbands);


      // prepare the Hamiltonian and Overlapping
      std::vector<double> potential(g.dim(2)*g.dim(1)*g.dim(0), 0.0); // flat effective local potential
      std::vector<atom_image::sho_atom_t> a(1); // sho_atoms
      std::vector<atom_image::atom_image_t> ai(1); // atom_images
      a[0]  = atom_image::sho_atom_t(3, 0.5, 999); // numax=3, sigma=0.5, atom_id=999
      ai[0] = atom_image::atom_image_t(g.dim(0)*g.h[0]/2, g.dim(1)*g.h[1]/2, g.dim(2)*g.h[2]/2, 999, 0); // image position at the center, index=0 maps into list of sho_atoms
      int const nprecond = control::get("conjugate_gradients.precond", 1.); // ToDo: implement preconditioner
      bool const use_precond = false; // (nprecond > 0);
      int const bc[] = {0, 0, 0}, nn[] = {8, 8, 8}, nn_precond[] = {nprecond, nprecond, nprecond};
      finite_difference::finite_difference_t<double> kinetic(g.h, bc, nn), precond(g.h, bc, nn_precond);
      for(int d = 0; d < 3; ++d) { precond.c2nd[d][1] = 1/12.; precond.c2nd[d][0] = 1/6.; } // preconditioner is a diffusion stencil
      
      int const cg0sd1 = control::get("conjugate_gradients.steepest.descent", 0.); // 1:steepest descent, 0:conjugate gradients
      bool const use_cg = (0 == cg0sd1);
      
      int const nv = g.all();
      view2D<real_t> s(eigenstates, nv); // set wrapper
      if (echo > 8) printf("# %s allocate %.3f MByte for %d adjoint wave functions\n", 
                              __func__, nbands*nv*1e-6*sizeof(real_t), nbands);
      view2D<real_t> Ss(nbands, nv, 0.0); // get memory

      int const northo[] = {2, cg0sd1, 1};

      // get memory for single vectors:
      int const nsinglevectors = 6 + (use_precond?1:0) + (use_cg?2:0);
      if (echo > 8) printf("# %s allocate %.3f MByte for %d single vectors functions\n", 
                              __func__, nsinglevectors*nv*1e-6*sizeof(real_t), nsinglevectors);
      view2D<real_t> mem(nsinglevectors, nv, 0.0);
      auto const Hs=mem[0], grd=mem[1], SPgrd=mem[2], dir=mem[3], Sdir=mem[4], Hcon=mem[5],
                Pgrd=use_precond ? mem[6] : grd,
                con =use_cg ? mem[nsinglevectors - 2] :  Pgrd,
                Scon=use_cg ? mem[nsinglevectors - 1] : SPgrd;
      if (echo > 9) printf("# CG vectors: Hs=%p, grd=%p, SPgrd=%p, dir=%p, \n#   Sdir=%p, Hcon=%p, Pgrd=%p, con=%p, Scon=%p\n",
          (void*)Hs, (void*)grd, (void*)SPgrd, (void*)dir, (void*)Sdir, (void*)Hcon, (void*)Pgrd, (void*)con, (void*)Scon);

      std::vector<double> energy(nbands, 0.0), residual(nbands, 9.9e99);

      int const num_outer_iterations = control::get("conjugate_gradients.num.outer.iterations", 1.);
      for(int outer_iteration = 0; outer_iteration < num_outer_iterations; ++outer_iteration) {
          if (echo > 4) printf("# CG start outer iteration #%i\n", outer_iteration);

          for(int ib_= 0; ib_< nbands; ++ib_) {  int const ib = ib_; // write protection to ib
              if (echo > 5) printf("# start CG for band #%i\n", ib);

              
              assert( 1 == D0 );
              // apply Hamiltonian and Overlap operator
//            stat += grid_operators::grid_Hamiltonian(Hs, s[ib], g, a, ai, kinetic, potential.data());
              stat += grid_operators::grid_Overlapping(Ss[ib], s[ib], g, a, ai);

              for(int iortho = 0; iortho < northo[0]; ++iortho) {
                  if (echo > 7 && ib > 0) printf("# orthogonalize band #%i against %d lower bands\n", ib, ib);
                  for(int jb = 0; jb < ib; ++jb) {
                      auto const a = inner_product(nv, s[ib], Ss[jb], g.dV());
                      if (echo > 8) printf("# orthogonalize band #%i against band #%i with %g\n", ib, jb, a);
                      add_product(Ss[ib], nv, Ss[jb], real_t(-a));
                      add_product( s[ib], nv,  s[jb], real_t(-a));
                  } // jb

                  auto const snorm = inner_product(nv, s[ib], Ss[ib], g.dV());
                  if (echo > 7) printf("# norm^2 of band #%i is %g\n", ib, snorm);
                  if (snorm < 1e-12) {
                      printf("\n# CG failed for band #%i\n\n", ib); // ToDo: convert to a warning!
                      return 1;
                  }
                  real_t const snrm = 1./std::sqrt(snorm);
                  scale(Ss[ib], nv, snrm);
                  scale( s[ib], nv, snrm);
              } // orthogonalize s[ib] against s[0...ib-1]
              
              double res_new{1}, res_old;
              double prev_energy{0};
              
              bool restart = true;

              int it_converged = 0;
              for(int iiter = 1; (iiter <= maxiter) && (0 == it_converged); ++iiter) {

                  if (0 == (iiter - 1) % restart_every_n_iterations) restart = true;
                  if (restart) {
                      // apply Hamiltonian and Overlap operator
                      stat += grid_operators::grid_Hamiltonian(Hs, s[ib], g, a, ai, kinetic, potential.data());
                      stat += grid_operators::grid_Overlapping(Ss[ib], s[ib], g, a, ai);

                      set(Sdir, nv, real_t(0)); // clear search direction
                      set( dir, nv, real_t(0)); // clear search direction
                    
                      restart = false;
                  } // (re)start
                
                  energy[ib] = inner_product(nv, s[ib], Hs, g.dV());
                  if (echo > 7) printf("# CG energy of band #%i is %g %s in iteration #%i\n", 
                                                              ib, energy[ib]*eV, _eV, iiter);
                  
                  // compute gradient = - ( H\psi - E*S\psi )
                  set(grd, nv, Hs, real_t(-1)); add_product(grd, nv, Ss[ib], real_t(energy[ib]));
                  if (echo > 7) printf("# norm^2 (no S) of gradient %.e\n", inner_product(nv, grd, grd, g.dV()));

                  if (use_precond) {
                      assert(0); // ToDo preconditioner
                  } else {
                      assert(Pgrd == grd); // Preconditioner=unity
                  } // use_precond

                  for(int iortho = 0; iortho < northo[1]; ++iortho) {
                      if (echo > 7) printf("# orthogonalize gradient against %d bands\n", ib + 1);
                      for(int jb = 0; jb <= ib; ++jb) {
                          auto const a = inner_product(nv, Pgrd, Ss[jb], g.dV());
                          if (echo > 8) printf("# orthogonalize gradient against band #%i with %g\n", jb, a);
                          add_product(Pgrd, nv, s[jb], real_t(-a));
                      } // jb
                  } // orthogonalize Pgrd against s[0...ib]

                  // apply Overlap operator
                  stat += grid_operators::grid_Overlapping(SPgrd, Pgrd, g, a, ai);
                  double const gnorm = inner_product(nv, SPgrd, Pgrd, g.dV());
                  if (echo > 7) printf("# norm^2 of%s gradient %.e\n", use_precond?" preconditioned":"", gnorm);

                  res_old = res_new; // pass
                  res_new = inner_product(nv, Pgrd, SPgrd, g.dV());
                  residual[ib] = std::abs(res_new); // store
                  if (echo > 6) printf("# CG band #%i energy change %g %s residual %.2e in iteration #%i\n", 
                                         ib, (energy[ib] - prev_energy)*eV, _eV, residual[ib], iiter);
                  prev_energy = energy[ib];

                  if (use_cg) {
                      // conjugate gradients method
                      real_t const f = ((std::abs(res_new) > tiny<real_t>()) ? res_old/res_new : 0);
                      if (echo > 7) printf("# CG step with old/new = %g/%g = %g\n", res_old, res_new, f);

                      scale( dir, nv, f); add_product( dir, nv,  Pgrd, real_t(1));
                      scale(Sdir, nv, f); add_product(Sdir, nv, SPgrd, real_t(1));

                      set(con,  nv,  dir); // copy
                      set(Scon, nv, Sdir); // copy
                      if (echo > 7) printf("# norm^2 of con %.e\n", inner_product(nv, con, Scon, g.dV()));
                      
                      for(int iortho = 0; iortho < northo[2]; ++iortho) {
                          if (echo > 7) printf("# orthogonalize conjugate direction against %d bands\n", ib + 1);
                          for(int jb = 0; jb <= ib; ++jb) {
                              auto const a = inner_product(nv, con, Ss[jb], g.dV());
                              if (echo > 8) printf("# orthogonalize conjugate direction against band #%i with %g\n", jb, a);
                              add_product(con,  nv,  s[jb], real_t(-a));
                              add_product(Scon, nv, Ss[jb], real_t(-a));
                          } // jb
                      } // orthogonalize con against s[0...ib]
                      
                  } else {
                      assert(Scon == SPgrd); assert(con == Pgrd); // steepest descent method
                  } // CG or SD


                  double const cnorm = inner_product(nv, con, Scon, g.dV());
                  if (echo > 7) printf("# norm^2 of conjugate direction is %.e\n", cnorm);
                  if (cnorm < 1e-12) {
                      it_converged = iiter; // converged, stop the iiter loop
                  } else { 
                      real_t const cf = 1./std::sqrt(cnorm);
                      scale( con, nv, cf);
                      if (1) {
                          scale(Scon, nv, cf); 
                      } else { // recompute it from S*con?
                          stat += grid_operators::grid_Overlapping(Scon, con, g, a, ai);
                      }
                      // apply Hamiltonian
                      stat += grid_operators::grid_Hamiltonian(Hcon, con, g, a, ai, kinetic, potential.data());
                      
                      double    const sHs = energy[ib];
                      complex_t const sHc = inner_product(nv, s[ib], Hcon, g.dV()); // must be complex_t
                      double    const cHc = inner_product(nv, con, Hcon, g.dV());
                      complex_t alpha{1};
                      double    beta{0};
                      submatrix2x2(sHs, sHc, cHc, alpha, beta);
                      
                      // update the band #ib
                      scale( s[ib], nv, alpha); add_product( s[ib], nv,  con, real_t(beta));
                      scale(Ss[ib], nv, alpha); add_product(Ss[ib], nv, Scon, real_t(beta));
                      scale(Hs    , nv, alpha); add_product(Hs    , nv, Hcon, real_t(beta));

                      if (res_new < threshold) it_converged = -iiter;
                  } // cnorm

              } // iiter
              
              if (maxiter < 1) {
                  energy[ib] = inner_product(nv, s[ib], Hs, g.dV());
              } // evaluate only the energy so the display is correct
              
              if (echo > 4) printf("# CG energy of band #%i is %.9g %s, residual %.1e %s in %d iterations\n", 
                                ib, energy[ib]*eV, _eV, residual[ib],
                                it_converged ? "converged" : "reached", 
                                it_converged ? std::abs(it_converged) : maxiter);
          } // ib
          
          if (echo > 4) show_matrix(energy.data(), nbands, 1, nbands, 
              print("Eigenvalues of outer iteration #%i", outer_iteration), '#', 0);
      } // outer iterations

      if (eigenvalues) {
          for(int ib = 0; ib < nbands; ++ib) {
              eigenvalues[ib] = energy[ib];
          } // ib
      } // export eigenvalues
      
      return stat;
  } // eigensolve

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS
  
  template <typename real_t>
  status_t test_solver(int const echo=9) {
      int const nbands = std::min(8, (int)control::get("conjugate_gradients.test.num.bands", 4));
      if (echo > 3) printf("\n# %s %s<%s> with %d bands\n", __FILE__, __func__, real_t_name<real_t>(), nbands);
      status_t stat{0};
      int constexpr D0 = 1; // 1: no vectorization
      int const dims[] = {8, 8, 8};
      // particle in a box: lowest mode: sin(xyz*pi/L)^3 --> k_x=k_y=k_z=pi/L
      // --> ground state energy = 3*(pi/9)**2 Rydberg = 0.182 Hartree
      //                           3*(pi/8.78)**2      = 0.192 (found)
      //                           3*(pi/8)**2         = 0.231
      // first excitation energies should be 2*(pi/9)**2 + (2*pi/9)**2 = 0.384 Hartree (3-fold degenerate)
      real_space_grid::grid_t<D0> g(dims);
      std::vector<real_t> psi(nbands*g.all(), 0.0);

      char const swm = *control::get("conjugate_gradients.test.start.waves", "a"); // 'a':good, 'r':random
      if ('a' == swm) { // scope: create good start wave functions
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
      } else if ('r' == swm) {
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

      int const nit = control::get("conjugate_gradients.test.max.iterations", 1.);
      for(int it = 0; it < nit; ++it) {
          stat += eigensolve(psi.data(), nbands, g, echo);
      } // it
      return stat;
  } // test_solver

  status_t all_tests(int const echo) {
    status_t stat{0};
    stat += test_solver<double>(echo);
    stat += test_solver<float>(echo); // test complation and convergence
    return stat;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace conjugate_gradients