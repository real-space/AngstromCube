#include <cstdio> // printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <algorithm> // std::swap<T>
#include <cmath> // std::sqrt

#include "iterative_poisson.hxx"

#include "real_space_grid.hxx" // ::grid_t
#include "data_view.hxx" // view2D<T>
#include "inline_math.hxx" // set, dot_product
#include "finite_difference.hxx" // ::finite_difference_t

#include "real_space_grid.hxx" // ::bessel_projection
#include "radial_grid.hxx" // ::radial_grid_t, ::create_equidistant_radial_grid
#include "constants.hxx" // ::pi
#include "bessel_transform.hxx" // ::transform_s_function

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

namespace iterative_poisson {
  // solve iteratively for the lowest eigenstates of an implicitly given Hamiltonian using the conjugate gradients method

  
  template<typename real_t>
  double norm2(real_t const v[], size_t const n) { double s{0}; for(size_t i{0}; i < n; ++i) { s += pow2(v[i]); } return s; }
//   { return dot_product(n, v, v); }

  template<typename real_t>
  double norm1(real_t const v[], size_t const n) { double s{0}; for(size_t i{0}; i < n; ++i) { s += v[i]; } return s; }

  template<typename real_t>
  double scalar_product(real_t const v[], real_t const w[], size_t const n) { return dot_product(n, v, w); }

  template<typename real_t>
  status_t solve(real_t x[] // result to -Laplace(x)/(4*pi) == b
                , real_t const b[] // inhomogeneous part b
                , real_space_grid::grid_t<1> g // grid descriptor
                , float const threshold=3e-8 // convergence criterion
                , double *residual=nullptr // residual that was reached
                , int const maxiter=999 // maximum number of iterations 
                , int const miniter=0   // minimum number of iterations
                , int restart=4096 // number of iterations before restrat, 1:steepest descent
                , int const echo=0 // log level
                ) {
        
    restart = std::max(1, restart);
    
    size_t const nall = size_t(g.dim(2))*size_t(g.dim(1))*size_t(g.dim(0));
    
    view2D<real_t> mem(5, nall, 0.0); // get memory
    auto const r=mem[0], p=mem[1], ax=mem[2], ap=mem[3], z=mem[4];    
    
    int const bc[] = {0, 0, 0}, nn[] = {8, 8, 8};
    finite_difference::finite_difference_t<real_t> fd(g.h, bc, nn);
    fd.scale_coefficients(-.25/constants::pi); // electrostatics prefactor
    
    double const cell_volume = nall*g.dV();
    double const threshold2 = cell_volume * pow2(threshold);

//     ! |Ax> = A|x>
//     ! |r> = |b> - |Ax>
//     ! |z> = P|r>
//     ! rz_old = <r|z>
//     ! |p> = |z>
//     ! it = 0
//     ! do while
//     !   |Ap> = A|p>
//     !   pAp = <p|Ap>
//     !   alpha = rz_old / pAp
//     !   |x> = |x> + alpha |p>
//     !   |r> = |r> - alpha |Ap>
//     !   res = <r|r>
//     !   |z> = P|r>
//     !   rz_new = <r|z>
//     !   beta = rz_new / rz_old
//     !   rz_old = rz_new
//     !   |p> = |z> + beta |p>
//     !   it = it+1

    status_t ist{0};

    // |Ax> := A|x>
    ist = finite_difference::Laplacian(ax, x, g, fd);
    if (ist) error("CG_solve: Laplacian failed!");

    dump_to_file("cg_start", nall, x, nullptr, 1, 1, "x", echo);
    
    // |r> = |b> - A|x> = |b> - |Ax>
    set(r, nall, b); add_product(r, nall, ax, real_t(-1));

    // res^2 = <r|r>
    double res2 = norm2(r, nall) * g.dV(); // g.comm
    double const res_start = std::sqrt(res2/cell_volume); // store staring residual
    if (echo > 8) printf("# %s start residual=%.1e\n", __FILE__, res_start);

    // |z> = |Pr> = P|r>
//     ist = Precon_Nf4( g, r, z )
//     if( ist /= 0 ) stop 'CG_solve: Precon_Nf4 failed!'
    set(z, nall, r); // Preconditioner == unity

    // rz_old = <r|z>
    double rz_old = scalar_product(r, z, nall) * g.dV(); // g.comm 

    // |p> = |z>
    set(p, nall, z);

    int it{0}; // init iteration counter

    // number of iterations is less then maxiter?
    bool run = (it < maxiter);
    while (run) {
      ++it;
//       !--------------------------------------
//       ! begin of the CG iteration
//       !--------------------------------------

      // |ap> = A|p>
      ist = finite_difference::Laplacian(ap, p, g, fd);
      if (ist) error("CG_solve: Laplacian failed!");

      double const pAp = scalar_product(p, ap, nall) * g.dV(); // g.comm

      double constexpr RZ_TINY = 1e-14;
      // alpha = rz_old / pAp
      double const alpha = (std::abs(pAp) < RZ_TINY) ? RZ_TINY : rz_old / pAp;

      // |x> = |x> + alpha |p>
      add_product(x, nall, p, real_t(alpha));

//       !============================================================
//       ! special treatment of completely periodic case
//       !============================================================
      if (fd.all_boundary_conditions_periodic()) {
          double const xnorm = norm1(x, nall)/nall; // g.comm
//        xnorm = xnorm/real( g%ng_all(1)*g%ng_all(2)*g%ng_all(3) )
          // subtract the average potential
          for(size_t i{0}; i < nall; ++i) x[i] -= xnorm;
      } // 3 periodic BCs
//       !============================================================

      if (0 == it % restart) {
          // |Ax> = A|x>
          ist = finite_difference::Laplacian(ax, x, g, fd);
          if (ist) error("CG_solve: Laplace failed!")
          // |r> = |b> - A|x> = |b> - |ax>
          set(r, nall, b); add_product(r, nall, ax, real_t(-1));
      } else {
          // |r> = |r> - alpha |ap>
          add_product(r, nall, ap, real_t(-alpha));
      } // restart?

      // res = <r|r>
      res2 = norm2(r, nall) * g.dV(); // g.comm

      // |z> = |Pr> = P|r>
//       ist = Precon_Nf4( g, r, z )
//       if( ist /= 0 ) stop 'CG_solve: Precon_Nf4 failed!'
      set(z, nall, r);

      // rz_new = <r|z>
      double const rz_new = scalar_product(r, z, nall) * g.dV(); // g.comm

      double constexpr RS_TINY = 1e-10;
      // beta = rz_new / rz_old
      double beta = rz_new / rz_old;
      if (rz_old < RS_TINY) {
          beta = 0;
          set(p, nall, z);
      } else {
          // |p> = |z> + beta |p>
          scale(p, nall, real_t(beta));
          add_product(p, nall, z, real_t(1));
      }

      if (echo > 7) printf("# %s it=%i alfa=%g beta=%g res=%.1e E=%.15f\n", __FILE__, 
        it, alpha, beta, std::sqrt(res2/cell_volume), scalar_product(x, b, nall));

      // rz_old = rz_new
      rz_old = rz_new;

      // decide if we continue to iterate
      run = (res2 > threshold2); // residual fell below threshold ?
      run = run && (it < maxiter); // maximum number of steps exceeded ?
      run = run || (it < miniter); // minimum number of steps not reached ?
//       !--------------------------------------
//       ! end of the CG iteration
//       !--------------------------------------
    } // while(run)

    double const res = std::sqrt(res2/cell_volume);
    if (residual) *residual = res; // export

    // show the result
    if (echo > 2) printf("# %s %.2e -> %.2e e/Bohr^3%s in %d%s iterations\n", __FILE__,
        res_start, res, (res < threshold)?" converged":"", it, (it < maxiter)?"":" (maximum)");
    
    return ist;
  } // solve


#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS
  
  template <typename real_t>
  status_t test_solver(int const echo=9, int const ng=24) {
      status_t stat{0};
      real_space_grid::grid_t<1> g(ng, ng, ng);
      view2D<real_t> xb(2, ng*ng*ng, 0.0);
      auto const x = xb[0], b = xb[1];
      double integral{0};
      double const cnt[3] = {.5*ng, .5*ng, .5*ng};
      for(int iz = 0; iz < ng; ++iz) {
      for(int iy = 0; iy < ng; ++iy) {
      for(int ix = 0; ix < ng; ++ix) {
          size_t const izyx = ix + ng*(iy + ng*iz);
          double const r2 = pow2(ix - cnt[0]) + pow2(iy - cnt[1]) + pow2(iz - cnt[2]);
          double const rho = std::exp(-.125*r2) - 8*std::exp(-.5*r2);
          b[izyx] = rho;
          integral += rho;
      }}} // ix iy iz
      if (echo > 2) printf("# %s integrated density %g\n", __FILE__, integral*g.dV());
      
      stat = solve(x, b, g);
      
      if (echo > 8) { // get a radial representation through Bessel transform
          float const dq = 1.f/16; int const nq = (int)(constants::pi/(g.smallest_grid_spacing()*dq));
          auto const rg = *radial_grid::create_equidistant_radial_grid(150, 15.);
          std::vector<real_t> q_coeff(nq, 0.0), f_radial(rg.n, 0.0);
          for(int i01 = 0; i01 < 2; ++i01) {
              real_space_grid::bessel_projection(q_coeff.data(), nq, dq, xb[i01], g, cnt);
              bessel_transform::transform_s_function(f_radial.data(), q_coeff.data(), rg, nq, dq, true); // transform back to real-space again
              printf("\n## r, %s (a.u.)\n", i01?"density":"V_electrostatic");
              for(int ir = 0; ir < rg.n; ++ir) {
                  printf("%g %g\n", rg.r[ir], f_radial[ir]);
              } // ir
          } // i01
      } // echo
      return stat;
  } // test_solver

  status_t all_tests(int const echo) {
    status_t stat{0};
    stat += test_solver<double>(echo);
    stat += test_solver<float>(echo); // test complation and convergence
    return stat;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace iterative_poisson
