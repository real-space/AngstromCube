#include <cstdio> // printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <algorithm> // std::swap<T>
#include <cmath> // std::sqrt
#include <type_traits> // std::is_same

#include "iterative_poisson.hxx"

#include "real_space.hxx" // ::grid_t
#include "data_view.hxx" // view2D<T>
#include "inline_math.hxx" // set, dot_product
#include "finite_difference.hxx" // ::stencil_t
#include "constants.hxx" // ::pi
#include "multi_grid.hxx" // ::restrict3D, ::interpolate3D, ::analyze_grid_sizes
#include "linear_algebra.hxx" // ::linear_solve

#ifndef NO_UNIT_TESTS
  #include "real_space.hxx" // ::bessel_projection
  #include "radial_grid.hxx" // ::radial_grid_t, ::create_equidistant_radial_grid
  #include "bessel_transform.hxx" // ::transform_s_function
#endif

namespace iterative_poisson {
  // solve the Poisson equation iteratively using the conjugate gradients method
  
  double constexpr m1over4pi = -.25/constants::pi; // -4*constants::pi is the electrostatics prefactor in Hartree atomic units

  template<typename real_t>
  double norm2(real_t const v[], size_t const n) { double s{0}; for(size_t i{0}; i < n; ++i) { s += pow2(v[i]); } return s; }
//   { return dot_product(n, v, v); }

  template<typename real_t>
  double norm1(real_t const v[], size_t const n) { double s{0}; for(size_t i{0}; i < n; ++i) { s += v[i]; } return s; }

  template<typename real_t>
  double scalar_product(real_t const v[], real_t const w[], size_t const n) { return dot_product(n, v, w); }

  
  template<typename real_t>
  double multi_grid_smoothen( real_t x[] // on entry x, on exit a slightly better to A*x == b
                            , real_t r[] // on exit r == b - A*x
                            , real_t const b[] // right hand side
                            , real_space::grid_t<1> const &g
                            , int const echo=0) {
      size_t const n = g.all();
      finite_difference::stencil_t<real_t> const A(g.h, 1, m1over4pi); // 1:lowest order FD stencil
      double const c0 = A.c2nd[0][0] + A.c2nd[1][0] + A.c2nd[2][0]; // diagonal element
      double const omega = 2/3.; // Jacobi damping

      // Jacobi update method
      finite_difference::apply(r, x, g, A); // r = A*x
      add_product(r, n, b, real_t(-1)); // now r = A*x - b
      add_product(x, n, r, real_t(-omega/c0));
      scale(r, n, real_t(-1)); // now r = b - A*x
      if (g.all_boundary_conditions_periodic()) {
          real_t const x_avg = norm1(x, n)/n;
          real_t const r_avg = norm1(r, n)/n;
          for(size_t i = 0; i < n; ++i) {
              x[i] -= x_avg;
              r[i] -= r_avg;
          } // i
      } // stabilize periodic
      return norm2(r, n);
  } // multi_grid_smoothen
  
  inline int multi_grid_level_number(real_space::grid_t<1> const &g) {
      return int(std::ceil(std::log2(std::max(std::max(g[0], g[1]), g[2]))));
  } // multi_grid_level_number
      
  void multi_grid_level_label(char *label, real_space::grid_t<1> const &g) {
      auto const lev = multi_grid_level_number(g);
      auto lab{label};
      for(int l = 0; l < 2*lev; ++l) *(lab++) = ' '; // indent two spaces per level to indicate the V-cyle visually in the log-output
      std::sprintf(lab, "level=%i", lev);
  } // multi_grid_level_label

  template<typename real_t>
  status_t multi_grid_exact(real_t x[] // on entry x, on exit a slightly better to A*x == b
                  , real_t const b[] // right hand side
                  , real_space::grid_t<1> const &g
                  , int const echo=0) {
      int const n = g.all();
      if (n > 8) return -1; // larger exact solutions not implemented
      if (n < 1) return -1; // strange
      int const peri = g.all_boundary_conditions_periodic() ? 1 : 0;

      // construct the right hand side
      double b8[8]; // get memory
      set(b8, n, b); // convert b into double (if real_t==float)
      if (peri) {
          double const b_avg = norm1(b8, n)/n;
          for(int i = 0; i < n; ++i) b8[i] -= b_avg; // make charge neutral
      } // periodic
#ifdef DEVEL
      double c8[8]; set(c8, 8, b8);
      view2D<double> c88(8, 8, 0.0);
#endif

      // construct the explicit matrix operator
      finite_difference::stencil_t<real_t> const A(g.h, 1, m1over4pi); // 1:lowest order FD stencil
      view2D<double> a88(8, 8, 0.0); // get memory
      for(int iz = 0; iz < g[2]; ++iz) {
          for(int iy = 0; iy < g[1]; ++iy) {
              for(int ix = 0; ix < g[0]; ++ix) {
                  int const ixyz = (iz*g[1] + iy)*g[0] + ix;
                  assert(ixyz < n);
                  int const ixiyiz[] = {ix, iy, iz};
                  for(int d = 0; d < 3; ++d) {
                      int const id = ixiyiz[d];
                      for(int ij = -1; ij <= 1; ++ij) {
                          double f{1};
                          int jd = id + ij;
                          if (0 == g.boundary_condition(d)) {
                              if (jd >= g[d]) f = 0;
                              if (jd < 0)     f = 0;
                          } // non-periodic boundaries
                          jd = (jd + 16*g[d]) % g[d]; // periodic wrap-around
                          int jxjyjz[] = {ix, iy, iz};
                          jxjyjz[d] = jd;
                          int const jxyz = (jxjyjz[2]*g[1] + jxjyjz[1])*g[0] + jxjyjz[0];
                          assert(jxyz < n);
                          a88(ixyz,jxyz) += f*A.c2nd[d][std::abs(ij)];
#ifdef DEVEL
                          c88(ixyz,jxyz) = a88(ixyz,jxyz); // copy
#endif
                      } // ij
                  } // d
              } // ix
          } // iy
      } // iz
      
      // solve
//       auto const info = linear_algebra::linear_solve(n - peri, a88.data(), a88.stride(), b8, 8); // if(peri) we solve for one degrees of freedom less
      auto const info = linear_algebra::linear_solve(n, a88.data(), a88.stride(), b8, 8);
      double *x8 = b8;
      // if (peri) x8[n - peri] = -norm1(x8, n - peri); // construct x8 neutral, i.e. such that norm1(x8) == 0

      double x_avg{0};
      if (peri) {
          x_avg = norm1(x8, n)/n;
          for(int i = 0; i < n; ++i) x8[i] -= x_avg; // make neutral
      } // periodic

#ifdef DEVEL
      if (echo > 9) {
          printf("\n# %s for %dx%dx%d=%d status=%i\n", __func__, g[0],g[1],g[2], n, info);
          for(int i = 0; i < n; ++i) {
              double b8i{0}; // check that A_ij*x_j == b_i
              printf("#%3i ", i);
              for(int j = 0; j < n; ++j) {
                  printf(" %g", c88(i,j));
                  b8i += c88(i,j) * x8[j];
              } // j
              printf(" \tx=%g \tb=%g \tAx=%g \tx-<x>=%g\n", x8[i] + x_avg, c8[i], b8i, x8[i]);
          } // i
          printf("\n");
      } // echo
#endif

      if (0 == info) {
          set(x, n, x8); // convert x into real_t
      } // solver success
      return info;
  } // multi_grid_exact
  
  template<typename real_t>
  status_t multi_grid_cycle(real_t x[] // to Laplace(x)/(-4*pi) == b
                            , real_t const b[] //
                            , real_space::grid_t<1> const &g
                            , int const echo=0
                            , bool const as_solver=true
                            , unsigned const n_pre=2
                            , unsigned const n_post=2
                            , float const tol2=1e-18) {
      // multi-grid V-cycle
      status_t stat(0);
      char label[96]; multi_grid_level_label(label, g);

      if (g.all() <= 8) return multi_grid_exact(x, b, g, echo);
      
      std::vector<real_t> residual_vector(g.all()); // get memory
      auto const r = residual_vector.data();
      real_t const *rd;

      if (as_solver) {
          float rn2{9e9}; // residual 2-norm
          for(int p = 0; p < n_pre; ++p) {
              rn2 = multi_grid_smoothen(x, r, b, g);
              if (echo > 19) printf("# %s %s  pre-smoothen step %i residual norm %.1e\n", __func__, label, p, rn2);
          } // p
          if (echo > 0) printf("# %s %s  pre-smoothen to residual norm %.1e\n", __func__, label, rn2);
          rd = r;
      } else {
          rd = b;
          if (echo > 0) printf("# %s %s            input residual norm %.1e\n", __func__, label, norm2(rd, g.all()));
      } // n pre > 0

      uint32_t ngc[3];
      stat += multi_grid::analyze_grid_sizes(g, ngc); // find the next coarser 2^k grid
      real_space::grid_t<1> gc(ngc);
      gc.set_boundary_conditions(g.boundary_conditions());
      gc.set_grid_spacing(g[0]*g.h[0]/ngc[0], g[1]*g.h[1]/ngc[1], g[2]*g.h[2]/ngc[2]);
      std::vector<real_t> rc(gc.all()), uc(gc.all(), real_t(0)); // two quantites on the coarse grid
      stat += multi_grid::restrict3D(rc.data(), gc, rd, g, 0); // mute
      stat += multi_grid_cycle(uc.data(), rc.data(), gc, echo - 1, true, n_pre, n_post, tol2); // recursive invokation
      stat += multi_grid::interpolate3D(r, g, uc.data(), gc, 0); // mute
      add_product(x, g.all(), r, real_t(1)); // correct x

      if (as_solver && n_post > 0) {
          float rn2{9e9}; // residual 2-norm
          for(int p = 0; p < n_post && rn2 > tol2; ++p) {
              rn2 = multi_grid_smoothen(x, r, b, g);
              if (echo > 19) printf("# %s %s post-smoothen step %i residual norm %.1e\n", __func__, label, p, rn2);
          } // p
          if (echo > 0) printf("# %s %s post-smoothen to residual norm %.1e\n", __func__, label, rn2);
      } // n_post > 0

      return stat;
  } // multi_grid_cycle
  
  template<typename real_t>
  status_t multi_grid_precond(real_t x[] // approx. solution to Laplace(x)/(-4*pi) == b
                            , real_t const b[] //
                            , real_space::grid_t<1> const &g
                            , int const echo=0) {
      return multi_grid_cycle(x, b, g, echo, false);
  } // multi_grid_precond
  
  template<typename real_t>
  status_t solve(real_t x[] // result to Laplace(x)/(-4*pi) == b
                , real_t const b[] // right hand side b
                , real_space::grid_t<1> const &g // grid descriptor
                , int const echo // =0 // log level
                , float const threshold // =3e-8 // convergence criterion
                , float *residual // =nullptr // residual that was reached
                , int const maxiter // =999 // maximum number of iterations 
                , int const miniter // =0   // minimum number of iterations
                , int restart // =4096 // number of iterations before restrat, 1:steepest descent
                , char const mixed_precision // use mixed precision as preconditioner
                ) {

    restart = std::max(1, restart);
    
    size_t const nall = size_t(g[2])*size_t(g[1])*size_t(g[0]);
    
    if (std::is_same<real_t, double>::value) {
        if ('m' == mixed_precision) {
            view2D<float> xb(2, nall, 0.0); // get memory
            auto const x32 = xb[0], b32 = xb[1];
            set(b32, nall, b); // convert to float
            set(x32, nall, x); // convert to float
            if (echo > 5) printf("# %s solve in <float> precision first\n", __FILE__);
            solve(x32, b32, g, echo + 1, threshold, residual, maxiter >> 4, miniter, restart);
            set(x, nall, x32); // convert to double
        } // mixed precision
    } // real_t==double

    int const nn_precond = -1; // 0:none, >0:stencil, <0:multi_grid
    bool const use_precond = (0 != nn_precond);
    
    view2D<real_t> mem(4 + use_precond, nall, 0.0); // get memory
    auto const r=mem[0], p=mem[1], ax=mem[2], ap=mem[3], z=use_precond?mem[4]:r;    
    
    finite_difference::stencil_t<real_t> const Laplacian(g.h, 8, m1over4pi); // 8: use a 17-point stencil
    
    finite_difference::stencil_t<real_t> precond;
    if (nn_precond > 0) {
        precond = finite_difference::stencil_t<real_t>(g.h, nn_precond);
        auto const nn = precond.nearest_neighbors();
        double nrm{0};
        for(int d = 0; d < 3; ++d) {
            for(int i = 0; i < nn[d]; ++i) {
                precond.c2nd[d][i] = std::abs(precond.c2nd[d][i]);
                nrm += precond.c2nd[d][i] * (1 + (i > 0));
            } // i
        } // d
        nrm = 1./nrm;
        for(int d = 0; d < 3; ++d) {
            for(int i = 0; i < nn[d]; ++i) {
                precond.c2nd[d][i] *= nrm;
            } // i
        } // d
        if (echo > 6) printf("# %s use a diffusion preconditioner with %d %d %d neighbors\n", 
                                __FILE__, nn[0], nn[1], nn[2]);
    } else if (use_precond) {
        if (echo > 4) printf("# %s use a multi-grid preconditioner starting at level %d\n", 
                                __FILE__, multi_grid_level_number(g));
    } // use_precond
    
    double const cell_volume = nall*g.dV();
    double const threshold2 = cell_volume * pow2(threshold);
    double constexpr RZ_TINY = 1e-14, RS_TINY = 1e-10;

    
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
    ist = finite_difference::apply(ax, x, g, Laplacian);
    if (ist) error("CG_solve: Laplacian failed!");

    // dump_to_file("cg_start", nall, x, nullptr, 1, 1, "x", echo);
    
    if (g.all_boundary_conditions_periodic()) {
        double const bnorm = norm1(b, nall)/nall * g.dV(); // g.comm
        if (echo > 8) printf("# %s all boundary conditions are periodic but system is charged with %g electrons\n", __FILE__, bnorm);
    } // all_boundary_conditions_periodic

    // |r> = |b> - A|x> = |b> - |Ax>
    set(r, nall, b); add_product(r, nall, ax, real_t(-1));

    // res^2 = <r|r>
    double res2 = norm2(r, nall) * g.dV(); // g.comm
    double const res_start = std::sqrt(res2/cell_volume); // store staring residual
    if (echo > 8) printf("# %s start residual=%.1e\n", __FILE__, res_start);

    // |z> = |Pr> = P|r>
    if (use_precond) {
        ist = (nn_precond > 0) ? finite_difference::apply(z, r, g, precond) 
                               : multi_grid_precond(z, r, g, echo);
        if (ist) error("CG_solve: Preconditioner failed!");
    } else assert(z == r);

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
        ist = finite_difference::apply(ap, p, g, Laplacian);
        if (ist) error("CG_solve: Laplacian failed!");

        double const pAp = scalar_product(p, ap, nall) * g.dV(); // g.comm

        // alpha = rz_old / pAp
        double const alpha = (std::abs(pAp) < RZ_TINY) ? RZ_TINY : rz_old / pAp;

        // |x> = |x> + alpha |p>
        add_product(x, nall, p, real_t(alpha));

//       !============================================================
//       ! special treatment of completely periodic case
//       !============================================================
        if (g.all_boundary_conditions_periodic()) {
            double const xnorm = norm1(x, nall)/nall; // g.comm
  //        xnorm = xnorm/real( g%ng_all(1)*g%ng_all(2)*g%ng_all(3) )
            // subtract the average potential
            for(size_t i{0}; i < nall; ++i) x[i] -= xnorm;
        } // 3 periodic BCs
//       !============================================================

        if (0 == it % restart) {
            // |Ax> = A|x>
            ist = finite_difference::apply(ax, x, g, Laplacian);
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
        if (use_precond) {
            ist = (nn_precond > 0) ? finite_difference::apply(z, r, g, precond) 
                                   : multi_grid_precond(z, r, g, echo);
            if (ist) error("CG_solve: Preconditioner failed!");
        } else assert(z == r);

        // rz_new = <r|z>
        double const rz_new = scalar_product(r, z, nall) * g.dV(); // g.comm

        // beta = rz_new / rz_old
        double beta = rz_new / rz_old;
        if (rz_old < RS_TINY) {
            beta = 0;
            set(p, nall, z);
        } else {
            // |p> = |z> + beta |p>
            scale(p, nall, real_t(beta));
            add_product(p, nall, z, real_t(1));
        } // rz_old < tiny

        if (echo > 9) printf("# %s it=%i alfa=%g beta=%g\n", __FILE__, it, alpha, beta);
        if (echo > 7) printf("# %s it=%i res=%.2e E=%.15f\n", __FILE__, it, 
            std::sqrt(res2/cell_volume), scalar_product(x, b, nall) * g.dV());

        // rz_old = rz_new
        rz_old = rz_new;

        // decide if we continue to iterate
        run = (res2 > threshold2); // residual fell below threshold ?
        run = run || (it < miniter); // minimum number of steps not reached ?
        run = run && (it < maxiter); // maximum number of steps exceeded ?
//       !--------------------------------------
//       ! end of the CG iteration
//       !--------------------------------------
    } // while(run)

    double const res = std::sqrt(res2/cell_volume);
    if (residual) *residual = res; // export

    // show the result
    if (echo > 2) printf("# %s %.2e -> %.2e e/Bohr^3%s in %d%s iterations\n", __FILE__,
        res_start, res, (res < threshold)?" converged":"", it, (it < maxiter)?"":" (maximum)");
    if (echo > 5) printf("# %s inner product <x|b> = %.15f\n", __FILE__, scalar_product(x, b, nall) * g.dV());

    return (res > threshold);
  } // solve

#ifdef  NO_UNIT_TESTS
  template // explicit template instantiation
  status_t solve(double x[], double const b[], real_space::grid_t<1> const &g // grid descriptor
                , int const echo, float const threshold, float *residual, int const maxiter
                , int const miniter, int restart, char const mixed_precision);

  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template <typename real_t>
  status_t test_solver(int const echo=9, int const ng=24) {
      status_t stat{0};
      real_space::grid_t<1> g(ng, ng, ng);
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

      stat = solve(x, b, g, echo);
      
      if (echo > 8) { // get a radial representation through Bessel transform
          float const dq = 1.f/16; int const nq = (int)(constants::pi/(g.smallest_grid_spacing()*dq));
          auto const rg = *radial_grid::create_equidistant_radial_grid(150, 15.);
          std::vector<real_t> q_coeff(nq, 0.0), f_radial(rg.n, 0.0);
          for(int i01 = 0; i01 < 2; ++i01) {
              real_space::bessel_projection(q_coeff.data(), nq, dq, xb[i01], g, cnt);
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
    stat += test_solver<float>(echo); // test compilation and convergence for float
    return stat;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace iterative_poisson
