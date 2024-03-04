// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, ::snprintf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <algorithm> // std::swap<T>
#include <cmath> // std::sqrt
#include <type_traits> // std::is_same

#include "iterative_poisson.hxx"

#include "real_space.hxx" // ::grid_t
#include "data_view.hxx" // view2D<T>
#include "inline_math.hxx" // set, dot_product, pow2, add_product
#include "finite_difference.hxx" // ::stencil_t, ::apply
#include "constants.hxx" // ::pi
#include "multi_grid.hxx" // ::restrict3D, ::interpolate3D, ::analyze_grid_sizes
#include "linear_algebra.hxx" // ::linear_solve
#include "boundary_condition.hxx" // Periodic_Boundary

#ifndef NO_UNIT_TESTS
  #include "real_space.hxx" // ::grid_t, ::Bessel_projection
  #include "radial_grid.h" // radial_grid_t
  #include "radial_grid.hxx" // ::create_radial_grid, ::equation_equidistant, ::destroy_radial_grid
  #include "bessel_transform.hxx" // ::transform_s_function
  #include "radial_potential.hxx" // ::Hartree_potential
  #include "fourier_poisson.hxx" // ::solve
  #include "control.hxx" // ::get
#endif

namespace iterative_poisson {
  // solve the Poisson equation iteratively using the conjugate gradients method + finite difference
  // or using a multi grid approach, all serial implementation

  double constexpr m1over4pi = -.25/constants::pi; // -4*constants::pi is the electrostatics prefactor in Hartree atomic units

  template <typename real_t>
  double norm2(real_t const v[], size_t const n) {
      double s{0};
      for (size_t i{0}; i < n; ++i) { 
          s += pow2(v[i]);
      } return s;
  } // norm2

  template <typename real_t>
  double norm1(real_t const v[], size_t const n) {
      double s{0};
      for (size_t i{0}; i < n; ++i) {
          s += v[i];
      }
      return s;
  } // norm1

  template <typename real_t>
  double scalar_product(real_t const v[], real_t const w[], size_t const n) 
      { return dot_product(n, v, w); }

  // Multi-grid method relies on a short range stencil [1 -2 1]
  // in order to speed this up, the analysis of the grid sizes could
  // run beforehand and be stored in a linked list (or vector) of grid descriptors
  // then, at the very bottom, the exactl solution use linear_solve (A*x == b),
  // which an explicit operator A (which currently is constructed every time again) 
  // this could be replaced by a matrix multiplication with A^{-1} --> faster than solving
  // A only depends on boundary_conditions, grid point numbers and grid spacings.
  // then, we can also in advance decide if we have to use the general restrict3D and interpolate3D
  // methods which are xyz --> Xyz --> XYz --> XYZ or specific interpolations routines
  // that do it in one step (e.g. for the case when all 3 grid dimensions are divisible by 2)

  template <typename real_t>
  double multi_grid_smoothen( 
        real_t x[] // on entry x, on exit slightly better to fullfil A*x == b
      , real_t r[] // on exit r == b - A*x residual vector
      , real_t const b[] // right hand side
      , real_space::grid_t const & g // Cartesian grid descriptor
      , int const echo=0 // log-level
  ) {
      size_t const n = g.all();
      finite_difference::stencil_t<real_t> const A(g.h, 1, m1over4pi); // 1:lowest order FD stencil
      double const c0 = A.c2nd[0][0] + A.c2nd[1][0] + A.c2nd[2][0]; // diagonal element
      double const omega = 2/3.; // Jacobi damping

      // Jacobi update method
      finite_difference::apply(r, x, g, A); // r = A*x
      add_product(r, n, b, real_t(-1)); // now r = A*x - b
      add_product(x, n, r, real_t(-omega/c0));
      scale(r, n, real_t(-1)); // now r = b - A*x
      if (3 == g.number_of_boundary_conditions(Periodic_Boundary)) {
          real_t const x_avg = norm1(x, n)/n;
          real_t const r_avg = norm1(r, n)/n;
          for (size_t i = 0; i < n; ++i) {
              x[i] -= x_avg;
              r[i] -= r_avg;
          } // i
      } // stabilize periodic
      return norm2(r, n);
  } // multi_grid_smoothen

  inline int multi_grid_level_number(real_space::grid_t const &g) {
      return int(std::ceil(std::log2(std::max(std::max(g[0], g[1]), g[2]))));
  } // multi_grid_level_number

  void multi_grid_level_label(char *label, real_space::grid_t const &g) {
      auto const lev = multi_grid_level_number(g);
      auto lab{label};
      for (int l = 0; l < 2*lev; ++l) *(lab++) = ' '; // indent two spaces per level to indicate the V-cyle visually in the log-output
      std::snprintf(lab, 16, "level=%i", lev);
  } // multi_grid_level_label

  template <typename real_t>
  status_t multi_grid_exact(real_t x[] // on entry x, on exit a slightly better to A*x == b
                  , real_t const b[] // right hand side
                  , real_space::grid_t const &g
                  , int const echo=0) {
      int const n = g.all();
      if (n > 8) return -1; // larger exact solutions not implemented
      if (n < 1) return -1; // strange
      int const peri = (g.number_of_boundary_conditions(Periodic_Boundary) == 3);

      // construct the right hand side
      double b8[8]; // get memory
      set(b8, n, b); // convert b into double (if real_t==float)
      if (peri) {
          double const b_avg = norm1(b8, n)/n;
          for (int i = 0; i < n; ++i) b8[i] -= b_avg; // make charge neutral
      } // periodic
#ifdef    DEVEL
      double c8[8]; set(c8, 8, b8);
      view2D<double> c88(8, 8, 0.0);
#endif // DEVEL

      // construct the explicit matrix operator
      finite_difference::stencil_t<real_t> const A(g.h, 1, m1over4pi); // 1:lowest order FD stencil
      view2D<double> a88(8, 8, 0.0); // get memory
      for (int iz = 0; iz < g[2]; ++iz) {
          for (int iy = 0; iy < g[1]; ++iy) {
              for (int ix = 0; ix < g[0]; ++ix) {
                  int const ixyz = (iz*g[1] + iy)*g[0] + ix;
                  assert(ixyz < n);
                  int const ixiyiz[] = {ix, iy, iz};
                  for (int d = 0; d < 3; ++d) {
                      int const id = ixiyiz[d];
                      for (int ij = -1; ij <= 1; ++ij) {
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
#ifdef    DEVEL
                          c88(ixyz,jxyz) = a88(ixyz,jxyz); // copy
#endif // DEVEL
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
          for (int i = 0; i < n; ++i) x8[i] -= x_avg; // make neutral
      } // periodic

#ifdef    DEVEL
      if (echo > 9) {
          std::printf("\n# %s for %dx%dx%d=%d status=%i\n", __func__, g[0],g[1],g[2], n, info);
          for (int i = 0; i < n; ++i) {
              double b8i{0}; // check that A_ij*x_j == b_i
              std::printf("#%3i ", i);
              for (int j = 0; j < n; ++j) {
                  std::printf(" %g", c88(i,j));
                  b8i += c88(i,j) * x8[j];
              } // j
              std::printf(" \tx=%g \tb=%g \tAx=%g \tx-<x>=%g\n", x8[i] + x_avg, c8[i], b8i, x8[i]);
          } // i
          std::printf("\n");
      } // echo
#endif // DEVEL

      if (0 == info) {
          set(x, n, x8); // convert x into real_t
      } // solver success
      return info;
  } // multi_grid_exact

  template <typename real_t>
  status_t multi_grid_cycle(real_t x[] // approx. solution to Laplace(x)/(-4*pi) == b
                            , real_t const b[] //
                            , real_space::grid_t const &g
                            , int const echo=0
                            , bool const as_solver=true
                            , float *residual=nullptr
                            , unsigned const n_pre=2
                            , unsigned const n_post=2
                            , float const tol2=1e-18) {
      // multi-grid V-cycle
      status_t stat(0);
      char label[96]; multi_grid_level_label(label, g);

      if (g.all() <= 8) return multi_grid_exact(x, b, g, echo);

      std::vector<real_t> residual_vector(g.all()); // get memory
      auto const r = residual_vector.data();
      real_t const *rd = r;

      if (as_solver) {
          float rn2{9e9}; // residual 2-norm
          for (int p = 0; p < n_pre; ++p) {
              rn2 = multi_grid_smoothen(x, r, b, g);
              if (echo > 19) std::printf("# %s %s  pre-smoothen step %i residual norm %.1e\n", __func__, label, p, rn2);
          } // p
          if (echo > 6) std::printf("# %s %s  pre-smoothen to residual norm %.1e\n", __func__, label, rn2);
      } else {
          rd = b;
          if (echo > 5) std::printf("# %s %s            input residual norm %.1e\n", __func__, label, norm2(rd, g.all()));
      } // n pre > 0

      uint32_t ngc[3];
      stat += multi_grid::analyze_grid_sizes(g, ngc); // find the next coarser 2^k grid
      real_space::grid_t gc(ngc);
      gc.set_boundary_conditions(g.boundary_conditions());
      gc.set_grid_spacing(g[0]*g.h[0]/ngc[0], g[1]*g.h[1]/ngc[1], g[2]*g.h[2]/ngc[2]);
      std::vector<real_t> rc(gc.all()), uc(gc.all(), real_t(0)); // two quantites on the coarse grid
      stat += multi_grid::restrict3D(rc.data(), gc, rd, g, 0); // mute
      stat += multi_grid_cycle(uc.data(), rc.data(), gc, echo - 2, true, nullptr, n_pre, n_post, tol2); // recursive invokation
      stat += multi_grid::interpolate3D(r, g, uc.data(), gc, 0); // mute
      add_product(x, g.all(), r, real_t(1)); // correct x

      if (as_solver && n_post > 0) {
          float rn2{9e9}; // residual 2-norm
          for (int p = 0; p < n_post && rn2 > tol2; ++p) {
              rn2 = multi_grid_smoothen(x, r, b, g);
              if (echo > 19) std::printf("# %s %s post-smoothen step %i residual norm %.1e\n", __func__, label, p, rn2);
          } // p
          if (echo > 5) std::printf("# %s %s post-smoothen to residual norm %.1e\n", __func__, label, rn2);
          if (residual) *residual = std::sqrt(rn2*g.dV()); // in units of the density
      } // n_post > 0

      return stat;
  } // multi_grid_cycle

  template <typename real_t>
  status_t multi_grid_solve(real_t x[] // one exit solution to Laplace(x)/(-4*pi) == b
                , real_t const b[] //
                , real_space::grid_t const &g
                , int const echo // =0 // log level
                , float const threshold=3e-8 // convergence criterion
                , float *residual=nullptr
                , int const maxiter=99 // maximum number of iterations
  ) {
      status_t stat(0);
      float res{9e9};
      for (int iter = 0; iter < maxiter && res > threshold; ++iter) {
//        if (echo > 0) std::printf("# %s iteration #%i\n", __func__, iter);
          stat += multi_grid_cycle(x, b, g, echo - 2, true, &res); // 'V'-cycle
          if (echo > 0) std::printf("# %s iteration #%i residual = %.1e a.u.\n", __func__, iter, res);
      } // iter
      if (residual) *residual = res;
      return stat;
  } // multi_grid_solve

  template <typename real_t>
  status_t multi_grid_precond(real_t x[] // approx. solution to Laplace(x)/(-4*pi) == b
                            , real_t const b[] //
                            , real_space::grid_t const &g
                            , int const echo=0) {
      return multi_grid_cycle(x, b, g, echo, false);
  } // multi_grid_precond

  template <typename real_t>
  status_t solve(real_t x[] // result to Laplace(x)/(-4*pi) == b
                , real_t const b[] // right hand side b
                , real_space::grid_t const &g // grid descriptor
                , char const method // use mixed precision as preconditioner
                , int const echo // =0 // log level
                , float const threshold // =3e-8 // convergence criterion
                , float *residual // =nullptr // residual that was reached
                , int const maxiter // =999 // maximum number of iterations 
                , int const miniter // =0   // minimum number of iterations
                , int restart // =4096 // number of iterations before restart, 1:steepest descent
                ) {

    size_t const nall = size_t(g[2])*size_t(g[1])*size_t(g[0]);

    status_t ist(0);
    
    restart = ('s' == method) ? 1 : std::max(1, restart);
    
    if (std::is_same<real_t, double>::value) {
        view2D<float> xb(2, nall, 0.0); // get memory
        auto const x32 = xb[0], b32 = xb[1];
        set(b32, nall, b); // convert to float
        set(x32, nall, x); // convert to float
        if (echo > 5) std::printf("# %s solve in <float> precision first\n", __FILE__);
        ist += solve(x32, b32, g, method, echo, threshold, residual, maxiter >> 4, miniter, restart);
        if (echo > 5) std::printf("# %s switch back to <double> precision\n", __FILE__);
        set(x, nall, x32); // convert to double
    } // real_t==double

    if ('m' == (method | 32)) {
        return ist + multi_grid_solve(x, b, g, echo, threshold, residual, maxiter);
    }
    
    int const nn_precond = 0; // -1; // 0:none, >0:stencil, <0:multi_grid (MG-precond does not work right...)
    bool const use_precond = (0 != nn_precond);

    view2D<real_t> mem(4 + use_precond, nall, 0.0); // get memory
    auto const r=mem[0], p=mem[1], ax=mem[2], ap=mem[3], z=use_precond?mem[4]:r;    

    finite_difference::stencil_t<real_t> const Laplacian(g.h, 8, m1over4pi); // 8: use a 17-point stencil
    
    finite_difference::stencil_t<real_t> precond;
    if (nn_precond > 0) {
        precond = finite_difference::stencil_t<real_t>(g.h, nn_precond);
        auto const nn = precond.nearest_neighbors();
        double nrm{0};
        for (int d = 0; d < 3; ++d) {
            for (int i = 0; i < nn[d]; ++i) {
                precond.c2nd[d][i] = std::abs(precond.c2nd[d][i]);
                nrm += precond.c2nd[d][i] * (1 + (i > 0));
            } // i
        } // d
        nrm = 1./nrm;
        for (int d = 0; d < 3; ++d) {
            for (int i = 0; i < nn[d]; ++i) {
                precond.c2nd[d][i] *= nrm;
            } // i
        } // d
        if (echo > 6) std::printf("# %s use a diffusion preconditioner with %d %d %d neighbors\n", 
                                __FILE__, nn[0], nn[1], nn[2]);
    } else if (use_precond) {
        if (echo > 4) std::printf("# %s use a multi-grid preconditioner starting at level %d\n", 
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


    // |Ax> := A|x>
    ist = finite_difference::apply(ax, x, g, Laplacian);
    if (ist) error("CG_solve: Laplacian failed with status %i", int(ist));

    // dump_to_file("cg_start", nall, x, nullptr, 1, 1, "x", echo);
    
    if (g.number_of_boundary_conditions(Periodic_Boundary) == 3) {
        double const bnorm = norm1(b, nall)/nall * g.dV(); // g.comm
        if (echo > 8) std::printf("# %s all boundary conditions are periodic but system is charged with %g electrons\n", __FILE__, bnorm);
    } // all_boundary_conditions_periodic

    // |r> = |b> - A|x> = |b> - |Ax>
    set(r, nall, b); add_product(r, nall, ax, real_t(-1));

    // res^2 = <r|r>
    double res2 = norm2(r, nall) * g.dV(); // g.comm
    double const res_start = std::sqrt(res2/cell_volume); // store starting residual
    if (echo > 8) std::printf("# %s start residual=%.1e\n", __FILE__, res_start);

    // |z> = |Pr> = P|r>
    if (use_precond) {
        ist = (nn_precond > 0) ? finite_difference::apply(z, r, g, precond) 
                               : multi_grid_precond(z, r, g, echo);
        if (ist) error("CG_solve: Preconditioner failed with status %i", int(ist));
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
        if (ist) error("CG_solve: Laplacian failed with status %i", int(ist));

        double const pAp = scalar_product(p, ap, nall) * g.dV(); // g.comm

        // alpha = rz_old / pAp
        double const alpha = (std::abs(pAp) < RZ_TINY) ? RZ_TINY : rz_old / pAp;

        // |x> = |x> + alpha |p>
        add_product(x, nall, p, real_t(alpha));

//       !============================================================
//       ! special treatment of completely periodic case
//       !============================================================
        if (g.number_of_boundary_conditions(Periodic_Boundary) == 3) {
            double const xnorm = norm1(x, nall)/nall; // g.comm
  //        xnorm = xnorm/real( g%ng_all(1)*g%ng_all(2)*g%ng_all(3) )
            // subtract the average potential
            for (size_t i{0}; i < nall; ++i) x[i] -= xnorm;
        } // 3 periodic BCs
//       !============================================================

        if (0 == it % restart) {
            // |Ax> = A|x>
            ist = finite_difference::apply(ax, x, g, Laplacian);
            if (ist) error("CG_solve: Laplace failed with status %i", int(ist))
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
            if (ist) error("CG_solve: Preconditioner failed with status %i", int(ist));
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

        if (echo > 9) std::printf("# %s it=%i alfa=%g beta=%g\n", __FILE__, it, alpha, beta);
        if (echo > 7) std::printf("# %s it=%i res=%.2e E=%.15f\n", __FILE__, it, 
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
    if (echo > 2) std::printf("# %s %.2e -> %.2e e/Bohr^3%s in %d%s iterations\n", __FILE__,
        res_start, res, (res < threshold)?" converged":"", it, (it < maxiter)?"":" (maximum)");
    if (echo > 5) std::printf("# %s inner product <x|b> = %.15f\n", __FILE__, scalar_product(x, b, nall) * g.dV());

    return (res > threshold);
  } // solve

#ifdef    NO_UNIT_TESTS
  template // explicit template instantiation for double
  status_t solve(double*, double const*, real_space::grid_t const &, char, int, float, float*, int, int, int);

  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

  template <typename real_t>
  status_t test_solver(int const echo=9, int const ng_default=32) {
      int const ng[] = {ng_default, ng_default, ng_default};
      real_space::grid_t g(ng[0], ng[1], ng[2]); // grid spacing = 1.0
      g.set_boundary_conditions(1); // all boundary conditions periodic, ToDo: fails for isolated BCs
      view2D<real_t> xb(3, ng[2]*ng[1]*ng[0], 0.0);
      int constexpr I_multigrid = 0, I_fft = 1, I_rho = 2;
      auto const x = xb[I_multigrid], x_fft = xb[I_fft], b = xb[I_rho];
      double constexpr c1 = 1, a1=.125, c2 = -8 + 1.284139e-7, a2=.5;
      double const cnt[] = {.5*ng[0], .5*ng[1], .5*ng[2]};
      double integral{0};
      for (int iz = 0; iz < ng[2]; ++iz) {
      for (int iy = 0; iy < ng[1]; ++iy) {
      for (int ix = 0; ix < ng[0]; ++ix) {
          size_t const izyx = (iz*ng[1] + iy)*ng[0] + ix;
          double const r2 = pow2(ix - cnt[0]) + pow2(iy - cnt[1]) + pow2(iz - cnt[2]);
          double const rho = c1*std::exp(-a1*r2) + c2*std::exp(-a2*r2);
          b[izyx] = rho;
          integral += rho;
      }}} // ix iy iz
      if (echo > 1) std::printf("# %s integrated density %g\n", __FILE__, integral*g.dV());

      float const threshold = (sizeof(real_t) > 4) ? 3e-8 : 5e-6;
      auto const method = control::get("parallel_poisson.test.method", "MultiGrid");
      auto const stat = solve(x, b, g, *method, echo, threshold); // method=M:multi_grid, 

      auto constexpr pi = constants::pi;
      double const mat[3][4] = {{2*pi/ng[0],0,0, 0},{0,2*pi/ng[1],0, 0}, {0,0,2*pi/ng[2], 0}};
      fourier_poisson::solve(x_fft, b, ng, mat);

      if (echo > 8) { // get a radial representation through Bessel transform
          auto & rg = *radial_grid::create_radial_grid(300, 15.f);

          int constexpr I_rho_radial = 3, I_hartree = 4, I_q2 = 5;
          view2D<double> fr(6, rg.n, 0.0); // radial functions
          for (int ir = 0; ir < rg.n; ++ir) {
              auto const r2 = pow2(rg.r[ir]);
              fr(I_rho_radial,ir) = c1*std::exp(-a1*r2) + c2*std::exp(-a2*r2);
          } // ir
          radial_potential::Hartree_potential(fr[I_hartree], rg, fr[I_rho_radial], rg.n, 0); // compute the Hartree potential by radial integration

          { // scope: Bessel-transform x and b
              float const dq = 1.f/128; // spacing in reciprocal space
              int const nq = int(constants::pi/(g.smallest_grid_spacing()*dq));
              std::vector<double> q_coeff(nq, 0.0);
              for (int xxb = 0; xxb <= 2; ++xxb) {
                  real_space::Bessel_projection(q_coeff.data(), nq, dq, xb[xxb], g, cnt, 1., 15.f);
                  bessel_transform::transform_s_function(fr[xxb], q_coeff.data(), rg, nq, dq, true); // transform back to real-space again
              } // x0b1
              // compute another reference solution by applying 4pi/q^2 using q_coeff of b
              for (int iq = 1; iq < nq; ++iq) {
                  q_coeff[iq] /= pow2(iq*dq);
              } // iq
              q_coeff[0] = 0;
              bessel_transform::transform_s_function(fr[I_q2], q_coeff.data(), rg, nq, dq, true); // transform back to real-space again
          } // scope

          double const f = 0.25/constants::pi; // 1/4pi
          std::printf("\n## r, V_mg, V_fft, V_rad, V_q^2, rho, rho_rad (all in a.u.)\n");
          for (int ir = 0; ir < rg.n; ++ir) {
              std::printf("%g %g %g %g %g %g %g\n", rg.r[ir], fr(I_multigrid,ir)*f, fr(I_fft,ir)*f, fr(I_hartree,ir), fr(I_q2,ir), fr(I_rho,ir)*f, fr(I_rho_radial,ir));
          } // ir
          std::printf("\n\n");
          radial_grid::destroy_radial_grid(&rg);

          std::printf("\n## r, V_mg, V_fft, rho (all in a.u.)\n"); // show all grid values (dots should not be connected by a line)
          for (int iz = 0; iz < ng[2]; ++iz) {
           for (int iy = 0; iy < ng[1]; ++iy) {
            for (int ix = 0; ix < ng[0]; ++ix) {
                size_t const izyx = (iz*ng[1] + iy)*ng[0] + ix;
                double const r2 = pow2(ix - cnt[0]) + pow2(iy - cnt[1]) + pow2(iz - cnt[2]);
                std::printf("%g %g %g %g\n", std::sqrt(r2), x[izyx], x_fft[izyx], b[izyx]); // point cloud
          }}} // ix iy iz

      } // echo
      return stat;
  } // test_solver

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_solver<double>(echo); // instantiation for both, double and float
 //   stat += test_solver<float>(echo);  // compilation and convergence tests
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace iterative_poisson
