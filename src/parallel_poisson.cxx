// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, ::snprintf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <algorithm> // std::swap<T>
#include <cmath> // std::sqrt
#include <type_traits> // std::is_same

#include "parallel_poisson.hxx"

#include "real_space.hxx" // ::grid_t
#include "data_view.hxx" // view2D<T>
#include "inline_math.hxx" // set, dot_product, pow2
#include "finite_difference.hxx" // ::stencil_t, ::apply
#include "constants.hxx" // ::pi
#include "linear_algebra.hxx" // ::linear_solve
#include "boundary_condition.hxx" // Periodic_Boundary

#ifndef NO_UNIT_TESTS
  #include "real_space.hxx" // ::Bessel_projection
  #include "radial_grid.hxx" // radial_grid_t, ::create_radial_grid, ::equation_equidistant, ::destroy_radial_grid
  #include "bessel_transform.hxx" // ::transform_s_function
  #include "radial_potential.hxx" // ::Hartree_potential
  #include "fourier_poisson.hxx" // ::solve
  #include "control.hxx" // ::get
#endif

namespace parallel_poisson {
  // solve the Poisson equation iteratively using the conjugate gradients method
  
  double constexpr m1over4pi = -.25/constants::pi; // -4*constants::pi is the electrostatics prefactor in Hartree atomic units

  template <typename real_t>
  double norm2(real_t const v[], size_t const n) {
      double s{0};
      for (size_t i{0}; i < n; ++i) { 
          s += pow2(v[i]);
      } // i 
      return s;
  } // norm2

  template <typename real_t>
  double norm1(real_t const v[], size_t const n) {
      double s{0};
      for (size_t i{0}; i < n; ++i) {
          s += v[i];
      } // i
      return s;
  } // norm1

  template <typename real_t>
  double scalar_product(real_t const v[], real_t const w[], size_t const n) {
      double dot{0};
      for (size_t i = 0; i < n; ++i) {
          dot += double(v[i])*double(w[i]); // conversion to double is different from dot_product define in inline_math.hxx
      } // i
      return dot;
  } // scalar_product
  
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
        ist += solve(x32, b32, g, method, echo, threshold, residual, maxiter, miniter, restart);
        if (echo > 5) std::printf("# %s switch back to <double> precision\n", __FILE__);
        set(x, nall, x32); // convert to double
    } // real_t==double

    // from here we use CG + order-16 FD

    int const nn_precond = 0; // 0:none, >0:range-1-stencil, <0:multi_grid (does not work properly yet)
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
    double const res_start = std::sqrt(res2/cell_volume); // store staring residual
    if (echo > 8) std::printf("# %s start residual=%.1e\n", __FILE__, res_start);

    // |z> = |Pr> = P|r>
    if (use_precond) {
        ist = finite_difference::apply(z, r, g, precond);
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
            if (ist) error("CG_solve: Laplacian failed with status %i", int(ist))
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
            ist = finite_difference::apply(z, r, g, precond);
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
            std::sqrt(res2/cell_volume), scalar_product(x, b, nall)*g.dV());

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
    if (echo > 5) std::printf("# %s inner product <x|b> = %.15f\n", __FILE__, scalar_product(x, b, nall)*g.dV());

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
      double constexpr c1 = 1, a1=.125, c2 = -8 + 1.284139e-7, a2=.5; // parameters for two Gaussians, in total close to neutral
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
      auto const method = control::get("parallel_poisson.test.method", "mix");
      int  const max_it = control::get("parallel_poisson.test.maxiter", 999.);
      float residual_reached{0};
      auto const stat = solve(x, b, g, *method, echo, threshold, &residual_reached, max_it);

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

} // namespace parallel_poisson
