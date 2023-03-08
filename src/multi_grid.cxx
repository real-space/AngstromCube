#include <cstdio> // std::printf
#include <cmath> // ::log10
#include <algorithm> // std::min, ::max
#include <vector> // std::vector<T>
#include <cstdint> // uint32_t

#include "multi_grid.hxx"

#ifndef   NO_UNIT_TESTS
  #include <cstdio> // std::fprintf, ::fopen, stdout, ::fclose
  #include <numeric> // std::iota
  #include "constants.hxx" // ::pi
  #include "simple_math.hxx" // ::random
  #include "debug_output.hxx" // dump_to_file
  #include "control.hxx" // ::get
#endif // NO_UNIT_TESTS

#include "real_space.hxx" // ::grid_t
#include "data_view.hxx" // view2D<T>
#include "inline_math.hxx" // set, add_product
#include "status.hxx" // status_t

namespace multi_grid {
  
  /* 
      Use multi-grid acceleration
      e.g. like this: Kohn-Sham is solved on a grid with ng grid points
          then usually the potential is constructed on 2*ng grid points.
          We need to find the largest integer k so that 2^k < ng.
          Then, we only have a single interface between coarse and fine grid
          where the refinement/coarsening factor is not 2.
          Suggestion: as simple as possible 
          (since the coarser levels are only preconditioners)
          compute in float and
          restrict from 2x finer level by c0=(f0 + f1)/2, c1=(f2 + f3)/2, ... (local operation, charge conserving)
          prolong  from 2x coarser level by f0=c0, f1=c0, f2=c1, f3=c1, ...   (local operation)
          and if we hit the non-factor-2-interface:
          example: restrict from 5 --> 3 (just for illustration, the coarser grid should actually be 2^k)
           |  0  |  1  |  2  |  3  |  4  |
           |    0    |    1    |    2    |
          with a banded (charge conserving) operator
          [.6 .4  0  0  0]
          [ 0 .2 .6 .2  0]
          [ 0  0  0 .4 .6]          
          and prolong with linear interpolation.
          Both operations need some data-exchange.
   */

  inline unsigned nearest_binary_power(unsigned const ng) {
      int k{0};
      while(ng > (1ull << (k + 1))) ++k;
      return k;
  } // nearest_binary_power


  status_t analyze_grid_sizes(
        real_space::grid_t const & g // coarse grid where the Kohn-Sham equation is typically solved
      , uint32_t *n_coarse
      , int const echo // =0
  ) {
      for (int d = 0; d < 3; ++d) {
          unsigned const ng = g[d];
          unsigned const k = nearest_binary_power(ng);
          unsigned const nb = 1ull << k;
          if (nb < ng) {
              if (echo > 2) std::printf("# grid in %c-direction can be coarsened from %d to %d = 2^%i grid points\n", 'x'+d, ng, nb, k);
              if (n_coarse) n_coarse[d] = nb;
          } else {
              // cannot be coarsened further
              if (n_coarse) n_coarse[d] = ng;
          }
      } // d
      return 0;
  } // analyze_grid_sizes

  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  int64_t grid_point_id(uint32_t const x, uint32_t const y, uint32_t const z) {
      assert(x < (1 << 21)); 
      assert(y < (1 << 21)); 
      assert(z < (1 << 21));
      int64_t id{0}; // 3D Morton number
      for (int b = 0; b < 21; ++b) {
          uint64_t const p = 1u << b; // probe bit #b
          id |= ((x & p) << (2*b + 0)); // bit #b of x becomes bit #3b   of id
          id |= ((y & p) << (2*b + 1)); // bit #b of y becomes bit #3b+1 of id
          id |= ((z & p) << (2*b + 2)); // bit #b of z becomes bit #3b+2 of id
      } // b
      return id;
  } // grid_point_id

  status_t grid_point_xyz(int32_t & x, int32_t & y, int32_t & z, int64_t const id) {
      if (id < 0) { x = -1; y = -1; z = -1; return -1; }
      x = 0; y = 0; z = 0;
      int64_t const one = 1;
      for (int b = 0; b < 21; ++b) {
          x |= ((id & (one << (3*b + 0))) >> (2*b + 0));
          y |= ((id & (one << (3*b + 1))) >> (2*b + 1));
          z |= ((id & (one << (3*b + 2))) >> (2*b + 2));
      } // b
      return 0;
  } // grid_point_xyz

  status_t test_Morton_indices(int const echo=0) {
      status_t stat(0);
      for (int n = 0; n < 99; ++n) {
          uint32_t xyz[3]; for (int d = 0; d < 3; ++d) xyz[d] = simple_math::random(0, 9999);
          int32_t abc[3];
          auto const id = grid_point_id(xyz[0], xyz[1], xyz[2]);
          stat += grid_point_xyz(abc[0], abc[1], abc[2], id);
          for (int d = 0; d < 3; ++d) stat += (abc[d] != xyz[d]);
          if (echo > 9) std::printf("# %s %5d %5d %5d --> %lld --> %5d %5d %5d\n",
              __func__, xyz[0], xyz[1], xyz[2], id, abc[0], abc[1], abc[2]);
      } // n
      int32_t abc[3];
      stat += (grid_point_xyz(abc[0], abc[1], abc[2], -9) != -1);
      stat +=   (abc[0] != -1) + (abc[1] != -1) + (abc[2] != -1);
      return stat;
  } // test_Morton_indices

  status_t check_general_restrict(unsigned const ng, int const echo=0, char const dir='?') {
      unsigned const k = nearest_binary_power(ng);
      unsigned const nb = 1ull << k;
      if (echo > 2) std::printf("# grid in %c-direction can be coarsened from %d to %d = 2^%i grid points\n", dir, ng, nb, k);
      // compute the restriction operator from ng to n2
      assert(nb < ng);
      double const rowref = 1; // each row sum must be equal to 1
      std::vector<double> colsum(ng, 0.0);
      double const b_over_g = nb/double(ng);
      double const colref = b_over_g; // each column sum must be nb/ng
      double rowdev{0};
      for (int ib = 0; ib < nb; ++ib) {
          if (echo > 8) std::printf("# %c row%4i  ", dir, ib);
          double rowsum{0};
          for (int ig = 0; ig < ng; ++ig) {
              // overlap between a grid compartment of length L/ng starting at ig*L/ng
              // with a grid compartment of length L/n2 starting at i2*L/n2
              //            |  0  |  1  |  2  |  3  |  4  |   g
              //            |    0    |    1    |    2    |   b
              // (the cell length L falls out of the equation)
              double const start = std::max((ig - 0)*b_over_g, double(ib - 0));
              double const end   = std::min((ig + 1)*b_over_g, double(ib + 1));
              double const ovl = std::max(0.0, end - start);
              if (ovl > 0) {
                  if (echo > 8) std::printf("  %i %g", ig, ovl); // show only non-zero matrix entries (with their column index) 
                  rowsum += ovl;
                  colsum[ig] += ovl;
              } // ovl non-zero
              // thid way to restrict requires no MPI communication across the outer boundaries, 
              // only communication between inner boundaries and onwership redistribution.
              // Ownership redistribution becomes necessary since at some point,
              // the communication times are larger than processing.
              // Or there might be a finit-difference stencil with nn > 1,
              // then, there is a natural limit of how few grid points per process elements
              // can be operatored. Then, we should redistribute to less process elements.
              // finally solving on the master only for the coarsest level.
          } // id
          if (echo > 8) std::printf(" sum=%g dev=%.1e\n", rowsum, rowsum - rowref);
          rowdev = std::max(rowdev, std::abs(rowsum - rowref));
      } // ib
      double coldev{0}; // column deviation
      for (int ig = 0; ig < ng; ++ig) {
          coldev = std::max(coldev, std::abs(colsum[ig] - colref));
      } // ig
      if (echo > 6) std::printf("# %c-direction largest deviation = %.1e (row) and %.1e (col)\n\n", dir, rowdev, coldev);
      return (coldev > 1e-14) + (rowdev > 1e-14);
      
      // in some cases only two coefficients per row are non-zero
      // (in particular the case of power of 2: coefficients are 0.5 0.5)
      // in general, there are up to three non-zero coefficients.
      
  } // check_general_restrict


  status_t test_transfer(int const echo=0) {
      if (echo < 1) return 0;
      status_t stat(0);
      int const ng = 32;
      std::vector<double> input_c(ng), result_c(ng), result_cdc(ng);
      
      std::vector<int> dense_grids{ng, (5*ng)/3, 2*ng}; // numbers of dense grid points
      for (auto mg : dense_grids) {
          std::vector<double> result_d(mg), input_d(mg);
          std::printf("\n## ik T (%s: interpolate from grid=%d to grid=%d)\n", __func__, ng, mg); // legend
          for (int ik = 1; ik < 3; ++ik) {
              double const k = (2*constants::pi*ik);
              for (int ig = 0; ig < ng; ++ig) {
                  double const x = (ig + 0.5)/ng;
                  input_c[ig] = std::cos(k*x);
              } // ig
              
              stat += linear_interpolation(result_d.data(), mg, input_c.data(), ng, 1, 1);
              
              std::printf("\n## x interpolate(cos(x)) cos(x) [n=%d m=%d k=%i]\n", ng, mg, ik);
              for (int jg = 0; jg < mg; ++jg) {
                  double const x = (jg + 0.5)/mg;
                  input_d[jg] = std::cos(k*x);
                  std::printf("%g %g %g\n", x, result_d[jg], input_d[jg]);
              } // jg
              
              stat += restrict_to_any_grid(result_c.data(), ng, input_d.data(), mg, 1, 1);
              stat += restrict_to_any_grid(result_cdc.data(), ng, result_d.data(), mg, 1, 1);

              std::printf("\n## x cos(x) restrict(interpolate(cos(x))) [n=%d m=%d k=%i]\n", ng, mg, ik);
              for (int ig = 0; ig < ng; ++ig) {
                  double const x = (ig + 0.5)/ng;
                  std::printf("%g %g %g %g\n", x, result_c[ig], result_cdc[ig], input_c[ig]);
              } // ig
              
          } // ik
      } // mg
      return stat;
  } // test_transfer
  
  status_t test_analysis(int const echo=0) {
      status_t stat(0);
      real_space::grid_t g(63, 64, 65);
      for (char dir = 'x'; dir <= 'z'; ++dir) {
          stat += check_general_restrict(g(dir), echo, dir);
      } // direction
      return stat + analyze_grid_sizes(g, nullptr, echo);
  } // test_analysis
  
  status_t test_restrict_interpolate(int const echo=0) {
      status_t stat(0);
      real_space::grid_t gi(15, 16, 17), go(8, 8, 16);
      std::vector<float> in(gi.all()), out(go.all());
      std::iota(in.begin(), in.end(), .5f);
      stat += restrict3D(out.data(), go, in.data(), gi);
      stat += interpolate3D(in.data(), gi, out.data(), go);
      return stat;
  } // test_restrict_interpolate

  
  
  
  
  
  template <typename real_t>
  double get_residual_norm(
        real_t const r[]
      , real_t const x[]
      , real_t const b[]
      , size_t const g
      , double const h
      , FILE* f=stdout
  ) {
      if (f) std::fprintf(f, "\n## %s i, x[i], residual[i], rhs[i]  for i < %ld\n", __func__, g);
      double const c0 = 2./(h*h), c1 = -.5*c0;
      double norm2{0}, normT{0};
      real_t r_prev = r[g - 1]; // init
      for (int i = 0; i < g; ++i) {
          if (f) std::fprintf(f, "%g %g %g %g\n", i*h, x[i], r[i], b[i]);
          real_t const ri = r[i];
          real_t const r_next = r[(i + 1) % g];
          double const Tr = c1*r_prev + c0*ri + c1*r_next;
          norm2 += ri*ri; // accumulate residual norm ||r||_2
          normT += Tr*Tr; // accumulate norm ||D^2 r||_2
      } // i
      if (f) std::fprintf(f, "%g %g %g %g\n\n", g*h, x[0], r[0], b[0]); // periodic
      if (0) std::printf("# residual norm (length %ld) is %g / %g = %g\n", g, normT, norm2, normT/norm2);
      return norm2/g;
  } // get_resiual_norm

  // toy model for testing the mg_cycle structure, operator A = stencil {1, -2, 1}/h^2
  template <typename real_t>
  double jacobi(
        real_t x[]
      , real_t r[]
      , real_t const b[]
      , size_t const g
      , double const h=1
      , int const echo=0
      , float const omega=.666666
  ) {
      // solve A*x == b on a g grid points
      // using a 2nd order finite-difference stencil: [-1/h^2  2/h^2  -1/h^2]
      double const c0inv = .5*h*h, c0 = 1./c0inv, c1 = -.5*c0;
      double const omega_c0inv = omega*c0inv;

      double norm2{0}, xavg{0}, ravg{0};
      real_t x_prev = x[g - 1]; // init
      for (size_t i = 0; i < g; ++i) {
          real_t const xi = x[i];
          real_t const x_next = x[(i + 1) % g];
          
          double const Ax = c1*x_prev + c0*xi + c1*x_next; // A*x
          double const ri = b[i] - Ax; // residual vector r = b - A*x
          x[i] = xi + omega_c0inv*ri; // new solution, Jacobi update formula
          r[i] = ri; 

          x_prev = xi; // loop-carried dependency!

          // accumulate
          norm2 += ri*ri; // residual norm ||r||_2 belonging to the input solution x (NOT TO THE OUTPUT)
          xavg  += x[i]; // the average of the new solution x
          ravg  += ri; // the average of the residual r
      } // i

      { // scope: subtract average (necessary if all boundary_conditions are periodic)
          xavg /= g; ravg /= g;
          for (size_t i = 0; i < g; ++i) {
              x[i] -= xavg;
              r[i] -= ravg; // in MG methods, the residual vector is used as rhs b on the next coarser level.
                            // since there, the same problem A*x==b is solved with the same boundary condition,
                            // we have to fulfill the physical constraint that b is charge neutral again.
          } // i
      } // scope: subtract average

      return norm2/g;
  } // jacobi

  template <typename real_t>
  status_t smoothen(real_t x[] // preliminary solution
                    , real_t r[] // residual
                    , real_t const b[] // right hand side
                    , size_t const g // number of grid points
                    , double const h // grid spacing
                    , unsigned const maxiter // max. number of iterations
                    , int const echo=0 // log-level
                    , char const *which=""
                    , float const tol=0
                    , int const level=0) { // for display
      unsigned iter{0};
      double rn{1}; // residual norm
      for (iter = 0; iter < maxiter && rn > tol; ++iter) {
          rn = jacobi(x, r, b, g, h);
          if (echo > 9) std::printf("# %s level=%i %ssmoothing step %i (of max %d) residual norm %g\n", __func__, level, which, iter, maxiter, std::log10(rn));
//           if (echo > 29) get_residual_norm(r, x, b, g, h, stdout);
      } // iter
//       rn = get_residual_norm(r, x, b, g, h, (echo > 19)?stdout:0);
      if (echo > 0) std::printf("# %s level=%i %d %ssmoothing steps residual norm %.1f\n", __func__, level, iter, which, std::log10(rn));
      return 0;
  } // smoothen

  template <typename real_t>
  status_t exact_solve(
        real_t x[] // preliminary solution
      , real_t const b[] // right hand side
      , size_t const g // number of grid points
      , double const h // grid spacing
      , int const echo=0 // log-level
  ) {
      if (1 == g) {
          x[0] = 0; // exact solution of a periodic electrostatic potential
      } else if (2 == g) {
          // solve A*x == b on a g grid points using a 2nd order finite-difference
          // stencil: [-1/h^2  2/h^2  -1/h^2] with folded back c1-coefficients:
          // double const c0 = 2./(h*h), c1 = -.5*c0;
          /*        / c0 2*c1 \ / x[0] \    / b[0] \ */
          /* solve  |         | |      | == |      | */
          /*        \ 2*c1 c0 / \ x[1] /    \ b[1] / */
          // under the constraint that x[0] = -x[1] (and b[0] = -b[1]) and knowing that 2*c1 = -c0
          // so :  c0*(x[0] - x[1]) == b[0]
          // and: -c0*(x[0] - x[1]) == b[1]
          double const c0inv = .5*(h*h);
          double const db = b[1] - b[0];
          double const dx = .5*c0inv*db;
          x[0] = -.5*dx; x[1] = .5*dx;
          if (echo > 5) std::printf("# %s h=%g b= %g %g, x= %g %g\n", __func__, h, b[0],b[1], x[0],x[1]);
      } else {
          return 1; // only these two exact solvers implemented!
      }
      return 0;
  } // exact_solve

  
  template <typename real_t>
  status_t mg_cycle(
        real_t x[] // solution
      , real_t const b[] // right hand side
      , unsigned const k // level
      , double const h=1 // grid spacing
      , char const scheme='V' // scheme in {'V','W'}
      , int const echo=0 // log-level
      , short const nu_pre=2  //  pre-smoothing Jacobi steps
      , short const nu_post=2 // post-smoothing Jacobi steps
  ) {
//    std::vector<char> tabs((echo > 0)*2*k + 1, ' '); tabs[(echo > 0)*2*k] = 0;
//    if (echo > 0) std::printf("# %s %s level = %i\n", __func__, tabs.data(), k);
      status_t stat(0);

      size_t const g = 1ul << k; // number of grid points
      std::vector<real_t> r_vec(g, 0.0); // get memory
      auto const r = r_vec.data(); // residual vector

      int const min_level = control::get("multi_grid.test.min.level", 1.); // 1:only 2 grid points

      if (k > min_level) {

          stat += smoothen(x, r, b, g, h, nu_pre, echo, " pre-", 0, k);

          size_t const gc = 1ul << (k - 1);
//        if (echo > 0) std::printf("# %s %s level = %i coarsen from %ld to %ld\n", __func__, tabs.data(), k, g, gc);
          std::vector<real_t> uc(gc, real_t(0)), rc(gc);
          double const hc = 2*h; // coarser grid spacing

          for (char vw = 'V'; vw <= scheme; ++vw) {

              stat += restrict_to_any_grid(rc.data(), gc, r, g, 1, 1, 0, true); // restriction

              stat += mg_cycle(uc.data(), rc.data(), k - 1, hc, scheme, echo, nu_pre, nu_post); // recursive invocation

              stat += linear_interpolation(r, g, uc.data(), gc, 1, 1, 0, true); // prolongation

              add_product(x, g, r, real_t(1));

              stat += smoothen(x, r, b, g, h, nu_post, echo, "post-", 0, k);

          } // vw

      } else { // k > min_level

          // exact solver can be plugged in here
          if (2 == g) {
              stat += exact_solve(x, b, g, h, echo);
          } else {
              stat += smoothen(x, r, b, g, h, 999, echo, "solve-", 1e-21);
          }

      } // k > min_level

//    if (echo > 0) std::printf("# %s %s level = %i status = %i\n", __func__, tabs.data(), k, int(stat));
      return stat;
  } // mg_cycle

  template <typename real_t>
  status_t test_mg_cycle(int const echo=0) {
      status_t stat(0);
      unsigned const k    = control::get("multi_grid.test.levels", 10.);
      short const nu_pre  = control::get("multi_grid.test.pre.jacobi",  2.);
      short const nu_post = control::get("multi_grid.test.post.jacobi", 2.);
      char const scheme  = *control::get("multi_grid.test.scheme", "V"); // {V, W, F}
      char const rhs     = *control::get("multi_grid.test.rhs", "random"); // r:random, s:sin(x*2pi/L)
      unsigned const nit  = control::get("multi_grid.test.iterations", 1.);

      size_t const g = 1ul << k;
      std::vector<real_t> x(g, 0.f), b(g, 0.f);
      double avg{0};
      for (int i = 0; i < g; ++i) {
          if ('s' == (rhs | 32)) {
              b[i] = std::sin(i*2*constants::pi/g);
          } else {
              b[i] = simple_math::random(-1.f, 1.f); // always get the random values in float, so they are the same for float and double versions
          }
          avg += b[i];
      } // i
      avg /= g;
      for (int i = 0; i < g; ++i) { 
          b[i] -= avg; // make b charge neutral
      } // i

      if (echo > 0) std::printf("\n\n\n# %s starting from 2^%i grid points\n", __func__, k);
      double const h = 1./g;
      
      for (int it = 0; it < nit; ++it) {
          stat += mg_cycle(x.data(), b.data(), k, h, scheme, echo, nu_pre, nu_post);
      } // it

      double norm2{0}; double const hm2 = 1./(h*h);
      auto f = (echo > 5) ? ( (echo > 9) ? stdout : std::fopen("multi_grid.out.mg_cycle.dat", "w") ) : nullptr;
      if (f) std::fprintf(f, "## index i, solution x[i], residual r[i], right hand side b[i], Ax[i]   for i < %ld\n", g);
      int const modulo_g = g - 1; // works since g is always a pure power of 2
      for (int i = 0; i < g; ++i) {
          int const i_plus_1  = (i + 1) & modulo_g; // == (i + 1) % g;
          int const i_minus_1 = (i - 1) & modulo_g; // == (i + 1 + g) % g;
          double const Ax = (-x[i_minus_1] + 2*x[i] -x[i_plus_1])*hm2; // simplest 1D finite-difference Laplacian
          double const res = b[i] - Ax;
          norm2 += res*res;
          if (f) std::fprintf(f, "%g %g %g %g %g\n", i*h, x[i], res, b[i], Ax);
      } // i
      if (f && (stdout != f)) std::fclose(f);
      if (echo > 1) std::printf("# %s residual %.3f\n", __func__, std::log10(norm2/g));
      return stat;
  } // test_mg_cycle

  status_t all_tests(int const echo) {
      status_t stat(0);
      int n{0}; int const t = control::get("multi_grid.select.test", -1.); // -1:all
      if (t & (1 << n++)) stat += test_mg_cycle<double>(echo);
      if (t & (1 << n++)) stat += test_transfer(echo);
      if (t & (1 << n++)) stat += test_analysis(echo);
      if (t & (1 << n++)) stat += test_restrict_interpolate(echo);
      if (t & (1 << n++)) stat += test_Morton_indices(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace multi_grid
