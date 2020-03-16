#pragma once

#include "real_space_grid.hxx" // ::grid_t
#include "inline_math.hxx"// set
#include "data_view.hxx"// view2D<T>

typedef int status_t;

namespace iterative_poisson {

  template<typename real_t>
  status_t solve(real_t x[] // result to Laplace(x)/(-4*pi) == b
                , real_t const b[] // right hand side b
                , real_space_grid::grid_t<1> g // grid descriptor
                , int const echo=0 // log level
                , float const threshold=3e-8 // convergence criterion
                , float *residual=nullptr // residual that was reached
                , int const maxiter=999 // maximum number of iterations 
                , int const miniter=3   // minimum number of iterations
                , int restart=4096 // number of iterations before restrat, 1:steepest descent
                );

  inline status_t solve_multiprecision(double x[] // result to Laplace(x)/(-4*pi) == b
                , double const b[] // right hand side b
                , real_space_grid::grid_t<1> g // grid descriptor
                , int const echo=0 // log level
                , float const threshold=3e-8 // convergence criterion
                , float *residual=nullptr // residual that was reached
                , int const maxiter=999 // maximum number of iterations 
                , int const miniter=3   // minimum number of iterations
                , int restart=4096 // number of iterations before restrat, 1:steepest descent
                ) {
      size_t const nall = g.dim(0) * g.dim(1) * g.dim(2);
      view2D<float> xb(2, nall, 0.0); // get memory
      auto const x32 = xb[0], b32 = xb[1];
      set(b32, nall, b); // convert to float
      set(x32, nall, x); // convert to float
      if (echo > 1) printf("# %s solve in <float> precision first\n", __FILE__);
      solve(x32, b32, g, echo + 1, threshold, residual, maxiter >> 4, miniter, restart);
      set(x, nall, x32); // convert to double
      return solve(x, b, g, echo, threshold, residual, maxiter, miniter, restart);
  } // solve_multiprecision
  
  status_t all_tests(int const echo=0);

} // namespace iterative_poisson
