#pragma once

#include "real_space_grid.hxx" // ::grid_t

typedef int status_t;

namespace iterative_poisson {

  template<typename real_t>
  status_t solve(real_t x[] // result to Laplace(x)/(-4*pi) == b
                , real_t const b[] // right hand side b
                , real_space_grid::grid_t<1> g // grid descriptor
                , int const echo=0 // log level
                , float const threshold=3e-8 // convergence criterion
                , float *residual=nullptr // residual that was reached
                , int const maxiter=199 // maximum number of iterations 
                , int const miniter=3  // minimum number of iterations
                , int restart=4096 // number of iterations before restrat, 1:steepest descent
                , char const mixed_precision='m'
                );
  
  status_t all_tests(int const echo=0);

} // namespace iterative_poisson
