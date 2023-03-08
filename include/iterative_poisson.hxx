#pragma once
// This file is part of AngstromCube under MIT License

#include "real_space.hxx" // ::grid_t

#include "status.hxx" // status_t

namespace iterative_poisson {

  template <typename real_t>
  status_t solve(
        real_t x[] // result to Laplace(x)/(-4*pi) == b
      , real_t const b[] // right hand side b
      , real_space::grid_t const & g // grid descriptor
      , char const method='M' // solver method M:multi-grid, c:conjugate-gradient, s:steepest-descent
      , int const echo=0 // log level
      , float const threshold=3e-8 // convergence criterion
      , float *residual=nullptr // residual that was reached
      , int const maxiter=199 // maximum number of iterations 
      , int const miniter=3  // minimum number of iterations
      , int restart=4096 // number of iterations before restart, 1:steepest descent
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace iterative_poisson
