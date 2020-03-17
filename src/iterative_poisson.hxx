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

  /* 
      ToDo: use multi-grid acceleration
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
  
} // namespace iterative_poisson
