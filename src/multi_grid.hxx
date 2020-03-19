#pragma once

#include <cmath> // std::min, std::max
#include <vector> // std::vector<T>

#include "real_space_grid.hxx" // ::grid_t

typedef int status_t;

namespace multi_grid {
  
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
  
  unsigned nearest_binary_power(unsigned const ng) {
      int k{0};
      while(ng > (1ull << (k + 1))) ++k;
      return k;
  } // nearest_binary_power
  
  template <int D0=1>
  inline status_t analyze_grid_sizes(
            real_space_grid::grid_t<D0> const & g // coarse grid where the Kohn-Sham equation is typically solved
          , int const echo=0) {
      status_t stat{0};
      for(int d = 0; d < 3; ++d) {
          unsigned const ng = g.dim(d);
          unsigned const k = nearest_binary_power(ng);
          unsigned const nb = 1ull << k;
          if (echo > 2) printf("# grid in %c-direction can be coarsened from %d to %d = 2^%i grid points\n", 'x'+d, ng, nb, k);
          // compute the restriction operator from ng to n2
          assert(nb < ng);
          double const rowref = 1; // each row sum must be equal to 1
          std::vector<double> colsum(ng, 0.0);
          double const b_over_g = nb/double(ng);
          double const colref = b_over_g; // each column sum must be nb/ng
          double rowdev{0};
          for(int ib = 0; ib < nb; ++ib) {
              if (echo > 8) printf("# %c row%4i  ", 'x'+d, ib);
              double rowsum{0};
              for(int ig = 0; ig < ng; ++ig) {
                  // overlap between a grid compartment of length L/ng starting at ig*L/ng
                  // with a grid compartment of length L/n2 starting at i2*L/n2
                  //            |  0  |  1  |  2  |  3  |  4  |   g
                  //            |    0    |    1    |    2    |   b
                  // (the cell length L falls out of the equation)
                  double const start = std::max((ig - 0)*b_over_g, double(ib - 0));
                  double const end   = std::min((ig + 1)*b_over_g, double(ib + 1));
                  double const ovl = std::max(0.0, end - start);
                  if (ovl > 0) {
                      if (echo > 8) printf("  %i %g", ig, ovl); // show only non-zero matrix entries (with their column index) 
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
              if (echo > 8) printf(" sum=%g dev=%.1e\n", rowsum, rowsum - rowref);
              rowdev = std::max(rowdev, std::abs(rowsum - rowref));
          } // ib
          double coldev{0}; // column deviation
          for(int ig = 0; ig < ng; ++ig) {
              coldev = std::max(coldev, std::abs(colsum[ig] - colref));
          } // ig
          if (echo > 6) printf("# %c-direction largest deviation = %.1e (row) and %.1e (col)\n\n", 'x'+d, rowdev, coldev);
          stat += (coldev > 1e-14) + (rowdev > 1e-14);
      } // d
      
      // in some cases only two coefficients per row are non-zero
      // (in particular the case of power of 2: coefficients are 0.5 0.5)
      // in general, there are up to three non-zero coefficients.
      
      return stat;
  } // analyze_grid_sizes
  

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  inline status_t test_analysis(int const echo=0) {
      real_space_grid::grid_t<1> g(63, 64, 65);
      return analyze_grid_sizes(g, echo);
  } // test_analysis

  inline status_t all_tests(int const echo=0) { 
      status_t stat{0};
      stat += test_analysis(echo);
      return stat; 
  } // all_tests
  
#endif // NO_UNIT_TESTS

} // namespace multi_grid
