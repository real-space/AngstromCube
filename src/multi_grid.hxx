#pragma once

#include <cmath> // std::min, std::max
#include <vector> // std::vector<T>

#include "real_space_grid.hxx" // ::grid_t
#include "data_view.hxx" // view2D<T>
#include "inline_math.hxx" // set, add_product, dot_product
#include "constants.hxx" // ::pi

#include "status.hxx" // status_t

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
  
  
  template <typename real_t, typename real_in_t>
  status_t restrict_to_any_grid(real_t out[], unsigned const go
                      , real_in_t const in[], unsigned const gi
                      , int const stride=1, int const bc=0) {
      if (go < 1) return 0; // early return
      if (gi < 1) return 0; // early return
      
      view2D<real_t> target(out, stride); // wrap
      view2D<real_in_t const> source(in, stride); // wrap
      
      double const ratio = go/double(gi);
      for(int io = 0; io < go; ++io) {
          real_t w8s[4] = {0,0,0,0};
          int        iis[4];
          int nw8 = 0;
          for(int ii = 0; ii < gi; ++ii) {
              double const start = std::max((ii - 0)*ratio, double(io - 0));
              double const end   = std::min((ii + 1)*ratio, double(io + 1));
              double const w8 = std::max(0.0, end - start);
              if (w8 > 0) {
                  assert(nw8 < 4);
                  w8s[nw8] = w8;
                  iis[nw8] = ii;
                  ++nw8;
              } // non-zero
          } // ii
          assert(std::abs((w8s[0] + w8s[1] + w8s[2] + w8s[3]) - 1) < 1e-6);
          set(target[io], stride, 0.0);
          for(int iw8 = 0; iw8 < nw8; ++iw8) {
              int const ii = iis[iw8];
              add_product(target[io], stride, source[ii], w8s[iw8]);
          } // iw8
      } // io
      
  } // restrict_to_any_grid
  
  template <typename real_t, typename real_in_t>
  status_t linear_interpolation(real_t out[], unsigned const go
                      , real_in_t const in[], unsigned const gi
                      , int const stride=1, int const bc=0) {
      if (go < 1) return 0; // early return
      if (gi < 1) return 0; // early return
      
      // linear interpolation between two grids with the same alignment.
      //            |  0  |  1  |  2  |  3  |  4  |    out
      //       -1   |    0    |    1    |    2    |    3    in
      
      view2D<real_t> target(out, stride); // wrap
      view2D<real_in_t const> source(in, stride); // wrap
      view2D<real_in_t> buffer(2, stride, 0.0); // lower and upper halo buffer extending in-array
      if (bc) {
          set(buffer[1], stride, source[0]);      // fill location #3  with element #0 
          set(buffer[0], stride, source[gi - 1]); // fill location #-1 with element #2 
      } // periodic boundary condition
      
      double const ratio = gi/double(go);
      for(int io = 0; io < go; ++io) {
          double const p = (io + 0.5)*ratio + 0.5;
          int const ii = int(p); // index of the upper grid point on the in-array
          double const wu = p - ii; // upper weight onto [ii]
          double const wl = 1 - wu; // lower weight onto [ii - 1]
          printf("# io=%i p=%g ii=%i w8s %g %g\n", io, p, ii, wl, wu); // DEBUG
          assert(std::abs((wl + wu) - 1) < 1e-6);
          if (0 == ii) {
              set(        target[io], stride, buffer[0], wl);
              add_product(target[io], stride, source[ii], wu);
          } else {
              set(        target[io], stride, source[ii - 1], wl);
              add_product(target[io], stride, (ii < gi) ? source[ii] : buffer[1], wu);
          } // ii > 0
      } // io
//       fflush(stdout); exit(42);
  } // linear_interpolation
  

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  inline status_t test_transfer(int const echo=0) {
      if (echo < 1) return 0;
      status_t stat(0);
      int const ng = 32;
      std::vector<double> input_c(ng), result_c(ng), result_cdc(ng);
      std::vector<int> dense_grids{ng, (5*ng)/3, 2*ng}; // numbers of dense grid points
      for(auto mg : dense_grids) {
          std::vector<double> result_d(mg), input_d(mg);
          printf("\n## ik T (%s: interpolate from grid=%d to grid=%d)\n", __func__, ng, mg); // legend
          for(int ik = 1; ik < 3; ++ik) {
              double const ki = (2*constants::pi*ik)/ng;
              for(int ig = 0; ig < ng; ++ig) {
                  input_c[ig] = std::cos(ki*(ig + 0.5));
              } // ig
              
              stat += linear_interpolation(result_d.data(), mg, input_c.data(), ng, 1, 1);
              
              printf("\n## x interpolate(cos(x)) cos(x) [n=%d m=%d k=%i]\n", ng, mg, ik);
              double const kj = (2*constants::pi*ik)/mg;
              for(int jg = 0; jg < mg; ++jg) {
                  input_d[jg] = std::cos(kj*(jg + 0.5));
                  printf("%g %g %g\n", (jg + 0.5)/mg, result_d[jg], input_d[jg]);
              } // jg
              
              stat += restrict_to_any_grid(result_c.data(), ng, input_d.data(), mg, 1, 1);
              stat += restrict_to_any_grid(result_cdc.data(), ng, result_d.data(), mg, 1, 1);

              printf("\n## x cos(x) restrict(interpolate(cos(x))) [n=%d m=%d k=%i]\n", ng, mg, ik);
              for(int ig = 0; ig < ng; ++ig) {
                  printf("%g %g %g %g\n", (ig + 0.5)/ng, result_c[ig], result_cdc[ig], input_c[ig]);
              } // ig
              
          } // ik
      } // gn
      return 0;
  } // test_transfer
  
  inline status_t test_analysis(int const echo=0) {
      real_space_grid::grid_t<1> g(63, 64, 65);
      return analyze_grid_sizes(g, echo);
  } // test_analysis

  inline status_t all_tests(int const echo=0) { 
      status_t stat{0};
//       stat += test_analysis(echo);
      stat += test_transfer(echo);
      return stat;
  } // all_tests
  
#endif // NO_UNIT_TESTS

} // namespace multi_grid
