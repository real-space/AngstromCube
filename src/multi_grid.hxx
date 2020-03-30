#pragma once

#include <cstdio> // printf
#include <cmath> // std::min, std::max
#include <vector> // std::vector<T>

#include <algorithm> // std::minmax_element

#ifndef NO_UNIT_TESTS
  #include <numeric> // std::iota
  #include "constants.hxx" // ::pi
  #include "simple_math.hxx" // ::random
  #include "debug_output.hxx" // dump_to_file
#endif

#include "real_space_grid.hxx" // ::grid_t
#include "data_view.hxx" // view2D<T>
#include "inline_math.hxx" // set, add_product, dot_product
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
  
  template <int D0=1>
  inline status_t analyze_grid_sizes(
            real_space_grid::grid_t<D0> const & g // coarse grid where the Kohn-Sham equation is typically solved
          , int const echo=0) {
      status_t stat{0};
      for(int d = 0; d < 3; ++d) {
          unsigned const ng = g[d];
          unsigned const k = nearest_binary_power(ng);
          unsigned const nb = 1ull << k;
          if (echo > 2) printf("# grid in %c-direction can be coarsened from %d to %d = 2^%i grid points\n", 'x'+d, ng, nb, k);
          assert(nb < ng);
      } // d
      return stat;
  } // analyze_grid_sizes
  
  template <typename real_t, typename real_in_t>
  status_t restrict_to_any_grid(real_t out[], unsigned const go
                      , real_in_t const in[], unsigned const gi
                      , size_t const stride=1, int const bc=0
                      , int const echo=0 // log-level
                      , bool const use_special_version_for_2x=true) {
      if (go < 1) return go; // early return
      if (gi < 1) return gi; // early return
      
      view2D<real_t> target(out, stride); // wrap
      view2D<real_in_t const> source(in, stride); // wrap
      
      if (use_special_version_for_2x && (2*go == gi)) {
          real_t const w8 = 0.5;
          for(int io = 0; io < go; ++io) {
              set(        target[io], stride, source[2*io + 0], w8);
              add_product(target[io], stride, source[2*io + 1], w8);
          } // io
          // printf("\n# %s specialized for m==2*n\n", __func__);

      } else { // use_special_version_for_2x
          assert(go <= gi); 
          if(go == gi) warn("restriction is a copy operation for %d grid points", go);
          
          // any other grid number ratio
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
              set(target[io], stride, real_t(0));
              for(int iw8 = 0; iw8 < nw8; ++iw8) {
                  int const ii = iis[iw8];
                  add_product(target[io], stride, source[ii], w8s[iw8]);
              } // iw8
          } // io
      } // 2x
      
      return 0;
  } // restrict_to_any_grid
  
  template <typename real_t, typename real_in_t>
  status_t linear_interpolation(real_t out[], unsigned const go
                      , real_in_t const in[], unsigned const gi
                      , size_t const stride=1, int const bc=0
                      , int const echo=0 // log-level
                      , bool const use_special_version_for_2x=true) {
      if (go < 1) return go; // early return
      if (gi < 1) return gi; // early return

      view2D<real_t> target(out, stride); // wrap
      view2D<real_in_t const> source(in, stride); // wrap
      view2D<real_in_t> buffer(2, stride, 0.0); // lower and upper halo buffer extending in-array
      if (bc) {
          set(buffer[1], stride, source[0]);      // fill location #3  with element #0 
          set(buffer[0], stride, source[gi - 1]); // fill location #-1 with element #2 
      } // periodic boundary condition

      if (echo > 2) printf("# %s from %d to %d grid points (inner dim %ld)\n", __func__, gi, go, stride);
      
      if ((go == 2*gi) && use_special_version_for_2x) {
          real_t const w14 = 0.25, w34 = 0.75;
          // linear interpolation from a grid to a 2x denser grid
          //            |  0  |  1  |  2  |  3  |  4  |  5  |         out
          //       -1   |     0     |     1     |     2     |     3    in
          set(        target[0], stride, buffer[0], w14); // first dense grid point
          add_product(target[0], stride, source[0], w34); // first dense grid point
          for(int ii = 1; ii < gi; ++ii) {
              set(        target[2*ii - 1], stride, source[ii - 1], w34); //  odd dense grid point
              add_product(target[2*ii - 1], stride, source[ii + 0], w14); //  odd dense grid point
              set(        target[2*ii + 0], stride, source[ii - 1], w14); // even dense grid point
              add_product(target[2*ii + 0], stride, source[ii + 0], w34); // even dense grid point
          } // ii
          set(        target[go - 1], stride, source[gi - 1], w34); // last dense grid point
          add_product(target[go - 1], stride, buffer[1],      w14); // last dense grid point
          // printf("\n# %s specialized for m==2*n\n", __func__);
        
      } else { // use_special_version_for_2x
        
          // any other grid number ratio
          double const ratio = gi/double(go);
          for(int io = 0; io < go; ++io) {
              // linear interpolation between two grids with the same alignment.
              //            |  0  |  1  |  2  |  3  |  4  |    out
              //       -1   |    0    |    1    |    2    |    3    in
              double const p = (io + 0.5)*ratio + 0.5;
              int const ii = int(p); // index of the upper grid point on the in-array
              real_t const wu = p - ii; // upper weight onto [ii]
              real_t const wl = 1 - wu; // lower weight onto [ii - 1]
    //           printf("# io=%i p=%g ii=%i w8s %g %g\n", io, p, ii, wl, wu); // DEBUG
              assert(std::abs((wl + wu) - 1) < 1e-6);
              set(        target[io], stride, (ii >  0) ? source[ii - 1] : buffer[0], wl);
              add_product(target[io], stride, (ii < gi) ? source[ii]     : buffer[1], wu);
          } // io
          
      } // 2x
      return 0;
  } // linear_interpolation
  
  template <typename real_t>
  void print_min_max(real_t const *const begin, real_t const *const end
        , int const echo=0, char const *title="<array>", char const *unit="") {
#ifdef DEVEL
      if (echo < 1) return;
      auto const mm = std::minmax_element(begin, end);
      printf("# %24s in range [%g, %g] %s\n", title, *mm.first, *mm.second, unit);
#endif
  } // print_min_max

  // now 3D functions:
  template <typename real_t, typename real_in_t=real_t, int D0=1>
  status_t restrict3D(real_t out[], real_space_grid::grid_t<D0> const & go
            , real_in_t const in[], real_space_grid::grid_t<D0> const & gi
            , int const echo=0) { // log-level
      status_t stat(0);

      { // scope: checks
          int bc[3];
          for(char d = 'x'; d  <= 'z'; ++d) {
              assert(go(d) < gi(d));
              bc[d-120] = gi.boundary_condition(d);
              assert(go.boundary_condition(d) == bc[d-120]);
              if (echo > 0) printf("# %s in %c-direction from %d to %d grid points, boundary condition %d\n", 
                                      __func__, d, gi(d), go(d), bc[d-120]);
          } // d
      } // scope

      print_min_max(in, in + gi.all(), echo, "input");

      size_t const nyx = gi('y')*gi('x');
      view2D<real_t> tz(go('z'), nyx); // get memory
      stat += restrict_to_any_grid(tz.data(), go('z'), 
                                          in, gi('z'), nyx, gi.boundary_condition('z'), echo);
      
      print_min_max(tz.data(), tz.data() + go('z')*nyx, echo, "z-restricted");

      size_t const nx = gi('x');
      view2D<real_t> ty(go('z'), go('y')*nx); // get memory
      for(int z = 0; z < go('z'); ++z) {
          stat += restrict_to_any_grid(ty[z], go('y'),
                                       tz[z], gi('y'), nx, gi.boundary_condition('y'), echo*(z == 0));
      } // z
      
      print_min_max(ty.data(), ty.data() + go('z')*go('y')*nx, echo, "zy-restricted");

      view3D<real_t> tx(out, go('y'), go('x')); // wrap
      for(int z = 0; z < go('z'); ++z) {
          for(int y = 0; y < go('y'); ++y) {
              stat += restrict_to_any_grid(tx(z,y), go('x'),
                                      ty[z] + y*nx, gi('x'), 1, gi.boundary_condition('x'), echo*(z == 0)*(y == 0));
          } // y
      } // z

      print_min_max(out, out + go.all(), echo, "output");

      return stat;
  } // restrict3D

  template <typename real_t, typename real_in_t=real_t, int D0=1>
  status_t interpolate3D(real_t out[], real_space_grid::grid_t<D0> const & go
               , real_in_t const in[], real_space_grid::grid_t<D0> const & gi, int const echo=0) {
      status_t stat(0);
      for(int d = 0; d < 3; ++d) {
          assert(go.boundary_condition(d) == gi.boundary_condition(d));
      } // d

      view3D<real_in_t const> ti(in, gi('y'), gi('x')); // wrap
      view3D<real_t>     tx(gi('z'), gi('y'), go('x')); // get memory
      for(int z = 0; z < gi('z'); ++z) {
          for(int y = 0; y < gi('y'); ++y) {
              stat += linear_interpolation(tx(z,y), go('x'),
                                           ti(z,y), gi('x'), 1, gi.boundary_condition('x'), echo*(z == 0)*(y == 0));
          } // y
      } // z

      size_t const nx = go('x');
      view2D<real_t> ty(go('z'), go('y')*nx); // get memory
      for(int z = 0; z < go('z'); ++z) {
          stat += linear_interpolation(ty[z], go('y'),
                                     tx(z,0), gi('y'), nx, gi.boundary_condition('y'), echo*(z == 0));
      } // z

      size_t const nyx = go('y')*go('x');
      stat += linear_interpolation(out, go('z'),
                             ty.data(), gi('z'), nyx, gi.boundary_condition('z'), echo);
      
      return stat;
  } // interpolate3D
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  inline status_t check_general_restrict(unsigned const ng, int const echo=0, char const dir='?') {
      unsigned const k = nearest_binary_power(ng);
      unsigned const nb = 1ull << k;
      if (echo > 2) printf("# grid in %c-direction can be coarsened from %d to %d = 2^%i grid points\n", dir, ng, nb, k);
      // compute the restriction operator from ng to n2
      assert(nb < ng);
      double const rowref = 1; // each row sum must be equal to 1
      std::vector<double> colsum(ng, 0.0);
      double const b_over_g = nb/double(ng);
      double const colref = b_over_g; // each column sum must be nb/ng
      double rowdev{0};
      for(int ib = 0; ib < nb; ++ib) {
          if (echo > 8) printf("# %c row%4i  ", dir, ib);
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
      if (echo > 6) printf("# %c-direction largest deviation = %.1e (row) and %.1e (col)\n\n", dir, rowdev, coldev);
      return (coldev > 1e-14) + (rowdev > 1e-14);
      
      // in some cases only two coefficients per row are non-zero
      // (in particular the case of power of 2: coefficients are 0.5 0.5)
      // in general, there are up to three non-zero coefficients.
      
  } // check_general_restrict


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
              double const k = (2*constants::pi*ik);
              for(int ig = 0; ig < ng; ++ig) {
                  double const x = (ig + 0.5)/ng;
                  input_c[ig] = std::cos(k*x);
              } // ig
              
              stat += linear_interpolation(result_d.data(), mg, input_c.data(), ng, 1, 1);
              
              printf("\n## x interpolate(cos(x)) cos(x) [n=%d m=%d k=%i]\n", ng, mg, ik);
              for(int jg = 0; jg < mg; ++jg) {
                  double const x = (jg + 0.5)/mg;
                  input_d[jg] = std::cos(k*x);
                  printf("%g %g %g\n", x, result_d[jg], input_d[jg]);
              } // jg
              
              stat += restrict_to_any_grid(result_c.data(), ng, input_d.data(), mg, 1, 1);
              stat += restrict_to_any_grid(result_cdc.data(), ng, result_d.data(), mg, 1, 1);

              printf("\n## x cos(x) restrict(interpolate(cos(x))) [n=%d m=%d k=%i]\n", ng, mg, ik);
              for(int ig = 0; ig < ng; ++ig) {
                  double const x = (ig + 0.5)/ng;
                  printf("%g %g %g %g\n", x, result_c[ig], result_cdc[ig], input_c[ig]);
              } // ig
              
          } // ik
      } // mg
      return stat;
  } // test_transfer
  
  inline status_t test_analysis(int const echo=0) {
      status_t stat(0);
      real_space_grid::grid_t<1> g(63, 64, 65);
      for(char dir = 'x'; dir <= 'z'; ++dir) {
          stat += check_general_restrict(g(dir), echo, dir);
      } // direction
      return stat + analyze_grid_sizes(g, echo);
  } // test_analysis
  
  inline status_t test_restrict_interpolate(int const echo=0) {
      status_t stat(0);
      real_space_grid::grid_t<1> gi(15, 16, 17), go(8, 8, 16);
      std::vector<float> in(gi.all()), out(go.all());
      std::iota(in.begin(), in.end(), .5f);
      stat += restrict3D(out.data(), go, in.data(), gi);
      stat += interpolate3D(in.data(), gi, out.data(), go);
      return stat;
  } // test_restrict_interpolate

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // toy model for testing the V_cycle structure, operator A = stencil{1, -2, 1}
  template <typename real_t>
  inline double jacobi(real_t x[], real_t r[], real_t const b[], unsigned const k, double const h) {
      // solve A*x == b on a 2^k-grid
      size_t const g = 1ul << k;
      double const c0inv = -.5*h*h, c0 = 1./c0inv, c1 = -.5*c0;
      
      double norm2{0};
      real_t x_prev = x[g - 1]; // init
      for(int i = 0; i < g; ++i) {
          real_t const xi = x[i];
          real_t const x_next = x[(i + 1) % g];
          
          real_t const Lx = b[i] - c1*x_prev - c1*x_next;
          x[i] = c0inv*Lx; // new solution, Jacobi update formula
          r[i] = c0*xi - Lx; // residual vector r = A*x - b
          
          norm2 += r[i]*r[i]; // accumulate residual norm ||r||_2
          x_prev = xi; // loop-carried dependency!
      } // i
      
      return norm2/g;
  } // jacobi

  template <typename real_t>
  inline status_t V_cycle(real_t x[], real_t const b[], unsigned const k, int const echo=0, double const h=1
                         , short const nu1=7, short const nu2=7) {
      std::vector<char> tabs((echo > 0)*4*k + 1, ' '); tabs[(echo > 0)*4*k] = 0;
      if (echo > 0) printf("# %s %s level = %i\n", __func__, tabs.data(), k);
      status_t stat(0);

      size_t const g = 1ul << k; // number of grid points
      std::vector<real_t> r_vec(g); // get memory
      auto const r = r_vec.data(); // residual vector

      {   double rn{-1};
          for(int i = 0; i < nu1; ++i) {
              rn = jacobi(x, r, b, k, h);
              if (echo > 9) printf("# %s %s level = %i  pre-smoothing step %i residual norm %g\n", __func__, tabs.data(), k, i, rn);
          } // pre-smoothing
          if (echo > 0) printf("# %s %s level = %i after %i  pre-smoothing steps residual norm %.2e\n", __func__, tabs.data(), k, nu1, rn);
      }
      
      if (k > 1) {
          size_t const gc = 1ul << (k - 1);
          if (echo > 0) printf("# %s %s level = %i coarsen from %ld to %ld\n", __func__, tabs.data(), k, g, gc);
          std::vector<real_t> uc(gc, 0.f), rc(gc);
          double const hc = 2*h; // coarser grid spacing

          stat += restrict_to_any_grid(rc.data(), gc, r, g, 1, 1, 0, true);

          stat += V_cycle(uc.data(), rc.data(), k-1, echo, hc); // recursive invocation
          
          stat += linear_interpolation(r, g, uc.data(), gc, 1, 1, 0, true);
          
          for(int i = 0; i < g; ++i) x[i] -= r[i];
      } // coarsen
      
      { // scope: subtract average (necessary if all boundary_conditions are periodic)
          double avg{0}; for(int i = 0; i < g; ++i) { avg += x[i]; }
          avg /= g;      for(int i = 0; i < g; ++i) { x[i] -= avg; }
      } // scope: subtract average

      {   double rn{-1};
          for(int i = 0; i < nu2; ++i) {
              rn = jacobi(x, r, b, k, h);
              if (echo > 9) printf("# %s %s level = %i post-smoothing step %i residual norm %g\n", __func__, tabs.data(), k, i, rn);
          } // post-smoothing
          if (echo > 0) printf("# %s %s level = %i after %i post-smoothing steps residual norm %.2e\n", __func__, tabs.data(), k, nu2, rn);
      }

      if (echo > 0) printf("# %s %s level = %i status = %i\n", __func__, tabs.data(), k, int(stat));
      return stat;
  } // V_cycle
  
  template <typename real_t>
  inline status_t test_V_cycle(int const echo=0) {
      unsigned const k = 12;
      size_t const g = 1ul << k;
      std::vector<real_t> x(g, 0.f), b(g, 0.f);
      double avg{0};
      for(int i = 0; i < g; ++i) {
          b[i] = simple_math::random(-1.f, 1.f); // always get the random values in float, so they are the same for float and double versions
          avg += b[i];
      } // i
      avg /= g; for(int i = 0; i < g; ++i) { b[i] -= avg; } // make b charge neutral

      if (echo > 0) printf("\n\n\n# %s starting from 2^%i grid points\n", __func__, k);
      auto const stat = V_cycle(x.data(), b.data(), k, echo);

      double norm2{0};
      FILE* f = (echo > 9) ? fopen("multi_grid.out.V_cycle.dat", "w") : nullptr;
      if (f) fprintf(f, "## index i, solution x[i], Ax[i], residual r[i], right hand side b[i] for i < %ld\n", g);
      for(int i = 0; i < g; ++i) {
          auto const Ax = x[(i - 1) % g] + x[(i + 1) % g] - 2*x[i]; // simplest 1D finite-difference Laplacian
          auto const res = Ax - b[i];
          norm2 += res*res;
          if (f) fprintf(f, "%i %g %g %g %g\n", i, x[i], Ax, res, b[i]);
      } // i
      if (f) fclose(f);
      if (echo > 1) printf("# %s residual %.2e\n", __func__, norm2/g);
      return stat;
  } // test_V_cycle
  
  inline status_t all_tests(int const echo=0) { 
      status_t stat{0};
//       stat += test_transfer(echo);
//       stat += test_analysis(echo);
//       stat += test_restrict_interpolate(echo);
      stat += test_V_cycle<double>(echo);
      return stat;
  } // all_tests
  
#endif // NO_UNIT_TESTS

} // namespace multi_grid
