#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cmath> // std::abs
#include <cstdint> // uint32_t
#include <algorithm> // std::min, ::max, ::minmax_element

#include "real_space.hxx" // ::grid_t
#include "data_view.hxx" // view2D<T>
#include "inline_math.hxx" // set, add_product
#include "status.hxx" // status_t

namespace multi_grid {

  status_t analyze_grid_sizes(
        real_space::grid_t const & g // coarse grid where the Kohn-Sham equation is typically solved
      , uint32_t *n_coarse
      , int const echo=0
  ); // declaration only

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

  template <typename real_t, typename real_in_t>
  status_t restrict_to_any_grid(
        real_t out[]
      , unsigned const go
      , real_in_t const in[]
      , unsigned const gi
      , size_t const stride=1
      , int const bc=0
      , int const echo=0 // log-level
      , bool const use_special_version_for_2x=true
  ) {
      if (go < 1) return go; // early return
      if (gi < 1) return gi; // early return
      
      view2D<real_t> target(out, stride); // wrap
      view2D<real_in_t const> source(in, stride); // wrap
      
      if (use_special_version_for_2x && (2*go == gi)) {
          real_t const w8 = 0.5;
          for (int io = 0; io < go; ++io) {
              set(        target[io], stride, source[2*io + 0], w8);
              add_product(target[io], stride, source[2*io + 1], w8);
          } // io
          // std::printf("\n# %s specialized for m==2*n\n", __func__);

      } else { // use_special_version_for_2x
          assert(go <= gi); 
          if (go == gi) warn("restriction is a copy operation for %d grid points", go);

          // any other grid number ratio
          double const ratio = go/double(gi);
          for (int io = 0; io < go; ++io) {
              real_t w8s[4] = {0,0,0,0};
              int iis[4];
              int nw8 = 0;
              for (int ii = 0; ii < gi; ++ii) {
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
              for (int iw8 = 0; iw8 < nw8; ++iw8) {
                  int const ii = iis[iw8];
                  add_product(target[io], stride, source[ii], w8s[iw8]);
              } // iw8
          } // io
      } // 2x
      
      return 0;
  } // restrict_to_any_grid


  template <typename real_t, typename real_in_t>
  status_t linear_interpolation(
        real_t out[]
      , unsigned const go
      , real_in_t const in[]
      , unsigned const gi
      , size_t const stride=1
      , int const periodic=0
      , int const echo=0 // log-level
      , bool const use_special_version_for_2x=true
  ) {
      if (go < 1) return go; // early return
      if (gi < 1) return gi; // early return

      view2D<real_t> target(out, stride); // wrap
      view2D<real_in_t const> source(in, stride); // wrap: source(a,b) --> in[a*stride + b]
      view2D<real_in_t> buffer(2, stride, 0.0); // get memory lower and upper halo buffer extending in-array
      if (periodic) {
          set(buffer[1], stride, source[0]);      // fill location #3  with element #0 
          set(buffer[0], stride, source[gi - 1]); // fill location #-1 with element #2 
      } // periodic boundary condition

      if (echo > 2) std::printf("# %s from %d to %d grid points (inner dim %ld)\n", __func__, gi, go, stride);
      
      if ((go == 2*gi) && use_special_version_for_2x) {
          real_t const w14 = 0.25, w34 = 0.75;
          // linear interpolation from a grid to a 2x denser grid
          //            |  0  |  1  |  2  |  3  |  4  |  5  |         out
          //       -1   |     0     |     1     |     2     |     3    in
          set(        target[0], stride, buffer[0], w14); // first dense grid point
          add_product(target[0], stride, source[0], w34); // first dense grid point
          for (int ii = 1; ii < gi; ++ii) {
              set(        target[2*ii - 1], stride, source[ii - 1], w34); //  odd dense grid point
              add_product(target[2*ii - 1], stride, source[ii + 0], w14); //  odd dense grid point
              set(        target[2*ii + 0], stride, source[ii - 1], w14); // even dense grid point
              add_product(target[2*ii + 0], stride, source[ii + 0], w34); // even dense grid point
          } // ii
          set(        target[go - 1], stride, source[gi - 1], w34); // last dense grid point
          add_product(target[go - 1], stride, buffer[1],      w14); // last dense grid point
          // std::printf("\n# %s specialized for m==2*n\n", __func__);
        
      } else { // use_special_version_for_2x
        
          // any other grid number ratio
          double const ratio = gi/double(go);
          for (int io = 0; io < go; ++io) {
              // linear interpolation between two grids with the same alignment.
              //            |  0  |  1  |  2  |  3  |  4  |    out
              //       -1   |    0    |    1    |    2    |    3    in
              double const p = (io + 0.5)*ratio + 0.5;
              int const ii = int(p); // index of the upper grid point on the in-array
              real_t const wu = p - ii; // upper weight onto [ii]
              real_t const wl = 1 - wu; // lower weight onto [ii - 1]
    //           std::printf("# io=%i p=%g ii=%i w8s %g %g\n", io, p, ii, wl, wu); // DEBUG
              assert(std::abs((wl + wu) - 1) < 1e-6);
              set(        target[io], stride, (ii >  0) ? source[ii - 1] : buffer[0], wl);
              add_product(target[io], stride, (ii < gi) ? source[ii]     : buffer[1], wu);
          } // io
          
      } // 2x
      return 0;
  } // linear_interpolation


  template <typename real_t>
  void print_min_max(
        real_t const *const begin
      , real_t const *const end
      , int const echo=0
      , char const *title="<array>"
      , char const *unit=""
  ) {
#ifdef DEVEL
      if (echo < 1) return;
      auto const mm = std::minmax_element(begin, end);
      std::printf("# %24s in range [%g, %g] %s\n", title, *mm.first, *mm.second, unit);
#endif
  } // print_min_max


  // now 3D functions:
  template <typename real_t, typename real_in_t=real_t>
  status_t restrict3D(
        real_t out[]
      , real_space::grid_t const & go
      , real_in_t const in[]
      , real_space::grid_t const & gi
      , int const echo=0 // log-level
  ) {
      status_t stat(0);

      { // scope: checks
          int bc[3];
          for (char d = 'x'; d  <= 'z'; ++d) {
              assert(go(d) <= gi(d));
              bc[d-120] = gi.boundary_condition(d);
              assert(go.boundary_condition(d) == bc[d-120]);
              if (echo > 0) std::printf("# %s in %c-direction from %d to %d grid points, boundary condition %d\n", 
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
      for (int z = 0; z < go('z'); ++z) {
          stat += restrict_to_any_grid(ty[z], go('y'),
                                       tz[z], gi('y'), nx, gi.boundary_condition('y'), echo*(z == 0));
      } // z
      
      print_min_max(ty.data(), ty.data() + go('z')*go('y')*nx, echo, "zy-restricted");

      view3D<real_t> tx(out, go('y'), go('x')); // wrap
      for (int z = 0; z < go('z'); ++z) {
          for (int y = 0; y < go('y'); ++y) {
              stat += restrict_to_any_grid(tx(z,y), go('x'),
                                      ty[z] + y*nx, gi('x'), 1, gi.boundary_condition('x'), echo*(z == 0)*(y == 0));
          } // y
      } // z

      print_min_max(out, out + go.all(), echo, "output");

      return stat;
  } // restrict3D


  template <typename real_t, typename real_in_t=real_t>
  status_t interpolate3D(
        real_t out[]
      , real_space::grid_t const & go
      , real_in_t const in[]
      , real_space::grid_t const & gi
      , int const echo=0 // log-level
  ) {
      status_t stat(0);
      for (int d = 0; d < 3; ++d) {
          assert(go.boundary_condition(d) == gi.boundary_condition(d));
      } // d

      view3D<real_in_t const> ti(in, gi('y'), gi('x')); // wrap
      view3D<real_t>     tx(gi('z'), gi('y'), go('x')); // get memory
      for (int z = 0; z < gi('z'); ++z) {
          for (int y = 0; y < gi('y'); ++y) {
              stat += linear_interpolation(tx(z,y), go('x'),
                                           ti(z,y), gi('x'), 1, gi.boundary_condition('x'), echo*(z == 0)*(y == 0));
          } // y
      } // z

      size_t const nx = go('x');
      view2D<real_t> ty(gi('z'), go('y')*nx); // get memory
      for (int z = 0; z < gi('z'); ++z) {
          stat += linear_interpolation(ty[z], go('y'),
                                     tx(z,0), gi('y'), nx, gi.boundary_condition('y'), echo*(z == 0));
      } // z

      size_t const nyx = go('y')*go('x');
      stat += linear_interpolation(out, go('z'),
                             ty.data(), gi('z'), nyx, gi.boundary_condition('z'), echo);
      
      return stat;
  } // interpolate3D


  status_t all_tests(int const echo=0); // declaration only

} // namespace multi_grid
