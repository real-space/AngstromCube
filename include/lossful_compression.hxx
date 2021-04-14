#pragma once

#include <cassert> // assert
#include <vector> // std::vector<T>
#include <cstdio> // std::fprintf, stdout, FILE

  template <typename real_t, typename y_real_t>
  int RamerDouglasPeucker(
      std::vector<bool> & active
    , real_t const x_list[]
    , y_real_t const y_list[]
    , int const end
    , int const begin=0
    , float const epsilon=1e-6
  ) {
    // Find the point with the maximum distance
    double d2max{0};
    int index{-1};
    // construct a straight line through the 1st and last point
    int const p0 = begin, pl = end - 1;
    assert(active[p0]); // start and end point:
    assert(active[pl]); // both points still have to be active
    double const x0 = x_list[p0], y0 = y_list[p0];
    double const x_dist = x_list[pl] - x0;
    double const y_dist = y_list[pl] - y0;
    double const det = x_dist*x_dist + y_dist*y_dist;
    if (det <= 0) return -1; // failed
    double const det_inv = 1./det;

    for(int i = p0 + 1; i < pl; ++i) { // loop over all points in between
        if (active[i]) {
            double const xi = x_list[i], yi = y_list[i];
            /*
            // solve
            //    p0x + s0*x_dist == pix + si*y_dist == 0
            //    p0y + s0*y_dist == piy - si*x_dist == 0
            // for s0, si:
            //    / x_dist   -y_dist \   / s0 \     / pix - p0x \
            //    |                  | * |    |  == |           |
            //    \ y_dist    x_dist /   \ si /     \ piy - p0y /
            //
            // solution: adjoint matrix divided by determinant
            //    / s0 \     /  x_dist   y_dist \   / pix - p0x \     1
            //    |    |  == |                  | * |           |  * ----
            //    \ si /     \ -y_dist   x_dist /   \ piy - p0y /    det
            */
//             double const s0 = (x_dist*(xi - x0) + y_dist*(yi - y0))*det_inv;
//             double const xf0 = x0 + s0*x_dist;
//             double const yf0 = y0 + s0*y_dist;

            double const si = (x_dist*(yi - y0) - y_dist*(xi - x0))*det_inv;
//             double const xf = xi + si*y_dist;
//             double const yf = yi - si*x_dist;
//             assert( pow2(xf - xf0) + pow2(yf - yf0) < 1e-12);

//             double const d2 = pow2(yf - yi) + pow2(xf - xi); // perpendicularDistance^2
            double const d2 = si*si*det;
    //         printf("# DouglasPeucker distance^2 = %g at point #%d\n", d2, i);
            if (d2 > d2max) {
                index = i;
                d2max = d2;
            } // find the maximum
        } // active?
    } // i

    // If max distance is greater than epsilon, recursively simplify
    if (d2max > epsilon*epsilon) {
        assert(index > p0); assert(index < pl);
//      printf("#      DouglasPeucker case N, distance^2 = %g at point #%d\n", d2max, index);
        // Recursive call
//      printf("# call DouglasPeucker on interval [%d, %d] and  [%d, %d] later \n", begin, index, index, end - 1);
        int const n0 = RamerDouglasPeucker(active, x_list, y_list, index + 1, begin, epsilon);
        int const n1 = RamerDouglasPeucker(active, x_list, y_list, end,       index, epsilon);
        return n0 + n1 - 1;
    } else {
        int n_off{0};
        for(int i = p0 + 1; i < pl; ++i) {
            n_off += (0 != active[i]);
            active[i] = false; // switch off every point in between since the curve is sufficiently smooth there
        } // i
//      if (n_off > 0) printf("# DouglasPeucker eleminate %d points, largest distance^2 is %g\n", n_off, d2max);
        return 2;
    } // if

  } // RamerDouglasPeucker

  template <typename real_t, typename y_real_t>
  std::vector<bool> RDP_lossful_compression(
        real_t const x[]
      , y_real_t const y[]
      , int const n
      , float const epsilon=1e-6
  ) {
      std::vector<bool> mask(n, true);
      RamerDouglasPeucker(mask, x, y, n, 0, epsilon);
      return mask;
  } // RDP_lossful_compression

  template <typename real_t, typename real_y_t>
  void print_compressed(
        real_t const x[]
      , real_y_t const y[]
      , int const n
      , float const epsilon=1e-6
      , FILE* os=stdout
  ) {
      auto const mask = RDP_lossful_compression(x, y, n, epsilon);
      for(int i = 0; i < n; ++i) {
          if (mask[i]) std::fprintf(os, "%g %g\n", x[i], y[i]);
      } // i
      std::fprintf(os, "\n");
  } // print_compressed
