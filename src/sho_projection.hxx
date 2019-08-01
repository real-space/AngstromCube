#pragma once

#include <cstdint> // uint32_t
#include <cstdio> // printf

typedef int status_t;

#include "sho_tools.hxx" // sho_tools::nSHO, sho_tools::get_nu, sho_tools::order_zyx
#include "real_space_grid.hxx" // real_space_grid::grid_t<real_t,D0>
#include "hermite_polynomial.hxx" // hermite_polys

namespace sho_projection {

  template<typename real_t, int D0, int PROJECT0_OR_ADD1>
  status_t sho_project_or_add(real_t coeff[] // result if projecting, coefficients are zyx-ordered
                     , int const numax // how many
                     , double const center[3] // where
                     , double const sigma
                     , real_t values[] // grid array, result if adding
                     , real_space_grid::grid_t<real_t,D0> const &g // grid descriptor, assume that g is a Cartesian grid
                     , int const echo=4
                      ) { //
      assert( 1 == D0 );
      double const rcut = 9*sigma;
      double const sigma_inv = 1./sigma;
      // determine the limitations of the projection domain
      int off[3], end[3], num[3];
      for(int dir = 0; dir < 3; ++dir) {
          off[dir] = std::ceil((center[dir] - rcut)*g.inv_h[dir]);
          end[dir] = std::ceil((center[dir] + rcut)*g.inv_h[dir]);
          if (echo > 3) printf("# prelim for %c-direction are [%d, %d)\n", 120+dir, off[dir], end[dir]);
          off[dir] = std::max(off[dir], 0); // lower
          end[dir] = std::min(end[dir], g.dim(dir)); // upper boundary
          if (echo > 3) printf("# limits for %c-direction are [%d, %d)\n", 120+dir, off[dir], end[dir]);
          num[dir] = std::max(0, end[dir] - off[dir]);
      } // dir
      if (echo > 2) printf("# rectangular sub-domain z:[%d, %d) y:[%d, %d) x:[%d, %d)\n", 
                           off[2], end[2], off[1], end[1], off[0], end[0]);
      long const nvolume = num[0] * num[1] * num[2];
      if (echo > 1) printf("# %s on rectangular sub-domain %d * %d * %d = %ld points\n", 
                 (0 == PROJECT0_OR_ADD1)?"project":"add", num[2], num[1], num[0], nvolume);
      if (nvolume < 1) return 0; // no range

      // ToDo: analyze if the grid spacing is small enough for this \sigma

      int const M = 1 + numax;
      real_t* H1d[3];
      for(int dir = 0; dir < 3; ++dir) {
          int const nd = num[dir];
          H1d[dir] = new real_t[nd*M]; // get memory

          real_t const grid_spacing = g.h[dir];
          if (echo > 5) printf("\n# Hermite polynomials for %c-direction:\n", 120+dir);
          for(int ii = 0; ii < nd; ++ii) {
              int const ix = ii + off[dir]; // offset
              real_t const x = (ix*grid_spacing - center[dir])*sigma_inv;
              hermite_polys(H1d[dir] + ii*M, x, numax);
              if (echo > 5) { printf("%g\t", x); for(int nu = 0; nu <= numax; ++nu) printf("%12.6f", H1d[dir][ii*M + nu]); printf("\n"); }
          } // i
      } // dir
   
      int const nSHO = sho_tools::nSHO(numax);
      if (0 == PROJECT0_OR_ADD1) set(coeff, nSHO, (real_t)0);
      
      // int ixyz = 0;
      // for(        int iz = 0; iz < g.dim('z'); ++iz) {
      //     for(    int iy = 0; iy < g.dim('y'); ++iy) {
      //         for(int ix = 0; ix < g.dim('x'); ++ix) {
      for(        int iz = off[2]; iz < end[2]; ++iz) {
          for(    int iy = off[1]; iy < end[1]; ++iy) {
              for(int ix = off[0]; ix < end[0]; ++ix) {
                  int const ixyz = (iz*g.dim('y') + iy)*g.dim('x') + ix;
                  auto const val0 = values[ixyz];
                  auto val = val0;
                  if ((0 != val) || (1 == PROJECT0_OR_ADD1))  {
//                    if (echo > 6) printf("%g %g\n", std::sqrt(vz*vz + vy*vy + vx*vx), val); // plot function value vs r
                      int iSHO = 0;
                      for(int nz = 0; nz <= numax; ++nz) {
                          for(int ny = 0; ny <= numax - nz; ++ny) {
                              for(int nx = 0; nx <= numax - nz - ny; ++nx) {
                                  auto const H3d = H1d[2][(iz - off[2])*M + nz]
                                                 * H1d[1][(iy - off[1])*M + ny]
                                                 * H1d[0][(ix - off[0])*M + nx];
                                  if (1 == PROJECT0_OR_ADD1) {
                                      val += coeff[iSHO] * H3d; // here, the addition happens                                          
                                  } else {
                                      coeff[iSHO] += val0 * H3d; // here, the projection happens                                          
                                  }
                                  ++iSHO;
                              } // nx
                          } // ny
                      } // nz
                      assert( nSHO == iSHO );
                  } // non-zero
                  if (1 == PROJECT0_OR_ADD1) values[ixyz] = val;
                  // ++ixyz;
              } // ix
          } // iy
      } // iz

      if (0 == PROJECT0_OR_ADD1) scale(coeff, nSHO, (real_t)g.dV()); // volume element of the grid

      for(int dir = 0; dir < 3; ++dir) {
          delete[] H1d[dir]; // free memory
      } // dir
      return 0; // success
  } // sho_project_or_add
  

  // wrapper function
  template<typename real_t, int D0>
  status_t sho_project(real_t coeff[] // result, coefficients are zyx-ordered
                     , int const numax // how many
                     , double const center[3] // where
                     , double const sigma
                     , real_space_grid::grid_t<real_t,D0> const &g // grid descriptor, assume that g is a Cartesian grid
                     , real_t const *values=nullptr // input, grid array, defaults to g.values
                     , int const echo=4) { //
      real_t const *val = (nullptr != values)? values : g.values;
      return sho_project_or_add<real_t,D0,0>(coeff, numax, center, sigma, (real_t*)val, g, echo);
  } // sho_project

  // wrapper function
  template<typename real_t, int D0>
  status_t sho_add(real_t values[] // result gets modified, grid array
                 , real_space_grid::grid_t<real_t,D0> const &g // grid descriptor, assume that g is a Cartesian grid
                 , real_t const coeff[] // input, coefficients are zyx-ordered
                 , int const numax // how many
                 , double const center[3] // where
                 , double const sigma
                 , int const echo=4) { //
      return sho_project_or_add<real_t,D0,1>((real_t*)coeff, numax, center, sigma, values, g, echo);
  } // sho_add
  
  status_t all_tests();

} // namespace sho_projection
