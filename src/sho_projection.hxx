#pragma once

#include <cstdint> // uint32_t
#include <cstdio> // printf

typedef int status_t;

#include "sho_tools.hxx" // sho_tools::nSHO, sho_tools::get_nu, sho_tools::order_zyx
#include "real_space_grid.hxx" // real_space_grid::grid_t<D0>
#include "hermite_polynomial.hxx" // hermite_polys
#include "inline_math.hxx" // pow3
#include "constants.hxx" // sqrtpi, factorial<>

namespace sho_projection {

  
  template<typename real_t>
  inline real_t truncation_radius(real_t const sigma, int const numax=-1) { return 9*sigma; }
  
  template<typename real_t, int D0, int PROJECT0_OR_ADD1>
  status_t _sho_project_or_add(real_t coeff[] // result if projecting, coefficients are zyx-ordered
                     , int const numax // how many
                     , double const center[3] // where
                     , double const sigma
                     , real_t values[] // grid array, result if adding
                     , real_space_grid::grid_t<D0> const &g // grid descriptor, assume that g is a Cartesian grid
                     , int const echo=4) { //
      assert( 1 == D0 );
      double const rcut = truncation_radius(sigma, numax);
      double const sigma_inv = 1./sigma;
      // determine the limitations of the projection domain
      int off[3], end[3], num[3];
      for(int dir = 0; dir < 3; ++dir) {
          off[dir] = std::ceil((center[dir] - rcut)*g.inv_h[dir]);
          end[dir] = std::ceil((center[dir] + rcut)*g.inv_h[dir]);
          if (echo > 9) printf("# prelim for %c-direction are [%d, %d)\n", 120+dir, off[dir], end[dir]);
          off[dir] = std::max(off[dir], 0); // lower
          end[dir] = std::min(end[dir], g.dim(dir)); // upper boundary
          if (echo > 9) printf("# limits for %c-direction are [%d, %d)\n", 120+dir, off[dir], end[dir]);
          num[dir] = std::max(0, end[dir] - off[dir]);
      } // dir
      long const nvolume = num[0] * num[1] * num[2];
      if ((nvolume < 1) && (echo < 7)) return 0; // no range
      if (echo > 2) printf("# %s on rectangular sub-domain x:[%d, %d) y:[%d, %d) y:[%d, %d) = %d * %d * %d = %ld points\n", 
                           (0 == PROJECT0_OR_ADD1)?"project":"add", off[0], end[0], off[1], end[1], off[2], end[2],
                           num[0], num[1], num[2], nvolume);
      if (nvolume < 1) return 0; // no range

      // ToDo: analyze if the grid spacing is small enough for this \sigma

      int const M = 1 + numax;
      real_t* H1d[3];
      for(int dir = 0; dir < 3; ++dir) {
          H1d[dir] = new real_t[num[dir]*M]; // get memory

          real_t const grid_spacing = g.h[dir];
          if (echo > 5) printf("\n# Hermite polynomials for %c-direction:\n", 120+dir);
          for(int ii = 0; ii < num[dir]; ++ii) {
              int const ix = ii + off[dir]; // offset
              real_t const x = (ix*grid_spacing - center[dir])*sigma_inv;
              hermite_polys(H1d[dir] + ii*M, x, numax);
              if (echo > 5) { printf("%g\t", x); for(int nu = 0; nu <= numax; ++nu) printf("%12.6f", H1d[dir][ii*M + nu]); printf("\n"); }
          } // i
      } // dir
   
      int const nSHO = sho_tools::nSHO(numax);
      if (0 == PROJECT0_OR_ADD1) set(coeff, nSHO*D0, real_t(0));
      
      // int ixyz = 0;
      for(        int iz = 0; iz < num[2]; ++iz) {
          for(    int iy = 0; iy < num[1]; ++iy) {
              for(int ix = 0; ix < num[0]; ++ix) {
                  int const ixyz = ((iz + off[2])*g.dim('y') + (iy + off[1]))*g.dim('x') + (ix + off[0]);
                  
                  real_t val[D0];
                  for(int i0 = 0; i0 < D0; ++i0) {
                      val[i0] = values[ixyz*D0 + i0];
                  } // i0 vectorization
                  if (true) {
//                    if (echo > 6) printf("%g %g\n", std::sqrt(vz*vz + vy*vy + vx*vx), val[0]); // plot function value vs r
                      int iSHO = 0;
                      for(int nz = 0; nz <= numax; ++nz) {
                          for(int ny = 0; ny <= numax - nz; ++ny) {
                              for(int nx = 0; nx <= numax - nz - ny; ++nx) {
                                  auto const H3d = H1d[2][iz*M + nz]
                                                 * H1d[1][iy*M + ny]
                                                 * H1d[0][ix*M + nx];
                                  for(int i0 = 0; i0 < D0; ++i0) {
                                      if (1 == PROJECT0_OR_ADD1) {
                                          val[i0] += coeff[iSHO*D0 + i0] * H3d; // here, the addition happens                                          
                                      } else {
                                          coeff[iSHO*D0 + i0] += val[i0] * H3d; // here, the projection happens                                          
                                      }
                                  } // i0 vectorization
                                  ++iSHO;
                              } // nx
                          } // ny
                      } // nz
                      assert( nSHO == iSHO );
                  } // true
                  if (1 == PROJECT0_OR_ADD1) {
                      for(int i0 = 0; i0 < D0; ++i0) {
                          values[ixyz*D0 + i0] = val[i0];
                      } // i0 vectorization
                  } // write back (add)
                  
              } // ix
          } // iy
      } // iz

      if (0 == PROJECT0_OR_ADD1) scale(coeff, nSHO, (real_t)g.dV()); // volume element of the grid

      for(int dir = 0; dir < 3; ++dir) {
          delete[] H1d[dir]; // free memory
      } // dir
      return 0; // success
  } // _sho_project_or_add
  

  // wrapper function
  template<typename real_t, int D0>
  status_t sho_project(real_t coeff[] // result, coefficients are zyx-ordered
                     , int const numax // how many
                     , double const center[3] // where
                     , double const sigma
                     , real_t const values[] // input, grid array
                     , real_space_grid::grid_t<D0> const &g // grid descriptor, assume that g is a Cartesian grid
                     , int const echo=0) { //
      return _sho_project_or_add<real_t,D0,0>(coeff, numax, center, sigma, (real_t*)values, g, echo); // un-const values pointer
  } // sho_project

  // wrapper function
  template<typename real_t, int D0>
  status_t sho_add(real_t values[] // result gets modified, grid array
                 , real_space_grid::grid_t<D0> const &g // grid descriptor, assume that g is a Cartesian grid
                 , real_t const coeff[] // input, coefficients are zyx-ordered
                 , int const numax // how many
                 , double const center[3] // where
                 , double const sigma
                 , int const echo=0) { //
      return _sho_project_or_add<real_t,D0,1>((real_t*)coeff, numax, center, sigma, values, g, echo); // un-const coeff pointer
  } // sho_add
  
  inline double sho_prefactor(int const nx, int const ny, int const nz, double const sigma) { 
      return std::sqrt( ( 1 << sho_tools::get_nu(nx, ny, nz) )
                       /(   pow3(constants::sqrtpi * sigma) 
                          * factorial<double>(nx) 
                          * factorial<double>(ny) 
                          * factorial<double>(nz) 
                        )
                      ); }

  status_t all_tests();

} // namespace sho_projection
