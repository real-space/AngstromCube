#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy, std::fill
#include <cmath> // std::floor, std::pow

#include "sho_projection.hxx"

#include "sho_tools.hxx" // sho_tools::nSHO, sho_tools::get_nu, sho_tools::order_zyx
#include "real_space_grid.hxx" // real_space_grid::grid_t<real_t,D0>
#include "hermite_polynomial.hxx" // hermite_polys
#include "constants.hxx" // constants::pi
#include "inline_math.hxx" // set, factorial

// #define FULL_DEBUG
// #define DEBUG

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
      
      if (0 == PROJECT0_OR_ADD1) {
          if (echo > 0) {
              int const nu_show = std::min(echo, numax);
              printf("# coefficients (up to nu = %d):\n", nu_show);
              int iSHO = 0;
              for(int nz = 0; nz <= numax; ++nz) {
                  for(int ny = 0; ny <= numax - nz; ++ny) {
                      for(int nx = 0; nx <= numax - nz - ny; ++nx) {
                          int const nu = sho_tools::get_nu(nx, ny, nz);
                          if (nu <= nu_show) printf("# %x%x%x nu=%d %16.9f\n", nz,ny,nx, nu, coeff[iSHO]);
                          ++iSHO;
                      } // nx
                  } // ny
              } // nz
              printf("\n\n");
          } // 1
          
          if (echo > 0) {
              auto const energy_ordered = new int[nSHO];
              auto const loop_ordered = new int[nSHO];
              sho_tools::construct_index_table(energy_ordered, numax, sho_tools::order_zyx, loop_ordered);

              int const nu_show = std::min(echo, numax);
              printf("# coefficients (ordered, up to nu = %d):\n", nu_show);
              for(int nzyx = 0; nzyx < sho_tools::nSHO(nu_show); ++nzyx) {
                  int const izyx = loop_ordered[nzyx];
                  assert( nzyx == energy_ordered[izyx] );
                  printf("# nu=%d %16.9f\n", sho_tools::get_nu(nzyx), coeff[izyx]);
              } // nzyx
              printf("\n\n");
              delete[] loop_ordered;
              delete[] energy_ordered;
          } // 1
      } // PROJECT

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


#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template<typename real_t>
  status_t test_orthogonality(int const echo=1) {
      if (echo > 0) printf("\n# %s\n", __func__);
      int const dims[] = {32, 31, 30};
      real_space_grid::grid_t<real_t,1> g(dims);
      std::fill(g.values, g.all() + g.values, 0.0);
      g.set_grid_spacing(0.45);
      double const pos[] = {g.dim('x')*.52*g.h[0], g.dim('y')*.51*g.h[1], g.dim('z')*.50*g.h[2]};
      int constexpr numax = 7;
      int const ncoeff = sho_tools::nSHO(numax);
      auto const icoeff = new uint8_t[ncoeff][4];
      {   int iSHO = 0;
          for(int nz = 0; nz <= numax; ++nz) {
              for(int ny = 0; ny <= numax - nz; ++ny) {
                  for(int nx = 0; nx <= numax - nz - ny; ++nx) {
                      icoeff[iSHO][0] = nx;
                      icoeff[iSHO][1] = ny;
                      icoeff[iSHO][2] = nz;
                      icoeff[iSHO][3] = sho_tools::get_nu(nx, ny, nz);
                      ++iSHO;
                  } // nx
              } // ny
          } // nz
      }
      auto fac = factorial<double>;
      double const pi_factor = std::pow(constants::pi, 1.5);
      auto const coeff = new real_t[ncoeff];
      status_t stat = 0;
      for(int i = 0; i < ncoeff; ++i) {
          double const diag = fac(icoeff[i][2]) * fac(icoeff[i][1]) * fac(icoeff[i][0]) 
                                  * pi_factor / (1 << icoeff[i][3]);
          set(coeff, ncoeff, (real_t)0);
          coeff[i] = 1;
          set(g.values, g.all(), (real_t)0);
          stat += sho_add(g.values, g, coeff, numax, pos, 1.0, 0);
          stat += sho_project(coeff, numax, pos, 1.0, g, g.values, 0);
          for(int j = 0; j < ncoeff; ++j) {
              if (i == j) {
                  if (echo > 2) printf("# diag i=%d  %g  nu=%d  zyx=%x%x%x\n", i, coeff[i] / diag
                                , icoeff[i][3], icoeff[i][2], icoeff[i][1], icoeff[i][0]);
                  if (std::abs(coeff[i] / diag - 1.0) > 1e-8) {
                      printf("# diagonal i=%d  %g\n", i, coeff[i] / diag - 1.0);
                      ++stat;
                  }
              } else {
                  if (std::abs(coeff[j]) > 1e-7) {
                      printf("# i=%d j=%d  %g\n", i, j, coeff[j]);
                      ++stat;
                  }
              } // i == j
          } // j
      } // i
      return stat;
  } // test_orthogonality

  status_t all_tests() {
    auto status = 0;
    status += test_orthogonality<double>();
//  status += test_orthogonality<float>();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace sho_projection
