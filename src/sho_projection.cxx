#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy, std::fill
#include <cmath> // std::floor

#include "sho_projection.hxx"

#include "sho_tools.hxx" // sho_tools::nSHO
#include "real_space_grid.hxx" // grid_t
#include "hermite_polynomials.hxx" // hermite_polys


#include "display_units.h" // eV, _eV, Ang, _Ang
#include "constants.hxx" // pi

// #define FULL_DEBUG
// #define DEBUG

namespace sho_projection {

  template<typename real_t, int D0>
  status_t sho_project(real_t coeff[] // result, coefficients are zyx-ordered
                     , int const numax // how many
                     , double const center[3] // where
                     , double const sigma
                     , real_space_grid::grid_t<real_t,D0> const g // grid values, assume that g is a Cartesian grid
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
      if (echo > 1) printf("# project on rectangular sub-domain %d x %d x %d = %ld points\n", 
                           num[0], num[1], num[2], nvolume);
      if (nvolume < 1) return 0; // no range

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
              hermite_polys(&(H1d[dir][ii*M]), x, numax);
              if (echo > 5) { printf("%g\t", x); for(int nu = 0; nu <= numax; ++nu) printf("%12.6f", H1d[dir][ii*M + nu]); printf("\n"); }
          } // i
      } // dir
   
      int const nSHO = sho_tools::nSHO(numax);
      for(int iSHO = 0; iSHO < nSHO; ++iSHO) coeff[iSHO] = 0; // clear
      
      int ixyz = 0;
      for(        int iz = 0; iz < g.dim('z'); ++iz) {
          for(    int iy = 0; iy < g.dim('y'); ++iy) {
              for(int ix = 0; ix < g.dim('x'); ++ix) {
                  auto const val = g.values[ixyz];
                  if (0 != val) {
//                    if (echo > 6) printf("%g %g\n", std::sqrt(vz*vz + vy*vy + vx*vx), val); // plot function value vs r
                      int iSHO = 0;
                      for(int nz = 0; nz <= numax; ++nz) {
                          for(int ny = 0; ny <= numax - nz; ++ny) {
                              for(int nx = 0; nx <= numax - nz - ny; ++nx) {
                                  auto const H3d = H1d[2][(iz - off[2])*M + nz]
                                                 * H1d[1][(iy - off[1])*M + ny]
                                                 * H1d[0][(ix - off[0])*M + nx];
                                  coeff[iSHO] += val * H3d; // here, the projection happens                                          
                                  ++iSHO;
                              } // nx
                          } // ny
                      } // nz
                      assert( nSHO == iSHO );
                  } // non-zero
                  ++ixyz;
              } // ix
          } // iy
      } // iz
      // g.dV(); // volume element of the grid has been ignored
      
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
          int energy_ordered[nSHO];
//           sho_tools::zyx_translation_table(energy_ordered, numax);
          int loop_ordered[nSHO];
//           sho_tools::zyx_translation_table(loop_ordered, numax, true);
          sho_tools::construct_index_table(energy_ordered, numax, sho_tools::order_zyx, loop_ordered);

          int const nu_show = std::min(echo, numax);
          printf("# coefficients (ordered, up to nu = %d):\n", nu_show);
          for(int nzyx = 0; nzyx < sho_tools::nSHO(nu_show); ++nzyx) {
              int const izyx = loop_ordered[nzyx];
              assert( nzyx == energy_ordered[izyx] );
              printf("# nu=%d %16.9f\n", sho_tools::get_nu(nzyx), coeff[izyx]);
          } // nzyx
          printf("\n\n");
      } // 1
      
      for(int dir = 0; dir < 3; ++dir) {
          delete[] H1d[dir]; // free memory
      } // dir
      return 0; // success
  } // sho_project
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template<typename real_t>
  status_t test_create_and_destroy(int echo=9) {
      if (echo > 0) printf("\n# %s\n", __func__);
      int const dims[] = {32, 31, 30};
      real_space_grid::grid_t<real_t,1> g(dims);
      std::fill(g.values, g.all() + g.values, 1.0);
      g.set_grid_spacing(0.5);
      double const pos[] = {g.dim('x')*.42*g.h[0], g.dim('y')*.51*g.h[1], g.dim('z')*.60*g.h[2]};
      int constexpr numax = 7;
      real_t coeff[sho_tools::nSHO(numax)];
      auto const stat = sho_project<real_t,1>(coeff, numax, pos, 1.0, g);
      return stat;
  } // test_create_and_destroy

  status_t all_tests() {
    auto status = 0;
    status += test_create_and_destroy<double>();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace sho_projection