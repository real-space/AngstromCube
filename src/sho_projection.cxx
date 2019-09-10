#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy, std::fill
#include <cmath> // std::floor, std::pow

#include "sho_projection.hxx"
#include "inline_math.hxx" // factorial<real_or_int_t>

// #define FULL_DEBUG
// #define DEBUG

namespace sho_projection {


#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template<typename real_t>
  status_t test_orthogonality(int const numax=7, int const echo=1) {
      if (echo > 0) printf("\n# %s\n", __func__);
      int const dims[] = {42, 41, 40};
      real_space_grid::grid_t<1> g(dims);
      double const sigma = 1.25;
      auto const values = new real_t[g.all()];
      set(values, g.all(), (real_t)0);
      g.set_grid_spacing(0.472432); // 0.25 Angstrom
      double const pos[] = {g.dim('x')*.52*g.h[0], g.dim('y')*.51*g.h[1], g.dim('z')*.50*g.h[2]};
      int const nSHO = sho_tools::nSHO(numax);
      auto const icoeff = new uint8_t[nSHO][4];
      { // scope
          int iSHO = 0;
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
      } // scope
//       auto fac = factorial<double>;
//       double const pi_factor = pow3(constants::sqrtpi); // == std::pow(constants::pi, 1.5);
      auto const coeff = new real_t[nSHO];
      double maxdev = 0;
      status_t stat = 0;
      for(int i = 0; i < nSHO; ++i) {
//           double const diag = fac(icoeff[i][0]) * fac(icoeff[i][1]) * fac(icoeff[i][2]) 
//                                   * pi_factor / (1 << icoeff[i][3]);
          double const diag = 1.0;
          double const prefactor = sho_prefactor(icoeff[i][0], icoeff[i][1], icoeff[i][2], sigma);
          set(coeff, nSHO, (real_t)0);
          coeff[i] = pow2(prefactor);
          set(values, g.all(), (real_t)0);
          stat += sho_add(values, g, coeff, numax, pos, sigma, 0);
          stat += sho_project(coeff, numax, pos, sigma, values, g, 0);

          if (echo > 8) {
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
          
          if (echo > 9) {
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

          for(int j = 0; j < nSHO; ++j) {
              if (i == j) {
                  if (echo > 2) printf("# diag i=%d  %g  nu=%d  zyx=%x%x%x\n", i, coeff[i] / diag
                                , icoeff[i][3], icoeff[i][2], icoeff[i][1], icoeff[i][0]);
                  double const dev = coeff[i] / diag - 1.0;
                  maxdev = std::max(maxdev, std::abs(dev));
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
      if (echo > 0) printf("# %s %s: max deviation from unity on the diagonal is %.1e\n", __FILE__, __func__, maxdev);
      delete[] values;
      delete[] icoeff;
      delete[] coeff;
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
