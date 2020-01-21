#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy, std::fill
#include <cmath> // std::floor, std::pow
#include <vector> // std::vector<T>

#include "sho_projection.hxx"
#include "inline_math.hxx" // factorial<real_or_int_t>

// #define FULL_DEBUG
// #define DEBUG

namespace sho_projection {


#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template<typename real_t>
  status_t test_orthogonality(int const numax=7, int const echo=2) {
      if (echo > 0) printf("\n# %s\n", __func__);
      int const dims[] = {42, 41, 40};
      real_space_grid::grid_t<1> g(dims);
      double const sigma = 1.25;
      std::vector<real_t> values(g.all(), 0);
      g.set_grid_spacing(0.472432); // 0.25 Angstrom
      if (echo > 1) printf("# %s %s: for sigma = %g with grid spacing %g\n", __FILE__, __func__, sigma, g.h[0]);
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
      std::vector<real_t> coeff(nSHO);
      double maxdev_diag = 0, maxdev_offd = 0;
      status_t stat = 0;
      for(int i = 0; i < nSHO; ++i) {
          double const diag = 1.0;
          double const prefactor = sho_prefactor(icoeff[i][0], icoeff[i][1], icoeff[i][2], sigma);
          set(coeff.data(), nSHO, (real_t)0);
          coeff[i] = pow2(prefactor);
          set(values.data(), g.all(), (real_t)0);
          stat += sho_add(values.data(), g, coeff.data(), numax, pos, sigma, 0);
          stat += sho_project(coeff.data(), numax, pos, sigma, values.data(), g, 0);

          if (echo > 8) {
              int const nu_show = std::min(echo, numax);
              printf("# coefficients (up to nu = %i):\n", nu_show);
              int iSHO = 0;
              for(int nz = 0; nz <= numax; ++nz) {
                  for(int ny = 0; ny <= numax - nz; ++ny) {
                      for(int nx = 0; nx <= numax - nz - ny; ++nx) {
                          int const nu = sho_tools::get_nu(nx, ny, nz);
                          if (nu <= nu_show) printf("# %x%x%x nu=%i %16.9f\n", nz,ny,nx, nu, coeff[iSHO]);
                          ++iSHO;
                      } // nx
                  } // ny
              } // nz
              printf("\n\n");
          } // 1
          
          if (echo > 9) {
              std::vector<int> energy_ordered(nSHO, 0);
              std::vector<int> loop_ordered(nSHO, 0);
              sho_tools::construct_index_table(energy_ordered.data(), numax, sho_tools::order_zyx, loop_ordered.data());

              int const nu_show = std::min(echo, numax);
              printf("# coefficients (ordered, up to nu = %i):\n", nu_show);
              for(int nzyx = 0; nzyx < sho_tools::nSHO(nu_show); ++nzyx) {
                  int const izyx = loop_ordered[nzyx];
                  assert( nzyx == energy_ordered[izyx] );
                  printf("# nu=%i %16.9f\n", sho_tools::get_nu(nzyx), coeff[izyx]);
              } // nzyx
              printf("\n\n");
          } // echo

          for(int j = 0; j < nSHO; ++j) {
              if (i == j) {
                  if (echo > 2) printf("# diag i=%i  %g  nu=%i  zyx=%x%x%x\n", i, coeff[i] / diag
                                , icoeff[i][3], icoeff[i][2], icoeff[i][1], icoeff[i][0]);
                  double const dev = coeff[i] / diag - 1.0; // relative deviation
                  maxdev_diag = std::max(maxdev_diag, std::abs(dev));
                  if (std::abs(dev) > 1e-8) {
                      printf("# diagonal i=%i  %g\n", i, dev);
                      ++stat;
                  }
              } else {
                  double const dev = coeff[j]; // absolute deviation
                  maxdev_offd = std::max(maxdev_offd, std::abs(dev)); // relative deviation
                  if (dev > 1e-7) {
                      printf("# i=%i j=%i  %g\n", i, j, coeff[j]);
                      ++stat;
                  }
              } // i == j
          } // j
      } // i
      if (echo > 0) printf("# %s %s: max deviation from unity on the diagonal is %.1e\n", __FILE__, __func__, maxdev_diag);
      if (echo > 0) printf("# %s %s: max deviation  of off-diagonal  elements is %.1e\n", __FILE__, __func__, maxdev_offd);
      delete[] icoeff;
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
