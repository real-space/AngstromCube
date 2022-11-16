#include <cstdio> // std::printf
#include <cassert> // assert
#include <algorithm> // std::copy, ::fill
#include <cmath> // std::floor
#include <vector> // std::vector<T>

#include "sho_projection.hxx"

#include "inline_math.hxx" // pow2
#include "sho_unitary.hxx" // ::Unitary_SHO_Transform<real_t>
#include "solid_harmonics.hxx" // ::rlXlm, ::cleanup<real_t>

namespace sho_projection {

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline constexpr char const * _diag(int const on) { return on ? "diagonal" : "off-diag"; }

  template <typename real_t>
  inline status_t _analyze_row(
        double maxdev[2]
      , int const i
      , real_t const out[]
      , int const nj
      , int const echo=0
  ) {
      status_t stat(0);
      double const threshold = (8 == sizeof(real_t)) ? 1e-8 : 2e-5;
      for (int j = 0; j < nj; ++j) {
          int const d = (i == j); // d==0:off-diagonal, d==1:diagonal
          double const dev = std::abs(out[j] - d);
          maxdev[d] = std::max(maxdev[d], dev);
          if (echo > 9) std::printf("# %s %s i=%i j=%i\t  %g %g\n", __func__, _diag(d), i, j, out[j], dev);
          if (dev > threshold) {
              if (echo > 7) std::printf("# %s %s i=%i j=%i\t  %g\tdev= %g\n", __func__, _diag(d), i, j, out[j], dev);
              ++stat;
          } // deviation large
      } // j
      return stat;
  } // _analyze_row

  template <typename real_t>
  status_t test_L2_orthogonality(int const echo=2, int const numax=5, double const sigma=1.05) {
      if (echo > 0) std::printf("\n# %s<%s>\n", __func__, (8 == sizeof(real_t))?"double":"float");
      int const dims[] = {42, 41, 40};
      real_space::grid_t g(dims);
      std::vector<real_t> values(g.all(), 0);
      g.set_grid_spacing(0.472432); // 0.25 Angstrom
      if (echo > 1) std::printf("# %s %s: for sigma= %g numax= %i with grid spacing %g\n", __FILE__, __func__, sigma, numax, g.h[0]);
      double const pos[] = {g[0]*.52*g.h[0], g[1]*.51*g.h[1], g[2]*.50*g.h[2]};
      int const nSHO = sho_tools::nSHO(numax);
      auto const icoeff = new uint8_t[nSHO][4];
      { // scope
          int iSHO{0};
          for (int nz = 0; nz <= numax; ++nz) {
              for (int ny = 0; ny <= numax - nz; ++ny) {
                  for (int nx = 0; nx <= numax - nz - ny; ++nx) {
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
      double maxdev[] = {0, 0}; // {off-diagonal, diagonal}
      status_t stat(0);
      for (int i = 0; i < nSHO; ++i) {
          double const prefactor = sho_prefactor(icoeff[i][0], icoeff[i][1], icoeff[i][2], sigma);
          set(coeff.data(), nSHO, (real_t)0);
          coeff[i] = pow2(prefactor);
          set(values.data(), g.all(), (real_t)0);
          stat += sho_add(values.data(), g, coeff.data(), numax, pos, sigma, 0);
          stat += sho_project(coeff.data(), numax, pos, sigma, values.data(), g, 0);

          if (echo > 8) {
              int const nu_show = std::min(echo, numax);
              std::printf("# coefficients (up to nu= %i):\n", nu_show);
              int iSHO = 0;
              for (int nz = 0; nz <= numax; ++nz) {
                  for (int ny = 0; ny <= numax - nz; ++ny) {
                      for (int nx = 0; nx <= numax - nz - ny; ++nx) {
                          int const nu = sho_tools::get_nu(nx, ny, nz);
                          if (nu <= nu_show) std::printf("# %x%x%x nu=%i %16.9f\n", nz,ny,nx, nu, coeff[iSHO]);
                          ++iSHO;
                      } // nx
                  } // ny
              } // nz
              std::printf("\n\n");
          } // 1

          if (echo > 9) {
              std::vector<int> energy_ordered(nSHO, 0);
              std::vector<int> loop_ordered(nSHO, 0);
              sho_tools::construct_index_table(energy_ordered.data(), numax, sho_tools::order_zyx, loop_ordered.data());

              int const nu_show = std::min(echo, numax);
              std::printf("# coefficients (ordered, up to nu= %i):\n", nu_show);
              for (int nzyx = 0; nzyx < sho_tools::nSHO(nu_show); ++nzyx) {
                  int const izyx = loop_ordered[nzyx];
                  assert( nzyx == energy_ordered[izyx] );
                  std::printf("# nu=%i %16.9f\n", sho_tools::get_nu(nzyx), coeff[izyx]);
              } // nzyx
              std::printf("\n\n");
          } // echo

          stat += _analyze_row(maxdev, i, coeff.data(), nSHO, echo);
      } // i
      if (echo > 0) {      
          for (int d = 0; d <= 1; ++d) {
              std::printf("# %s %s: max deviation of %s elements is %.1e\n", __FILE__, __func__, _diag(d), maxdev[d]);
          } // d
      } // echo
      delete[] icoeff;
      return stat;
  } // test_L2_orthogonality

  template <typename real_t=double>
  status_t test_renormalize_electrostatics(int const echo=2, int const numax=2) {
      double const sigma = 1.0;
      if (echo > 0) std::printf("\n# %s with sigma= %g\n", __func__, sigma);
      int const dims[] = {42, 41, 40};
      real_space::grid_t g(dims);
      std::vector<real_t> values(g.all(), 0);
      g.set_grid_spacing(0.472432); // 0.25 Angstrom
      if (echo > 1) std::printf("# %s %s: for sigma= %g numax= %i with grid spacing %g\n", __FILE__, __func__, sigma, numax, g.h[0]);
      double const pos[] = {g[0]*.52*g.h[0], g[1]*.51*g.h[1], g[2]*.50*g.h[2]};
      int const nSHO = sho_tools::nSHO(numax);

      std::vector<int> energy_ordered(nSHO, 0);
      std::vector<int> loop_ordered(nSHO, 0);
      sho_tools::construct_index_table(energy_ordered.data(), numax, sho_tools::order_zyx, loop_ordered.data());
//       std::vector<char> sho_label(nSHO*8);
//       sho_tools::construct_label_table(sho_label.data(), numax, sho_tools::order_Ezyx);

      sho_unitary::Unitary_SHO_Transform<real_t> u(numax);

//       solid_harmonics::cleanup();

      std::vector<real_t> coeff(nSHO);
      double maxdev[] = {0, 0}; // {off-diagonal, diagonal}
      status_t stat(0);
      for (int ell = 0; ell <= numax; ++ell) { // angular momentum quantum number
          for (int emm = -ell; emm <= ell; ++emm) { // magnetic quantum number
              int const lm = sho_tools::lm_index(ell, emm);

              set(values.data(), g.all(), real_t(0)); // clear all values on the grid

              { // scope: construct non-decaying solid harmonics on the grid r^ell*X_{ell m}
                  double v[3], xlm[64]; assert( numax < 8 ); // sufficient up to lmax=7
                  for (        int iz = 0; iz < g('z'); ++iz) { v[2] = iz*g.h[2] - pos[2];
                      for (    int iy = 0; iy < g('y'); ++iy) { v[1] = iy*g.h[1] - pos[1];
                          for (int ix = 0; ix < g('x'); ++ix) { v[0] = ix*g.h[0] - pos[0];
                              int const ixyz = (iz*g('y') + iy)*g('x') + ix;
                              solid_harmonics::rlXlm(xlm, numax, v);
                              values[ixyz] = xlm[lm];
//                            std::printf("# %s %s: values[%d][%d]= %g\n", __FILE__, __func__, lm,ixyz, values[ixyz]);
    //                        if (00 == lm) std::printf(" %g", values[ixyz]); // found 0.282095 == 1/sqrt(4*pi)
                          } // ix
                      } // iy
                  } // iz
              } // scope

              stat += sho_project(coeff.data(), numax, pos, sigma, values.data(), g, 0);
//               for (int izyx = 0; izyx < nSHO; ++izyx) {
//                   std::printf("# %s %s: coefficients[%d][%d]= %g\n", __FILE__, __func__, lm,izyx, coeff[izyx]);
//               } // izyx

              int const nlm = pow2(1 + numax);
              std::vector<double> vlm(nlm, 0.0);
              stat += renormalize_electrostatics(vlm.data(), coeff.data(), numax, sigma, u, echo);

              for (int ilm = 0; ilm < nlm; ++ilm) {
//                std::printf("# %s %s: element[%d][%d]= %g\n", __FILE__, __func__, lm,ilm, vlm[ilm]);
                  vlm[ilm] = std::abs(vlm[ilm]); // signs may differ
              } // ilm
              // analyze the matrix row, warning: the matrix is rectangular for numax > 1
              stat += _analyze_row(maxdev, lm, vlm.data(), nlm, echo);
          } // emm
      } // ell
      if (echo > 0) {
          for (int d = 0; d <= 1; ++d) {
              std::printf("# %s %s: max deviation of %s elements is %.1e\n", __FILE__, __func__, _diag(d), maxdev[d]);
          } // d
      } // echo
      return stat;
  } // test_renormalize_electrostatics

  status_t test_L2_prefactors(int const echo=0, int const numax=9) {
      double dev{0};
      for (double sigma = 0.25; sigma < 5; sigma *= 2) {
          for (int z = 0; z <= numax; ++z) {           auto const rz = sho_1D_prefactor(z, sigma);
          for (int y = 0; y <= numax - z; ++y) {       auto const ry = sho_1D_prefactor(y, sigma);
          for (int x = 0; x <= numax - z - y; ++x) {   auto const rx = sho_1D_prefactor(x, sigma);
              dev = std::max(dev, std::abs(sho_prefactor(x, y, z, sigma) - rx*ry*rz));
          }}} // zyx
      } // sigma
      if (echo > 3) std::printf("\n# %s: deviation is %.1e\n\n", __func__, dev);
      return (dev > 1e-12);
  } // test_L2_prefactors

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_L2_prefactors(echo);
      stat += test_renormalize_electrostatics(echo);
      stat += test_L2_orthogonality<double>(echo); // takes a while
//    stat += test_L2_orthogonality<float>(echo);
      if (stat) warn("status= %i", int(stat));
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace sho_projection
