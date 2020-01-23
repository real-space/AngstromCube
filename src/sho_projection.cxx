#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy, std::fill
#include <cmath> // std::floor, std::pow
#include <vector> // std::vector<T>

#include "sho_projection.hxx"
#include "inline_math.hxx" // factorial<real_or_int_t>
#include "solid_harmonics.hxx" // ::rlXlm
#include "sho_unitary.hxx" // ::Unitary_SHO_Transform<real_t>

// #define FULL_DEBUG
// #define DEBUG

namespace sho_projection {


#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template<typename real_t>
  status_t test_L2_orthogonality(int const numax=7, int const echo=2) {
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
  } // test_L2_orthogonality

  status_t test_electrostatic_normalization(int const numax=2, int const echo=2) {
      if (echo > 0) printf("\n# %s\n", __func__);
      int const dims[] = {42, 41, 40};
      real_space_grid::grid_t<1> g(dims);
      double const sigma = 1.25;
      typedef double real_t;
      std::vector<real_t> values(g.all(), 0);
      g.set_grid_spacing(0.472432); // 0.25 Angstrom
      if (echo > 1) printf("# %s %s: for sigma = %g with grid spacing %g\n", __FILE__, __func__, sigma, g.h[0]);
      double const pos[] = {g.dim('x')*.52*g.h[0], g.dim('y')*.51*g.h[1], g.dim('z')*.50*g.h[2]};
      int const nSHO = sho_tools::nSHO(numax);
      
      std::vector<int> energy_ordered(nSHO, 0);
      std::vector<int> loop_ordered(nSHO, 0);
      sho_tools::construct_index_table(energy_ordered.data(), numax, sho_tools::order_zyx, loop_ordered.data());
      std::vector<char> sho_label(nSHO*8, '\0');
      sho_tools::construct_label_table(sho_label.data(), numax, sho_tools::order_Ezyx);
      
      sho_unitary::Unitary_SHO_Transform<real_t> u(numax);
      
      std::vector<real_t> coeff(nSHO);
      double maxdev_diag = 0, maxdev_offd = 0;
      status_t stat = 0;
      for(int ell = 0; ell <= numax; ++ell) {
        for(int emm = -ell; emm <= ell; ++emm) {
          int const lm = sho_tools::lm_index(ell, emm);
          
          double const diag = 1.0;
          set(values.data(), g.all(), (real_t)0); // clear all values on the grid

          { // scope: construct non-decaying solid harmonics on the grid r^ell*X_{ell m}
              double v[3], xlm[64]; assert( numax < 8 ); // sufficient up to lmax=7
              for(        int iz = 0; iz < g.dim('z'); ++iz) { v[2] = iz*g.h[2] - pos[2];
                  for(    int iy = 0; iy < g.dim('y'); ++iy) { v[1] = iy*g.h[1] - pos[1];
                      for(int ix = 0; ix < g.dim('x'); ++ix) { v[0] = ix*g.h[0] - pos[0];
                          int const ixyz = (iz*g.dim('y') + iy)*g.dim('x') + ix;
                          solid_harmonics::rlXlm(xlm, numax, v);
                          values[ixyz] = xlm[lm];
                      } // ix
                  } // iy
              } // iz
          } // scope

          stat += sho_project(coeff.data(), numax, pos, sigma, values.data(), g, 0);

          { // scope
              int iSHO = 0;
              for(int nz = 0; nz <= numax; ++nz) {
                  for(int ny = 0; ny <= numax - nz; ++ny) {
                      for(int nx = 0; nx <= numax - nz - ny; ++nx) {
                          coeff[iSHO] *= sho_prefactor(nx, ny, nz, sigma);
                          ++iSHO;
                      } // nx
                  } // ny
              } // nz
          } // scope

          if (echo > 8) {
              int const nu_show = std::min(echo, numax);
              printf("# coefficients (ordered, up to nu = %i):\n", nu_show);
              for(int nzyx = 0; nzyx < sho_tools::nSHO(nu_show); ++nzyx) {
                  int const izyx = loop_ordered[nzyx];
                  assert( nzyx == energy_ordered[izyx] );
                  printf("# nu=%i %s %16.9f\n", sho_tools::get_nu(nzyx), &sho_label[nzyx*8], coeff[izyx]);
              } // nzyx
              printf("\n\n");
          } // echo

          std::vector<double> vnlm(nSHO, 0.0);
          u.transform_vector(vnlm.data(), sho_tools::order_nlm, 
                             coeff.data(), sho_tools::order_zyx, numax, 0);

          { // scope
              for(int l = 0; l <= numax; ++l) {
                  auto f = electrostatics_prefactor(l, sigma);
                  f *= radial_L2_prefactor(l, sigma) / electrostatic_L1_prefactor(l, sigma); // correction
                  for(int m = -l; m <= l; ++m) {
                      int const jlm = sho_tools::lm_index(l, m);
                      vnlm[jlm] *= f; // scale
                  } // m
              } // l
              for(int iSHO = pow2(1 + numax); iSHO < nSHO; ++iSHO) {
                  vnlm[iSHO] = 0; // clear out all nrn > 0 contributions
              } // iSHO
          } // scope
          
          printf("# vnlm for l=%i m=%i\t ", ell,emm);
          for(int j = 0; j < nSHO; ++j) {
              printf(" %11.6f", vnlm[j]);
              
              if (lm == j) {
                  double const dev = vnlm[j] / diag - 1.0; // relative deviation
                  maxdev_diag = std::max(maxdev_diag, std::abs(dev));
                  if (std::abs(dev) > 1e-8) {
                      if (echo > 9) printf("# diagonal l=%i m=%i\t  %g\n", ell,emm, dev);
                      ++stat;
                  }
              } else {
                  double const dev = vnlm[j]; // absolute deviation
                  maxdev_offd = std::max(maxdev_offd, std::abs(dev)); // relative deviation
                  if (dev > 1e-7) {
                      if (echo > 9) printf("# l=%i m=%i  j=%i\t  %g\n", ell,emm, j, vnlm[j]);
                      ++stat;
                  }
              } // lm == j
              
          } // j
          printf("\n");
      }} // lm
      if (echo > 0) printf("# %s %s: max deviation from unity on the diagonal is %.1e\n", __FILE__, __func__, maxdev_diag);
      if (echo > 0) printf("# %s %s: max deviation  of off-diagonal  elements is %.1e\n", __FILE__, __func__, maxdev_offd);
      return stat;
  } // test_electrostatic_normalization
  
  status_t all_tests() {
    auto status = 0;
    status += test_electrostatic_normalization();
    status += test_L2_orthogonality<double>(); // takes a while
//  status += test_L2_orthogonality<float>();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace sho_projection
