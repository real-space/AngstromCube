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

  inline constexpr char const * _diag(bool const diagonal) { return diagonal ? "diagonal" : "off-diag"; }

  template<typename real_t>
  inline status_t _analyze_row(int const i, real_t const out[], int const nj, double maxdev[2], int const echo=0) {
      status_t stat = 0;
      double const threshold = (8 == sizeof(real_t)) ? 1e-8 : 2e-5;
      for(int j = 0; j < nj; ++j) {
          int const d = (i == j); // d==0:off-diagonal, d==1:diagonal
          double const dev = std::abs(out[j] - d);
          maxdev[d] = std::max(maxdev[d], std::abs(dev));
          if (echo > 9) printf("# %s %s i=%i j=%i\t  %g %g\n", __func__, _diag(d), i, j, out[j], dev);
          if (std::abs(dev) > threshold) {
              if (echo > 7) printf("# %s %s i=%i j=%i\t  %g %g\n", __func__, _diag(d), i, j, out[j], dev);
              ++stat;
          } // deviation large
      } // j
      return stat;
  } // _analyze_row
  
  template<typename real_t>
  status_t test_L2_orthogonality(int const numax=7, int const echo=2) {
      if (echo > 0) printf("\n# %s<%s>\n", __func__, (8 == sizeof(real_t))?"double":"float");
      int const dims[] = {42, 41, 40};
      real_space_grid::grid_t<1> g(dims);
      double const sigma = 1.25;
      std::vector<real_t> values(g.all(), 0);
      g.set_grid_spacing(0.472432); // 0.25 Angstrom
      if (echo > 1) printf("# %s %s: for sigma = %g numax = %i with grid spacing %g\n", __FILE__, __func__, sigma, numax, g.h[0]);
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
      double maxdev[] = {0, 0}; // {off-diagonal, diagonal}
      status_t stat = 0;
      for(int i = 0; i < nSHO; ++i) {
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

          stat += _analyze_row(i, coeff.data(), nSHO, maxdev, echo);
      } // i
      for(int d = 0; d <= 1; ++d) {
          if (echo > 0) printf("# %s %s: max deviation of %s elements is %.1e\n", __FILE__, __func__, _diag(d), maxdev[d]);
      } // d
      delete[] icoeff;
      return stat;
  } // test_L2_orthogonality

  status_t test_electrostatic_normalization(int const numax=2, int const echo=2) {
      if (echo > 0) printf("\n# %s\n", __func__);
      int const dims[] = {42, 41, 40};
      real_space_grid::grid_t<1> g(dims);
      double const sigma = 0.95;
      typedef double real_t;
      std::vector<real_t> values(g.all(), 0);
      g.set_grid_spacing(0.472432); // 0.25 Angstrom
      if (echo > 1) printf("# %s %s: for sigma = %g numax = %i with grid spacing %g\n", __FILE__, __func__, sigma, numax, g.h[0]);
      double const pos[] = {g.dim('x')*.52*g.h[0], g.dim('y')*.51*g.h[1], g.dim('z')*.50*g.h[2]};
      int const nSHO = sho_tools::nSHO(numax);
      
      std::vector<int> energy_ordered(nSHO, 0);
      std::vector<int> loop_ordered(nSHO, 0);
      sho_tools::construct_index_table(energy_ordered.data(), numax, sho_tools::order_zyx, loop_ordered.data());
      std::vector<char> sho_label(nSHO*8, '\0');
      sho_tools::construct_label_table(sho_label.data(), numax, sho_tools::order_Ezyx);
      
      sho_unitary::Unitary_SHO_Transform<real_t> u(numax);
      
      std::vector<real_t> coeff(nSHO);
      double maxdev[] = {0, 0}; // {off-diagonal, diagonal}
      status_t stat = 0;
      for(int ell = 0; ell <= numax; ++ell) { // angular momentum quantum number
        for(int emm = -ell; emm <= ell; ++emm) { // magnetic quantum number
          int const lm = sho_tools::lm_index(ell, emm);
          
          set(values.data(), g.all(), (real_t)0); // clear all values on the grid

          { // scope: construct non-decaying solid harmonics on the grid r^ell*X_{ell m}
              double v[3], xlm[64]; assert( numax < 8 ); // sufficient up to lmax=7
              for(        int iz = 0; iz < g.dim('z'); ++iz) { v[2] = iz*g.h[2] - pos[2];
                  for(    int iy = 0; iy < g.dim('y'); ++iy) { v[1] = iy*g.h[1] - pos[1];
                      for(int ix = 0; ix < g.dim('x'); ++ix) { v[0] = ix*g.h[0] - pos[0];
                          int const ixyz = (iz*g.dim('y') + iy)*g.dim('x') + ix;
                          solid_harmonics::rlXlm(xlm, numax, v);
                          values[ixyz] = xlm[lm];
//                        if (00 == lm) printf(" %g", values[ixyz]); // found 0.282095 == 1/sqrt(4*pi)
                      } // ix
                  } // iy
              } // iz
          } // scope

          stat += sho_project(coeff.data(), numax, pos, sigma, values.data(), g, 0);

          { // scope: rescale with Cartesian L2 normalization factor, order_zyx unchanged
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

          // Ensure that coeff is L2 normalized before SHO unitary transform is applied!
          std::vector<double> vnlm(nSHO, 0.0);
          u.transform_vector(vnlm.data(), sho_tools::order_nlm, 
                            coeff.data(), sho_tools::order_zyx, numax, 0);
          // now vnlm is L2 normalized

          { // scope: renormalize
              for(int l = 0; l <= numax; ++l) { // angular momentum quantum number
                  auto const f = radial_L1_prefactor(l, sigma)
                               / radial_L2_prefactor(l, sigma); // correction
                  for(int m = -l; m <= l; ++m) { // magnetic quantum number
                      int const jlm = sho_tools::lm_index(l, m);
                      vnlm[jlm] *= f; // scale
                  } // m
              } // l

              // for numax > 1, there are contributions from nrn > 0
              for(int iSHO = pow2(1 + numax); iSHO < nSHO; ++iSHO) {
                  vnlm[iSHO] = 0; // clear out all nrn > 0 contributions
              } // iSHO
          } // scope
          
          if (echo > 0) {
              printf("# vnlm for l=%i m=%i \t", ell,emm);
              for(int j = 0; j < nSHO; ++j) {
                  printf(" %11.6f", vnlm[j]); // show the matrix row
              } // j
              printf("\n");
          } // echo

          for(int j = 0; j < nSHO; ++j) vnlm[j] = std::abs(vnlm[j]); // some signs do not agree
          
          // analyze the matrix row, warning: the matrix is rectangular for numax > 1
          stat += _analyze_row(lm, vnlm.data(), nSHO, maxdev, echo);
      }} // lm
      for(int d = 0; d <= 1; ++d) {
          if (echo > 0) printf("# %s %s: max deviation of %s elements is %.1e\n", __FILE__, __func__, _diag(d), maxdev[d]);
      } // d
      return stat;
  } // test_electrostatic_normalization

  status_t test_renormalize_electrostatics(int const numax=2, int const echo=2) {
      double const sigma = 1.0;
      if (echo > 0) printf("\n# %s with sigma = %g\n", __func__, sigma);
      int const dims[] = {42, 41, 40};
      real_space_grid::grid_t<1> g(dims);
      typedef double real_t;
      std::vector<real_t> values(g.all(), 0);
      g.set_grid_spacing(0.472432); // 0.25 Angstrom
      if (echo > 1) printf("# %s %s: for sigma = %g numax = %i with grid spacing %g\n", __FILE__, __func__, sigma, numax, g.h[0]);
      double const pos[] = {g.dim('x')*.52*g.h[0], g.dim('y')*.51*g.h[1], g.dim('z')*.50*g.h[2]};
      int const nSHO = sho_tools::nSHO(numax);
      
      std::vector<int> energy_ordered(nSHO, 0);
      std::vector<int> loop_ordered(nSHO, 0);
      sho_tools::construct_index_table(energy_ordered.data(), numax, sho_tools::order_zyx, loop_ordered.data());
      std::vector<char> sho_label(nSHO*8, '\0');
      sho_tools::construct_label_table(sho_label.data(), numax, sho_tools::order_Ezyx);
      
      sho_unitary::Unitary_SHO_Transform<real_t> u(numax);
      
      std::vector<real_t> coeff(nSHO);
      double maxdev[] = {0, 0}; // {off-diagonal, diagonal}
      status_t stat = 0;
      for(int ell = 0; ell <= numax; ++ell) { // angular momentum quantum number
        for(int emm = -ell; emm <= ell; ++emm) { // magnetic quantum number
          int const lm = sho_tools::lm_index(ell, emm);
          
          set(values.data(), g.all(), (real_t)0); // clear all values on the grid

          { // scope: construct non-decaying solid harmonics on the grid r^ell*X_{ell m}
              double v[3], xlm[64]; assert( numax < 8 ); // sufficient up to lmax=7
              for(        int iz = 0; iz < g.dim('z'); ++iz) { v[2] = iz*g.h[2] - pos[2];
                  for(    int iy = 0; iy < g.dim('y'); ++iy) { v[1] = iy*g.h[1] - pos[1];
                      for(int ix = 0; ix < g.dim('x'); ++ix) { v[0] = ix*g.h[0] - pos[0];
                          int const ixyz = (iz*g.dim('y') + iy)*g.dim('x') + ix;
                          solid_harmonics::rlXlm(xlm, numax, v);
                          values[ixyz] = xlm[lm];
//                        if (00 == lm) printf(" %g", values[ixyz]); // found 0.282095 == 1/sqrt(4*pi)
                      } // ix
                  } // iy
              } // iz
          } // scope

          stat += sho_project(coeff.data(), numax, pos, sigma, values.data(), g, 0);

          int const nlm = pow2(1 + numax);
          std::vector<double> vlm(nlm, 0.0);
          stat += renormalize_electrostatics(vlm.data(), coeff.data(), numax, sigma, u, echo);

          for(int ilm = 0; ilm < nlm; ++ilm) vlm[ilm] = std::abs(vlm[ilm]); // signs may differ
          // analyze the matrix row, warning: the matrix is rectangular for numax > 1
          stat += _analyze_row(lm, vlm.data(), nlm, maxdev, echo);
      }} // lm
      for(int d = 0; d <= 1; ++d) {
          if (echo > 0) printf("# %s %s: max deviation of %s elements is %.1e\n", __FILE__, __func__, _diag(d), maxdev[d]);
      } // d
      return stat;
  } // test_renormalize_electrostatics

  status_t test_normalize_electrostatics(int const numax=2, int const echo=9) {
      double const sigma = 1.456; // any
      if (echo > 0) printf("\n# %s with sigma = %g\n", __func__, sigma);
      
      sho_unitary::Unitary_SHO_Transform<double> u(numax);
      status_t stat = (u.test_unitarity(echo/2) > 3e-16);

      int const nSHO = sho_tools::nSHO(numax);
      int const nlm = pow2(1 + numax); // lm is a combined 2D angular momentum and magnetic quantum number
      for(int inverse = 0; inverse <= 1; ++inverse) {
          int const n = inverse ? nlm : nSHO;
          int const m = inverse ? nSHO : nlm;
          if (echo > 3) printf("# %s %s: n=%i m=%i n=%i\n", __FILE__, __func__, n, m, n);
          std::vector<double> inp(n), tmp(m), out(n);
          double maxdev[] = {0, 0}; // {off-diagonal, diagonal}
          for(int i = 0; i < n; ++i) {

              set(out.data(), n, 0.0); // clear
              set(tmp.data(), m, 0.0); // clear
              set(inp.data(), n, 0.0); // clear
              inp[i] = 1.;
  
              if (inverse) {
                  assert( m >= n );
                  stat += denormalize_electrostatics(tmp.data(), inp.data(), numax, sigma, u, echo);
                  if ((0 == i) && (echo > 1)) printf("# %s after denormalization tmp[000] = %g\n", __func__, tmp[0]); 
                  stat += renormalize_electrostatics(out.data(), tmp.data(), numax, sigma, u, echo);
              } else {
                  assert( n >= m );
                  stat += renormalize_electrostatics(tmp.data(), inp.data(), numax, sigma, u, echo);
                  if ((0 == i) && (echo > 1)) printf("# %s after renormalization tmp[00] = %g\n", __func__, tmp[0]); 
                  stat += denormalize_electrostatics(out.data(), tmp.data(), numax, sigma, u, echo);
              } // inverse

              stat += _analyze_row(i, out.data(), std::min(n, m), maxdev, echo);
          } // i
          if ((n > m) && (numax > 1) && (stat >=  6)) stat -=  6;
          if ((n > m) && (numax > 2) && (stat >= 18)) stat -= 18;
          if ((n > m) && (numax > 3) && (stat >= 36)) stat -= 36;
          for(int d = 0; d <= 1; ++d) {
              if (echo > 0) printf("# %s %s: max deviation of %s elements is %.1e\n", __FILE__, __func__, _diag(d), maxdev[d]);
          } // d
      } // inverse
      return stat;
  } // test_normalize_electrostatics
  
  
  status_t all_tests() {
    auto status = 0;
    status += test_normalize_electrostatics();
    status += test_renormalize_electrostatics();
    status += test_electrostatic_normalization();
    status += test_L2_orthogonality<double>(); // takes a while
    status += test_L2_orthogonality<float>();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace sho_projection
