#include <vector> // std::vector
#include <cstdio> // printf
// #include <math.h> // C_PI
#include <cassert> // assert

#include "radial_potential.hxx"

#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // create_exponential_radial_grid
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
#include "output_units.h" // eV, _eV, Ang, _Ang
#include "solid_harmonics.hxx" // lm_index
#include "constants.hxx" // pi

// #define FULL_DEBUG
// #define DEBUG

namespace radial_potential {
  
  double Hartree_potential(
            double rV[],            // r*Hartree-potential(r)
            radial_grid_t const &g, // radial grid descriptor
            double const rho4pi[]) { // 4*pi*density(r)

      double vH1 = 0;
      for(int ir = 0; ir < g.n; ++ir) {
          vH1 += rho4pi[ir]*g.rdr[ir];
      } // ir
      double const Coulomb = vH1; // return value

      double vH2 = 0;
      rV[0] = 0;
      for(int ir = 1; ir < g.n; ++ir) { // start from 1 since for some radial grids r[ir=0] == 0
          vH2 += rho4pi[ir]*g.r2dr[ir];
          vH1 -= rho4pi[ir]*g.rdr[ir];
          rV[ir] = vH2 + vH1*g.r[ir];
      } // ir

      return Coulomb; // the integral 4 \pi rho(r) * r dr is needed for the Coulomb energy
  } // Hartree_potential _spherical

  
  void Hartree_potential(
            double vHt[], // Hartree-potential_lm(r)
            radial_grid_t const &g, // radial grid descriptor
            double const rho[],  // density_lm(r)
            int const stride, // stride between differen lm-compenents in rho and V
            ell_QN_t const ellmax,
            double const qlm[], // defaults to nullptr
            double const q0) { // defaults to zero
    
      auto const rl = new double[g.n];
      auto const rm = new double[g.n];

      for(int ell = 0; ell <= ellmax; ++ell) { // run forward and serial
          double const f = (4*constants::pi)/(2.*ell + 1.);

          if (0 == ell) {
              rl[0] = 0; rm[0] = 0;
              for(int ir = 1; ir < g.n; ++ir) {
                  rl[ir] = 1; // r^{0}
                  rm[ir] = g.rinv[ir]; // r^{-1}
              } // ir
          } else {
              for(int ir = 1; ir < g.n; ++ir) {
                  rl[ir] *= g.r[ir]; // prepare r^{\ell} for the next iteration
                  rm[ir] *= g.rinv[ir]; // prepare r^{-1-\ell} for the next iteration
              } // ir
          } // 0 == ell

          for(int emm = -ell; emm <= ell; ++emm) {
              int const lm = solid_harmonics::lm_index(ell, emm);

              double charge2 = ell ? 0 : q0;
              for(int ir = 0; ir < g.n; ++ir) {
                  charge2 += rho[lm*stride + ir]*rl[ir]*g.r2dr[ir];
                  vHt[lm*stride + ir]  = charge2*rm[ir];
              } // ir

              double charge1 = qlm ? qlm[lm] : 0;
              for(int ir = g.n - 1; ir > 0; --ir) {
                  vHt[lm*stride + ir] += charge1*rl[ir];
                  charge1 += rho[lm*stride + ir]*rm[ir]*g.r2dr[ir];
                  vHt[lm*stride + ir] *= f;
              } // ir

          } // emm
          
      } // ell

  } // Hartree_potential_ell
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_radial_potential(
    radial_grid_t const g) { // radial grid descriptor
    auto const rho = std::vector<double>(g.n, 1);
    auto       rVH = std::vector<double>(g.n);
    auto       vHt = std::vector<double>(g.n);
    Hartree_potential(rVH.data(), g, rho.data());
    Hartree_potential(vHt.data(), g, rho.data(), 0, 0);
    double const R = g.rmax;
    auto const V0 = .50754*R*R; // small correction by 1.5%
    printf("# Rmax = %g V0 = %g \n", R, V0);
    for(int ir = 0; ir < g.n; ++ir) {
        auto const r = g.r[ir];
        double const model = V0 - r*r/6;
        printf("%g %g %g %g\n", r, rVH[ir], r*model, r*vHt[ir]/(4*constants::pi));
    } // ir
    return 0;
  } // test_radial_potential

  status_t all_tests() {
    auto status = 0;
    status += test_radial_potential(*radial_grid::create_exponential_radial_grid(512));
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace radial_potential
