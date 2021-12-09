#include <vector> // std::vector
#include <cstdio> // std::printf
#include <cassert> // assert

#include "radial_potential.hxx"

#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // ::create_radial_grid, ::destroy_radial_grid
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
#include "display_units.h" // eV, _eV, Ang, _Ang
#include "solid_harmonics.hxx" // ::lm_index
#include "constants.hxx" // ::pi

namespace radial_potential {
  
  double Hartree_potential( // returns the Coulomb integral
        double rV[] // result: r*Hartree-potential(r)
      , radial_grid_t const & g // radial grid descriptor
      , double const rho4pi[] // 4*\pi*density(r)
  ) {
      double vH1{0};
      for (int ir = 0; ir < g.n; ++ir) {
          vH1 += rho4pi[ir]*g.rdr[ir];
      } // ir
      double const Coulomb = vH1; // the value to be returned

      double vH2{0};
      rV[0] = 0;
      for (int ir = 1; ir < g.n; ++ir) { // start from 1 since for some radial grids r[ir=0] == 0
          vH2 += rho4pi[ir]*g.r2dr[ir];
          vH1 -= rho4pi[ir]*g.rdr[ir];
          rV[ir] = vH2 + vH1*g.r[ir];
      } // ir

      return Coulomb; // the integral 4 \pi rho(r) * r dr is needed for the Coulomb energy
  } // Hartree_potential (spherical)

  
  void Hartree_potential(
        double vHt[] // result: Hartree-potential_lm(r)
      , radial_grid_t const & g // radial grid descriptor
      , double const rho[] // density_lm(r)
      , int const stride // stride between differen lm-compenents in rho and V
      , ell_QN_t const ellmax // largest angular momentum quantum number
      , double const q0 // =0 // singularity
      , double const qlm[] // =nullptr // external boundary conditions, functionality redundant
  ) {
      std::vector<double> rm(g.n), rl(g.n, 1.0);
      for (int ell = 0; ell <= ellmax; ++ell) { // run forward and serial

          if (0 == ell) {
              for (int ir = 0; ir < g.n; ++ir) {
//                rl[ir] = 1; // r^{0} --> already done at initialization
                  rm[ir] = g.rinv[ir]; // r^{-1}
              } // ir
          } else {
              for (int ir = 0; ir < g.n; ++ir) {
                  rl[ir] *= g.r[ir]; // prepare r^{\ell} for the next iteration
                  rm[ir] *= g.rinv[ir]; // prepare r^{-1-\ell} for the next iteration
              } // ir
          } // 0 == ell

          double const factor = (4*constants::pi)/(2*ell + 1);
          for (int emm = -ell; emm <= ell; ++emm) {
              int const lm = solid_harmonics::lm_index(ell, emm);

              double charge2 = ell ? 0.0 : q0;
              for (int ir = 0; ir < g.n; ++ir) {
                  charge2 += rho[lm*stride + ir]*rl[ir]*g.r2dr[ir];
                  vHt[lm*stride + ir]  = charge2*rm[ir];
              } // ir

              double charge1 = qlm ? qlm[lm] : 0.0;
              for (int ir = g.n - 1; ir >= 0; --ir) {
                  vHt[lm*stride + ir] += charge1*rl[ir]; // beware: it makes a difference if we put this line before the next
                  charge1 += rho[lm*stride + ir]*rm[ir]*g.r2dr[ir];
                  vHt[lm*stride + ir] *= factor;
              } // ir

          } // emm
      } // ell

  } // Hartree_potential (ell-emm resolved)
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_radial_Hartree_potential(int const echo=0) {
      auto & g = *radial_grid::create_radial_grid(512);
      std::vector<double> const rho(g.n, 1); // a constant density rho(r) == 1
      std::vector<double> rVH(g.n), vHt(g.n);
      Hartree_potential(rVH.data(), g, rho.data()); // spherical version
      Hartree_potential(vHt.data(), g, rho.data(), 0, 0); //  lm-version
      double const R = g.rmax, V0 = .50754*R*R; // a small correction by 1.5% is needed here
      if (echo > 3) std::printf("\n# %s: Rmax = %g Bohr, V0 = %g Ha\n", __func__, R, V0);
      if (echo > 4) std::printf("## r, r*V_spherical(r), r*V_00(r), r*analytical(r) in a.u.:\n");
      double devsph{0}, devana{0};
      for (int ir = 0; ir < g.n; ++ir) {
          auto const r = g.r[ir];
          double const analytical = V0 - r*r/6;
          double const rV00 = r*vHt[ir]/(4*constants::pi);
          devsph = std::max(devsph, std::abs(rVH[ir] - rV00));
          devana = std::max(devana, std::abs(rVH[ir] - r*analytical));
          if (echo > 4) std::printf("%g %g %g %g\n", r, rVH[ir], rV00, r*analytical);
      } // ir
      radial_grid::destroy_radial_grid(&g);
      if (echo > 3) std::printf("\n# %s deviates %.1e from spherical and %.1e from analytical\n", __func__, devsph, devana);
      return (devsph > 1e-12);
  } // test_radial_Hartree_potential

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_radial_Hartree_potential(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace radial_potential
