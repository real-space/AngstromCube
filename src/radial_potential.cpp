#include <vector> // std::vector
#include <cstdio> // printf
// #include <cstdlib> // abs
#include <cmath> // sqrt, pow, log
#include <math.h> // C_PI
#include <cassert> // assert

#include "radial_potential.hxx"

#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // create_exponential_radial_grid
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
#include "output_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // pow2

// #define FULL_DEBUG
// #define DEBUG

namespace radial_potential {
  
  double Hartree_potential(
            double rV[],            // r*Hartree-potential(r)
            radial_grid_t const &g, // radial grid descriptor
            double const rho4pi[],  // 4*pi*density(r)
            ell_QN_t const ell) {  // default ell=0
    assert(0 == ell);

    double vH1 = 0;
    for(int ir = 0; ir < g.n; ++ir) {
        vH1 += rho4pi[ir]*g.rdr[ir];
    } // ir
    double const Coulomb = vH1; // return value

    double vH2 = 0;
    rV[0] = 0;
    for(int ir = 1; ir < g.n; ++ir) {
        vH2 += rho4pi[ir]*g.r2dr[ir];
        vH1 -= rho4pi[ir]*g.rdr[ir];
        rV[ir] = vH2 + vH1*g.r[ir];
    } // ir
    // printf("# vH1 = %g after integration\n", vH1); // check that vH1 is small

    return Coulomb; // the integral 4 \pi rho(r) * r dr is needed for the Coulomb energy
  } // Hartree_potential

  
  double lda_PZ81_kernel(double const rho, double &Vup, 
                       double const mag, double *Vdn) { // optional
    
    double constexpr C_PI = 3.14159265358979323846; // pi
    double constexpr THIRD  = 1./3., TINYDEN = 1e-20;
    double const tpt5 = .6108870577108572; // (2.25/(C_PI*C_PI))**THIRD

    if (nullptr == Vdn) { //  LDA

      if (rho < TINYDEN) {
        Vup = 0;
        return 0;
      } else { // negligible
        double Exc;
        auto const rs = pow(3.0/(4.0*C_PI*rho), THIRD);
        if (rs > 1.) {
          auto const srs = sqrt(rs);
          Exc = -0.1423/(1. + 1.0529*srs + 0.3334*rs);
          Vup = Exc - rs/3.0*(0.1423*(0.3334 + 0.52645/srs)/pow2(1. + 1.0529*srs + 0.3334*rs));
        } else {
          auto const lrs = log(rs);
          Exc = -0.048 + 0.0311*lrs - 0.0116*rs + 0.002*rs*lrs;
          Vup = Exc - rs/3.0*(0.0311/rs - 0.0096 + 0.002*lrs);
        } //
        Exc = Exc - 0.75*tpt5/rs;
        Vup = Vup - 0.75*tpt5/rs - rs/3.0*(0.75*tpt5/(rs*rs));
        return Exc;
      } // negligible

    } else { // spin resolved

      if (rho < TINYDEN) {
        *Vdn = 0;
         Vup = 0;
        return 0;
      } else { // negligible
        double Exc = 0;
         Vup = 0;
        *Vdn = 0;
        return Exc;
      } // negligible
    } // Magnetization

  } // lda_PZ81_kernel
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_radial_potential(
    radial_grid_t const g) { // radial grid descriptor
    auto const rho = std::vector<double>(g.n, 1);
    auto       rVH = std::vector<double>(g.n);
    Hartree_potential(rVH.data(), g, rho.data());
    double const R = g.rmax;
    auto const V0 = .50754*R*R; // small correction by 1.5%
    printf("# Rmax = %g V0 = %g \n", R, V0);
    for(int ir = 0; ir < g.n; ++ir) {
        auto const r = g.r[ir];
        double const model = V0 - r*r/6;
        printf("%g %g %g %g\n", r, rVH[ir], r*model, rho[ir]);
    } // ir
    return 0;
  } // test_radial_potential

  status_t test_radial_PZ81_potential() {
    for(double rho = .625e-20; rho < 1e9; rho *= 2) {
        double V; double const E = lda_PZ81_kernel(rho, V);
        printf("%g %g %g\n", rho, E, V);
    } // rho
    return 0;
  } // test_radial_PZ81_potential

  status_t all_tests() {
    auto status = 0;
//  status += test_radial_potential(*create_exponential_radial_grid(512));
    status += test_radial_PZ81_potential();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace radial_potential
