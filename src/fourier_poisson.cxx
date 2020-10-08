#include <cstdio> // printf
#include <cmath> // std::cos, std::abs, std::sqrt
#include <cassert> // assert
#include <vector> // std::vector<T>

#include "fourier_poisson.hxx" // solve

#include "inline_math.hxx" // pow2
#include "constants.hxx" // pi

namespace fourier_poisson {

#ifndef NO_UNIT_TESTS

  status_t test_FFT_Poisson_solver(int const echo=3) {
      if (echo > 1) printf("\n# %s:\n", __func__);
      auto const pi = constants::pi;
      status_t stat(0);
      int const ng[3] = {32, 32, 32}, ngall = ng[2]*ng[1]*ng[0];
      double const mat[3][4] = {{2*pi/ng[0],0,0, 0},{0,2*pi/ng[1],0, 0}, {0,0,2*pi/ng[2], 0}};
      double const alpha = 1./pow2(8.); // in units of h^-2
      std::vector<double> rho(ngall), V(ngall);
      double charge{0};
      for(int i01 = 0; i01 <= 1; ++i01) {
          double q{0}; // count charge
          for(int z = 0; z < ng[2]; ++z) {
          for(int y = 0; y < ng[1]; ++y) {
          for(int x = 0; x < ng[0]; ++x) {
                      double const r2 = pow2(x - .5*ng[0]) + pow2(y - .5*ng[1]) + pow2(z - .5*ng[0]);
                      int const i = (z*ng[1] + y)*ng[0] + x;
                      rho[i] = std::exp(-alpha*r2) - charge;
                      q += rho[i];
                      if (i01 && (echo > 6)) printf("%g %g %g\n", std::sqrt(r2), rho[i], V[i]);
          }}} // zyx
          if (0 == i01) {
              stat += solve(V.data(), rho.data(), ng, mat);
              charge = q/ngall;
          } // first time
          if (echo > 2) printf("# charge in cell %g %g\n", q, charge);
      } // i01
      if (echo > 4) printf("\n# radial density and 1/r Coulomb potential\n");
      double const dr = 1./8., pi4dr = 4*pi*dr;
      double V_rad{0};
      for(int i01 = 0; i01 <= 1; ++i01) {
          double q_rad{0};
          for(int ir = 0; ir < ng[0]/dr; ++ir) {
              auto const r = (ir + .125)*dr, r2 = r*r;
              auto const rho_rad = std::exp(-alpha*r2) - charge;
              // show density, (shifted) potential, integrated charge up to r
              if (i01 && (echo > 4)) printf("%g %g %g %g\n", r, rho_rad, V_rad + q_rad/r, q_rad); 
              q_rad += rho_rad * r2 * pi4dr;
              V_rad -= rho_rad * r  * pi4dr * (2*i01 - 1); // sum up in 1st iteration and subtract in second
          } // ir
          if (echo > 3) printf("\n# radial integrated charge %g, V_rad %g\n", q_rad, V_rad);
      } // i01
      if (echo > 1) printf("# %s: status = %i\n\n", __func__, stat);
      return stat;
  } // test_FFT_Poisson_solver

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_FFT_Poisson_solver(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS
  
} // namespace fourier_poisson
