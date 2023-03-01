#pragma once

#include <cstdio> // std::printf

#include "status.hxx" // status_t

#include "fourier_transform.hxx" // ::fft
#include "vector_math.hxx" // vec<n,T>
#include "constants.hxx" // ::pi
#include "inline_math.hxx" // pow2, align<nBits>

namespace fourier_poisson {

  double constexpr epsilon0 = 4*constants::pi; // in atomic units

  template <typename real_t>
  status_t solve(
        real_t x[] // (out) solution of Laplacian*x == b
      , real_t const b[] // right hand side b
      , int const ng[3] // grid numbers
      , double const reci[3][4] // shape of the reciprocal space
      , double const factor=-epsilon0
      , int const echo=0
  ) {

      size_t const ng_all = size_t(ng[0]) * size_t(ng[1]) * size_t(ng[2]);
      auto const mg_all = align<3>(ng_all); // aligned to 8 real_t numbers
      std::vector<real_t> mem(3*mg_all, real_t(0)); // get memory
      auto const x_Re = mem.data(),
                 x_Im = mem.data() + mg_all, // point to the second half of that array
                 b_Im = mem.data() + 2*mg_all,
                 neglect = b_Im;

      status_t stat(0);
      stat += fourier_transform::fft(x_Re, x_Im, b, b_Im, ng, true); // transform b into reciprocal space
      if (0 != stat) error("fourier transform failed with status ", stat);

      if (echo > 0) std::printf("# %s charge neutrality = %g %g\n", __func__, x_Re[0], x_Im[0]);
      x_Re[0] = 0; x_Im[0] = 0; // charge neutrality, clear the k=[0 0 0]-component

      real_t const scale = -factor/ng_all;

      typedef vector_math::vec<3,double> vec3;
      vec3 rec[3]; for (int d = 0; d < 3; ++d) rec[d] = reci[d];

      int const nh[] = {ng[0]/2, ng[1]/2, ng[2]/2};

      for (        int j2 = 0; j2 < ng[2]; ++j2) { int const k2 = j2 - (j2 > nh[2])*ng[2]; vec3 const vec2   = rec[2]*k2;
          for (    int j1 = 0; j1 < ng[1]; ++j1) { int const k1 = j1 - (j1 > nh[1])*ng[1]; vec3 const vec21  = rec[1]*k1 + vec2;
              for (int j0 = 0; j0 < ng[0]; ++j0) { int const k0 = j0 - (j0 > nh[0])*ng[0]; vec3 const vec210 = rec[0]*k0 + vec21;
                  int const i = (j2*ng[1] + j1)*ng[0] + j0;

                  int const kk = k0*k0 + k1*k1 + k2*k2;
                  if (kk > 0) {
                      real_t const invLaplacian = scale/norm(vec210);
                      // modify x_Re and x_Im in-place
                      x_Re[i] *= invLaplacian;
                      x_Im[i] *= invLaplacian;
                  } // kk > 0

              } // j0
          } // j1
      } // j2

      stat += fourier_transform::fft(x, neglect, x_Re, x_Im, ng, false); // transform solution x back into real-space

      return stat;
  } // solve








#ifdef NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_FFT_Poisson_solver(int const echo=3) {
      if (echo > 1) std::printf("\n# %s:\n", __func__);
      auto const pi = constants::pi;
      status_t stat(0);
      int const ng[3] = {32, 32, 32}, ngall = ng[2]*ng[1]*ng[0];
      double const mat[3][4] = {{2*pi/ng[0],0,0, 0},{0,2*pi/ng[1],0, 0}, {0,0,2*pi/ng[2], 0}};
      double const alpha = 1./pow2(8.); // in units of h^-2
      std::vector<double> rho(ngall), V(ngall);
      double charge{0};
      for (int i01 = 0; i01 <= 1; ++i01) {
          double q{0}; // count charge
          for (int z = 0; z < ng[2]; ++z) {
          for (int y = 0; y < ng[1]; ++y) {
          for (int x = 0; x < ng[0]; ++x) {
                      double const r2 = pow2(x - .5*ng[0]) + pow2(y - .5*ng[1]) + pow2(z - .5*ng[0]);
                      int const i = (z*ng[1] + y)*ng[0] + x;
                      rho[i] = std::exp(-alpha*r2) - charge;
                      q += rho[i];
                      if (i01 && (echo > 6)) std::printf("%g %g %g\n", std::sqrt(r2), rho[i], V[i]);
          }}} // zyx
          if (0 == i01) {
              stat += solve(V.data(), rho.data(), ng, mat);
              charge = q/ngall;
          } // first time
          if (echo > 2) std::printf("# charge in cell %g %g\n", q, charge);
      } // i01
      if (echo > 4) std::printf("\n# radial density and 1/r Coulomb potential\n");
      double const dr = 1./8., pi4dr = 4*pi*dr;
      double V_rad{0};
      for (int i01 = 0; i01 <= 1; ++i01) {
          double q_rad{0};
          for (int ir = 0; ir < ng[0]/dr; ++ir) {
              auto const r = (ir + .125)*dr, r2 = r*r;
              auto const rho_rad = std::exp(-alpha*r2) - charge;
              // show density, (shifted) potential, integrated charge up to r
              if (i01 && (echo > 4)) std::printf("%g %g %g %g\n", r, rho_rad, V_rad + q_rad/r, q_rad);
              q_rad += rho_rad * r2 * pi4dr;
              V_rad -= rho_rad * r  * pi4dr * (2*i01 - 1); // sum up in 1st iteration and subtract in second
          } // ir
          if (echo > 3) std::printf("\n# radial integrated charge %g, V_rad %g\n", q_rad, V_rad);
      } // i01
      if (echo > 1) std::printf("# %s: status = %i\n\n", __func__, stat);
      return stat;
  } // test_FFT_Poisson_solver

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_FFT_Poisson_solver(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace fourier_poisson
