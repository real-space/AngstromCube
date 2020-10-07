#pragma once

#include "status.hxx" // status_t

#include "fourier_transform.hxx" // ::fft
#include "inline_tools.hxx" // align
#include "vector_math.hxx" // vec<n,T>
#include "constants.hxx" // ::pi

namespace fourier_poisson {

  double constexpr epsilon0 = 4*constants::pi; // in atomic units

  template<typename real_t>
  status_t solve(real_t x[] // (out) solution of Laplacian*x == b
              , real_t const b[] // right hand side b
              , int const ng[3] // grid numbers
              , double const reci[3][4] // shape of the reciprocal space
              , double const factor=-epsilon0
              , int const echo=0) {

      size_t const ng_all = size_t(ng[0]) * ng[1] * ng[2];
      auto const mg_all = align<3>(ng_all); // aligned to 8 real_t numbers
      std::vector<real_t> mem(2*mg_all); // get memory     
      auto const x_Re = mem.data(), x_Im = mem.data() + mg_all; // point to the second half of that array

      status_t stat = 0;
      stat += fourier_transform::fft<real_t, true>(x_Re, x_Im, b, ng); // transform b into reciprocal space

      if (echo > 0) printf("# %s charge neutrality = %g %g\n", __func__, x_Re[0], x_Im[0]);
      x_Re[0] = 0; x_Im[0] = 0; // charge neutrality, clear the k=[0 0 0]-component

      real_t const scale = -factor/ng_all;

      typedef vector_math::vec<3,double> vec3;
      vec3 rec[3]; for(int d = 0; d < 3; ++d) rec[d] = reci[d];

      int const nh[] = {ng[0]/2, ng[1]/2, ng[2]/2};

      for(        int j2 = 0; j2 < ng[2]; ++j2) { int const k2 = j2 - (j2 > nh[2])*ng[2]; vec3 const vec2   = rec[2]*k2;
          for(    int j1 = 0; j1 < ng[1]; ++j1) { int const k1 = j1 - (j1 > nh[1])*ng[1]; vec3 const vec21  = rec[1]*k1 + vec2;
              for(int j0 = 0; j0 < ng[0]; ++j0) { int const k0 = j0 - (j0 > nh[0])*ng[0]; vec3 const vec210 = rec[0]*k0 + vec21;
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

      stat += fourier_transform::fft<real_t, false>(x, x_Im, x_Re, ng); // transform solution x back into real-space

      return stat;
  } // solve

#ifdef NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS
  status_t all_tests(int const echo=0); // declaration only
#endif // NO_UNIT_TESTS

} // namespace fourier_poisson
