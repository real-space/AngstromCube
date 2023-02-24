#pragma once

#include <cstdint> // uint32_t
#include <cstdio> // std::printf
#include <complex> // std::real

#include "status.hxx" // status_t

#include "sho_tools.hxx" // ::nSHO, ::get_nu, ::order_zyx
#include "real_space.hxx" // ::grid_t<D0>
#include "hermite_polynomial.hxx" // Gauss_Hermite_polynomials
#include "inline_math.hxx" // pow3, factorial<T>
#include "constants.hxx" // ::sqrtpi
#include "sho_unitary.hxx" // Unitary_SHO_Transform


namespace sho_projection {
  // Projection and expansion (aka addition) operations using factorizable
  // Hermite Gauss functions, here called Spherical Harmonic Oscillator basis.

  template <typename real_t>
  inline real_t truncation_radius(real_t const sigma, int const numax=-1) { return 9*sigma; }

  template <typename complex_t, int PROJECT0_OR_ADD1> inline
  status_t _sho_project_or_add(
        complex_t coeff[] // result if projecting, coefficients are zyx-ordered
      , int const numax // how many
      , double const center[3] // where
      , double const sigma
      , complex_t values[] // grid array, result if adding
      , real_space::grid_t const &g // grid descriptor, assume that g is a Cartesian grid
      , int const echo=0 // log-level
  ) {
      using real_t = decltype(std::real(complex_t(1))); // base type

      auto const rcut = truncation_radius(sigma, numax);
      assert(sigma > 0);
      double const sigma_inv = 1./sigma;
      // determine the limitations of the projection domain
      int off[3], end[3], num[3];
      for (int d = 0; d < 3; ++d) {
          if (echo > 9) std::printf("# %c-direction: center= %g, rcut= %g", 120+d, center[d]*g.inv_h[d], rcut*g.inv_h[d]);
          off[d] = std::ceil((center[d] - rcut)*g.inv_h[d]);
          end[d] = std::ceil((center[d] + rcut)*g.inv_h[d]);
          if (echo > 9) std::printf(", prelim limits [%d, %d)", off[d], end[d]);
          off[d] = std::max(off[d], 0); // lower
          end[d] = std::min(end[d], g[d]); // upper boundary
          if (echo > 9) std::printf(", limits [%d, %d)\n", off[d], end[d]);
          num[d] = std::max(0, end[d] - off[d]);
      } // d
      auto const nvolume = (size_t(num[0]) * num[1]) * num[2];
      if ((nvolume < 1) && (echo < 7)) return 0; // no range
      if (echo > 2) std::printf("# %s on rectangular sub-domain x:[%d, %d) y:[%d, %d) y:[%d, %d) = %d * %d * %d = %ld points\n",
                           (0 == PROJECT0_OR_ADD1)?"project":"add", off[0], end[0], off[1], end[1], off[2], end[2],
                           num[0], num[1], num[2], nvolume);

      int const nSHO = sho_tools::nSHO(numax);
      if (0 == PROJECT0_OR_ADD1) set(coeff, nSHO, complex_t(0));

      if (nvolume < 1) return 0; // no range

      // ToDo: analyze if the grid spacing is small enough for this \sigma

      int const M = sho_tools::n1HO(numax);
      std::vector<real_t> H1d[3];
      for (int dir = 0; dir < 3; ++dir) {
          H1d[dir] = std::vector<real_t>(num[dir]*M); // get memory
          auto const h1d = H1d[dir].data();

          double const grid_spacing = g.h[dir];
          if (echo > 55) std::printf("\n# Hermite polynomials for %c-direction:\n", 'x' + dir);
          for (int ii = 0; ii < num[dir]; ++ii) {
              int const ix = ii + off[dir]; // offset
              real_t const x = (ix*grid_spacing - center[dir])*sigma_inv;
              Gauss_Hermite_polynomials(h1d + ii*M, x, numax);
#ifdef DEVEL
              if (echo > 55) {
                  std::printf("%g\t", x);
                  for (int nu = 0; nu <= numax; ++nu) {
                      std::printf(" %11.6f", h1d[ii*M + nu]);
                  } // nu
                  std::printf("\n");
              } // echo
#endif // DEVEL
          } // i
      } // dir


#ifdef DEVEL
      if (1 == PROJECT0_OR_ADD1) {
          if (echo > 6) {
              std::printf("# addition coefficients ");
              for (int iSHO = 0; iSHO < nSHO; ++iSHO) {
                  std::printf(" %g", std::real(coeff[iSHO]));
              } // iSHO
              std::printf("\n\n");
          } // echo
      } // ADD
#endif // DEVEL


      for (        int iz = 0; iz < num[2]; ++iz) {
          for (    int iy = 0; iy < num[1]; ++iy) {
              for (int ix = 0; ix < num[0]; ++ix) {
                  int const ixyz = ((iz + off[2])*g('y') + (iy + off[1]))*g('x') + (ix + off[0]);

                  complex_t val(0);
                  if (true) {
//                    if (echo > 6) std::printf("%g %g\n", std::sqrt(vz*vz + vy*vy + vx*vx), val); // plot function value vs r
                      int iSHO{0};
                      for (int nz = 0; nz <= numax; ++nz) {                    auto const H1d_z = H1d[2][iz*M + nz];
                          for (int ny = 0; ny <= numax - nz; ++ny) {           auto const H1d_y = H1d[1][iy*M + ny];
                              for (int nx = 0; nx <= numax - nz - ny; ++nx) {  auto const H1d_x = H1d[0][ix*M + nx];
                                  auto const H3d = H1d_z * H1d_y * H1d_x;
                                  if (1 == PROJECT0_OR_ADD1) {
                                      val += coeff[iSHO] * H3d; // here, the addition happens
                                  } else {
                                      coeff[iSHO] += val * H3d; // here, the projection happens
                                  }
                                  ++iSHO; // in sho_tools::zyx_order
                              } // nx
                          } // ny
                      } // nz
                      assert( nSHO == iSHO );
                  } // true
                  if (1 == PROJECT0_OR_ADD1) {
                      values[ixyz] += val; // load-modify-store, must be atomic if threads are involved
                  } // write back (add)

              } // ix
          } // iy
      } // iz

      if (0 == PROJECT0_OR_ADD1) scale(coeff, nSHO, complex_t(g.dV())); // volume element of the grid

#ifdef DEVEL
      if (0 == PROJECT0_OR_ADD1) {
          if (echo > 6) {
              std::printf("# projection coefficients ");
              for (int iSHO = 0; iSHO < nSHO; ++iSHO) {
                  std::printf(" %g", std::real(coeff[iSHO]));
              } // iSHO
              std::printf("\n\n");
          } // echo
      } // PROJECT
#endif // DEVEL

      return 0; // success
  } // _sho_project_or_add


  template <typename complex_t>
  status_t sho_project( // wrapper function
        complex_t coeff[] // result, coefficients are zyx-ordered
      , int const numax // SHO basis size
      , double const center[3] // where
      , double const sigma // SHO basis spread
      , complex_t const values[] // input, grid array
      , real_space::grid_t const &g // grid descriptor, assume that g is a Cartesian grid
      , int const echo=0 //
  ) {
      return _sho_project_or_add<complex_t,0>(coeff, numax, center, sigma, (complex_t*)values, g, echo); // un-const values pointer
  } // sho_project

  template <typename complex_t>
  status_t sho_add( // wrapper function
        complex_t values[] // result gets modified, grid array
      , real_space::grid_t const &g // grid descriptor, assume that g is a Cartesian grid
      , complex_t const coeff[] // input, coefficients are zyx-ordered
      , int const numax // SHO basis size
      , double const center[3] // where
      , double const sigma // SHO basis spread
      , int const echo=0 // log-level
  ) {
      return _sho_project_or_add<complex_t,1>((complex_t*)coeff, numax, center, sigma, values, g, echo); // un-const coeff pointer
  } // sho_add


  inline double sho_1D_prefactor(int const nu, double const sigma) {
      return std::sqrt( ( 1 << nu ) / ( constants::sqrtpi * sigma * factorial(nu) ) ); // 1 << nu == 2^nu
  } // sho_1D_prefactor (L2-normalization)

  inline double sho_prefactor(int const nx, int const ny, int const nz, double const sigma) {
//    return sho_1D_prefactor(nx, sigma) * sho_1D_prefactor(ny, sigma) * sho_1D_prefactor(nz, sigma);
      return std::sqrt(   double(1 << nx)*double(1 << ny)*double(1 << nz) /
                        (  factorial(nx) * factorial(ny) * factorial(nz) * pow3(constants::sqrtpi * sigma) ) );
  } // sho_prefactor (L2-normalization)

  template <typename real_t>
  std::vector<real_t> get_sho_prefactors(
        int const numax
      , double const sigma
  ) {
      int const nSHO = sho_tools::nSHO(numax);
      std::vector<real_t> f(nSHO);
      int iSHO{0};
      for (int nz = 0; nz <= numax; ++nz) {                    auto const fz = sho_1D_prefactor(nz, sigma);
          for (int ny = 0; ny <= numax - nz; ++ny) {           auto const fy = sho_1D_prefactor(ny, sigma);
              for (int nx = 0; nx <= numax - nz - ny; ++nx) {  auto const fx = sho_1D_prefactor(nx, sigma);
                  f[iSHO] = fx * fy * fz;
                  ++iSHO;
              } // nx
          } // ny
      } // nz
      assert(nSHO == iSHO);
      return f;
  } // get_sho_prefactors


  template <typename real_t> inline
  status_t renormalize_coefficients(
        real_t out[] // normalized with sho_prefactor [and energy ordered], de-normalized if inverse
      , real_t const in[] // zyx_ordered, input unnormalized, if inverse input is assumed normalized
      , int const numax
      , double const sigma
  ) {
      int iSHO{0};
      for (int nz = 0; nz <= numax; ++nz) {                    auto const fz = sho_1D_prefactor(nz, sigma);
          for (int ny = 0; ny <= numax - nz; ++ny) {           auto const fy = sho_1D_prefactor(ny, sigma);
              for (int nx = 0; nx <= numax - nz - ny; ++nx) {  auto const fx = sho_1D_prefactor(nx, sigma);
                  auto const f = fx * fy * fz;
                  out[iSHO] = in[iSHO] * f;
                  ++iSHO;
              } // nx
          } // ny
      } // nz
      assert( sho_tools::nSHO(numax) == iSHO );
      return 0;
  } // renormalize_coefficients


  inline double radial_L1_prefactor(int const ell, double const sigma) {
      return std::sqrt(2) / ( constants::sqrtpi * std::pow(sigma, 2*ell + 3) * factorial<2>(2*ell + 1) );
  } // radial_L1_prefactor (L1-normalization)

  inline double radial_L2_prefactor(int const ell, double const sigma, int const nrn=0) {
      assert( 0 == nrn ); // derived only for the node-less radial SHO eigenfunctions
      double const fm2 = std::pow(sigma, 2*ell + 3) * factorial<2>(2*ell + 1) * constants::sqrtpi * std::pow(0.5, 2 + ell);
      return 1/std::sqrt(fm2);
  } // radial_L2_prefactor (L2-normalization)

  template <typename real_t> inline
  status_t renormalize_radial_coeff(
        real_t out[] // out[(1 + ellmax)^2]
      , real_t const in[] // in[(1 + ellmax)^2], assert all nrn > 0 contributions vanish
      , int const ellmax
      , double const sigma
  ) {
      for (int ell = 0; ell <= ellmax; ++ell) { // angular momentum quantum number
          auto const pfc = radial_L1_prefactor(ell, sigma) / radial_L2_prefactor(ell, sigma); // prefactor correction
          for (int emm = -ell; emm <= ell; ++emm) { // magnetic quantum number
              int const lm = sho_tools::lm_index(ell, emm);
              out[lm] = pfc * in[lm]; // scale and set only nrn==0 contributions
          } // emm
      } // ell
      return 0;
  } // renormalize_radial_coeff

  inline
  status_t renormalize_electrostatics(
        double vlm[] // result vlm[(1 + ellmax)^2]
      , double const vzyx[] // input zyx_ordered, unnormalized, [nSHO(ellmax)]
      , int const ellmax
      , double const sigma
      , sho_unitary::Unitary_SHO_Transform<double> const & u
      , int const echo=0
  ) {
      status_t stat(0);

      int const nSHO = sho_tools::nSHO(ellmax);
      std::vector<double> vzyx_nrm(nSHO, 0.0); // L2-normalized order_zyx

      // rescale with Cartesian L2 normalization factor, order_zyx unchanged
      stat += renormalize_coefficients(vzyx_nrm.data(), vzyx, ellmax, sigma);

      std::vector<double> vnlm(nSHO, 0.0); // L2-normalized order_nlm
      stat += u.transform_vector(vnlm.data(), sho_tools::order_nlm,
                             vzyx_nrm.data(), sho_tools::order_zyx, ellmax, 0);

      stat += renormalize_radial_coeff(vlm, vnlm.data(), ellmax, sigma);

      return stat;
  } // renormalize_electrostatics

  inline
  status_t denormalize_electrostatics(
        double qzyx[] // result qzyx[nSHO(ellmax)], denormalized
      , double const qlm[] // input charge moments, normalized, [pow2(1 + ellmax)]
      , int const ellmax
      , double const sigma
      , sho_unitary::Unitary_SHO_Transform<double> const & u
      , int const echo=0
  ) {
      status_t stat(0);

      int const nSHO = sho_tools::nSHO(ellmax);
      std::vector<double> qnlm(nSHO, 0.0); // L2-normalized order_nlm

      stat += renormalize_radial_coeff(qnlm.data(), qlm, ellmax, sigma);

      std::vector<double> qzyx_nrm(nSHO, 0.0); // L2-normalized order_zyx
      stat += u.transform_vector(qzyx_nrm.data(), sho_tools::order_zyx,
                                     qnlm.data(), sho_tools::order_nlm, ellmax, 0);

      // rescale with Cartesian L2 normalization factor, order_zyx unchanged
      stat += renormalize_coefficients(qzyx, qzyx_nrm.data(), ellmax, sigma);

      return stat;
  } // denormalize_electrostatics

  status_t all_tests(int const echo=0); // declaration only

} // namespace sho_projection
