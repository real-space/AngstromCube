#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // uint32_t, int16_t
#include <vector> // std::vector<T>

#include "real_space.hxx" // ::grid_t
#include "boundary_condition.hxx" // *_Boundary
#include "recorded_warnings.hxx" // warn
#include "inline_math.hxx" // intpow, pow2

#include "status.hxx" // status_t
#include "uniform_laplacian.hxx" // ::get

namespace finite_difference {

    int constexpr nnArraySize = 16;

    template <typename real_t>
    int set_Laplacian_coefficients(
          real_t c[] // c[1+nn]
        , int const nn=1 // returns nn on success
        , double const grid_spacing=1
        , char const direction='?'
    ) {
            return uniform_laplacian::get(c, nn, grid_spacing, direction);
    } // set_Laplacian_coefficients


  template <typename real_t> // real_t may be float or double
  class stencil_t {
    public:
      real_t c2nd[3][nnArraySize]; // coefficients for the 2nd derivative
    private:
      int8_t _nn[3]; // number of FD neighbors

      void _constructor(double const grid_spacing[3], int const nneighbors[3], double const scale_factor) {
          for (int d = 0; d < 3; ++d) {
              for (int i = 0; i < nnArraySize; ++i) { c2nd[d][i] = 0; } // clear
              _nn[d] = set_Laplacian_coefficients(c2nd[d], nneighbors[d], grid_spacing[d], 'x' + d);
              if (_nn[d] < nneighbors[d]) {
                  warn("In stencil_t requested nn=%i but use nn=%i for %c-direction",
                                  nneighbors[d], _nn[d], 'x' + d);
              }
          } // d spatial direction
          scale_coefficients(scale_factor);
      } // _constructor

    public:

      stencil_t(double const grid_spacing[3], int const nneighbors[3], double const scale_factor=1) {
          _constructor(grid_spacing, nneighbors, scale_factor);
      } // preferred constructor

      stencil_t(double const grid_spacing[3], int const nn=4, double const scale_factor=1) {
          int const nns[3] = {nn, nn, nn};
          _constructor(grid_spacing, nns, scale_factor);
      } // isotropic nn constructor

      stencil_t(double const h=1, int const nn=4, double const scale_factor=1) {
          double const hgs[3] = {h, h, h};
          int const nns[3] = {nn, nn, nn};
          _constructor(hgs, nns, scale_factor);
      } // isotropic constructor, default constructor

      double clear_diagonal_elements() { // modifies the coefficients c2nd[][]
          double diag{0};
          for (int d = 0; d < 3; ++d) {
              diag += c2nd[d][0];
              c2nd[d][0] = 0; // clear diagonal elements
          } // d
          return diag;
      } // clear_diagonal_elements

    public:
      void scale_coefficients(double const f[3]) {
          for (int d = 0; d < 3; ++d) {
              for (int i = 0; i < nnArraySize; ++i) {
                  c2nd[d][i] *= f[d];
              } // i
          } // d
      } // scale_coefficients

      void scale_coefficients(double const f) { double const f3[] = {f, f, f}; scale_coefficients(f3); }
    public:
      int8_t const * nearest_neighbors() const { return _nn; }
      int nearest_neighbors(int const d) const { assert(d >= 0); assert(d < 3); return _nn[d]; }
      real_t * Laplace_coefficients(int const d) { assert(d >= 0); assert(d < 3); return c2nd[d]; }

  }; // class stencil_t




  template <typename complex_out_t // result is stored in this precision
           ,typename complex_in_t // input comes in this precision
           ,typename real_fd_t> // computations are executed in this precision
  status_t apply(
        complex_out_t out[]
      , complex_in_t const in[]
      , real_space::grid_t const & g
      , stencil_t<real_fd_t> const & fd
      , double const factor=1
      , complex_in_t const boundary_phase[3][2]=nullptr
  ) {

      int const n16 = nnArraySize; // max number of finite difference neighbors, typically 16
      typedef int16_t list_integ_t;
      std::vector<list_integ_t> list[3]; // can be of type int16_t
      std::vector<complex_in_t> phas[3];
      for (int d = 0; d < 3; ++d) {
          int const n = g[d];
          assert(n >= 0);
          assert(n <= std::numeric_limits<list_integ_t>::max());
          // ToDo: check that n is smaller than the upper limit of int
          int const bc = g.boundary_condition(d);
          int const nf = fd.nearest_neighbors(d);
          if (nf > n) error("finite-difference range (%d) in %c-direction is larger than grid (%d grid points)", nf, 'x' + d, n);
          assert(nf <= n);
          assert(nf <= n16);
          int const nh = n16 + n + n16; // number including largest halos
          list[d] = std::vector<list_integ_t>(nh, -1); // get memory, init as -1:non-existing
          phas[d] = std::vector<complex_in_t>(nh, 0); // get memory, init neutral

          // core region
          for (int j = 0; j < n; ++j) {
              list[d][n16 + j] = j;
              phas[d][n16 + j] = 1;
          } // j

          complex_in_t const phase_low = boundary_phase ? boundary_phase[d][0] : 1;
          complex_in_t const phase_upp = boundary_phase ? boundary_phase[d][1] : 1;

          // lower boundary
          if (Periodic_Boundary == bc) { // periodic BC
              for (int j = -nf; j < 0; ++j) {
                  list[d][n16 + j] = (n + j) % n; // wrap around
                  phas[d][n16 + j] = phase_low; // incorrect if nf > n
              } // j
          } else if (Mirrored_Boundary == bc) { // mirror BC
              for (int j = -nf; j < 0; ++j) {
                  list[d][n16 + j] = - 1 - j; // mirror at -1 | 0
                  phas[d][n16 + j] = phase_low; // incorrect if nf > n
              } // j
          } // else open BC, list[:] = -1

          // upper boundary
          if (Periodic_Boundary == bc) { // periodic BC
              for (int j = 0; j < nf; ++j) {
                  list[d][n16 + n + j] = (n + j) % n; // wrap around
                  phas[d][n16 + n + j] = phase_upp; // incorrect if nf > n
              } // j
          } else if (Mirrored_Boundary == bc) { // mirror BC
              for (int j = 0; j < nf; ++j) {
                  list[d][n16 + n + j] = n - 1 - j; // mirror at n-1 | n
                  phas[d][n16 + n + j] = phase_upp; // incorrect if nf > n
              } // j
          } // else open BC, list[:] = -1

          if (0) { // DEBUG: show indirection list and phase factors
              std::printf("# indirection list for %c-direction ", 'x'+d);
              for (int j = -nf; j < n + nf; ++j) {
                  if (0 == j || n == j) std::printf(" |");
                  std::printf(" %i", list[d][n16 + j]);
              } // j
              std::printf("\n");
              std::printf("# phase factor list for %c-direction ", 'x'+d);
              for (int j = -nf; j < n + nf; ++j) {
                  if (0 == j || n == j) std::printf(" |");
                  std::printf("  %g %g", std::real(phas[d][n16 + j]), std::imag(phas[d][n16 + j]));
              } // j
              std::printf("\n");
          } // show indirection list

      } // spatial direction d

      real_fd_t const scale_factor = factor;
      for (int z = 0; z < g('z'); ++z) {
          for (int y = 0; y < g('y'); ++y) {
              for (int x = 0; x < g('x'); ++x) {

                  complex_out_t t(0); // init result

                  for (int d = 0; d < 3; ++d) {
                      int const nf = fd.nearest_neighbors(d);
                      int zyx[3] = {x, y, z};
                      int const i_center = zyx[d];
                      for (int jmi = -nf; jmi <= nf; ++jmi) {
                          int const j = i_center + jmi;
                          int const index = list[d][n16 + j];
                          if (index >= 0) {
                              zyx[d] = index;
#ifdef    GENERAL_CELL
                              // allow shift-rectangular cells from lower triangular cell matrices
                              if (1 == d) { // derive in y-direction
                                  auto const jy = int(j < 0) - int(j >= g('y'));
                                  zyx[0] = (zyx[0] + jy*g.shift_yx + 9*g('x')) % g('x');
                              } else // 'y'
                              if (2 == d) { // derive in z-direction
                                  auto const jz = int(j < 0) - int(j >= g('z'));
                                  zyx[0] = (zyx[0] + jz*g.shift_zx + 9*g('x')) % g('x');
                                  zyx[1] = (zyx[1] + jz*g.shift_zy + 9*g('y')) % g('y');
                              } // 'z'
                              assert(zyx[0] >= 0); assert(zyx[0] < g('x'));
                              assert(zyx[1] >= 0); assert(zyx[1] < g('y'));
                              assert(zyx[2] >= 0); assert(zyx[2] < g('z'));
#endif // GENERAL_CELL
                              int const jzyx = (zyx[2]*g('y') + zyx[1])*g('x') + zyx[0];
                              auto const coeff = fd.c2nd[d][std::abs(jmi)];
                              t += (phas[d][n16 + j] * in[jzyx]) * coeff;
                          } // index exists
                      } // jmi
                  } // d direction of the derivative

                  int const izyx = (z*g('y') + y)*g('x') + x;
                  out[izyx] = t * scale_factor; // store

              } // x
          } // y
      } // z

      return 0; // success
  } // apply









#ifdef    NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

  inline status_t test_create_and_destroy(int const echo=9) {
      auto const f = new stencil_t<float>();
      f->~stencil_t();
      stencil_t<double> d;
      return 0;
  } // test_create_and_destroy

  template <typename real_t>
  inline status_t test_Laplacian(int const echo=3) {
      status_t stat(0);
      double const h[3] = {1, 1, 1}; // unit grid spacings
      for (int dir = 0; dir < 3; ++dir) {
          int nn[3] = {0,0,0}; nn[dir] = 12; // switch FD off for the two perpendicular directions
          stencil_t<real_t> Laplacian(h, nn);
          int dims[] = {1,1,1}; dims[dir] = 127 + dir;
          real_space::grid_t g(dims);
          g.set_boundary_conditions(Periodic_Boundary);
          double const k = (1 + dir)*2*constants::pi/g[dir]; // wave vector of a single plane wave
          std::vector<real_t> values(g.all()), result(g.all());
          for (size_t i = 0; i < g.all(); ++i) values[i] = std::cos(k*i); // fill with some non-zero values
          stat += finite_difference::apply(result.data(), values.data(), g, Laplacian);
          if (echo > 5) std::printf("\n# in, result, ref values:\n");
          double dev{0};
          for (size_t i = 0; i < g.all(); ++i) {
              auto const ref = -k*k*values[i]; // analytic solution to the Laplacian operator applied to a plane wave
              if (echo > 5) std::printf("%ld %g %g %g\n", i, values[i], result[i], ref);
              // compare in the middle range result and ref values
              dev += std::abs(result[i] - ref);
          } // i
          if (echo > 2) std::printf("# %s %c-direction: dev = %g\n", __func__, 'x'+dir, dev);
      } // direction
      return stat;
  } // test_Laplacian

  template <typename real_t>
  inline status_t test_Bloch_wave(int const echo=3) {
      status_t stat(0);
      double const h[3] = {1, 1, 1}; // unit grid spacings
      std::complex<real_t> boundary_phase[3][2] = {{-1,-1}, {-1,-1}, {-1,-1}};
      double maxdev{0};
      for (int dir = 0; dir < 3; ++dir) {
          int nn[3] = {0,0,0}; nn[dir] = 12; // switch FD off for the two perpendicular directions
          stencil_t<real_t> Laplacian(h, nn);
          int dims[] = {1,1,1}; dims[dir] = 127 + dir;
          real_space::grid_t g(dims);
          g.set_boundary_conditions(Periodic_Boundary);
          std::vector<std::complex<real_t>> values(g.all()), result(g.all());
          for (int iphase = 0; iphase <= 180; iphase += 20) {
              double const k = (1 + dir + iphase/360.)*2*constants::pi/g[dir]; // wave vector of a single plane wave
              double const arg = iphase*constants::pi/180.;
              boundary_phase[dir][0] = std::complex<real_t>(std::cos(arg), -std::sin(arg));
              boundary_phase[dir][1] = real_t(1)/boundary_phase[dir][0];
              for (size_t i = 0; i < g.all(); ++i) {
                  values[i] = std::complex<real_t>(std::cos(k*i), std::sin(k*i));
              } // i
              stat += finite_difference::apply(result.data(), values.data(), g, Laplacian, 1, boundary_phase);
              double dev{0};
              for (size_t i = 0; i < g.all(); ++i) {
                  auto const ref = -k*k*values[i]; // analytic solution to the Laplacian operator applied to a plane wave
                  // compare in the middle range result and ref values
                  dev += std::abs(result[i] - ref);
              } // i
              if (echo > 3) std::printf("# %s %c-direction: dev = %g\n", __func__, 'x'+dir, dev);
              maxdev = std::max(maxdev, std::abs(dev));
          } // iphase
      } // direction
      if (echo > 1) std::printf("\n# %s largest deviation is %.1e\n", __func__, maxdev);
      return stat;
  } // test_Bloch_wave

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_create_and_destroy(echo);
      stat += test_Laplacian<double>(echo);
      stat += test_Bloch_wave<double>(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace finite_difference
