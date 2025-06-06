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

  template <typename real_t>
  double polar(std::complex<real_t> const c) { return std::atan2(c.imag(), c.real())*(180./constants::pi); }
  template <typename real_t> double polar(real_t const r) { return (r < 0)*180; }

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
          // check that n is smaller than the upper limit of int
          assert(0 < n); assert(n <= std::numeric_limits<list_integ_t>::max());
          int const bc = g.boundary_condition(d);
          int const nf = fd.nearest_neighbors(d);
          if (nf > n) error("finite-difference range (%d) in %c-direction is larger than grid (%d grid points)", nf, 'x' + d, n);
          assert(nf <= n);
          assert(nf <= n16);
          int const nh = n16 + n + n16; // number including largest halos
          list[d] = std::vector<list_integ_t>(nh, -1); // get memory, init as -1:non-existing
          phas[d] = std::vector<complex_in_t>(nh, 0); // get memory, init zero

          // core region
          for (int j = 0; j < n; ++j) {
              list[d][n16 + j] = j; // identity
              phas[d][n16 + j] = 1; // neutral
          } // j

          complex_in_t const phase_low = boundary_phase ? boundary_phase[d][0] : 1;
          complex_in_t const phase_upp = boundary_phase ? boundary_phase[d][1] : 1;

          // lower boundary
          if (Periodic_Boundary == bc || Shifted_Boundary == bc) { // periodic BC
              for (int j = -nf; j < 0; ++j) {
                  list[d][n16 + j] = (n + j) % n; // wrap around
                  phas[d][n16 + j] = phase_low; // incorrect if nf > n
              } // j
          } else if (Mirrored_Boundary == bc) { // mirror BC
              for (int j = -nf; j < 0; ++j) {
                  list[d][n16 + j] = - 1 - j; // mirror at -1 | 0
                  phas[d][n16 + j] = phase_low; // incorrect if nf > n
              } // j
          } // else open BC, list[:] = -1, phas[:] = 0

          // upper boundary
          if (Periodic_Boundary == bc || Shifted_Boundary == bc) { // periodic BC
              for (int j = 0; j < nf; ++j) {
                  list[d][n16 + n + j] = (n + j) % n; // wrap around
                  phas[d][n16 + n + j] = phase_upp; // incorrect if nf > n
              } // j
          } else if (Mirrored_Boundary == bc) { // mirror BC
              for (int j = 0; j < nf; ++j) {
                  list[d][n16 + n + j] = n - 1 - j; // mirror at n-1 | n
                  phas[d][n16 + n + j] = phase_upp; // incorrect if nf > n
              } // j
          } // else open BC, list[:] = -1, phas[:] = 0

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

      auto const nx = g('x'), ny = g('y'), nz = g('z');
#ifdef    GENERAL_CELL
      complex_in_t const phase_xy_low = std::pow(boundary_phase ? boundary_phase[0][0] : 1, complex_in_t(g.shift_yx/double(nx))),
                         phase_xy_upp = std::pow(boundary_phase ? boundary_phase[0][1] : 1, complex_in_t(g.shift_yx/double(nx)));
      if (0) { std::printf("# %s: low= %g %g = %g degrees, upp= %g %g = %g degrees\n", __func__,
                            std::real(phase_xy_low), std::imag(phase_xy_low), polar(phase_xy_low),
                            std::real(phase_xy_upp), std::imag(phase_xy_upp), polar(phase_xy_upp)); }
      assert(0 <= g.shift_yx); assert(g.shift_yx < nx);

      complex_in_t const phase_xz_low = std::pow(boundary_phase ? boundary_phase[0][0] : 1, complex_in_t(g.shift_zx/double(nx))),
                         phase_xz_upp = std::pow(boundary_phase ? boundary_phase[0][1] : 1, complex_in_t(g.shift_zx/double(nx)));
      complex_in_t const phase_yz_low = std::pow(boundary_phase ? boundary_phase[1][0] : 1, complex_in_t(g.shift_zy/double(ny))),
                         phase_yz_upp = std::pow(boundary_phase ? boundary_phase[1][1] : 1, complex_in_t(g.shift_zy/double(ny)));
#endif // GENERAL_CELL

      real_fd_t const scale_factor = factor;
      for (int z = 0; z < nz; ++z) {
          for (int y = 0; y < ny; ++y) {
              for (int x = 0; x < nx; ++x) {

                  complex_out_t t(0); // init result

                  for (int d = 0; d < 3; ++d) {
                      int const nf = fd.nearest_neighbors(d);
                      int zyx[3] = {x, y, z};
                      int const i_center = zyx[d];
                      for (int jmi = -nf; jmi <= nf; ++jmi) {
                          int const j = i_center + jmi;
                          int const index = list[d][n16 + j]; // indirection to capture boundary
                          if (index >= 0) {
                              zyx[d] = index;
                              auto phase = phas[d][n16 + j]; // phase factor when crossing a periodic boundary
#ifdef    GENERAL_CELL
                              // ToDo: put this part into a section which is taken out in regular versions by a boolean template parameter

                              // allow shift-rectangular cells from lower triangular cell matrices
                              if (1 == d) { // derive in y-direction

                                    if (j >= ny) {
                                        zyx[0] = (x - g.shift_yx + nx) % nx; // modify the x-coordinate of the source
                                        phase *= phase_xy_upp;
                                        if (zyx[0] > x) { phase *= phas[0][n16 - 1]; }
                                    } else
                                    if (j < 0) {
                                        zyx[0] = (x + g.shift_yx + nx) % nx; // modify the x-coordinate of the source
                                        phase *= phase_xy_low;
                                        if (zyx[0] < x) { phase *= phas[0][n16 + nx]; } 
                                    } else {
                                        zyx[0] = x;
                                    }

                              } else // 'y'
                              if (2 == d) { // derive in z-direction

                                    if (j >= nz) {
                                        zyx[0] = (x - g.shift_zx + nx) % nx; // modify the x-coordinate of the source
                                        phase *= phase_xz_upp;
                                        if (zyx[0] > x) { phase *= phas[0][n16 - 1]; }
                                    } else
                                    if (j < 0) {
                                        zyx[0] = (x + g.shift_zx + nx) % nx; // modify the x-coordinate of the source
                                        phase *= phase_xz_low;
                                        if (zyx[0] < x) { phase *= phas[0][n16 + nx]; } 
                                    } else {
                                        zyx[0] = x;
                                    }

                                    if (j >= nz) {
                                        zyx[1] = (y - g.shift_zy + ny) % ny; // modify the y-coordinate of the source
                                        phase *= phase_yz_upp;
                                        if (zyx[1] > y) { phase *= phas[1][n16 - 1]; }
                                    } else
                                    if (j < 0) {
                                        zyx[1] = (y + g.shift_zy + ny) % ny; // modify the y-coordinate of the source
                                        phase *= phase_yz_low;
                                        if (zyx[1] < y) { phase *= phas[1][n16 + ny]; } 
                                    } else {
                                        zyx[1] = y;
                                    }

                              } // 'z'
                              assert(0 <= zyx[0]); assert(zyx[0] < nx);
                              assert(0 <= zyx[1]); assert(zyx[1] < ny);
                              assert(0 <= zyx[2]); assert(zyx[2] < nz);
#endif // GENERAL_CELL
                              int const jzyx = (zyx[2]*ny + zyx[1])*nx + zyx[0]; // source index
                              auto const coeff = fd.c2nd[d][std::abs(jmi)];
                              auto const contrib = phase * in[jzyx];
                            //   if (29 == x && 0 == y) { std::printf("# FD add %6.1f + %6.1f [%9.3f] = %6.1f \t from [%i %i %i]\n",
                            //                 polar(phase), polar(in[jzyx]), coeff, polar(contrib), zyx[2], zyx[1], zyx[0]); }
                              t += contrib*coeff;
                          } // index exists
                      } // jmi
                  } // d direction of the derivative

                  int const izyx = (z*ny + y)*nx + x;
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

#ifdef    GENERAL_CELL
    template <typename real_t>
    inline status_t test_general_cell(int const echo=3) {
        // ToDo: not covered: more than one shift non-zero at the same time
        status_t stat(0);
        // create a shifted Cartesian unit cell with angles 60 degree
        int const dims[] = {15, 13, 12};
        double const h[3] = {1./dims[0], std::sqrt(3/4.)/dims[1], std::sqrt(2/3.)/dims[2]}; // grid spacings
        int const nn[3] = {12, 12, 12}; // FD-order
        stencil_t<real_t> Laplacian(h, nn);
        real_space::grid_t g(dims);
        std::vector<std::complex<real_t>> values(g.all()), result(g.all());
        std::complex<real_t> boundary_phase[3][2];
        g.set_boundary_conditions(Periodic_Boundary, Shifted_Boundary, Shifted_Boundary);
        if (echo > 1) std::printf("\n# %s start\n", __func__);
                                    //  000        001   010     011    100    101      110     111
        char const shift_name[8][8] = {"noshift", "xy", "xz", "xz+yz", "yz", "xy+yz", "xy+xz", "xyxzyz"};
        double maxdevall{0};
        for (int is{0}; is < 7; ++is) { // is==7 needs a very long time, test only 0--6 to be faster
            int const shift_dims[] = {((is >> 0) & 0x1)*(dims[0] - 1),      // max shift along x-direction on crossing a y-boundary
                                      ((is >> 1) & 0x1)*(dims[0] - 1),      // max shift along x-direction on crossing a z-boundary
                                      ((is >> 2) & 0x1)*(dims[1] - 1)};     // max shift along y-direction on crossing a z-boundary
            if (echo > 3) { std::printf("# %s: test xy_shift in [0, %d], xz_shift in [0, %d], yz_shift in [0, %d], name= %s\n",
                                __func__, shift_dims[0], shift_dims[1], shift_dims[2], shift_name[is]); }
            double maxdev_is{0};
            for (int idirection{0}; idirection <= 90; idirection += 15) { // angle
                // prepare plane wave vector (this could be moved outside the shift loops to save time)
                double const k = 1.6; // sqRy
                auto constexpr arc = constants::pi/180.;
                auto const k_cos = k*std::cos(idirection*arc), k_sin = k*std::sin(idirection*arc);
                double const kv_dir[] = {k_cos, k_sin, .1, k_cos, k_sin, .1, k_cos, k_sin, .1, k_cos};
                double const *const kv = kv_dir + is; // only reference *kv as kv[3]
                // prepare wave values
                for (int iz{0}; iz < g('z'); ++iz) {
                    for (int iy{0}; iy < g('y'); ++iy) {
                        for (int ix{0}; ix < g('x'); ++ix) {
                            auto const arg = kv[0]*ix*h[0] + kv[1]*iy*h[1] + kv[2]*iz*h[2]; // prepare an arbitrary plane wave
                            values[(iz*g('y') + iy)*g('x') + ix] = std::complex<real_t>(std::cos(arg), std::sin(arg));
                        } // ix
                    } // iy
                } // iz
                // prepare matching boundary phases
                for (int dir{0}; dir < 3; ++dir) {
                    auto const arg = kv[dir]*g[dir]*h[dir];
                    boundary_phase[dir][1] = std::complex<real_t>(std::cos(arg), std::sin(arg));
                    boundary_phase[dir][0] = real_t(1)/boundary_phase[dir][1];
                    // if (echo > 11) { std::printf("# %s(%d deg) %c-phase= %g %g\n", __func__, idirection, 'x'+dir, boundary_phase[dir][1].real(), boundary_phase[dir][1].imag()); }
                } // dir

        double maxdev{0};
        for (int xy_shift{0}; xy_shift <= shift_dims[0]; ++xy_shift) {          // shift along x-direction on crossing a y-boundary
        for (int xz_shift{0}; xz_shift <= shift_dims[1]; ++xz_shift) {          // shift along x-direction on crossing a z-boundary
        for (int yz_shift{0}; yz_shift <= shift_dims[2]; ++yz_shift) {          // shift along y-direction on crossing a z-boundary
            if (echo > 13) { std::printf("# %s: test %s-shifts\n", __func__, shift_name[is]); }
            double const cell_shape[3][4] = {{h[0]*dims[0],             0,             0, 0},
                                             {h[0]*xy_shift, h[1]*dims[1],             0, 0},  // lower triangular matrix
                                             {h[0]*xz_shift, h[1]*yz_shift, h[2]*dims[2], 0}};
            g.set_cell_shape(cell_shape, echo/4);

            // apply
            stat += finite_difference::apply(result.data(), values.data(), g, Laplacian, 1, boundary_phase);

            // compare to anaytical Laplacian
            auto const k2 = pow2(kv[0]) + pow2(kv[1]) + pow2(kv[2]);
            double dev{0};
            for (size_t i{0}; i < g.all(); ++i) {
                auto const val = values[i], res = result[i], ref = -k2*val; // reference is the analytic solution to the Laplacian operator applied to a plane wave
                auto const dev_add = std::abs(res - ref);
                dev += dev_add;
            } // i
            if (echo > 9) { std::printf("# %s direction=%4d degrees xy= %d/%d, xz= %d/%d, yz= %d/%d, dev= %g\n", __func__,
                                  idirection,  xy_shift, dims[0],  xz_shift, dims[0],  yz_shift, dims[1],  dev/g.all()); }
            maxdev = std::max(maxdev, std::abs(dev/g.all()));
        }}} // *_shift

        stat += (maxdev > 1e-12);
        if (echo > 4) std::printf("# %s direction=%4d degrees, dev= %g\n", __func__, idirection, maxdev);
        maxdev_is = std::max(maxdev_is, maxdev);
    } // idirection

        if (echo > 0) { std::printf("# %s largest deviation tests \'%s\' is %.1e\n", __func__, shift_name[is], maxdev_is); }
        maxdevall = std::max(maxdevall, maxdev_is);
        } // is
        if (echo > 0) { std::printf("# %s largest deviation of all tests is %.1e\n", __func__, maxdevall); }
        return stat;
    } // test_general_cell
#endif // GENERAL_CELL

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
#ifdef    GENERAL_CELL
      stat += test_general_cell<double>(echo);
#endif // GENERAL_CELL
      stat += test_create_and_destroy(echo);
      stat += test_Laplacian<double>(echo);
      stat += test_Bloch_wave<double>(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace finite_difference
