#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // uint32_t, int8_t
#include <algorithm> // std::min, ::max
#include <cmath> // std::floor, ::ceil, ::sqrt, ::abs
#include <cassert> // assert

#include "status.hxx" // status_t

#include "inline_math.hxx" // set, scale
#include "simple_math.hxx" // ::determinant
#include "constants.hxx" // ::pi
#include "display_units.h" // Ang, _Ang
#include "recorded_warnings.hxx" // warn
#include "boundary_condition.hxx" // Periodic_Boundary, Isolated_Boundary, Mirrored_Boundary, Invalid_Boundary

namespace real_space {

  int constexpr debug = 0;

  double inline length2(double const v[3]) { return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; }
  double inline length(double const v[3]) { return std::sqrt(length2(v)); }
  double inline angle(double const v[3], double const w[3]) {
      double const ll = std::sqrt(length2(v)*length2(w));
      double const dot = v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
      return (ll > 0) ? std::acos(std::min(std::max(-1., dot/ll), 1.)) : 0;
  } // angle

  class grid_t {
  private:
      uint32_t dims[4]; // 0,1,2:real-space grid dimensions, 3:outer dim
      int8_t bc[3]; // boundary conditions
  public:
      double h[3], inv_h[3]; // grid spacings and their inverse
      double cell[3][4]; // cell shape (only [3][3] used, only diagonal elements {cell[0][0], cell[1][1], cell[2][2]} if is_Cartesian)
#ifdef    GENERAL_CELL
      int32_t shift_yx, shift_zx, shift_zy;

      status_t correct_shift_cell_parameters(int32_t & n_shift_yx, int const y, int const x, int const echo=0) {
          double constexpr threshold = 1e-6;
          double const old_cell_param = cell[y][x];
          n_shift_yx = std::round(old_cell_param*inv_h[x]);
          double const new_cell_param = n_shift_yx*h[x];
          if (echo > 6) std::printf("# shift_%c%c=%d or %6.3f %%\n", 'x'+y, 'x'+x, n_shift_yx, n_shift_yx/(dims[x]*.01));

          double const dev = old_cell_param - new_cell_param;
          status_t stat(0);
          if (std::abs(dev) > threshold*cell[x][x]) {
              warn("inaccurate shift_%c%c: %g - %d*%g = %g %s", 'x'+y, 'x'+x, old_cell_param*Ang, n_shift_yx, h[x]*Ang, dev*Ang, _Ang);
              ++stat;
          } else if (echo > 8) std::printf("# shift_%c%c=%d, relative deviation is %.1e\n", 'x'+y, 'x'+x, n_shift_yx, dev/cell[x][x]);
          if (std::abs(n_shift_yx) >= dims[x]) {
              error("May not shift more than one cell on perpendicular translation in %c-direction!", 'x'+y); // avoid problems with periodic images
          }
          cell[y][x] = new_cell_param;
          if (Periodic_Boundary != bc[x]|| Periodic_Boundary != bc[y]) {
              warn("for shift_%c%c=%d boundary conditions must be periodic, found bc= %d and %d", 'x'+y, 'x'+x, n_shift_yx, bc[x], bc[y]);
              ++stat;
          } // boundary is not periodic
          return stat;
      } // correct_shift_cell_parameters
#endif // GENERAL_CELL

      grid_t(void) : dims{0,0,0,1}, bc{0,0,0}, h{1,1,1}, inv_h{1,1,1} {
          set(cell[0], 12, 0.0);
#ifdef    GENERAL_CELL
          shift_yx = shift_zx = shift_zy = 0;
#endif // GENERAL_CELL
      } // default constructor

      grid_t(int const d0, int const d1, int const d2, int const dim_outer=1)
        : bc{0,0,0}, h{1,1,1}, inv_h{1,1,1} {
          dims[0] = std::max(1, d0); // x
          dims[1] = std::max(1, d1); // y
          dims[2] = std::max(1, d2); // z
          dims[3] = std::max(1, dim_outer);
          long const nnumbers = dims[3] * dims[2] * dims[1] * dims[0];
          if (nnumbers > 0) {
              if (debug) std::printf("# grid with %d x %d x %d * %d = %.6f M numbers\n", 
                  dims[0], dims[1], dims[2], dims[3], nnumbers*1e-6);
          } else {
              if (debug) std::printf("# grid invalid: dims={%d, %d, %d,  %d}\n",
                  dims[0], dims[1], dims[2], dims[3]);
          }
          set(cell[0], 12, 0.0);
          for (int d = 0; d < 3; ++d) {
              cell[d][d] = dims[d]*h[d];
          } // d
#ifdef    GENERAL_CELL
          shift_yx = shift_zx = shift_zy = 0;
#endif // GENERAL_CELL
      } // constructor

      template <typename int_t>
      grid_t(int_t const dim[3], int const dim_outer=1) : grid_t(dim[0], dim[1], dim[2], dim_outer) {} // delegating contructor

      ~grid_t() {
          if (debug) std::printf("# release a grid with %d x %d x %d * %d = %.6f M numbers\n",
              dims[0], dims[1], dims[2], dims[3], 1e-6*dims[3]*dims[2]*dims[1]*dims[0]);
      } // destructor

      status_t set_cell_shape(double const general_cell[3][4], int const echo=0, double const factor=1) {
          if (nullptr == general_cell) return -1;
          status_t stat(0);
          for (int d = 0; d < 3; ++d) {
              set(cell[d], 3, general_cell[d], factor); 
              cell[d][3] = 0; // not used, only for memory alignment
          } // d
          auto const c0 = cell[0], c1 = cell[1], c2 = cell[2];
          if (dims[0] > 0 && dims[1] > 0 && dims[2] > 0) {
              if (is_Cartesian()) {
                  stat += set_grid_spacing(cell[0][0]/dims[0], cell[1][1]/dims[1], cell[2][2]/dims[2], echo);
              } else 
#ifdef    GENERAL_CELL
              if (is_shifted()) { // shifted
                  if (echo > 3) std::printf("# create shifted cell with  %d %d %d  grid points\n", dims[0], dims[1], dims[2]);
                  stat += set_grid_spacing(cell[0][0]/dims[0], cell[1][1]/dims[1], cell[2][2]/dims[2], echo);
                  if (echo > 3) std::printf("# grid spacings  %g %g %g %s\n", h[0]*Ang, h[1]*Ang, h[2]*Ang, _Ang);
                  stat += correct_shift_cell_parameters(shift_yx, 1, 0, echo);
                  stat += correct_shift_cell_parameters(shift_zx, 2, 0, echo);
                  stat += correct_shift_cell_parameters(shift_zy, 2, 1, echo);
                  // show the lengths and angles of unit vectors
                  if (echo > 3) std::printf("# shifted cell vector lengths  %g %g %g  %s\n",
                                    length(c0)*Ang, length(c1)*Ang, length(c2)*Ang, _Ang);
                  double constexpr deg = 180/constants::pi;
                  if (echo > 3) std::printf("# cell vector angles  %g %g %g  degrees\n",
                                    angle(c1, c2)*deg, angle(c2, c0)*deg, angle(c0, c1)*deg);
              } else
#endif // GENERAL_CELL
              {
                  stat += set_grid_spacing(length(c0)/dims[0], length(c1)/dims[1], length(c2)/dims[2], echo);
              }
          } else {
              if (echo > 2) std::printf("# cannot set grid spacing as grid dims are %d %d %d\n", dims[0], dims[1], dims[2]);
          }
          if (echo > 4) std::printf("# cell shape %g %g %g  %g %g %g  %g %g %g %s, type=%s\n",
                            c0[0]*Ang, c0[1]*Ang, c0[2]*Ang,
                            c1[0]*Ang, c1[1]*Ang, c1[2]*Ang,
                            c2[0]*Ang, c2[1]*Ang, c2[2]*Ang, _Ang,
                            is_Cartesian()?"Cartesian":(has_upper_elements()?"general":"shifted"));
          return stat;
      } // set_cell_shape

      status_t set_grid_spacing(double const hx, double const hy=-1, double const hz=-1, int const echo=0) {
          status_t stat(0);
          double const h3[3] = {hx, (hy<0)?hx:hy, (hz<0)?hx:hz};
          for (int i3 = 0; i3 < 3; ++i3) {
              h[i3] = h3[i3]; // convert to double
              if (h[i3] > 0) {
                  inv_h[i3] = 1./h[i3]; // invert only here
              } else {
                  ++stat;
              } // h > 0
          } // i3
          return stat;
      } // set

      status_t set_boundary_conditions(int8_t const bc3[3]) { set(bc, 3, bc3); return 0; }
      status_t set_boundary_conditions(int8_t const bcx
                                     , int8_t const bcy=Invalid_Boundary 
                                     , int8_t const bcz=Invalid_Boundary) {
          bc[0] = bcx;
          bc[1] = (bcy == Invalid_Boundary) ? bcx : bcy;
          bc[2] = (bcz == Invalid_Boundary) ? bcx : bcz;
          return  (bcx == Invalid_Boundary);
      } // set

      inline int has_upper_elements() const {
            return int(0 != cell[0][1]) + int(0 != cell[0][2]) + int(0 != cell[1][2]);
      } // has_upper_elements

      inline int has_lower_elements() const {
            return int(0 != cell[1][0]) + int(0 != cell[2][0]) + int(0 != cell[2][1]);
      } // has_lower_elements

      inline int is_Cartesian() const { // diagonal elements must be positive, off-diagonals zero
            return (cell[0][0] > 0) && (cell[1][1] > 0) && (cell[2][2] > 0) &&
                (0 == has_lower_elements()) && (0 == has_upper_elements());
      } // is_Cartesian

      inline int is_shifted(int const including_Cartesian=0) const {
            return (cell[0][0] > 0) && (cell[1][1] > 0) && (cell[2][2] > 0) &&
                (has_lower_elements() >= including_Cartesian) && (0 == has_upper_elements());
      } // is_shifted
      // is_shifted(1) --> shifted but not Cartesian
      // is_shifted(0) --> shifted, can be Cartesian

      inline double volume() const { return simple_math::determinant(
          cell[0][0], cell[0][1], cell[0][2],
          cell[1][0], cell[1][1], cell[1][2],
          cell[2][0], cell[2][1], cell[2][2]); }
      inline int operator[] (int const d) const { assert(0 <= d); assert(d < 4); return dims[d]; }
      inline int operator() (char const c) const { assert('x' <= (c|32)); assert((c|32) <= 'z'); return dims[(c|32) - 120]; }
      inline double dV() const { return is_shifted() ? h[0]*h[1]*h[2] : volume()/std::max(1u, dims[0]*dims[1]*dims[2]); } // volume element for integration
      inline double grid_spacing(int const d) const { assert(0 >= d); assert(d < 3); return h[d]; } // so far not used
      inline double const * grid_spacings() const { return h; } // so far not used
      inline double smallest_grid_spacing() const { return std::min(std::min(h[0], h[1]), h[2]); }
      inline size_t all() const { return ((size_t(dims[3]) * dims[2]) * dims[1]) * dims[0]; }
      inline int8_t boundary_condition(int  const d) const { assert(0 <= d); assert(d < 3); return bc[d]; }
      inline int8_t boundary_condition(char const c) const { return boundary_condition((c|32) - 120); }
      inline int8_t const * boundary_conditions() const { return bc; }
      inline uint32_t const * grid_points() const { return dims; }
      inline int number_of_boundary_conditions(int const bc_ref=Periodic_Boundary) const { 
                    return (bc_ref == bc[0]) + (bc_ref == bc[1]) + (bc_ref == bc[2]); };
  }; // class grid_t



  template <typename real_t>
  status_t add_function_general(
        real_t values[] // grid values which are modified
      , grid_t const & g // grid descriptor
      , double const r2coeff[] // coefficients of the radial function on r^2-grid
      , int const ncoeff // number of coefficients on the r^2-grid
      , float const hcoeff // r^2-grid parameter
      , double *added=nullptr // optional result: how much (e.g. charge) was added
      , double const center[3]=nullptr // spherical center w.r.t. the position of grid point (0,0,0), internal coords!
      , double const factor=1 // optional scaling
      , float const r_cut=-1 // radial truncation, -1:use the max radius of the r^2-grid
  ) {
      // Add a spherically symmetric regular function to a general grid.
      // The function is tabulated as r2coeff[0 <= hcoeff*r^2 < ncoeff]
      status_t stat(0);
      double c[3] = {0,0,0}; if (center) set(c, 3, center);
      double const r_max = std::sqrt((ncoeff - 1.)/hcoeff); // largest radius of the r^2-grid
      double const rcut = (-1 == r_cut) ? r_max : std::min(double(r_cut), r_max);
      double const r2cut = rcut*rcut;
      assert(hcoeff*r2cut < ncoeff);
      assert(hcoeff > 0);
      double added_charge{0}; // clear
      double const denom[] = {1./g[0], 1./g[1], 1./g[2]}; // grid denominators
      for (        int iz = 0; iz < g[2]; ++iz) {
          for (    int iy = 0; iy < g[1]; ++iy) {
              for (int ix = 0; ix < g[0]; ++ix) {
                  double const iv[] = {ix*denom[0], iy*denom[1], iz*denom[2]};
                  double rv[3];
                  for (int d = 0; d < 3; ++d) {
                      rv[d] = iv[0]*g.cell[0][d] + iv[1]*g.cell[1][d] + iv[2]*g.cell[2][d] - c[d];
                  } // d
                  double const r2 = pow2(rv[0]) + pow2(rv[1]) + pow2(rv[2]);
                  if (r2 < r2cut) {
                      int const ir2 = int(hcoeff*r2);
                      if (ir2 < ncoeff) {
                          int const izyx = (iz*g('y') + iy)*g('x') + ix;
                          double const w8 = hcoeff*r2 - ir2; // linear interpolation weight
                          int const ir2p1 = ir2 + 1;
                          auto const value_to_add = (r2coeff[ir2] * (1 - w8)
                               + ((ir2p1 < ncoeff) ? r2coeff[ir2p1] : 0)*w8);
                          values[izyx] += factor*value_to_add;
                          added_charge += factor*value_to_add;
                      } // ir2 < ncoeff
                  } // inside rcut
              } // ix
          } // iy
      } // iz
      if (added) *added = added_charge * g.dV(); // volume integral
      return stat;
  } // add_function_general




  template <typename real_t>
  status_t add_function(
        real_t values[] // grid values which are modified
      , grid_t const & g // grid descriptor
      , double const r2coeff[] // coefficients of the radial function on r^2-grid
      , int const ncoeff // number of coefficients on the r^2-grid
      , float const hcoeff // r^2-grid parameter
      , double *added=nullptr // optional result: how much (e.g. charge) was added
      , double const center[3]=nullptr // spherical center w.r.t. the position of grid point (0,0,0)
      , double const factor=1 // optional scaling
      , float const r_cut=-1 // radial truncation, -1:use the max radius of the r^2-grid
  ) {
      // Add a spherically symmetric regular function to the grid.
      // The function is tabulated as r2coeff[0 <= hcoeff*r^2 < ncoeff]
      status_t stat(0);
      double c[3] = {0,0,0}; if (center) set(c, 3, center);
      double const r_max = std::sqrt((ncoeff - 1.)/hcoeff); // largest radius of the r^2-grid
      double const rcut = (-1 == r_cut) ? r_max : std::min(double(r_cut), r_max);
      double const r2cut = rcut*rcut;
      assert(hcoeff*r2cut < ncoeff);
      if (!g.is_Cartesian()) {
          return add_function_general(values, g, r2coeff, ncoeff, hcoeff, added, center, factor, r_cut);
      } // not Cartesian
      int imn[3], imx[3];
#ifdef    DEBUG
      size_t nwindow{1};
      size_t modified{0};
#endif // DEBUG
      for (int d = 0; d < 3; ++d) {
          imn[d] = std::max(0, int(std::floor((c[d] - rcut)*g.inv_h[d])));
          imx[d] = std::min(   int(std::ceil ((c[d] + rcut)*g.inv_h[d])), g[d] - 1);
#ifdef    DEBUG
          std::printf("# %s window %c = %d elements from %d to %d\n", __func__, 'x'+d, imx[d] + 1 - imn[d], imn[d], imx[d]);
          nwindow *= std::max(0, imx[d] + 1 - imn[d]);
#endif // DEBUG
      } // d
      assert(hcoeff > 0);
      double added_charge{0}; // clear
      size_t out_of_range{0};
      for (            int iz = imn[2]; iz <= imx[2]; ++iz) {  double const vz = iz*g.h[2] - c[2], vz2 = vz*vz;
          for (        int iy = imn[1]; iy <= imx[1]; ++iy) {  double const vy = iy*g.h[1] - c[1], vy2 = vy*vy;
              if (vz2 + vy2 < r2cut) {
                  for (int ix = imn[0]; ix <= imx[0]; ++ix) {  double const vx = ix*g.h[0] - c[0], vx2 = vx*vx;
                      double const r2 = vz2 + vy2 + vx2;
                      if (r2 < r2cut) {
                          int const ir2 = int(hcoeff*r2);
                          if (ir2 < ncoeff) {
                              int const izyx = (iz*g('y') + iy)*g('x') + ix;
                              double const w8 = hcoeff*r2 - ir2; // linear interpolation weight
                              int const ir2p1 = ir2 + 1;
                              auto const value_to_add = (r2coeff[ir2] * (1 - w8)
                                   + ((ir2p1 < ncoeff) ? r2coeff[ir2p1] : 0)*w8);
                              values[izyx] += factor*value_to_add;
                              added_charge += factor*value_to_add;
#ifdef    DEBUG
                              ++modified;
#endif // DEBUG

#if 0
//        std::printf("#rs %g %g\n", std::sqrt(r2), value_to_add);
//        std::printf("#rs %.1f %.1f %.1f %.12f\n", vx*g.inv_h[0], vy*g.inv_h[1], vz*g.inv_h[2], value_to_add);
#endif // 0
                          } else {
                              ++out_of_range;
                          } // ir2 < ncoeff
                      } // inside rcut
                  } // ix
              } // rcut for (y,z)
          } // iy
      } // iz
      if (added) *added = added_charge * g.dV(); // volume integral
#ifdef    DEBUG
      std::printf("# %s modified %.3f k inside a window of %.3f k on a grid of %.3f k grid values.\n", 
              __func__, modified*1e-3, nwindow*1e-3, g('x')*g('y')*g('z')*1e-3); // show stats
#endif // DEBUG
      if (out_of_range > 0) {
          stat += 0 < warn("Found %ld entries out of range of the radial function!", out_of_range);
      } // out of range of the radial function
      return stat;
  } // add_function

  template <typename real_t>
  status_t Bessel_projection(
        double q_coeff[] // result Bessel coefficients
      , int const nq // number of coefficients
      , float const dq // wave number spacing in reciprocal space
      , real_t const values[] // grid array
      , grid_t const & g // grid descriptor
      , double const center[3]=nullptr // spherical center
      , double const factor=1 // scaling factor
      , float const r_cut=10 // radial truncation, default 10 Bohr
  ) {
      double c[3] = {0,0,0}; if (center) set(c, 3, center);
      double const rcut = r_cut, r2cut = rcut*rcut; // stop at 10 Bohr
      int imn[3], imx[3];
      for (int d = 0; d < 3; ++d) {
          imn[d] = std::max(0, int(std::floor((c[d] - rcut)*g.inv_h[d])));
          imx[d] = std::min(   int(std::ceil ((c[d] + rcut)*g.inv_h[d])), g[d] - 1);
// #ifdef    DEBUG
//           std::printf("# %s window %c = %d elements from %d to %d\n", __func__, 'x'+d, imx[d] + 1 - imn[d], imn[d], imx[d]);
// #endif // DEBUG
      } // d
      set(q_coeff, nq, 0.0); // clear
      for (            int iz = imn[2]; iz <= imx[2]; ++iz) {  double const vz = iz*g.h[2] - c[2], vz2 = vz*vz;
          for (        int iy = imn[1]; iy <= imx[1]; ++iy) {  double const vy = iy*g.h[1] - c[1], vy2 = vy*vy;
              if (vz2 + vy2 < r2cut) {
                  for (int ix = imn[0]; ix <= imx[0]; ++ix) {  double const vx = ix*g.h[0] - c[0], vx2 = vx*vx;
                      double const r2 = vz2 + vy2 + vx2;
                      if (r2 < r2cut) {
                          int const ixyz = (iz*g('y') + iy)*g('x') + ix;
                          double const r = std::sqrt(r2);
                          double const val = double(values[ixyz]);
//                        std::printf("%g %g\n", r, val); // DEBUG
                          for (int iq = 0; iq < nq; ++iq) {
                              double const q = iq*dq;
                              double const x = q*r;
                         //   double const j0 = bessel_transform::Bessel_j0(x);
                              double const j0 = (x*x < 1e-16) ? (1. - x*x*(1/6.)) : (std::sin(x)/x);
                              q_coeff[iq] += val * j0;
                          } // iq
                      } // inside rcut
                  } // ix
//                std::printf("\n"); // DEBUG
              } // rcut for (y,z)
          } // iy
      } // iz
      double const sqrt2pi = std::sqrt(2./constants::pi); // this makes the transform symmetric
      scale(q_coeff, nq, g.dV()*factor*sqrt2pi); // volume element, external factor, Bessel transform factor
      return 0; // success
  } // Bessel_projection






#ifdef    NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

  inline status_t test_create_and_destroy(int const echo=9) {
      int const dims[] = {10, 20, 30};
      auto gp = new grid_t(dims); // construct
      gp->~grid_t(); // explicity call the destructor
      return 0;
  } // test_create_and_destroy

  inline status_t test_add_function(int const echo=9) {
      if (echo > 0) std::printf("\n# %s\n", __func__);
      int const dims[] = {32, 31, 30};
      grid_t g(dims);
      g.set_grid_spacing(0.333);
      double const cnt[] = {g[0]*.42*g.h[0],
                            g[1]*.51*g.h[1],
                            g[2]*.60*g.h[2]}; // center is slightly shifted from exact grid point positions
      int const nr2 = 1 << 11;
      float const rcut = 4, inv_hr2 = nr2/(rcut*rcut);
      double const hr2 = 1./inv_hr2;
      double r2c[nr2], rad_integral{0};
      if (echo > 4) std::printf("\n# values on the radial grid\n");
      for (int ir2 = 0; ir2 < nr2; ++ir2) { // sample r^2
          double const r2 = ir2*hr2, r = std::sqrt(r2);
          r2c[ir2] = std::exp(-r2); // function evaluation here
          if (echo > 4) std::printf("%g %g\n", r, r2c[ir2]); // plot function value vs radius r
          rad_integral += r2c[ir2] * r;
      } // ir2
      rad_integral *= 2*constants::pi/inv_hr2;

      if (echo > 2) std::printf("\n# add_function()\n\n");
      double added{0};
      std::vector<double> values(g.all(), 0.0);
      add_function(values.data(), g, r2c, nr2, inv_hr2, &added, cnt);
      if (echo > 6) std::printf("\n# non-zero values on the Cartesian grid (sum = %g)\n", added);
      double xyz_integral{0};
      for (        int iz = 0; iz < g('z'); ++iz) {  double const vz = iz*g.h[2] - cnt[2];
          for (    int iy = 0; iy < g('y'); ++iy) {  double const vy = iy*g.h[1] - cnt[1];
              for (int ix = 0; ix < g('x'); ++ix) {  double const vx = ix*g.h[0] - cnt[0];
                  auto const ixyz = (iz*g('y') + iy)*g('x') + ix;
                  auto const val = values[ixyz];
                  if (0 != val) {
                      if (echo > 6) std::printf("%g %g\n", std::sqrt(vz*vz + vy*vy + vx*vx), val); // plot function value vs radius r
                      xyz_integral += val;
                  } // non-zero
              } // ix
          } // iy
      } // iz
      xyz_integral *= g.dV(); // volume element
      auto const diff = xyz_integral - rad_integral;
      if (echo > 1) std::printf("# grid integral = %g  radial integral = %g  difference = %.1e (%.3f %%)\n",
                                  xyz_integral, rad_integral, diff, 100*diff/rad_integral);
      return std::abs(diff/rad_integral) > 4e-4;
  } // test_add_function

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_create_and_destroy(echo);
      stat += test_add_function(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace real_space
