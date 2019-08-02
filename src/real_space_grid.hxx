#pragma once

#include <cstdint> // uint32_t
#include <cstdio> // printf
#include "inline_math.hxx" // set, scale
#include "bessel_transform.hxx" // Bessel_j0

typedef int status_t;

namespace real_space_grid {

  template<typename real_t, int D0> // D0: inner dimension, vector length
  class grid_t {
  public:
      real_t* values;
  private:
      uint32_t dims[4]; // 0,1,2:real-space grid dimensions, 3:outer dim
  public:
      double h[3], inv_h[3]; // grid spacings and their inverse

      grid_t(int const dim[3], int const dim_outer=1, real_t* const v=nullptr)
       : h{1,1,1}, inv_h{1,1,1} {
          dims[0] = std::max(1, dim[0]); // x
          dims[1] = std::max(1, dim[1]); // y
          dims[2] = std::max(1, dim[2]); // z
          dims[3] = std::max(1, dim_outer);
          long const nnumbers = dims[3] * dims[2] * dims[1] * dims[0] * D0;
          if (nnumbers > 0) {
              if (nullptr == v) {
                  printf("# allocate a grid with %d * %d x %d x %d * %d real_%ld = %.6f GByte\n", 
                      dims[3], dims[2], dims[1], dims[0], D0, sizeof(real_t), nnumbers*1e-9*sizeof(real_t));
                  values = new real_t[nnumbers];
              } else {
                  printf("# use %p as grid with %d * %d x %d x %d * %d real_%ld\n", 
                      (void*)v, dims[3], dims[2], dims[1], dims[0], D0, sizeof(real_t));
                  values = v;
              }
          } else {
              printf("# grid invalid: <D0=%d> dims={%d, %d, %d,  %d}\n", 
                      D0, dims[0], dims[1], dims[2], dims[3]);
              values = nullptr;
          }
      } // constructor

      ~grid_t() {
          if (nullptr != values) {
              long const nnumbers = dims[3] * dims[2] * dims[1] * dims[0] * D0;
              printf("# release a grid with %d * %d x %d x %d * %d real_%ld = %.6f GByte, values=%p\n",
                  dims[3], dims[2], dims[1], dims[0], D0, sizeof(real_t), nnumbers*1e-9*sizeof(real_t), (void*)values);
//               delete[] values; // leads to a segfault
          } // free memory
      } // destructor
      
      status_t set_grid_spacing(float const hx, float const hy=-1, float const hz=-1) {
          status_t stat = 0;
          float const h3[] = {hx, (hy<0)?hx:hy, (hz<0)?hx:hy};
          for(int i3 = 0; i3 < 3; ++i3) {
              h[i3] = h3[i3]; // convert to double
              if ((0.0 != h[i3]) && (h[i3] == h[i3])) {
                  inv_h[i3] = 1./h[i3]; // invert only here
              } else ++stat;
          } // i3
          return stat;
      } // set
      
      inline int dim(char const xyz) const { return ('w' == (xyz | 32)) ? D0 : dims[(xyz | 32) - 120]; }
      inline int dim(int const x0y1z2) const { return dims[x0y1z2]; }
      inline double dV(bool const Cartesian=true) const { return h[0]*h[1]*h[2]; } // volume element, assuming a Cartesian grid

      inline size_t all() const { return dims[3] * dims[2] * dims[1] * dims[0] * D0; }
      
  }; // class grid_t
  
  
  
  template<typename real_t, int D0> // inner dimension
  status_t add_function(grid_t<real_t,D0> &g, real_t added[], // grid values are modified
                        double const r2coeff[], int const ncoeff, float const hcoeff,
                        double const center[3]=nullptr, float const rcut=-1, double const factor=1) {
  // Add a spherically symmetric regular function to the grid.
  // The function is tabulated as r2coeff[0 <= hcoeff*r^2 < ncoeff][D0]
      assert(D0 == g.dim('w'));
      double c[3] = {0,0,0}; if (center) set(c, 3, center);
      bool const cutoff = (rcut >= 0); // negative values mean that we do not need a cutoff
      double const r2cut = cutoff ? rcut*rcut : (ncoeff - 1)/hcoeff;
      int imn[3], imx[3];
      size_t nwindow = 1;
      for(int i3 = 0; i3 < 3; ++i3) {
          int const M = g.dim(i3) - 1; // highest index
          if (cutoff) {
              imn[i3] = std::max(0, (int)std::floor((c[i3] - rcut)*g.inv_h[i3]));
              imx[i3] = std::min(M, (int)std::ceil ((c[i3] + rcut)*g.inv_h[i3]));
          } else {
              imn[i3] = 0;
              imx[i3] = M;
          } // cutoff
#ifdef  DEBUG
          printf("# %s window %c = %d elements from %d to %d\n", __func__, 'x'+i3, imx[i3] + 1 - imn[i3], imn[i3], imx[i3]);
#endif
          nwindow *= std::max(0, imx[i3] + 1 - imn[i3]);
      } // i3
      assert(hcoeff > 0);
      set(added, D0, (real_t)0); // clear
      size_t modified = 0, out_of_range = 0;
      for(            int iz = imn[2]; iz <= imx[2]; ++iz) {  double const vz = iz*g.h[2] - c[2], vz2 = vz*vz;
          for(        int iy = imn[1]; iy <= imx[1]; ++iy) {  double const vy = iy*g.h[1] - c[1], vy2 = vy*vy;
              if (vz2 + vy2 < r2cut) {
                  for(int ix = imn[0]; ix <= imx[0]; ++ix) {  double const vx = ix*g.h[0] - c[0], vx2 = vx*vx;
                      double const r2 = vz2 + vy2 + vx2;
                      if (r2 < r2cut) {
                          int const ixyz = (iz*g.dim('y') + iy)*g.dim('x') + ix;
                          int const ir2 = (int)(hcoeff*r2);
                          if (ir2 < ncoeff) {
                              double const w8 = hcoeff*r2 - ir2; // linear interpolation weight
                              int const ir2p1 = ir2 + 1;
                              for(int i0 = 0; i0 < D0; ++i0) { // vectorize
                                  auto const value_to_add = (r2coeff[ir2*D0 + i0]*(1 - w8)
                                     + ((ir2p1 < ncoeff) ? r2coeff[ir2p1*D0 + i0] : 0)*w8);
                                  g.values[ixyz*D0 + i0] += factor*value_to_add;
                                  added[i0]              += factor*value_to_add;
                              } // i0
                              ++modified;
                          } else ++out_of_range;
                      } // inside rcut
                  } // ix
              } // rcut for (y,z)
          } // iy
      } // iz
      scale(added, D0, g.dV()); // volume integral
#ifdef  DEBUG
      printf("# %s modified %.3f k inside a window of %.3f k on a grid of %.3f k grid values.\n", 
              __func__, modified*1e-3, nwindow*1e-3, g.dim('x')*g.dim('y')*g.dim('z')*1e-3); // show stats
#endif
      if (out_of_range)
          printf("# Warning! %s modified %ld points and found %ld entries out of range of the radial function!\n", 
              __func__, modified, out_of_range); // show stats
      return 0; // success
  } // add_function

  template<typename real_t, int D0> // inner dimension
  status_t bessel_projection(real_t q_coeff[], int const nq, float const dq,
                grid_t<real_t,D0> const &g, double const center[3]=nullptr, 
                float const rcut=-1, double const factor=1) {
      assert(D0 == g.dim('w'));
      double c[3] = {0,0,0}; if (center) set(c, 3, center);
      bool const cutoff = (rcut >= 0); // negative values mean that we do not need a cutoff
      double const r2cut = cutoff ? rcut*rcut : 100.; // stop at 10 Bohr
      int imn[3], imx[3];
      size_t nwindow = 1;
      for(int i3 = 0; i3 < 3; ++i3) {
          int const M = g.dim(i3) - 1; // highest index
          if (cutoff) {
              imn[i3] = std::max(0, (int)std::floor((c[i3] - rcut)*g.inv_h[i3]));
              imx[i3] = std::min(M, (int)std::ceil ((c[i3] + rcut)*g.inv_h[i3]));
          } else {
              imn[i3] = 0;
              imx[i3] = M;
          } // cutoff
#ifdef  DEBUG
          printf("# %s window %c = %d elements from %d to %d\n", __func__, 'x'+i3, imx[i3] + 1 - imn[i3], imn[i3], imx[i3]);
#endif
          nwindow *= std::max(0, imx[i3] + 1 - imn[i3]);
      } // i3
      set(q_coeff, nq*D0, (real_t)0); // clear
      for(            int iz = imn[2]; iz <= imx[2]; ++iz) {  double const vz = iz*g.h[2] - c[2], vz2 = vz*vz;
          for(        int iy = imn[1]; iy <= imx[1]; ++iy) {  double const vy = iy*g.h[1] - c[1], vy2 = vy*vy;
              if (vz2 + vy2 < r2cut) {
                  for(int ix = imn[0]; ix <= imx[0]; ++ix) {  double const vx = ix*g.h[0] - c[0], vx2 = vx*vx;
                      double const r2 = vz2 + vy2 + vx2;
                      if (r2 < r2cut) {
                          int const ixyz = (iz*g.dim('y') + iy)*g.dim('x') + ix;
                          double const r = std::sqrt(r2);
//                           printf("%g %g\n", r, g.values[ixyz*D0 + 0]); // DEBUG
                          for(int iq = 0; iq < nq; ++iq) {
                              double const q = iq*dq;
                              for(int i0 = 0; i0 < D0; ++i0) { // vectorize
                                  q_coeff[iq*D0 + i0] += g.values[ixyz*D0 + i0] * bessel_transform::Bessel_j0(q*r);
                              } // i0
                          } // iq
                      } // inside rcut
                  } // ix
//                   printf("\n"); // DEBUG
              } // rcut for (y,z)
          } // iy
      } // iz
      double const sqrt2pi = std::sqrt(2./constants::pi); // this makes the transform symmetric
      scale(q_coeff, nq*D0, (real_t)(g.dV()*factor*sqrt2pi)); // volume element, external factor, Bessel transform factor
      return 0; // success
  } // bessel_projection
  
  status_t all_tests();

} // namespace real_space_grid
