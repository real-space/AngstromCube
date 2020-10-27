#pragma once

#include <cstdio> // printf
#include <algorithm> // std::max
#include <cmath> // std::round

#include "status.hxx" // status_t
// #include "complex_tools.hxx" // complex_name, is_complex, conjugate, to_complex_t
#include "data_view.hxx" // view2D, view4D
#include "inline_math.hxx" // set, pow2, product
#include "print_tools.hxx" // printf_vector

namespace brillouin_zone {

  inline
  int get_kpoint_mesh(
        view2D<double> & mesh
      , unsigned const nv[3]
      , int const echo=0 // log-level
  ) {
      int constexpr WEIGHT=3, COUNT=0, WRITE=1;
      unsigned const n[] = {std::max(1u, nv[0]), std::max(1u, nv[1]), std::max(1u, nv[2])};
      size_t const nfull = n[2]*n[1]*n[0];
      view4D<double> full(n[2], n[1], n[0], 4); // get temporary memory
      int const ishift[] = {int(n[0] - 1), int(n[1] - 1), int(n[2] - 1)}; // could also be double
      double const denom[] = {.5/n[0], .5/n[1], .5/n[2]};
      double xyzw[] = {0, 0, 0, 1};
      for(int iz = 0; iz < n[2]; ++iz) {  xyzw[2] = (2*iz - ishift[2])*denom[2];
      for(int iy = 0; iy < n[1]; ++iy) {  xyzw[1] = (2*iy - ishift[1])*denom[1];
      for(int ix = 0; ix < n[0]; ++ix) {  xyzw[0] = (2*ix - ishift[0])*denom[0];
          set(full(iz, iy, ix), 4, xyzw);
          if (echo > 18) printf("# kpoint mesh entry %9.6f %9.6f %9.6f weight= %g\n", xyzw[0],xyzw[1],xyzw[2], xyzw[3]);
      }}} // iz iy iz

//       int const nsymmetries = 1;
//       view3D<int8_t> symmetry(nsymmetries, 4, 4, 0); // use only the 3x3 sub matrices of the 4x4 memory chunks
//       for(int d = 0; d < 3; ++d) symmetry(0,d,d) = -1; // time reversal

      for(int iz = 0; iz < n[2]; ++iz) {
      for(int iy = 0; iy < n[1]; ++iy) {
      for(int ix = 0; ix < n[0]; ++ix) {
          // apply time reversal symmetry (k and -k produce the same density)
          double xyz[3];
          set(xyz, 3, full(iz,iy,ix), -1.0); // go from k --> -k, ToDo: implement symmetry operation
          // convert to nearest integers
          int jxyz[3];
          for(int d = 0; d < 3; ++d) {
              jxyz[d] = int(std::round(n[d]*xyz[d] + 0.5*ishift[d]));
          } // d
          int const jx = jxyz[0], jy = jxyz[1], jz = jxyz[2];
          if (echo > 16) {
              printf("# symmetry #%i maps kpoint", 0);
              printf_vector(" %9.6f", full(iz,iy,ix), 3, " to");
              printf_vector(" %9.6f", full(jz,jy,jx), 3);
          } // echo

          auto const w8 = full(jz,jy,jx,WEIGHT);
          if (w8 > 0) {
              double diff2{0};
              for(int d = 0; d < 3; ++d) diff2 += pow2(xyz[d] - full(jz,jy,jx,d));
              if (diff2 < 1e-12) {
                  // k + (-k) == nullvector, transfer weight to (iz,iy,ix)
                  full(iz,iy,ix,WEIGHT) += w8;
                  full(jz,jy,jx,WEIGHT) -= w8;
                  if (echo > 14) {
                      printf("# transfer kpoint weight %g from", w8);
                      printf_vector(" %9.6f", full(jz,jy,jx), 3, " to");
                      printf_vector(" %9.6f", full(iz,iy,ix), 3);
                  } // echo
              } // d2 is zero
          } // w8 > 0
      }}} // iz iy iz

      // copy kpoints with positive weight into result list
      // up to here, all weights were integer, so we have to divide by the denominator nfull
      double const wfull = 1./nfull;
      double const weighted[] = {1, 1, 1, wfull};
      int nmesh{0};
      for(int i01 = COUNT; i01 <= WRITE; i01 += (WRITE - COUNT)) {
          double w8sum{0};
          int imesh{0};
          for(int iz = 0; iz < n[2]; ++iz) {
          for(int iy = 0; iy < n[1]; ++iy) {
          for(int ix = 0; ix < n[0]; ++ix) {
              auto const w8 = full(iz,iy,ix,WEIGHT);
              if (w8 > 0) {
                  if (WRITE == i01) product(mesh[imesh], 4, full(iz,iy,ix), weighted); // copy and apply weight denominator
                  w8sum += w8;
                  ++imesh;
              } // w8 > 0
          }}} // iz iy iz
          if (WRITE == i01) {
              assert(nmesh == imesh && "1st and 2nd time counting did not agree");
          } else {
              nmesh = imesh;
              if (echo > 8) printf("# %d kpoints have positive weight, sum = 1 + %.1e\n", nmesh, w8sum*wfull - 1);
              mesh = view2D<double>(nmesh, 4); // get memory
          } // COUNT or WRITE
      } // twice

      if (echo > 3) printf("# kpoint mesh with %d x %d x %d has %d points\n", n[0],n[1],n[2], nmesh);
      return nmesh;
  } // get_kpoint_mesh

  inline
  int get_kpoint_mesh(
        view2D<double> & mesh
      , unsigned const n
      , int const echo=0 // log-level
  ) {
      unsigned const nv[] = {n, n, n};
      return get_kpoint_mesh(mesh, nv, echo);
  } // get_kpoint_mesh


#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_mesh(int const echo=0, unsigned const n=8) {
      view2D<double> mesh;
      auto const nm = get_kpoint_mesh(mesh, n, echo);
      return (nm < 1);
  } // test_mesh

  inline status_t all_tests(int const echo=0) {
      status_t status(0);
      for(int n = 0; n < 9; ++n) {
          status += test_mesh(echo, n);
      } // n
      return status;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace brillouin_zone
