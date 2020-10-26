#pragma once

#include <cstdio> // printf
#include <algorithm> // std::max

#include "status.hxx" // status_t
// #include "complex_tools.hxx" // complex_name, is_complex, conjugate, to_complex_t
#include "data_view.hxx" // view2D, view4D
#include "inline_math.hxx" // set

namespace brillouin_zone {

  inline
  int get_kpoint_mesh(
        view2D<double> & mesh
      , unsigned const nv[3]
      , int const echo=0 // log-level
  ) {
      unsigned const n[] = {std::max(1u, nv[0]), std::max(1u, nv[1]), std::max(1u, nv[2])};
      size_t const nfull = n[2]*n[1]*n[0];
      view4D<double> full(n[2], n[1], n[0], 4); // get temporary memory
      int const ishift[] = {int(n[0] - 1), int(n[1] - 1), int(n[2] - 1)};
      double const denom[] = {.5/n[0], .5/n[1], .5/n[2]};
      double const wfull = 1./nfull;
      double xyzw[] = {0, 0, 0, wfull};
      for(int iz = 0; iz < n[2]; ++iz) {  xyzw[2] = (2*iz - ishift[2])*denom[2];
      for(int iy = 0; iy < n[1]; ++iy) {  xyzw[1] = (2*iy - ishift[1])*denom[1];
      for(int ix = 0; ix < n[0]; ++ix) {  xyzw[0] = (2*ix - ishift[0])*denom[0];
          set(full(iz, iy, ix), 4, xyzw);
          if (echo > 8) printf("# kpoint mesh entry %9.6f %9.6f %9.6f weight= %g\n", xyzw[0],xyzw[1],xyzw[2], xyzw[3]);
      }}} // iz iy iz
      int const nmesh = nfull; // unreduced
      mesh = view2D<double>(nmesh, 4); // get memory
      set(mesh.data(), nmesh*4, full.data()); // copy
      if (echo > 3) printf("# kpoint mesh with %d x %d x %d has %d points\n", nv[0],nv[1],nv[2], nmesh);
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
      return nm - n*n*n;
  } // test_mesh

  inline status_t all_tests(int const echo=0) {
      status_t status(0);
      status += test_mesh(echo);
      return status;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace brillouin_zone
