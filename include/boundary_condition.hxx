#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cmath> // std::ceil, ::sqrt, ::abs
#include <cstdint> // int8_t

#include "inline_math.hxx"
#include "data_view.hxx" // view2D<T>
#include "status.hxx" // status_t
#ifndef   STANDALONE_TEST
  #include "recorded_warnings.hxx" // warn
#else  // STANDALONE_TEST
  #define warn std::printf
#endif // STANDALONE_TEST

  int8_t constexpr Periodic_Boundary =  1;
  int8_t constexpr Isolated_Boundary =  0;
  int8_t constexpr Mirrored_Boundary = -1;
  int8_t constexpr Invalid_Boundary  = -2;

  int8_t constexpr Vacuum_Boundary = 2;
  // The vacuum boundary condition is an addition to Isolated_Boundary and Periodic_Boundary from boundary_condition.hxx
  // Vacuum_Boundary means that the Green function extends from its source coordinate up to the truncation radius
  // even if this exceeds the isolated boundary at which potential values stop to be defined.
  // The potential is continued as zero beyond the isolated boundary but the tail of the Green function is allowed to fade there.
  // k-points are not relevant.

  // We could think of a useful new boundary condition
  int8_t constexpr Repeat_Boundary = 3;
  // The repeat boundary allows to compute with small unit cells but with realistic truncation radii.
  // The local potential is repeated periodically (like Periodic_Boundary).
  // Sources are only relevant inside the small unit cell but again, the Green function tail exceeds the cell boundaries
  // and extends up to the truncation radius.
  // For the dyadic potential, we may not reduce over periodic atom images but create copies of the atoms.
  // k-points are not relevant.

  // The wrap boundary condition is an addition to Periodic_Boundary from boundary_condition.hxx
  int8_t constexpr Wrap_Boundary = 5;
  // Wrap_Boundary means that the truncation sphere fits into the cell, so k-points have no effect.
  // Nevertheless, it cannot be treated like Isolated_Boundary since target block coordinates may need to be wrapped.
  // However, it could be viewed as Repeat_Boundary...


namespace boundary_condition {

  inline int periodic_images( // returns ni, the number of images found
        view2D<double> & ipos   // array of periodic positions (ni,4)
      , double const cell[3][4] // lower triangular cell matrix
      , int8_t const bc[3]      // boundary condition selectors
      , float  const rcut       // truncation radius
      , int    const echo=0     // log-level
      , view2D<int8_t> *iidx=nullptr // optional: pointer to array of indices (ni,4)
  ) {
      double const cell_diagonal2 = pow2(rcut)
                                  + pow2(cell[0][0]) + pow2(cell[1][1]) + pow2(cell[2][2]); // ToDo: needs adjustment?
      int ni_xyz[3], ni_max{1};
      if (rcut < 0) warn("A negative cutoff radius leads to only one image! rcut = %g a.u.", rcut);
      for (int d = 0; d < 3; ++d) {
          if (Periodic_Boundary == bc[d]) {
              ni_xyz[d] = std::max(0, int(std::ceil(rcut/std::abs(cell[d][d]))));
              assert( ni_xyz[d] <= 127 ); // warning: int8_t has range [-128, 127]
              ni_max *= (ni_xyz[d]*2 + 1);
          } else {
              // other boundary conditions, e.g. isolated
              ni_xyz[d] = 0;
          } // periodic
      } // d
      int const nx = ni_xyz[0], ny = ni_xyz[1], nz = ni_xyz[2];
      if (echo > 5) std::printf("# %s: check %d x %d x %d = %d images max.\n",
              __func__, nx+1+nx, ny+1+ny, nz+1+nz, ni_max);

#ifndef   GENERAL_CELL
      assert((0 == cell[0][1]) && (0 == cell[1][2]) && (0 == cell[0][2]) && "the cell is not a lower triangular matrix");
#endif // GENERAL_CELL

#ifdef    DEVEL
      view3D<char> mark(2*ny + 1 + 1, 2*nz + 1, 2*nx + 1 + 3, ' '); // data layout y,z,x for plotting
#endif // DEVEL
      view2D<double> pos(ni_max, 4, 0.0); // get memory
      view2D<int8_t> idx(ni_max, 4, 0);   // get memory
      int ni{1}; // at least one periodic images is always there: (0,0,0)
      for         (int iz = -nz; iz <= nz; ++iz) {
#ifndef   GENERAL_CELL
          double const pz[3]  = {iz*cell[2][0], iz*cell[2][1], iz*cell[2][2]}; // can deal with a lower triangular cell matrix
#endif // GENERAL_CELL
          for     (int iy = -ny; iy <= ny; ++iy) {
#ifndef   GENERAL_CELL
              double const pyz[3] = {iy*cell[1][0] + pz[0], iy*cell[1][1] + pz[1], pz[2]}; // can deal with a lower triangular cell matrix
#endif // GENERAL_CELL
              for (int ix = -nx; ix <= nx; ++ix) {
#ifdef    GENERAL_CELL
                  double p[3];
                  for (int d{0}; d < 3; ++d) {
                      p[d] = ix*cell[0][d] + iy*cell[1][d] + iz*cell[2][d];
                  } // d
#else  // GENERAL_CELL
                  double const px = ix*cell[0][0]; // can deal with a lower triangular cell matrix
                  double const p[3] = {pyz[0] + px, pyz[1], pyz[2]};
#endif // GENERAL_CELL
                  double const d2 = pow2(p[0]) + pow2(p[1]) + pow2(p[2]);
                  if (d2 < cell_diagonal2) {
                      if (d2 > 0) { // exclude the origin (that is already index #0)
                          pos(ni,0) = p[0];
                          pos(ni,1) = p[1];
                          pos(ni,2) = p[2];
                          pos(ni,3) = d2; // distance^2 - no use
                          // the indices are important for the Bloch-phases
                          idx(ni,0) = ix;
                          idx(ni,1) = iy;
                          idx(ni,2) = iz;
                          idx(ni,3) =  0; // no use
                          ++ni; // count the number of images inside
#ifdef    DEVEL
                          mark(iy + ny,iz + nz,ix + nx) = 'o'; // data layout y,z,x
                      } else {
                          mark(iy + ny,iz + nz,ix + nx) = 'x'; // data layout y,z,x
#endif // DEVEL
                      } // d2 > 0
                  } // d2 < cell_diagonal2
              } // ix
          } // iy
      } // iz
      if (echo > 1) std::printf("# %s: found %d of %d images\n", __func__, ni, ni_max);
#ifdef    DEVEL
      if (echo > 6) { // display in terminal
          int const mz = std::max(1, std::min(2*nz, 96/(2*nx + 1 + 3))); // maximum number of z-slices so the line does not get too long
          std::printf("# %s x=%d...%d, z=%d...%d(max %d), total=%d:\n", __func__, -nx, nx, -nz, mz - nz, nz, ni);
          for     (int iz{0}; iz <= mz;   ++iz) {
              for (int ix{0}; ix <= 2*nx; ++ix) {
                  mark(2*ny + 1,iz,ix) = char('0' + (std::abs(ix - nx) % 10)); // make an x-legend
              } // ix
          } // iz
          mark(2*ny + 1,mz,2*nx + 3) = '\0'; // zero-termination so the legend does not get too long
          std::printf("#     x=   %s\n", mark(2*ny + 1,0)); // the x-legend
          for     (int iy{0}; iy <= 2*ny; ++iy) {
              for (int iz{0}; iz <= 2*nz; ++iz) {
                  mark(iy,iz,2*nx + 2) = '|';
              } // iz
              mark(iy,mz,2*nx + 3) = '\0'; // zero-termination so the line does not get too long
              std::printf("# y=%4i | %s\n", iy - ny, mark(iy,0)); // plot
          } // iy
      } // echo
#endif // DEVEL

      // export array of periodic positions
      ipos = view2D<double>(ni, 4); // get memory
      set(ipos.data(), ni*4, pos.data()); // deep copy

      if (nullptr != iidx) {
          *iidx = view2D<int8_t>(ni, 4); // get memory
          set(iidx->data(), ni*4, idx.data()); // deep copy
      } // export indices

      return ni;
  } // periodic_images

  inline int8_t fromString(
        char const *string
      , int const echo=0
      , char const dir='?'
  ) {
      int8_t bc{Invalid_Boundary};
      if (nullptr != string) {
          char const first = *string;
          switch (first | 32) { // ignore case with | 32
              case 'p': case '1': bc = Periodic_Boundary; break;
              case 'i': case '0': bc = Isolated_Boundary; break;
              case 'm': case '-': bc = Mirrored_Boundary; break; // experimental
          } // switch
      } // nullptr != string
      if (echo > 0) {
          char const bc_names[][12] = {"isolated", "periodic", "invalid", "mirror"};
          std::printf("# interpret \"%s\" as %s boundary condition in %c-direction\n",
                      string, bc_names[bc & 0x3], dir);
      } // echo
      return bc;
  } // fromString

#ifdef    NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

  inline status_t test_periodic_images(int const echo=0) {
      if (echo > 2) std::printf("\n# %s %s \n", __FILE__, __func__);
      double const cell[3][4] = {{1,0,0,0}, {0,2,0,0}, {0,0,3,0}};
      float  const rcut = 6;
      int8_t const bc[] = {Periodic_Boundary, Periodic_Boundary, Isolated_Boundary};
      view2D<double> ipos;
      view2D<int8_t> iidx;
      auto const nai = periodic_images(ipos, cell, bc, rcut, echo, &iidx);
      if (echo > 2) std::printf("# found %d periodic images\n", nai);
      auto const nai2 = periodic_images(ipos, cell, bc, rcut);
      return (nai2 - nai);
  } // test_periodic_images

  inline status_t test_fromString_single(char const bc_strings[][16], int const echo=0) {
      if (echo > 2) std::printf("\n# %s %s \n", __FILE__, __func__);
      status_t stat(0);
      for (int8_t bc = Invalid_Boundary; bc <= Periodic_Boundary; ++bc) {
          stat += (bc != fromString(bc_strings[bc & 0x3], echo));
      } // bc
      return stat;
  } // test_fromString_single

  inline status_t test_fromString(int const echo=0) {
      // test the parser with different strings
      if (echo > 2) std::printf("\n# %s %s \n", __FILE__, __func__);
      status_t stat(0);
      {   char const bc_strings[][16] = {"isolated", "periodic", "?invalid", "mirror"}; // {0, 1, -2, -1}
          stat += test_fromString_single(bc_strings, echo);   }
      {   char const bc_strings[][16] = {"i", "p", "_", "m"}; // {0, 1, -2, -1}
          stat += test_fromString_single(bc_strings, echo);   }
      {   char const bc_strings[][16] = {"I", "P", "#", "M"}; // {0, 1, -2, -1}
          stat += test_fromString_single(bc_strings, echo);   }
      {   char const bc_strings[][16] = {"0", "1", "*", "-"}; // {0, 1, -2, -1}
          stat += test_fromString_single(bc_strings, echo);   }
      return stat;
  } // test_fromString

  inline status_t all_tests(int const echo=0) {
      if (echo > 0) std::printf("\n# %s %s\n", __FILE__, __func__);
      status_t stat(0);
      stat += test_periodic_images(echo);
      stat += test_fromString(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace boundary_condition
