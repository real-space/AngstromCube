#pragma once

#include <cstdio> // std::printf
#include <cmath> // std::ceil, std::sqrt
#include <cstdint> // int8_t

#include "inline_math.hxx"
#include "data_view.hxx" // view2D<T>
#include "status.hxx" // status_t
#ifndef STANDALONE_TEST
  #include "recorded_warnings.hxx" // warn
#else
  #define warn std::printf
#endif

  int8_t constexpr Periodic_Boundary =  1;
  int8_t constexpr Isolated_Boundary =  0;
  int8_t constexpr Mirrored_Boundary = -1;
  int8_t constexpr Invalid_Boundary  = -2;

namespace boundary_condition {

  inline int periodic_images( // returns the number of images found
        view2D<double> & ipos // array of periodic positions (n,4)
      , double const cell[3]  // orthorhombic cell parameters
      , int8_t const bc[3]       // boundary condition selectors
      , float  const rcut      // truncation radius
      , int    const echo=0      // log-level
      , view2D<int8_t> *iidx=nullptr // optional: pointer to array of indices (n,4)
  ) {
      double const cell_diagonal2 = pow2(rcut)
                                  + pow2(cell[0]) + pow2(cell[1]) + pow2(cell[2]);
      int ni_xyz[3], ni_max{1};
      if (rcut < 0) warn("A negative cutoff radius leads to only one image! rcut = %g a.u.", rcut);
      for (int d = 0; d < 3; ++d) {
          if (Periodic_Boundary == bc[d]) {
              ni_xyz[d] = std::max(0, int(std::ceil(rcut/cell[d])));
              assert( ni_xyz[d] <= 127 ); // warning: int8_t has range [-128, 127]
              ni_max *= (ni_xyz[d]*2 + 1);
          } else {
              // other boundary conditions, e.g. isolated
              ni_xyz[d] = 0;
          } // periodic
      } // d
      if (echo > 5) std::printf("# %s: check %d x %d x %d = %d images max.\n",
              __func__, 1+2*ni_xyz[0], 1+2*ni_xyz[1], 1+2*ni_xyz[2], ni_max);
      view2D<double> pos(ni_max, 4, 0.0); // get memory
      view2D<int8_t> idx(ni_max, 4, 0); // get memory
      int ni{1}; // at least one periodic images is always there: (0,0,0)
      for         (int iz = -ni_xyz[2]; iz <= ni_xyz[2]; ++iz) {  auto const pz = iz*cell[2];
          for     (int iy = -ni_xyz[1]; iy <= ni_xyz[1]; ++iy) {  auto const py = iy*cell[1];
              for (int ix = -ni_xyz[0]; ix <= ni_xyz[0]; ++ix) {  auto const px = ix*cell[0];
                  auto const d2 = pow2(px) + pow2(py) + pow2(pz);
#ifdef DEVEL
                  char mark{' '};
#endif // DEVEL
                  if (d2 < cell_diagonal2) {
                      if (d2 > 0) { // exclude the origin (that is already index #0)
                          pos(ni,0) = px;
                          pos(ni,1) = py;
                          pos(ni,2) = pz;
                          pos(ni,3) = d2; // distance^2 - no use
                          // the indices are important for the Bloch-phases
                          idx(ni,0) = ix;
                          idx(ni,1) = iy;
                          idx(ni,2) = iz;
                          idx(ni,3) =  0; // no use
                          ++ni; // count the number of images inside
#ifdef DEVEL
                          mark = 'o';
                      } else {
                          mark = 'x';
#endif // DEVEL
                      } // d2 > 0
                  } // d2 < cell_diagonal2
#ifdef DEVEL
                  if (echo > 6) {
                      if (ix == -ni_xyz[0]) {
                          if (iy == -ni_xyz[1]) std::printf("# %s z=%i\n", __func__, iz);
                          std::printf("#%4i  | ", iy); // before first x
                      } // first x
                      std::printf("%c", mark);
                      if (ix == ni_xyz[0]) std::printf(" |\n"); // after last x
                  } // echo
#endif // DEVEL
              } // ix
          } // iy
      } // iz
      if (echo > 1) std::printf("# %s: found %d of %d images\n", __func__, ni, ni_max);

      // export array of periodic positions
      ipos = view2D<double>(ni, 4); // get memory
      set(ipos.data(), ni*4, pos.data()); // copy

      if (iidx) { 
          *iidx = view2D<int8_t>(ni, 4); // get memory
          set(iidx->data(), ni*4, idx.data()); // copy
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
              case 'm': case '-': bc = Mirrored_Boundary; break;
          } // switch
      } // nullptr != string
      if (echo > 0) {
          char const bc_names[][12] = {"isolated", "periodic", "invalid", "mirror"};
          std::printf("# interpret \"%s\" as %s boundary condition in %c-direction\n", 
                      string, bc_names[bc & 0x3], dir);
      } // echo
      return bc;
  } // fromString

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_periodic_images(int const echo=0) {
      if (echo > 2) std::printf("\n# %s %s \n", __FILE__, __func__);
      double const cell[] = {1,2,3};
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
