#pragma once

#include <cstdio> // printf
#include <cmath> // std::ceil, std::sqrt
#include <cstdint> // int8_t

#include "inline_math.hxx"
#include "recorded_warnings.hxx" // warn
#include "data_view.hxx" // view2D<T>
#include "status.hxx" // status_t

  int constexpr Periodic_Boundary =  1;
  int constexpr Isolated_Boundary =  0;
  int constexpr Mirrored_Boundary = -1;
  int constexpr Invalid_Boundary  = -9;
  
namespace boundary_condition {

  inline int periodic_images(view2D<double> & ipos // pointer to array of periodic positions (n,4)
                           , double const cell[3] // orthorhombic cell parameters
                           , int const bc[3] // boundary condition selectors
                           , float const rcut // truncation radius
                           , int const echo=0 // log-level
                           , view2D<int8_t> *iidx=nullptr // pointer to array of indices (n,4)
                            ) { // log-level
      auto const cell_diagonal2 = pow2(rcut)
                                + pow2(cell[0]) + pow2(cell[1]) + pow2(cell[2]);
      int ni_xyz[3], ni_max{1};
      if (rcut < 0) warn("A negative cutoff radius leads to only one image! rcut = %g a.u.", rcut);
      for(int d = 0; d < 3; ++d) {
          if (Periodic_Boundary == bc[d]) {
              ni_xyz[d] = std::max(0, int(std::ceil(rcut/cell[d])));
              assert( ni_xyz[d] <= 127 ); // warning: int8_t has range [-128, 127]
              ni_max *= (ni_xyz[d] + 1 + ni_xyz[d]);
          } else {
              // other boundary conditions, e.g. isolated
              ni_xyz[d] = 0;
          } // periodic
      } // d
      if (echo > 5) printf("# %s: check %d x %d x %d = %d images max.\n",
          __func__, 1+2*ni_xyz[0], 1+2*ni_xyz[1], 1+2*ni_xyz[2], ni_max);
      view2D<double> pos(ni_max, 4, 0.0); // get memory
      view2D<int8_t> idx(ni_max, 4, 0); // get memory
      int ni = 1; // at least one periodic images is always there: (0,0,0)
      for         (int iz = -ni_xyz[2]; iz <= ni_xyz[2]; ++iz) {  auto const pz = iz*cell[2];
          for     (int iy = -ni_xyz[1]; iy <= ni_xyz[1]; ++iy) {  auto const py = iy*cell[1];
              for (int ix = -ni_xyz[0]; ix <= ni_xyz[0]; ++ix) {  auto const px = ix*cell[0];
                  auto const d2 = pow2(px) + pow2(py) + pow2(pz);
#ifdef DEVEL
                  char mark{' '};
#endif                  
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
#endif
                      } // d2 > 0
                  } // d2 < cell_diagonal2
#ifdef DEVEL
                  if (echo > 6) {
                      if (ix == -ni_xyz[0]) {
                          if (iy == -ni_xyz[1]) printf("# %s z=%i\n", __func__, iz);
                          printf("#%4i  | ", iy); // before first x
                      } // first x
                      printf("%c", mark);
                      if (ix == ni_xyz[0]) printf(" |\n"); // after last x
                  } // echo
#endif
              } // ix
          } // iy
      } // iz
      if (echo > 1) printf("# %s: found %d of %d images\n", __func__, ni, ni_max);
      
      // export array of periodic positions
      ipos = view2D<double>(ni, 4); // get memory
      set(ipos.data(), ni*4, pos.data()); // copy

      if (iidx) { 
          *iidx = view2D<int8_t>(ni, 4); // get memory
          set(iidx->data(), ni*4, idx.data()); // copy
      } // export indices
      
      return ni;
  } // periodic_images
  
  
  inline int periodic_images_old(double **ipos // pointer to array of periodic positions [4*n]
                           , double const cell[3] // orthorhombic cell parameters
                           , int const bc[3] // boundary condition selectors
                           , float const rcut // truncation radius
                           , int const echo=0 // log-level
                           , int8_t **iidx=nullptr // pointer to array of indices [4*n]
                            ) { // log-level
      auto const cell_diagonal2 = pow2(rcut)
                                + pow2(cell[0]) + pow2(cell[1]) + pow2(cell[2]);
      int ni_xyz[3], ni_max{1};
      if (rcut < 0) warn("A negative cutoff radius leads to only one image! rcut = %g a.u.", rcut);
      for(int d = 0; d < 3; ++d) {
          if (Periodic_Boundary == bc[d]) {
              ni_xyz[d] = std::max(0, int(std::ceil(rcut/cell[d])));
              assert(ni_xyz[d] <= 127); // warning: int8_t has range [-128, 127]
              ni_max *= (ni_xyz[d] + 1 + ni_xyz[d]);
          } else {
              // other boundary conditions, e.g. isolated
              ni_xyz[d] = 0;
          } // periodic
      } // d
      if (echo > 5) printf("# %s: check %d x %d x %d = %d images max.\n",
          __func__, 1+2*ni_xyz[0], 1+2*ni_xyz[1], 1+2*ni_xyz[2], ni_max);
      auto const pos = new double[ni_max*4];
      auto const idx = new int8_t[ni_max*4];
      for(int d = 0; d < 4; ++d) {
          pos[0*4 + d] = 0; // zero-th image/original position
      } // d
      int ni = 1;
      for         (int iz = -ni_xyz[2]; iz <= ni_xyz[2]; ++iz) {
          for     (int iy = -ni_xyz[1]; iy <= ni_xyz[1]; ++iy) {
              for (int ix = -ni_xyz[0]; ix <= ni_xyz[0]; ++ix) {
                  auto const px = ix*cell[0];
                  auto const py = iy*cell[1];
                  auto const pz = iz*cell[2];
                  auto const d2 = pow2(px) + pow2(py) + pow2(pz);
#ifdef DEVEL
                  char mark{' '};
#endif                  
                  if (d2 < cell_diagonal2) {
                      if (d2 > 0) { // exclude the origin (because we want the original image as index #0)
                          pos[ni*4 + 0] = px;
                          pos[ni*4 + 1] = py;
                          pos[ni*4 + 2] = pz;
                          pos[ni*4 + 3] = d2; // distance^2 - no use
                          // the indices are important for the Bloch-phases
                          idx[ni*4 + 0] = ix;
                          idx[ni*4 + 1] = iy;
                          idx[ni*4 + 2] = iz;
                          idx[ni*4 + 3] =  0; // no use
                          ++ni; // count the number of images inside
#ifdef DEVEL
                          mark = 'o';
                      } else {
                          mark = 'x';
#endif
                      } // d2 > 0
                  } // d2 < cell_diagonal2
#ifdef DEVEL
                  if (echo > 6) {
                      if (ix == -ni_xyz[0]) {
                          if (iy == -ni_xyz[1]) printf("# %s z=%i\n", __func__, iz);
                          printf("#%4i  | ", iy); // before first x
                      } // first x
                      printf("%c", mark);
                      if (ix == ni_xyz[0]) printf(" |\n"); // after last x
                  } // echo
#endif
              } // ix
          } // iy
      } // iz
      *ipos = pos; // export array of periodic positions
      if (nullptr != iidx) { *iidx = idx; } else { delete[] idx; } // export or release
      if (echo > 1) printf("# %s: found %d of %d images\n", __func__, ni, ni_max);
      return ni;
  } // periodic_images
  
  inline int fromString(char const *string, int const echo=0, char const dir='?') {
      char const first = *string;
      switch (first | 32) { // ignore case with | 32
        case 'p': case '1':
            if (echo > 0) printf("# interpret \"%s\" as periodic boundary condition in %c-direction\n", string, dir);
            return Periodic_Boundary;
        case 'i': case '0':
            if (echo > 0) printf("# interpret \"%s\" as isolated boundary condition in %c-direction\n", string, dir);
            return Isolated_Boundary;
        case 'm': case '-':
            if (echo > 0) printf("# interpret \"%s\" as mirror boundary condition in %c-direction\n", string, dir);
            return Mirrored_Boundary;
        default :
            if (echo > 0) printf("# cannot interpret \"%s\" as boundary condition in %c-direction\n", string, dir);
            return Invalid_Boundary;
      } // switch
  } // fromString
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  inline status_t test_periodic_images(int const echo=7) {
      if (echo > 2) printf("\n# %s %s \n", __FILE__, __func__);
      double const cell[] = {1,2,3}, rcut = 6.f;
      int const bc[] = {Periodic_Boundary, Periodic_Boundary, Isolated_Boundary};
      view2D<double> ipos;
      view2D<int8_t> iidx;
      auto const nai = periodic_images(ipos, cell, bc, rcut, echo, &iidx);
      if (echo > 2) printf("# found %d images\n", nai);
      return 0;
  } // test_periodic_images

  inline status_t all_tests(int const echo=7) {
      if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
      status_t status(0);
      status += test_periodic_images(echo);
      return status;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace boundary_condition
