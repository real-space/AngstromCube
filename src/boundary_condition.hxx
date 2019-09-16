#pragma once

#include <cstdio> // printf
#include <cmath> // std::ceil
#include "inline_math.hxx"

typedef int status_t;

  int constexpr Periodic_Boundary =  1;
  int constexpr Isolated_Boundary =  0;
  int constexpr Mirrored_Boundary = -1;
  int constexpr Invalid_Boundary = -9;

namespace boundary_condition {

  inline int periodic_images(double **ipos, double const cell[3], // orthorhombic cell parameters
                             int const bc[3], double const rcut, int const echo=0) {
      auto const cell_diagonal2 = pow2(cell[0]) + pow2(cell[1]) + pow2(cell[2]) + pow2(rcut);
      int ni_xyz[3], ni_max = 1;
      for(int d = 0; d < 3; ++d) {
          ni_xyz[d] = 0;
          if (Periodic_Boundary == bc[d]) {
              ni_xyz[d] = (int)std::ceil(rcut/cell[d]);
          } // periodic
          ni_max *= (ni_xyz[d] + 1 + ni_xyz[d]);
      } // d
      if (echo > 5) printf("# check %d x %d x %d = %d images\n", 1+2*ni_xyz[0], 1+2*ni_xyz[1], 1+2*ni_xyz[2], ni_max);
      auto const pos = new double[ni_max*4];
      for(int d = 0; d < 4; ++d) pos[d] = 0; // zero image
      int ni = 1;
      for         (int iz = -ni_xyz[2]; iz <= ni_xyz[2]; ++iz) {
          for     (int iy = -ni_xyz[1]; iy <= ni_xyz[1]; ++iy) {
              if (echo > 6) printf("#   ");
              for (int ix = -ni_xyz[0]; ix <= ni_xyz[0]; ++ix) {
                  auto const px = ix*cell[0];
                  auto const py = iy*cell[1];
                  auto const pz = iz*cell[2];
                  auto const d2 = pow2(px) + pow2(py) + pow2(pz);
                  char mark;
                  if (d2 > 0) { // exclude the origin (because we want the original image as index #0)
                      if (d2 < cell_diagonal2) {
                          pos[ni*4 + 0] = px;
                          pos[ni*4 + 1] = py;
                          pos[ni*4 + 2] = pz;
                          pos[ni*4 + 3] = 0;
                          ++ni;
                              mark = 'x';
                      } else  mark = ' ';
                  } else      mark = 'o';
                  if (echo > 6) printf("%c", mark);
              } // ix
              if (echo > 6) printf("\n");
          } // iy
      } // iz
      *ipos = pos;
      if (echo > 1) printf("# %s: found %d images\n", __func__, ni);
      return ni;
  } // periodic_images
  
  inline int fromString(char const *string, int const echo=0) {
      char const first = *string;
      switch (first | 32) { // ignore case with | 32
        case 'p': case '1':
            if (echo > 0) printf("# interpret \"%s\" as periodic boundary condition\n", string);
            return Periodic_Boundary; break;
        case 'i': case '0': 
            if (echo > 0) printf("# interpret \"%s\" as isolated boundary condition\n", string);
            return Isolated_Boundary; break;
        default :
            if (echo > 0) printf("# cannot interpret \"%s\" as boundary condition\n", string);
            return Invalid_Boundary;
      } // switch
  } // fromString
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  inline status_t test_periodic_images(int const echo=7) {
      if (echo > 2) printf("\n# %s %s \n", __FILE__, __func__);
      double const cell[] = {1,2,3}, rcut = 6;
      int const bc[] = {Periodic_Boundary, Periodic_Boundary, Isolated_Boundary};
      double *ipos = nullptr;
      int const nai = periodic_images(&ipos, cell, bc, rcut, echo);
      delete[] ipos;
      if (echo > 2) printf("# found %d images\n", nai);
      return 0;
  } // test_periodic_images

  inline status_t all_tests(int const echo=7) {
    if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
    auto status = 0;
    status += test_periodic_images(echo);
    return status;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace boundary_condition
