#pragma once

#ifndef NO_UNIT_TESTS
  #include <cstdio> // if (echo > 3) printf
  #include <cstdint> // int64_t, std::sif (echo > 3) printf, uint8_t
  #include <string> // std::string
  #include <vector> // std::vector<T>

#endif // NO_UNIT_TESTS

#include "constants.hxx" // ::pi
#include "display_units.h" // Ang, _Ang

typedef int status_t;

namespace shift_boundary {
  // A shift boundary condition is helpful for Cartesian methods:
  // If we force the (periodic) unit cell to be orthorhombic
  // in order to be able to use factorizable finite-difference stencils
  // for the kinetic energy operator we suffer from the drawback
  // that we can only compute simple cubic crystal structure with a
  // single atom in the unit cell.
  // In order to test the code, calculations with a moderate number
  // of degrees of freedom in both, the density and the wave functions
  // or equivalent Green function are important.
  // A single atom per unit cell in a highly symmetric
  // setup might require more k-points in the Brillouin zone but
  // fixes the symmetries which should give a fast SCF convergence.
  // With shift boundary conditions, the unit cell is still rectangular
  // but not periodic in a traditional sense:
  // We allow non-zero (an we restrict ourselfes to positive) entries
  // in the upper triangular Bravais matrix:
  //        xx  xy  xz
  //        0   yy  yz
  //        0   0   zz
  // For simplicity, we look at the 2D case now
  //        xx  xy
  //        0   yy
  // Then, a unit cell is connected to its periodic images as
  //
  //                  |
  //         --+------+-----------------+--
  //           |                        |
  //           |                        |
  //           |                        |
  //           |                        |
  //     <-xy->|                        |    --> x-direction
  //  --+------+-----------------+------+--
  //  ^ |                        |
  //  | |                        |
  // yy |                        |
  //  | |                        |             ^
  //  v |                        |             | y-direction
  //  --+-----------------+------+--           |
  //     <----------xx---------->
  //
  // so the x-direction is periodic in the usual sense,
  // however, when we translate by one unit cell in the y-direction,
  // we have to accound for a shift of magnitude xy in x-direction too.
  //
  // This allows a Cartesian code to represent more crystal lattices,
  // however, there is still a restriction to the generality:
  // unless we accept an interpolation operator implied with the boundary
  // operations, we need to match the lattices, i.e. xy and xz need to be an
  // integer multiple of hx, the grid spacing in x-direction.
  // Similarly, zy needs to be in integer of hy.
  //
  // However, the most important high symmetry cases are FCC, HCP and BCC.
  // For FCC and BCC, it is sufficient that the number of grid points is even:
  //    - FCC: the cubic unit cell (a,a,a) with 4 atoms is reduced to 1 atom
  //      in a cell (a,a/2,a/2) with yz=a/2 and xz=a/2 or xy=a/2 (but not both)
  //    - BCC: the cubic unit cell (a,a,a) with 2 atoms is reduced to 1 atom
  //      in a cell (a,a,a/2) with shifts xz=a/2 and yz=a/2 and xy=0
  //    - HCP: The HEXagonal 2D basis of HCP can be brought into [1, sqrt(3/4)]
  //      shape with an xy-shift of 1/2 (in units of the nn-distance).
  //      Also here, the unit cell will stay 2-atomic. If we permute the indices,
  //      the could also be the xz

  // Note: we could start by implementing xz and yz!

  // There is an option to treat FCC lattices in an orthorhombic 2-atomic cell:
  // (a*sqrt(1/2), a*sqrt(1/2), a), however, this artifically breaks isotropy as
  // the x,y-directions cannot have the same grid spacing as the z-direction.
  // Furthermore, to reduce it to a 1-atomic cell, the rules for BCC apply.

  // Comment: a bit more systematically would be to allow xy, yz, zx
  // to be non-zero as these form a cyclic sequence. However, this
  // opens an undescribed gap of volume xy*yz*zx, so we see that it
  // must be an upper triangular matrix and entries below the diagonal
  // must vanish. It becomes clear that xy, yz, zx are all on the same
  // side-diagonal of the matrix and therefore, their product enters the
  // determinant and, hence, changes the volume.

  // How to get there systematically?

  // Example for FCC, Bravais matrix in units of a/2
  //      1  1  0
  //      0  1  1
  //      1  0  1
  // now step 1) rotate (here, we do Gauss elemination by row-wise operations)
  //      1  1  0
  //      0  1  1
  //      0 -1  1
  // now step 2) and identify:
  //      1  1  0              zz  yz  xz           yy  yz  xy
  //      0  1  1               0  yy  xy    or      0  zz  xz
  //      0  0  2               0   0  xx            0   0  xx

  // Example for BCC, Bravais matrix in units of a/2
  //      1  1 -1
  //     -1  1  1
  //      1 -1  1
  // now step 1) rotate (here, we do Gauss elemination by row-wise operations)
  //      1  1 -1
  //      0  2  0
  //      0 -2  2
  // now step 2) and identify:
  //      1  1 -1              zz  yz  xz
  //      0  2  0               0  yy  xy
  //      0  0  2               0   0  xx

  // Alternative: (more complicated)
  // Step 1) Rotate a given Bravais matrix such that one cell vector is (xx,0,0).
  // Step 2) Rotate around that vector until the second cell vector is (xy,yy,0).
  // Done. We found (xz,yz,zz).

  // For reasons of the performance, having always real periodicity in x-direction
  // is good if the spatial indices ixyz = (iz*Ny + iy)*Nx + ix
  // e.g. in the finite-difference Laplacian of the Poisson solver.


  // For density and potentials, we can work in all three cases FCC, BCC and SC
  // with (a/2,a/2,a/2) with mirror boundaries everywhere if the atoms are positioned
  // on the cell borders.

  // ToDo: how to treat k-points?

    inline status_t test_plane_wave(int const echo=9, int const structure=4) {
      status_t stat{0};
      char const structure_name[][4] = {"sc ","bcc","hcp","fcc"};
      double amat[3][4], bmat[3][4];
      for(int ij = 0; ij < 3*4; ++ij) { amat[0][ij] = 0; bmat[0][ij] = 0; }
      double const alat = 4.1741; // e.g. Gold in hcp or fcc
      double const ahalf = 0.5 * alat;
      if (echo > 3) printf("\n# structure = %s  lattice constant = %g %s\n", structure_name[structure - 1], alat*Ang, _Ang);
      if (4 == structure) { // fcc
          amat[0][0] = 2*ahalf; amat[0][1] = ahalf;  amat[0][2] = 0;
          amat[1][0] = 0;       amat[1][1] = ahalf;  amat[1][2] = ahalf;
          amat[2][0] = 0;       amat[2][1] = 0;      amat[2][2] = ahalf;
      } else
      if (3 == structure) { // hex in xy-direction, c/a for hcp
          double const s34 = std::sqrt(.75), s83=std::sqrt(8/3.);
          double const ann = ahalf*std::sqrt(2.); // nearest neighbor bond length as in fcc
          amat[0][0] = ann; amat[0][1] = ann*0.5;
          amat[1][0] = 0;   amat[1][1] = ann*s34;
                                                      amat[2][2] = ann*s83;
      } else
      if (2 == structure) { // bcc
          amat[0][0] = 2*ahalf; amat[0][1] = ahalf;    amat[0][2] = 0;
          amat[1][0] = 0;       amat[1][1] = 2*ahalf;  amat[1][2] = ahalf;
          amat[2][0] = 0;       amat[2][1] = 0;        amat[2][2] = ahalf;
      } else
      if (1 == structure) { // sc
          for(int d = 0; d < 3; ++d) amat[d][d] = 2*ahalf;
      } // fcc

      // invert amat to find bmat
      double const detinv = 1./(amat[0][0]*amat[1][1]*amat[2][2]);
      for(int i = 0; i < 3; ++i) {     int const i1 = (i + 1)%3, i2 = (i + 2)%3;
          for(int j = 0; j < 3; ++j) { int const j1 = (j + 1)%3, j2 = (j + 2)%3;
              bmat[j][i] = ( amat[i1][j1] * amat[i2][j2]
                           - amat[i1][j2] * amat[i2][j1] )*detinv;
          } // j
      } // i
      
      if (echo > 4) {
          // show both matrices
          for(int i = 0; i < 3; ++i) {
              printf("#  bmat %c %8.3f%8.3f%8.3f   amat %8.3f%8.3f%8.3f\n", i+'x',
                  bmat[i][0],bmat[i][1],bmat[i][2],   amat[i][0],amat[i][1],amat[i][2]);
          } // i
      } // echo
      
      // check if product a*b and b*a are 3x3 unit matrices
      double maxdev[] = {0, 0};
      for(int i = 0; i < 3; ++i) {
        if (echo > 6) printf("# i=%i ", i);
        for(int j = 0; j < 3; ++j) {
          double uij{0}, uji{0};
          for(int k = 0; k < 3; ++k) {
            uij += amat[i][k] * bmat[k][j];
            uji += bmat[i][k] * amat[k][j];
          } // k - contraction index
          maxdev[0] = std::max(maxdev[0], uij);
          maxdev[1] = std::max(maxdev[1], uji);
          if (echo > 6) printf("%8.3f%8.3f ", uji, uij);
          if (echo > 8) printf("%.1e %.1e ", uji - (i == j), uij - (i == j));
        } // j
        if (echo > 6) printf("\n");
      } // i
      if (echo > 3) printf("# %s after inversion largest deviation is %.1e (a*b) and %.1e (b*a)\n", 
                               __func__, maxdev[0], maxdev[1]);

      // test: set up periodic+shifted BC and diagonalize the free electron Hamiltonian
      // check that

      //                  |
      //         --+------+-----------------+--
      //           |                        |
      //           |                        |
      //           |    phi=e^{i*Ly*kyy}    |     phi=e^{i*(Lx*kxx+Ly*kyy)}
      //           |                        |
      //           |                        |    --> x-direction
      //  --+------+-----------------+------+--
      //    |                        |
      //    |                        |
      //    |        phi=e^0         |   phi=e^{i*Lx*kxx}
      //    |                        |             ^
      //    |                        |             | y-direction
      //  --+-----------------+------+----------   |
      //                      |

      // With Cartesian real space lattices, it is even simple
      // to implement non-symmorphic symmetries like glide reflections
      // as needed in HCP:
// ### HCP: Layer structure ABABAB (while FCC has a layer structure ABCABC)
// #
// #  A--- ---B---A
// #  |           |
// #  | B   A   _ | a
// #  |           |
// #  A--- ---B---A
// #    a sqrt 3
// #        
   // cell    1 sqrt3 c/a
   // # from initially 4 atoms at the fractional coordinates
   // Cd      0   0   0     A
   // Cd      1:2 3:6 0     A
   // Cd      1:2 1:6 1:2   B
   // Cd      0   4:6 1:2   B

   // # we come to 2 atoms with an xy=0.5 shift boundary
// #  A--- -+ B   A
// #  |     |     |
// #  | B   A-----+
// #  |     |     |
// #  A--- -+ B   A  
// #    a sqrt 3

   // cell    1 sqrt.75 c/a        
   // Cd      0   0   0     A
   // Cd      1:2 1:6 1:2   B

   // and with another symmetry enabler, we can reduce to one atom:
   // cell    1 sqrt.75 c/a/2        
   // if we have a glide reflection boundary (shift mirror) by 1/6:

  //
  //    |                        |
  //  --+-------+----------------+-------+--
  //            | 25  24  23  22  21  20 |
  //            B 15  14  13  12  11  10 |
  //            | 05  04  03  02  01  00 |       --> y-direction
  //  --+-------+-------()-------+-------+--
  //    | 00  01  02  03  04  05 |      
  //    | 10  11  12  13  14  15 A
  //    | 20  21  22  23  24  25 |             ^
  //  --+-------+----------------+-------+--   | z-direction
  //            |                        |
  // so we reflect at the point marked by () in order to map A onto B    


   // What do we need to simulate diamond structure with 1 atoms?
        // diamond = FCC + basis, has inversion symmetry
        // the FCC approach with Cartesian cells is
        // cell (a, a/2, a/2)
        // which can, with a mirror plane be reduced to (a/2, a/2, a/2), [also for the wave functions??]
        // this reduction cannot be done for diamond, as 
        // in one of the cubes, there is the second base atom body centered.
        // however, if we reflect by a plane through atom #2 and rotate by pi around it,
        // so the cell may be shaped (a, a/2, a/4)

    // Probably, it is simpler to work on the lattice periodic part of
    // a Bloch wave with a k-dependent operator instead of having only
    // the boundary k-dependent and the Hamiltonian inside the bulk real.

      return stat;
    } // test_plane_wave

#ifdef  NO_UNIT_TESTS
    inline status_t all_tests(int const echo=0) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

    inline status_t all_tests(int const echo=3) {
      if (echo > 1) printf("\n# %s: %s\n\n", __FILE__, __func__);
      status_t stat{0};
      for(int structure = 1; structure <= 4; ++structure) {
          stat += test_plane_wave(echo, structure);
      } // structure
      return stat;
    } // all_tests
#endif // NO_UNIT_TESTS

} // namespace shift_boundary
