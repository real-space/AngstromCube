#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, std::snprintf
#include <cstdint> // uint32_t
#include <vector> // std::vector<T>

#include "status.hxx" // status_t
#include "geometry_analysis.hxx" // ::read_xyz_file
#include "data_view.hxx" // view2D<>
#include "display_units.h" // Ang, _Ang
#include "inline_math.hxx" // set
#include "simple_math.hxx" // ::invert3x3

namespace symmetry_group {

  double constexpr s30 = 0.866025403784438597;

  double const _rotation_matrices[32*3*3] = {
    1.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  1.0,
   -1.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  1.0,
   -1.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0, -1.0,
    1.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0, -1.0,
    0.0,  1.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0, -1.0,
    0.0, -1.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,
    0.0, -1.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  1.0,
    0.0,  1.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0,
    0.0,  0.0,  1.0,  0.0, -1.0,  0.0,  1.0,  0.0,  0.0,
    0.0,  0.0, -1.0,  0.0, -1.0,  0.0, -1.0,  0.0,  0.0,
    0.0,  0.0, -1.0,  0.0,  1.0,  0.0,  1.0,  0.0,  0.0,
    0.0,  0.0,  1.0,  0.0,  1.0,  0.0, -1.0,  0.0,  0.0,
   -1.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  1.0,  0.0,
   -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0, -1.0,  0.0,
    1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  1.0,  0.0,
    1.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0, -1.0,  0.0,
    0.0,  0.0,  1.0,  1.0,  0.0,  0.0,  0.0,  1.0,  0.0,
    0.0,  0.0, -1.0, -1.0,  0.0,  0.0,  0.0,  1.0,  0.0,
    0.0,  0.0, -1.0,  1.0,  0.0,  0.0,  0.0, -1.0,  0.0,
    0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  0.0, -1.0,  0.0,
    0.0,  1.0,  0.0,  0.0,  0.0,  1.0,  1.0,  0.0,  0.0,
    0.0, -1.0,  0.0,  0.0,  0.0, -1.0,  1.0,  0.0,  0.0,
    0.0, -1.0,  0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0,
    0.0,  1.0,  0.0,  0.0,  0.0, -1.0, -1.0,  0.0,  0.0,
    0.5,  s30,  0.0, -s30,  0.5,  0.0,  0.0,  0.0,  1.0,
    0.5, -s30,  0.0,  s30,  0.5,  0.0,  0.0,  0.0,  1.0,
   -0.5,  s30,  0.0, -s30, -0.5,  0.0,  0.0,  0.0,  1.0,
   -0.5, -s30,  0.0,  s30, -0.5,  0.0,  0.0,  0.0,  1.0,
    0.5, -s30,  0.0, -s30, -0.5,  0.0,  0.0,  0.0, -1.0,
    0.5,  s30,  0.0,  s30, -0.5,  0.0,  0.0,  0.0, -1.0,
   -0.5, -s30,  0.0, -s30,  0.5,  0.0,  0.0,  0.0, -1.0,
   -0.5,  s30,  0.0,  s30,  0.5,  0.0,  0.0,  0.0, -1.0};


  inline double _check_set_s(double const overlap[3][3], double const rot[3][3], int64_t s[3][4]=nullptr) {
    // The inverse of the overlap matrix is applied
    int constexpr echo = 0;
    if (echo > 7) std::printf("# %s  ", __func__);
    double max_eps{0};
    for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) {
            auto const value = overlap[0][j]*rot[k][0]
                             + overlap[1][j]*rot[k][1]
                             + overlap[2][j]*rot[k][2];
            if (echo > 7) std::printf(" %g", value);
            double const eps = std::abs(std::round(value) - value);
            // if (eps > 1e-6) {
            //     // If a noninteger is obtained, this implies that this operation
            //     // is not a symmetry operation for the given lattice
            //     if (echo > 7) std::printf(" @(%d,%d) --> no\n", j, k);
            //     return 1 + j*3 + k;
            // } // is_integer
            max_eps = std::max(max_eps, eps);
            if (nullptr != s) s[j][k] = std::round(value);
        } // k
    } // j
    if (echo > 7) std::printf("  --> yes\n");
    return max_eps;
  } // _check_set_s


  inline status_t find_symm_bravais_latt(double const cell[3][4], int const echo=0) {

      double rot[3][3];
      for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
              rot[i][j] = cell[i][0]*cell[j][0]
                        + cell[i][1]*cell[j][1]
                        + cell[i][2]*cell[j][2];
          } // j
          if (echo > 8) std::printf("# lattice vector overlap:  %15.9f %15.9f %15.9f Bohr^2\n", rot[i][0], rot[i][1], rot[i][2]);
      } // i

      double overlap[3][3];
      auto const det_rot = simple_math::invert3x3(overlap[0], 3, rot[0], 3);
      if (echo > 0) std::printf("# %s: determinant is %g Bohr^6\n", __func__, det_rot);
      for (int d = 0; d < 3*(echo > 7); ++d) {
          std::printf("# overlap inverse:  %15.9f %15.9f %15.9f Bohr^{-2}\n", overlap[d][0], overlap[d][1], overlap[d][2]);
      } // d

      auto const rotation_matrices = (double const (*)[3][3])_rotation_matrices;

      int64_t s[48][3][4]; // integer representations

      char yes[33]; yes[32] = '\0'; // string summary which ops are included

      double max_eps_accepted{0},
             max_eps_rejected{0};
      int nrot{0};
      for (int irot = 0; irot < 32; ++irot) {
          auto const rm = rotation_matrices[irot];
          // for each possible symmetry
          for (int j = 0; j < 3; ++j) {
              double rat[3]; // Cartesian coordinates of the rotated vector
              for (int i = 0; i < 3; ++i) {
                  rat[i] = cell[j][0]*rm[0][i] + cell[j][1]*rm[1][i] + cell[j][2]*rm[2][i];
              } // i
              if (echo > 19) std::printf("# irot=%d, j=%d rat:  %15.9f %15.9f %15.9f\n", irot, j, rat[0], rat[1], rat[2]);

              // The rotated vector is projected onto the direct lattice
              for (int k = 0; k < 3; ++k) {
                  rot[j][k] = cell[k][0]*rat[0] + cell[k][1]*rat[1] + cell[k][2]*rat[2];
              } // k
              if (echo > 17) std::printf("# irot=%d, j=%d rot:  %15.9f %15.9f %15.9f\n", irot, j, rot[j][0], rot[j][1], rot[j][2]);
          } // j

          auto const eps = _check_set_s(overlap, rot, s[irot]);
          if (eps < 1e-6) {
              if (echo > 13) std::printf("# symmetry op#%d is included, epsilon= %.1e\n", irot, eps);
              ++nrot;
              yes[irot] = '0' + ((irot + 1) % 10); // mark each included symmetry operation
              max_eps_accepted = std::max(max_eps_accepted, eps);
          } else {
              if (echo > 15) std::printf("# symmetry op#%d is not included, epsilon= %g\n", irot, eps);
              yes[irot] = ' ';
              max_eps_rejected = std::max(max_eps_rejected, eps);
          }

      } // irot
      if (echo > 3) std::printf("# %d symmetry operations \'%s\'\n", nrot, yes);
      if (echo > 4) std::printf("# %d symmetry operations included (max.eps= %.1e), %d not included (max.eps= %.1g)\n",
                                    nrot, max_eps_accepted, 32 - nrot, max_eps_rejected);

      switch (nrot) {
          case 1: case 2: case 4: case 6: case 8: case 12: case 24: break; // ok
          default:
            warn("Bravais lattice has wrong number of symmetries: ", nrot);
            if (echo > 1) std::printf("# cannot have %d symmetries in a Bravais lattice, symmetries are disabled\n", nrot);
            nrot = 1;
      } // switch nrot

      // Set the inversion symmetry (Bravais lattices always have inversion symmetry)
      for (int irot = 0; irot < nrot; ++irot) {
          set(s[irot + nrot][0], 12, s[irot][0], -1ll);
      } // irot
      nrot *= 2;
      if (echo > 2) std::printf("# %d symmetry operations included\n", nrot);

      return 0;
  } // find_symm_bravais_latt


#ifdef NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_symmetry_group(int const echo=0) {
      status_t stat(0);
      double constexpr a = 1.25; // in Bohr
      double const cubic[3][3][4] = {{{ a,0,0,0},{0, a,0,0},{0,0, a,0}}   // simple cubic
                                    ,{{-a,a,a,0},{a,-a,a,0},{a,a,-a,0}}   // bcc
                                    ,{{ 0,a,a,0},{a, 0,a,0},{a,a, 0,0}}}; // fcc
      double cell_in_file[3][4];
      char const test_name[4][32] = {"simple cubic", "body-centered cubic", "face-centered cubic", "atoms.xyz"};
      for (int t = 0; t < 4; ++t) {
          if (echo > 0) std::printf("\n# %s: test %s\n", __func__, test_name[t]);
          auto lattice_vectors = cubic[t];
          if (t >= 3) {
              int n_atoms;
              view2D<double> xyzZ;
              stat += geometry_analysis::read_xyz_file(xyzZ, n_atoms, cell_in_file, nullptr, "atoms.xyz", echo);
              lattice_vectors = cell_in_file;
          }
          for (int d = 0; d < 3*(echo > 4); ++d) {
              std::printf("# lattice vector #%d:  %15.9f %15.9f %15.9f %s\n", d+1,
                lattice_vectors[d][0]*Ang, lattice_vectors[d][1]*Ang, lattice_vectors[d][2]*Ang, _Ang);
          } // d
          stat += find_symm_bravais_latt(lattice_vectors, echo);
      } // t
      return stat;
  } // test_symmetry_group

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_symmetry_group(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace symmetry_group
