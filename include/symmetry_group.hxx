#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, std::snprintf
#include <cstdint> // uint32_t, int32_t
#include <vector> // std::vector<T>

#include "status.hxx" // status_t
#include "geometry_input.hxx" // ::read_xyz_file
#include "data_view.hxx" // view2D<>, view3D<>
#include "display_units.h" // Ang, _Ang
#include "inline_math.hxx" // set
#include "print_tools.hxx" // printf_vector(format, vec, n, final, scale, add)
#include "simple_math.hxx" // ::invert3x3, ::determinant
#include "recorded_warnings.hxx" // warn, error

namespace symmetry_group {

  template <typename T>
  inline T determinant(T const mat[3][3]) {
      return simple_math::determinant(3, mat[0], 3);
  } // determinant


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
          default: {
              warn("Bravais lattice has wrong number of symmetries: ", nrot);
              if (echo > 1) std::printf("# cannot have %d symmetries in a Bravais lattice, symmetries are disabled\n", nrot);
              nrot = 1;
          }
      } // switch nrot

      // Set the inversion symmetry (Bravais lattices always have inversion symmetry)
      for (int irot = 0; irot < nrot; ++irot) {
          set(s[irot + nrot][0], 12, s[irot][0], -1ll);
      } // irot
      nrot *= 2;
      if (echo > 2) std::printf("# %d symmetry operations included\n", nrot);

      return 0;
  } // find_symm_bravais_latt

  template <typename T>
  void get_string(char str[3], T const mat[3][3]) {
      for (int i = 0; i < 3; ++i) {
          char c{'?'};
          for (int j = 0; j < 3; ++j) {
              auto const rm = mat[i][j];
              if (0 != rm) {
                //   assert('?' == c);
                  c = j + ((rm > 0) ? 'x' : 'X');
              }
          } // j
        //   assert('?' != c);
          str[i] = c;
      } // i
  } // get_string

  template <typename real_t, int Stride=3>
  status_t generate_cubic_symmetry_matrix(real_t mat[][Stride], int const i48) {
      if (i48 <   0) return -1;
      if (i48 >= 48) return  1;
      int8_t constexpr permutation[6][3] = {{0,1,2},{1,2,0},{2,0,1},{2,1,0},{0,2,1},{1,0,2}};
      int const i8 = i48 & 0x7; // in [0,7]
      int const i6 = i48 >> 3;  // in [0,5]
      assert(i6*8 + i8 == i48);
      for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
              mat[i][j] = 0; // clear
          } // j
          int const sign = (i8 >> i) & 0x1;
          int const j3 = permutation[i6][i];
          mat[i][j3] = real_t(1 - 2*sign);
      } // i
      return 0; // success
  } // generate_cubic_symmetry_matrix


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
              stat += geometry_input::read_xyz_file(xyzZ, n_atoms, cell_in_file, nullptr, "atoms.xyz", echo);
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

  inline status_t test_generate_group(int const echo=0) {
      status_t stat(0);
      // generate all operations of 3x3 matrices with only {-1,0,1} entries:
      int n_det_1{0};
      for (int i19683 = 0; i19683 < 19683; ++i19683) {
          int mat[3][3];
          {
              int j{i19683};
              for (int ij = 0; ij < 3*3; ++ij) {
                  mat[0][ij] = (j % 3) - 1;
                  j /= 3;
              } // ij
          }
          auto const det = determinant(mat);
          if (1 == det*det) {
              char str[4] = "???";
              get_string(str, mat);
              auto const *const m = mat[0]; // flat array int[9] for display
              if (echo > 11) std::printf("# i19683=%6i string=%s mat= %d %d %d  %d %d %d  %d %d %d\n",
                                            i19683, str, m[0],m[1],m[2], m[3],m[4],m[5], m[6],m[7],m[8]);
              ++n_det_1;
          } // |det| == 1
      } // i19683
      if (echo > 3) std::printf("# %s: %d of 19683 have |det| == 1\n", __func__, n_det_1);
      return stat;
  } // test_generate_group


  typedef int32_t hast18_t; // 18 bits needed

  template <typename T>
  inline hast18_t mat_hash18(T const mat[3][3]) {
      hast18_t hash{0};
      for (int ij = 0; ij < 3*3; ++ij) {
          auto const rm = int(mat[0][ij]);
          assert(std::abs(rm) < 2); // allowed matrix entries are {-1, 0, 1}
          hash <<= 2; // *= 4
          hash |= (rm & 0x3); // use bits 00(0.0), 01(1.0) or 11(-1.0), never 10
      } // ij
      return hash;
  } // mat_hash18

  template <typename T>
  inline int hash2mat(T mat[3][3], hast18_t const hash) {
      if (hash < 0)       return -1;
      if (hash >= 262144) return  1;
      int8_t constexpr translate[4] = {0, 1, 99, -1};
      auto h{hash};
      for (int ij = 8; ij >= 0; --ij) {
          int const i3 = h & 0x3;
          mat[0][ij] = translate[i3 & 0x3];
          assert(99 != mat[0][ij] && "invalid hash?");
          h >>= 2; // divide by 4
      } // ij
      return 0; // success
  } // hash2mat

  // a more compact hash would be int16_t which would require that we encode into 19683 = 3^9
  typedef int16_t hash15_t; // 15 bits needed

  template <typename T>
  inline hash15_t mat_hash15(T const mat[3][3]) {
      hash15_t hash{0}; // 15bit needed
      for (int ij = 0; ij < 3*3; ++ij) {
          auto const rm = int(mat[0][ij]);
          assert(std::abs(rm) < 2); // allowed matrix entries are {-1, 0, 1}
          hash = hash*3 + rm + 3*(rm < 0); // translated tokens   { 2, 0, 1}
      } // ij
      return hash;
  } // mat_hash15

  template <typename T>
  inline int hash2mat(T mat[3][3], hash15_t const hash) {
      if (hash < 0)      return -1;
      if (hash >= 19683) return  1;
      int8_t constexpr translate[4] = {0, 1, -1, 99};
      auto h{hash};
      for (int ij = 8; ij >= 0; --ij) {
          int const i3 = h % 3;
          mat[0][ij] = translate[i3 & 0x3];
          assert(99 != mat[0][ij] && "invalid hash?");
          h /= 3; // divide by 3
      } // ij
      return 0; // success
  } // hash2mat


  template <typename T, int N=3>
  inline void matmul(T axb[N][N], T const a[N][N], T const b[N][N]) {
      for (int i = 0; i < N; ++i) {
          for (int j = 0; j < N; ++j) {
              T t(0);
              for (int k = 0; k < N; ++k) {
                  t += a[i][k] * b[k][j];
              } // k
              axb[i][j] = t;
          } // j
      } // i
  } // matmul

  inline status_t test_check_group24(int const echo=0) {
      status_t stat(0);
      int mat[24][3][3];
      hash15_t hash15s[24];
      hast18_t hash18s[24];
      for (int i24 = 0; i24 < 24; ++i24) {
          for (int ij = 0; ij < 3*3; ++ij) {
              int const rm = std::round(_rotation_matrices[i24*3*3 + ij]);
              assert(std::abs(rm) <= 1); // rm may only be in {-1, 0, 1}
              mat[i24][0][ij] = rm;
          } // ij
          hash15s[i24] = mat_hash15(mat[i24]);
          hash18s[i24] = mat_hash18(mat[i24]);
      } // i24
      // multiplication table
      int not_found{0};
      int multab[24][24];
      for (int i24 = 0; i24 < 24; ++i24) {
          for (int j24 = 0; j24 < 24; ++j24) {
              int prod[3][3];
              matmul(prod, mat[i24], mat[j24]);
              auto const hash = mat_hash15(prod);
              // now search hash in hashes
              int k24{-1}; // -1:not_found
              for (int k = 0; k < 24; ++k) {
                  if (hash == hash15s[k]) { k24 = k; };
              } // k
              multab[i24][j24] = k24;
              not_found += (-1 == k24); // count how many times it has not been found in hashes
          } // j24
          if (echo > 5) { std::printf("# multiplication[%2d]:  ", i24); printf_vector(" %d", multab[i24], 24); }
      } // i24
      if (not_found > 0) {
          warn("%d symmetry matrix products are not part of the group", not_found);
          ++stat;
      } // not_found
      return stat;
  } // test_check_group24

  inline status_t test_check_group48(int const echo=0) {
      status_t stat(0);
      int mat[48][3][3]; // generate all 48 cubic symmetry matrices
      hash15_t hash15s[48]; // and their hash codes
      hast18_t hash18s[48]; // and their hash codes
      for (int i48 = 0; i48 < 48; ++i48) {
          generate_cubic_symmetry_matrix(mat[i48], i48);
          hash15s[i48] = mat_hash15(mat[i48]);
          hash18s[i48] = mat_hash18(mat[i48]);
          { // scope: check the hash2mat functions
              int m[3][3];
              for (int ibits = 15; ibits <= 18; ibits += 3) {
                auto const works = (15 == ibits) ? hash2mat(m, hash15s[i48])
                                                 : hash2mat(m, hash18s[i48]);
                assert(0 == works);
                for (int ij = 0; ij < 3*3; ++ij) {
                    if (m[0][ij] != mat[i48][0][ij]) error("i48=%i ij=%i m=%i /= %i =mat", i48, ij, m[0][ij], mat[i48][0][ij]);
                } // ij
              } // ibits in {15, 18}
          } // scope
      } // i48
      int8_t multab[48][48]; // multiplication table
      int not_found{0};
      for (int i48 = 0; i48 < 48; ++i48) {
          for (int j48 = 0; j48 < 48; ++j48) {
              // make sure that symmetries are unique
              assert((hash15s[i48] == hash15s[j48]) == (i48 == j48));
              // check product
              int prod[3][3];
              matmul(prod, mat[i48], mat[j48]);
              auto const prod_mat_hash15 = mat_hash15(prod);
              auto const prod_mat_hash18 = mat_hash18(prod);
              // now search hash in hashes
              int k48{-1}; // -1:not_found
              for (int k = 0; k < 48; ++k) {
                  if (prod_mat_hash15 == hash15s[k]) {
                      k48 = k; // found
                      assert(prod_mat_hash18 == hash18s[k]);
                  }
              } // k
              multab[i48][j48] = k48;
              not_found += (-1 == k48); // count how many times it has not been found in hashes
          } // j48
          if (echo > 5) { std::printf("# multiplication[%2d]:  ", i48); printf_vector(" %d", multab[i48], 48); }
      } // i48
      if (not_found > 0) {
          warn("%d symmetry matrix products are not part of the group", not_found);
          ++stat;
      } // not_found
      return stat;
  } // test_check_group48

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_symmetry_group(echo);
      stat += test_generate_group(echo);
      stat += test_check_group24(echo);
      stat += test_check_group48(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace symmetry_group
