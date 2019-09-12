#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::floor
#include <cstdint> // uint8_t
#include <fstream>
#include <sstream>
#include <string>
#include <vector> // std::vector<T>

#include "geometry_analysis.hxx"

#include "boundary_condition.hxx" // Periodic_Boundary, Isolated_Boundary, periodic_images
#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set, pow2
#include "vector_math.hxx" // vec<n,T>
#include "chemical_symbol.h" // element_symbols

// #define FULL_DEBUG
// #define DEBUG

namespace geometry_analysis {
  
  double constexpr Bohr2Angstrom = 0.52917724924;
  double constexpr Angstrom2Bohr = 1./Bohr2Angstrom;
  
  status_t read_xyz_file(double **xyzz, int *n_atoms, char const *filename, int const echo=5) {
      
      std::ifstream infile(filename, std::ifstream::in);
      int natoms = 0, linenumber = 2;
      infile >> natoms; // read the number of atoms
      std::string line;
      std::getline(infile, line);
      std::getline(infile, line);
      if (echo > 2) printf("# expect %d atoms, comment line in file: %s\n", natoms, line.c_str());
      auto const xyzZ = (natoms > 0) ? new double[natoms*4] : nullptr;
      int na = 0;
      while (std::getline(infile, line)) {
          ++linenumber;
          std::istringstream iss(line);
          std::string Symbol;
          double px, py, pz; // positions
          if (!(iss >> Symbol >> px >> py >> pz)) { 
              printf("# failed parsing in %s:%d reads \"%s\", stop\n", filename, linenumber, line.c_str());
              break; // error
          }
          char const *Sy = Symbol.c_str();
          char const S = Sy[0], y = (Sy[1])? Sy[1] : ' ';
          if (echo > 7) printf("# %c%c  %16.9f %16.9f %16.9f\n", S,y, px, py, pz);
          int iZ = -1;
          { // scope: search matching iZ
              for(int Z = 0; Z < 128; ++Z) {
                  if ((S == element_symbols[2*Z + 0]) &&
                      (y == element_symbols[2*Z + 1])) iZ = Z;
              } // Z
          } // scope
          xyzZ[na*4 + 0] = px*Angstrom2Bohr;
          xyzZ[na*4 + 1] = py*Angstrom2Bohr;
          xyzZ[na*4 + 2] = pz*Angstrom2Bohr;
          xyzZ[na*4 + 3] = iZ;
          ++na;
          assert(na <= natoms);
          // process pair 
      } // parse file line by line
      
//    auto const cell[] = {63.872738414, 45.353423726, 45.353423726}; // DNA-cell in Bohr, bc={periodic,isolated,isolated}
      *xyzz = (natoms == na) ? xyzZ : nullptr;
      *n_atoms = na;
      return na - natoms; // returns 0 if the exact number of atoms has been found
  } // read_xyz_file
  
  

  float default_bond_length(int const Z1, int const Z2=0) {
    // data from http://chemwiki.ucdavis.edu/Theoretical_Chemistry/Chemical_Bonding/Bond_Order_and_Lengths
    uint8_t const bl[128] = { 0 // 0 // artificial
    ,31   //  H     1
    ,28   //  He    2
    ,128  //  Li    3
    ,96   //  Be    4
    ,84   //  B     5
    ,76   //  C     6
    ,73   //  N     7
    ,69   //  O     8
    ,71   //  F     9
    ,66   //  Ne    10
    ,57   //  Na    11
    ,141  //  Mg    12
    ,121  //  Al    13
    ,111  //  Si    14
    ,107  //  P     15
    ,105  //  S     16
    ,102  //  Cl    17
    ,106  //  Ar    18
    ,203  //  K     19
    ,176  //  Ca    20
    ,170  //  Sc    21
    ,160  //  Ti    22
    ,153  //  V     23
    ,139  //  Cr    24
    ,132  //  Mn    25
    ,110  //  Fe    26  110pm is good for bcc bulk
    ,139  //  Co    27
    ,124  //  Ni    28
    ,132  //  Cu    29
    ,122  //  Zn    30
    ,122  //  Ga    31
    ,120  //  Ge    32
    ,119  //  As    33
    ,120  //  Se    34
    ,120  //  Br    35
    ,116  //  Kr    36
    ,220  //  Rb    37
    ,195  //  Sr    38
    ,190  //  Y     39
    ,175  //  Zr    40
    ,164  //  Nb    41
    ,154  //  Mo    42
    ,147  //  Tc    43
    ,146  //  Ru    44
    ,142  //  Rh    45
    ,139  //  Pd    46
    ,145  //  Ag    47
    ,144  //  Cd    48
    ,142  //  In    49
    ,139  //  Sn    50
    ,139  //  Sb    51
    ,138  //  Te    52
    ,139  //  I     53
    ,140  //  Xe    54
    ,244  //  Cs    55
    ,215  //  Ba    56
    ,207  //  La    57
    ,204  //  Ce    58
    ,203  //  Pr    59
    ,201  //  Nd    60
    ,199  //  Pm    61
    ,198  //  Sm    62
    ,198  //  Eu    63
    ,196  //  Gd    64
    ,194  //  Tb    65
    ,192  //  Dy    66
    ,192  //  Ho    67
    ,189  //  Er    68
    ,190  //  Tm    69
    ,187  //  Yb    70
    ,187  //  Lu    71
    ,175  //  Hf    72
    ,170  //  Ta    73
    ,162  //  W     74
    ,151  //  Re    75
    ,144  //  Os    76
    ,141  //  Ir    77
    ,136  //  Pt    78
    ,136  //  Au    79
    ,132  //  Hg    80
    ,145  //  Tl    81
    ,146  //  Pb    82
    ,148  //  Bi    83
    ,140  //  Po    84
    ,150  //  At    85
    ,150  //  Rn    86
    ,255  //  Fr    87  // Fr reduced from 260
    ,221  //  Ra    88
    ,215  //  Ac    89
    ,206  //  Th    90
    ,200  //  Pa    91
    ,196  //  U     92
    ,190  //  Np    93
    ,187  //  Pu    94
    ,180  //  Am    95
    ,169  //  Cm    96
    ,166  //  Bk    97
    ,168  //  Cf    98
    ,165  //  Es    99
    ,167  //  Fm    100
    ,173  //  Md    101
    ,176  //  No    102
    ,161  //  Lr    103
    ,157  //  Rf    104
    ,149  //  Dr    105
    ,143  //  Sg    106
    ,141  //  Bh    107
    ,134  //  Hs    108
    ,129  //  Mt    109
    ,128  //  Ds    110
    ,121  //  Rg    111
    ,122  //  Cn    112
    ,136  //  ut    113
    ,143  //  uq    114
    ,162  //  up    115
    ,175  //  uh    116
    ,165  //  us    117
    ,157  //  uo    118
    ,159 ,160 ,161 ,162 ,163 ,164 ,165 ,166, 167}; // 119--127 invented numbers
    // data from http://chemwiki.ucdavis.edu/Theoretical_Chemistry/Chemical_Bonding/Bond_Order_and_Lengths
    float const picometer = .01889726;
    return (bl[Z1] + bl[Z2]) * picometer;
  } // default_bond_length
  
  
  status_t analysis(double const xyzZ[], int const natoms, 
                    double const cell[3], int const bc[3], int const echo=6) {
      if (echo > 1) printf("\n# %s:%s\n", __FILE__, __func__);
      double *image_pos = nullptr;
      
      float const elongation = 1.25f; // a bond elongated by 25% over his default length is still counted
      
      double const rcut = 5.11*Angstrom2Bohr; // maximum analysis range is 5 Angstrom
//       int const num_bins = 0; // no histogram
//       int const num_bins = 1 << 9; // 512 bins for 5.12 Angstrom
//       double const rcut = 4*5.11*Angstrom2Bohr; // maximum analysis range is 20 Angstrom
//       int const num_bins = 1 << 11; // 2048 bins for 20.48 Angstrom
//       double const bin_width = 0.01*Angstrom2Bohr;

      double const bin_width = 0.01*Angstrom2Bohr;
      int const num_bins = std::ceil(rcut/bin_width);
      
      if (echo > 4) {
          printf("# Bond search within interaction radius %.3f %s\n", rcut*Ang,_Ang);
          if (num_bins > 0) printf("# Distance histogram bin width is %.6f %s\n", bin_width*Ang,_Ang);
      } // echo
      
      int const nimages = boundary_condition::periodic_images(&image_pos, cell, bc, rcut, echo);
      
      auto ispecies = std::vector<int8_t>(natoms, 0); 
      auto occurrence = std::vector<int>(128, 0); 
      int8_t species_of_Z[128];
      int8_t Z_of_species[128];
      char  Sy_of_species[128][4];
      int nspecies = 0; // number of different species
      { // scope: fill ispecies and species_of_Z
          for(int first = 1; first >= 0; --first) {
              for(int ia = 0; ia < natoms; ++ia) {
                  int const Z_ia = std::round(xyzZ[ia*4 + 3]);
                  int const Z_ia_mod = Z_ia & 127;
                  if (first) {
                      ++occurrence[Z_ia_mod];
                  } else {
                      ispecies[ia] = species_of_Z[Z_ia_mod];
                  } // first
              } // ia
              if (first) {
                  // evaluate the histogram only in the first iteration
                  for(int Z = 0; Z < 128; ++Z) {
                      if (occurrence[Z] > 0) {
                          species_of_Z[Z] = nspecies;
                          Z_of_species[nspecies] = Z;
                          Sy_of_species[nspecies][0] = element_symbols[2*Z + 0];
                          Sy_of_species[nspecies][1] = element_symbols[2*Z + 1];
                          Sy_of_species[nspecies][2] = 0;
                          Sy_of_species[nspecies][3] = 0;
                          ++nspecies;
                      } // non-zero count
                  } // Z
              } // i10
          } // run twice

      } // scope
      
      if (echo > 2) {
          printf("# Found %d different elements: ", nspecies);
          for(int is = 0; is < nspecies; ++is) {
              printf("  %dx %s", occurrence[Z_of_species[is]], Sy_of_species[is]); 
          } // is
          printf("\n");
      } // echo
      
      double const inv_bin_width = 1./bin_width;
      auto const dist_hist = (num_bins > 0) ? new int[nspecies*nspecies][num_bins] : nullptr;
      set((int*) dist_hist, nspecies*nspecies*num_bins, 0); // clear
      auto bond_hist = std::vector<int>(nspecies*nspecies, 0);

      auto coordination_number = std::vector<uint8_t>(natoms, 0);
      float const too_large = 188.973;
      auto smallest_distance = std::vector<float>(nspecies*nspecies, too_large);

      typedef vector_math::vec<3,double> vec3;
      double const rcut2 = pow2(rcut);
      int npairs = 0, nbonds = 0;
      int nfar = 0, near = 0;
      for(int ia = 0; ia < natoms; ++ia) {
          vec3 const pos_ia = &xyzZ[ia*4];
          int const Z_ia = std::round(xyzZ[ia*4 + 3]);
          int const isi = ispecies[ia];
          assert(Z_ia == Z_of_species[isi]);
          if (echo > 6) printf("# [ia=%d] pos_ia = %g %g %g Z=%d\n", ia, pos_ia[0], pos_ia[1], pos_ia[2], Z_ia);
          for(int ja = 0; ja <= ia; ++ja) {
              vec3 const pos_ja = &xyzZ[ja*4];
              int const isj = ispecies[ja];
              int const Z_ja = Z_of_species[isj];
              auto const expected_bond_length = default_bond_length(Z_ia, Z_ja);
              if (echo > 7) printf("# [ia=%d, ja=%d] pos_ja = %g %g %g Z=%d\n", ia, ja, pos_ja[0], pos_ja[1], pos_ja[2], Z_ja);
              for(int ii = (ia == ja); ii < nimages; ++ii) { // start from 0 unless self-interaction
                  vec3 const pos_ii = &image_pos[ii*4];
                  vec3 const diff = pos_ja + pos_ii - pos_ia;
                  auto const d2 = norm(diff);
                  if (d2 > rcut2) {
                      ++nfar; // too far to be analyzed
                  } else {
                      ++near;
                      if (echo > 8) printf("# [%d,%d,%d] %g\n", ia, ja, ii, d2);
                      if (echo > 8) printf("%g\n", d2);
                      auto const dist = std::sqrt(d2);
                      int const ijs = isi*nspecies + isj;
                      int const jis = isj*nspecies + isi;
                      if (nullptr != dist_hist) {
                          int const ibin = (int)(dist*inv_bin_width);
                          assert(ibin < num_bins);
                          ++dist_hist[ijs][ibin];
                          ++dist_hist[jis][ibin];
                      }
                      if (dist < expected_bond_length*elongation) {
                          ++nbonds;
                          ++coordination_number[ia];
                          ++coordination_number[ja];
                          ++bond_hist[ijs];
                          ++bond_hist[jis];
                      } // atoms are close enough to assume a chemical bond
                      smallest_distance[ijs] = std::min(smallest_distance[ijs], (float)dist);
                      smallest_distance[jis] = std::min(smallest_distance[jis], (float)dist);
                  }
                  ++npairs;
              } // ii
          } // ja
      } // ia
      if (echo > 0) printf("# checked %d atom-atom pairs, %d near and %d far\n", npairs, near, nfar);
      assert((natoms*(natoms + 1))/2 * nimages - natoms == npairs);

      if (echo > 2) {
          if (num_bins > 0) {
              printf("\n## bond histogram (in %s)\n", _Ang);
              for(int ibin = 0; ibin < num_bins; ++ibin) {
                  float const dist = ibin*bin_width;
                  printf("%.3f ", dist*Ang);
                  for(int ijs = 0; ijs < nspecies*nspecies; ++ijs) {
                      printf(" %d", dist_hist[ijs][ibin]);
                  } // ijs
                  printf("\n");
              } // ibin
          } // num_bins > 0

          printf("\n# bond counts  ");
          for(int js = 0; js < nspecies; ++js) {
              printf("     %s ", Sy_of_species[js]); // create legend
          } // js
          printf("total = %d\n", nbonds);
          int check_nbonds = 0;
          for(int is = 0; is < nspecies; ++is) {
              printf("# bond count ");
              for(int js = 0; js < nspecies; ++js) {
                  check_nbonds += bond_hist[is*nspecies + js];
                  if (js >= is) {
                      printf("%8d", bond_hist[is*nspecies + js] >> (js == is)); // divide diagonal elements by 2
                  } else {
                      printf("        "); // do not show elements below the diagonal
                  }
              } // js
              printf("  %s\n", Sy_of_species[is]);
          } // is
          printf("# check total = %d vs %d\n", check_nbonds, 2*nbonds);
          printf("\n");
          assert(check_nbonds == 2*nbonds);
          
          printf("# shortest distances");
          for(int js = 0; js < nspecies; ++js) {
              printf("     %s ", Sy_of_species[js]); // create legend
          } // js
          printf("      in %s\n", _Ang);
          for(int is = 0; is < nspecies; ++is) {
              printf("# shortest distance ");
              for(int js = 0; js < nspecies; ++js) {
                  if (js >= is) {
                      if (smallest_distance[is*nspecies + js] < too_large) {
                          printf("%8.3f", smallest_distance[is*nspecies + js]*Ang);
                      } else {
                          printf("   n/a  "); // no distance below rcut found
                      }
                  } else {
                      printf("        "); // do not show elements below the diagonal
                  }
              } // js
              printf("  %s  in %s\n", Sy_of_species[is], _Ang);
          } // is
      } // echo
      
      // analyze coordination numbers
      if (echo > 3) {
          int cn_exceeds = 0;
          int const max_cn = 16;
          auto const cn_hist = new int[nspecies][max_cn];
          set((int*)cn_hist, nspecies*max_cn, 0);
          for(int ia = 0; ia < natoms; ++ia) {
              int const isi = ispecies[ia];
              int const cni = coordination_number[ia];
              if (cni < max_cn) ++cn_hist[isi][cni]; else ++cn_exceeds;
          } // ia
          printf("\n# coordination numbers (radius in %s)\n", _Ang);
          for(int is = 0; is < nspecies; ++is) {
              printf("# coordination number for %s (%.3f)", Sy_of_species[is], default_bond_length(Z_of_species[is])*Ang);
              for(int cn = 0; cn < max_cn; ++cn) {
                  if (cn_hist[is][cn] > 0) {
                      printf("  %dx%d", cn_hist[is][cn], cn);
                  } // histogram count non-zero
              } // cn
              printf("\n");
          } // is
          printf("\n");
          delete[] cn_hist;
      } // echo
      
      if (nullptr != dist_hist) delete[] dist_hist;
      delete[] image_pos;
      return 0;
  } // analysis
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_analysis(char const *filename="dna.xyz", int const echo=9) {

    double *xyzZ = nullptr;
    int natoms;
    status_t stat = read_xyz_file(&xyzZ, &natoms, filename, echo);

    double const cell[] = {63.872738414, 45.353423726, 45.353423726}; // DNA-cell in Bohr
    int const bc[] = {Periodic_Boundary, Isolated_Boundary, Isolated_Boundary};

    stat += analysis(xyzZ, natoms, cell, bc);
    
    delete[] xyzZ;
    return stat;
  } // test_analysis

  status_t all_tests() {
    auto status = 0;
    status += test_analysis();
    return status;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace geometry_analysis
