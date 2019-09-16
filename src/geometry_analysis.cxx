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
#include "recorded_warnings.hxx" // warn

// #define FULL_DEBUG
// #define DEBUG

namespace geometry_analysis {
  
  template<typename int_t>
  class BoxStructure {
  private:
    int static constexpr X=0, Y=1, Z=2;
    std::vector<std::vector<int_t>> atom_index; // list of atomic indices
    std::vector<int> indirection;
    std::vector<int16_t> image_index;
    int64_t n_all_atoms; // number of all atoms
    int nboxes[3];
    int nhalo[3];
    
    template<int D> int inline get_center_box(int const j) { return (j + 999*nboxes[D]) % nboxes[D]; }
    int inline get_center_index(int const ix, int const iy, int const iz) {
        assert(ix >= 0); assert(ix < nboxes[X]);
        assert(iy >= 0); assert(iy < nboxes[Y]);
        assert(iz >= 0); assert(iz < nboxes[Z]);
        return (iz*nboxes[Y] + iy)*nboxes[X] + ix; }

  public:
      BoxStructure(double const cell[3], int const bc[3], 
                   double const radius, size_t const natoms, double const xyzZ[], int const echo=9) {
          assert(radius > 0);
          double box_size[3], inv_box_size[3];
          int hnh[3];
          int nbx = 1; // number of all boxes, including halos
          int mbx = 1; // number of central boxes, no halos
          for(int d = 0; d < 3; ++d) {
              assert(cell[d] > 0);
              nboxes[d] = std::max(1, (int)(cell[d]/radius));
              box_size[d] = cell[d]/nboxes[d];
              inv_box_size[d] = 1./box_size[d];
              nhalo[d] = (Periodic_Boundary == bc[d]) ? (int)std::ceil(radius*inv_box_size[d]) : 1;
              hnh[d] = nboxes[d] + 2*nhalo[d];
              nbx *= hnh[d]; // halo-enlarged
              mbx *= nboxes[d]; // only central boxes
          } // d
          if (echo > 2) printf("# %s: divide cell %.3f x %.3f x %.3f %s^3 into %d x %d x %d boxes\n", __func__, 
               cell[X]*Ang, cell[Y]*Ang, cell[Z]*Ang, _Ang, nboxes[X], nboxes[Y], nboxes[Z]);
          if (echo > 2) printf("# %s: box size %.3f x %.3f x %.3f %s^3 for interaction radius %.3f %s\n", __func__, 
               box_size[X]*Ang, box_size[Y]*Ang, box_size[Z]*Ang, _Ang, radius*Ang, _Ang);
          if (echo > 2) printf("# %s: use %d %d %d halo boxes on each side\n", __func__, nhalo[X], nhalo[Y], nhalo[Z]);

          indirection.resize(nbx);
          image_index.resize(4*nbx);
          { // scope: fill indirection list
              for        (int jz = 0; jz < hnh[Z]; ++jz) {  int const iz = get_center_box<Z>(jz - nhalo[Z]);
                  for    (int jy = 0; jy < hnh[Y]; ++jy) {  int const iy = get_center_box<Y>(jy - nhalo[Y]);
                      if (echo > 9) printf("# %s: indirection(iz=%2d, iy=%2d):", __func__, jz - nhalo[Z], jy - nhalo[Y]);                        
                      for(int jx = 0; jx < hnh[X]; ++jx) {  int const ix = get_center_box<X>(jx - nhalo[X]);
                          int ib = get_center_index(ix, iy, iz);
                          int const jb = (jz*hnh[Y] + jy)*hnh[X] + jx;
                          image_index[4*jb + X] = (int)std::floor((jx - nhalo[X])/((double)nboxes[X]));
                          image_index[4*jb + Y] = (int)std::floor((jy - nhalo[Y])/((double)nboxes[Y]));
                          image_index[4*jb + Z] = (int)std::floor((jz - nhalo[Z])/((double)nboxes[Z]));
                          if ((Isolated_Boundary == bc[X]) && (0 != image_index[4*jb + X])) ib = -1;
                          if ((Isolated_Boundary == bc[Y]) && (0 != image_index[4*jb + Y])) ib = -1;
                          if ((Isolated_Boundary == bc[Z]) && (0 != image_index[4*jb + Z])) ib = -1;
                          indirection[jb] = ib;
                          if (echo > 9) printf(" %2d", ib);
                      } // jx
                      if (echo > 9) printf("\n");                        
                  } // jy
              } // jz
          } // scope
          
          atom_index.resize(mbx);
          size_t const estimate = std::ceil(natoms/(mbx*.875));
          if (echo > 2) printf("# %s: use %ld atoms as box balance estimate\n", __func__, estimate);
          for(int ib = 0; ib < mbx; ++ib) {
              atom_index[ib].reserve(estimate);
          } // ib
          
          if (mbx > 1) {
              // start to sort the atomic positions into the boxes
              if (echo > 8) printf("# %s: inverse box size is %g %g %g\n", __func__, inv_box_size[X],inv_box_size[Y],inv_box_size[Z]);
              double min_coords[3] = {9e99, 9e99, 9e99};
              for(size_t ia = 0; ia < natoms; ++ia) {
                  for(int d = 0; d < 3; ++d) {
                      min_coords[d] = std::min(min_coords[d], xyzZ[ia*4 + d]);
                  } // d
              } // ia
              double const offset[3] = {-min_coords[X], -min_coords[Y], -min_coords[Z]};
              
              for(size_t ia = 0; ia < natoms; ++ia) {
                  int const iz = get_center_box<Z>((int)std::floor((xyzZ[ia*4 + Z] + offset[Z])*inv_box_size[Z]));
                  int const iy = get_center_box<Y>((int)std::floor((xyzZ[ia*4 + Y] + offset[Y])*inv_box_size[Y]));
                  int const ix = get_center_box<X>((int)std::floor((xyzZ[ia*4 + X] + offset[X])*inv_box_size[X]));
                  if (echo > 9) printf("# %s: atom #%ld Z=%.1f at %g %g %g goes into box %d %d %d\n", __func__,
                                ia, xyzZ[ia*4+3], xyzZ[ia*4+X],xyzZ[ia*4+Y],xyzZ[ia*4+Z], ix, iy, iz);
                  int const ib = get_center_index(ix, iy, iz);
                  atom_index[ib].push_back(ia);
              } // ia
          } else {
              if (echo > 2) printf("# %s: only 1 box, list is trivial\n", __func__);
              for(size_t ia = 0; ia < natoms; ++ia) {
                  atom_index[0].push_back(ia);
              } // ia
          } // mbx > 1
          
          // report box balance
          { // scope: get some stats
              double den = 0, sum = 0, squ = 0, mn = 9e99, mx = -mn;
              for(size_t ib = 0; ib < atom_index.size(); ++ib) {
                  double const val = atom_index[ib].size();
                  mn = std::min(mn, val); mx = std::max(mx, val);
                  den += 1; sum += val; squ += pow2(val);
              } // ib
              auto const avg = sum / den;
              auto const rms = std::sqrt(std::max(0., squ / den - avg*avg));
              if (echo > 2) printf("# %s: box balance min %d avg %.2f rms %.2f max %d\n", __func__, (int)mn, avg, rms, (int)mx);
              assert(natoms == sum);
          } // scope

      } // constructor
      
      size_t get_atom_list(int_t const **list, int jx, int jy, int jz, int *ii=nullptr) const {
          int const hnh[2] = {nboxes[X] + 2*nhalo[X], nboxes[Y] + 2*nhalo[Y]};
          assert(jx >= -nhalo[X]); assert(jx < nboxes[X] + nhalo[X]);
          assert(jy >= -nhalo[Y]); assert(jy < nboxes[Y] + nhalo[Y]);
          assert(jz >= -nhalo[Z]); assert(jz < nboxes[Z] + nhalo[Z]);
          int const jb = ((jz + nhalo[Z])*hnh[Y] + (jy + nhalo[Y]))*hnh[X] + (jx + nhalo[X]);
          for(int d = 0; d < 3; ++d) {
              if (ii != nullptr) ii[d] = image_index[4*jb + d];
          } // d
          int const ib = indirection[jb];
          *list = (ib >= 0) ? atom_index[ib].data() : nullptr;
          return  (ib >= 0) ? atom_index[ib].size() : 0;
      } // get_atom_list

      int get_number_of_boxes() const { return nboxes[X]*nboxes[Y]*nboxes[Z]; }
      int get_number_of_boxes(uint const D) const { assert(D < 3); return nboxes[D]; }
      inline int get_halo_thickness(uint const D) const { assert(D < 3); return nhalo[D]; }

  }; // class BoxStructure
  
  
  double constexpr Bohr2Angstrom = 0.52917724924;
  double constexpr Angstrom2Bohr = 1./Bohr2Angstrom;
  
  status_t read_xyz_file(double **xyzz, int *n_atoms, char const *filename, 
                         double *cell, int *bc, int const echo) { // optionals

      std::ifstream infile(filename, std::ifstream::in);
      int natoms = 0, linenumber = 2;
      infile >> natoms; // read the number of atoms
      std::string line;
      std::getline(infile, line);
      std::getline(infile, line);
      if (echo > 2) printf("# expect %d atoms, comment line in file: %s\n", natoms, line.c_str());
      { // scope parse comment line
          std::istringstream iss(line);
          std::string Cell, B[3];
          double L[3]; // positions
          iss >> Cell >> L[0] >> L[1] >> L[2] >> B[0] >> B[1] >> B[2];
          for(int d = 0; d < 3; ++d) {
              if (nullptr != cell) cell[d] = L[d] * Angstrom2Bohr;
              if (nullptr != bc) bc[d] = boundary_condition::fromString(B[d].c_str(), echo);
          } // d
      } // scope
      auto const xyzZ = (natoms > 0) ? new double[natoms*4] : nullptr;
      int na = 0;
      while (std::getline(infile, line)) {
          ++linenumber;
          std::istringstream iss(line);
          std::string Symbol;
          double px, py, pz; // positions
          if (!(iss >> Symbol >> px >> py >> pz)) {
              if (na < natoms) // we expect more atom lines
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

  
  status_t analysis(double const xyzZ[], int64_t const natoms, 
                    double const cell[3], int const bc[3], int const echo=6) {
      status_t stat = 0;
      if (echo > 1) printf("\n# %s:%s\n", __FILE__, __func__);
      double *image_pos = nullptr;
      
      float const elongation = 1.25f; // a bond elongated by 25% over his default length is still counted
      
      double const rcut = 5.11*Angstrom2Bohr; // maximum analysis range is 5 Angstrom
//       int const num_bins = 0; // no histogram
      int const num_bins = 1 << 9; // 512 bins for 5.12 Angstrom
//       double const rcut = 4*5.11*Angstrom2Bohr; // maximum analysis range is 20 Angstrom
//       int const num_bins = 1 << 11; // 2048 bins for 20.48 Angstrom
      double const bin_width = 0.01*Angstrom2Bohr;

//       double const bin_width = 0.02*Angstrom2Bohr;
//       int const num_bins = std::ceil(rcut/bin_width);
      
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
          printf("# Found %d different elements for %ld atoms: ", nspecies, natoms);
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
      int64_t nzero = 0, nstrange = 0;
      int64_t npairs = 0, nbonds = 0;
      int64_t nfar = 0, near = 0;
      for(int ia = 0; ia < natoms; ++ia) { // without BoxStructure
          vec3 const pos_ia = &xyzZ[ia*4];
          int const Z_ia = std::round(xyzZ[ia*4 + 3]);
          int const isi = ispecies[ia];
          assert(Z_ia == Z_of_species[isi]);
          if (echo > 6) printf("# [ia=%d] pos_ia = %g %g %g Z=%d\n", ia, pos_ia[0], pos_ia[1], pos_ia[2], Z_ia);
//           for(int ja = 0; ja <= ia; ++ja) { // without BoxStructure
          for(int ja = 0; ja < natoms; ++ja) { // without BoxStructure
              vec3 const pos_ja = &xyzZ[ja*4];
              int const isj = ispecies[ja];
              int const Z_ja = Z_of_species[isj];
              auto const expected_bond_length = default_bond_length(Z_ia, Z_ja);
              if (echo > 7) printf("# [ia=%d, ja=%d] pos_ja = %g %g %g Z=%d\n", ia, ja, pos_ja[0], pos_ja[1], pos_ja[2], Z_ja);
//               for(int ii = (ia == ja); ii < nimages; ++ii) { // start from 0 unless self-interaction
              for(int ii = 0; ii < nimages; ++ii) { // include self-interaction
                  vec3 const pos_ii = &image_pos[ii*4];
                  vec3 const diff = pos_ja + pos_ii - pos_ia;
                  auto const d2 = norm(diff);
                  if (d2 > rcut2) {
                      ++nfar; // too far to be analyzed
                  } else if (d2 < 1e-6) {
                      if (ia == ja) {
                          ++nzero; // ok - self interaction
                      } else ++nstrange; // ?
                      ++near;
                  } else {
                      ++near;
                      if (echo > 8) printf("# [%d,%d,%d] %g\n", ia, ja, ii, d2); // very verbose!!
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
      if (echo > 0) printf("# checked %.6f M atom-atom pairs, %.3f k near and %.6f M far\n", 1e-6*npairs, 1e-3*near, 1e-6*nfar);
//       assert((natoms*(natoms + 1))/2 * nimages == npairs + natoms);
      assert(natoms == nzero);
      assert(0 == nstrange);

      float minimum_distance = 9e9;
      for(int ijs = 0; ijs < nspecies*nspecies; ++ijs) {
          minimum_distance = std::min(minimum_distance, smallest_distance[ijs]);
      } // ijs
      if (minimum_distance < 1) { // 1 Bohr is reasonable to launch a warning
          ++stat;
          sprintf(warn, "Minimum distance between two atoms is %.1f %s", minimum_distance*Ang,_Ang);
      }
      
      if (echo > 2) {
          if (num_bins > 0) {
              printf("\n## bond histogram (in %s)\n", _Ang);
              for(int ibin = 0; ibin < num_bins; ++ibin) {
                  float const dist = ibin*bin_width;
                  printf("%.3f ", dist*Ang);
                  for(int ijs = 0; ijs < nspecies*nspecies; ++ijs) {
                      printf(" %d", dist_hist[ijs][ibin] >> 1); // divide by 2
                  } // ijs
                  printf("\n");
              } // ibin
          } // num_bins > 0

          printf("\n# bond counts  ");
          for(int js = 0; js < nspecies; ++js) {
              printf("     %s ", Sy_of_species[js]); // create legend
          } // js
          printf("total = %ld\n", nbonds);
          int64_t check_nbonds = 0;
          for(int is = 0; is < nspecies; ++is) {
              printf("# bond count ");
              for(int js = 0; js < nspecies; ++js) {
                  check_nbonds += bond_hist[is*nspecies + js];
                  if (js >= is) {
                      printf("%8d", bond_hist[is*nspecies + js] >> 1); // divide by 2
                  } else {
                      printf("        "); // do not show elements below the diagonal
                  }
              } // js
              printf("  %s\n", Sy_of_species[is]);
          } // is
          printf("# check total = %ld vs %ld\n", check_nbonds, 2*nbonds);
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
              int const cni = coordination_number[ia] >> 1;
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
          if (cn_exceeds > 0) {
              sprintf(warn, "In %d cases, the max. coordination (%d) was exceeded", cn_exceeds, max_cn);
              ++stat;
          }
          delete[] cn_hist;
      } // echo
      
      if (nullptr != dist_hist) delete[] dist_hist;
      delete[] image_pos;
      return stat;
  } // analysis

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  status_t analysis_box(double const xyzZ[], int64_t const natoms, 
                    double const cell[3], int const bc[3], int const echo=6) {
      status_t stat = 0;
      if (echo > 1) printf("\n# %s:%s\n", __FILE__, __func__);
      
      float const elongation = 1.25f; // a bond elongated by 25% over his default length is still counted
      
      double const rcut = 5.11*Angstrom2Bohr; // maximum analysis range is 5 Angstrom
      int const num_bins = 0; // no histogram
//       int const num_bins = 1 << 9; // 512 bins for 5.12 Angstrom
//       double const rcut = 4*5.11*Angstrom2Bohr; // maximum analysis range is 20 Angstrom
//       int const num_bins = 1 << 11; // 2048 bins for 20.48 Angstrom
      double const bin_width = 0.01*Angstrom2Bohr;

//       double const bin_width = 0.02*Angstrom2Bohr;
//       int const num_bins = std::ceil(rcut/bin_width);
      
      if (echo > 4) {
          printf("# Bond search within interaction radius %.3f %s\n", rcut*Ang,_Ang);
          if (num_bins > 0) printf("# Distance histogram bin width is %.6f %s\n", bin_width*Ang,_Ang);
      } // echo
      
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
          printf("# Found %d different elements for %ld atoms: ", nspecies, natoms);
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

      
      BoxStructure<int> box(cell, bc, 1*rcut, natoms, xyzZ);

      typedef vector_math::vec<3,double> vec3;
      double const rcut2 = pow2(rcut);
      int64_t nzero = 0, nstrange = 0;
      int64_t npairs = 0, nbonds = 0;
      int64_t nfar = 0, near = 0;
      for(int ibz = 0; ibz < box.get_number_of_boxes(2); ++ibz) {
      for(int iby = 0; iby < box.get_number_of_boxes(1); ++iby) {
      for(int ibx = 0; ibx < box.get_number_of_boxes(0); ++ibx) {
        int const *list_i;
        int const na_i = box.get_atom_list(&list_i, ibx, iby, ibz);
        for(int jbz = ibz - box.get_halo_thickness(2); jbz <= ibz + box.get_halo_thickness(2); ++jbz) {
        for(int jby = iby - box.get_halo_thickness(1); jby <= iby + box.get_halo_thickness(1); ++jby) {
        for(int jbx = ibx - box.get_halo_thickness(0); jbx <= ibx + box.get_halo_thickness(0); ++jbx) {
            int const *list_j;
            int shift[3];
            int const na_j = box.get_atom_list(&list_j, jbx, jby, jbz, shift);
            double const ppos_i[3] = {shift[0]*cell[0], shift[1]*cell[1], shift[2]*cell[2]}; // periodic shift if applicable
            vec3 const pos_ii = ppos_i; // convert to vector

            for(int iia = 0; iia < na_i; ++iia) { int const ia = list_i[iia];
              //========================================================================================================
              vec3 const pos_ia = &xyzZ[ia*4];
              int const Z_ia = std::round(xyzZ[ia*4 + 3]);
              int const isi = ispecies[ia];
              assert(Z_ia == Z_of_species[isi]);
              if (echo > 6) printf("# [ia=%d] pos_ia = %g %g %g Z=%d\n", ia, pos_ia[0], pos_ia[1], pos_ia[2], Z_ia);
              //========================================================================================================
              for(int ija = 0; ija < na_j; ++ija) { int const ja = list_j[ija];
                  //========================================================================================================
                  vec3 const pos_ja = &xyzZ[ja*4];
                  int const isj = ispecies[ja];
                  int const Z_ja = Z_of_species[isj];
                  auto const expected_bond_length = default_bond_length(Z_ia, Z_ja);
                  if (echo > 7) printf("# [ia=%d, ja=%d] pos_ja = %g %g %g Z=%d\n", ia, ja, pos_ja[0], pos_ja[1], pos_ja[2], Z_ja);
                  //========================================================================================================
                  vec3 const diff = pos_ja + pos_ii - pos_ia;
                  auto const d2 = norm(diff);
                  if (d2 > rcut2) {
                      ++nfar; // too far to be analyzed
                  } else if (d2 < 1e-6) {
                      if (ia == ja) {
                          ++nzero; // ok - self interaction
                      } else ++nstrange; // ?
                      ++near;
                  } else {
                      ++near;
                      if (echo > 8) printf("# [%d,%d,%d %d %d] %g\n", ia, ja, shift[0],shift[1],shift[2], d2); // very verbose!!
                      auto const dist = std::sqrt(d2);
                      int const ijs = isi*nspecies + isj;
                      if (nullptr != dist_hist) {
                          int const ibin = (int)(dist*inv_bin_width);
                          assert(ibin < num_bins);
                          ++dist_hist[ijs][ibin];
                      }
                      if (dist < expected_bond_length*elongation) {
                          ++nbonds;
                          ++coordination_number[ia];
                          ++bond_hist[ijs];
//                        if ((0 == isi) || (0 == isj)) printf("# hydrogen bond i#%d j#%d\n", (int)ia, (int)ja);
                      } // atoms are close enough to assume a chemical bond
                      smallest_distance[ijs] = std::min(smallest_distance[ijs], (float)dist);
                  }
                  ++npairs;

          } // ija
        } // iia
        
      } // jbx
      } // jby
      } // jbz
          
      } // ibx
      } // iby
      } // ibz
//       assert(natoms == nzero);
      if (natoms != nzero) {
          sprintf(warn, "Should find %ld exact zero distances but found %ld", natoms, nzero);
          ++stat;
      } // warn
      assert(0 == nstrange);
                  
      if (echo > 0) printf("# checked %.6f M atom-atom pairs, %.3f k near and %.6f M far\n", 1e-6*npairs, 1e-3*near, 1e-6*nfar);
//       assert((natoms*(natoms + 1))/2 * nimages == npairs + natoms);

      float minimum_distance = 9e9;
      for(int ijs = 0; ijs < nspecies*nspecies; ++ijs) {
          minimum_distance = std::min(minimum_distance, smallest_distance[ijs]);
      } // ijs
      if (minimum_distance < 1) { // 1 Bohr is reasonable to launch a warning
          ++stat;
          sprintf(warn, "Minimum distance between two atoms is %.1f %s", minimum_distance*Ang,_Ang);
      }
      
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
          printf("total = %ld\n", nbonds);
          int64_t check_nbonds = 0;
          for(int is = 0; is < nspecies; ++is) {
              printf("# bond count ");
              for(int js = 0; js < nspecies; ++js) {
                  check_nbonds += bond_hist[is*nspecies + js];
                  if (js >= is) {
                      printf("%8d", bond_hist[is*nspecies + js]);
                  } else {
                      printf("        "); // do not show elements below the diagonal
                  }
              } // js
              printf("  %s\n", Sy_of_species[is]);
          } // is
          printf("# check total = %ld vs %ld\n", 2*check_nbonds, 2*nbonds);
          printf("\n");
          assert(check_nbonds == nbonds);
          
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
          if (cn_exceeds > 0) {
              sprintf(warn, "In %d cases, the max. coordination (%d) was exceeded", cn_exceeds, max_cn);
              ++stat;
          }
          delete[] cn_hist;
      } // echo
      
      if (nullptr != dist_hist) delete[] dist_hist;
      return stat;
  } // analysis
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_analysis(int const echo=9) {

    double *xyzZ = nullptr;
    int natoms = 0;
    auto const filename = "gst.xyz";
    double cell[3] = {0, 0, 0}; 
    int bc[3] = {-7, -7, -7};
    status_t stat = read_xyz_file(&xyzZ, &natoms, filename, cell, bc, 0);
    
    if (echo > 2) printf("# found %d atoms in file \"%s\" with cell=[%.3f %.3f %.3f] %s and bc=[%d %d %d]\n",
                             natoms, filename, cell[0]*Ang, cell[1]*Ang, cell[2]*Ang, _Ang, bc[0], bc[1], bc[2]);

    stat += analysis_box(xyzZ, natoms, cell, bc);

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
