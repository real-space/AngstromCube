#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::floor
#include <cstdint> // uint8_t
#include <fstream> // std::fstream
#include <sstream> // std::sstream
#include <string> // std::string
#include <vector> // std::vector<T>

#include "geometry_analysis.hxx"

#include "boundary_condition.hxx" // Periodic_Boundary, Isolated_Boundary, ::periodic_images
#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set, pow2
#include "constants.hxx" // ::pi
#include "data_view.hxx" // view2D<T>
#include "inline_tools.hxx" // required_bits
#include "vector_math.hxx" // ::vec<n,T>
#include "chemical_symbol.h" // element_symbols
#include "recorded_warnings.hxx" // warn, error
#include "simple_timer.hxx" // SimpleTimer
#include "control.hxx" // ::get

#define FULL_DEBUG
#define DEBUG

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
                   double const radius, size_t const natoms, view2D<double> const & xyzZ, int const echo=0) {
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
              if (echo > 8) printf("# %s: inverse box size is %g %g %g\n", 
                                      __func__, inv_box_size[X], inv_box_size[Y], inv_box_size[Z]);
              double min_coords[3] = {9e99, 9e99, 9e99};
              for(size_t ia = 0; ia < natoms; ++ia) {
                  for(int d = 0; d < 3; ++d) {
                      min_coords[d] = std::min(min_coords[d], xyzZ[ia][d]);
                  } // d
              } // ia
              double const offset[3] = {-min_coords[X], -min_coords[Y], -min_coords[Z]};
              
              for(size_t ia = 0; ia < natoms; ++ia) {
                  int const iz = get_center_box<Z>((int)std::floor((xyzZ[ia][Z] + offset[Z])*inv_box_size[Z]));
                  int const iy = get_center_box<Y>((int)std::floor((xyzZ[ia][Y] + offset[Y])*inv_box_size[Y]));
                  int const ix = get_center_box<X>((int)std::floor((xyzZ[ia][X] + offset[X])*inv_box_size[X]));
                  if (echo > 9) printf("# %s: atom #%ld Z=%.1f at %g %g %g goes into box %d %d %d\n", __func__,
                                ia, xyzZ[ia][3], xyzZ[ia][X], xyzZ[ia][Y], xyzZ[ia][Z], ix, iy, iz);
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
  
  status_t read_xyz_file(view2D<double> & xyzZ, int & n_atoms, char const *filename, 
                         double *cell, int *bc, int const echo) { // optionals

      std::ifstream infile(filename, std::ifstream::in);
      if (infile.fail()) error("Unable to open file '%s' for reading coordinates", filename);

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
              if (nullptr != cell) {
                  cell[d] = L[d] * Angstrom2Bohr;
                  assert(cell[d] > 0);
              }
              if (nullptr != bc) bc[d] = boundary_condition::fromString(B[d].c_str(), echo, 120+d);
          } // d
      } // scope
      xyzZ = view2D<double>(natoms, 4);
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
          xyzZ(na,0) = px*Angstrom2Bohr;
          xyzZ(na,1) = py*Angstrom2Bohr;
          xyzZ(na,2) = pz*Angstrom2Bohr;
          xyzZ(na,3) = iZ;
          ++na;
          assert(na <= natoms);
          // process pair 
      } // parse file line by line
      
//    auto const cell[] = {63.872738414, 45.353423726, 45.353423726}; // DNA-cell in Bohr, bc={periodic,isolated,isolated}
      n_atoms = na;
      return na - natoms; // returns 0 if the exact number of atoms has been found
  } // read_xyz_file
  
  float default_bond_length(int const Z1, int const Z2=0) {
    // data originally from http://chemwiki.ucdavis.edu/Theoretical_Chemistry/Chemical_Bonding/Bond_Order_and_Lengths
    // can now be found for single bonds at
    // https://chem.libretexts.org/Ancillary_Materials/Reference/
    //              Reference_Tables/Atomic_and_Molecular_Properties/A3%3A_Covalent_Radii
    // which references https://doi.org/10.1002/chem.200901472 (Pekka Pyykk\"o and Michiko Atsumi)
    // Fe 26 with 110pm is good for bcc bulk
    // Fr 87 reduced from 260 to 255 to fit data format
    // Sc 21 reduced from 170 to 130 to make no Sc-Sc bonds inside Sc2O3
    uint8_t const bl[128] = {0, // now following Z=1..118 in Aco Z. Muradjan-ordering
      31,  28,                                                             // H  He
     128,  96,                                                             // Li Be
      84,  76,  73,  69,  71,  66,                                         // B  C  N  O  F  Ne
      57, 141,                                                             // Na Mg
     121, 111, 107, 105, 102, 106,                                         // Al Si P  S  Cl Ar
     203, 176,                                                             // K  Ca
     130, 160, 153, 139, 132, 110, 139, 124, 132, 122,                     // Sc Ti V  Cr Mn Fe Co Ni Cu Zn
     122, 120, 119, 120, 120, 116,                                         // Ga Ge As Se Br Kr
     220, 195,                                                             // Rb Sr
     190, 175, 164, 154, 147, 146, 142, 139, 145, 144,                     // Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd
     142, 139, 139, 138, 139, 140,                                         // In Sn Sb Te I  Xe
     244, 215,                                                             // Cs Ba
     207, 204, 203, 201, 199, 198, 198, 196, 194, 192, 192, 189, 190, 187, // La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb
     187, 175, 170, 162, 151, 144, 141, 136, 136, 132,                     // Lu Hf Ta W  Re Os Ir Pt Au Hg
     145, 146, 148, 140, 150, 150,                                         // Tl Pb Bi Po At Rn
     255, 221,                                                             // Fr Ra
     215, 206, 200, 196, 190, 187, 180, 169, 166, 168, 165, 167, 173, 176, // Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No
     161, 157, 149, 143, 141, 134, 129, 128, 121, 122,                     // Lr Rf Db Sg Bh Hs Mt Ds Rg Cn
     136, 143, 162, 175, 165, 157,                                         // ut uq up uh us uo
     159, 160,                                                             // un ud
     161, 162, 163, 164, 165, 166, 167};                                   // Z=121 .. Z=127 invented numbers
    float const picometer2Bohr = .01889726;
    return (bl[Z1] + bl[Z2]) * picometer2Bohr;
  } // default_bond_length
  
  
  typedef uint32_t index_t;
  class atom_image_index_t {
  public:
      index_t ia; // global index of the atom in the coordinate list
      int8_t is; // species index
      int8_t ix, iy, iz; // periodic shifts
      atom_image_index_t(index_t const atom_index=0, int const species_index=0, 
                         int const iix=0, int const iiy=0, int const iiz=0) 
      : ia(atom_index), is(species_index), ix(iix), iy(iiy), iz(iiz) {} // constructor
  }; // class atom_image_index_t

  
  template<typename real_t, typename int_t=short>
  int print_summary(char string[], // result
          real_t values[], int_t const na, // intput array (will be sorted on exit)
          double const grid_factor=1, double const display_factor=1, 
          char const mult_char='x', char const *fmt1=" %g",
          int const MaxBufferLen=-1, int const MaxEntryLen=24) {
      auto const string_start = string;
      bool const no_length_check = (MaxBufferLen < 0);

      std::sort(values, values + na);

      char fmtn[16]; { int const nfc = std::sprintf(fmtn, "%s%c%%d", fmt1, mult_char); assert(nfc < 16); }

      std::vector<int_t> ivalues(na);
      std::vector<int_t> occurrence(na, 1); // init all counters as 1, will be modified later
      for(int_t ia = 0; ia < na; ++ia) {
//           ivalues[ia] = int_t(values[ia]*grid_factor + 0.5); // integerize value
          ivalues[ia] = std::round(values[ia]*grid_factor); // integerize value
          if (ia > 0) {
              if (ivalues[ia - 1] == ivalues[ia]) { // values are the same on the integer grid
                  // transfer counts to the last one in a sequence of same values  
                  occurrence[ia] += occurrence[ia - 1];
                  occurrence[ia - 1] = 0;
              } // values are the same
          } // not for value #0
      } // ia

      int nbytes{0}; // init return value
      for(int_t ia = 0; ia < na; ++ia) {
          int_t  const count = occurrence[ia];
          double const value = ivalues[ia]*display_factor;
          if (count > 0) {
              if (no_length_check || ((string + MaxEntryLen) < (string_start + MaxBufferLen))) {
                  // write into the string
                  int const nchars = (count < 2) ? std::sprintf(string, fmt1, value):
                                            std::sprintf(string, fmtn, value, count);
                  assert( nchars <= MaxEntryLen );
                  nbytes += nchars;
                  string += nchars;
              } // there is still space in the string
          } // count
      } // ia
      return nbytes;
  } // print_summary
  
  
  template<typename real_t>
  void analyze_bond_structure(char* string, int const nb, real_t const bond_vectors[], float const Z) {
      string[0] = '\0'; if (nb < 1) return;
      view2D<real_t const> const bv(bond_vectors, 4); // wrap
//    string += sprintf(string, " coordination=%d", nb);
//       char const multiplicity_b = '^';
//       char const multiplicity_a = '*';
//       real_t max_len2 = 0, min_len2 = 9e99;

      std::vector<real_t> bond_length(nb);
      { // scope: compute all bond lengths
          for(int ib = 0; ib < nb; ++ib) {    auto const bvi = bv[ib];
              auto const d2i = pow2(bvi[0]) + pow2(bvi[1]) + pow2(bvi[2]);
//               max_len2 = std::max(max_len2, d2i);
//               min_len2 = std::min(min_len2, d2i);
              bond_length[ib] = std::sqrt(d2i);
//            printf(" %.3f", bond_length[ib]*Ang);
          } // ib
      } // scope

      // show the minimum and maximum bond length
//    string += sprintf(string, " [%.3f, %.3f]", std::sqrt(min_len2)*Ang, std::sqrt(max_len2)*Ang);

      int const na = (nb*(nb - 1))/2; // number of bond angles
      std::vector<real_t> bond_angle(na);
      { // scope: compute all possible bond angles
          int ia{0};
          for(int ib = 0; ib < nb; ++ib) {      auto const bvi = bv[ib];
              for(int jb = 0; jb < ib; ++jb) {  auto const bvj = bv[jb];
                  auto const dot = bvj[0]*bvi[0] + bvj[1]*bvi[1] + bvj[2]*bvi[2];
                  auto const bliblj = bond_length[ib]*bond_length[jb];
                  bond_angle[ia] = (-dot < bliblj) ? std::acos(dot/bliblj) 
                                                   : constants::pi;
                  ++ia;
//                printf(" %d", (int)(bond_angle[ia]*Deg));
              } // jb
          } // ib
          assert(na == ia);
      } // scope

      real_t constexpr Deg = 180/constants::pi; // Degrees
      string += print_summary(string, bond_length.data(), nb, Ang/.01, .01, '^', " %.2f"); // bin width 0.01 Ang
      string += std::sprintf(string, "  "); // some separator
      string += print_summary(string, bond_angle.data(), na, Deg/2., 2., '*'); // bin width 2 degrees
      
  } // analyze_bond_structure
  
  
  status_t analysis(view2D<double> const & xyzZ, index_t const natoms, 
                    double const cell[3], int const bc[3], int const echo=6) {
      status_t stat = 0;
      if (echo > 1) printf("\n# %s:%s\n", __FILE__, __func__);
      
      float const elongation = 1.25f; // a bond elongated by 25% over his default length is still counted
      if (echo > 3) printf("# Count interatomic distances up to %.1f%% of the default bond length as bond\n", (elongation - 1)*100);
      
      double const rcut = 6.03*Angstrom2Bohr; // maximum analysis range is 6 Angstrom
      double const bin_width = 0.02*Angstrom2Bohr;
      double const inv_bin_width = 1./bin_width;
      int const num_bins = (int)std::ceil(rcut*inv_bin_width); // 0:no histogram

      int constexpr MaxBP = 24; // maximum number of bond partners stored for later detailed analysis, 0:inactive
      std::vector<std::vector<atom_image_index_t>> bond_partner((MaxBP > 0)*natoms);
      index_t const natoms_BP = std::min(natoms, (index_t)2000); // limit the number of atoms for which the bonds are analyzed

      if (echo > 4) {
          printf("# Bond search within interaction radius %.3f %s\n", rcut*Ang,_Ang);
          if (num_bins > 0) printf("# Distance histogram bin width is %.6f %s\n", bin_width*Ang,_Ang);
      } // echo
      
      std::vector<int8_t> ispecies(natoms, 0); 
      std::vector<int> occurrence(128, 0); 
      int8_t species_of_Z[128];
      int8_t Z_of_species[128];
      char  Sy_of_species[128][4];    // examples: "Au\0", "H \0"
      char  Sy_of_species_right[128][4]; // right alignment: examples: "Au\0", " H\0"
      int nspecies{0}; // number of different species

      { // scope: fill ispecies and species_of_Z
          for(int first = 1; first >= 0; --first) {
              for(index_t ia = 0; ia < natoms; ++ia) {
                  int const Z_ia = std::round(xyzZ[ia][3]);
                  int const Z_ia_mod = Z_ia & 127; // modulo 128
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
                          char const S = element_symbols[2*Z + 0],
                                     y = element_symbols[2*Z + 1];
                          set(Sy_of_species[nspecies], 4, char(0));
                          Sy_of_species[nspecies][0] = S;
                          Sy_of_species[nspecies][1] = y;
                          set(Sy_of_species_right[nspecies], 4, char(0));
                          Sy_of_species_right[nspecies][0] = (' ' == y) ? ' ' : S;
                          Sy_of_species_right[nspecies][1] = (' ' == y) ?  S  : y;
                          ++nspecies;
                      } // non-zero count
                  } // Z
              } // i10
          } // run twice

      } // scope
      
      if (echo > 2) {
          printf("# Found %d different elements for %d atoms: ", nspecies, natoms);
          for(int is = 0; is < nspecies; ++is) {
              printf("  %dx %s", occurrence[Z_of_species[is]], Sy_of_species[is]); 
          } // is
          printf("\n");
      } // echo
      
      auto dist_hist = std::vector<uint16_t>(num_bins*nspecies*nspecies, 0);
      auto bond_hist = std::vector<int>(nspecies*nspecies, 0);

      auto coordination_number = std::vector<uint8_t>(natoms, 0);
      int constexpr MAX_coordination_number = std::numeric_limits<uint8_t>::max();
      float const too_large = 188.973;
      auto smallest_distance = std::vector<float>(nspecies*nspecies, too_large);

      view2D<double> image_pos;
      view2D<int8_t> image_shift;
      typedef vector_math::vec<3,double> vec3;
      double const rcut2 = pow2(rcut);
      int64_t nzero = 0, nstrange = 0, npairs = 0, nbonds = 0, nfar = 0, near = 0; // init counters
      int64_t bp_exceeded = 0, bp_truncated = 0; // init counters

// #define GEO_ORDER_N2
#ifdef  GEO_ORDER_N2
      int const nimages = boundary_condition::periodic_images(image_pos, cell, bc, rcut, echo, &image_shift);
      if (echo > 2) printf("# use N^2-algorithm, expect to visit %.1e atom pairs\n", nimages*pow2(double(natoms)));

      {{{{{ // open 5 scopes because the box structure has 5 loops
      for(int ii = 0; ii < nimages; ++ii) { // includes self-interaction
          vec3 const pos_ii = image_pos[ii];
          auto const shift = image_shift[ii]; // not tested
          for(int ia = 0; ia < natoms; ++ia) {
              //========================================================================================================
#else
      BoxStructure<index_t> box(cell, bc, rcut, natoms, xyzZ);

      for(int ibz = 0; ibz < box.get_number_of_boxes(2); ++ibz) {
      for(int iby = 0; iby < box.get_number_of_boxes(1); ++iby) {
      for(int ibx = 0; ibx < box.get_number_of_boxes(0); ++ibx) {
        index_t const *list_i;
        int const na_i = box.get_atom_list(&list_i, ibx, iby, ibz);
        for(int jbz = ibz - box.get_halo_thickness(2); jbz <= ibz + box.get_halo_thickness(2); ++jbz) {
        for(int jby = iby - box.get_halo_thickness(1); jby <= iby + box.get_halo_thickness(1); ++jby) {
        for(int jbx = ibx - box.get_halo_thickness(0); jbx <= ibx + box.get_halo_thickness(0); ++jbx) {
          index_t const *list_j;
          int shift[3];
          int const na_j = box.get_atom_list(&list_j, jbx, jby, jbz, shift);
          double const ppos_i[3] = {shift[0]*cell[0], shift[1]*cell[1], shift[2]*cell[2]}; // periodic shift if applicable
          vec3 const pos_ii = ppos_i; // convert to vector

          for(int iia = 0; iia < na_i; ++iia) { index_t const ia = list_i[iia];
              //========================================================================================================
#endif
              //========================================================================================================
              vec3 const pos_ia = xyzZ[ia];
              int const isi = ispecies[ia];
              int const Z_ia = Z_of_species[isi];
              if (echo > 8) printf("# [ia=%d] pos_ia = %g %g %g Z=%d\n", ia, pos_ia[0],pos_ia[1],pos_ia[2], Z_ia);
              //========================================================================================================
              vec3 const pos_ii_minus_ia = pos_ii - pos_ia;
#ifdef  GEO_ORDER_N2
              for(int ja = 0; ja < natoms; ++ja) {
#else
              for(int ija = 0; ija < na_j; ++ija) { index_t const ja = list_j[ija];
#endif
                  //========================================================================================================
                  vec3 const pos_ja = xyzZ[ja];
                  int const isj = ispecies[ja];
                  int const Z_ja = Z_of_species[isj];
                  if (echo > 9) printf("# [ia=%d, ja=%d] pos_ja = %g %g %g Z=%d\n", ia, ja, pos_ja[0],pos_ja[1],pos_ja[2], Z_ja);
                  //========================================================================================================
                  vec3 const diff = pos_ja + pos_ii_minus_ia;
                  auto const d2 = norm(diff);
                  if (d2 >= rcut2) {
                      ++nfar; // too far to be analyzed
                  } else if (d2 < 1e-6) {
                      if (ia == ja) {
                          ++nzero; // ok - self interaction
                      } else ++nstrange; // ?
                      ++near;
                  } else {
                      ++near;
//                    if (echo > 8) printf("# [%d,%d,%d] %g\n", ia, ja, ii, d2); // very verbose!!
//                    if (echo > 8) printf("# [%d,%d,%d %d %d] %g\n", ia, ja, shift[0],shift[1],shift[2], d2); // very verbose!!
                      auto const dist = std::sqrt(d2);
                      int const ijs = isi*nspecies + isj;
                      if (num_bins > 0) {
                          int const ibin = (int)(dist*inv_bin_width);
                          assert(ibin < num_bins);
                          ++dist_hist[ibin*nspecies*nspecies + ijs];
                      }
                      if (dist < elongation*default_bond_length(Z_ia, Z_ja)) {
                          ++nbonds;
                          if (MaxBP > 0) { // storing of bond partners active
                              int const cn = coordination_number[ia];
                              if (0 == cn) bond_partner[ia].clear();
                              if (cn < MaxBP) {
                                  if (ia < natoms_BP) {
                                      bond_partner[ia].push_back(atom_image_index_t(ja, isj, shift[0], shift[1], shift[2]));
                                  } else ++bp_truncated;
                              } else ++bp_exceeded;
                          } // MaxBP > 0
                          ++coordination_number[ia];
                          assert( coordination_number[ia] <= MAX_coordination_number );
                          ++bond_hist[ijs];
//                           if (echo > 2) printf("# bond between a#%d %s-%s a#%d  %g %s\n", 
//                             ia, Sy_of_species[isi], Sy_of_species[isj], ja, dist*Ang,_Ang);
                      } // atoms are close enough to assume a chemical bond
                      smallest_distance[ijs] = std::min(smallest_distance[ijs], (float)dist);
                  }
                  ++npairs;
                  //========================================================================================================
              } // ja
          } // ia
      } // ii
      }}}}} // close 5 loops for the box structure
      
      if (echo > 2) printf("# checked %.6f M atom-atom pairs, %.3f k near and %.6f M far\n", 1e-6*npairs, 1e-3*near, 1e-6*nfar);
      if (natoms != nzero) {
          warn("Should find %d exact zero distances but found %ld", natoms, nzero);
          ++stat;
      } // warn
      assert(natoms == nzero);
      assert(0 == nstrange);

      float minimum_distance = 9e9;
      for(int ijs = 0; ijs < nspecies*nspecies; ++ijs) {
          minimum_distance = std::min(minimum_distance, smallest_distance[ijs]);
      } // ijs
      if (minimum_distance < 1) { // 1 Bohr is reasonable to launch a warning
          ++stat;
          warn("Minimum distance between two atoms is %.1f %s", minimum_distance*Ang,_Ang);
      }
      
      if (echo > 5) {
        
          if (num_bins > 0) {
              printf("\n## distance histogram (in %s)\n", _Ang);
              int last_bins[4] = {0, -1, -1, -1}; // show the first bin always, -1: do not show
              for(int ibin = 0; ibin < num_bins; ++ibin) {
#if 1
                  bool nonzero = false; // sparsify the plotting of the distance histogram
                  for(int ijs = 0; ijs < nspecies*nspecies; ++ijs) {
                      nonzero = nonzero || (dist_hist[ibin*nspecies*nspecies + ijs] > 0); // analyze dist_hist[ibin]
                  } // ijs
                  if (nonzero) {
                      last_bins[(ibin + 1) & 3] = ibin - 1;
                      last_bins[(ibin + 2) & 3] = ibin; // also plot bins to the left and right of this nonzero line
                      last_bins[(ibin + 3) & 3] = ibin + 1;
                  } // nonzero
                  int const jbin = last_bins[ibin & 3];
                  last_bins[ibin & 3] = -1; // clear
#else
                  int const jbin = ibin;
#endif
                  if (jbin >= 0) {
                      float const dist = (jbin + 0.5)*bin_width; // display the center of the bin
                      printf("%.3f ", dist*Ang);
                      for(int ijs = 0; ijs < nspecies*nspecies; ++ijs) {
                          printf(" %d", dist_hist[jbin*nspecies*nspecies + ijs]);
                      } // ijs
                      printf("\n");
                  } // non-zero or before non-zero or after non-zero
              } // ibin
          } // num_bins > 0

      } // echo

      if (echo > 2) {
      
          printf("\n# bond counts ");
          for(int js = 0; js < nspecies; ++js) {
              printf("     %s ", Sy_of_species_right[js]); // create legend
          } // js
          std::vector<int> spec_sum(nspecies, 0);
          printf("\n");
          int64_t check_nbonds = 0;
          for(int is = 0; is < nspecies; ++is) {
              int row_sum{0};
              printf("# bond count ");
              for(int js = 0; js < nspecies; ++js) {
                  check_nbonds += bond_hist[is*nspecies + js];
                  spec_sum[js] += bond_hist[is*nspecies + js];
                  row_sum      += bond_hist[is*nspecies + js];
                  if (js >= is) {
                      printf("%8d", bond_hist[is*nspecies + js]);
                  } else {
                      printf("        "); // do not show elements below the diagonal
                  } // show only the upper triangular matrix
              } // js
              printf("  %s sum= %d\n", Sy_of_species[is], row_sum);
          } // is
          printf("# bond counts");
          for(int js = 0; js < nspecies; ++js) {
              printf("%8d", spec_sum[js]);
          } // js
          printf("   total= %ld\n", nbonds);
          if (check_nbonds != nbonds) error("Checksum for bonds does not agree: %d vs %d", check_nbonds, nbonds);
          printf("\n");

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
          int cn_exceeds{0};
          int const max_cn = 24;
          view2D<int> cn_hist(nspecies, max_cn, 0);
          for(index_t ia = 0; ia < natoms; ++ia) {
              int const isi = ispecies[ia];
              int const cni = coordination_number[ia];
              if (cni < max_cn) ++cn_hist(isi,cni); else ++cn_exceeds;
          } // ia
          printf("\n# coordination numbers (radius in %s)\n", _Ang);
          for(int is = 0; is < nspecies; ++is) {
              printf("# coordination number for %s (%.3f)", Sy_of_species[is], default_bond_length(Z_of_species[is])*Ang);
              for(int cn = 0; cn < max_cn; ++cn) {
                  if (cn_hist(is,cn) > 0) {
                      printf("  %dx%d", cn_hist(is,cn), cn);
                  } // histogram count non-zero
              } // cn
              printf("\n");
          } // is
          printf("\n");
          if (cn_exceeds > 0) {
              warn("In %d cases, the max. coordination (%d) was exceeded", cn_exceeds, max_cn);
              ++stat;
          } // cn_exceeds
      } // echo
     

      if (echo > 5) {
          // analyze local bond structure
          if (bond_partner.size() > 0) {
              if (echo > 2) printf("# show a bond structure analysis:\n"
                                   "# bond lengths are in %s, angles in degree\n\n",_Ang);
// #pragma omp parallel for
              for(index_t ia = 0; ia < natoms_BP; ++ia) {
                  int const cn = std::min((int)coordination_number[ia], MaxBP);
                  double const xyz_ia[3] = {xyzZ[ia][0], xyzZ[ia][1], xyzZ[ia][2]}; // load center coordinates
                  assert(cn == bond_partner[ia].size());
                  view2D<float> coords(cn, 4, 0.0); // get memory, decide float or double for the bond structure analysis
                  for(int ip = 0; ip < cn; ++ip) {
                      auto const &partner = bond_partner[ia][ip];
                      auto const ja = partner.ia; assert(ja >= 0); 
                      coords[ip][0] = xyzZ[ja][0] + partner.ix*cell[0] - xyz_ia[0];
                      coords[ip][1] = xyzZ[ja][1] + partner.iy*cell[1] - xyz_ia[1];
                      coords[ip][2] = xyzZ[ja][2] + partner.iz*cell[2] - xyz_ia[2];
                      coords[ip][3] = xyzZ[ja][3]; // Z of the bond partner, if needed
                  } // ip
                  int const isi = ispecies[ia];
                  bond_partner[ia].clear(); // free memory
                  char string_buffer[2048];
                  analyze_bond_structure(string_buffer, cn, coords.data(), xyzZ[4][3]);
// #pragma omp critical
                  if (echo > 4) printf("# a#%i %s %s\n", ia, Sy_of_species[isi], string_buffer); // no new line
              } // ia
          } // bond_partner
          if (bp_exceeded > 0 && MaxBP > 0) {
              stat += 0 < warn("In %ld cases, the max. number of bond partners (MaxBP=%d) was exceeded", bp_exceeded, MaxBP);
          } // maximum number of bond partners was exceeded
          if (bp_truncated > 0) {
              stat += 0 < warn("Bond partner analysis is performed only for the first %d atoms", natoms_BP);
          } // the number of atoms is larger than the max. number of atoms for which a bond structure analysis is done
      } // echo
      
      return stat;
  } // analysis
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_analysis(int const echo=9) {
    
    view2D<double> xyzZ;
    int natoms = 0;
    auto const geo_file = control::get("geometry.file", "atoms.xyz");
    double cell[3] = {0, 0, 0}; 
    int bc[3] = {-7, -7, -7};
    status_t stat = read_xyz_file(xyzZ, natoms, geo_file, cell, bc, 0);
    
    if (echo > 2) printf("# found %d atoms in file \"%s\" with cell=[%.3f %.3f %.3f] %s and bc=[%d %d %d]\n",
                             natoms, geo_file, cell[0]*Ang, cell[1]*Ang, cell[2]*Ang, _Ang, bc[0], bc[1], bc[2]);
    { // SimpleTimer timer(__FILE__, __LINE__, "analysis");
        stat += analysis(xyzZ, natoms, cell, bc, echo);
    } // timer

    return stat;
  } // test_analysis

  status_t all_tests(int const echo) {
    status_t status(0);
    status += test_analysis(echo);
    return status;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace geometry_analysis
