#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::floor
#include <cstdint> // uint8_t
#include <fstream>
#include <sstream>
#include <string>

#include "geometry_analysis.hxx"

#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set, pow2
#include "chemical_symbol.h" // element_symbols

// #define FULL_DEBUG
// #define DEBUG

namespace geometry_analysis {
  
  
  

  float default_bond_length(int const Z1, int const Z2) {
    // data from http://chemwiki.ucdavis.edu/Theoretical_Chemistry/Chemical_Bonding/Bond_Order_and_Lengths
    uint8_t const bl[128] = {
     24   //        0
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
    ,255  //  Fr    87  // used to be 260
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
  } // def_bond_length
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_analysis(int const echo=5, char const *filename="dna.xyz") {

      double const Anstrom2Bohr = 1./0.52917724924;
      
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
          char Sy[2];
          double px, py, pz; // positions
          if (!(iss >> Sy >> px >> py >> pz)) { 
              printf("# failed parsing in %s:%d reads \"%s\"\n", filename, linenumber, line.c_str());
              break; // error
          }
          if (echo > 7) printf("# %c%c  %16.9f %16.9f %16.9f\n", Sy[0],Sy[1], px, py, pz);
          int iZ = -1; for(int Z = 0; Z < 128; ++Z) if ((Sy[0] == element_symbols[2*Z]) && (Sy[1] == element_symbols[2*Z + 1])) iZ = Z;
          xyzZ[na*4 + 0] = px*Anstrom2Bohr;
          xyzZ[na*4 + 1] = py*Anstrom2Bohr;
          xyzZ[na*4 + 2] = pz*Anstrom2Bohr;
          xyzZ[na*4 + 3] = iZ;
          ++na;
          assert(na <= natoms);
          // process pair 
      } // parse file line by line
            
//    auto const cell[] = {63.872738414, 45.353423726, 45.353423726}; // DNA-cell in Bohr, bc={periodic,isolated,isolated}

      return na - natoms; // returns 0 if the exact number of atoms has been found
  } // test_analysis

  status_t all_tests() {
    auto status = 0;
    status += test_analysis();
    return status;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace geometry_analysis
