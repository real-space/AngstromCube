#include <cstdio> // printf
#include <cassert> // assert
#include <cmath> // std::abs
#include <vector> // std::vector<T>

#include "chemical_symbol.hxx"

#include "chemical_symbol.h" // char element_symbols[128*2]

namespace chemical_symbol {

  template<typename int_t>
  inline int8_t mod128(int_t const Z) { return static_cast<int8_t>(Z & 127); } // Z % (2^7)

  // public:
  status_t get(int const Z, char* Sy) {
      if (nullptr == Sy) return 2; // error code 2: result pointer is null
      auto const z7 = mod128(Z);
      Sy[0] = element_symbols[2*z7 + 0];
      char const y = element_symbols[2*z7 + 0];
      if (' ' == y) { Sy[1] = y; Sy[2] = 0; } else { Sy[1] = 0; }
      return ((Z < 0) || (Z > 127)); // error code 1: out of bounds error
  } // get
  
  int8_t map(char const y, char const ys[], int8_t const Zs[]) {
      int i = 0;
      for( ; 0 != ys[i]; ++i) {
          if (ys[i] == y) return Zs[i];
      } // i
      return Zs[i]; // error code must be the last entry
  } // map

  // public:
  int8_t decode(char const S, char const y) {
      switch (S) {
        // from here on, the code is auto-generated by test_decoding, see below
        case 'A': { int8_t const Zs[] = {13, 18, 33, 47, 79, 85, 89, 95, -65};
                    return map(y, "lrsgutcm", Zs); }
        case 'B': { int8_t const Zs[] = {4, 5, 35, 56, 83, 97, 107, -66};
                    return map(y, "e raikh", Zs); }
        case 'C': { int8_t const Zs[] = {6, 17, 20, 24, 27, 29, 48, 55, 58, 96, 98, 112, -67};
                    return map(y, " laroudsemfn", Zs); }
        case 'D': { int8_t const Zs[] = {66, 105, 110, -68};
                    return map(y, "ybs", Zs); }
        case 'E': { int8_t const Zs[] = {63, 68, 99, -69};
                    return map(y, "urs", Zs); }
        case 'F': { int8_t const Zs[] = {9, 26, 87, 100, -70};
                    return map(y, " erm", Zs); }
        case 'G': { int8_t const Zs[] = {31, 32, 64, -71};
                    return map(y, "aed", Zs); }
        case 'H': { int8_t const Zs[] = {1, 2, 67, 72, 80, 108, -72};
                    return map(y, " eofgs", Zs); }
        case 'I': { int8_t const Zs[] = {49, 53, 77, -73};
                    return map(y, "n r", Zs); }
        case 'K': { int8_t const Zs[] = {19, 36, -75};
                    return map(y, " r", Zs); }
        case 'L': { int8_t const Zs[] = {3, 57, 71, 103, -76};
                    return map(y, "iaur", Zs); }
        case 'M': { int8_t const Zs[] = {12, 25, 42, 101, 109, -77};
                    return map(y, "gnodt", Zs); }
        case 'N': { int8_t const Zs[] = {7, 10, 11, 28, 41, 60, 93, 102, -78};
                    return map(y, " eaibdpo", Zs); }
        case 'O': { int8_t const Zs[] = {8, 76, -79};
                    return map(y, " s", Zs); }
        case 'P': { int8_t const Zs[] = {15, 46, 59, 61, 78, 82, 84, 91, 94, -80};
                    return map(y, " drmtboau", Zs); }
        case 'R': { int8_t const Zs[] = {37, 44, 45, 75, 86, 88, 104, 111, -82};
                    return map(y, "buhenafg", Zs); }
        case 'S': { int8_t const Zs[] = {14, 16, 21, 34, 38, 50, 51, 62, 106, -83};
                    return map(y, "i cernbmg", Zs); }
        case 'T': { int8_t const Zs[] = {22, 43, 52, 65, 69, 73, 81, 90, -84};
                    return map(y, "icebmalh", Zs); }
        case 'U': { int8_t const Zs[] = {92, -85};
                    return map(y, " ", Zs); }
        case 'V': { int8_t const Zs[] = {23, -86};
                    return map(y, " ", Zs); }
        case 'W': { int8_t const Zs[] = {74, -87};
                    return map(y, " ", Zs); }
        case 'X': { int8_t const Zs[] = {54, -88};
                    return map(y, "e", Zs); }
        case 'Y': { int8_t const Zs[] = {39, 70, -89};
                    return map(y, " b", Zs); }
        case 'Z': { int8_t const Zs[] = {30, 40, -90};
                    return map(y, "nr", Zs); }
        case '_': { int8_t const Zs[] = {0, -95};
                    return map(y, "_", Zs); }
        case 'e': { int8_t const Zs[] = {127, -101};
                    return map(y, " ", Zs); }
        case 'u': { int8_t const Zs[] = {113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, -117};
                    return map(y, "tqphsond123456", Zs); }
        //   up to here, the code is auto-generated by test_decoding, see below
      } // switch (code)
      return -63; // not found in list, display error element '?'
      // in priciple, the error return values could be different depending of there is
      // no element with that char S as first character or if the second char y does not exist for that S.
  } // decode

// #define _WITH_14BIT_ENCODING
#ifdef  _WITH_14BIT_ENCODING
  // strategy: combine both printable ASCII char S,y into one 14bitnumber
  inline int16_t encode_14bit(char const S, char const y) { return (S - 48) * 128 + (y - 48); }

  int16_t encode_14bit(int const Z) {
      int8_t const z7 = mod128(Z);
      return encode_14bit(element_symbols[2*z7 + 0], element_symbols[2*z7 + 1]);
  } // encode_14bit

  int8_t decode_14bit(char const S, char const y) {
      switch (encode_14bit(S, y)) {
        // from here on, the code is auto-generated by test_decoding_14bit, see below
        case 6063: return 0; // "__"
        case 3056: return 1; // "H "
        case 3125: return 2; // "He"
        case 3641: return 3; // "Li"
        case 2357: return 4; // "Be"
        case 2288: return 5; // "B "
        case 2416: return 6; // "C "
        case 3824: return 7; // "N "
        case 3952: return 8; // "O "
        case 2800: return 9; // "F "
        case 3893: return 10; // "Ne"
        case 3889: return 11; // "Na"
        case 3767: return 12; // "Mg"
        case 2236: return 13; // "Al"
        case 4537: return 14; // "Si"
        case 4080: return 15; // "P "
        case 4464: return 16; // "S "
        case 2492: return 17; // "Cl"
        case 2242: return 18; // "Ar"
        case 3440: return 19; // "K "
        case 2481: return 20; // "Ca"
        case 4531: return 21; // "Sc"
        case 4665: return 22; // "Ti"
        case 4848: return 23; // "V "
        case 2498: return 24; // "Cr"
        case 3774: return 25; // "Mn"
        case 2869: return 26; // "Fe"
        case 2495: return 27; // "Co"
        case 3897: return 28; // "Ni"
        case 2501: return 29; // "Cu"
        case 5438: return 30; // "Zn"
        case 2993: return 31; // "Ga"
        case 2997: return 32; // "Ge"
        case 2243: return 33; // "As"
        case 4533: return 34; // "Se"
        case 2370: return 35; // "Br"
        case 3522: return 36; // "Kr"
        case 4402: return 37; // "Rb"
        case 4546: return 38; // "Sr"
        case 5232: return 39; // "Y "
        case 5442: return 40; // "Zr"
        case 3890: return 41; // "Nb"
        case 3775: return 42; // "Mo"
        case 4659: return 43; // "Tc"
        case 4421: return 44; // "Ru"
        case 4408: return 45; // "Rh"
        case 4148: return 46; // "Pd"
        case 2231: return 47; // "Ag"
        case 2484: return 48; // "Cd"
        case 3262: return 49; // "In"
        case 4542: return 50; // "Sn"
        case 4530: return 51; // "Sb"
        case 4661: return 52; // "Te"
        case 3184: return 53; // "I "
        case 5173: return 54; // "Xe"
        case 2499: return 55; // "Cs"
        case 2353: return 56; // "Ba"
        case 3633: return 57; // "La"
        case 2485: return 58; // "Ce"
        case 4162: return 59; // "Pr"
        case 3892: return 60; // "Nd"
        case 4157: return 61; // "Pm"
        case 4541: return 62; // "Sm"
        case 2757: return 63; // "Eu"
        case 2996: return 64; // "Gd"
        case 4658: return 65; // "Tb"
        case 2633: return 66; // "Dy"
        case 3135: return 67; // "Ho"
        case 2754: return 68; // "Er"
        case 4669: return 69; // "Tm"
        case 5298: return 70; // "Yb"
        case 3653: return 71; // "Lu"
        case 3126: return 72; // "Hf"
        case 4657: return 73; // "Ta"
        case 4976: return 74; // "W "
        case 4405: return 75; // "Re"
        case 4035: return 76; // "Os"
        case 3266: return 77; // "Ir"
        case 4164: return 78; // "Pt"
        case 2245: return 79; // "Au"
        case 3127: return 80; // "Hg"
        case 4668: return 81; // "Tl"
        case 4146: return 82; // "Pb"
        case 2361: return 83; // "Bi"
        case 4159: return 84; // "Po"
        case 2244: return 85; // "At"
        case 4414: return 86; // "Rn"
        case 2882: return 87; // "Fr"
        case 4401: return 88; // "Ra"
        case 2227: return 89; // "Ac"
        case 4664: return 90; // "Th"
        case 4145: return 91; // "Pa"
        case 4720: return 92; // "U "
        case 3904: return 93; // "Np"
        case 4165: return 94; // "Pu"
        case 2237: return 95; // "Am"
        case 2493: return 96; // "Cm"
        case 2363: return 97; // "Bk"
        case 2486: return 98; // "Cf"
        case 2755: return 99; // "Es"
        case 2877: return 100; // "Fm"
        case 3764: return 101; // "Md"
        case 3903: return 102; // "No"
        case 3650: return 103; // "Lr"
        case 4406: return 104; // "Rf"
        case 2610: return 105; // "Db"
        case 4535: return 106; // "Sg"
        case 2360: return 107; // "Bh"
        case 3139: return 108; // "Hs"
        case 3780: return 109; // "Mt"
        case 2627: return 110; // "Ds"
        case 4407: return 111; // "Rg"
        case 2494: return 112; // "Cn"
        case 8900: return 113; // "ut"
        case 8897: return 114; // "uq"
        case 8896: return 115; // "up"
        case 8888: return 116; // "uh"
        case 8899: return 117; // "us"
        case 8895: return 118; // "uo"
        case 8894: return 119; // "un"
        case 8884: return 120; // "ud"
        case 8833: return 121; // "u1"
        case 8834: return 122; // "u2"
        case 8835: return 123; // "u3"
        case 8836: return 124; // "u4"
        case 8837: return 125; // "u5"
        case 8838: return 126; // "u6"
        case 6768: return 127; // "e "
        //   up to here, the code is auto-generated by test_decoding_14bit, see below
      } // switch (code)
      return -1; // not found in list.
      // in priciple, the error return values could be different depending of there is
      // no element with that char S as first character or if the second char y does not exist for that S.
  } // decode_14bit
#endif
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_decoding(int const echo=1) { // echo>4: run code generation
      // code generator, not used
      if (echo < 1) return 0;
      printf("\n# %s %s\n", __FILE__, __func__);
      if (echo < 5) return 0;

      int constexpr mx = 16;
      std::vector<int8_t> iZ(128*mx, -1);
      std::vector<uint8_t> nZ(128, 0);

      for(int m = 0; m < 128; ++m) {
          char const S = element_symbols[2*m + 0];
          iZ[S*mx + nZ[S]] = m;
          ++nZ[S];
      } // m

      uint8_t max_nZ = 0;
      for(char S = 'A'; S <= 'z'; ++S) {
          if (nZ[S] > 0) {
              max_nZ = std::max(max_nZ, nZ[S]);
// code generation for using std::map<char,int8_t>
//               printf("# \'%c\' --> {", S);
//               for(int n = 0; n < nZ[S]; ++n) {
//                   auto const Z = iZ[S*mx + n];
//                   char const y = element_symbols[2*Z + 1];
//                   printf("\'%c\':%i, ", y, Z);
//               } // n
              printf("        case \'%c\': { int8_t const Zs[] = {", S);
              for(int n = 0; n < nZ[S]; ++n) {
                  auto const Z = iZ[S*mx + n];
                  printf("%i, ", Z);
              } // n
              printf("-%i};\n                    return map(y, \"", S);
              for(int n = 0; n < nZ[S]; ++n) {
                  auto const Z = iZ[S*mx + n];
                  char const y = element_symbols[2*Z + 1];
                  printf("%c", y);
              } // n
              printf("\", Zs); }\n");
          } // nZ[S] > 0
      } // S
      if (echo > 6) printf("\n# %s %s: largest number of y-chars is %i\n", __FILE__, __func__, max_nZ);
      return (mx < max_nZ);
  } // test_decoding

  status_t test_decoding_14bit(int const echo=2) {
      // check that encode_14bit(decode_14bit(x)) == x
      status_t stat = 0;
#ifdef  _WITH_14BIT_ENCODING
      if (echo > 1) printf("\n# %s %s\n", __FILE__, __func__);
      for(int Z = 0; Z < 128; ++Z) {
          char const S = element_symbols[2*Z + 0], y = element_symbols[2*Z + 1];
          if (echo > 5) printf("      case %i: return %i; // \"%c%c\"\n", encode_14bit(Z), Z, S,y);
          stat += (decode_14bit(S, y) != Z);
      } // Z
      if ((echo > 0) && (stat > 0)) printf("# %s %s failed for %i cases, run code generation again!\n", __FILE__, __func__, stat);
#endif
      return stat;
  } // test_decoding_14bit

  status_t test_consistency(int const echo=5) {
      // check that decode(ChemSymbol[Z]) == Z
      if (echo > 1) printf("\n# %s %s\n", __FILE__, __func__);
      status_t status = 0;
      for(int Z = 0; Z < 128; ++Z) {
          auto const stat = (decode(element_symbols[2*Z + 0], element_symbols[2*Z + 1]) != Z);
          if ((echo > 3) && (stat > 0)) printf("# %s %s failed for case Z=%i, run code generation again!\n", __FILE__, __func__, Z);
          status += (0 != stat);
      } // Z
      if ((echo > 0) && (status > 0)) printf("# %s %s failed for %i cases, run code generation again!\n", __FILE__, __func__, status);
      return status;
  } // test_consistency
  
  status_t all_tests(int echo) {
    auto status = 0;
    status += std::abs(test_decoding_14bit());
    status += std::abs(test_decoding());
    status += std::abs(test_consistency());
    return status;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace chemical_symbol
