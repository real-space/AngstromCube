#include <cstdio> // std::printf, std::snprintf
#include <cmath> // std::floor
#include <cstdint> // int8_t
#include <vector> // std::vector<T>
#include <cassert> // assert

#include "sigma_config.hxx" // ::element_t, ::get, ::all_tests

#include "status.hxx" // status_t
#include "control.hxx" // ::get
#include "inline_math.hxx" // set
#include "chemical_symbol.hxx" // ::get
#include "recorded_warnings.hxx" // warn, error
#include "spherical_state.hxx" // core, semicore, valence, csv_undefined, csv_name

namespace sigma_config {

  char const ellchar[12] = "spdfghijkl?";

#define EXPERIMENTAL

  char const * default_config(unsigned const iZ) { // compiled into the code

        switch (iZ) {
            case  1: return "1s* 1 0 2p | 0.9 sigma .247";                              // H
            case  2: return "1s* 2 2p | 1.5 sigma .48";                                 // He
            case  3: return "2s* 1 0 2p 2e-99 | 2.0 sigma 1";                           // Li
            case  4: return "2s* 2 2p 2e-99 0 | 1.5 sigma .45";                         // Be
            case  5: return "2s* 2 2p* 1 0 3d | 1.2 sigma .45";                         // B
            case  6: return "2s* 2 2p* 2 0 3d | 1.2 sigma .43";                         // C
            case  7: return "2s* 2 2p* 3 0 3d | 1.0 sigma .33";                         // N
            case  8: return "2s* 2 2p* 3 1 3d | 1.13 sigma .322";                       // O
            case  9: return "2s* 2 2p* 3 2 3d | 1.2 sigma .323";                        // F
            case 10: return "2s* 2 2p* 6 3d | 1.8 sigma .564";                          // Ne
            case 11: return "2s 2 3s 1 0 2p 6 3p 2e-99 3d | 2.27 sigma .69";            // Na
//          case 12: return "3s* 2 3p 2e-99 3d | 1.96 sigma .41";                       // Mg (2 e)
            case 12: return "2s 2 3s 2 2p 6 3p 2e-99 3d | 1.96 sigma .41";              // Mg (10 e)
            case 13: return "3s* 2 3p* 1 0 3d | 2.05 sigma .645";                       // Al
            case 14: return "3s* 2 3p* 2 0 3d | 2.0 sigma .643";                        // Si
            case 15: return "3s* 2 3p* 3 0 3d | 1.8 sigma .512";                        // P
            case 16: return "3s* 2 3p* 3 1 3d | 1.7 sigma .535";                        // S
            case 17: return "3s* 2 3p* 3 2 3d | 1.5 sigma .503";                        // Cl
            case 18: return "3s* 2 3p* 6 3d | 1.6 sigma .546";                          // Ar
            case 19: return "3s 2 4s 1 0 3p* 6 3d | 1.77 sigma .466";                   // K
            case 20: return "3s 2 4s 2 3p* 6 3d | 1.77 sigma .487";                     // Ca
            case 21: return "3s 2 4s 2 3p 6 4p 2e-99 3d* 1 0 | 2.32 sigma .58";         // Sc
            case 22: return "3s 2 4s 2 3p 6 4p 2e-99 3d* 2 0 | 2.0 sigma .58";          // Ti
            case 23: return "3s 2 4s 2 3p 6 4p 2e-99 3d* 3 0 | 2.4 sigma .56";          // V
            case 24: return "4s* 1 0 4p* 2e-99 3d* 5 0 | 2.1 sigma .667";               // Cr
            case 25: return "3s 2 4s 2 3p 6 4p 2e-99 3d* 5 0 | 2.41 sigma .554";        // Mn
            case 26: return "4s* 2 4p* 2e-99 3d* 5 1 | 2.0 sigma .65";                  // Fe
            case 27: return "4s* 2 4p* 2e-99 3d* 5 2 | 1.9 sigma .608";                 // Co
            case 28: return "4s* 2 3p 6 4p 2e-99 3d* 5 3 | 2.15 sigma .48";             // Ni
            case 29: return "4s* 1 0 4p* 2e-99 3d* 10 | 2.0 sigma .61";                 // Cu
            case 30: return "4s* 2 4p* 2e-99 3d 10 | 2.23 sigma .596 V=parabola";       // Zn
            case 31: return "4s* 2 4p* 1 0 | 2.2 sigma .686";                           // Ga
            case 32: return "4s* 2 4p* 2 0 4d | 1.9 sigma .606";                        // Ge
            case 33: return "4s* 2 4p* 3 0 4d | 2.0 sigma .62";                         // As
            case 34: return "4s* 2 4p* 3 1 4d | 1.6 sigma .521";                        // Se
            case 35: return "4s* 2 4p* 3 2 4d | 2.1 sigma .6";                          // Br
            case 36: return "4s* 2 4p* 6 4d | 2.2 sigma .61";                           // Kr
            case 37: return "4s 2 5s 1 0 4p* 6 4d | 2.3 sigma .78";                     // Rb
            case 38: return "4s 2 5s 2 4p* 6 4d | 2.37 sigma .666";                     // Sr
            case 39: return "4s 2 5s 2 4p 6 5p 2e-99 4d 1 0 | 2.43 sigma .9";           // Y
            case 40: return "4s 2 5s 2 4p 6 5p 2e-99 4d* 2 0 | 2.35 sigma .58";         // Zr
            case 41: return "4s 2 5s 1 0 4p 6 5p 2e-99 4d* 4 0 | 2.35 sigma .59";       // Nb
            case 42: return "4s 2 5s 1 0 4p 6 5p 2e-99 4d* 5 0 | 2.34 sigma .585";      // Mo
            case 43: return "4s 2 5s 1 0 4p 6 5p 2e-99 4d* 5 1 | 2.4 sigma .58";        // Tc
            case 44: return "4s 2 5s 1 0 4p 6 5p 2e-99 4d* 5 2 | 2.37 sigma .571";      // Ru
            case 45: return "5s* 1 0 4p 6 5p 2e-99 4d* 5 3 | 2.35 sigma .91";           // Rh
            case 46: return "4s 2 5s 2e-99 0 4p* 6 4d* 10 | 2.32 sigma .904";           // Pd
            case 47: return "5s* 1 0 4p* 6 4d* 10 | 2.23 sigma .608";                   // Ag
            case 48: return "5s* 2 5p* 2e-99 4d* 10 | 2.2 sigma .6";                    // Cd
            case 49: return "5s* 2 5p* 1 0 4d 10 | 2.17 sigma .7 V=parabola";           // In
            case 50: return "5s* 2 5p* 2 0 4d 10 | 2.24 sigma .755 V=parabola";         // Sn
            case 51: return "5s* 2 5p* 3 0 4d 10 | 2.18 sigma .73 V=parabola";          // Sb
            case 52: return "5s* 2 5p* 3 1 5d | 2.23 sigma .654";                       // Te
            case 53: return "5s* 2 5p* 3 2 5d | 2.2 sigma .68";                         // I
            case 54: return "5s* 2 5p* 6 5d | 2.24 sigma .62";                          // Xe
            case 55: return "5s 2 6s 1 0 5p* 6 5d | 2.2 sigma .61";                     // Cs
            case 56: return "5s 2 6s 2 5p* 6 5d | 2.2 sigma .645";                      // Ba
            case 57: return "5s 2 6s 2 5p* 6 5d 1 0 4f | 2.45 sigma .59";               // La
#ifdef EXPERIMENTAL
            case 58: return "5s 2 6s 2 5p* 6 5d 0 4f 2 0 | 2.335 sigma .65";            // Ce
            case 59: return "5s 2 6s 2 5p* 6 5d 0 4f 3 0 | 2.49 sigma .686";            // Pr
            case 60: return "5s 2 6s 2 5p* 6 5d 0 4f 4 0 | 2.45 sigma .7";              // Nd
            case 61: return "5s 2 6s 2 5p* 6 5d 0 4f 5 0 | 2.45 sigma .72";             // Pm
            case 62: return "5s 2 6s 2 5p* 6 5d 0 4f 6 0 | 2.516 sigma .78";            // Sm
            case 63: return "5s 2 6s 2 5p* 6 5d 0 4f 7 0 | 2.559 sigma .823";           // Eu
            case 64: return "5s 2 6s 2 5p* 6 5d 0 4f 7 1 | 2.417 sigma .75";            // Gd
            case 65: return "5s 2 6s 2 5p* 6 5d 0 4f 7 2 | 2.37 sigma .736";            // Tb
            case 66: return "5s 2 6s 2 5p* 6 5d 0 4f 7 3 | 2.392 sigma .757";           // Dy
            case 67: return "5s 2 6s 2 5p* 6 5d 0 4f 7 4 | 2.37 sigma .75";             // Ho
            case 68: return "5s 2 6s 2 5p* 6 5d 0 4f 7 5 | 2.375 sigma .753";           // Er
            case 69: return "5s 2 6s 2 5p* 6 5d 0 4f 7 6 | 2.393 sigma .762";           // Tm
            case 70: return "5s 2 6s 2 5p* 6 5d 0 4f 7 7 | 2.4 sigma .761";             // Yb
#endif // EXPERIMENTAL
            case 71: return "5s 2 6s 2 5p* 6 5d 1 0 4f 14 | 2.25 sigma .6";             // Lu
            case 72: return "5s 2 6s 2 5p* 6 5d 2 0 4f 14 | 2.47 sigma .826";           // Hf
            case 73: return "5s 2 6s 2 5p* 6 5d* 3 0 4f 14 | 2.47 sigma .6";            // Ta
            case 74: return "5s 2 6s 2 5p 6 6p 2e-99 5d* 4 0 4f 14 | 2.32 sigma .62";   // W
            case 75: return "5s 2 6s 2 5p 6 6p 2e-99 5d* 5 0 4f 14 | 2.47 sigma .63";   // Re
            case 76: return "6s* 2 5p 6 6p 2e-99 5d* 5 1 | 2.35 sigma .58";             // Os
            case 77: return "6s* 2 5p 6 6p 2e-99 5d* 5 2 | 2.43 sigma .62";             // Ir
            case 78: return "6s* 1 0 5p 6 6p 2e-99 5d* 5 4 | 2.47 sigma .59";           // Pt
            case 79: return "6s* 1 0 6p* 2e-99 5d* 10 | 2.5 sigma .667";                // Au
            case 80: return "6s* 2 5p 6 6p 2e-99 5d* 10 | 2.44 sigma .59";              // Hg
            case 81: return "6s* 2 6p* 1 0 5d 10 5f | 2.25 sigma .718";                 // Tl
            case 82: return "6s* 2 6p* 2 0 5d* 10 | 2.3 sigma .59";                     // Pb
            case 83: return "6s* 2 6p* 3 0 5d* 10 | 2.41 sigma .605";                   // Bi
            case 84: return "6s* 2 6p* 3 1 5d* 10 | 2.3 sigma .54";                     // Po
            case 85: return "6s* 2 6p* 3 2 5d* 10 | 2.3 sigma .54";                     // At
            case 86: return "6s* 2 6p* 6 6d | 2.29 sigma .54";                          // Rn
#ifdef EXPERIMENTAL
            case 87: return "6s 2 7s 1 0 6p 6 7p 2e-99 | 2. sigma .5";                  // Fr
            case 88: return "7s 2 7p | 2. sigma .5";                                    // Ra
            case 89: return "7s 2 7p 6d 5f 1 0 | 2. sigma .5";                          // Ac
            case 90: return "7s 2 7p 6d 5f 2 0 | 2. sigma .5";                          // Th
            case 91: return "7s 2 7p 6d 5f 3 0 | 2. sigma .5";                          // Pa
            case 92: return "7s 2 7p 6d 5f 4 0 | 2. sigma .5";                          // U
            case 93: return "7s 2 7p 6d 5f 5 0 | 2. sigma .5";                          // Np
            case 94: return "7s 2 7p 6d 5f 6 0 | 2. sigma .5";                          // Pu
            case 95: return "7s 2 7p 6d 5f 7 0 | 2. sigma .5";                          // Am
//          case 96: return "6s 2 7s 2 6p 6 7p 2e-99 6d 5f 7 1 | 2. sigma .6";          // Cm
            case 96: return "6s 2 7s 2 6p* 6 6d 5f 7 1 | 2. sigma .6";                  // Cm
            case 97: return "7s 2 7p 6d 5f 7 2 | 2. sigma .5";                          // Bk
            case 98: return "6s 2 7s 2 6p 6 7p 2e-99 6d 2e-99 5f 7 3 | 2. sigma .5";    // Cf
            case 99: return "7s 2 7p 6d 5f 7 4 | 2. sigma .5";                          // Es
            case 100: return "7s 2 7p 6d 5f 7 5 | 2. sigma .5";                         // Fm
            case 101: return "7s 2 7p 6d 5f 7 6 | 2. sigma .5";                         // Md
            case 102: return "7s 2 7p 6d 5f 14 | 2. sigma .5";                          // No
            case 103: return "7s 2 7p 6d 1 0 | 2. sigma .5";                            // Lr
            case 104: return "7s 2 7p 6d 2 0 | 2. sigma .5";                            // Rf
            case 105: return "7s 2 7p 6d 3 0 | 2. sigma .5";                            // Db
            case 106: return "7s 2 7p 6d 4 0 | 2. sigma .5";                            // Sg
            case 107: return "7s 2 7p 6d 5 0 | 2. sigma .5";                            // Bh
            case 108: return "7s 2 7p 6d 5 1 | 2. sigma .5";                            // Hs
            case 109: return "7s 2 7p 6d 5 2 | 2. sigma .5";                            // Mt
            case 110: return "7s 2 7p 6d 5 3 | 2. sigma .5";                            // Ds
            case 111: return "7s 2 7p 6d 5 4 | 2. sigma .5";                            // Rg
            case 112: return "7s 2 7p 6d 10 | 2. sigma .5";                             // Cn
            case 113: return "7s 2 7p 1 0 7d | 2. sigma .5";                            // Nh
            case 114: return "7s 2 7p 2 0 7d | 2. sigma .5";                            // Fl
            case 115: return "7s 2 7p 3 0 7d | 2. sigma .5";                            // Mc
            case 116: return "7s 2 7p 3 1 7d | 2. sigma .5";                            // Lv
            case 117: return "7s 2 7p 3 2 7d | 3. sigma .9";                            // Ts
            case 118: return "7s 2 7p 6 7d | 3. sigma .8";                              // Og
            case 119: return "7s 2 8s 1 0 7p 6 6d 10 6f | 3. sigma .9";                 // ue
            case 120: return "7s 2 8s 2 7p 6 6d 10 6f | 3. sigma .9";                   // u0
            case 121: /* ================================== */                          // u1
            case 122: /*  Elements u1 (Z=121) through       */                          // u2
            case 123: /*           u6 (Z=126) are reserved  */                          // u3
            case 124: /*           for user configuration   */                          // u4
            case 125: /* ================================== */                          // u5
            case 126: return "warning_not_configured";                                  // u6
            case 127: return "1s -1 | 1.0 sigma .25 Z= -1";                             // e
            case   0: return "1s | 1.0 sigma .5";                                       // __
#endif // EXPERIMENTAL
            default:  warn("no default element configuration given for Z=%d", iZ);
                      return "";
        } // switch
  } // default config

//
// Some elements differ from the automatic choice of occupation numbers: Z= 24, 29, 41--47, 57, 78, 79
//

    int8_t constexpr KeyIgnore = 0, KeyRcut = -1, KeySigma = -2, KeyZcore = -3,
                     KeyMethod = -4, KeyWarn = -5, KeyNumax = -6, KeyNumeric = -7, KeyUndef = -8;
    char const Key2String[][8] = {"ignored", "|", "sigma", "Z=", "V", "warn", "numax", "numeric", "?"};

    inline int8_t char2ell(char const c) {
        switch (c) {
            case 's': return 0;
            case 'p': return 1;
            case 'd': return 2;
            case 'f': return 3;
            case 'g': return 4;
            case 'h': return 5;
            case 'i': return 6; // ToDo: according to chemistry standards, character 'i' should not be used!
            case 'j': return 7;
            case 'k': return 8;
            default : return -1; // invalid ell-character
        } // switch
    } // char2ell

    inline int8_t char2key(char const c) {
        // interpret the leading character of a word
        switch (c) { // case sensitive
            case ' ': case '\0': case '\t': case '\n': return KeyIgnore; // ignore whitespaces
            case 'r': case 'R': case '|': return KeyRcut;
            case 's': case 'S': return KeySigma;
            case 'Z': case 'z': return KeyZcore;
            case 'N': case 'n': return KeyNumax;
            case 'V': case 'v': return KeyMethod;
            case 'W': case 'w': return KeyWarn;
            case '0': case '.': case '+': case '-': return KeyNumeric; // numeric reading
            case '1': case '2': case '3': case '4': case '5':
            case '6': case '7': case '8': case '9': return c - '0'; // enn quantum number of an orbital
            default : return KeyUndef; // other characters, e.g. ':'
        } // switch
    } // char2key

    // same definition as atom_core::nl_index
    inline int nl_index(int const enn, int const ell) { 
        assert(ell < enn);
        return (enn*(enn - 1))/2 + ell;
    } // nl_index

    void set_default_core_shells(int ncmx[4], double const Z) {
        ncmx[0] = 0 + (Z >= 120) + (Z >= 88) + (Z >= 56) + (Z >= 38) + (Z >= 20) + (Z >= 12) + (Z >= 4) + (Z >= 2);
        ncmx[1] = 1 + (Z >= 118) + (Z >= 86) + (Z >= 54) + (Z >= 36) + (Z >= 18) + (Z >= 10);
        ncmx[2] = 2 + (Z >= 112) + (Z >= 80) + (Z >= 48) + (Z >= 30);
        ncmx[3] = 3 + (Z >= 102) + (Z >= 70);
    } // set_default_core_shells

    struct parsed_word_t {
        double value;
        int8_t key;
        int8_t enn;
        int8_t ell;
        int8_t mrn;
    }; // parsed_word_t

    
    element_t const & get(
          double const Zcore // nuclear charge
        , int const echo // =0 log-level
        , char const **configuration // =nullptr string that has been parsed
    ) {
        char symbol[4];
        int const iZ = chemical_symbol::get(symbol, Zcore);
        char element_Sy[16];
        std::snprintf(element_Sy, 15, "element_%s", symbol);
        auto const default_config_iZ = default_config(iZ);
        auto const config = control::get(element_Sy, default_config_iZ);
        if (configuration) *configuration = config;

        if (echo > 3) std::printf("# for Z=%g use %sconfiguration +%s=\"%s\"\n",
            Zcore, std::strcmp(default_config_iZ, config)?"":"default ", element_Sy, config);
        // now convert config into an element_t
        auto & e = *(new element_t);

        e.Z = ((iZ + 1) & 127) - 1; // preliminary integer number of protons
        e.rcut = 2.;
        e.sigma = .5;
        e.numax = -1; // automatic
        set(e.method, 16, '\0');
        set(e.nn, 8, uint8_t(0));

        if (nullptr == config) {
            warn("occupation numbers for Z=%g not set", e.Z);
            return e;
        }

        // start to parse the config string
        std::vector<parsed_word_t> words;
        words.reserve(32);
        int iword{0};
        char local_potential_method[32];

        char const * string{config + ('"' == config[0])}; // drop first char if it is '"'
        char c0{*string};
        while(c0) {
            if (echo > 11) std::printf("# start from '%s'\n", string);

            words.push_back(parsed_word_t());
            auto & w = words[iword];
            w.value = 0; // init neutral
            w.enn =  0; // init invalid
            w.ell = -1; // init invalid
            w.mrn = -3; // init invalid

            bool try_numeric{false};
            char const cn = *string;
            w.key = char2key(cn);
            if (echo > 10) std::printf("# found key=%i (%s) in '%s'\n", w.key, (w.key > 0)?"enn":Key2String[-w.key], string);
            if (w.key > 0) {
                // the first non-blank char was a digit larger than 0
                w.enn = w.key;
                if (w.enn > 9) error("enn=%i > 9 not supported in '%s'!", w.enn, string);
                char const cl = *(string + 1); // expect an ellchar here
                w.ell = char2ell(cl); // should we use (cl | 32) to allow also upper case ell-characters SPDFGH?
                if (w.ell >= 0) {
                    if (w.ell >= w.enn) error("unphysical ell=%i >= enn=%i in '%s'!", w.ell, w.enn, string);

                    char const *asterisk{string + 2};
                    w.mrn = 0; while ('*' == *asterisk) { ++w.mrn; ++asterisk; } // count the '*' chars

                    char const cm = *(string + 2);
                    if (0 == w.mrn) { // has no asterisks, check for 's'/'S':semicore or 'c'/'C':core specification
                        if ('s' == (cm | 32)) w.mrn = -1; // semicore state
                        if ('c' == (cm | 32)) w.mrn = -2; // core state (used to modify core state occupations)
                        if ('h' == (cm | 32)) error("config string might have \'hole\'-syntax: \"%s\"", config);
                    }
                    if (echo > 9) std::printf("# found enn=%i ell=%i mrn=%i in '%c%c%c'\n", w.enn, w.ell, w.mrn, cn,cl,cm);

                } else { // ell is a valid angular momentum quantum number
                    try_numeric = true;
                }
            } else {
                // leading character is not in {'1', ..., '9'}, i.e. a valid principal quantum number
                if (KeyNumeric == w.key) {
                    try_numeric = true;
                } else if (KeyMethod == w.key) {
                    std::strncpy(local_potential_method, string, 31);
                    if (echo > 7) std::printf("# found local potential method '%s'\n", local_potential_method);
                } else {
                    if (echo > 8) std::printf("# found special '%s'\n", string);
                }
            }

            if (try_numeric) {
                w.value = std::atof(string);
                if (echo > 8) std::printf("# found numeric value %g in '%s'\n", w.value, string);
                w.key = KeyNumeric; // -9:numeric
            } // try_numeric

            ++iword;
 
            auto const next_blank = std::strchr(string, ' '); // forward to the next w, ToDo: how does it react to \t?
            if (next_blank) {
                string = next_blank;
                while(*string == ' ') ++string; // forward to the next non-blank
                c0 = *string;
            } else {
                c0 = 0; // stop the while loop
            }
        } // while;
        int const nwords = iword; // how many words were in the string
        if (echo > 8) std::printf("# process %d words\n", nwords);

        if (echo > 7) {
            // repeat what was just parsed
            std::printf("# repeat config string '");
            char const mrn2string[][4] = {"", "*", "**", "***"};
            for (int iword = 0; iword < nwords; ++iword) {
                auto const & w = words[iword];
                if (w.key > 0) {
                    assert(w.enn == w.key && w.ell >= 0 && w.ell < w.enn && w.mrn > -3); // orbital
                    std::printf("%i%c%s ", w.enn, ellchar[w.ell], (w.mrn >= 0) ? mrn2string[w.mrn] : csv_name(w.mrn + 2));
                } else if (KeyNumeric == w.key) {
                    std::printf("%g ", w.value);
                } else if (KeyIgnore == w.key) {
                    // do not print anything
                } else {
                    assert(w.key <= 0);
                    std::printf("%s ", Key2String[-w.key]);
                }
            } // iword
            std::printf("'\n");
        } // echo

        int constexpr max_enn = 9, max_inl = (max_enn*(max_enn + 1))/2;
        double occ[max_inl][2]; // spin resolved occupation numbers
        set(occ[0], max_inl*2, 0.0); // clear
        std::vector<int8_t> csv(max_inl, csv_undefined);

        int ncmx[8]; set(ncmx, 8, 0);
        set_default_core_shells(ncmx, e.Z);      

        int constexpr max_nstack = 4;
        double stack[max_nstack] = {0, 0, 0, 0};
        int nstack{0};

        // process the parsed contents backwards
        for (int iword = nwords - 1; iword >= 0; --iword) {
            auto const & w = words[iword];
            if (KeyNumeric == w.key) {
                if (echo > 21) std::printf("# nstack=%i push\n", nstack);
                assert(nstack < max_nstack);
                stack[nstack++] = w.value; // push to the stack
                if (echo > 21) std::printf("# nstack=%i stack[%i]=%g\n", nstack, nstack - 1, stack[nstack - 1]);
            } else if (KeyMethod == w.key) {
//              warn("method specifier in config string for %s ignored : %s", symbol, config);
                std::strncpy(e.method, local_potential_method + 2, 15); // +2 to strip the "v=" or "V="
                if (echo > 9) std::printf("# found local potential method = \'%s\'\n", e.method);
            } else if (KeyWarn == w.key) {
                warn("config string for %s may be experimental: %s", symbol, config);
            } else if (KeyIgnore == w.key) {
                // do not modify the stack!
                if (echo > 6) std::printf("# white spaces are ignored for Z= %g\n", e.Z);
            } else {
                double value{0};
                if (nstack > 0) {
                    if (echo > 21) std::printf("# nstack=%i stack[%i]=%g\n", nstack, nstack - 1, stack[nstack - 1]);
                    value = stack[--nstack]; // pop the stack
                    if (echo > 21) std::printf("# nstack=%i popped\n", nstack);
                }
                if (w.key > 0) { // orbital

                    // read one or two occupation numbers from the stack
                    double occs[2] = {0, 0};
                    if (nstack > 0) {
                        occs[0] = value;
                        if (echo > 21) std::printf("# nstack=%i stack[%i]=%g\n", nstack, nstack - 1, stack[nstack - 1]);
                        value = stack[--nstack]; // pop the stack a 2nd time
                        if (echo > 21) std::printf("# nstack=%i popped\n", nstack);
                        occs[1] = value;
                    } else {
                        set(occs, 2, .5*value); // split the value between dn-spin and up-spin
                    }
                    // checks for the boundaries of occupation numbers
                    auto const max_occ = 2*w.ell + 1;
                    for (int spin = 0; spin < 2; ++spin) {
                        if (occs[spin] < 0 && e.Z > 0) {
                            warn("found a negative occupation number %g in the %s-%i%c orbital", occs[spin], symbol, w.enn, ellchar[w.ell]);
                            occs[spin] = 0;
                        }
                        if (occs[spin] > max_occ) {
                            warn("occupation number %g is too large for a %s-%i%c orbital", occs[spin], symbol, w.enn, ellchar[w.ell]);
                            occs[spin] = max_occ;
                        }
                    } // spin

                    int const inl = nl_index(w.enn, w.ell);
                    csv[inl] = std::min(valence, w.mrn + 2); // classify: mrn == -2:core, -1:semicore, 0:valence, 1:valence, ...

                    if (echo > 9) std::printf("# found orbital %i%c occ= %g %g inl=%i\n", w.enn,ellchar[w.ell], occs[0],occs[1], inl);
                    assert(w.key > 0); // this is a valence orbital specification
                    set(occ[inl], 2, occs); // set valence occupation numbers
                    if ((valence == csv[inl]) && (w.ell < 8)) {
                        e.nn[w.ell] += 1 + w.mrn; // number of partial waves
                        ncmx[w.ell] = w.enn - 1; // highest enn of a core state
                    }

                } else if (KeyRcut == w.key) {
                    e.rcut = value;
                    if (echo > 9) std::printf("# found cutoff radius rcut= %g\n", e.rcut);
                    if (e.rcut <= 0) warn("rcut must be positive but found rcut=%g", e.rcut);
                } else if (KeySigma == w.key) {
                    e.sigma = value;
                    if (echo > 9) std::printf("# found projector spread sigma= %g\n", e.sigma);
                    if (e.sigma <= 0) warn("sigma must be positive but found sigma=%g", e.sigma);
                } else if (KeyNumax == w.key) {
                    e.numax = int(value);
                    if (echo > 9) std::printf("# found SHO projector cutoff numax= %d\n", e.numax);
                    if (std::abs(e.numax - value) > 1e-6) warn("numax must be a positive integer found %g --> %d", value, e.numax);
                } else if (KeyZcore == w.key) {
                    e.Z = value;
                    set_default_core_shells(ncmx, e.Z); // adjust default core shells
                    if (echo > 9) std::printf("# found core charge Z= %g for %s\n", e.Z, symbol);
                    if (e.Z >= 120) warn("some routine may not be prepared for Z= %g >= 120", e.Z);
                } else {
                    warn("%s unknown key= %d (%s) in word #%i", symbol, w.key, Key2String[-w.key], iword);
                }
            }
        } // iword
        if (nstack > 0) warn("after parsing, some value %g was left on the stack", stack[0]);

        // fill the core
        if (echo > 6) {
            std::printf("# fill the core up to the principal quantum number n=");
            for (int ell = 0; ell < 4; ++ell) {
                std::printf(" %d", ncmx[ell]*(ncmx[ell] > ell)); // show zeros if there are no core electrons
            } // ell
            std::printf(" for s,p,d,f\n");
        } // echo

        // set all core states
        for (int ell = 0; ell < 4; ++ell) {
            double const max_occ = 2*ell + 1;
            for (int enn = ell + 1; enn <= ncmx[ell]; ++enn) {
                int const inl = nl_index(enn, ell);
                if (csv_undefined == csv[inl]) {
                    csv[inl] = core; // classify as core
                    for (int spin = 0; spin < 2; ++spin) {
                        occ[inl][spin] = max_occ; // occupy fully
                    } // spin
                } else if (valence != csv[inl]) {
                    // inform about not fully occupied core or semicore states, "core holes"
                    if (occ[inl][0] + occ[inl][1] < 2*max_occ) {
                        if (echo > 0) std::printf("# %s-%d%c %s orbital is only occupied with %g + %g electrons\n",
                                           symbol, enn,ellchar[ell], csv_name(csv[inl]), occ[inl][0], occ[inl][1]);
                    }
                }
            } // enn
        } // ell

        double ne[4] = {0, 0, 0, 0}; // number of valence, semicore and core electrons, others?
        for (int inl = 0; inl < max_inl; ++inl) {
            if (inl < 32) {
                ne[csv[inl]] += occ[inl][0] + occ[inl][1];
                set(e.occ[inl], 2, occ[inl]); // copy occupation numbers
                e.csv[inl] = csv[inl];
            } else {
                ne[csv_undefined] += occ[inl][0] + occ[inl][1]; // electrons lost
            } // inl < 32
        } // inl

        if (echo > 11) { // show all occupation numbers
            std::printf("# Z=%3i occ=", iZ);
            for (int inl = 0; inl < 30; ++inl) {
                std::printf("%3i", int(std::round(std::abs(occ[inl][0]) + std::abs(occ[inl][1]))));
            } // inl
            std::printf("\n");
        } // echo

        double const nelectrons = ne[core] + ne[semicore] + ne[valence];
        if (echo > 4) std::printf("# found %g electrons = %g core + %g semicore + %g valence electrons for Z= %g protons\n",
                                           nelectrons, ne[core], ne[semicore], ne[valence], e.Z);
        if (ne[csv_undefined] != 0) warn("lost %g electrons for Z= %g", ne[csv_undefined], e.Z);

        if (echo > 4) {
            std::printf("# PAW setup for %s (Z=%g) suggests", symbol, e.Z);
            for (int ell = 0; ell < 8; ++ell) {
                std::printf(" %d", e.nn[ell]);
            } // ell
            std::printf(" partial waves for s,p,d,...\n");
        } // echo

        if ((std::abs(nelectrons - e.Z) > 1e-12) && (e.Z > 0)) {
            warn("PAW setup for %s (Z=%g) is charged with %g electrons", symbol, e.Z, nelectrons - e.Z);
        }

        return e;
    } // get

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_parsing(int const echo=0) {

      int const iZ_show = control::get("sigma_config.test.Z", -120); // use -127 to check also the custom configuration strings and vacuum
      if (iZ_show >= 0) {

          // parse one selected element and show the configuration string used
          if (echo > 0) std::printf("\n# running with +sigma_config.test.Z=%i\n", iZ_show);
          char Sy[4];
          auto const iZ = chemical_symbol::get(Sy, iZ_show);
          char const *actual_config{nullptr};
          auto const e = get(iZ, echo, &actual_config);
          if (echo > 0) std::printf("\n+element_%s=\"%s\"\n\n", Sy, actual_config);
          if (echo > 0) std::printf("# found Z= %g\n", e.Z);

      } else { // iZ_show >= 0

          int const Z_max = std::min(127, -iZ_show);
          // parse all elements in [1, Z_max]
#ifdef EXPERIMENTAL
          if (echo > 8) std::printf("\n\n# sizeof(element_t) = %ld Byte\n", sizeof(element_t));

          if (echo > 2) std::printf("\n\n# parse EXPERIMENTAL elements 58--70, 87--%d\n\n", Z_max);
          for (int iZ = 58; iZ <= Z_max; ++iZ) {
              if (std::abs(iZ - 64) < 7 || iZ > 86) {
                  if (echo > 4) std::printf("\n");
                  auto const e = get(iZ, echo);
                  if (echo > 4) std::printf("# Z=%g rcut=%g sigma=%g Bohr\n", e.Z, e.rcut, e.sigma);
              } // with 58 through 70 and more
          } // iZ
#endif // EXPERIMENTAL

          if (echo > 2) std::printf("\n\n# parse configuration strings for elements 86--71, 57--1\n\n");
          for (int iZ = std::min(86, Z_max); iZ > 0; --iZ) {
              if (std::abs(iZ - 64) >= 7) {
                  if (echo > 4) std::printf("\n");
                  auto const e = get(iZ, echo);
                  if (echo > 4) std::printf("# Z=%g rcut=%g sigma=%g Bohr\n", e.Z, e.rcut, e.sigma);
              } // without 58 through 70
          } // iZ

      } // iZ_show >= 0
      return 0; // success (all parsing errors are fatal)
  } // test_parsing
  
  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_parsing(echo);
      return stat;
  } // all_tests

#endif

} // namespace sigma_config
