#include <cstdio> // printf, std::snprintf
#include <cmath> // std::floor
#include <cstdint> // int8_t
#include <vector> // std::vector<T>
#include <cassert> // assert

#include "sigma_config.hxx"

#include "inline_math.hxx" // set
#include "status.hxx" // status_t
#include "chemical_symbol.hxx" // ::get
#include "control.hxx" // ::get
#include "recorded_warnings.hxx" // warn, error

namespace sigma_config {

    char const ellchar[12] = "spdfghijkl?";

#define EXPERIMENTAL

    char const * default_config(unsigned const iZ) {

        switch (iZ) {
            case  1: return "1s* 1 0 2p | 0.9 sigma .247";                              // H   
            case  2: return "1s* 2 2p | 1.5 sigma .48";                                 // He  
            case  3: return "2s* 1 0 2p 2e-99 | 2.0 sigma .6";                          // Li  
            case  4: return "2s* 2 2p 2e-99 0 | 1.5 sigma .45";                         // Be  
            case  5: return "2s* 2 2p* 1 0 3d | 1.2 sigma .45";                         // B   
            case  6: return "2s* 2 2p* 2 0 3d | 1.2 sigma .43";                         // C   
            case  7: return "2s* 2 2p* 3 0 3d | 1.0 sigma .33";                         // N   
            case  8: return "2s* 2 2p* 3 1 3d | 1.13 sigma .297";                       // O   
            case  9: return "2s* 2 2p* 3 2 3d | 1.2 sigma .323";                        // F   
            case 10: return "2s* 2 2p* 6 3d | 1.8 sigma .564";                          // Ne  
            case 11: return "3s* 1 0 2p 6 3p 2e-99 3d | 2.27 sigma .69";                // Na
            case 12: return "3s* 2 3p 2e-99 3d | 1.96 sigma .41";                       // Mg (2 e)
//          case 12: return "2s 2 3s 2 2p 6 3p 2e-99 3d | 1.96 sigma .41";              // Mg (10 e)
            case 13: return "3s* 2 3p* 1 0 3d | 2.05 sigma .645";                       // Al  
            case 14: return "3s* 2 3p* 2 0 3d | 2.0 sigma .643";                        // Si  
            case 15: return "3s* 2 3p* 3 0 3d | 1.8 sigma .512";                        // P   
            case 16: return "3s* 2 3p* 3 1 3d | 1.7 sigma .535";                        // S   
            case 17: return "3s* 2 3p* 3 2 3d | 1.5 sigma .503";                        // Cl  
            case 18: return "3s* 2 3p* 6 3d | 1.6 sigma .546";                          // Ar  
            case 19: return "3s 2 4s 1 0 3p* 6 3d | 1.77 sigma .47";                    // K   
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
            case 30: return "4s* 2 4p* 2e-99 3d* 10 | 2.23 sigma .577";                 // Zn  
            case 31: return "4s* 2 4p* 1 0 4d | 2.2 sigma .686";                        // Ga  
            case 32: return "4s* 2 4p* 2 0 4d | 1.9 sigma .606";                        // Ge  
            case 33: return "4s* 2 4p* 3 0 4d | 2.0 sigma .62";                         // As  
            case 34: return "4s* 2 4p* 3 1 4d | 1.6 sigma .521";                        // Se  
            case 35: return "4s* 2 4p* 3 2 4d | 2.1 sigma .6";                          // Br  
            case 36: return "4s* 2 4p* 6 4d | 2.2 sigma .61";                           // Kr 
            case 37: return "4s 2 5s 1 0 4p* 6 4d | 2.3 sigma .78";                     // Rb 
            case 38: return "4s 2 5s 2 4p* 6 4d | 2.37 sigma .666";                     // Sr
            case 39: return "4s 2 5s 2 4p 6 5p 2e-99 4d* 1 0 | 2.43 sigma .6";          // Y 
            case 40: return "4s 2 5s 2 4p 6 5p 2e-99 4d* 2 0 | 2.35 sigma .58";         // Zr
            case 41: return "4s 2 5s 1 0 4p 6 5p 2e-99 4d* 4 0 | 2.35 sigma .59";       // Nb
            case 42: return "4s 2 5s 1 0 4p 6 5p 2e-99 4d* 5 0 | 2.34 sigma .585";      // Mo
            case 43: return "4s 2 5s 1 0 4p 6 5p 2e-99 4d* 5 1 | 2.4 sigma .58";        // Tc
            case 44: return "4s 2 5s 1 0 4p 6 5p 2e-99 4d* 5 2 | 2.37 sigma .571";      // Ru
            case 45: return "5s* 1 0 4p 6 5p 2e-99 4d* 5 3 | 2.35 sigma .58";           // Rh
            case 46: return "5s* 2e-99 0 4p 6 5p 2e-99 4d* 10 | 2.32 sigma .585";       // Pd
            case 47: return "5s* 1 0 4p 6 5p 2e-99 4d* 10 | 2.23 sigma .57";            // Ag
            case 48: return "5s* 2 5p* 2e-99 4d* 10 | 2.2 sigma .563";                  // Cd 
            case 49: return "5s* 2 5p* 1 0 4d* 10 | 2.17 sigma .565";                   // In 
            case 50: return "5s* 2 5p* 2 0 4d* 10 | 2.24 sigma .585";                   // Sn
            case 51: return "5s* 2 5p* 3 0 4d* 10 | 2.18 sigma .57";                    // Sb
            case 52: return "5s* 2 5p* 3 1 5d | 2.23 sigma .555";                       // Te
            case 53: return "5s* 2 5p* 3 2 5d | 2.2 sigma .68";                         // I 
            case 54: return "5s* 2 5p* 6 5d | 2.24 sigma .62";                          // Xe 
            case 55: return "5s 2 6s 1 0 5p* 6 5d | 2.0 sigma .61";                     // Cs 
            case 56: return "5s 2 6s 2 5p* 6 5d* | 2.2 sigma .645";                     // Ba  
            case 57: return "5s 2 6s 2 5p* 6 5d* 1 0 | 1.9 sigma .59";                  // La
#ifdef EXPERIMENTAL
            case 58: return "5s 2 6s 2 5p* 6 5d* 0 4f 2 0 | 2. sigma .6";               // Ce
            case 59: return "5s 2 6s 2 5p* 6 5d* 0 4f 3 0 | 2. sigma .6";               // Pr
            case 60: return "5s 2 6s 2 5p* 6 5d* 0 4f 4 0 | 2. sigma .6";               // Nd
            case 61: return "5s 2 6s 2 5p* 6 5d* 0 4f 5 0 | 2. sigma .6";               // Pm
            case 62: return "5s 2 6s 2 5p* 6 5d* 0 4f 6 0 | 2. sigma .6";               // Sm
            case 63: return "5s 2 6s 2 5p* 6 5d* 0 4f 7 0 | 2. sigma .6";               // Eu
            case 64: return "5s 2 6s 2 5p* 6 5d* 0 4f 7 1 | 2. sigma .6";               // Gd
            case 65: return "5s 2 6s 2 5p* 6 5d* 0 4f 7 2 | 2. sigma .6";               // Tb
            case 66: return "5s 2 6s 2 5p* 6 5d* 0 4f 7 3 | 2. sigma .6";               // Dy
            case 67: return "5s 2 6s 2 5p* 6 5d* 0 4f 7 4 | 2. sigma .6";               // Ho
            case 68: return "5s 2 6s 2 5p* 6 5d* 0 4f 7 5 | 2. sigma .6";               // Er
            case 69: return "5s 2 6s 2 5p* 6 5d* 0 4f 7 6 | 2. sigma .6";               // Tm
            case 70: return "5s 2 6s 2 5p* 6 5d* 0 4f 7 7 | 2. sigma .6";               // Yb
#endif
            case 71: return "5s 2 6s 2 5p 6 6p 2e-99 5d* 1 0 | 2.4 sigma .6";           // Lu  
            case 72: return "5s 2 6s 2 5p 6 6p 2e-99 5d* 2 0 | 2.47 sigma .6077";       // Hf  
            case 73: return "5s 2 6s 2 5p 6 6p 2e-99 5d* 3 0 | 2.47 sigma .6";          // Ta  
            case 74: return "5s 2 6s 2 5p 6 6p 2e-99 5d* 4 0 | 2.32 sigma .62";         // W   
            case 75: return "6s* 2 5p 6 6p 2e-99 5d* 5 0 | 2.47 sigma .63";             // Re  
            case 76: return "6s* 2 5p 6 6p 2e-99 5d* 5 1 | 2.35 sigma .58";             // Os  
            case 77: return "6s* 2 5p 6 6p 2e-99 5d* 5 2 | 2.43 sigma .62";             // Ir  
            case 78: return "6s* 1 0 5p 6 6p 2e-99 5d* 5 4 | 2.47 sigma .59";           // Pt  
            case 79: return "6s* 1 0 6p* 2e-99 5d* 10 | 2.5 sigma .667";                // Au  
            case 80: return "6s* 2 5p 6 6p 2e-99 5d* 10 | 2.44 sigma .59";              // Hg  
            case 81: return "6s* 2 6p* 1 0 5d* 10 | 2.25 sigma .567";                   // Tl  
            case 82: return "6s* 2 6p* 2 0 5d* 10 | 2.3 sigma .59";                     // Pb  
            case 83: return "6s* 2 6p* 3 0 5d* 10 | 2.41 sigma .605";                   // Bi  
            case 84: return "6s* 2 6p* 3 1 5d* 10 | 2.3 sigma .54";                     // Po  
            case 85: return "6s* 2 6p* 3 2 5d* 10 | 2.3 sigma .54";                     // At  
            case 86: return "6s* 2 6p* 6 6d | 2.29 sigma .54";                          // Rn  
#ifdef EXPERIMENTAL
            case 87: return "7s 1 0 7p | 2. sigma .5";                                  // Fr
            case 88: return "7s 2 7p | 2. sigma .5";                                    // Ra
            case 89: return "7s 2 7p 6d 5f 1 0 | 2. sigma .5";                          // Ac
            case 90: return "7s 2 7p 6d 5f 2 0 | 2. sigma .5";                          // Th
            case 91: return "7s 2 7p 6d 5f 3 0 | 2. sigma .5";                          // Pa
            case 92: return "7s 2 7p 6d 5f 4 0 | 2. sigma .5";                          // U
            case 93: return "7s 2 7p 6d 5f 5 0 | 2. sigma .5";                          // Np
            case 94: return "7s 2 7p 6d 5f 6 0 | 2. sigma .5";                          // Pu
            case 95: return "7s 2 7p 6d 5f 7 0 | 2. sigma .5";                          // Am
            case 96: return "7s 2 7p 6d 5f 7 1 | 2. sigma .5";                          // Cm
            case 97: return "7s 2 7p 6d 5f 7 2 | 2. sigma .5";                          // Bk
            case 98: return "7s 2 7p 6d 5f 7 3 | 2. sigma .5";                          // Cf
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
            case 118: return "7s 2 7p 6 7d | 3. sigma .9";                              // Og
            case 119: return "7s 2 8s 1 0 7p 6 6d 10 6f | 3. sigma .9";                 // ue
            case 120: return "7s 2 8s 2 7p 6 6d 10 6f | 3. sigma .9";                   // u0
            case 121: /* ================================== */                          // u1
            case 122: /*  Elements u1 (Z=121) through       */                          // u2
            case 123: /*           u6 (Z=126) are reserved  */                          // u3
            case 124: /*           for user configuration   */                          // u4
            case 125: /* ================================== */                          // u5
            case 126: return "warning_not_configured";                                  // u6
            case 127: return "1s -1 | 1.0 sigma .25 Z= -1";                             // e
            case   0: return "1s | 1.0 sigma .25";                                      // __
#endif
            default:
                     warn("no default element configuration given for Z=%d", iZ);
                     return "";
        } // switch
    } // default config

    
    int8_t constexpr KeyIgnore = 0, KeyRcut = -1, KeySigma = -2, KeyZcore = -3,
      KeyMethod = -4, KeyHole = -5, KeyWarn = -6, KeyUndef = -8, KeyNumeric = -9;
    char constexpr Key2Char[] = "_|sZVhW?Un"; // negative keys needed: [-key]
    
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
        switch (c) { // case sensitive
            case ' ': case '\0': case '\t': case '\n': return KeyIgnore;
            case 'r': case 'R': case '|': return KeyRcut;
            case 's': case 'S': return KeySigma;
            case 'Z': case 'z': return KeyZcore;
            case 'V': case 'v': return KeyMethod;
            case 'W': case 'w': return KeyWarn;
            case '0': case '.': case '-': return KeyNumeric; // numeric reading
            default : return c - '0'; // enn quantum number of an orbital
        } // switch
    } // char2key
    

    element_t& get(double const Zcore, int const echo) {
        
        char symbol[4], element_Sy[16];
        int const iZ = chemical_symbol::get(symbol, Zcore);
        std::snprintf(element_Sy, 15, "element_%s", symbol);
        auto const config = control::get(element_Sy, default_config(iZ));

        if (echo > 3) printf("# for Z=%g use configuration %s=\"%s\"\n", Zcore, element_Sy, config);
        // now convert config into an element_t
        auto e = new element_t;
        
        double const Z = ((iZ + 1) & 127) - 1;
        e->Z = Z;
        e->rcut = 2.;
        e->sigma = .5;
        e->q_core_hole[0] = 0;
        e->q_core_hole[1] = 0;
        e->inl_core_hole = -1;
        set(e->nn, 8, uint8_t(0));
        e->ncmx[0] = (Z >= 2) + (Z >= 4) + (Z >= 12) + (Z >= 20) + (Z >= 38) + (Z >= 56) + (Z >= 88) + (Z >= 120);
        e->ncmx[1] = 1 + (Z >= 10) + (Z >= 18) + (Z >= 36) + (Z >= 54) + (Z >= 86) + (Z >= 118);
        e->ncmx[2] = 2 + (Z >= 30) + (Z >= 48) + (Z >= 80) + (Z >= 112);
        e->ncmx[3] = 3 + (Z >= 70) + (Z >= 102);

        if (nullptr == config) return *e;

        int constexpr mwords = 32;
        double values[mwords];
        int8_t keys[mwords];
        int8_t enns[mwords];
        int8_t ells[mwords];
        int8_t mrns[mwords];

        int iword{0};
        
        char const * string{config};
        char c0{*string};
        while(c0) {
//          if (echo > 0) printf("# start from '%s'\n", string);
          
            assert(iword < mwords);
          
            values[iword] = 0.0;
            enns[iword] =  0; // undefined
            ells[iword] = -1; // undefined
            mrns[iword] = -1; // undefined
            keys[iword] = KeyUndef; //
          
            bool try_numeric{false};
            char const cn = *string;
            auto key = char2key(cn);
//                 if (echo > 0) printf("# found key=%i in '%s'\n", key, string);
            if (key > KeyIgnore) {
                int const enn = key;
                if (enn > 9) error("enn=%i > 9 not supported in '%s'!", enn, string);
                char const cl = *(string + 1);
                auto const ell = char2ell(cl);
                if (ell >= 0) {
                    if (ell >= enn) error("unphysical ell=%i >= enn=%i in '%s'!", ell, enn, string);

                    char const cm = *(string + 2);
                    if ((cm | 32) == 'h') key = KeyHole; // core hole, modify key
                    int const mrn = ('*' == cm); // max radial nodes
                    if (echo > 9) printf("# found enn=%i ell=%i mrn=%i in '%c%c%c'\n", enn, ell, mrn, cn,cl,cm);
                    enns[iword] = enn; // store
                    ells[iword] = ell; // store
                    mrns[iword] = mrn; // store

                } else { // ell is a valid angular momentum quantum number
                    try_numeric = true;
                }
            } else { // enn is valid principal quantum number
                try_numeric = (KeyNumeric == key);
                if (!try_numeric && echo > 8) printf("# found special '%s'\n", string);
            }
            keys[iword] = key; // store

            if (try_numeric) {
                double const value = std::atof(string);
                if (echo > 8) printf("# found numeric value %g in '%s'\n", value, string);
                values[iword] = value; // store
                keys[iword] = KeyNumeric; // -9:numeric
            } // try_numeric

            ++iword;
            
            auto const next_blank = std::strchr(string, ' '); // forward to the next word, ToDo: how does it react to \t?
            if (next_blank) {
                string = next_blank;
                while(*string == ' ') ++string; // forward to the next non-blank
                c0 = *string;
            } else {
                c0 = 0; // stop the while loop
            }
        } // while;
        int const nwords = iword; // how many words were in the string
        if (echo > 8) printf("# process %d words\n", nwords);
        
        if (echo > 8) {
            // repeat what was just parsed
            printf("# repeat config string '");
            for (int iword = 0; iword < nwords; ++iword) {
                int8_t const key = keys[iword];
                if (key > 0) { auto const enn = enns[iword]; assert(enn == key); // orbital
                    printf("%i%c%s ", enn, ellchar[ells[iword]], mrns[iword]?"*":"");
                } else if (KeyNumeric == key) {
                    printf("%g ", values[iword]);
                } else if (KeyIgnore == key) {
                    // do not print anything
                } else {
                    printf("%c ", Key2Char[-key]);
                }
            } // iword
            printf("'\n");
        } // echo

        int constexpr max_enn = 9, max_inl = (max_enn*(max_enn + 1))/2;
        double occ[max_inl][2];
        for (int inl = 0; inl < max_inl; ++inl) { occ[inl][0] = occ[inl][1] = 0; }
        
        double stack[4] = {0, 0, 0, 0};
        int nstack{0};

        // now put into context:
        for(int iword = nwords - 1; iword >= 0; --iword) {
            int8_t const key = keys[iword];
            if (KeyNumeric == key) {
//                 if (echo > 0) printf("# nstack=%i push\n", nstack);
                stack[nstack++] = values[iword]; // push to the stack
//                 if (echo > 0) printf("# nstack=%i stack[%i]=%g\n", nstack, nstack - 1, stack[nstack - 1]);
            } else if (KeyMethod == key) {
                warn("method specifier in config string for %s ignored : %s", symbol, config);
            } else if (KeyWarn == key) {
                warn("config string for %s may be experimental: %s", symbol, config);
            } else if (KeyIgnore == key) {
                // do not modify the stack!
                if (echo > 6) printf("# white spaces are ignored for Z = %g\n", e->Z);
            } else {
                double value{0};
                if (nstack > 0) {
//                     if (echo > 0) printf("# nstack=%i stack[%i]=%g\n", nstack, nstack - 1, stack[nstack - 1]);
                    value = stack[--nstack]; // pop the stack
//                     if (echo > 0) printf("# nstack=%i popped\n", nstack);
                }
                if (key > KeyIgnore || key == KeyHole) { // orbital
                    int const ell = ells[iword];
                    int const enn = enns[iword];
                    int const inl = ell + (enn*(enn-1))/2; assert(ell < enn);
                    double occs[2] = {0, 0};
                    if (nstack > 0) {
                        occs[0] = value;
//                         if (echo > 0) printf("# nstack=%i stack[%i]=%g\n", nstack, nstack - 1, stack[nstack - 1]);
                        value = stack[--nstack]; // pop the stack a 2nd time
//                         if (echo > 0) printf("# nstack=%i popped\n", nstack);
                        occs[1] = value;
                    } else {
                        set(occs, 2, .5*value);
                    }
                    // checks for the boundaries of occupation numbers
                    for (int spin = 0; spin < 2; ++spin) {
                        if (occs[spin] < 0 && e->Z > 0) {
                            warn("found a negative occupation number %g in the %s-%i%c orbital", occs[spin], symbol, enn, ellchar[ell]);
                            occs[spin] = 0;
                        }
                        if (occs[spin] > 2*ell + 1) {
                            warn("occupation number %g is too large for a %s-%i%c orbital", occs[spin], symbol, enn, ellchar[ell]);
                            occs[spin] = 2*ell + 1;
                        }
                    } // spin
                    if (echo > 9) printf("# found orbital %i%c occ= %g %g inl=%i\n", enn,ellchar[ell], occs[0], occs[1], inl);
                    if (key > KeyIgnore) { // orbital
                        set(occ[inl], 2, occs); // set valence occupation numbers positive
                        if (ell < 8) e->nn[ell] += 1 + mrns[iword];
                        if (ell < 4) e->ncmx[ell] = enn - 1;
                    } else if (key == KeyHole) {
                        set(e->q_core_hole, 2, occs);
                        e->inl_core_hole = inl;
                        if (echo > 2) printf("# found a %i%c-core hole of charges = %g %g in %s\n", enn,ellchar[ell], occs[0],occs[1], symbol);
                    } else assert(false); // should not occur
                } else if (KeyRcut == key) {
                    e->rcut = value;     
                    if (echo > 9) printf("# found cutoff radius rcut = %g\n", e->rcut);
                    if(e->rcut <= 0) warn("rcut must be positive but found rcut=%g", e->rcut);
                } else if (KeySigma == key) {  
                    e->sigma = value;
                    if (echo > 9) printf("# found projector spread sigma = %g\n", e->sigma);
                    if(e->sigma <= 0) warn("sigma must be positive but found sigma=%g", e->sigma);
                } else if (KeyZcore == key) {  
                    e->Z = value;
                    if (echo > 9) printf("# found core charge Z = %g for %s\n", e->Z, symbol);
                    if(e->Z >= 120) warn("some routine may not be prepared for Z = %g >= 120", e->Z);
                }
            }
        } // iword
        if (nstack > 0) warn("after parsing, some value %g was left on the stack", stack[0]);
        
        
        // fill the core
        if (echo > 6) {
            printf("# fill the core up to enn =");
            for(int ell = 0; ell < 4; ++ell) {
                printf(" %d", e->ncmx[ell]*(e->ncmx[ell] > ell)); // show zeros if there are no core electrons
            } // ell
            printf(" for s,p,d,f\n");
        } // echo
        for (int ell = 0; ell < 4; ++ell) {
            for(int enn = ell + 1; enn <= e->ncmx[ell]; ++enn) {
                int const inl = ell + (enn*(enn-1))/2; assert(ell < enn);
                for(int spin = 0; spin < 2; ++spin) {
                    if (occ[inl][spin] <= 0) occ[inl][spin] = -(2*ell + 1); // negative occupation numbers indicate core electrons
                } // spin
            } // enn
        } // ell

        { // scope: introduce a core hole
            int const inl = e->inl_core_hole;
            if (inl >= 0) {
                if (echo > 2) printf("# introduce a core hole in inl=%i with charge %g %g electrons\n", inl, e->q_core_hole[0], e->q_core_hole[1]);
                occ[inl][0] += e->q_core_hole[0]; // add because core occupation numbers are negative
                occ[inl][1] += e->q_core_hole[1];
            } // inl valid
        } // scope

        double nve{0}, nce{0}, ncv{0}; // number of valence and core electrons, others?
        for(int inl = 0; inl < max_inl; ++inl) {
            nve += std::max(0.0, occ[inl][0]) + std::max(0.0, occ[inl][1]);
            nce -= std::min(0.0, occ[inl][0]) + std::min(0.0, occ[inl][1]);
            if (inl < 32) {
                set(e->occ[inl], 2, occ[inl]);
            } else {
                ncv += std::abs(occ[inl][0]) + std::abs(occ[inl][1]);
            } // inl < 32
        } // inl
        double const nelectrons = nve + nce;
        if (echo > 4) printf("# found %g = %g core + %g valence electrons and Z = %g\n", nelectrons, nce, nve, e->Z);
        if (ncv) warn("lost %g electrons for Z = %g\n", ncv, e->Z);

        if (echo > 4) {
            printf("# PAW setup for %s (Z=%g) suggests", symbol, e->Z);
            for(int ell = 0; ell < 8; ++ell) {
                printf(" %d", e->nn[ell]);
            } // ell
            printf(" partial waves\n");
        } // echo
        
        if (std::abs(nelectrons - e->Z) > 1e-12 && e->Z > 0) warn("PAW setup for %s (Z=%g) is charged with %g electrons", symbol, e->Z, nelectrons - e->Z);
        
        return *e;
    } // get

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  inline status_t test_86(int const echo=0) {
      for (int iZ = 86; iZ > 0; --iZ) {
          if (std::abs(iZ - 64) >= 7) {
              if (echo > 4) printf("\n");
              auto const e = get(iZ, echo);
              if (echo > 0) printf("# Z=%g rcut=%g sigma=%g Bohr\n", e.Z, e.rcut, e.sigma);
          } // without 58 through 70
      } // iZ
#ifdef EXPERIMENTAL
      if (echo > 0) printf("\n\n# EXPERIMENTAL elements 58--70, 87--128\n\n");
      for (int iZ = 58; iZ <= 120; ++iZ) {
          if (std::abs(iZ - 64) < 7 || iZ > 86) {
              if (echo > 4) printf("\n");
              auto const e = get(iZ & 127, echo);
              if (echo > 0) printf("# Z=%g rcut=%g sigma=%g Bohr\n", e.Z, e.rcut, e.sigma);
          } // with 58 through 70 and more
      } // iZ
#endif
      return 0;
  } // test_86

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_86(echo);
      return stat;
  } // all_tests

#endif

} // namespace sigma_config
