#pragma once

#include <cstdio> // printf, std::snprintf
#include <cmath> // std::floor
#include <cstdint> // int8_t
#include <vector> // std::vector<T>

#include "chemical_symbol.hxx" // ::get
#include "control.hxx" // ::get
#include "recorded_warnings.hxx" // warn, error

#ifndef NO_UNIT_TESTS
    #include <cassert> // assert
#endif // NO_UNIT_TESTS

namespace sigma_config {
    
    char const ellchar[12] = "spdfghijkl?";
    
    typedef struct {
        double  rcut;
        double  sigma;
        double  Z;
        double  q_core_hole;
        uint8_t nn[8];
        int8_t  ncmx[4];
        int8_t  iln_core_hole;
    } element_t;

    inline char const * default_config(unsigned const iZ) {
        switch (iZ) {
            case  1: return "1s* 1 0 2p | 0.9 sigma .247";                               // H   
            case  2: return "1s* 2 2p | 1.5 sigma .48";                                  // He  
            case  3: return "2s* 1 0 2p 1e-99 | 2.0 sigma .6";                           // Li  
            case  4: return "2s* 2 2p 1e-99 0 | 1.5 sigma .45";                          // Be  
            case  5: return "2s* 2 2p* 1 0 3d | 1.2 sigma .45";                          // B   
            case  6: return "2s* 2 2p* 2 0 3d | 1.2 sigma .43";                          // C   
            case  7: return "2s* 2 2p* 3 0 3d | 1.0 sigma .33";                          // N   
            case  8: return "2s* 2 2p* 3 1 3d | 1.13 sigma .297";                        // O   
            case  9: return "2s* 2 2p* 3 2 3d | 1.2 sigma .323";                         // F   
            case 10: return "2s* 2 2p* 6 3d | 1.8 sigma .564";                           // Ne  
            case 11: return "3s* 1 0 2p 6 3p 1e-99 3d | 2.27 sigma .69";                 // Na  
            case 12: return "2s 2 3s 2 2p 6 3p 1e-99 3d | 1.96 sigma .41";               // Mg  
            case 13: return "3s* 2 3p* 1 0 3d | 2.05 sigma .645";                        // Al  
            case 14: return "3s* 2 3p* 2 0 3d | 2.0 sigma .643";                         // Si  
            case 15: return "3s* 2 3p* 3 0 3d | 1.8 sigma .512";                         // P   
            case 16: return "3s* 2 3p* 3 1 3d | 1.7 sigma .535";                         // S   
            case 17: return "3s* 2 3p* 3 2 3d | 1.5 sigma .503";                         // Cl  
            case 18: return "3s* 2 3p* 6 3d | 1.6 sigma .546";                           // Ar  
            case 19: return "3s 2 4s 1 0 3p* 6 3d | 1.77 sigma .47";                     // K   
            case 20: return "3s 2 4s 2 3p* 6 3d | 1.77 sigma .487";                      // Ca  
            case 21: return "3s 2 4s 2 3p 6 4p 1e-99 3d* 1 0 | 2.32 sigma .58";          // Sc  
            case 22: return "3s 2 4s 2 3p 6 4p 1e-99 3d* 2 0 | 2.0 sigma .58";           // Ti  
            case 23: return "3s 2 4s 2 3p 6 4p 1e-99 3d* 3 0 | 2.4 sigma .56";           // V   
            case 24: return "4s* 1 0 4p* 1e-99 3d* 5 0 | 2.1 sigma .667";                // Cr  
            case 25: return "3s 2 4s 2 3p 6 4p 1e-99 3d* 5 0 | 2.41 sigma .554";         // Mn  
            case 26: return "4s* 2 4p* 1e-99 3d* 5 1 | 2.0 sigma .65";                   // Fe  
            case 27: return "4s* 2 4p* 1e-99 3d* 5 2 | 1.9 sigma .608";                  // Co  
            case 28: return "4s* 2 3p 6 4p 1e-99 3d* 5 3 | 2.15 sigma .48";              // Ni  
            case 29: return "4s* 1 0 4p* 1e-99 3d* 10 | 2.0 sigma .61";                  // Cu  
            case 30: return "4s* 1 1 4p* 1e-99 3d* 10 | 2.23 sigma .577";                // Zn  
            case 31: return "4s* 2 4p* 1 0 4d | 2.2 sigma .686";                         // Ga  
            case 32: return "4s* 2 4p* 2 0 4d | 1.9 sigma .606";                         // Ge  
            case 33: return "4s* 2 4p* 3 0 4d | 2.0 sigma .62";                          // As  
            case 34: return "4s* 2 4p* 3 1 4d | 1.6 sigma .521";                         // Se  
            case 35: return "4s* 2 4p* 3 2 4d | 2.1 sigma .6";                           // Br  
            case 36: return "4s* 2 4p* 6 4d | 2.2 sigma .61";                            // Kr 
            case 37: return "4s 2 5s 1 0 4p* 6 4d | 2.3 sigma .78";                      // Rb 
            case 38: return "4s 2 5s 2 4p* 6 4d | 2.37 sigma .666";                      // Sr
            case 39: return "4s 2 5s 2 4p 6 5p 1e-99 4d* 1 0 | 2.43 sigma .6";           // Y 
            case 40: return "4s 2 5s 2 4p 6 5p 1e-99 4d* 2 0 | 2.35 sigma .58";          // Zr
            case 41: return "4s 2 5s 1 0 4p 6 5p 1e-99 4d* 4 0 | 2.35 sigma .59";        // Nb
            case 42: return "4s 2 5s 1 0 4p 6 5p 1e-99 4d* 5 0 | 2.34 sigma .585";       // Mo
            case 43: return "4s 2 5s 1 0 4p 6 5p 1e-99 4d* 5 1 | 2.4 sigma .58";         // Tc
            case 44: return "4s 2 5s 1 0 4p 6 5p 1e-99 4d* 5 2 | 2.37 sigma .571";       // Ru
            case 45: return "5s* 1 0 4p 6 5p 1e-99 4d* 5 3 | 2.35 sigma .58";            // Rh
            case 46: return "5s* 1e-99 0 4p 6 5p 1e-99 4d* 10 | 2.32 sigma .585";        // Pd
            case 47: return "5s* 1 0 4p 6 5p 1e-99 4d* 10 | 2.23 sigma .57";             // Ag
            case 48: return "5s* 1 1 5p* 1e-99 4d* 10 | 2.2 sigma .563";                 // Cd 
            case 49: return "5s* 2 5p* 1 0 4d* 10 | 2.17 sigma .565";                    // In 
            case 50: return "5s* 2 5p* 2 0 4d* 10 | 2.24 sigma .585";                    // Sn
            case 51: return "5s* 2 5p* 3 0 4d* 10 | 2.18 sigma .57";                     // Sb
            case 52: return "5s* 2 5p* 3 1 5d | 2.23 sigma .555";                        // Te
            case 53: return "5s* 2 5p* 3 2 5d | 2.2 sigma .68";                          // I 
            case 54: return "5s* 2 5p* 6 5d | 2.24 sigma .62";                           // Xe 
            case 55: return "5s 2 6s 1 0 5p* 6 5d | 2.0 sigma .61";                      // Cs 
            case 56: return "5s 2 6s 2 5p* 6 5d* | 2.2 sigma .645";                      // Ba  
            case 57: return "5s 2 6s 2 5p* 6 5d* 1 0 | 1.9 sigma .59";                   // La  
            // elements Ce_58 through Yb_70 missing
            case 71: return "5s 2 6s 2 5p 6 6p 1e-99 5d* 1 0 | 2.4 sigma .6";            // Lu  
            case 72: return "5s 2 6s 2 5p 6 6p 1e-99 5d* 2 0 | 2.47 sigma .6077";        // Hf  
            case 73: return "5s 2 6s 2 5p 6 6p 1e-99 5d* 3 0 | 2.47 sigma .6";           // Ta  
            case 74: return "5s 2 6s 2 5p 6 6p 1e-99 5d* 4 0 | 2.32 sigma .62";          // W   
            case 75: return "6s* 2 5p 6 6p 1e-99 5d* 5 0 | 2.47 sigma .63";              // Re  
            case 76: return "6s* 2 5p 6 6p 1e-99 5d* 5 1 | 2.35 sigma .58";              // Os  
            case 77: return "6s* 2 5p 6 6p 1e-99 5d* 5 2 | 2.43 sigma .62";              // Ir  
            case 78: return "6s* 1 0 5p 6 6p 1e-99 5d* 5 4 | 2.47 sigma .59";            // Pt  
            case 79: return "6s* 1 0 6p* 1e-99 5d* 10 | 2.5 sigma .667";                 // Au  
            case 80: return "6s* 2 5p 6 6p 1e-99 5d* 10 | 2.44 sigma .59";               // Hg  
            case 81: return "6s* 2 6p* 1 0 5d* 10 | 2.25 sigma .567";                    // Tl  
            case 82: return "6s* 2 6p* 2 0 5d* 10 | 2.3 sigma .59";                      // Pb  
            case 83: return "6s* 2 6p* 3 0 5d* 10 | 2.41 sigma .605";                    // Bi  
            case 84: return "6s* 2 6p* 3 1 5d* 10 | 2.3 sigma .54";                      // Po  
            case 85: return "6s* 2 6p* 3 2 5d* 10 | 2.3 sigma .54";                      // At  
            case 86: return "6s* 2 6p* 6 6d | 2.29 sigma .54";                           // Rn  
            default:
                     warn("no default element configuration given for Z=%d", iZ);
                     return "";
        } // switch
    } // default config

    inline int8_t char2ell(char const c) {
        switch (c | 32) { // ignore if uppercase or lowercase
            case ' ': case '\0': case '\t': case '\n': return -1;
            case 'e': return -9; // used for exponential numeric values, otherwise "1e-99" would be an orbital with unphysical ell
            case 's': return 0;
            case 'p': return 1;
            case 'd': return 2;
            default : return (c | 32) - 'c'; // 'f' -> 3, 'g' -> 4, ...
        } // switch
    } // char2ell

    inline int8_t char2enn(char const c) {
        switch (c) { // case sensitive
            case ' ': case '\0': case '\t': case '\n': return 0;
            case 'r': case 'R': case '|': return -1;
            case 's': case 'S': return -2;
            case 'Z': case 'z': return -3;
            case 'V': case 'v': return -4;
            case '0': case '.': return -9; // numeric reading
            default : return c - '0';
        } // switch
    } // char2enn
    
    
    typedef struct {
        double value;
        int8_t enn, ell, nrn, key;
    } orbital_t;

    inline element_t& get(double const Z, int const echo=0) {
        
        char symbol[4], element_Sy[16];
        int const iZ = chemical_symbol::get(symbol, Z);
        std::snprintf(element_Sy, 15, "element_%s");
        auto const config = control::get(element_Sy, default_config(iZ));

        if (echo > 0) printf("# for Z=%g use %s=\"%s\"\n", Z, element_Sy, config);
        // now convert config into an element_t
        auto e = new element_t;
        
        e->rcut = 2.;
        e->sigma = .5;
        e->Z = iZ;
        e->q_core_hole = 0;
        set(e->nn, 8, uint8_t(0));
        e->ncmx[0] = (Z > 2) + (Z > 4) + (Z > 12) + (Z > 20) + (Z > 38) + (Z > 56) + (Z > 88) + (Z > 120);
        e->ncmx[1] = 1 + (Z > 10) + (Z > 18) + (Z > 36) + (Z > 54) + (Z > 86) + (Z > 118);
        e->ncmx[2] = 2 + (Z > 30) + (Z > 48) + (Z > 80) + (Z > 112);
        e->ncmx[3] = 3 + (Z > 70) + (Z > 102);
        e->iln_core_hole = -1;

        if (nullptr == config) return *e;

        int constexpr mwords = 32;
        double values[mwords];
        int8_t enns[mwords];
        int8_t ells[mwords];
        int8_t mrns[mwords];

        int iword{0};
        
        char const * string{config};
        char c0{*string};
        while(c0) {
//          if (echo > 0) printf("# start from '%s'\n", string); fflush(stdout);
          
            assert(iword < mwords);
          
            values[iword] = 0.0;
            enns[iword] =  0; // undefined
            ells[iword] = -1; // undefined
            mrns[iword] = -1; // undefined
          
            bool try_numeric{false};
            char const cn = *string;
            auto const enn = char2enn(cn);
//                 if (echo > 0) printf("# found enn=%i in '%s'\n", enn, string); fflush(stdout);
            if (enn > 0 && enn <= 9) {
                char const cl = *(string + 1);
                auto const ell = char2ell(cl);
                if (ell >= 0) {
                    if (ell >= enn) error("unphysical ell=%i >= enn=%i in '%s'!", ell, enn, string);
                    ells[iword] = ell;
                    
                    char const cm = *(string + 2);
                    int const mrn = ('*'== cm); // max radial nodes
                    if (echo > 0) printf("# found enn=%i ell=%i mrn=%i in '%c%c%c'\n", enn, ell, mrn, cn,cl,cm); fflush(stdout);
                    mrns[iword] = mrn;
// later
//                         e->nn[ell] += (1 + mrn);
//                         if (ell < 4) e->ncmx[ell] = enn - 1;

                } else { // ell is a valid angular momentum quantum number
                    try_numeric = true;
                }
            } else { // enn is valid principal quantum number
                try_numeric = (-9 == enn);
                if (!try_numeric && echo > 0) printf("# found special '%s'\n", string); fflush(stdout);
            }
            enns[iword] = enn;

            if (try_numeric) {
                double const value = std::atof(string);
                if (echo > 0) printf("# found numeric value %g in '%s'\n", value, string); fflush(stdout);
                values[iword] = value; // store
                enns[iword] = -9; // -9:numeric
            } // try_numeric

            ++iword;
            
            auto const next_blank = std::strchr(string, ' '); // forward to the next word
            if (next_blank) {
                string = next_blank + 1;
                c0 = *string;
            } else {
                c0 = 0; // stop the while loop
            }
        } // while;
        int const nwords = iword; // how many words were in the string
        if (echo > 0) printf("# process %d words\n", nwords); fflush(stdout);
        
        int constexpr max_enn = 9, max_inl = (max_enn*(max_enn + 1))/2;
        double occ[max_inl][2];
        for (int inl = 0; inl < max_inl; ++inl) { occ[inl][0] = occ[inl][1] = 0; }

        if (echo > 0) {
            // repeat what was just parsed
            char const special[] = "_|sZV";
            printf("# repeat config string '");
            for (int iword = 0; iword < nwords; ++iword) {
                int8_t const enn = enns[iword];
                if (enn > 0) {
                    printf("%i%c%s ", enn, ellchar[ells[iword]], mrns[iword]?"*":"");
                } else if (-9 == enn) {
                    printf("%g ", values[iword]);
                } else {
                    printf("%c ", special[-enn]);
                }
            } // iword
            printf("'\n");
        } // echo
        
        
        typedef enum {
            ExpectOrbital = 1,
            ExpectDownOcc = 2,
            ExpectUpOcc = 3,
            ExpectRadius = -1,
            ExpectSigma = -2,
            ExpectZcore = -3,
            ExpectMethod = -4 // not supported
        } parse_state_t;
        
        parse_state_t state{ExpectOrbital};
       
//         int8_t enn{0}, ell{-1}, nrn{-1};
  
        double stack[4] = {0, 0, 0, 0};
        int nstack{0};

        // now put into context:
        for(int iword = nwords - 1; iword >= 0; --iword) {
            int8_t const enn = enns[iword];
            if (-9 == enn) {
//                 if (echo > 0) printf("# nstack=%i push\n", nstack); fflush(stdout);
                stack[nstack++] = values[iword]; // push to the stack
//                 if (echo > 0) printf("# nstack=%i stack[%i]=%g\n", nstack, nstack - 1, stack[nstack - 1]); fflush(stdout);
            } else {
                double value{0};
                if (nstack > 0) {
//                     if (echo > 0) printf("# nstack=%i stack[%i]=%g\n", nstack, nstack - 1, stack[nstack - 1]); fflush(stdout);
                    value = stack[--nstack]; // pop the stack
//                     if (echo > 0) printf("# nstack=%i popped\n", nstack); fflush(stdout);
                }
                if (enn > 0) { // orbital
                    int const ell = ells[iword];
                    int const inl = ell + (enn*(enn-1))/2; assert(ell < enn);
                    if (nstack > 0) {
                        occ[inl][0] = value;
//                         if (echo > 0) printf("# nstack=%i stack[%i]=%g\n", nstack, nstack - 1, stack[nstack - 1]); fflush(stdout);
                        value = stack[--nstack]; // pop the stack a 2nd time
//                         if (echo > 0) printf("# nstack=%i popped\n", nstack); fflush(stdout);
                        occ[inl][1] = value;
                    } else {
                        occ[inl][0] = .5*value;
                        occ[inl][1] = .5*value;
                    }
                    for (int spin = 0; spin < 2; ++spin) {
                        if (occ[inl][spin] < 0) {
                            warn("found a negative occupation number %g in the %s-%i%c orbital", occ[inl][spin], symbol, enn, ellchar[ell]);
                            occ[inl][spin] = 0;
                        }
                        if (occ[inl][spin] > 2*ell + 1) {
                            warn("occupation number %g is too large for a %s-%i%c orbital", occ[inl][spin], symbol, enn, ellchar[ell]);
                            occ[inl][spin] = 2*ell + 1;
                        }
                    } // spin
                    if (echo > 0) printf("# found orbital %i%c occ= %g %g inl=%i\n", enn,ellchar[ell], occ[inl][0], occ[inl][1], inl); fflush(stdout);
                    if (ell < 8) e->nn[ell] += 1 + mrns[iword];
                    if (ell < 4) e->ncmx[ell] = enn - 1;
                } else if (-1 == enn) {
                    e->rcut = value;     
                    if (echo > 0) printf("# found cutoff radius rcut = %g\n", e->rcut); fflush(stdout);
                    if(e->rcut <= 0) warn("rcut must be positive but found rcut=%g", e->rcut);
                } else if (-2 == enn) {  
                    e->sigma = value;
                    if (echo > 0) printf("# found projector spread sigma = %g\n", e->sigma); fflush(stdout);
                    if(e->sigma <= 0) warn("sigma must be positive but found sigma=%g", e->sigma);
                } else if (-3 == enn) {  
                    e->Z = value;
                    if (echo > 0) printf("# found core charge Z = %g\n", e->Z); fflush(stdout);
                    if(e->Z >= 128) warn("some routine may not be prepared for Z = %g >= 128", e->Z);
                }
            }
        } // iword
        if (nstack > 0) warn("after parsing, some value %g was left on the stack", stack[0]);
        
        
        // fill the core
        if (echo > 0) printf("# fill the core up to enn = %d %d %d %d for s,p,d,f\n", e->ncmx[0], e->ncmx[1], e->ncmx[2], e->ncmx[3]); fflush(stdout);
        for (int ell = 0; ell < 4; ++ell) {
            for(int enn = ell + 1; enn <= e->ncmx[ell]; ++enn) {
                int const inl = ell + (enn*(enn-1))/2; assert(ell < enn);
                occ[inl][0] = occ[inl][1] = -(2*ell + 1); // negative occupation numbers indicate core electrons
            } // enn
        } // ell

        double nve{0}, nce{0}; // number of valence and core electrons
        for(int inl = 0; inl < max_inl; ++inl) {
            nve += std::max(0.0, occ[inl][0]) + std::max(0.0, occ[inl][1]);
            nce -= std::min(0.0, occ[inl][0]) + std::min(0.0, occ[inl][1]);
        } // inl
        double const nelectrons = nve + nce;
        if (echo > 0) printf("# found %g = %g core + %g valence electrons and Z = %g\n", nelectrons, nce, nve, e->Z); fflush(stdout);

        if (echo > 0) {
            printf("# PAW setup for %s (Z=%g) uses", symbol, e->Z);
            for(int ell = 0; ell < 8; ++ell) {
                printf(" %d", e->nn[ell]);
            } // ell
            printf(" partial waves\n");
        } // echo
        
        if (std::abs(nelectrons - e->Z) > 1e-12) warn("PAW setup for %s (Z=%g) is charged with %g electrons", symbol, e->Z, nelectrons - e->Z);
        
        return *e;
    } // get

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  inline status_t test_86(int const echo=0) {
      for (int iZ = 86; iZ > 0; --iZ) {
          if (std::abs(iZ - 64) >= 7) {
              if (echo > 0) printf("\n");
              auto const e = get(iZ, echo);
          } // without 58 through 70
      } // iZ
      return 0;
  } // test_86

  inline status_t all_tests(int const echo=0) {  
      status_t stat(0);
      stat += test_86(echo);
      return stat;
  } // all_tests

#endif

} // namespace sigma_config
