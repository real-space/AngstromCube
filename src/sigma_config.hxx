#pragma once

#include <cstdio> // printf, std::snprintf
#include <cmath> // std::floor
#include <cstdint> // int8_t

#include "chemical_symbol.hxx" // ::get
#include "control.hxx" // ::get
#include "recorded_warnings.hxx" // warn

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
            case 23: return "3s 2 4s 2 3p 6 4p 1d-99 3d* 3 0 | 2.4 sigma .56";           // V   
            case 24: return "4s* 1 0 4p* 1d-99 3d* 5 0 | 2.1 sigma .667";                // Cr  
            case 25: return "3s 2 4s 2 3p 6 4p 1d-99 3d* 5 0 | 2.41 sigma .554";         // Mn  
            case 26: return "4s* 2 4p* 1d-99 3d* 5 1 | 2.0 sigma .65";                   // Fe  
            case 27: return "4s* 2 4p* 1d-99 3d* 5 2 | 1.9 sigma .608";                  // Co  
            case 28: return "4s* 2 3p 6 4p 1d-99 3d* 5 3 | 2.15 sigma .48";              // Ni  
            case 29: return "4s* 1 0 4p* 1d-99 3d* 10 | 2.0 sigma .61";                  // Cu  
            case 30: return "4s* 1 1 4p* 1d-99 3d* 10 | 2.23 sigma .577";                // Zn  
            case 31: return "4s* 2 4p* 1 0 4d | 2.2 sigma .686";                         // Ga  
            case 32: return "4s* 2 4p* 2 0 4d | 1.9 sigma .606";                         // Ge  
            case 33: return "4s* 2 4p* 3 0 4d | 2.0 sigma .62";                          // As  
            case 34: return "4s* 2 4p* 3 1 4d | 1.6 sigma .521";                         // Se  
            case 35: return "4s* 2 4p* 3 2 4d | 2.1 sigma .6";                           // Br  
            case 36: return "4s* 2 4p* 6 4d | 2.2 sigma .61";                            // Kr  
            case 37: return "5s* 2 5p* 1 0 4d* 10 | 2.17 sigma .565";                    // In  
            case 38: return "5s* 2 5p* 2 0 4d* 10 | 2.24 sigma .585";                    // Sn  
            case 39: return "5s* 2 5p* 3 0 4d* 10 | 2.18 sigma .57";                     // Sb  
            case 40: return "5s* 2 5p* 3 1 5d | 2.23 sigma .555";                        // Te  
            case 41: return "5s* 2 5p* 3 2 5d | 2.2 sigma .68";                          // I   
            case 42: return "5s* 2 5p* 6 5d | 2.24 sigma .62";                           // Xe  
            case 43: return "4s 2 5s 1 0 4p* 6 4d | 2.3 sigma .78";                      // Rb  
            case 44: return "4s 2 5s 2 4p* 6 4d | 2.37 sigma .666";                      // Sr  
            case 45: return "4s 2 5s 2 4p 6 5p 1d-99 4d* 1 0 | 2.43 sigma .6";           // Y   
            case 46: return "4s 2 5s 2 4p 6 5p 1d-99 4d* 2 0 | 2.35 sigma .58";          // Zr  
            case 47: return "4s 2 5s 1 0 4p 6 5p 1d-99 4d* 4 0 | 2.35 sigma .59";        // Nb  
            case 48: return "4s 2 5s 1 0 4p 6 5p 1d-99 4d* 5 0 | 2.34 sigma .585";       // Mo  
            case 49: return "4s 2 5s 1 0 4p 6 5p 1d-99 4d* 5 1 | 2.4 sigma .58";         // Tc  
            case 50: return "4s 2 5s 1 0 4p 6 5p 1d-99 4d* 5 2 | 2.37 sigma .571";       // Ru  
            case 51: return "5s* 1 0 4p 6 5p 1d-99 4d* 5 3 | 2.35 sigma .58";            // Rh  
            case 52: return "5s* 1e-99 0 4p 6 5p 1e-99 4d* 10 | 2.32 sigma .585";        // Pd  
            case 53: return "5s* 1 0 4p 6 5p 1d-99 4d* 10 | 2.23 sigma .57";             // Ag  
            case 54: return "5s* 1 1 5p* 1d-99 4d* 10 | 2.2 sigma .563";                 // Cd  
            case 55: return "5s 2 6s 1 0 5p* 6 5d | 2.0 sigma .61";                      // Cs  
            case 56: return "5s 2 6s 2 5p* 6 5d* | 2.2 sigma .645";                      // Ba  
            case 57: return "5s 2 6s 2 5p* 6 5d* 1 0 | 1.9 sigma .59";                   // La  
            // elements Ce_58 through Yb_70 missing
            case 71: return "5s 2 6s 2 5p 6 6p 1d-99 5d* 1 0 | 2.4 sigma .6";            // Lu  
            case 72: return "5s 2 6s 2 5p 6 6p 1d-99 5d* 2 0 | 2.47 sigma .6077";        // Hf  
            case 73: return "5s 2 6s 2 5p 6 6p 1d-99 5d* 3 0 | 2.47 sigma .6";           // Ta  
            case 74: return "5s 2 6s 2 5p 6 6p 1d-99 5d* 4 0 | 2.32 sigma .62";          // W   
            case 75: return "6s* 2 5p 6 6p 1d-99 5d* 5 0 | 2.47 sigma .63";              // Re  
            case 76: return "6s* 2 5p 6 6p 1d-99 5d* 5 1 | 2.35 sigma .58";              // Os  
            case 77: return "6s* 2 5p 6 6p 1d-99 5d* 5 2 | 2.43 sigma .62";              // Ir  
            case 78: return "6s* 1 0 5p 6 6p 1d-99 5d* 5 4 | 2.47 sigma .59";            // Pt  
            case 79: return "6s* 1 0 6p* 1d-99 5d* 10 | 2.5 sigma .667";                 // Au  
            case 80: return "6s* 2 5p 6 6p 1d-99 5d* 10 | 2.44 sigma .59";               // Hg  
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

    inline element_t get(double const Z, int const echo=0) {
        
        char symbol[4], element_Sy[16];
        int const iZ = chemical_symbol::get(symbol, Z);
        std::snprintf(element_Sy, 15, "element_%s");
        auto const config = control::get(element_Sy, default_config(iZ));

        if (echo > 0) printf("# for Z=%g use %s=\"%s\"\n", Z, element_Sy, config);
        // now convert config into an element_t
        element_t e;
        
        return e;
    } // get

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  inline status_t test_86(int const echo=0) {
      for (int iZ = 1; iZ <= 86; ++iZ) {
          if (std::abs(iZ - 64) >= 7) {
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
