#include <cstdio> // printf
#include <cstdlib> 
#include <cmath>

#include "element_configuration.hxx"

#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
#include "output_units.h" // eV, _eV, Ang, _Ang
#include "min_and_max.h" // max
#include "element_symbols.h" // element_symbols

// #define FULL_DEBUG
// #define DEBUG

typedef struct {
  char A;
  char b;
} element_symbol_t;

namespace element_configuration {
  
  element_symbol_t element_symbol(int const Z) {
      int const z128 = Z & 127; // == Z % 128
      element_symbol_t El;
      El.A = element_symbols[2*z128];
      El.b = element_symbols[2*z128 + 1];
      return El;
  } // element_symbol

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_show_all_elements() {
    float iocc[32];
    char const ellchar[] = "spdfghi";
    int constexpr echo = 2; // 0: no output
    
    for(int spin = 1; spin >= 0; --spin) {
      if (echo) {
          if (spin > 0) printf("# including spin-orbit\n");
          printf("\n   nl    j    occ  index\n");
      } // echo
      for(int inl = 0; inl < 32; ++inl) iocc[inl] = 0; // clear array
      int iZ = 0; // atomic number
      int ishell = 0; // shell index
      for(int m = 0; m < 8; ++m) { // auxiliary number
        int enn = (m + 1)/2;
        for(int ell = m/2; ell >= 0; --ell) { // angular momentum character
          ++enn; // principal quantum number
          for(int jj = 2*ell + spin; jj >= max(0, 2*ell - spin); jj -= 2) { // quantum number j = l + s, jj=2j
            if(echo) 
                printf("%4d%c%6.1f%4d%9d    %c", enn, ellchar[ell], jj*.5f, (2 - spin)*(jj + 1), ishell, (1==echo || spin)?'\n':' ');
            for(int mm = -jj; mm <= jj; mm += 2) { // angular momentum z-component emm, mm=2*emm
              for(int s = 0; s <= 1 - spin; ++s) { // when 0==spin, this loop fills each orbital with 2 electrons
                iocc[ishell] += 1;
                ++iZ; // next element
                if (echo) {
                    auto const El = element_symbol(iZ);
                    if (1 == echo) {
                        printf("%c%c occupation ", El.A, El.b);
                        for(int inl = 0; inl < ishell; ++inl)
                            printf("%.1f ", iocc[inl]);
                        printf("\n");
                    } else if ((echo > 1) && (0 == spin)) {
                        printf("%c%c ", El.A, El.b);
                    }
                } // echo
              } // s
            } // mm
            ++ishell; // next shell
            if ((echo > 1) && (0 == spin)) {
                printf("\n");
            } // echo
          } // jj
        } // ell
      } // m
      if(echo) {
        double sum_occ = 0;
        for(int inl = 0; inl < ishell; ++inl) {
            sum_occ += iocc[inl];
        } // inl
        printf("# sum(occ) = %.3f\n", sum_occ);
      } // echo
    } // spin
      printf("# table of elements according to Aco Z. Muradjan\n\n");
      return 0;
  } // test_show_all_elements

//  core electron configurations for predefined cores:
// 
//  rare gases        [frozen d]        [frozen f]
//  __ '    '   0
//  He '1   '   2
//  Ne '22  '  10
//  Ar '33  '  18     3d '333 '  28
//  Kr '443 '  36     4d '444 '  46
//  Xe '554 '  54     5d '5554'  78     4f '5544'  68
//  Rn '6654'  86     6d '6665' 110     6d '6655' 100
//  uo '7765' 118
// 
//   1s
//   2s 2p
//   3s 3p 3d
//   4s 4p 4d 4f
//   5s 5p 5d 5f
//   6s 6p 6d
//   7s 7p
//   8s

  status_t all_tests() {
    return test_show_all_elements();
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace element_configuration
