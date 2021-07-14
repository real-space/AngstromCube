#pragma once

#include <cstdio> // std::printf
#include <algorithm> // std::max

#include "chemical_symbol.h" // element_symbols
#include "status.hxx" // status_t

namespace element_config {
  
  inline char const* element_symbol(int const Z) {
      int const z128 = Z & 127; // == Z % 128
      return &(element_symbols[2*z128]); // this and the following char belong to this element
  } // element_symbol

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_show_all_elements(int const echo=2) {
      if (echo > 0) std::printf("# %s:  %s\n", __FILE__, __func__);
      float iocc[32];
      char const ellchar[] = "spdfgh+";

      for (int spin = 1; spin >= 0; --spin) {
          auto const echo_occ = 9; // how high must the log-level be to display also the occupation numbers?
          if (echo > 2 + 3*spin) { // how high must the log-level be to display also the spin-orbit table?
              if (spin > 0) std::printf("# including spin-orbit\n");
              std::printf("\n#  nl    j   occ    index\n");
              for (int inl = 0; inl < 32; ++inl) iocc[inl] = 0; // clear array
              int iZ{0}; // atomic number
              int ishell{0}; // shell index
              for (int m = 0; m < 8; ++m) { // auxiliary number
                  int enn{(m + 1)/2};
                  for (int ell = m/2; ell >= 0; --ell) { // angular momentum character
                      ++enn; // principal quantum number
                      for (int jj = 2*ell + spin; jj >= std::max(0, 2*ell - spin); jj -= 2) { // quantum number j = l + s, jj=2j
                          std::printf("%4d%c%6.1f%4d%9d    %c", enn, ellchar[ell], jj*.5f, 
                                      (2 - spin)*(jj + 1), ishell, (echo > echo_occ)?'\n':' ');
                          for (int mm = -jj; mm <= jj; mm += 2) { // angular momentum z-component emm, mm=2*emm
                              for (int s = 0; s <= 1 - spin; ++s) { // when 0==spin, this loop fills each orbital with 2 electrons
                                  iocc[ishell] += 1;
                                  ++iZ; // next element
                                  auto const El = element_symbol(iZ);
                                  if (echo > echo_occ) {
                                      std::printf("# %c%c occupation", El[0], El[1]);
                                      for (int inl = 0; inl <= ishell; ++inl) {
                                          std::printf(" %g", iocc[inl]);
                                      } // inl
                                      std::printf("\n");
                                  } else {
                                      std::printf(" %c%c", El[0], El[1]);
                                  } // echo_occ
                              } // s
                          } // mm
                          ++ishell; // next shell
                          std::printf("\n");
                      } // jj
                  } // ell
              } // m
              double sum_occ{0};
              for (int inl = 0; inl < ishell; ++inl) {
                  sum_occ += iocc[inl];
              } // inl
              std::printf("# sum(occ) = %.3f\n", sum_occ);
              assert(120 == sum_occ); // sanity check
              assert(20 + 12*spin == ishell); // sanity check
          } // echo > echos
      } // spin
      if (echo > 0) std::printf("# table of elements according to Aco Z. Muradjan\n\n");
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

  inline status_t all_tests(int const echo=0) {
      return test_show_all_elements(echo);
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace element_config

#ifdef DEBUG
  #undef DEBUG
#endif
#ifdef FULL_DEBUG
  #undef FULL_DEBUG
#endif
