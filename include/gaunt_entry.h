#pragma once

#include <stdint.h> // int16_t, int32_t

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
      double G; // Gaunt coefficient
      int32_t lm; // combined index ell0 emm0
      int16_t lm1, lm2; // combined indices ell1 emm1 and ell2 emm2
  } gaunt_entry_t;

#ifdef __cplusplus
} // extern "C"
#endif
