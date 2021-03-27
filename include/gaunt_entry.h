#pragma once

#include <stdint.h> // int16_t, int32_t

#ifdef __cplusplus
extern "C" {
#endif


typedef struct { 
  double G;
  int32_t lm;
  int16_t lm1, lm2;
} gaunt_entry_t;

#ifdef __cplusplus
} // extern "C"
#endif
