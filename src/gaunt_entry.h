#pragma once

#include <stdint.h> // int16_t, int32_t

typedef struct { 
  double G;
  int32_t lm;
  int16_t lm1, lm2;
} gaunt_entry_t;
