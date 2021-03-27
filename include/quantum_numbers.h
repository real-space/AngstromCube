#pragma once

#include <stdint.h> // uint8_t, int8_t

/* Quantum Numbers */
typedef uint8_t enn_QN_t; // energy quantum number
typedef uint8_t ell_QN_t; // angular momentum quantum number
typedef int8_t  emm_QN_t; // angular momentum L_z component
typedef int8_t  spin_QN_t; // spin quantum number

emm_QN_t constexpr emm_Degenerate = -99;
spin_QN_t constexpr spin_Degenerate = -9;
