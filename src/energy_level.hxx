#pragma once

#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t, emm_Degenerate, spin_QN_t, spin_Degenerate

  int constexpr TRU=0, SMT=1, TRU_AND_SMT=2, TRU_ONLY=1;

  // ToDo: write a class with constructors and destructors to handle the memory for *wave
  template<int Pseudo> // Pseudo=1: core states, Pseudo=2: valence states
  struct energy_level {
      double* wave[Pseudo]; // for valence states points to the true and smooth partial waves
      double* wKin[Pseudo]; // kinetic energy operator onto r*wave
      double energy; // energy level in Hartree atomic units
      double occupation; // occupation number
      enn_QN_t enn; // main quantum_number
      ell_QN_t ell; // angular momentum quantum_number
      emm_QN_t emm; // usually emm == emm_Degenerate
      spin_QN_t spin; // usually spin == spin_Degenerate
      enn_QN_t nrn[Pseudo]; // number of radial nodes
  }; // struct energy_level

  typedef struct energy_level<TRU_ONLY> core_level_t;
  typedef struct energy_level<TRU_AND_SMT> valence_level_t;
