#pragma once

#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t, emm_Degenerate, spin_QN_t, spin_Degenerate

  int constexpr TRU=0, SMT=1, TRU_AND_SMT=2;

  template<int Pseudo> // Pseudo: see description at instanciation
  struct energy_level {
      double* wave[Pseudo]; // for valence states points to the true and smooth partial waves
      double* wKin[Pseudo]; // kinetic energy operator onto wave, r*T*wave
      double energy; // energy level in Hartree atomic units
      double occupation; // occupation number
      char tag[8]; // label
      enn_QN_t enn; // main quantum_number
      ell_QN_t ell; // angular momentum quantum_number
      emm_QN_t emm; // usually emm == emm_Degenerate
      spin_QN_t spin; // usually spin == spin_Degenerate
      enn_QN_t nrn[Pseudo]; // number of radial nodes
      int8_t csv; // 0:core, 1:semicore, 2:valence, 3:undefined
  }; // struct energy_level

  typedef struct energy_level<1> spherical_orbital_t; // Pseudo=1: spherical states, e.g. core states, only a TRU wave is stored
  typedef struct energy_level<TRU_AND_SMT> partial_wave_t; // Pseudo=2: states with a smooth conterpart, e.g. partial waves describing valence states
