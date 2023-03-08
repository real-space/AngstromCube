#pragma once
// This file is part of AngstromCube under MIT License

#include "quantum_numbers.h" // enn_QN_t, ell_QN_t

  int constexpr TRU=0, SMT=1;
  int constexpr TRU_AND_SMT=2, TRU_ONLY=1;

  template <int Pseudo> // Pseudo: see description at instanciation
  struct energy_level_t {
      double* wave[Pseudo]; // for valence states points to the true and smooth partial waves
      double* wKin[Pseudo]; // kinetic energy operator onto wave, r*T*wave
      double energy; // energy level in Hartree atomic units
      double kinetic_energy; // true kinetic energy
      double occupation; // occupation number
      char tag[8]; // label
      enn_QN_t nrn[Pseudo]; // number of radial nodes
      enn_QN_t enn; // principal quantum_number
      ell_QN_t ell; // angular momentum quantum_number
      int8_t csv; // 0:core, 1:semicore, 2:valence, 3:undefined
  }; // struct energy_level_t

// now defined in spherical_state.hxx
//   typedef struct energy_level_t<TRU_ONLY> spherical_state_t; // Pseudo=1: spherical states, e.g. core states, only a TRU wave is stored

  typedef struct energy_level_t<TRU_AND_SMT> partial_wave_t; // Pseudo=2: states with a smooth counterpart, e.g. partial waves describing valence states
