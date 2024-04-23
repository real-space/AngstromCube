#pragma once
// This file is part of AngstromCube under MIT License

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _Output_Units_Fixed
  double constexpr eV  = 1; char const _eV[] = "Ha"; // Hartree atomic energy unit
  double constexpr Ang = 1; char const _Ang[] = "Bohr"; // Bohr atomic length unit
#else
  // these global variables are only to be modified in unit_system.hxx
  extern double eV;  extern char const *_eV;  // dynamic energy unit
  extern double Ang; extern char const *_Ang; // dynamic length unit
#endif
  double constexpr Kelvin = 315773.244215; char const _Kelvin[] = "Kelvin";
  double constexpr GByte = 1e-9; char const *const _GByte = "GByte";
  double constexpr MByte = 1e-6; char const *const _MByte = "MByte";

#ifdef __cplusplus
} // extern "C"
#endif
