#pragma once

#include <cassert> // assert

#include "display_units.h" // eV, _eV, Ang, _Ang, Kelvin, _Kelvin
#include "status.hxx" // status_t
#include "recorded_warnings.hxx" // warn

namespace unit_system {
  
  double constexpr Hartree2electronVolt = 27.210282768626; // conversion factor from Hartree to eVolt

  char constexpr _electron_Volt[] = "eV"; // electron Volt energy unit
  char constexpr _Rydberg[] = "Ry"; // Rydberg atomic energy unit
  char constexpr _Hartree[] = "Ha"; // Hartree atomic energy unit
  
  inline double energy_unit(char const *which, char const **const symbol) {
      char const w = which[0] | 32; // take first char and ignore case
      if ('e' == w) {
          *symbol = _electron_Volt;   return Hartree2electronVolt; // "eV"
      } else if ('r' == w) {
          *symbol = _Rydberg;         return 2; // "Ry"
      } else if ('k' == w) {
          *symbol = _Kelvin;          return Kelvin; // "Kelvin"
      } else {
          if ('h' != w) warn("unknown energy unit \"%s\", default to Ha (Hartree)", which);
          *symbol = _Hartree;         return 1; // "Ha"
      } // w
  } // energy_unit

  double constexpr Bohr2nanometer = .052917724924; // conversion factor from Bohr to nm

  char constexpr _nano_meter[] = "nm"; // nanometer length unit
  char constexpr _Angstrom[] = "Ang"; // Angstrom length unit
  char constexpr _Bohr[] = "Bohr"; // Bohr atomic length unit
  
  inline double length_unit(char const *which, char const **const symbol) {
      char const w = which[0] | 32; // take first char and ignore case
      if ('a' == w) {
          *symbol = _Angstrom;    return 10*Bohr2nanometer; // "Ang"
      } else if ('n' == w) {
          *symbol = _nano_meter;  return Bohr2nanometer; // "nm"
      } else {
          if ('b' != w) warn("unknown length unit \"%s\", default to Bohr", which);
          *symbol = _Bohr;        return 1; // "Bohr"
      } // w
  } // length_unit

  inline status_t set_output_units(char const *energy, char const *length) {
#ifdef _Output_Units_Fixed
      if ('H' != *energy || 'B' != *length) {
          warn("output units cannot be changed to {%s, %s} at runtime", energy, length);
      } // units differ from "Ha" and "Bohr"
      return -1; // cannot modify
#else
      eV  = energy_unit(energy, &_eV);
      Ang = length_unit(length, &_Ang);
      assert( eV  > 0 );
      assert( Ang > 0 );
      return 0;
#endif
  } // set_output_units

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_all_combinations(int const echo=0) {
      status_t stat(0);
      char const sy[2][3][8] = {{"Ha", "Ry", "eV"}, {"Bohr", "Ang", "nm"}};
      char const el_name[][8] = {"energy", "length"};
      for(int el = 0; el < 2; ++el) {
          for(int iu = 0; iu < 3; ++iu) {
              char const *_si;
              auto const f = el ? length_unit(sy[el][iu], &_si) : energy_unit(sy[el][iu], &_si);
              if (echo > 2) printf("# %s unit factor %.12f for %s\n", el_name[el], f, _si);
              auto const fi = 1./f;
              for(int ou = 0; ou < 3; ++ou) {
                  char const *_so;
                  auto const fo = el ? length_unit(sy[el][ou], &_so) : energy_unit(sy[el][ou], &_so);
                  if (echo > 3) printf("# %s unit factors %.9f * %.9f = %.9f %s/%s\n", el_name[el], fo, fi, fo*fi, _si, _so);
                  stat += (iu == ou)*((fo*fi - 1)*(fo*fi - 1) > 4e-32); // check that inverse factors produce 1.0 exactly
              } // ou
          } // iu input unit
          if (echo > 2) printf("#\n");
      } // el {energy, length}
      return stat;
  } // test_all_combinations

  inline status_t all_tests(int const echo=0) {
      if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
      status_t stat(0);
      stat += test_all_combinations(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  
  
} // namespace unit_system
