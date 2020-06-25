#pragma once

namespace constants {
  
    double constexpr pi = 3.14159265358979323846; // pi
    double constexpr sqrt2 = 1.4142135623730950488; // sqrt(2)
    double constexpr sqrtpi = 1.77245385090551602729816748334115;

//  char constexpr ignore_case = 32; // use c | ignore_case to convert all upper case to lower case letters

  double constexpr electron_Volt = 27.210282768626; 
  double constexpr nano_meter = .052917724924; 
  double constexpr Angstrom = 10*nano_meter;
  
  char constexpr _electron_Volt[] = "eV";
  char constexpr _nano_meter[] = "nm";
  char constexpr _Rydberg[] = "Ha"; // Rydberg atomic energy unit
  char constexpr _Hartree[] = "Ha"; // Hartree atomic energy unit
  char constexpr _Bohr[] = "Bohr"; // Bohr atomic length unit
    
} // namespace constants
