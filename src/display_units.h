#pragma once

#if 1
double const eV = 27.210282768626; char const _eV[] = "eV"; // eVolt energy unit
double const Ang = 0.52917724924; char const _Ang[] = "Ang"; // Angstrom length unit
// double const eV = 2; char const _eV[] = "Ryd"; // Rydberg atomic energy unit
// double const Ang = .052917724924; char const _Ang[] = "nm"; // nanometer length unit
#else
double const eV = 1; char const _eV[] = "Ha"; // Hartree atomic energy unit
double const Ang = 1; char const _Ang[] = "Bohr"; // Bohr atomic length unit
#endif
