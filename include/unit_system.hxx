#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert

#include "status.hxx" // status_t

namespace unit_system {

    char constexpr _Rydberg[] = "Ry"; // Rydberg atomic energy unit

    double energy_unit(char const *which, char const **const symbol); // declaration only

    double length_unit(char const *which, char const **const symbol); // declaration only

    status_t set(char const *length, char const *energy, int const echo=0); // declaration only

    status_t all_tests(int const echo=0); // declaration only

} // namespace unit_system
