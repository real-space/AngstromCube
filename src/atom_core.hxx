#pragma once

#include "radial_grid.h" // radial_grid_t

typedef int status_t;

namespace atom_core {

  status_t scf_atom(
      radial_grid_t const &g, // radial grid descriptor
      float const Z, // atomic number
      int const echo); // log output level
  
  status_t all_tests();

} // namespace atom_core
