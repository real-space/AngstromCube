#pragma once

#include "radial_grid.h" // radial_grid_t

typedef int status_t;

namespace atom_core {

  int scf_atom(radial_grid_t const &g, float const Z);
  
  status_t all_tests();

} // namespace atom_core
