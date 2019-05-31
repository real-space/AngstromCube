#pragma once

#include "radial_grid.h" // radial_grid_t

namespace radial_grid {

  radial_grid_t* create_exponential_radial_grid( // returns a pointer to a new radial grid descriptor
      int const n, // number of grid points
      float const r=9.45, // [optional] largest radius
      float const a=.015); // [optional] anisotropy parameter

  radial_grid_t* create_pseudo_radial_grid(radial_grid_t const &tru, double const r_min=.001);

  void destroy_radial_grid(radial_grid_t* g); // radial grid descriptor

  int all_tests();
  
} // namespace radial_grid
