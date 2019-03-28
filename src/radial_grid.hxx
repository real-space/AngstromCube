#pragma once

#include "radial_grid.h" // radial_grid_t

namespace radial_grid {

  radial_grid_t* create_exponential_radial_grid(int n, float r=9.45, float a=.015);

  void destroy_radial_grid(radial_grid_t* g);

  int all_tests();
  
} // namespace radial_grid
