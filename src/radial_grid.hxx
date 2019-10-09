#pragma once

#include <cmath> // std::sqrt

#include "radial_grid.h" // radial_grid_t

typedef int status_t;

namespace radial_grid {
  
  radial_grid_t* create_exponential_radial_grid( // returns a pointer to a new radial grid descriptor
      int const npoints, // number of grid points
      float const Rmax=9.45, // [optional] largest radius
      float const anisotropy=.015); // [optional] anisotropy parameter

  radial_grid_t* create_equidistant_radial_grid(int const npoints, float const Rmax=9.45);
  
  inline radial_grid_t* create_default_radial_grid(float const Z_nucleons)
    { return create_exponential_radial_grid(250*std::sqrt(Z_nucleons + 9.)+.5); }

  radial_grid_t* create_pseudo_radial_grid(radial_grid_t const &tru, double const r_min=1e-3, int const echo=0);

  void destroy_radial_grid(radial_grid_t* g); // radial grid descriptor

  status_t all_tests();
  
} // namespace radial_grid
