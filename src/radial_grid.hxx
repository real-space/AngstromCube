#pragma once

#include "radial_grid.h" // radial_grid_t

namespace radial_grid {

  radial_grid_t* create_exponential_radial_grid( // returns a pointer to a new radial grid descriptor
      int const n, // number of grid points
      float const r=9.45, // [optional] largest radius
      float const a=.015); // [optional] anisotropy parameter

  radial_grid_t* create_pseudo_radial_grid(radial_grid_t const &tru, double const r_min=1e-3);

  void destroy_radial_grid(radial_grid_t* g); // radial grid descriptor

  template<typename real_t>
  real_t inline dot_product(int const n, real_t const bra[], real_t const ket[]) {
      real_t dot = 0;
      for(int i = 0; i < n; ++i) {
          dot += bra[i]*ket[i];
      } // i
      return dot;
  } // dot_product
  
  int all_tests();
  
} // namespace radial_grid
