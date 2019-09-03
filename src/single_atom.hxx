#pragma once

#include "radial_grid.h" // radial_grid_t

typedef int status_t;

namespace single_atom {
  
  status_t update(float const Za[], int const na, 
                  double **rho=nullptr, radial_grid_t **rg=nullptr, double *sigma_cmp=nullptr);

  status_t all_tests();

} // namespace single_atom
