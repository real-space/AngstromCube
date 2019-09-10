#pragma once

#include "radial_grid.h" // radial_grid_t

typedef int status_t;

namespace single_atom {

  status_t update(int const na, float const Za[], float const ion[], 
                  radial_grid_t **rg=nullptr, double *sigma_cmp=nullptr, 
                  double **rho=nullptr, double **qlm=nullptr, double **vlm=nullptr);

  status_t all_tests();

} // namespace single_atom
