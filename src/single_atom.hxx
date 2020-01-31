#pragma once

#include "radial_grid.h" // radial_grid_t

typedef int status_t;

namespace single_atom {

  status_t update(int const na, float const Za[], float const ion[],
                  radial_grid_t **rg=nullptr, double *sigma_cmp=nullptr,
                  double **rho=nullptr, double **qlm=nullptr, double **vlm=nullptr,
                  int *lmax=nullptr, int *lmax_cmp=nullptr, double **zero_pot=nullptr);

  status_t all_tests(int const echo=0);

} // namespace single_atom
