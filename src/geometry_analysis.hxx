#pragma once

#include "status.hxx" // status_t

namespace geometry_analysis {

  status_t read_xyz_file(double **xyzZ, int *n_atoms, char const *filename, 
                         double *cell=nullptr, int *bc=nullptr, int const echo=5);

  status_t all_tests(int const echo=0);

} // namespace geometry_analysis
