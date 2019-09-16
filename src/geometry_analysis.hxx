#pragma once

typedef int status_t;

namespace geometry_analysis {

  status_t read_xyz_file(double **xyzz, int *n_atoms, char const *filename, 
                         double *cell=nullptr, int *bc=nullptr, int const echo=5);

  status_t all_tests();

} // namespace geometry_analysis
