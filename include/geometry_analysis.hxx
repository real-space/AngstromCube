#pragma once

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>

namespace geometry_analysis {

  status_t read_xyz_file(
        view2D<double> & xyzZ
      , int & n_atoms
      , char const *filename="atoms.xyz"
      , double cell[]=nullptr
      , int bc[]=nullptr
      , int const echo=5
  );

  inline double fold_back(double const position, double const cell_extend) { 
      double x{position};
      while (x >= 0.5*cell_extend) x -= cell_extend;
      while (x < -0.5*cell_extend) x += cell_extend;
      return x;
  } // fold_back

  status_t all_tests(int const echo=0); // declaration only

} // namespace geometry_analysis