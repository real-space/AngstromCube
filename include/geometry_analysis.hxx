#pragma once

#include <cstdint> // int8_t

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>

namespace geometry_analysis {

  status_t read_xyz_file(
        view2D<double> & xyzZ // result: atom info, x,y,z,Z,...
      , int & n_atoms // result: number of atoms
      , char const *filename="atoms.xyz" // filename
      , double cell[]=nullptr // result: cell extent
      , int8_t bc[]=nullptr // result: boundary conditions
      , int const echo=5 // log level
  ); // declaration only

  inline double fold_back(double const position, double const cell_extend) { 
      double x{position};
      while (x >= 0.5*cell_extend) x -= cell_extend;
      while (x < -0.5*cell_extend) x += cell_extend;
      return x;
  } // fold_back

  status_t all_tests(int const echo=0); // declaration only

} // namespace geometry_analysis
