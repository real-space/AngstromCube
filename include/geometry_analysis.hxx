#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdint> // int8_t

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>

  inline double length_squared(double x, double y, double z) { return x*x + y*y + z*z; }
  inline double length(double x, double y, double z) { return std::sqrt(length_squared(x, y, z)); }
  inline double length(double const xyz[3]) { return length(xyz[0], xyz[1], xyz[2]); }

namespace geometry_analysis {

  status_t read_xyz_file(
        view2D<double> & xyzZ // result: atom info, x,y,z,Z,...
      , int & n_atoms // result: number of atoms
      , double cell[3][4] // result: cell shape
      , int8_t bc[3]=nullptr // result: boundary conditions
      , char const *filename="atoms.xyz" // filename
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
