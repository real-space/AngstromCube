#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdint> // int32_t

#include "status.hxx" // status_t
#include "real_space.hxx" // ::grid_t
#include "data_view.hxx" // view2D

namespace geometry_input {

    status_t read_xyz_file(
          view2D<double> & xyzZ
        , int32_t & n_atoms
        , double cell[3][4]
        , int8_t bc[3]=nullptr
        , char const *const filename="atoms.xyz"
        , int const echo=5 // log-level
    ); // declaration only

    status_t init_geometry_and_grid(
          real_space::grid_t & g // output grid descriptor
        , view2D<double> & xyzZ // output atom coordinates and core charges Z
        , int32_t & natoms // output number of atoms found
        , unsigned const n_even=2 // make sure the number of grid points can be divided by n_even
        , int const echo=0 // log-level
    ); // declaration only

    double get_temperature(int const echo, double const def=1e-3); // declaration only

    status_t all_tests(int const echo=0); // declaration only

} // namespace geometry_input
