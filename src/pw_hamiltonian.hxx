#pragma once

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>
#include "real_space.hxx" // ::grid_t

namespace pw_hamiltonian {

  status_t solve(int const natoms_prj // number of PAW atoms
          , view2D<double const> const & xyzZ // (natoms, 4)
          , real_space::grid_t const & g // Cartesian grid descriptor for vtot
          , double const *const vtot // total effective potential on grid
          , double const *const sigma_prj // =nullptr
          , int    const *const numax_prj // =nullptr
          , double *const *const atom_mat=nullptr // PAW charge-deficit and Hamiltonian correction matrices
          , int const echo=0); // log-level

  status_t all_tests(int const echo=0);

} // namespace pw_hamiltonian
