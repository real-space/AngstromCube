#pragma once

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>
#include "real_space.hxx" // ::grid_t

namespace sho_hamiltonian {

  status_t solve(
        int const natoms // number of SHO basis centers
      , view2D<double> const & xyzZ // (natoms, 4)
      , real_space::grid_t const & g // Cartesian grid descriptor for vtot
      , double const *const vtot // total effective potential on grid
      , int const nkpoints
      , view2D<double> const & kmesh
      , int const natoms_prj=-1 // number of PAW centers
      , double const *const sigma_prj=nullptr // spreads of the SHO-type PAW projectors
      , int    const *const numax_prj=nullptr // cutoffs of the SHO-type PAW projectors
      , double *const *const atom_mat=nullptr // PAW charge-deficit and Hamiltonian correction matrices
      , int const echo=0 // log-level
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace sho_hamiltonian
