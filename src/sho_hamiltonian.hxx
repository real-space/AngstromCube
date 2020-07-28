#pragma once

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>
#include "real_space.hxx" // ::grid_t

namespace sho_hamiltonian {

  status_t solve(int const natoms
          , view2D<double const> const & xyzZ // (natoms, 4)
          , real_space::grid_t const & g
          , double const *const vtot
          , int const natoms_prj=-1
          , double const *const sigma_prj=nullptr
          , int    const *const numax_prj=nullptr
          , double *const *const atom_mat=nullptr
          , int const echo=0);
  
  status_t all_tests(int const echo=0);

} // namespace sho_hamiltonian
