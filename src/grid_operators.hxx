#pragma once

typedef int status_t;

#include "grid_operators.hxx"

#include "atom_image.hxx" // atom_image_t, sho_atom_t
#include "real_space_grid.hxx" // grid_t
#include "finite_difference.hxx" // ::finite_difference_t<real_t>

namespace grid_operators {

  // Hamiltonian
  template<typename real_t, int D0>
  status_t grid_Hamiltonian(real_t Hpsi[] // result
                          , real_t const psi[] // input wave functions
                          , real_space_grid::grid_t<D0> const &g // 3D Cartesian grid descriptor
                          , std::vector<atom_image::sho_atom_t> const &a
                          , std::vector<atom_image::atom_image_t> const &ai
                          , finite_difference::finite_difference_t<real_t> const &fd // finite difference
                          , double const potential[] // diagonal potential operator
                          , double const *boundary_phase=nullptr); // phase shifts at the boundary [optional]

  // Overlap operator
  template<typename real_t, int D0>
  status_t grid_Overlapping(real_t Spsi[] // result
                          , real_t const psi[] // input wave functions
                          , real_space_grid::grid_t<D0> const &g // 3D Cartesian grid descriptor
                          , std::vector<atom_image::sho_atom_t> const &a
                          , std::vector<atom_image::atom_image_t> const &ai
                          , double const *boundary_phase=nullptr);

  status_t all_tests(int const echo=0);

} // namespace grid_operators
