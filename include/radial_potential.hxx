#pragma once
// This file is part of AngstromCube under MIT License

#include "radial_grid.h" // radial_grid_t
#include "quantum_numbers.h" // ell_QN_t

#include "status.hxx" // status_t

namespace radial_potential {

  double Hartree_potential( // returns the Coulomb integral
        double rV[] // result: r*Hartree-potential(r)
      , radial_grid_t const & g // radial grid descriptor
      , double const rho4pi[] // 4*pi*density(r)
  ); // declaration only

  void Hartree_potential(
        double vHt[] // result: Hartree-potential_lm(r)
      , radial_grid_t const & g // radial grid descriptor
      , double const rho[] // density_lm(r)
      , int const stride // stride between different lm-compenents in rho and V
      , ell_QN_t const ellmax // largest angular momentum quantum number
      , double const q0=0 // singularity
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace radial_potential
