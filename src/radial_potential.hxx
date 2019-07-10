#pragma once

#include "radial_grid.h" // radial_grid_t
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t

typedef int status_t;

namespace radial_potential {
  
  double Hartree_potential(   // returns the Coulomb integral
      double rV[],            // r*Hartree-potential(r)
      radial_grid_t const &g, // radial grid descriptor
      double const rho4pi[],  // 4*pi*density(r)
      ell_QN_t const ell=0);  // [optional] angular momentum channel

  status_t all_tests();

} // namespace radial_potential
