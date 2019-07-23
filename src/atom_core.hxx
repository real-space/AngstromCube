#pragma once

#include "radial_grid.h" // radial_grid_t

typedef int status_t;

namespace atom_core {

  status_t scf_atom(
      radial_grid_t const &g, // radial grid descriptor
      float const Z, // atomic number
      int const echo); // log output level

  double initial_density(double r2rho[], radial_grid_t const &g, double const Z, double const charged=0);

  void rad_pot(double rV[], radial_grid_t const &g, double const rho4pi[], double const Z=0, double *energies=nullptr);
  
//   double dot_product(int const n, double const bra[], double const ket[]); // moved to radial_grid.hxx
  
  status_t all_tests();

} // namespace atom_core
