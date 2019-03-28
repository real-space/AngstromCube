#pragma once

#include "radial_grid.h" // radial_grid_t
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t

typedef int status_t;

namespace radial_potential {
  
  double Hartree_potential(double vH[], radial_grid_t const &g, double const rho4pi[], ell_QN_t const ell=0);

  double lda_PZ81_kernel(double const rho, double &Vup, double const mag=0, double *Vdn=nullptr);
  
  status_t all_tests();

} // namespace radial_potential
