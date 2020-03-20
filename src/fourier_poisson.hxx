#pragma once

#include "status.hxx" // status_t

#include "constants.hxx" // pi

namespace fourier_poisson {

  double constexpr epsilon0 = 4*constants::pi; // in atomic units
  
  template<typename real_t>
  status_t fourier_solve(real_t x[], real_t const b[], int const ng[3], double const bravais[3][4],
                         double const factor=-epsilon0, int const echo=0);

  status_t all_tests(int const echo=0);
  
} // namespace fourier_poisson
