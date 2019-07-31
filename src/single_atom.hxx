#pragma once

typedef int status_t;

namespace single_atom {
  
  status_t update(float const Za[], int const na, double **rho=nullptr);
  
  status_t all_tests();

} // namespace single_atom
