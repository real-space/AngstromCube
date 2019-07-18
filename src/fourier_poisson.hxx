#pragma once

typedef int status_t;

namespace fourier_poisson {

  status_t fourier_poisson(double x[], double const b[], int const ng[3], double const bravais[3][4]);

  int all_tests();
  
} // namespace fourier_poisson
