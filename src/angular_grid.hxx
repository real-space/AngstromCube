#pragma once

#include "angular_grid.h" // angular_grid_t

typedef int status_t;

namespace angular_grid {

  angular_grid_t* get_grid(int const ellmax, int const echo=0);
  
  template<typename real_t>
  void transform(real_t *out, real_t const *in, int const nrad, int const ellmax, bool const back=false, int echo=0);

  int Lebedev_grid_size(int const ellmax, int echo=0);

  template<typename real_t>
  status_t create_Lebedev_grid(int const ellmax, real_t xyzw[][4], int echo=9);

  status_t all_tests();
  
} // namespace angular_grid
