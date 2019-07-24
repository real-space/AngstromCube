#pragma once

#include "angular_grid.h" // angular_grid_t
#include "gaunt_entry.h" // gaunt_entry_t
#include <vector> // std::vector

typedef int status_t;

namespace angular_grid {

  angular_grid_t* get_grid(int const ellmax, int const echo=0);

  template<typename real_t>
  void transform(real_t *out, real_t const *in, int const nrad, int const ellmax, bool const back=false, int const echo=0);

  int Lebedev_grid_size(int const ellmax, int echo=0);

  template<typename real_t>
  status_t create_Lebedev_grid(int const ellmax, real_t xyzw[][4], int echo=9);

  status_t all_tests();

  template<int lmax>
  status_t create_numerical_Gaunt(std::vector<gaunt_entry_t>* gaunt=nullptr, int const echo=0);

} // namespace angular_grid
