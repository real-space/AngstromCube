#pragma once

#include <vector> // std::vector

#include "angular_grid.h" // angular_grid_t
#include "gaunt_entry.h" // gaunt_entry_t

#include "status.hxx" // status_t

namespace angular_grid {

  angular_grid_t* get_grid(int const ellmax, int const echo=0);

// status_t get_weights(double weights[], int const ellmax, int const echo=0); // needed?

  template<typename real_t>
  status_t transform(real_t out[], real_t const in[], int const stride, int const ellmax, bool const back=false, int const echo=0);

  int Lebedev_grid_size(int const ellmax, int const echo=0);

  template<typename real_t>
  int create_Lebedev_grid(int const ellmax, real_t xyzw[][4], int const echo=0);

  template<int lmax>
  status_t create_numerical_Gaunt(std::vector<gaunt_entry_t> & gaunt, int const echo=0);

  status_t all_tests(int const echo=0); // declaration only

} // namespace angular_grid
