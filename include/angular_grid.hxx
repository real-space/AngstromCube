#pragma once

#include <vector> // std::vector

#include "status.hxx" // status_t
#include "gaunt_entry.h" // gaunt_entry_t

namespace angular_grid {

  int get_grid_size(int const ellmax, int const echo=0); // declaration only

  template <typename real_t>
  status_t transform(real_t out[], real_t const in[], int const stride, int const ellmax, bool const back=false, int const echo=0);

  template <int lmax>
  std::vector<gaunt_entry_t> create_numerical_Gaunt(int const echo=0);

  void cleanup(int const echo=0); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace angular_grid
