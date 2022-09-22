#pragma once

#include <vector> // std::vector

#include "status.hxx" // status_t
#include "gaunt_entry.h" // gaunt_entry_t

namespace angular_grid {

  int get_grid_size(int const ellmax, int const echo=0); // declaration only

  status_t transform(
        double out[] // resulting matrix with stride
      , double const in[] // input matrix with stride
      , int const stride // stride for in and out
      , int const ellmax // largest angular momentum
      , bool const back=false // backtransform?
      , int const echo=0 // log-level
  ); // declaration only

  std::vector<gaunt_entry_t> create_numerical_Gaunt(int const ellmax, int const echo=0); // declaration only

  void cleanup(int const echo=0); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace angular_grid
