#pragma once

#include <vector> // std::vector<T>

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>

namespace sho_basis {

  status_t load(
        std::vector<view2D<double>> & basis
      , std::vector<int> & indirection
      , int const natoms // number of SHO basis centers
      , double const Z_core[] // Z_core[atoms]
      , int const echo=0 // log-level
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace sho_basis
