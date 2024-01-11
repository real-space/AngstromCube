#pragma once
// This file is part of AngstromCube under MIT License

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>

namespace sho_basis {

  template <typename complex_t>
  status_t generate(
        view2D<complex_t> & matrix // result: reduction matrix[nsho,nbasis]
      , double & sigma // result: SHO spread
      , double const Z // input core charge
      , int const numax=-1 // input SHO basis size, -1=smallest possible
      , int const echo=0 // log-level verbosity
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace sho_basis
