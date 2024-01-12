#pragma once
// This file is part of AngstromCube under MIT License

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>

namespace sho_basis {

  status_t get( // check if radial functions can be loaded from pseudo_basis.xml
        double & sigma // result: SHO spread
      , int & numax // result: basis size, input: -1 for minimum given in file
      , int & nbasis // result: size of the loaded basis (including emm-multiplicity)
      , double const Z // input: number of protons in the core
      , int const echo=0 // log-level
  ); // declaration only

  template <typename complex_t>
  status_t generate(
        view2D<complex_t> & matrix // result: reduction matrix[nsho,nbasis]
      , double & sigma // result: SHO spread
      , int & numax // input/output SHO basis size, -1=smallest possible
      , double const Z_core // input core charge
      , int const echo=0 // log-level verbosity
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace sho_basis
