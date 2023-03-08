#pragma once
// This file is part of AngstromCube under MIT License

#include <vector> // std::vector<T>

#include "status.hxx" // status_t

namespace green_input {

  status_t load_Hamiltonian(
        uint32_t ng[3] // numbers of grid points
      , int8_t bc[3] // boundary conditions
      , double hg[3] // grid spacings
      , std::vector<double> & Veff
      , int & natoms
      , std::vector<double> & xyzZinso
      , std::vector<std::vector<double>> & atom_mat
      , char const *const filename="Hmt.json" // input
      , int const echo=0 // log-level
  ); // declaration only

  status_t all_tests(int echo=0); // declaration only

} // namespace green_input
