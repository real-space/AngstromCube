#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdint>    // int64_t, int32_t, uint32_t, int16_t, uint16_t, int8_t, uint8_t
#include <vector>     // std::vector<T>
#include <complex>    // std::complex

#include "status.hxx" // status_t
#include "green_action.hxx" // ::plan_t
#include "green_dyadic.hxx" // ::dyadic_plan_t


 /*
  *  Future plan:
  *   Support also density matrix purification scheme (McWeeney filter: x^2(3-2x)
  *   or with a non-trivial overlap operator S (3xSx - 2xSxSx from Phys. Rev. B 50, 17622)
  *   maybe better with norm-conserving PAW formalism --> S == 1
  *   in particular suitable for the real Hamiltonian action (which is not supported by tfQMRgpu)
  */

namespace green_function {

  char const boundary_condition_name[][16] = {"isolated", "periodic", "vacuum"};
  char const boundary_condition_shortname[][8] = {"iso", "peri", "vacu"};

  status_t construct_Green_function(
        green_action::plan_t & p // result, create a plan how to apply the SHO-PAW Hamiltonian to a block-sparse truncated Green function
      , uint32_t const ng[3] // numbers of grid points of the unit cell in with the potential is defined
      , int8_t const boundary_condition[3] // boundary conditions
      , double const hg[3] // grid spacings
      , std::vector<double> const & Veff // [ng[2]*ng[1]*ng[0]]
      , std::vector<double> const & xyzZinso // [natoms*8]
      , std::vector<std::vector<double>> const & AtomMatrices // atomic hamiltonian and overlap matrix, [natoms][2*nsho^2]
      , int const echo=0 // log-level
      , std::complex<double> const *energy_parameter=nullptr // E in G = (H - E*S)^{-1}
      , int const Noco=2
  ); // declaration only

  status_t update_atom_matrices(
        green_dyadic::dyadic_plan_t & p
      , std::complex<double> E_param
      , std::vector<std::vector<double>> const & AtomMatrices
      , double const dVol // volume element of the grid
      , int const Noco // =1
      , double const scale_H // =1
      , int const echo // =0
  ); // declaration only

  status_t update_phases(
        green_action::plan_t & p
      , double const k_point[3]
      , int const Noco=1
      , int const echo=0 // verbosity
  ); // declaration only

  inline status_t update_energy_parameter(
        green_action::plan_t & p
      , std::complex<double> E_param
      , std::vector<std::vector<double>> const & AtomMatrices
      , double const dVol
      , int const Noco=1
      , double const scale_H=1
      , int const echo=0
  ) {
      p.E_param = E_param;
      return update_atom_matrices(p.dyadic_plan, E_param, AtomMatrices, dVol, Noco, scale_H, echo);
  } // update_energy_parameter

  status_t all_tests(int const echo=0); // declaration only

} // namespace green_function
