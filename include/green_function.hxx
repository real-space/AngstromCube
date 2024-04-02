#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdint> // uint32_t, int8_t
#include <vector>  // std::vector<T>
#include <complex> // std::complex<real_t>

#include "status.hxx" // status_t
#include "action_plan.hxx" // action_plan_t
#include "green_parallel.hxx" // ::RequestList_t
#include "data_view.hxx" // view3D<T>

 /*
  *  Future plan:
  *   Support also density matrix purification scheme (McWeeney filter: x^2(3-2x)
  *   or with a non-trivial overlap operator S (3xSx - 2xSxSx from Phys. Rev. B 50, 17622)
  *   maybe better with norm-conserving PAW formalism --> S == 1
  *   in particular suitable for the real Hamiltonian action (which is not supported by tfQMRgpu)
  */

namespace green_function {

    status_t construct_Green_function(
          action_plan_t & p // result, create a plan how to apply the SHO-PAW Hamiltonian to a block-sparse truncated Green function
        , uint32_t const ng[3] // numbers of grid points of the unit cell in with the potential is defined
        , int8_t const boundary_condition[3] // boundary conditions
        , double const hg[3] // grid spacings
        , std::vector<double> const & xyzZinso // [natoms*8]
        , int const echo=0 // verbosity
        , int const Noco=2 // 1:collinear spins, 2:Non-collinear
    ); // declaration only

    status_t update_potential(
          action_plan_t & p // modify
        , uint32_t const nb[3] // numbers of 4*4*4 grid blocks of the unit cell in with the potential is defined
        , std::vector<double> const & Veff // [nb[2]*4*nb[1]*4*nb[0]*4]
        , view3D<uint16_t> const & owner_rank // [nb[2],nb[1],nb[0]]
        , std::vector<std::vector<double>> const & AtomMatrices
        , int const echo=0 // verbosity
        , int const Noco=2 // 1:collinear spins, 2:Non-collinear
    ); // declaration only

    status_t update_energy_parameter(
          action_plan_t & p // modify
        , std::complex<double> E_param
     // , std::vector<std::vector<double>> const & AtomMatrices
        , double const dVol // volume element of the grid
        , int const echo=0 // verbosity
        , int const Noco=1 // 1:collinear spins, 2:Non-collinear
        , double const scale_H=1
     // , green_parallel::RequestList_t const *requests=nullptr
    ); // declaration only

    status_t update_phases(
          action_plan_t & p // modify
        , double const k_point[3]
        , int const echo=0 // verbosity
        , int const Noco=1 // // 1:collinear spins, 2:Non-collinear
    ); // declaration only

    status_t all_tests(int const echo=0); // declaration only

} // namespace green_function
