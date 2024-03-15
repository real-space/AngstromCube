#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdint>    // int64_t, int32_t, uint32_t, int16_t, uint16_t, int8_t, uint8_t
#include <vector>     // std::vector<T>
#include <complex>    // std::complex

#include "status.hxx" // status_t
// #include "green_action.hxx" // ::plan_t
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

//     char const boundary_condition_name[][16] = {"isolated", "periodic", "vacuum"};
//     char const boundary_condition_shortname[][8] = {"iso", "peri", "vacu"};

    status_t construct_Green_function(
          action_plan_t & p // result, create a plan how to apply the SHO-PAW Hamiltonian to a block-sparse truncated Green function
        , uint32_t const ng[3] // numbers of grid points of the unit cell in with the potential is defined
        , int8_t const boundary_condition[3] // boundary conditions
        , double const hg[3] // grid spacings
  //    , std::vector<double> const & Veff // [ng[2]*ng[1]*ng[0]]
        , std::vector<double> const & xyzZinso // [natoms*8]
  //    , std::vector<std::vector<double>> const & AtomMatrices // atomic hamiltonian and overlap matrix, [natoms][2*nsho^2]
        , int const echo=0 // verbosity
     // , std::complex<double> const *energy_parameter=nullptr // E in G = (H - E*S)^{-1}
        , int const Noco=2 // 1:collinear spins, 2:Non-collinear
    ); // declaration only

    status_t update_potential(
          action_plan_t & p // modify
        , uint32_t const nb[3] // numbers of 4*4*4 grid blocks of the unit cell in with the potential is defined
        , std::vector<double> const & Veff // [nb[2]*4*nb[1]*4*nb[0]*4]
        , view3D<uint16_t> const & owner_rank // [nb[2],nb[1],nb[0]]
        , int const echo=0 // verbosity
        , int const Noco=2 // 1:collinear spins, 2:Non-collinear
    ); // declaration only

    status_t update_energy_parameter(
          action_plan_t & p // modify
        , std::complex<double> E_param
        , std::vector<std::vector<double>> const & AtomMatrices
        , double const dVol // volume element of the grid
        , double const scale_H // =1
        , int const echo=0 // verbosity
        , int const Noco=1 // 1:collinear spins, 2:Non-collinear
        , green_parallel::RequestList_t const *requests=nullptr
    ); // declaration only

    status_t update_phases(
          action_plan_t & p // modify
        , double const k_point[3]
        , int const echo=0 // verbosity
        , int const Noco=1 // // 1:collinear spins, 2:Non-collinear
    ); // declaration only

    status_t all_tests(int const echo=0); // declaration only

} // namespace green_function
