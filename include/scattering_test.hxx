#pragma once
// This file is part of AngstromCube under MIT License

#include "radial_grid.h" // radial_grid_t
#include "energy_level.hxx" // TRU_AND_SMT
#include "status.hxx" // status_t

namespace scattering_test {

  status_t expand_sho_projectors( // pure SHO basis functions
        double prj[]  // output projectors [nln*stride]
      , int const stride // stride >= rg.n
      , radial_grid_t const & rg // radial grid descriptor
      , double const sigma // SHO basis spread
      , int const numax // SHO basis size
      , int const rpow=0 // return r^rpow * p(r)
      , int const echo=0 // log-level
      , double dprj[]=nullptr // derivative of projectors w.r.t. sigma
  ); // declaration only

  status_t logarithmic_derivative(
        radial_grid_t const rg[TRU_AND_SMT] // radial grid descriptors for Vtru, Vsmt
      , double const *const rV[TRU_AND_SMT] // true and smooth potential given on the radial grid *r
      , double const sigma // sigma spread of SHO projectors
      , int const ellmax // ellmax up to which the analysis should go
      , int const numax
      , double const aHm[] // non-local Hamiltonian elements in ln_basis
      , double const aSm[] // non-local overlap matrix elements
      , double const energy_range[3] // {lower, step, upper}
      , char const *label="" // log-prefix
      , int const echo=0 // log-level
      , float const Rlog_over_sigma=6.f
  ); // declaration only

  status_t eigenstate_analysis(
        radial_grid_t const & gV // radial grid descriptor for Vsmt
      , double const Vsmt[] // smooth potential given on radial grid
      , double const sigma // sigma spread of SHO projectors
      , int const ellmax // ellmax up to which the analysis should go
      , int const numax // SHO basis size
      , double const aHm[] // non-local Hamiltonian elements in ln_basis, assume stride nln
      , double const aSm[] // non-local overlap matrix elements, assume stride nln
      , int const nr=384 // number of radial grid points in equidistance mesh
      , double const Vshift=0 // potential shift
      , char const *label="" // log-prefix
      , int const echo=0 // log-level
      , float const reference[3][4]=nullptr // up to 3 reference eigenvalues for s,p,d,f
      , float const warning_threshold=3e-3
  ); // declaration only

  status_t emm_average(
        double Mln[] // output emm-averaged or emm-summed array[nln*nln]
      , double const Mlmn[] // input array[nlmn*nlmn]
      , int const numax // SHO basis size
      , int const avg2sum0=2 // 2:average over emm-values withing each ell-channel, 0:sum only
      , int const echo=0 // log-level
      , int const stride=-1 // optional stride for Mlmn, defaults to nlmn
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace scattering_test
