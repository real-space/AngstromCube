#pragma once

#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
#include "radial_grid.h" // radial_grid_t

#include "status.hxx" // status_t

namespace radial_eigensolver {

  status_t shooting_method(
        int const sra // 1:scalar relativistic approximation, 0:Schroedinger equation
      , radial_grid_t const & g // radial grid descriptor
      , double const rV[] // radial potential r*V_Hxc(r) - e^2*Z
      , enn_QN_t const enn // energy quantum number
      , ell_QN_t const ell // angular momentum quantum number
      , double & E // energy (eigen-)value in Hartree
      , double* rf=nullptr // radial wave function*r
      , double* r2rho=nullptr // density of that wave function*r^2
      , int const maxiter=999 // maximum number of iterations
      , float const threshold=1e-15 // threshold for eigenvalue convergence
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace radial_eigensolver
