#pragma once
// This file is part of AngstromCube under MIT License

#include "radial_grid.h" // radial_grid_t
#include "quantum_numbers.h" // ell_QN_t

#include "status.hxx" // status_t

namespace radial_integrator {

  double shoot( // returns the kink
        int const sra // 0:Schroedinger equation, 1:scalar relativistic approximation, 2:appromimated mass
      , radial_grid_t const & g // radial grid descriptor
      , double const rV[] // radial potential r*V_Hxc(r) - e^2*Z
      , ell_QN_t const ell // angular momentum quantum number
      , double const E // energy (eigen-)value in Hartree
      , int & nnodes // inout: number of nodes
      , double* rf=nullptr // optional result: radial wave function*r
      , double* r2rho=nullptr // optional result: density of that wave function*r^2
  ); // declaration only

  template <int SRA> // 1:scalar relativistic approximation, 0:Schroedinger equation
  int integrate_inwards( // return the number of nodes
        double gg[] // result: large component*r
      , double ff[] // result: small component*r
      , radial_grid_t const & g // radial grid descriptor
      , double const rV[] // radial potential r*V_Hxc(r) - e^2*Z
      , ell_QN_t const ell // angular momentum quantum number
      , double const E // energy (eigen-)value in Hartree
      , double const valder[2]=nullptr // optional: value and derivative
      , double *dg=nullptr // optional result: derivative at end point
      , int *ir_stopped=nullptr // optional result: index at which the inwards integration stopped
      , int const ir_start=-1 // start index, -1:(g.n - 1)
      , int const ir_stop=4 // latest stop index
  ); // declaration only

  template <int SRA> // 1:scalar relativistic approximation, 0:Schroedinger equation
  int integrate_outwards( // return the number of nodes
        double gg[] // result: large component*r
      , double ff[] // result: small component*r
      , radial_grid_t const & g // radial grid descriptor
      , double const rV[] // radial potential r*V_Hxc(r) - e^2*Z
      , ell_QN_t const ell // angular momentum quantum number
      , double const E // energy (eigen-)value in Hartree
      , int const ir_stop=-1 // latest stop index, -1:(g.n - 1)
      , double *dg=nullptr // optional result: derivative at end point
      , double const *rp=nullptr // optional input: inhomogeneity*r
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace radial_integrator
