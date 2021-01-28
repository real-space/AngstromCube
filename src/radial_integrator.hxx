#pragma once

#include "radial_grid.h" // radial_grid_t
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t

#include "status.hxx" // status_t

namespace radial_integrator {

  double shoot( // returns the kink
      int const sra, // 1:scalar relativistic approximation, 0:Schroedinger equation
      radial_grid_t const &g, // radial grid descriptor
      double const rV[], // radial potential r*V_Hxc(r) - e^2*Z
      ell_QN_t const ell, // angular momentum quantum number
      double const E, // energy (eigen-)value in Hartree
      int &nnodes, // number of nodes
      double* rf=nullptr, // radial wave function*r
      double* r2rho=nullptr); // density of that wave function*r^2

  template <int SRA> // 1:scalar relativistic approximation, 0:Schroedinger equation
  int integrate_inwards( // return the number of nodes
      radial_grid_t const &g, // radial grid descriptor
      double const rV[], // radial potential r*V_Hxc(r) - e^2*Z
      ell_QN_t const ell, // angular momentum quantum number
      double const E, // energy (eigen-)value in Hartree
      double gg[], // large component*r
      double ff[], // small component*r
      double const valder[2]=nullptr,
      double *dg=nullptr, // derivative at end point
      int *ir_stopped=nullptr, // index at which the inwards integration stopped
      int const ir_start=-1, // start index, -1:(g.n - 1)
      int const ir_stop=4); // latest stop index

  template <int SRA> // 1:scalar relativistic approximation, 0:Schroedinger equation
  int integrate_outwards( // return the number of nodes
      radial_grid_t const &g, // radial grid descriptor
      double const rV[], // radial potential r*V_Hxc(r) - e^2*Z
      ell_QN_t const ell, // angular momentum quantum number
      double const E, // energy (eigen-)value in Hartree
      double gg[], // large component*r
      double ff[], // small component*r
      int const ir_stop=-1, // latest stop index, -1:(g.n - 1)
      double *dg=nullptr, // derivative at end point
      double const *rp=nullptr); // inhomogeneity*r, only outward

  status_t all_tests(int const echo=0); // declaration only

} // namespace radial_integrator
