#pragma once

#include <cstdint> // int8_t, uint8_t

#include "status.hxx" // status_t

namespace sigma_config {

  typedef struct {
      double  occ[32][2]; // occupation numbers up to nl-index 31 (8f-orbital)
      int8_t  csv[32]; // 0:core, 1:semicore, 2:valence, see spherical_state.hxx
      double  Z; // number of protons in the core
      double  rcut; // cutoff-radius
      double  sigma; // Gaussian spread parameter for SHO-type projectors
      uint8_t nn[8]; // number of radial projectors per ell for ell=0..7
      char    method[15]; // method for the local potential
      int8_t  numax; // user specified numax, -1:auto
  } element_t;

  element_t const & get(
        double const Zcore // nuclear charge
      , int const echo=0 // log-level
      , char const **configuration=nullptr // export the string that has been parsed
  ); // declaration only

  void set_default_core_shells(int ncmx[4], double const Z);

  status_t all_tests(int const echo=0); // declaration only

} // namespace sigma_config
