#pragma once

#include <cstdint> // int8_t, uint8_t
#include <cassert> // assert

#include "status.hxx" // status_t

namespace sigma_config {

  typedef struct {
      double  occ[32][2]; // occupation numbers up to nl-index 31 (8f-orbital)
      int8_t  csv[32]; // 0:core, 1:semicore, 2:valence
      double  Z; // number of protons in the core
      double  rcut; // cutoff-radius
      double  sigma; // Gaussian spread parameter for SHO-type projectors
      double  q_core_hole[2]; // charge of the core hole, down and up-spin
      char    method[16]; // method for the local potential
      uint8_t nn[8]; // number of radial projectors per ell for ell=0..7
      int8_t  ncmx[4]; // highest enn quantum number of core states for ell=0..3
      int8_t  numax; // user specified numax, -1:auto
      int8_t  inl_core_hole; // state index where to insert the core hole
  } element_t;

  element_t const & get(
        double const Zcore // nuclear charge
      , int const echo=0 // log-level
      , char const **configuration=nullptr // string that has been parsed
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace sigma_config
