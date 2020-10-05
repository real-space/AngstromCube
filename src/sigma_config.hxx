#pragma once

#include <cstdint> // int8_t, uint8_t

#include "status.hxx" // status_t

namespace sigma_config {

    typedef struct {
        double  occ[32][2]; // occupation numbers up nl-index 31 (8f-orbital)
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

    element_t & get(double const Zcore, int const echo=0);

    inline int nl_index(int const enn, int const ell) { assert(ell < enn); return (enn*(enn - 1))/2 + ell; }
    
    status_t all_tests(int const echo=0);

} // namespace sigma_config
