#pragma once

#include <cstdint> // int8_t

#include "status.hxx" // status_t

namespace sigma_config {

    typedef struct {
        double  Z;
        double  rcut;
        double  sigma;
        double  q_core_hole;
        uint8_t nn[8];
        int8_t  ncmx[4];
        int8_t  iln_core_hole;
    } element_t;

    element_t & get(double const Zcore, int const echo=0);

    status_t all_tests(int const echo=0);

} // namespace sigma_config
