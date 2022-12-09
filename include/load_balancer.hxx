#pragma once

#include <cstdint> // uint32_t, uint16_t

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

namespace load_balancer {

  double get(
        uint32_t const comm_size // number of MPI processes in this communicator
      , int32_t  const comm_rank // rank of this MPI process
      , uint32_t const nb[3] // number of blocks in X/Y/Z direction
      , int const echo=0 // log level
      , double rank_center[4]=nullptr // export the rank center [0/1/2] and number of items [3]
      , uint16_t *owner_rank=nullptr // export the owner rank of each task, [nb[Z]*nb[Y]*nb[X]]
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace load_balancer
