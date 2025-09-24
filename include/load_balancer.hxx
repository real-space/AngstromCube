#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdint> // uint32_t, uint16_t
#include <vector> // std::vector


#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

namespace load_balancer {

  typedef uint16_t rank_int_t;
  rank_int_t constexpr no_owner = (1ull << 16) - 1; // 65535

  struct WeightInfo
  {
    uint32_t weightContributionForKinetic = 0; // Homogeneous
    uint32_t weightContributionForSHOadd  = 0; // Basisfunctions
    uint32_t weightContributionForSHOprj  = 0; // AtomImages * number of corners * 8, Note map 9 to 8
  };

  std::vector<float> calculate_weights(std::vector<int64_t>& global_source_indices, std::vector<load_balancer::WeightInfo>& weight_infos, size_t grid_size);

  double get(
        uint32_t const comm_size // number of MPI processes in this communicator
      , int32_t  const comm_rank // rank of this MPI process
      , uint32_t const nb[3] // number of cubes in X/Y/Z direction
      , float const *const block_weights=nullptr // stores the weight of each block, [nb[Z]*nb[Y]*nb[X]] 
      , int const echo=0 // log level
      , double rank_center[4]=nullptr // export the rank center [0/1/2] and number of items [3]
      , rank_int_t *owner_rank=nullptr // export the owner rank of each task, [nb[Z]*nb[Y]*nb[X]], needs an MPI_MAX
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace load_balancer
