#pragma once
// This file is part of AngstromCube under MIT License

#include "status.hxx" // status_t
#include "simple_stats.hxx" // ::Stats<>

namespace green_parallel {

  typedef uint16_t rank_int_t;

  class RequestList_t {
  public:

      RequestList_t() : window_size{0} {}
      RequestList_t( // constructor
            std::vector<int64_t> const & requests
          , std::vector<int64_t> const & offerings
          , rank_int_t const owner_rank[] // where to find it, [nb[Z]*nb[Y]*nb[X]]
          , uint32_t const nb[3] // global bounding box or {natoms,0,0}
          , int const echo=0 // log-level
      ); // declaration only

  public:
      std::vector<int32_t> owner; // owner rank of the requested data item
      std::vector<int32_t> index; // local index in owning process
      std::vector<int64_t> requested_id; // original identifyer (for debug only)
      std::size_t size() const { return owner.size(); }
      std::size_t window() const { return window_size; }
  private:
      uint32_t window_size;
  }; // RequestList_t

  template <typename real_t=double>
  status_t exchange(
        real_t       *const data_out // output data, data layout data_out[nrequests*count]
      , real_t const *const data_inp //  input data, data layout data_inp[nowned   *count]
      , RequestList_t const & requests // contains info which ids are pulled from where
      , uint32_t const count=1 // how many real_t per package
      , int const echo=0 // log-level
      , char const *what=nullptr // quantity
  ); // declaration only

  status_t potential_exchange(
        double    (*const Veff[4])[64]  // output effective potentials,  data layout Veff[Noco^2][nrows][64]
      , double const (*const Vinp)[64]  //  input effective potentials,  data layout Vinp[ncols*Noco^2 ][64]
      , RequestList_t const & requests // contains info which ids are pulled from where
      , int const Noco=1 // 1:no spin, 2: (non-collinear) spin
      , int const echo=0 // log-level
  ); // declaration only

  status_t all_tests(int echo=0); // declaration only

} // namespace green_parallel

