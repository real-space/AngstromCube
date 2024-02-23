#pragma once
// This file is part of AngstromCube under MIT License

#include "status.hxx" // status_t
#include "simple_stats.hxx" // ::Stats<>

namespace green_parallel {

  int init(int argc=0, char **argv=nullptr); // declaration only

  unsigned size(void); // declaration only

  int rank(void); // declaration only

  int finalize(void); // declaration only

  int max(uint16_t data[], size_t const n); // declaration only

  int allreduce(simple_stats::Stats<double> & stats); // declaration only

  typedef uint16_t rank_int_t;

  status_t dyadic_exchange(                                    // here, nSHO is the global maximum of nsho[ia]
        double       *const mat_out // output effective atom matrices, data layout mat_out[nrows*Noco^2*2*nSHO^2]
      , std::vector<int64_t> const & requests  // indices requested by this MPI process,  [nrows]
      , double const *const mat_inp //  input effective atom matrices, data layout mat_inp[ncols*Noco^2*2*nSHO^2]
      , std::vector<int64_t> const & offerings // indices offered   by this MPI process,  [ncols]
      , rank_int_t const owner_rank[] // where to find it, [nall]
      , uint32_t const nall // number of all atoms
      , int const count
      , bool const debug=true
      , int const echo=0 // log-level
  ); // declaration only

  status_t potential_exchange(
        double    (*const Veff[4])[64]  // output effective potentials,  data layout Veff[Noco^2][nrows][64]
      , std::vector<int64_t> const & requests  // indices requested by this MPI process,         [nrows]
      , double const (*const Vinp)[64]  //  input effective potentials,  data layout Vinp[ncols*Noco^2 ][64]
      , std::vector<int64_t> const & offerings //  indices offered  by this MPI process, [ncols]
      , rank_int_t const owner_rank[] // where to find it, [nb[Z]*nb[Y]*nb[X]]
      , uint32_t const nb[3] // global bounding box, maybe group this with owner_rank into a translator object global_index --> owner_rank
      , int const Noco=1 // 1:no spin, 2: (non-collinear) spin
      , bool const debug=true
      , int const echo=0 // log-level
  ); // declaration only

  class RequestList_t {
  public:

      RequestList_t() : window_size{0} {}
      RequestList_t(
            std::vector<int64_t> const & requests
          , std::vector<int64_t> const & offerings
          , rank_int_t const owner_rank[] // where to find it, [nb[Z]*nb[Y]*nb[X]]
          , uint32_t const nb[3] // global bounding box or {natoms,0,0}
          , int const echo=0 // log-level
      ); // declaration only

  public:
      std::vector<int32_t>  owner; // owner rank of the requested data item
      std::vector<uint32_t> index; // local index in owning process
      std::vector<int64_t>  requested_id; // original identifyer (for debug only)
      std::size_t size() const { return owner.size(); }
      std::size_t window() const { return window_size; }
  private:
      uint32_t window_size;
  }; // RequestList_t

   status_t exchange(
        double       *const data_out // output data, data layout data_out[nrequests*count]
      , double const *const data_inp //  input data, data layout data_inp[nowned   *count]
      , RequestList_t const & requests // contains info which ids are pulled from where
      , int const count=1 // how many doubles per package
      , int const echo=0 // log-level
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

