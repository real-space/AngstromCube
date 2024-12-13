#pragma once
// This file is part of AngstromCube under MIT License

#include <vector> // std::vector<T>

#define HAS_ONESIDED_MPI

#include "status.hxx" // status_t
#include "simple_stats.hxx" // ::Stats<>
#include "mpi_parallel.hxx" // MPI_Comm, MPI_COMM_WORLD
#include "load_balancer.hxx" // rank_int_t



namespace green_parallel {

    typedef load_balancer::rank_int_t rank_int_t;
    auto constexpr no_owner = load_balancer::no_owner;

    class RequestList_t {
    public:

        RequestList_t() : comm_{MPI_COMM_WORLD} {} // default constructor
        RequestList_t( // constructor
              std::vector<int64_t> const & requests
            , std::vector<int64_t> const & offerings
            , rank_int_t const owner_rank[] // where to find it, [nb[Z]*nb[Y]*nb[X]]
            , uint32_t const nb[3] // global bounding box or {natoms,0,0}
            , MPI_Comm const comm=MPI_COMM_WORLD
            , int const echo=0 // log-level
            , char const *const what="?"
        ); // declaration only

    public:
        std::size_t size()   const { return owner.size(); }
        std::size_t window() const { return window_size; }
        MPI_Comm comm() const { return comm_; }

    template <typename real_t=double>
    status_t exchange(
          real_t       *const data_out // output data, data layout data_out[nrequests*count]
        , real_t const *const data_inp //  input data, data layout data_inp[nowned   *count]
        , uint32_t const count=1 // how many real_t per package
        , int const echo=0 // log-level
        , char const *what=nullptr // quantity
    ) const ; // declaration only

#ifdef    HAS_ONESIDED_MPI
    template <typename real_t=double>
    status_t exchange_onesided(
          real_t       *const data_out // output data, data layout data_out[nrequests*count]
        , real_t const *const data_inp //  input data, data layout data_inp[nowned   *count]
        , uint32_t const count=1 // how many real_t per package
        , int const echo=0 // log-level
        , char const *what=nullptr // quantity
    ) const ; // declaration only
#endif // HAS_ONESIDED_MPI

    status_t potential_exchange(
          double    (*const Veff[4])[64]  // output effective potentials,  data layout Veff[Noco^2][nrows][64]
        , double const (*const Vinp)[64]  //  input effective potentials,  data layout Vinp[ncols*Noco^2 ][64]
        , int const Noco=1 // 1:no spin, 2: non-collinear spin
        , int const echo=0 // log-level
    ) const ; // declaration only

    status_t self_test(int const echo=0) const ; // declaration only

    std::vector<int32_t> const & owners() const { return owner; }

#ifdef    HAS_ONESIDED_MPI
    bool get_use1sided() const { return use1sided_; }
#else  // HAS_ONESIDED_MPI
    bool constexpr get_use1sided() const { return false; }
#endif // HAS_ONESIDED_MPI

    private:
        // for 1-sided or 2-sided communication (could be private if we only used 2-sided)
        std::vector<int32_t> owner; // owner rank of the requested data item
        std::vector<int32_t> local_indices; // local index in owning process
        std::vector<int64_t> requested_id; // original identifyer
        std::vector<int64_t> offered_id;   // original identifyer (for self-test)
        uint32_t window_size = 0;
        MPI_Comm comm_;
        // for 2-sided communication only
        std::vector<int32_t> send_packages_to_ranks; // ranks to send data to
        std::vector<std::vector<uint32_t>> send_package_index;
        std::vector<int32_t> recv_packages_from_ranks; // ranks to recveive data from
        std::vector<std::vector<uint32_t>> recv_package_index;
        std::vector<rank_int_t> recv_buffer_index; // if the request is remote, in which recv-buffer is it?
        std::vector<uint32_t> index_in_recv_buffer; // if the request is remote, where in the recv-buffer is it?
#ifdef    HAS_ONESIDED_MPI
        bool use1sided_ = false;
#endif // HAS_ONESIDED_MPI
    }; // class RequestList_t

    status_t all_tests(int const echo=0); // declaration only

} // namespace green_parallel

