// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // uint16_t
#include <vector> // std::vector<T>

#include "green_parallel.hxx" // rank_int_t

#include "status.hxx" // status_t
#include "inline_math.hxx" // set
#include "mpi_parallel.hxx" // ::init, ::size, ::rank, ::finalize, ::min, ::max, ::sum, ::allreduce, ::barrier, MPI_Comm
#include "global_coordinates.hxx" // ::get
#include "print_tools.hxx" // printf_vector
#include "recorded_warnings.hxx" // warn, error
#include "data_view.hxx" // view3D
#include "control.hxx" // ::get

#define DEBUG

#include "debug_output.hxx" // here

namespace green_parallel {

    inline char const * spin_name(int const Noco, int const spin) {
        if (1 == Noco && 0 == spin) return "";
        if (2 == Noco) {
            switch (spin) {
                case 0: return " V_down";
                case 1: return " V_up";
                case 2: return " V_x";
                case 3: return " V_y";
            } // spin
        } // 2 == Noco
        return " ???";
    } // spin_name


    size_t translate(int64_t const global_id, uint32_t const nb[3]) {
        uint32_t xyz[3];
        global_coordinates::get(xyz, global_id);
        for (int d{0}; d < 3; ++d) {
            assert(xyz[d] < nb[d] && "requested coordinates exceed grid box");
        } // d
        return (xyz[2]*size_t(nb[1]) + xyz[1])*size_t(nb[0]) + xyz[0];
    } // translate


    RequestList_t::RequestList_t( // constructor implementation
        std::vector<int64_t> const & requests
      , std::vector<int64_t> const & offerings // do we need this at all? the offerings should match owner_rank[]==me anyway
      , rank_int_t const owner_rank[] // where to find it, [nb[Z]*nb[Y]*nb[X]]
      , uint32_t const nb[3] // global bounding box{nb[X],nb[Y],nb[Z]} or {natoms,0,0}
      , MPI_Comm const comm // =MPI_COMM_WORLD
      , int const echo // =0 // log-level
      , char const *const what // ="?"
    ) {
        comm_ = comm; // copy the MPI communicator used in this context
        auto const nprocs = mpi_parallel::size(comm);
        auto const me     = mpi_parallel::rank(comm, nprocs);

        if (echo > 9) { std::printf("# rank#%i waits in barrier at %s:%d nb=%d %d %d, what=%s\n",
                        me, __FILE__, __LINE__, nb[0], nb[1], nb[2], what); std::fflush(stdout); }
        mpi_parallel::barrier(comm);

        int constexpr X=0, Y=1, Z=2;
        auto const grid = size_t(nb[Z])*size_t(nb[Y])*size_t(nb[X]);
        auto const nall = grid ? grid : nb[X] + nb[Y] + nb[Z];
        auto const nown = offerings.size(); // number of offerings
        auto const nreq = requests.size(); // number of requests
        if (echo > 7) { std::printf("# rank#%i \tRequestList_t [%d %d %d], nall= %ld, offered= %ld, requested= %ld\n",
                                                      me, nb[X],nb[Y],nb[Z], nall, nown, nreq); std::fflush(stdout); }

        int const user_wants_onesided_mpi = control::get("green_parallel.onesided", 0.);
#ifdef    HAS_ONESIDED_MPI
        // use one-sided MPI communication routines or not?
        use1sided_ = (1 == user_wants_onesided_mpi);
#else  // HAS_ONESIDED_MPI
        if (1 == user_wants_onesided_mpi) { warn("found +green_parallel.onesided=%i not compiled with -DHAS_ONESIDED_MPI", user_wants_onesided_mpi); }
#endif // HAS_ONESIDED_MPI
        if (echo > 7) { std::printf("# use %s-sided MPI communication\n", get_use1sided() ? "one" : "two"); std::fflush(stdout); }

        // create a debug aid: nloc_rank
        std::vector<uint32_t> nloc_rank(nprocs, 0); // init all entries as zero
        nloc_rank.at(me) = nown;
        mpi_parallel::sum(nloc_rank.data(), nprocs, comm);
        assert(nown == nloc_rank[me] && "no other rank may add to my contribution!");

#ifndef   HAS_NO_MPI
        bool const debug = 1;

        std::vector<uint32_t> local_index(nall, 0);
        std::vector<rank_int_t> local_check(nall*unsigned(debug), 0);

        for (size_t iloc = 0; iloc < nown; ++iloc) {
            auto const global_id = offerings[iloc];
            assert(global_id > -1);

            // translate global_id into an index iall
            size_t const iall = grid ? translate(global_id, nb) : global_id;
            assert(iall < nall);

            if (me != owner_rank[iall]) { error("rank#%i offers %s id %li but owned by rank#%i", me, what, global_id, owner_rank[iall]); }
            assert(me == owner_rank[iall] && "all offerings must be owned");

            local_index[iall] = iloc;
            assert(iloc == local_index[iall] && "uint32_t insufficient");
            if (debug) { ++local_check[iall]; }
        } // iloc

        if (debug) {
            if (echo > 7) { std::printf("# rank#%i local_check before ", me); printf_vector("%i", local_check); }
            for (size_t iall = 0; iall < nall; ++iall) {
                assert(local_check[iall] <= 1 && "duplicates found");
            } // iall

            auto const stat = mpi_parallel::sum(local_check.data(), nall, comm);
            if (stat) warn("MPI_Allreduce(local_check) failed with status= %i", int(stat));

            if (echo > 7) { std::printf("# rank#%i local_check after  ", me); printf_vector("%i", local_check); }
            for (size_t iall = 0; iall < nall; ++iall) {
                assert(1 == local_check[iall] && "not all covered");
            } // iall
            local_check.resize(0);

            std::vector<rank_int_t> owner_check(nall, 0);
            mpi_parallel::allreduce(owner_check.data(), MPI_MAX, comm, nall, owner_rank);
            if (stat) warn("MPI_Allmax(owner_rank) failed with status= %i", int(stat));
            for (size_t iall = 0; iall < nall; ++iall) {
                if (owner_check[iall] != owner_rank[iall]) {
                    error("rank#%i owner_rank[%li] differs: MPI-maximum is %d but expected %d",
                                me, iall, owner_check[iall], owner_rank[iall]);
                }
                assert(owner_check[iall] == owner_rank[iall] && "owner differs after MPI_MAX");
            } // iall
            mpi_parallel::allreduce(owner_check.data(), MPI_MIN, comm, nall, owner_rank);
            for (size_t iall = 0; iall < nall; ++iall) {
                if (owner_check[iall] != owner_rank[iall]) {
                    error("rank#%i owner_rank[%li] differs: MPI-minimum is %d but expected %d",
                                me, iall, owner_check[iall], owner_rank[iall]);
                }
                assert(owner_check[iall] == owner_rank[iall] && "owner differs after MPI_MIN");
            } // iall
        } // debug

        // get a global list of which local index is where
        {
            // auto const stat = MPI_Allreduce(MPI_IN_PLACE, local_index.data(), nall, MPI_UINT32_T, MPI_MAX, comm);
            auto const stat = mpi_parallel::allreduce(local_index.data(), MPI_MAX, comm, nall);
            if (stat) warn("MPI_Allreduce(local_index) failed with status= %i", int(stat));
        }
        // if this is too expensive see ALTERNATIVE
        //
        // ALTERNATIVE:
        // use MPI_Send and MPI_Recv to distribute the necessary info about local_index 
        // see below
        //

#endif // HAS_NO_MPI

        // initialize member fields
        this->owner = std::vector<int32_t>(nreq, 0); // initialize with master rank for the serial version
        this->local_indices = std::vector<int32_t>(nreq, -1);
        this->requested_id = requests; // deep copy
        this->offered_id = offerings; // deep copy
        this->window_size = nown; // number of owned data items

        size_t not_found{0};
        int64_t id_not_found_1st{-1}, id_not_found_last{-1};

        size_t stats[] = {0, 0, 0}; // {clear, local, remote}

        for (size_t ireq = 0; ireq < nreq; ++ireq) {
            auto const global_id = requests.at(ireq);
            assert(global_id == this->requested_id.at(ireq));
            if (global_id > -1) {
#ifndef   HAS_NO_MPI

                // translate global_id into an index iall
                size_t const iall = grid ? translate(global_id, nb) : global_id;
                assert(iall < nall && "internal index exceeded");

                auto const rank = owner_rank[iall];
                auto const iloc = local_index[iall];

#else  // HAS_NO_MPI

                // without MPI search global_id in offerings
                int64_t iloc_{-1};
                for (size_t iown = 0; iown < nown && iloc_ < 0; ++iown) {
                    if (global_id == offerings.at(iown)) { iloc_ = iown; }
                } // iown
                int32_t const iloc = iloc_;
                if (-1 == iloc) {
                    id_not_found_last = global_id; if (0 == not_found) { id_not_found_1st = global_id; }
                    ++not_found;
                } // not_found
                auto const rank = me;

#endif // HAS_NO_MPI

                if (iloc >= nloc_rank.at(rank)) {
                    error("rank#%i request#%i has owner rank#%i and remote local index %i but maximum is %d",
                                me, ireq, rank, int(iloc), nloc_rank[rank]);
                } // index larger than offered by remote process

                this->owner.at(ireq) = rank;
                this->local_indices.at(ireq) = iloc;
                ++stats[1 + (me != rank)]; // local or remote
            } else { // global_id > -1
                ++stats[0]; // clear
                this->local_indices.at(ireq) = -1;
                this->owner.at(ireq) = no_owner;
            } // global_id > -1
        } // ireq

        if (not_found > 0) {
#ifndef   HAS_NO_MPI
            bool const not_found_is_error = (0 == control::get("mpi.fake.size", 0.));
#else  // HAS_NO_MPI
            auto constexpr not_found_is_error = true;
#endif // HAS_NO_MPI
            if (not_found_is_error) {
                error("rank #%i failed to find %ld global_ids in offerings of \'%s\', 1st id= %li, last id= %li",
                             me, not_found, what, id_not_found_1st, id_not_found_last);
            } else {
                warn( "rank #%i failed to find %ld global_ids in offerings of \'%s\', 1st id= %li, last id= %li",
                             me, not_found, what, id_not_found_1st, id_not_found_last);
            }
        } // not_found

        if (echo > 6) { std::printf("# rank#%i \tRequestList_t expect %.3f k clear, %.3f k copies, %.3f k exchanges\n",
                                          me, stats[0]*1e-3, stats[1]*1e-3, stats[2]*1e-3); std::fflush(stdout); }

        if (echo > 5) { std::printf( "# prepare two-sided communication pattern\n"); std::fflush(stdout); }
        mpi_parallel::barrier(comm);        

        // find out how many packages we need to sendrecv from remote ranks
        std::vector<uint32_t> n_packages_from_rank(nprocs, 0);
        for (size_t ireq{0}; ireq < nreq; ++ireq) {
            if (local_indices.at(ireq) > -1) {
                auto const rank = this->owner.at(ireq);
                if (me != rank) {
                    ++n_packages_from_rank.at(rank);
                }
            }
        } // ireq
        assert(0 == n_packages_from_rank.at(me));

        std::vector<uint32_t> n_packages_to_recv(0); // number of packages to recv
        std::vector<int32_t> rank_index(nprocs, -1); // translation table
        this->recv_packages_from_ranks.resize(0);
        uint32_t ri{0};
        for (int rank{0}; rank < nprocs; ++rank) {
            auto const n_packages = n_packages_from_rank.at(rank);
            if (n_packages > 0) {
                n_packages_to_recv.push_back(n_packages);
                this->recv_packages_from_ranks.push_back(rank);
                assert(-1 == rank_index.at(rank));
                rank_index.at(rank) = ri;
                ++ri;
            }
        } // rank
        auto const n_recv_partners = ri;
        assert(n_recv_partners == this->recv_packages_from_ranks.size());
        assert(n_recv_partners == n_packages_to_recv.size());

        std::vector<std::vector<int64_t>> recv_package_global_id(n_recv_partners);
        this->recv_package_index.resize(n_recv_partners);
        for (uint32_t ri{0}; ri < n_recv_partners; ++ri) {
            this->recv_package_index.at(ri).resize(0);
            recv_package_global_id.at(ri).resize(0);
        } // ri

        for (size_t ireq{0}; ireq < nreq; ++ireq) {
            auto const global_id = requests.at(ireq);
            auto const iloc = local_indices.at(ireq);
            if (iloc > -1) {
                auto const rank = this->owner.at(ireq);
                if (me != rank) {
                    auto const ri = rank_index.at(rank);
                    assert(0 <= ri); assert(ri < n_recv_partners);
                    this->recv_package_index.at(ri).push_back(iloc);
                    recv_package_global_id.at(ri).push_back(global_id);
                }
            }
        } // ireq

        if (echo > 3) { std::printf( "# rank#%i receives from %d other ranks in 2-sided MPI communication\n", me, n_recv_partners); }
        if (echo > 7) { std::printf( "# rank#%i receives from these %d ranks: ", me, n_recv_partners); printf_vector(" %i", this->recv_packages_from_ranks); }
        for (uint32_t ri{0}; ri < n_recv_partners; ++ri) {
            if (echo > 8) { std::printf( "# rank#%i receives these %d local elements from rank#%i : ", me, n_packages_to_recv.at(ri),
                this->recv_packages_from_ranks.at(ri)); printf_vector(" %i", this->recv_package_index.at(ri)); std::fflush(stdout); }
            // consistency check
            assert(this->recv_package_index.at(ri).size() == n_packages_to_recv.at(ri));
        } // ri

        mpi_parallel::barrier(comm);
        // now set up the sender side

        std::vector<uint32_t> n_packages_to_rank(nprocs, 0);
#ifndef   HAS_NO_MPI
        MPI_Alltoall(n_packages_from_rank.data(), 1, MPI_UINT32_T,
                     n_packages_to_rank.data(),   1, MPI_UINT32_T, comm);
#endif // HAS_NO_MPI
        n_packages_from_rank.resize(0);

        std::vector<uint32_t> n_packages_to_send(0);
        this->send_packages_to_ranks.resize(0);

        for (int rank{0}; rank < nprocs; ++rank) {
            auto const n_packages = n_packages_to_rank.at(rank);
            if (n_packages > 0) {
                n_packages_to_send.push_back(n_packages);
                this->send_packages_to_ranks.push_back(rank);
            }
        } // rank
        uint32_t const n_send_partners = this->send_packages_to_ranks.size();
        assert(n_send_partners == n_packages_to_send.size());

        std::vector<std::vector<int64_t>> send_package_global_id(n_send_partners);
        this->send_package_index.resize(n_send_partners);
        for (uint32_t rj{0}; rj < n_send_partners; ++rj) {
            auto const n_packages = n_packages_to_send.at(rj);
            this->send_package_index.at(rj).resize(n_packages);
            send_package_global_id.at(rj).resize(n_packages);
        } // rj

        // now send the request indices to the sender ranks
#ifndef   HAS_NO_MPI
        { // scope: exchange indices, slightly confusing in terms of naming ....
            //  ... but yes, we send the list of indices that we want to receive ...
            //  ... and we receive the list of indices we need to send.
            int const tag = __LINE__;
            auto const nr = n_recv_partners + n_send_partners;
            std::vector<MPI_Request> mpi_req(nr);

            for (uint32_t ri{0}; ri < n_recv_partners; ++ri) {
                auto const rank = recv_packages_from_ranks.at(ri);
                MPI_Isend(this->recv_package_index.at(ri).data(), this->recv_package_index.at(ri).size(), 
                            MPI_UINT32_T, rank, tag, comm, &mpi_req.at(ri));
            } // ri

            for (uint32_t rj{0}; rj < n_send_partners; ++rj) {
                auto const rank = send_packages_to_ranks.at(rj);
                MPI_Irecv(this->send_package_index.at(rj).data(), this->send_package_index.at(rj).size(),
                            MPI_UINT32_T, rank, tag, comm, &mpi_req.at(n_recv_partners + rj));
            } // rj

            MPI_Waitall(nr, mpi_req.data(), MPI_STATUSES_IGNORE);


            // now also communicate the global indices

            for (uint32_t ri{0}; ri < n_recv_partners; ++ri) {
                auto const rank = recv_packages_from_ranks.at(ri);
                MPI_Isend(recv_package_global_id.at(ri).data(), recv_package_global_id.at(ri).size(), 
                            MPI_INT64_T, rank, tag, comm, &mpi_req.at(ri));
            } // ri

            for (uint32_t rj{0}; rj < n_send_partners; ++rj) {
                auto const rank = send_packages_to_ranks.at(rj);
                MPI_Irecv(send_package_global_id.at(rj).data(), send_package_global_id.at(rj).size(),
                            MPI_INT64_T, rank, tag, comm, &mpi_req.at(n_recv_partners + rj));
            } // rj

            MPI_Waitall(nr, mpi_req.data(), MPI_STATUSES_IGNORE);

            // check that all global_ids requested from this rank are owned locally
            for (uint32_t rj{0}; rj < n_send_partners; ++rj) {
                auto const n_packages = send_package_global_id.at(rj).size();
                for (uint32_t ip{0}; ip < n_packages; ++ip) {
                    auto const global_id = send_package_global_id.at(rj).at(ip);
                    assert(global_id > -1); // vacuum cells do not to be communicated
                    size_t const iall = grid ? translate(global_id, nb) : global_id;
                    assert(iall < nall && "internal index exceeded");
                    assert(me == owner_rank[iall] && "all packages marked for sending must be owned locally!");
                } // ip
            } // rj

        } // scope
#endif // HAS_NO_MPI

        assert(me == rank_int_t(me));
        this->recv_buffer_index = std::vector<rank_int_t>(nreq, rank_int_t(me)); // if the request is remote, in which recv-buffer is it?
        this->index_in_recv_buffer.resize(nreq, 0) ; // if the request is remote, where in the recv-buffer is it?

        size_t new_stats[] = {0, 0, 0}; // get element from {0:clear, 1:local 2:remote, 2:clear}
        for (size_t ireq = 0; ireq < nreq; ++ireq) {
            auto const global_id = this->requested_id.at(ireq);
            auto const rank      = this->owner.at(ireq);
            auto const iloc      = this->local_indices.at(ireq);

            if (no_owner == rank) {
                ++new_stats[0]; // clear
                assert(-1 == iloc);
                assert(-1 == global_id);
            } else if (me == rank) {
                ++new_stats[1]; // local copy
                if (echo > 18) std::printf("# exchange: rank#%i get data of item#%lli  copy local element %i\n", me, global_id, iloc);
                assert(iloc < nown);
            } else { // me == rank
                ++new_stats[2]; // remote access or sendrecv
                if (echo > 17) std::printf("# exchange: rank#%i get data of item#%lli from rank#%i element %i\n", me, global_id, rank, iloc);
#ifndef   HAS_NO_MPI
                if (rank >= nprocs) { error("rank#%i tries to MPI_Get from rank#%i but only %d processes running, global_id=%li",
                                            me, rank, nprocs, global_id); }
                assert(0 <= rank); assert(rank < nprocs);
                assert(iloc < nloc_rank.at(rank));
                auto const ri = rank_index.at(rank);
                assert(ri >= 0 && "did not expect elements from this rank, error in RequestList_t constructor");
                // now which package in the buffer belongs to iloc?
                int ibuf{-1};
                auto const n_packages = this->recv_package_index.at(ri).size();
                for (uint32_t ip{0}; ip < n_packages; ++ip) {
                    if (iloc == this->recv_package_index.at(ri).at(ip)) { ibuf = ip; }
                } // ip
                if (-1 == ibuf) {
                    error("rank#%i failed to find local index %i in buffer received from rank#%i of %d, global_id=%li",
                                                me, iloc, rank, nprocs, global_id);
                } // not found
                if (echo > 27) { std::printf("# rank#%i found item#%lli in buffer[%i] from rank#%i\n",
                                                me, global_id, ibuf, rank); std::fflush(stdout); }
                index_in_recv_buffer.at(ireq) = ibuf;
                recv_buffer_index.at(ireq) = ri;
#else  // HAS_NO_MPI
                error("Without MPI all entries must reside in the same process, me=%i, owner=%i", me, rank);
#endif // HAS_NO_MPI
            } // me == rank
        } // ireq

        for (int i3{0}; i3 < 3; ++i3) { assert(stats[i3] == new_stats[i3]); }

        mpi_parallel::sum(stats, 3, comm);
        if (echo > 5) { std::printf( "# total  \tRequestList_t expect %.3f k clear, %.3f k copies, %.3f k exchanges\n",
                                              stats[0]*1e-3, stats[1]*1e-3, stats[2]*1e-3); std::fflush(stdout); }

        auto const stat = this->self_test(echo);
        if (0 != stat) { warn("RequestList_t constructor failed in self_check with status= %i", int(stat)); }

    } // constructor implementation


#ifdef    HAS_ONESIDED_MPI

    template <typename real_t>
    status_t RequestList_t::exchange_onesided(
          real_t       *const data_out // output data, data layout data_out[nrequests*count]
        , real_t const *const data_inp //  input data, data layout data_inp[nowned   *count]
        , uint32_t const count // number of real_t per package
        , int const echo // =0, log-level
        , char const *what // =nullptr // quantity
    ) const {
        status_t status(0);
        what = what ? what : "?";
        auto const comm = this->comm();
        auto const nprocs = mpi_parallel::size(comm); // number of processes
        auto const me = mpi_parallel::rank(comm, nprocs);

        // The number of local atoms is limited to 2^16 == 65536
        if (echo > 5) std::printf("# exchange using MPI one-sided communication, packages of %d numbers, %.3f kByte %s\n",
                                                                                  count, count*sizeof(real_t)*.001, what);
        auto const nreq = this->size(); // number of requests
        auto const nwin = this->window(); // number of offerings
        if (nullptr == data_out) assert(0 == nreq && "may not be called with a nullptr for output");
        if (nullptr == data_inp) assert(0 == nwin && "may not be called with a nullptr for input");

#ifndef   HAS_NO_MPI
        // set up a memory window to read from
        MPI_Win window;
        uint const disp_unit = count*sizeof(real_t); // in Bytes
        size_t const win_size = nwin*disp_unit; // in Bytes
        int const assertions = MPI_MODE_NOPUT; // use bitwise or, e.g. MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE | MPI_MODE_NOSUCCEED;
        auto const data_type = mpi_parallel::get(real_t(0));
        // since the window descriptor depends on disp_unit, we cannot create the window in the constructor
        status += MPI_Win_create((void*)data_inp, win_size, disp_unit, MPI_INFO_NULL, comm, &window);
        // synchronize processes
        status += MPI_Win_fence(assertions, window);
#endif // HAS_NO_MPI

        for (size_t ireq = 0; ireq < nreq; ++ireq) {
            auto const global_id = this->requested_id.at(ireq);
            auto const rank      = this->owner.at(ireq);
            auto const iloc      = this->local_indices.at(ireq);
            if (no_owner == rank) {
                assert(-1 == global_id);
                set(&data_out[ireq*count], count, real_t(0)); // clear package
            } else if (me == rank) {
                if (echo > 18) std::printf("# exchange: rank#%i get data of item#%lli  copy local element %i\n", me, global_id, iloc);
                assert(iloc < nwin);
                set(&data_out[ireq*count], count, &data_inp[iloc*count]); // copy package
            } else { // me == rank
                if (echo > 17) std::printf("# exchange: rank#%i get data of item#%lli from rank#%i element %i\n", me, global_id, rank, iloc);
#ifndef   HAS_NO_MPI
                assert(rank >= 0); assert(rank < nprocs);
                // get package from remote process via RDMA
                status += MPI_Get(&data_out[ireq*count], count, data_type, rank, iloc, count, data_type, window);
#else  // HAS_NO_MPI
                ++status; // cannot do that without MPI
#endif // HAS_NO_MPI
            } // me == rank
        } // ireq

#ifndef   HAS_NO_MPI
        // synchronize processes
        status += MPI_Barrier(comm); // { std::printf("# rank#%d hit MPI_Barrier at line %d\n", me, __LINE__); std::fflush(stdout); }
        status += MPI_Win_fence(assertions, window);
        status += MPI_Win_free(&window);
#endif // HAS_NO_MPI
        return status;
    } // RequestList_t::exchange_onesided

#endif // HAS_ONESIDED_MPI


    template <typename real_t>
    status_t RequestList_t::exchange(
          real_t       *const data_out // output data, data layout data_out[nrequests*count]
        , real_t const *const data_inp //  input data, data layout data_inp[nowned   *count]
        , uint32_t const count // number of real_t per package
        , int const echo // =0, log-level
        , char const *what // =nullptr // quantity
    ) const {

        what = what ? what : "?";
        auto const comm = this->comm();
        auto const nprocs = mpi_parallel::size(comm); // number of processes
        auto const me = mpi_parallel::rank(comm, nprocs);

        auto const nreq = this->size(); // number of requests
        auto const nwin = this->window(); // number of offerings
        if (nullptr == data_out) assert(0 == nreq && "may not be called with a nullptr for output");
        if (nullptr == data_inp) assert(0 == nwin && "may not be called with a nullptr for input");

#ifdef    HAS_ONESIDED_MPI
        if (get_use1sided()) {
            return this->exchange_onesided(data_out, data_inp, count, echo, what);
        } // use one-sided MPI communication routines
#endif // HAS_ONESIDED_MPI
        if (echo > 5) std::printf("# exchange using MPI two-sided communication, packages of %d numbers, %.3f kByte %s\n",
                                                                                  count, count*sizeof(real_t)*.001, what);
        status_t status(0);

#ifndef   HAS_NO_MPI
        auto const data_type = mpi_parallel::get(real_t(0));

        auto const ns = this->send_packages_to_ranks.size();
        auto const nr = this->recv_packages_from_ranks.size();
        std::vector<MPI_Request> mpi_req(ns + nr);

        int const tag = sizeof(real_t);

        std::vector<std::vector<real_t>> send_buff(ns);
        for (uint32_t rj{0}; rj < ns; ++rj) {
            auto & buffer = send_buff.at(rj);
            auto const n_packages = this->send_package_index.at(rj).size();
            buffer.resize(n_packages*count, real_t(0));
            for (uint32_t ip{0}; ip < n_packages; ++ip) {
                auto const iloc = this->send_package_index.at(rj).at(ip);
                set(&buffer[ip*count], count, &data_inp[iloc*count]);
            } // ip
            auto const rank = this->send_packages_to_ranks.at(rj);
            // int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
            MPI_Isend(buffer.data(), buffer.size(), data_type, rank, tag, comm, &mpi_req.at(rj));
        } // rj

        std::vector<std::vector<real_t>> recv_buff(nr);
        for (uint32_t ri{0}; ri < nr; ++ri) {
            auto & buffer = recv_buff.at(ri);
            auto const n_packages = this->recv_package_index.at(ri).size();
            buffer.resize(n_packages*count, real_t(0));
            auto const rank = this->recv_packages_from_ranks.at(ri);
            // int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
            MPI_Irecv(buffer.data(), buffer.size(), data_type, rank, tag, comm, &mpi_req.at(ns + ri));
        } // ri

        // wait for all messages to be sent and to have arrived
        MPI_Waitall(mpi_req.size(), mpi_req.data(), MPI_STATUSES_IGNORE);

        send_buff.resize(0); // release memory of senders

#endif // HAS_NO_MPI

        for (size_t ireq = 0; ireq < nreq; ++ireq) {
            auto const rank = this->owner.at(ireq);
            if (no_owner == rank) {
                assert(-1 == this->requested_id.at(ireq));
                set(&data_out[ireq*count], count, real_t(0)); // clear package
            } else if (me == rank) {
                auto const iloc = this->local_indices[ireq];
                if (echo > 18) std::printf("# exchange: rank#%i get data of item#%lli  copy local element %i\n", me, this->requested_id.at(ireq), iloc);
                assert(iloc < nwin);
                set(&data_out[ireq*count], count, &data_inp[iloc*count]); // copy package
            } else { // me == rank
                auto const ibuf = this->index_in_recv_buffer.at(ireq);
                if (echo > 17) std::printf("# exchange: rank#%i get data of item#%lli from rank#%i buffer[%i]\n", me, this->requested_id.at(ireq), rank, ibuf);
#ifndef   HAS_NO_MPI
                assert(0 <= rank); assert(rank < nprocs);
                auto const ri = this->recv_buffer_index.at(ireq);
                assert(ri >= 0 && "did not expect elements from this rank, error in RequestList_t constructor");
                auto const & buffer = recv_buff.at(ri);
                set(&data_out[ireq*count], count, &buffer[ibuf*count]); // copy package from receive buffer
#else  // HAS_NO_MPI
                ++status; // cannot do that without MPI
#endif // HAS_NO_MPI
            } // me == rank
        } // ireq

        mpi_parallel::barrier(comm); // synchronize processes
        return status;
    } // RequestList_t::exchange


    status_t RequestList_t::potential_exchange(
          double    (*const Veff[4])[64]  // output effective potentials,  data layout Veff[Noco*Noco][nreq][64]
        , double const (*const Vinp)[64]  //  input effective potentials,  data layout Vinp[ncols*Noco*Noco][64]
        , int const Noco // =1, 1:no spin, 2: (non-collinear) spin
        , int const echo // =0, log-level
    ) const {
        if (echo > 0) std::printf("# MPI data exchange of potential, Noco=%d\n", Noco);
        assert(1 == Noco || 2 == Noco);

        assert(Veff && "may not be called with a nullptr for output");
        for (int spin = 0; spin < Noco*Noco; ++spin) {
            assert(Veff[spin] && "may not be called with a nullptr for spin output");
        } // spin
        assert(Vinp && "may not be called with a nullptr for input");

        auto const nreq = this->size();
        view3D<double> Vout(nreq,Noco*Noco,64, 0.0); // get temporary CPU memory in [nreq][Noco*Noco][64] layout

        auto const status = this->exchange(Vout.data(), Vinp[0], Noco*Noco*64, echo, "potential");

        // convert Vout[nreq][Noco*Noco][64] into special data layout of Veff[Noco*Noco][nreq][64] (in GPU memory) 
        for (size_t ireq = 0; ireq < nreq; ++ireq) {
            for (int spin = 0; spin < Noco*Noco; ++spin) {
                set(Veff[spin][ireq], 64, Vout(ireq,spin)); // copy blocks of 4*4*4 grid points
            } // spin
        } // ireq

        return status;
    } // RequestList_t::potential_exchange 


    status_t RequestList_t::self_test(int const echo) const {
        // sanity check routine testing exchange with 1 global_id per package
        uint32_t const nr = this->size();   // number of requested elements
        uint32_t const ns = this->window(); // number of offered elements
        typedef float real_t; // if real_t == float, this makes the explicit template instantiation below redundant
        std::vector<real_t> inp(ns, real_t(0));
        for (uint32_t is{0}; is < ns; ++is) {
            inp.at(is) = real_t(this->offered_id.at(is)); // input are the global ids offered, converted to real_t
        } // is
        std::vector<real_t> out(nr, real_t(0));
        auto stat = this->exchange(out.data(), inp.data(), 1, echo, "global_ids in self_test");
        if (0 != stat) {
            warn("exchange of global_ids as self_test failed with status= %i", int(stat));
            stat = 0;
        }
        // now check if the requested ids have been transmitted
        for (uint32_t ir{0}; ir < nr; ++ir) {
            auto const reference_id = this->requested_id.at(ir);
            stat += (out.at(ir) != real_t(reference_id)) * (reference_id >= 0);
            // reference_id == -1 will be mapped to 0.0f;
        } // ir
        return stat;
    } // RequestList_t::self_test

    // template // explicit template instantiation for real_t=double
    // status_t RequestList_t::exchange(double*, double const*, uint32_t, int, char const*) const;

    // template // explicit template instantiation for real_t=float
    // status_t RequestList_t::exchange(float* , float  const*, uint32_t, int, char const*) const;





















#ifdef    NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    status_t test_exchanges(int echo=0, int const nSHO=20) {
        status_t stat(0);
        uint32_t const nb[] = {2, 2, 2}; // grid box is 2x2x2 or 8 atoms
        std::vector<int64_t> requests = {0,7,6,1,5,2,4,3}; // every process requests all 8 of these ids

        auto const nrows = requests.size();
        auto const comm = mpi_parallel::comm();
        auto const nprocs = mpi_parallel::size(comm);       assert(nprocs > 0);
        auto const me = mpi_parallel::rank(comm, nprocs);

        std::vector<int64_t> offerings(0);
        uint32_t const nall = nb[2]*nb[1]*nb[0];
        std::vector<uint16_t> owner_rank(nall, 0);
        for(int id{0}; id < nall; ++id) {
            int const rank = id % nprocs; // block-cyclic distribution
            owner_rank[id] = rank;
            if (me == rank) offerings.push_back(id);
        } // id
        uint32_t const na[] = {nall, 0, 0};
        RequestList_t rlV(requests, offerings, owner_rank.data(), nb, comm, echo, "test_V");
        RequestList_t rlD(requests, offerings, owner_rank.data(), na, comm, echo, "test_D");

        view3D<double> pot_out_memory(2*2,nrows,64, 0.0);
        double (*pot_out[2*2])[64];
        for (int spin{0}; spin < 2*2; ++spin) {
            pot_out[spin] = (double(*)[64]) pot_out_memory(spin,0); // prepare for the Noco=2 test case
        } // spin

        int const ncols = offerings.size();
        int nerr{0};
        for (int Noco = 1; Noco <= 2; ++Noco) {
            mpi_parallel::barrier(comm);
            {
                auto const pot_inp = new double[ncols*Noco*Noco][64];
                for (int col{0}; col < ncols; ++col) { pot_inp[col*Noco*Noco][0] = 0.5 + me; } // ear-mark with owner rank
                stat += rlV.potential_exchange(pot_out, pot_inp, Noco, echo);
                for (int row{0}; row < nrows; ++row) { nerr += (pot_out[0][row][0] != (0.5 + owner_rank[requests[row]])); }
                if (echo > 0) { std::printf("# potential_exchange Noco= %d status= %d errors= %d\n\n", Noco, int(stat), nerr); std::fflush(stdout); }
                delete[] pot_inp;
            }
            mpi_parallel::barrier(comm);
            {
                char const *const what = (Noco > 1) ? "testMatNoco=2" : "testMatNoco=1";
                int const count = Noco*Noco*2*nSHO*nSHO; // number of doubles per package
                std::vector<double> mat_out(nrows*count), mat_inp(ncols*count);
                for (int col{0}; col < ncols; ++col) { mat_inp[col*count] = 0.5 + me; } // ear-mark with owner rank
                stat += rlD.exchange(mat_out.data(), mat_inp.data(), count, echo, what);
                for (int row{0}; row < nrows; ++row) { nerr += (mat_out[row*count] != (0.5 + owner_rank[requests[row]])); }
                if (echo > 0) { std::printf("# rank#%i exchange Noco= %d status= %d errors= %d\n\n", me, Noco, int(stat), nerr); std::fflush(stdout); }
            }
            mpi_parallel::barrier(comm);
        } // Noco

        return stat;
    } // test_exchanges

    status_t all_tests(int const echo) {
        status_t stat(0);
        bool const already_initialized = mpi_parallel::init();
        stat += test_exchanges(echo);
        if (!already_initialized) mpi_parallel::finalize();
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_parallel
