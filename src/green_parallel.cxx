// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // uint16_t
#include <vector> // std::vector<T>

#include "green_parallel.hxx" // rank_int_t

#include "status.hxx" // status_t
#include "inline_math.hxx" // set
#include "mpi_parallel.hxx" // ::init, ::size, ::rank, ::finalize, ::min, ::max, ::sum, ::allreduce, ::comm, ::barrier
#include "global_coordinates.hxx" // ::get
#include "print_tools.hxx" // printf_vector
#include "recorded_warnings.hxx" // warn
#include "data_view.hxx" // view3D

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


    RequestList_t::RequestList_t( // constructor implementation
        std::vector<int64_t> const & requests
      , std::vector<int64_t> const & offerings // do we need this at all? the offerings should match owner_rank[]==me anyway
      , rank_int_t const owner_rank[] // where to find it, [nb[Z]*nb[Y]*nb[X]]
      , uint32_t const nb[3] // global bounding box or {natoms,0,0}
      , int const echo // =0 // log-level
      , char const *const what // ="?"
    ) {
        auto const comm = mpi_parallel::comm(); // MPI_COMM_WORLD
        auto const me   = mpi_parallel::rank(comm);

        if (echo > 9) { std::printf("# rank#%i waits in barrier at %s:%d nb=%d %d %d\n", me, __FILE__, __LINE__, nb[0], nb[1], nb[2]); std::fflush(stdout); }
        mpi_parallel::barrier(comm);


        int constexpr X=0, Y=1, Z=2;
        auto const grid = size_t(nb[Z])*size_t(nb[Y])*size_t(nb[X]);
        auto const nall = grid ? grid : nb[X] + nb[Y] + nb[Z];
        auto const nown = offerings.size(); // number of offerings
        auto const nreq = requests.size(); // number of requests
        if (echo > 7) std::printf("# rank#%i \tRequestList_t [%d %d %d], nall= %ld, offered= %ld, requested= %ld\n",
                                          me,         nb[X],nb[Y],nb[Z], nall,              nown,           nreq);

#ifndef   HAS_NO_MPI
        bool const debug = 1;

        std::vector<uint16_t> local_index(nall, 0);
        std::vector<uint16_t> local_check(nall*unsigned(debug), 0);

        assert(nown <= (1ul << 16) && "each process can hold max 2^16 locally owned items");
        for (size_t iown = 0; iown < nown; ++iown) {
            auto const global_id = offerings[iown];

            // translate global_id into an index iall
            assert(global_id > -1);
            size_t iall = global_id;
            if (grid) {
                uint32_t xyz[3]; global_coordinates::get(xyz, global_id);
                for (int d = 0; d < 3; ++d) assert(xyz[d] < nb[d] && "requested coordinates exceed box");
                iall = (xyz[Z]*size_t(nb[Y]) + xyz[Y])*nb[X] + xyz[X];
            } // grid
            assert(iall < nall);

            assert(me == owner_rank[iall] && "all offerings must be owned");

            local_index[iall] = iown;
            if (debug) ++local_check[iall];
        } // iown

        if (debug) {
            if (echo > 7) { std::printf("# rank#%i local_check before ", me); printf_vector(" %i", local_check); }
            for (size_t iall = 0; iall < nall; ++iall) {
                assert(local_check[iall] <= 1 && "duplicates found");
            } // iall

            auto const stat = mpi_parallel::sum(local_check.data(), nall, comm);
            if (stat) warn("MPI_Allreduce(local_check) failed with status= %d", int(stat));

            if (echo > 7) { std::printf("# rank#%i local_check after  ", me); printf_vector(" %i", local_check); }
            for (size_t iall = 0; iall < nall; ++iall) {
                assert(1 == local_check[iall] && "not all covered");
            } // iall
            local_check.resize(0);

            std::vector<rank_int_t> owner_check(nall, 0);
            mpi_parallel::allreduce(owner_check.data(), MPI_MAX, comm, nall, owner_rank);
            if (stat) warn("MPI_Allmax(owner_rank) failed with status= %d", int(stat));
            for (size_t iall = 0; iall < nall; ++iall) {
                assert(owner_check[iall] == owner_rank[iall] && "owner differs after MPI_MAX");
            } // iall
            mpi_parallel::allreduce(owner_check.data(), MPI_MIN, comm, nall, owner_rank);
            for (size_t iall = 0; iall < nall; ++iall) {
                assert(owner_check[iall] == owner_rank[iall] && "owner differs after MPI_MIN");
            } // iall
        } // debug

        // get a global list of which local index is where
        auto const stat = MPI_Allreduce(MPI_IN_PLACE, local_index.data(), nall, MPI_UINT16_T, MPI_MAX, comm);
        if (stat) warn("MPI_Allreduce(local_index) failed with status= %d", int(stat));
        // if this is too expensive see ALTERNATIVE
        //
        // ALTERNATIVE:
        // use MPI_Send and MPI_Recv to distribute the necessary info about local_index 
        //

#endif // HAS_NO_MPI

        // initialize member fields
        owner = std::vector<int32_t>(nreq, 0); // initialize with master rank for the serial version
        index = std::vector<int32_t>(nreq, -1);
        requested_id = std::vector<int64_t>(nreq, 0);
        window_size = nown; // number of owned data items

        size_t stats[] = {0, 0, 0}; // {local, remote, clear}

        for (size_t ireq = 0; ireq < nreq; ++ireq) {
            auto const global_id = requests[ireq];

            if (global_id > -1) {
#ifndef   HAS_NO_MPI

                // translate global_id into an index iall
                size_t iall = global_id;
                if (grid) {
                    uint32_t xyz[3]; global_coordinates::get(xyz, global_id);
                    for (int d = 0; d < 3; ++d) assert(xyz[d] < nb[d] && "requested coordinates exceed box");
                    iall = (xyz[Z]*size_t(nb[Y]) + xyz[Y])*nb[X] + xyz[X];
                } // grid
                assert(iall < nall && "internal index exceeded");

                owner[ireq]     =  owner_rank[iall];
             // std::printf("# rank#%i request#%i owner=rank#%i\n", me, ireq, owner[ireq]);
                auto const iloc = local_index[iall];

#else  // HAS_NO_MPI

                // without MPI search global_id in offerings
                int64_t iloc{-1};
                for (size_t iown = 0; iown < nown && iloc < 0; ++iown) {
                    if (global_id == offerings[iown]) iloc = iown;
                } // iown
                if (-1 == iloc) error("failed to find global_id= %li in offerings of %s", global_id, what);
                assert(iloc > -1 && "index not found in local offerings");
                owner[ireq] = me;

#endif // HAS_NO_MPI
                ++stats[me != owner[ireq]];

                index[ireq] = iloc;
            } else { // global_id > -1
                index[ireq] = -1;
                owner[ireq] = me;
                ++stats[2]; // clear
            } // global_id > -1

            requested_id[ireq] = global_id; // copy the array of requests
        } // ireq

        if (echo > 6) std::printf("# rank#%i \tRequestList_t expect %.3f k copies, %.3f k exchanges, %.3f k initializations\n",
                                          me, stats[0]*1e-3, stats[1]*1e-3, stats[2]*1e-3);
        mpi_parallel::sum(stats, 3, comm);
        if (echo > 5) std::printf( "# total  \tRequestList_t expect %.3f k copies, %.3f k exchanges, %.3f k initializations\n",
                                              stats[0]*1e-3, stats[1]*1e-3, stats[2]*1e-3);
    } // constructor implementation


    template <typename real_t>
    status_t exchange(
            real_t       *const data_out // output data, data layout data_out[nrequests*count]
        , real_t const *const data_inp //  input data, data layout data_inp[nowned   *count]
        , RequestList_t const & requests
        , uint32_t const count // number of real_t per package
        , int const echo // =0, log-level
        , char const *what // =nullptr // quantity
    ) {
        what = what ? what : "?";
        auto const comm = mpi_parallel::comm(); // MPI_COMM_WORLD
        auto const me = mpi_parallel::rank(comm);

        // The number of local atoms is limited to 2^16 == 65536
        if (echo > 0) std::printf("# exchange using MPI one-sided communication, packages of %.3f k numbers, %.3f kByte %s\n",
                                                                                count*.001, count*sizeof(real_t)*.001, what);
        auto const nreq = requests.size(); // number of requests
        auto const nwin = requests.window(); // number of offerings
        if (nullptr == data_out) assert(0 == nreq && "may not be called with a nullptr for output");
        if (nullptr == data_inp) assert(0 == nwin && "may not be called with a nullptr for input");

        status_t status(0);

    #ifndef   HAS_NO_MPI

        auto const np = mpi_parallel::size(comm); // number of processes
        // set up a memory window to read from
        MPI_Win window;
        uint const disp_unit = count*sizeof(real_t); // in Bytes
        size_t const win_size = nwin*disp_unit; // in Bytes
        status += MPI_Win_create((void*)data_inp, win_size, disp_unit, MPI_INFO_NULL, comm, &window);
        auto const data_type = (8 == sizeof(real_t)) ? MPI_DOUBLE : MPI_FLOAT;

        // synchronize processes
        int const assertions = MPI_MODE_NOPUT; // use bitwise or, e.g. MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE | MPI_MODE_NOSUCCEED;
        status += MPI_Win_fence(assertions, window);

    #endif // HAS_NO_MPI

        size_t stats[] = {0, 0, 0}; // get element from {0:local, 1:remote, 2:clear}
        for (size_t ireq = 0; ireq < nreq; ++ireq) {
            auto const global_id = requests.requested_id[ireq];
            auto const owner     = requests.owner[ireq];
            auto const iloc      = requests.index[ireq];

            if (iloc >= 0) {
                if (me == owner) {
                    if (echo > 7) std::printf("# exchange: rank#%i get data of item#%lli  copy local element %i\n", me, global_id, iloc);
                    set(&data_out[ireq*count], count, &data_inp[iloc*count]); // copy
                    ++stats[0]; // local
                } else { // me == owner
                    ++stats[1]; // remote
#ifndef   HAS_NO_MPI
                    if (echo > 7) std::printf("# exchange: rank#%i get data of item#%lli from rank#%i element %i\n", me, global_id, owner, iloc);
                    assert(owner >= 0);
                    if (owner >= np) error("rank#%i tries to MPI_Get %.3f kByte from rank#%i but only %d processes running, global_id=%li", me, count*sizeof(real_t)*.001, owner, np, global_id);
                    status += MPI_Get(&data_out[ireq*count], count, data_type, owner, iloc, count, data_type, window);
#else  // HAS_NO_MPI
                    error("Without MPI all atom matrices must reside in the same process, me=%i, owner=%i", me, owner);
#endif // HAS_NO_MPI
                } // me == owner
            } else {
                ++stats[2]; // clear
                assert(-1 == iloc);
                assert(-1 == global_id);
                set(&data_out[ireq*count], count, real_t(0)); // clear
            }
        } // ireq

#ifndef   HAS_NO_MPI
        // synchronize processes
        status += MPI_Barrier(comm); // { std::printf("# rank#%d hit MPI_Barrier at line %d\n", me, __LINE__); std::fflush(stdout); }
        status += MPI_Win_fence(assertions, window);
        status += MPI_Win_free(&window);
#endif // HAS_NO_MPI

        if (echo > 6) std::printf("# rank#%i \tcopied %.3f k, pulled %.3f k and cleared %.3f k elements\n",
                                            me, stats[0]*.001, stats[1]*.001, stats[2]*.001);
        assert(stats[0] + stats[1] + stats[2] == nreq);
        mpi_parallel::sum(stats, 3, comm);
        if (echo > 5) std::printf("# total  \tcopied %.3f k, pulled %.3f k and cleared %.3f k elements\n",
                                                stats[0]*.001, stats[1]*.001, stats[2]*.001);

        return status;
    } // exchange

    template // explicit template instantiation for real_t=double
    status_t exchange(double*, double const*, RequestList_t const &, uint32_t, int, char const*);

    template // explicit template instantiation for real_t=float
    status_t exchange(float* , float  const*, RequestList_t const &, uint32_t, int, char const*);

    status_t potential_exchange(
            double    (*const Veff[4])[64]  // output effective potentials,  data layout Veff[Noco^2][nreq][64]
        , double const (*const Vinp)[64]  //  input effective potentials,  data layout Vinp[ncols*Noco^2][64]
        , RequestList_t const & requests
        , int const Noco // =1, 1:no spin, 2: (non-collinear) spin
        , int const echo // =0, log-level
    ) {
        if (echo > 0) std::printf("# new MPI exchange of potential, MPI one-sided communication, Noco=%d\n", Noco);
        assert(1 == Noco || 2 == Noco);

        assert(Veff && "may not be called with a nullptr for output");
        for (int spin = 0; spin < Noco*Noco; ++spin) {
            assert(Veff[spin] && "may not be called with a nullptr for spin output");
        } // spin
        assert(Vinp && "may not be called with a nullptr for input");

        auto const nreq = requests.size();
        view3D<double> Vout(nreq,Noco*Noco,64, 0.0); // get temporary CPU memory in [Noco^2][64] layout

        auto const status = exchange(Vout(0,0), Vinp[0], requests, Noco*Noco*64, echo, "potential");

        // convert Vout into special data layout of Veff (GPU memory) 
        for (size_t ireq = 0; ireq < nreq; ++ireq) {
            for (int spin = 0; spin < Noco*Noco; ++spin) {
                set(Veff[spin][ireq], 64, Vout(ireq,spin));
            } // spin
        } // ireq

        return status;
    } // potential_exchange with RequestList_t


















#ifdef    NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    status_t test_exchanges(int echo=0, int const nSHO=20) {
        status_t stat(0);
        uint32_t const nb[] = {2, 2, 2}; // grid box is 2x2x2 or 8 atoms
        std::vector<int64_t> requests = {0,7,6,1,5,2,4,3}; // every process requests all of these ids

        auto const nrows = requests.size();
        int const me = mpi_parallel::rank(),
                    np = mpi_parallel::size(); assert(np > 0);

        std::vector<int64_t> offerings(0);
        uint32_t const nall[] = {nb[2]*nb[1]*nb[0], 0, 0};
        std::vector<uint16_t> owner_rank(nall[0], 0);
        for(int id{0}; id < nall[0]; ++id) {
            int const rank = id % np; // block-cyclic distribution
            owner_rank[id] = rank;
            if (me == rank) offerings.push_back(id);
        } // id
        int const ncols = offerings.size();
        RequestList_t rlV(requests, offerings, owner_rank.data(), nb, echo);
        RequestList_t rlD(requests, offerings, owner_rank.data(), nall, echo);

        double (*pot_out[2*2])[64]; view3D<double> ptr(4,nrows,64);
        for (int spin = 0; spin < 2*2; ++spin) { pot_out[spin] = (double(*)[64]) ptr(spin,0); }

        int nerr{0};
        for (int Noco = 1; Noco <= 2; ++Noco) {
            {
                auto const pot_inp = new double[nrows*Noco*Noco][64];
                for (int col = 0; col < ncols; ++col) pot_inp[col*Noco*Noco][0] = 0.5 + me; // ear-mark with owner rank
                stat += green_parallel::potential_exchange(pot_out, pot_inp, rlV, Noco, echo);
                for (int row = 0; row < nrows; ++row) nerr += (pot_out[0][row][0] != (0.5 + owner_rank[requests[row]]));
                if (echo > 0) std::printf("# potential_exchange Noco= %d status= %d errors= %d\n\n", Noco, int(stat), nerr);
                delete[] pot_inp;
            }
            {
                int const count = Noco*Noco*2*nSHO*nSHO; // number of doubles per package
                std::vector<double> mat_out(nrows*count), mat_inp(ncols*count);
                for (int col = 0; col < ncols; ++col) mat_inp[col*count] = 0.5 + me; // ear-mark with owner rank
                stat += green_parallel::exchange(mat_out.data(), mat_inp.data(), rlD, count, echo);
                for (int row = 0; row < nrows; ++row) { nerr += (mat_out[row*count] != (0.5 + owner_rank[requests[row]])); }
                if (echo > 0) std::printf("# rank#%i exchange Noco= %d status= %d errors= %d\n\n", me, Noco, int(stat), nerr);
            }
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
