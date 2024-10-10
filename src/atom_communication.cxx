// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <cstdint> // int32_t, uint32_t

#include "status.hxx" // status_t
#include "control.hxx" // ::get
#include "mpi_parallel.hxx" // ::init, ::finalize, MPI_COMM_WORLD, ::barrier
#include "data_view.hxx" // view2D<>, view3D<>
#include "simple_stats.hxx" // ::Stats<>
#include "simple_timer.hxx" // SimpleTimer
#include "recorded_warnings.hxx" // warn
#include "inline_math.hxx" // set, add_product
#include "data_list.hxx" // data_list<T>
#include "atom_communication.hxx" // ::AtomCommList_t

namespace atom_communication {

    AtomCommList_t::AtomCommList_t(
              size_t const n_all_atoms
            , std::vector<uint32_t> const & global_atom_ids // [natoms], global ids of contributing atoms on this rank
            , MPI_Comm const comm // =MPI_COMM_WORLD
            , int const echo // =0
        ) {
            auto const nprocs = mpi_parallel::size(comm); assert(nprocs > 0);
            auto const me     = mpi_parallel::rank(comm); assert(me < nprocs);
            int const na = (n_all_atoms + nprocs - 1 - me)/nprocs; // simple model, owner rank = global_atom_id % nprocs
            int const na_max =  (n_all_atoms + nprocs - 1)/nprocs; // largest  of all numbers of local atoms
            int const na_min =               (n_all_atoms)/nprocs; // smallest of all numbers of local atoms
            if (echo > 9) std::printf("# rank#%i has %d (min %d max %d) owned atoms\n", me, na, na_min, na_max);
            int const na_max8 = (na_max + 7) >> 3; // number of Bytes needed to store at least na_max bits
            if (echo > 5) std::printf("# %s: use MPI_Alltoall with %d--%d bits in %d Byte\n", __func__, na_min, na_max, na_max8);
            comm_ = comm;
            natoms_ = global_atom_ids.size(); // number of contributing atoms
            list_.resize(na);
#ifdef    HAS_NO_MPI
            for (int ia{0}; ia < na; ++ia) { list_[ia].resize(1, 0); } // only rank zero contributes
#else  // HAS_NO_MPI
            for (int ia{0}; ia < na; ++ia) { list_[ia].resize(0); } // init

            view3D<uint8_t> bits(2, nprocs, na_max8, uint8_t(0)); // group 8 atom-process pairings into 1 Byte, total size: ~ n_all_atoms/4 Byte
            auto bits_send = bits[0], bits_recv = bits[1];
            for (auto const & global_atom_id : global_atom_ids) {
                auto const atom_owner = global_atom_id % nprocs; // atom owner rank
                auto const ia =         global_atom_id / nprocs; // local index in atom owner process
                bits_send(atom_owner,ia >> 3) |= (uint8_t(1) << (ia & 7)); // set bit #ia
                //   bits(0, atom_owner, ia) = 1; (if we used a view3D<bool>(2, nprocs, na_max) array)
            } // global_atom_id

            auto const stat = MPI_Alltoall(bits_send[0], na_max8, MPI_UINT8_T, bits_recv[0], na_max8, MPI_UINT8_T, comm);
            if (stat != 0) warn("MPI_Alltoall failed with status= %i", int(stat));
            // Alternative: first use Alltoall for the number of atoms, then Alltoallv for the atom_ids or ias

            for (uint32_t rank{0}; rank < nprocs; ++rank) {
                for (int ia{0}; ia < na; ++ia) {
                    bool const rank_contributes = (bits_recv(rank,ia >> 3) >> (ia & 7)) & 1;
                    if (rank_contributes) {
                        list_[ia].push_back(rank);
                    }
                } // ia
                for (int ia{na}; ia < na_max8*8; ++ia) {
                    bool const rank_contributes = (bits_recv(rank,ia >> 3) >> (ia & 7)) & 1;
                    assert(!rank_contributes && "unused bits must be unset!");
                } // ia
            } // rank

            if (echo > 9) {
                for (int ia{0}; ia < na; ++ia) {
                    auto const global_atom_id = ia*nprocs + me;
                    std::printf("# rank#%i communicates with %ld ranks for atom#%i global#%i\n", 
                                        me, list_[ia].size(), ia, global_atom_id);
                } // ia
            } // echo
#endif // HAS_NO_MPI

        } // constructor

    status_t atom_data_broadcast(
          data_list<double> & atom_data // result [natoms]
        , data_list<double> const & owner_data // input [na], only accessed in atom owner rank
        , char const *const what
        , AtomCommList_t const & atom_comm_list
        , std::vector<uint32_t> const & global_atom_ids
        , int const echo // =0 // log level
    ) {
        // send atom data updated by the atom owner to contributing MPI ranks
        status_t stat(0);
        auto const comm = atom_comm_list.comm();
        auto const nprocs = mpi_parallel::size(comm); assert(nprocs > 0);
        auto const me     = mpi_parallel::rank(comm); assert(me < nprocs);

        mpi_parallel::barrier(comm);

        // atom owners send their data
        auto const na = owner_data.nrows();
        auto const list = atom_comm_list.list();
        assert(na == list.size());
#ifdef    HAS_NO_MPI
        bool const remote_atom_is_error = (0 == control::get("mpi.fake.size", 0.));
#else  // HAS_NO_MPI
        for (int ia{0}; ia < na; ++ia) { // loop over owned atoms
            auto const list_ia = list.at(ia);
            auto const count = owner_data.ncols(ia);
            for (auto const rank : list_ia) {
                if (rank != me) {
                    if (echo > 13) std::printf("# rank#%i %s: send %s, %d doubles for my owned atom#%i, global#%i to contributing rank#%i\n",
                                                       me, __func__, what, count, ia, ia*nprocs + me, rank);
                    MPI_Request send_request;
                    stat += MPI_Isend(owner_data[ia], count, MPI_DOUBLE, rank, ia, comm, &send_request);
                } // remote
            } // rank
        } // ia
        uint32_t irequest{0};
#endif // HAS_NO_MPI


        // contributing atoms receive the data
        uint32_t const natoms = global_atom_ids.size();
        if (echo > 8) std::printf("# %s of %s, %d owned atoms to %d atoms\n", __func__, what, na, natoms);
        assert(natoms == atom_comm_list.natoms());
        std::vector<MPI_Request> recv_requests(natoms, MPI_REQUEST_NULL);
        for (uint32_t iatom{0}; iatom < natoms; ++iatom) { // loop over contributing atoms
            auto const global_atom_id = global_atom_ids[iatom];
            auto const atom_owner = global_atom_id % nprocs;
            auto const ia         = global_atom_id / nprocs;
            int const count = atom_data.ncols(iatom);
            if (atom_owner == me) {
                assert(count == owner_data.ncols(ia));
                if (echo > 11) std::printf("# rank#%i %s: copy %s, %d doubles for owned atom#%i to contributing atom#%i, global#%i, owner rank#%i\n",
                                                   me, __func__, what, count, ia, iatom, global_atom_id, atom_owner);
                set(atom_data[iatom], count, owner_data[ia]); // local copy
            } else {
#ifdef    HAS_NO_MPI
                if (remote_atom_is_error) error("cannot operate remote atoms without MPI, iatom= %i", iatom);
#else  // HAS_NO_MPI
                if (echo > 11) std::printf("# rank#%i %s: recv %s, %d doubles for owned atom#%i from owner rank#%i to contributing atom#%i, global#%i\n",
                                                   me, __func__, what, count, ia, atom_owner, iatom, global_atom_id);
                stat += MPI_Irecv(atom_data[iatom], count, MPI_DOUBLE, atom_owner, ia, comm, &recv_requests[irequest]);
                ++irequest;
#endif // HAS_NO_MPI
            }
        } // iatom

#ifndef   HAS_NO_MPI
        auto const nrequests = irequest;
        std::vector<MPI_Status> statuses(nrequests);
        stat += MPI_Waitall(nrequests, recv_requests.data(), statuses.data());
#endif // HAS_NO_MPI

        mpi_parallel::barrier(comm);

        if (stat) warn("failed for %s with status= %i", what, int(stat));
        return stat;
    } // atom_data_broadcast




    status_t atom_data_allreduce(
          data_list<double> & owner_data // result [na], only correct in atom owner rank
        , data_list<double> const & atom_data // input [natoms]
        , char const *const what
        , AtomCommList_t const & atom_comm_list
        , std::vector<uint32_t> const & global_atom_ids
        , double const factor // =1 // scaling factor, usually g.dV(), the grid volume element
        , int const echo // =0 // log level
    ) {
        // collect the unrenormalized atom_vlm data from contributing MPI ranks
        status_t stat(0);
        auto const comm = atom_comm_list.comm();
        auto const nprocs = mpi_parallel::size(comm);
        auto const me     = mpi_parallel::rank(comm);

        mpi_parallel::barrier(comm);

        // initialize the accumulators
        auto const na = owner_data.nrows();
        for (int ia{0}; ia < na; ++ia) {
            set(owner_data[ia], owner_data.ncols(ia), 0.0); // clear
        } // ia

        // contributing atoms send data
        uint32_t const natoms = global_atom_ids.size();
        if (echo > 8) std::printf("# %s of %s, %d atoms to %d owned atoms\n", __func__, what, natoms, na);
        assert(natoms == atom_comm_list.natoms());
        for (uint32_t iatom{0}; iatom < natoms; ++iatom) { // loop over contributing atoms
            auto const global_atom_id = global_atom_ids[iatom];
            auto const atom_owner = global_atom_id % nprocs;
            auto const ia         = global_atom_id / nprocs;
            int const count = atom_data.ncols(iatom);
            if (atom_owner == me) {
                assert(count == owner_data.ncols(ia));
                if (echo > 11) std::printf("# rank#%i %s:  add %s, %d doubles for owned atom#%i to contributing atom#%i, global#%i, owner rank#%i\n",
                                                   me, __func__, what, count, ia, iatom, global_atom_id, atom_owner);
                add_product(owner_data[ia], count, atom_data[iatom], factor); // local accumulation
            } else {
#ifdef    HAS_NO_MPI
                error("cannot operate remote atoms without MPI, iatom= %i", iatom);
#else  // HAS_NO_MPI
                if (echo > 11) std::printf("# rank#%i %s: send %s, %d doubles for contributing atom#%i, global#%i to owned atom#%i at owner rank#%i\n",
                                                   me, __func__, what, count, iatom, global_atom_id, ia, atom_owner);
                MPI_Request send_request;
                stat += MPI_Isend(atom_data[iatom], count, MPI_DOUBLE, atom_owner, ia, comm, &send_request);
#endif // HAS_NO_MPI
            }
        } // iatom

        // atom owners receive and collect the data
        auto const list = atom_comm_list.list();
        assert(na == list.size());
#ifndef   HAS_NO_MPI
        std::vector<MPI_Request> recv_requests(1);
        for (int ia{0}; ia < na; ++ia) { // loop over owned atoms
            auto const list_ia = list.at(ia);
            auto const count = owner_data.ncols(ia);
            std::vector<double> contrib(count);
            for (auto const rank : list_ia) {
                if (rank != me) {
                    if (echo > 13) std::printf("# rank#%i %s: recv %s, %d doubles for my owned atom#%i, global#%i from contributing rank#%i\n",
                                                       me, __func__, what, count, ia, ia*nprocs + me, rank);
                    MPI_Status status;     // Mind that this is a blocking communication routine
                    stat += MPI_Recv(contrib.data(), count, MPI_DOUBLE, rank, ia, comm, &status);
                    add_product(owner_data[ia], count, contrib.data(), factor); // accumulation
                } // remote
            } // rank
        } // ia
#endif // HAS_NO_MPI

        mpi_parallel::barrier(comm);

        if (stat) warn("failed with status= %i", int(stat));
        return stat;
    } // atom_data_allreduce


#ifdef  NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

    status_t all_tests(int const echo) {
        status_t stat(0);
        auto const already_initialized = mpi_parallel::init();
        if (echo > 0) std::printf("# no tests implemented yet!\n");
        if (!already_initialized) mpi_parallel::finalize();
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace atom_communication
