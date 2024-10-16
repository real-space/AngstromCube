// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, ::fflush, stdout
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <algorithm> // std::min, ::max
#include <cmath> // std::sqrt, ::abs, ::exp, ::erf
#include <type_traits> // std::is_same

#include "parallel_poisson.hxx"

#include "real_space.hxx" // ::grid_t
#include "data_view.hxx" // view2D<T>, view3D<T>
#include "inline_math.hxx" // set, dot_product, pow2, add_product
#include "constants.hxx" // ::pi, ::sqrtpi
#include "boundary_condition.hxx" // Periodic_Boundary
#include "mpi_parallel.hxx" // MPI_Comm, MPI_COMM_WORLD, MPI_COMM_NULL, ::sum, ::rank, ::size, ::min, ::barrier
#include "load_balancer.hxx" // ::get, ::no_owner
#include "simple_timer.hxx" // strip_path
#include "global_coordinates.hxx" // ::get
#include "boundary_condition.hxx" // Isolated_Boundary, Periodic_Boundary
#include "green_parallel.hxx" // ::exchange, ::RequestList_t
#include "print_tools.hxx" // printf_vector
#include "control.hxx" // ::get
#include "recorded_warnings.hxx" // error

#ifndef   NO_UNIT_TESTS
  #include <array> // std::array<T,N>
  #include <algorithm> // std::stable_sort

  #include "fourier_poisson.hxx" // ::solve
  #include "lossful_compression.hxx" // print_compressed
#endif // NO_UNIT_TESTS

namespace parallel_poisson {
    // solve the Poisson equation iteratively using the conjugate gradients method,
    // 16th-order finite differences, 8x8x8 grid point grouping and MPI data exchange

    double constexpr m1over4pi = -.25/constants::pi; // -4*constants::pi is the electrostatics prefactor in Hartree atomic units

    template <typename real_t>
    double norm2(real_t const v[], size_t const n, MPI_Comm const comm=MPI_COMM_NULL) {
        double s{0};
        for (size_t i{0}; i < n; ++i) { 
            s += pow2(v[i]);
        } // i 
        if (MPI_COMM_NULL != comm) mpi_parallel::sum(&s, 1, comm);
        return s;
    } // norm2

    template <typename real_t>
    double norm1(real_t const v[], size_t const n, MPI_Comm const comm=MPI_COMM_NULL, int const echo=0) {
        double s{0};
        for (size_t i{0}; i < n; ++i) {
            s += v[i];
        } // i
        auto const norm1_local = s;
        if (MPI_COMM_NULL != comm) mpi_parallel::sum(&s, 1, comm);
        if (echo > 0) std::printf("# norm1_local= %g, norm1= %g\n", norm1_local, s);
        return s;
    } // norm1

    template <typename real_t>
    double scalar_product(real_t const v[], real_t const w[], size_t const n, MPI_Comm const comm=MPI_COMM_NULL) {
        double dot{0};
        for (size_t i = 0; i < n; ++i) {
            dot += double(v[i])*double(w[i]); // conversion to double is different from dot_product define in inline_math.hxx
        } // i
        if (MPI_COMM_NULL != comm) mpi_parallel::sum(&dot, 1, comm);
        return dot;
    } // scalar_product


    load_balancing_t::load_balancing_t(
        real_space::grid_t const & g // grid descriptor of the entire grid
      , MPI_Comm const comm // MPI communicator
      , unsigned const n8 // number of grid points per block edge
      , int const echo // =0 log-level
    ) { // constructor

        comm_ = comm; // copy the communicator
        int32_t const me = mpi_parallel::rank(comm);
        uint32_t const np = control::get("parallel_poisson.nprocs", double(mpi_parallel::size(comm))); // listens to fake MPI setups

        auto nb = nb_;
        auto const ng = g.grid_points();
        if (echo > 8) std::printf("# %s(%d x %d x %d grid points in blocks of %dx%dx%d)\n", __func__, ng[0], ng[1], ng[2], n8,n8,n8);
        assert(n8 > 0);
        for (int d{0}; d < 3; ++d) {
            nb[d] = ng[d]/n8; // divide by n8
            assert(nb[d]*n8 == ng[d] && "grid numbers must be a positive multiple of the block edge");
            assert(nb[d] > 0 && "at least one block of grid points needed");
        } // d
        if (echo > 3) std::printf("# %s(%d x %d x %d grid points in %d x %d x %d blocks of %dx%dx%d)\n", 
                                __func__, ng[0], ng[1], ng[2], nb[0], nb[1], nb[2],        n8,n8,n8);

        owner_rank_ = view3D<green_parallel::rank_int_t>(nb[2], nb[1], nb[0], load_balancer::no_owner);

        double rank_center[4] = {0,0,0,  0};

        load_ = load_balancer::get(np, me, nb, echo, rank_center, owner_rank_.data());
        n_local_cubes_ = rank_center[3]; // the 4th component contains the number of items
        if (echo > 7) std::printf("# rank#%i rank center %g %g %g\n", me, rank_center[0], rank_center[1], rank_center[2]);

        if (true) {
            simple_stats::Stats<double> load_stats(0);
            load_stats.add(load_);
            mpi_parallel::allreduce(load_stats, comm);
            if (echo > 2) std::printf("# load_balancer distribution over %d ranks is %s\n", np, load_stats.interval().c_str());
        } // true

        // load_balancer has the mission to produce coherent domains with not too lengthy extents in space
        auto const nall = size_t(nb[2])*size_t(nb[1])*size_t(nb[0]);
        if (MPI_COMM_NULL != comm) {
            if (echo > 19) {
                std::printf("# rank#%i owner_rank before MPI_MIN ", me);
                printf_vector(" %i", owner_rank_.data(), nall);
            } // echo
            mpi_parallel::min(owner_rank_.data(), nall, comm);
            if (echo > 19) {
                std::printf("# rank#%i owner_rank after  MPI_MIN ", me);
                printf_vector(" %i", owner_rank_.data(), nall);
            } // echo
            if (echo > 2) {
                simple_stats::Stats<> ors;
                auto const owner_rank_data = owner_rank_.data();
                for (size_t iall{0}; iall < nall; ++iall) {
                    ors.add(owner_rank_data[iall]);
                } // iall
                std::printf("# rank#%i owner ranks in %s\n", me, ors.interval().c_str());
            } // echo
        } // not comm_null
        mpi_parallel::barrier(comm);

        set(min_domain_, 3, int32_t((1ull << 31) - 1)); 
        set(max_domain_, 3, -1);
        double dom_center[] = {0, 0, 0};
        size_t nown{0};
        for (uint32_t iz = 0; iz < nb[2]; ++iz) {
        for (uint32_t iy = 0; iy < nb[1]; ++iy) {  // this triple loop does not scale well as it is the same range for all processes
        for (uint32_t ix = 0; ix < nb[0]; ++ix) {
            auto const owner_rank_xyz = owner_rank_(iz,iy,ix);
            if (load_balancer::no_owner == owner_rank_xyz) {
                if (echo > 9) {
                    std::printf("# rank#%i owner_rank after  MPI_MIN ", me);
                    printf_vector(" %i", owner_rank_.data(), nall);
                    std::fflush(stdout);
                } // echo
                error("rank#%i entry[ix=%d,iy=%d,iz=%d] has no owner rank", me, ix,iy,iz);
            } else
            if (me == owner_rank_xyz) {
                ++nown;
                int32_t const ixyz[] = {int32_t(ix), int32_t(iy), int32_t(iz)};
                for (int d = 0; d < 3; ++d) {
                    min_domain_[d] = std::min(min_domain_[d], ixyz[d]);
                    max_domain_[d] = std::max(max_domain_[d], ixyz[d]);
                    dom_center[d] += (ixyz[d] + 0.5);
                } // d
            } // my
        }}} // iz iy ix
        if (nown != n_local_cubes_) {
            warn("expected match between n_local_blocks= %d and count(owner_rank[]==me)= %ld", n_local_cubes_, nown);
            n_local_cubes_ = nown;
        }
        if (echo > 5) std::printf("# rank#%i %s: load_balancer::get = %g, %g items, %d local blocks\n",
                                            me, __func__, load_, rank_center[3], n_local_cubes_);
        auto const by_nown = nown ? 1./nown : 0;
        for (int d = 0; d < 3; ++d) {
            dom_center_[d] = dom_center[d] * by_nown;
            assert(max_domain_[d] >= min_domain_[d]);
        } // d
        if (echo > 7) std::printf("# rank#%i domain center %g %g %g\n", me, dom_center_[0], dom_center_[1], dom_center_[2]);
        if (echo > 7) std::printf("# rank#%i   rank center %g %g %g\n", me, rank_center[0], rank_center[1], rank_center[2]);

    } // load_balancing_t constructor






    parallel_grid_t::parallel_grid_t(
        real_space::grid_t const & g // grid descriptor of the entire grid
      , load_balancing_t const & lb
      , int const echo // =0 log-level
      , char const *const what // ="FD1"
    ) { // constructor

        comm_ = lb.comm(); // copy the communicator
        int32_t const me = mpi_parallel::rank(comm_);
        auto nb = nb_;
        set(nb, 3, lb.grid_cubes());

        auto const bc = g.boundary_conditions();
        if (echo > 3) std::printf("# %s(nb= [%d %d %d], nall= %ld, bc=[%d %d %d], what=%s)\n", __func__,
            nb[0], nb[1], nb[2], nb[2]*size_t(nb[1])*size_t(nb[0]), bc[0], bc[1], bc[2], what);

        for (int d{0}; d < 3; ++d) {
            bc_[d] = bc[d]; // copy
            assert(0 != g.h[d]);
            h2_[d] = 1./pow2(g.h[d]);
        } // d
        dVol_ = g.dV();
        nperiodic_ = g.number_of_boundary_conditions(Periodic_Boundary);

        int nstencil_; // non-const version
        int8_t const (*stencil_)[4] {nullptr}; // non-const version

        int8_t const stencil_FD1[6][4] = { {-1,0,0,  1}, {1,0,0,  1},
                                           {0,-1,0,  1}, {0,1,0,  1},
                                           {0,0,-1,  1}, {0,0,1,  1} };

        int8_t stencil_333[27][4];
        if (*what == 'I') { // 3x3x3 Interpolation stencil
            nstencil_ = 27;
            stencil_ = stencil_333;
            for (int z = -1; z <= 1; ++z) {
            for (int y = -1; y <= 1; ++y) {
            for (int x = -1; x <= 1; ++x) {
                auto const k = ((z + 1)*3 + (y + 1))*3 + (x + 1);
                stencil_333[k][0] = x;
                stencil_333[k][1] = y;
                stencil_333[k][2] = z;
                stencil_333[k][3] = x*x + y*y + z*z; // distance^2
            }}} // x y z
        } else {
            // 3D-Finite Difference Star-shaped stencil
            stencil_ = stencil_FD1;
            nstencil_ = 6;
        }
        int8_t const (*const stencil)[4] = stencil_;
        int const nstencil = nstencil_;
        assert(nstencil > 0);


        uint32_t const n_local_blocks = lb.n_local();
        auto const & owner_rank = lb.owner_rank();
        local_global_ids_.resize(0);

        if (n_local_blocks > 0) { // scope: setup of star and remote_global_ids, determination of n_remote_blocks

            auto const *const max_domain = lb.max_domain();
            auto const *const min_domain = lb.min_domain();

            int32_t constexpr HALO=1;
            int32_t const ndom[] = {max_domain[0] - min_domain[0] + 1 + 2*HALO,
                                    max_domain[1] - min_domain[1] + 1 + 2*HALO,
                                    max_domain[2] - min_domain[2] + 1 + 2*HALO};
            int32_t const ioff[] = {min_domain[0] - HALO, min_domain[1] - HALO, min_domain[2] - HALO};
            if (echo > 7) std::printf("# rank#%i domain(%d %d %d), halo= %d\n", me, ndom[0], ndom[1], ndom[2], HALO);
            assert(ndom[0] > 0); assert(ndom[1] > 0); assert(ndom[2] > 0);

            uint8_t constexpr OUTSIDE=0, BORDER=1, INSIDE=2, VACUUM=7;
            view3D<uint8_t> domain(ndom[2],ndom[1],ndom[0], OUTSIDE);

            if (echo > 7) { std::printf("# rank#%i here %s:%d\n", me, strip_path(__FILE__), __LINE__); std::fflush(stdout); }

            view3D<int32_t> domain_index(ndom[2],ndom[1],ndom[0], -1); // -1: not assigned, <n_local_blocks: INSIDE

            if (echo > 7) { std::printf("# rank#%i here %s:%d\n", me, strip_path(__FILE__), __LINE__); std::fflush(stdout); }
            if (echo > 7) std::printf("# rank#%i n_local_blocks=%d\n", me, n_local_blocks);

            local_global_ids_.resize(n_local_blocks, int64_t(-1));

         // if (echo > 7) { std::printf("# rank#%i here %s:%d\n", me, strip_path(__FILE__), __LINE__); std::fflush(stdout); }

            uint32_t ilb{0}; // index of the local block
            for (int32_t iz = HALO; iz < ndom[2] - HALO; ++iz) {
            for (int32_t iy = HALO; iy < ndom[1] - HALO; ++iy) { // 1st domain loop, serial
            for (int32_t ix = HALO; ix < ndom[0] - HALO; ++ix) {
                // translate domain indices {ix,iy,iz} into box indices
                int32_t const ixyz[] = {ix + ioff[0], iy + ioff[1], iz + ioff[2]};
                for (int d = 0; d < 3; ++d) { assert(ixyz[d] >= 0); assert(ixyz[d] < nb[d]); }
                if (me == owner_rank(ixyz[2],ixyz[1],ixyz[0])) {
                    domain(iz,iy,ix) |= INSIDE;
                    domain_index(iz,iy,ix) = ilb; // domain index for inner elements

                    // mark stencil coverage as border regions
                    for (int k{0}; k < nstencil; ++k) {
                        domain(iz + stencil[k][2], iy + stencil[k][1], ix + stencil[k][0]) |= BORDER;
                    } // k

                    local_global_ids_[ilb] = global_coordinates::get(ixyz);

                    ++ilb; // count inside elements
                } // me == owner
            }}} // iz iy ix
            assert(ilb <= (1ull << 31));
            if (ilb != n_local_blocks) error("expected match between n_local_blocks=%d and count(owner_rank[]==me)=%d\n", n_local_blocks, ilb);

            if (echo > 7) { std::printf("# rank#%i here %s:%d\n", me, strip_path(__FILE__), __LINE__); std::fflush(stdout); }


            // loop over the domain again, this time including the halos
            uint32_t jrb{0}; // index for border-only blocks
            size_t st[4] = {0, 0, 0}; // statistics for display
            for (int32_t iz = 0; iz < ndom[2]; ++iz) {
            for (int32_t iy = 0; iy < ndom[1]; ++iy) { // 2nd domain loop, serial
            for (int32_t ix = 0; ix < ndom[0]; ++ix) {
                auto const dom = domain(iz,iy,ix);
                if (BORDER == dom) { // is border-only
                    domain_index(iz,iy,ix) = n_local_blocks + jrb; // domain index for border elements
                    ++jrb; // count remote elements
                } // border
                ++st[dom & 0x3]; // dom should be in [0, 3] anyway but better safe than sorry
            }}} // iz iy ix
            if (echo > 5) std::printf("# rank#%i has %ld outside, %ld border, %ld inside, %ld inside+border elements\n",
                                                me, st[OUTSIDE], st[BORDER], st[INSIDE], st[INSIDE+BORDER]);
            uint32_t const n_remote_blocks = jrb;
            if (echo > 5) std::printf("# rank#%i has %d local and %d remote blocks, %d in total\n",
                                                me, n_local_blocks, n_remote_blocks, n_local_blocks + n_remote_blocks);
            
            remote_global_ids_.resize(n_remote_blocks, -1); // init remote element request lists
            star_ = view2D<uint32_t>(n_local_blocks, nstencil, uint32_t(-1)); // init finite-difference neighborhood lists

            inner_cell_.resize(n_local_blocks, false);

            size_t vacuum_assigned{0}, inner_cell_found{0};
            // loop over domain again (3rd time), with halos
            for (int32_t iz = HALO; iz < ndom[2] - HALO; ++iz) {
            for (int32_t iy = HALO; iy < ndom[1] - HALO; ++iy) { // 3rd domain loop, could be parallel
            for (int32_t ix = HALO; ix < ndom[0] - HALO; ++ix) {
                // translate domain indices {ix,iy,iz} into box indices
                int32_t const ixyz[] = {ix + ioff[0], iy + ioff[1], iz + ioff[2]};
                for (int d = 0; d < 3; ++d) { assert(ixyz[d] >= 0); assert(ixyz[d] < nb[d]); } // should still hold...
                if (domain(iz,iy,ix) & INSIDE) {
                    auto const id0 = domain_index(iz,iy,ix);
                    assert(id0 >= 0); assert(id0 < n_local_blocks);

                    for (int k{0}; k < nstencil; ++k) {
                        int32_t       jxyz[] = {ixyz[0] + stencil[k][0], ixyz[1] + stencil[k][1], ixyz[2] + stencil[k][2]}; // global coordinates
                        int32_t const jdom[] = {ix      + stencil[k][0], iy      + stencil[k][1], iz      + stencil[k][2]}; // domain coordinates
                        for (int d = 0; d < 3; ++d) { assert(jdom[d] >= 0); assert(jdom[d] < ndom[d]); }

                        assert( BORDER & domain(jdom[2],jdom[1],jdom[0]) ); // consistent with 1st domain loop
                        auto const klb = domain_index(jdom[2],jdom[1],jdom[0]);
                        assert(klb >= 0);
                        star_(id0,k) = klb;

                        int is_vacuum{0};
                        for (int d = 0; d < 3; ++d) {
                            auto const bc_shift = int(jxyz[d] < 0) - int(jxyz[d] >= int32_t(nb[d]));
                            jxyz[d] += bc_shift*int32_t(nb[d]); // fold back into [0, nb)
                            is_vacuum += int((0 != bc_shift) && (Isolated_Boundary == bc[d]));
                        }
                        if (is_vacuum) {
                            // id = -1; // not existing
                            auto & dom = domain(jdom[2],jdom[1],jdom[0]);
                            // mark as vacuum
                            assert(dom != INSIDE);
                            if (VACUUM == dom) {
                                assert(27 == nstencil); // in the case of an interpolation stencil it can already be marked as vacuum before, 
                            } else {                    // in the case of a finite-difference stencil this should not happen
                                assert(BORDER == dom);  // but the it must be a border-only cell
                                dom = VACUUM; // mark in the mask (dom is a reference)
                                ++vacuum_assigned;
                            }
                        } else {
                            if (klb >= n_local_blocks) {
                                assert(BORDER == domain(jdom[2],jdom[1],jdom[0])); // must be a border-only element
                                auto const irb = klb - n_local_blocks;
                                assert(irb >= 0);
                                for (int d = 0; d < 3; ++d) { assert(jxyz[d] >= 0); assert(jxyz[d] < nb[d]); }
                                auto const gid = global_coordinates::get(jxyz);
                                assert((gid == remote_global_ids_[irb]) || (-1 == remote_global_ids_[irb])); // either unassigned(-1) or the same value
                                remote_global_ids_[irb] = gid;
                            } // add to remote list
                        } // is_vacuum
                    } // k
                    if (6 == nstencil && echo > 19) std::printf("# rank#%i star[%i,:] = {%i %i %i %i %i %i}\n", me, id0,
                                    star_(id0,0), star_(id0,1), star_(id0,2), star_(id0,3), star_(id0,4), star_(id0,5));
                    
                    // to be an inner_cell_ all blocks hit by the stencil must be local blocks
                    int is_inner_cell{0};
                    for (int k{0}; k < nstencil; ++k) {
                        is_inner_cell += (star_(id0,k) < n_local_blocks);
                    } // k
                    inner_cell_[id0]  = (is_inner_cell == nstencil);
                    inner_cell_found += (is_inner_cell == nstencil);

                } // is INSIDE
            }}} // iz iy ix

            if (echo > 8) std::printf("# rank#%i star list generated, %ld inner cells found\n", me, inner_cell_found); std::fflush(stdout);

            size_t vacuum_requested{0};
            for (auto id : remote_global_ids_) {
                vacuum_requested += (-1 == id);
            } // id
            if (vacuum_requested && echo > 3) std::printf("# rank#%i assigned %ld, request %ld vacuum cells\n", me, vacuum_assigned, vacuum_requested);
            assert(vacuum_assigned == vacuum_requested);

        } else { // n_local_blocks > 0
            star_ = view2D<uint32_t>(nullptr, 6); // dummy
        } // n_local_blocks > 0

        if (echo > 8) {
            std::printf("# rank#%i %s: requests={", me, __func__);
            for (auto rq : remote_global_ids_) {
                std::printf(" %lli", rq);
            } // rq
            std::printf(" }, %ld items\n", remote_global_ids_.size());

            std::printf("# rank#%i %s: offering={", me, __func__);
            for (auto of : local_global_ids_) {
                std::printf(" %lli", of);
            } // of
            std::printf(" }, %ld items\n", local_global_ids_.size());
        } // echo

        if (echo > 9) { std::printf("# rank#%i waits in barrier at %s:%d nb=%d %d %d\n", me, strip_path(__FILE__), __LINE__, nb[0], nb[1], nb[2]); std::fflush(stdout); }
        mpi_parallel::barrier(comm_);

        requests_ = green_parallel::RequestList_t(remote_global_ids_, local_global_ids_, owner_rank.data(), nb, echo, what);

        if (echo > 8) {
            std::printf("# rank#%i %s: RequestList.owner={", me, __func__);
            for (auto ow : requests_.owner) {
                std::printf(" %i", ow);
            } // ow
            std::printf(" }, %ld items\n", requests_.owner.size());
        } // echo

    } // parallel_grid_t constructor







    template <typename real_t>
    status_t data_exchange(
          real_t *v  // input and result array, data layout v[n_local_remote][count]
        , parallel_grid_t const & pg // descriptor
        , size_t const count // number of real_t per package
        , int const echo=0 // log-level
        , char const *const what="??"
    ) {
        auto const n_local = pg.n_local();
        auto const stat = green_parallel::exchange(v + count*n_local, v, pg.requests(), count, echo, what); // ToDo: inser communicator
        return stat;
    } // data_exchange



    template <typename real_t> // =double
    status_t block_interpolation(
          real_t       *const v888 // result array, data layout v888[n_local_blocks][8*8*8]
        , real_t const *const v444 // input  array, data layout v444[n_local_blocks][4*4*4]
        , parallel_grid_t const & pg // descriptor, must be prepared with "3x3x3"
        , int const echo // =0 // log level
        , double const factor // =1
        , char const *const what // ="!"
    ) {

        auto const nlb = pg.n_local();
        auto const nrb = pg.n_remote();
        if (echo > 9) std::printf("\n# %s start what=%s %d local blocks, %d remote blocks\n", __func__, what, nlb, nrb);
        view2D<real_t> v4(nlb + nrb, 4*4*4, real_t(0));
        set(v4[0], nlb*64, v444); // copy in
        auto const stat = data_exchange(v4.data(), pg, 4*4*4, echo, __func__); // fill remote blocks

        if (echo > 9) std::printf("\n# %s %d remote blocks exchanged\n", __func__, nrb);

        assert(27 == pg.star_dim());
        auto const star = (uint32_t const(*)[27])pg.star();

        double const f = factor/(4*4*4); // interpolation weights are [0.25 0.75] expressed as [1 3]/4 --> denominator 4 per dimension

        for (uint32_t ilb = 0; ilb < nlb; ++ilb) { // loop over local blocks --> CUDA block-parallel
            auto const *const nn = star[ilb]; // 27 nearest-neighbor blocks of block ilb, load into GPU shared memory
            if (echo > 39) { std::printf("# star[%i,:] = ", ilb); printf_vector(" %i", nn, 27); }

            // copy data into a halo=1-enlarged block v666
            real_t v666[6][6][6]; // real_t=float 864 Byte, real_t=double 1.728 kByte
            for (int z = -1; z < 5; ++z) { int const z3 = (z + 2) >> 2;
            for (int y = -1; y < 5; ++y) { int const y3 = (y + 2) >> 2; // mapping of y:[-1,0,1,2,3,4] --> y3:[0,1,1,1,1,2]
            for (int x = -1; x < 5; ++x) { int const x3 = (x + 2) >> 2;
                auto const i3zyx = (z3*3 + y3)*3 + x3;
                assert(i3zyx >= 0); assert(i3zyx < 27);
                auto const i64 = nn[i3zyx];
                assert(i64 < nlb + nrb);
                auto const i4zyx = ((z & 0x3)*4 + (y & 0x3))*4 + (x & 0x3); // & 0x3 is the same as % 4
                v666[z + 1][y + 1][x + 1] = v4(i64,i4zyx)*f;
            }}} // x y z

            // interpolate linearly in x-direction, weights are {1,3}
            double v668[6][6][8]; // 2.304 kByte
            for (int z = 0; z < 6; ++z) {
            for (int y = 0; y < 6; ++y) {
            for (int x = 0; x < 4; ++x) {
                v668[z][y][2*x+0] = v666[z][y][x+0] + 3.*v666[z][y][x+1];
                v668[z][y][2*x+1] = v666[z][y][x+1]*3. + v666[z][y][x+2];
            }}} // x y z

            // interpolate linearly in y-direction, weights are {1,3}
            double v688[6][8][8]; // 3.072 kByte
            for (int z = 0; z < 6; ++z) {
            for (int y = 0; y < 4; ++y) {
            for (int x = 0; x < 8; ++x) {
                v688[z][2*y+0][x] = v668[z][y+0][x] + 3*v668[z][y+1][x];
                v688[z][2*y+1][x] = v668[z][y+1][x]*3 + v668[z][y+2][x];
            }}} // x y z

            // interpolate linearly in z-direction, weights are {1,3}
            auto const i512 = size_t(ilb) << 9; // block offset in v8 blocks, write blocks of 4.096 kByte
            for (int z = 0; z < 4; ++z) {
            for (int y = 0; y < 8; ++y) {
            for (int x = 0; x < 8; ++x) {
                v888[i512 + ((2*z+0)*8 + y)*8 + x] = v688[z+0][y][x] + 3*v688[z+1][y][x];
                v888[i512 + ((2*z+1)*8 + y)*8 + x] = v688[z+1][y][x]*3 + v688[z+2][y][x];
            }}} // x y z

        } // ilb

        if (echo > 9) std::printf("# %s done\n\n", __func__);
        return stat;
    } // block_interpolation

    template // explicit template instantiation for real_t=double
    status_t block_interpolation(double*, double const*, parallel_grid_t const &, int, double, char const*);


    template <typename real_t, typename double_t=double>
    status_t Laplace16th(
          real_t *Av // result array, data layout Av[n_local_blocks][8*8*8]
        , real_t *v  // input  array, data layout  v[n_local_remote][8*8*8], cannot be const due to call data_exchange onto v
        , parallel_grid_t const & pg // descriptor
        , int const echo=0 // log level
        , double const prefactor=1
    ) {
        if (echo > 9) std::printf("\n# %s start\n", __func__);

        // to reduce the latencies, we could start to apply the stencil to inner cells that do not depend on remote data

        // mode must be in {0, 1, 2, 3}, 0: only communicate, 1:only do inner cells, 2: only do outer cells, 3: do both

        auto const stat = data_exchange(v, pg, 8*8*8, echo, __func__);
        double const *const h2 = pg.get_prefactors();

        // prepare finite-difference coefficients (isotropic)
        //            c_0        c_1         c_2       c_3        c_4      c_5       c_6     c_7     c_8
        // FD16th = [-924708642, 538137600, -94174080, 22830080, -5350800, 1053696, -156800, 15360, -735] / 302702400
        double_t const norm = prefactor/302702400.;
        double_t const cFD[9] = {-924708642*norm,
                                  538137600*norm,
                                  -94174080*norm,
                                   22830080*norm,
                                   -5350800*norm,
                                    1053696*norm,
                                    -156800*norm,
                                      15360*norm,
                                       -735*norm};
        int const nlb = pg.n_local();
        assert(6 == pg.star_dim());
        auto const star = (uint32_t const(*)[6])pg.star();
        for (uint32_t ilb = 0; ilb < nlb; ++ilb) { // loop over local blocks --> CUDA block-parallel
            auto const i512 = size_t(ilb) << 9; // block offset
            auto const *const nn = star[ilb]; // nearest-neighbor blocks of block ilb, load into GPU shared memory
            if (echo > 11) std::printf("# Laplace16th: for ilb= %i take from neighbors{%i %i, %i %i, %i %i}\n",
                                                           ilb, nn[0], nn[1],  nn[2], nn[3],  nn[4], nn[5]);
            for (int iz = 0; iz < 8; ++iz) {
            for (int iy = 0; iy < 8; ++iy) { // loops over block elements --> CUDA thread-parallel
            for (int ix = 0; ix < 8; ++ix) {
                auto const izyx = iz*64 + iy*8 + ix;
                auto const i0 = i512 + izyx;
                double_t const av = cFD[0]*double_t(v[i0]);
                double_t ax{av}, ay{av}, az{av}; // accumulators
                // if (echo > 9) std::printf("# Av[%i][%3.3o] init as %g\n", ilb, izyx, av);
                for (int ifd = 1; ifd <= 8; ++ifd) {
                    // as long as ifd is small enough, we take from the central block of v, otherwise from neighbor blocks
                    auto const ixm = (ix >=    ifd) ? i0 - ifd : (nn[0] << 9) + izyx + 8 - ifd;
                    auto const ixp = (ix + ifd < 8) ? i0 + ifd : (nn[1] << 9) + izyx - 8 + ifd;
                    ax += cFD[ifd]*(double_t(v[ixm]) + double_t(v[ixp]));
                    auto const iym = (iy >=    ifd) ? i0 - 8*ifd : (nn[2] << 9) + izyx + 64 - 8*ifd;
                    auto const iyp = (iy + ifd < 8) ? i0 + 8*ifd : (nn[3] << 9) + izyx - 64 + 8*ifd;
                    ay += cFD[ifd]*(double_t(v[iym]) + double_t(v[iyp]));
                    auto const izm = (iz >=    ifd) ? i0 - 64*ifd : (nn[4] << 9) + izyx + 512 - 64*ifd;
                    auto const izp = (iz + ifd < 8) ? i0 + 64*ifd : (nn[5] << 9) + izyx - 512 + 64*ifd;
                    az += cFD[ifd]*(double_t(v[izm]) + double_t(v[izp]));
                    // if (echo > 9) std::printf("# %d += x[%i][%3.3o] + x[%i][%3.3o] + y[%i][%3.3o] + y[%i][%3.3o] + z[%i][%3.3o] + z[%i][%3.3o]\n", ifd,
                    //                  ixm>>9, ixm&511, ixp>>9, ixp&511, iym>>9, iym&511, iyp>>9, iyp&511, izm>>9, izm&511, izp>>9, izp&511);
                } // ifd
                Av[i0] = ax*h2[0] + ay*h2[1] + az*h2[2]; // store, possible conversion from double_t to real_t
            }}} // ix iy iz
        } // ilb

        if (echo > 9) std::printf("# %s done\n\n", __func__);
        return stat;
    } // Laplace16th


    template <typename real_t>
    status_t solve(
          real_t xx[] // result to Laplace(x)/(-4*pi) == b, only rank-local cubes, data layout x[][8*8*8]
        , real_t const bb[] // right hand side b          , only rank-local cubes, data layout b[][8*8*8]
        , parallel_grid_t const & pg // parallel grid descriptor
        , char const method // ='c' // use mixed precision as preconditioner
        , int const echo // =0 // log level
        , float const threshold // =3e-8 // convergence criterion
        , float *residual // =nullptr // residual that was reached
        , int const maxiter // =999 // maximum number of iterations 
        , int const miniter // =0   // minimum number of iterations
        , int restart // =4096 // number of iterations before restart, 1:steepest descent
        , double *inner_xx_bb // = nullptr
    ) {

        auto const comm = pg.comm();
        int const echo_L = echo >> 3; // verbosity of Lapacian16th

        auto const nb = pg.grid_cubes();
        size_t const n_all_grid_points = size_t(nb[2]*8)*size_t(nb[1]*8)*size_t(nb[0]*8);
        auto const nall = pg.n_local()* size_t(512),
                   nrem = pg.n_remote()*size_t(512);

        status_t ist(0);

        restart = ('s' == method) ? 1 : std::max(1, restart);

        if (std::is_same<real_t,double>::value) {
            view2D<float> xb(2, nall); // get memory
            auto const x32=xb[0], b32=xb[1];
            set(b32, nall, bb); // convert to float
            set(x32, nall, xx); // convert to float
            if (echo > 5) std::printf("# %s solve in <float> precision first\n", strip_path(__FILE__));
            ist += solve(x32, b32, pg, method, echo, threshold, residual, maxiter, miniter, restart);
            if (echo > 5) std::printf("# %s switch back to <double> precision\n", strip_path(__FILE__));
            set(xx, nall, x32); // convert to double
        } // real_t == double

        // we use CG + order-16 FD Laplace operator, no preconditioning
        bool constexpr use_precond = false;

        // find memory aligned nloc
        view2D<real_t> mem(6 + use_precond, nall + nrem, 0.0); // get memory
        auto const x=mem[0], r=mem[1], p=mem[2], ax=mem[3], ap=mem[4], b=mem[5], z=use_precond?mem[6]:r; 
        // ToDo: need 3x (nall+nrem) for x, p, r
        //       need 4x (nall)      for b, ax, ap, z
        
        set(x, nall, xx);
        set(b, nall, bb);

        double const cell_volume = n_all_grid_points * pg.dV();
        double const threshold2 = cell_volume * pow2(threshold);
        double constexpr RZ_TINY = 1e-14, RS_TINY = 1e-10;

    
//     ! |Ax> = A|x>
//     ! |r> = |b> - |Ax>
//     ! |z> = P|r>
//     ! rz_old = <r|z>
//     ! |p> = |z>
//     ! it = 0
//     ! do while
//     !   |Ap> = A|p>
//     !   pAp = <p|Ap>
//     !   alpha = rz_old / pAp
//     !   |x> = |x> + alpha |p>
//     !   |r> = |r> - alpha |Ap>
//     !   res = <r|r>
//     !   |z> = P|r>
//     !   rz_new = <r|z>
//     !   beta = rz_new / rz_old
//     !   rz_old = rz_new
//     !   |p> = |z> + beta |p>
//     !   it = it+1


        // |Ax> := A|x>
        ist = Laplace16th(ax, x, pg, echo_L, m1over4pi);
        if (ist) error("CG_solve: Laplacian failed with status %i", int(ist));
        
        if (pg.all_periodic_boundary_conditions()) {
            double const bnorm = norm1(b, nall, comm) * pg.dV();
            if (echo > 8) std::printf("# %s all boundary conditions are periodic but system is charged with %g electrons\n", strip_path(__FILE__), bnorm);
        } // all_boundary_conditions_periodic

        // |r> = |b> - A|x> = |b> - |Ax>
        set(r, nall, b); add_product(r, nall, ax, real_t(-1));

        // res^2 = <r|r>
        double res2 = norm2(r, nall, comm) * pg.dV();
        double const res_start = std::sqrt(res2/cell_volume); // store starting residual
        if (echo > 8) std::printf("# %s start residual=%.1e\n", strip_path(__FILE__), res_start);

        // |z> = |Pr> = P|r>
        if (use_precond) {
            error("CG_solve: Preconditioner deactivated in line %i", __LINE__);
        } else assert(z == r);

        // rz_old = <r|z>
        double rz_old = scalar_product(r, z, nall, comm) * pg.dV();

        // |p> = |z>
        set(p, nall, z);

        int it{0}; // init iteration counter

        // number of iterations is less then maxiter?
        bool run = (it < maxiter);
        while (run) {
            ++it;
//       !--------------------------------------
//       ! begin of the CG iteration
//       !--------------------------------------

            // |ap> = A|p>
            ist = Laplace16th(ap, p, pg, echo_L, m1over4pi);
            if (ist) error("CG_solve: Laplacian failed with status %i", int(ist));

            double const pAp = scalar_product(p, ap, nall, comm) * pg.dV();

            // alpha = rz_old / pAp
            double const alpha = (std::abs(pAp) < RZ_TINY) ? RZ_TINY : rz_old / pAp;

            // |x> = |x> + alpha |p>
            add_product(x, nall, p, real_t(alpha));

//       !============================================================
//       ! special treatment of completely periodic case
//       !============================================================
            if (pg.all_periodic_boundary_conditions()) {
                real_t const xnorm = norm1(x, nall, comm)/n_all_grid_points;
                // subtract the average potential
                for (size_t i{0}; i < nall; ++i) { x[i] -= xnorm; }
            } // 3 periodic BCs
//       !============================================================

            if (0 == (it % restart)) {
                // |Ax> = A|x> for restart
                ist = Laplace16th(ax, x, pg, echo_L, m1over4pi);
                if (ist) error("CG_solve: Laplacian failed with status %i", int(ist))
                // |r> = |b> - A|x> = |b> - |ax>
                set(r, nall, b);
                add_product(r, nall, ax, real_t(-1));
            } else {
                // |r> = |r> - alpha |ap>
                add_product(r, nall, ap, real_t(-alpha));
            } // restart

            // res = <r|r>
            res2 = norm2(r, nall, comm) * pg.dV();

            // |z> = |Pr> = P|r>
            if (use_precond) {
                error("CG_solve: Preconditioner deactivated in line %i", __LINE__);
            } else assert(z == r);

            // rz_new = <r|z>
            double const rz_new = scalar_product(r, z, nall, comm) * pg.dV();

            // beta = rz_new / rz_old
            double beta{0};
            if (rz_old < RS_TINY) {
                set(p, nall, z); // beta == 0
            } else {
                beta = rz_new / rz_old;
                // |p> = |z> + beta |p>
                scale(p, nall, real_t(beta));
                add_product(p, nall, z, real_t(1));
            } // rz_old < tiny

            if (echo > 9) std::printf("# %s it=%i alfa=%g beta=%g\n", strip_path(__FILE__), it, alpha, beta);
            auto const inner = scalar_product(x, b, nall, comm) * pg.dV(); // this synchronization point is for display only
            if (echo > 7) std::printf("# %s it=%i res=%.2e E=%.15f\n", strip_path(__FILE__), it, std::sqrt(res2/cell_volume), inner);

            // rz_old = rz_new
            rz_old = rz_new;

            // decide if we continue to iterate
            run = (res2 > threshold2); // residual fell below threshold ?
            run = run || (it < miniter); // minimum number of steps not reached ?
            run = run && (it < maxiter); // maximum number of steps exceeded ?
//       !--------------------------------------
//       ! end of the CG iteration
//       !--------------------------------------
        } // while(run)

        double const res = std::sqrt(res2/cell_volume); // the residual has the unit of a density
        if (residual) *residual = res; // export

        // show the result
        if (echo > 2) std::printf("# %s %.2e -> %.2e e/Bohr^3%s in %d%s iterations\n", strip_path(__FILE__),
            res_start, res, (res < threshold)?" converged":"", it, (it < maxiter)?"":" (maximum)");

        auto const inner = scalar_product(x, b, nall, comm) * pg.dV();
        if (echo > 5) std::printf("# %s inner product <x|b> = %.15f\n", strip_path(__FILE__), inner);

        if (nullptr != inner_xx_bb) *inner_xx_bb = inner;

        set(xx, nall, x);

        return (res > threshold); // returns 0 when converged
    } // solve



















#ifdef    NO_UNIT_TESTS
    template // explicit template instantiation for double (and float implicitly)
    status_t solve(double*, double const*, parallel_grid_t const &, char, int, float, float*, int, int, int, double*);

    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    inline size_t parallel_grid_index(uint32_t const nb[3], int const ix, int const iy, int const iz) {
        uint32_t const i888 = ((iz & 0x7)*8 + (iy & 0x7))*8 + (ix & 0x7);
        uint64_t const i512 = ((iz >> 3)*nb[1] + (iy >> 3))*nb[0] + (ix >> 3);
        return (i512 << 9) + i888;
    } // parallel_grid_index

    template <typename real_t>
    status_t test_solver(int const echo=9, uint32_t const nb_default=4) {
        uint32_t const nb[] = {nb_default, nb_default, nb_default}; // number of 8*8*8 blocks
        real_space::grid_t g(nb[0]*8, nb[1]*8, nb[2]*8); // grid spacing == 1.0
        if (echo > 2) std::printf("\n# %s<%s> ng=[%d %d %d]\n", __func__, (8 == sizeof(real_t))?"double":"float", g[0], g[1], g[2]);
        g.set_boundary_conditions(1); // all boundary conditions periodic, ToDo: fails for isolated BCs
        auto const ng_all = size_t(g[2])*size_t(g[1])*size_t(g[0]);
        view2D<real_t> xb(4, ng_all, real_t(0)); // get memory
        auto const x = xb[0], x_fft = xb[1], b = xb[2], b_fft = xb[3];
        double constexpr c1 = 1, a1=.125, c2 = -8 + 1.284139e-7, a2=.5; // parameters for two Gaussians, in total close to neutral
        double const cnt[] = {.5*g[0], .5*g[1], .5*g[2]};
        { // scope: prepare the charge density (right-hand-side) rho
            double integral{0};
            for (int iz{0}; iz < g[2]; ++iz) {
            for (int iy{0}; iy < g[1]; ++iy) {
            for (int ix{0}; ix < g[0]; ++ix) {
                double const r2 = pow2(ix - cnt[0]) + pow2(iy - cnt[1]) + pow2(iz - cnt[2]);
                double const rho = c1*std::exp(-a1*r2) + c2*std::exp(-a2*r2);
                integral += rho;
                size_t const izyx = parallel_grid_index(nb, ix, iy, iz);
                b[izyx] = rho;
                size_t const ifft = (iz*g[1] + iy)*g[0] + ix;
                b_fft[ifft] = rho;
            }}} // ix iy iz
            if (echo > 3) std::printf("# %s integrated density %g\n", strip_path(__FILE__), integral*g.dV());
        } // scope

        load_balancing_t const lb(g, MPI_COMM_WORLD, 8, echo);
        parallel_grid_t const pg(g, lb, echo);

        view3D<real_t> xb_local(2, std::max(pg.n_local(), 1u), 512, real_t(0)); // create parallelized memory load
        { // scope: copy in
            auto const local_ids = pg.local_ids();
            for (int ilb{0}; ilb < pg.n_local(); ++ilb) {
                uint32_t ixyz[3]; global_coordinates::get(ixyz, local_ids[ilb]);
                size_t const j512 = (ixyz[2]*nb[1] + ixyz[1])*nb[0] + ixyz[0];
                set(xb_local(1,ilb), 512, b + j512*512); // copy one block of b
            } // ilb
        } // scope

        float const threshold = (sizeof(real_t) > 4) ? 3e-8 : 5e-6;
        auto const method = control::get("parallel_poisson.test.method", "mix");
        int  const max_it = control::get("parallel_poisson.test.maxiter", 999.);
        float residual_reached{0};

        auto const stat = solve(xb_local(0,0), xb_local(1,0), pg, *method, echo, threshold, &residual_reached, max_it);

        { // scope: copy out
            auto const local_ids = pg.local_ids();
            for (int ilb{0}; ilb < pg.n_local(); ++ilb) {
                uint32_t ixyz[3]; global_coordinates::get(ixyz, local_ids[ilb]);
                size_t const j512 = (ixyz[2]*nb[1] + ixyz[1])*nb[0] + ixyz[0];
                set(x + j512*512, 512, xb_local(0,ilb)); // copy one block of x
            } // ilb
            if (mpi_parallel::size() > 1) mpi_parallel::sum(x, ng_all);
        } // scope

        auto constexpr pi = constants::pi;
        { // scope: create a reference solution by FFT (not MPI parallel)
            int const ng[] = {g[0], g[1], g[2]};
            double const mat[3][4] = {{2*pi/g[0],0,0, 0},{0,2*pi/g[1],0, 0}, {0,0,2*pi/g[2], 0}};
            auto const stat_fft = fourier_poisson::solve(x_fft, b_fft, ng, mat);
            if (0 != stat_fft) warn("fourier_poisson::solve returned status= %i", int(stat_fft));
        } // scope

        if (0 == stat && echo > 7) { // get a radial representation from a point cloud plot
            float const compressed = control::get("parallel_poisson.plot.compressed", 1e-9); // 0: do not even sort, <0: plot all points, >0: use RDP compression
            int  const sorted = (0 != compressed);
            auto const ng_all = size_t(g[2])*size_t(g[1])*size_t(g[0]);
            std::vector<std::array<float,4>> vec(sorted*ng_all);
            if (0 == compressed) std::printf("\n\n## r, V_fd, V_fft, rho (all in a.u.)\n"); // show all grid values
            for (int iz{0}; iz < g[2]; ++iz) {
            for (int iy{0}; iy < g[1]; ++iy) {
            for (int ix{0}; ix < g[0]; ++ix) {
                    size_t const ifft = (iz*g[1] + iy)*g[0] + ix;
                    size_t const izyx = parallel_grid_index(nb, ix, iy, iz);
                    double const r2 = pow2(ix - cnt[0]) + pow2(iy - cnt[1]) + pow2(iz - cnt[2]), r = std::sqrt(r2);
                    if (0 == sorted) {
                        std::printf("%g %g %g %g\n", r, x[izyx], x_fft[ifft], b[izyx]); // point cloud (dots should not be connected by a line)
                    } else {
                        vec[izyx] = {float(r), float(x[izyx]), float(x_fft[ifft]), float(b[izyx])}; // store
                    }
            }}} // ix iy iz
            if (sorted) {
                auto lambda = [](std::array<float,4> const & left, std::array<float,4> const & right) { return left[0] < right[0]; };
                std::stable_sort(vec.begin(), vec.end(), lambda);
                view2D<float> columns(4*(compressed > 0), ng_all, 0.f);
                for (size_t izyx{0}; izyx < ng_all; ++izyx) {
                    auto const & v = vec[izyx];
                    if (0 == compressed) {
                        std::printf("%g %g %g %g\n", v[0], v[1], v[2], v[3]); // V_fd, V_fft and rho as lines
                    } else {
                        for (int ic{0}; ic < 4; ++ic) { columns(ic,izyx) = v[ic]; } // store
                    }
                } // izyx
                for (int ic{1}; ic < 4*(compressed > 0); ++ic) {
                    if (echo > 6 + ic) {
                        auto const *const label = (1==ic)?"V_fd":((2==ic)?"V_fft":"rho");
                        std::printf("\n# r, %s (in a.u.):\n", label);
                        auto const n_printed = print_compressed(columns[0], columns[ic], ng_all, compressed);
                        std::printf("# %s lossfully compressed from %ld to %ld points at threshold= %.1e\n\n", label, ng_all, n_printed, compressed);
                    } // echo
                } // ic
            } // sorted
            std::printf("\n"); // empty line to separate plots from other data
            if (echo > 7) {
                std::printf("\n\n# r, V, rho\n"); // also plot the radial function of V_analytical and rho
                auto const sa1 = std::sqrt(a1), sa2 = std::sqrt(a2), sc1 = c1/(4*a1*sa1), sc2 = c2/(4*a2*sa2);
                for (int ir{0}; ir <= 80*nb_default; ++ir) {
                    auto const r = 0.1*ir, r2 = r*r;
                    auto const rho = c1*std::exp(-a1*r2) + c2*std::exp(-a2*r2);
                    auto const V = (r < 1e-6) ? (sc1*2*sa1 + sc2*2*sa2)*4*pi :
                                   (sc1*std::erf(sa1*r) + sc2*std::erf(sa2*r))*(constants::sqrtpi*4*pi)/r;
                    std::printf("%g %g %g\n", r, V, rho);
                } // ir
            } // echo
            std::fflush(stdout);
        } // echo
        mpi_parallel::barrier();
        return stat;
    } // test_solver

    status_t test_parallel_grid(int const echo=0) {
        // test all combinations of isolated and periodic boundary conditions
        uint32_t const gm = control::get("parallel_poisson.grid.max", 9.); // and grids up to this number^3
        int8_t constexpr nBCs = 2; // can be used to limit it to one
        int8_t const BCs[] = {Isolated_Boundary, Periodic_Boundary};
            if (echo > 7) std::printf("\n#\n");
            // test various combinations of grids
        char what[] = "???";
        for (char w{'F'}; w <= 'I'; w += 'I' - 'F') { what[0] = w;
        for (uint32_t gz{1}; gz <= 1+0*gm; ++gz) {
        for (uint32_t gy{3}; gy <= 3+0*gm; ++gy) {
        for (uint32_t gx{1}; gx <= 1+0*gm; ++gx) {
            if (echo > 9) std::printf("\n\n\n\n\n\n\n\n\n\n\n\n\n");
            real_space::grid_t g(8*gx, 8*gy, 8*gz);
            if (echo > 7) std::printf("\n#\n# %s with grid [%d %d %d]\n", __func__, g[0], g[1], g[2]);
            load_balancing_t const lb(g, MPI_COMM_WORLD, 8, echo); // reacts to +parallel_poisson.nprocs (fake MPI processes)
        for (int8_t bz{0}; bz < nBCs; ++bz) {
        for (int8_t by{0}; by < nBCs; ++by) {
        for (int8_t bx{0}; bx < nBCs; ++bx) {
            int8_t const bc[] = {BCs[bx], BCs[by], BCs[bz]};
            if (echo > 3) std::printf("# %s with boundary conditions [%d %d %d]\n", __func__, bc[0], bc[1], bc[2]);
            std::fflush(stdout);
            mpi_parallel::barrier(); 

            g.set_boundary_conditions(bc);
            parallel_grid_t pg(g, lb, echo/8, what);

            }}} // bx by bz
        }}}} // gx gy gz w
        return 0;
    } // test_parallel_grid

    template <typename real_t>
    status_t test_Laplace16th(int8_t const bc[3], int const echo=9) {
        if (echo > 4) std::printf("\n# %s<%s>(bc=[%d %d %d])\n", __func__, (8 == sizeof(real_t))?"double":"float", bc[0], bc[1], bc[2]);
        status_t stat(0);
        real_space::grid_t g(4*8, 4*8, 4*8);
        load_balancing_t const lb(g, MPI_COMM_WORLD, 8, echo);
        g.set_boundary_conditions(bc);
        parallel_grid_t pg(g, lb, echo);
        auto const nl = pg.n_local(), nr = pg.n_remote();
        view3D<real_t> xAx(2, std::max(1, int(nl + nr)), 512, real_t(0));
        auto const  x = (real_t (*)[512]) xAx(0,0);
        auto const Ax = (real_t (*)[512]) xAx(1,0);
        for (int il{0}; il < nl; ++il) { 
            for (int i512{0}; i512 < 512; ++i512) {
                x[il][i512] = 1;
            } // i512
        } // il
        if (echo > 3) std::printf("# %s array prepared: x[%d][512]\n", __func__, std::max(1, int(nl + nr)));

        stat += Laplace16th(Ax[0], x[0], pg, echo);

        if (pg.n_local() > 0) {
            double diff{0};
            for (int i512{0}; i512 < 512; ++i512) { diff += std::abs(Ax[0][i512]); }
            if (echo > 3) std::printf("# %s<%s>: diff= %g\n\n", __func__, (8 == sizeof(real_t))?"double":"float", diff/512);
        } // n_local > 0
        // with double_t=float, we find diff=5.5e-7 per grid point,
        // with double_t=double         diff=6.3e-16 
        return stat;
    } // test_Laplace16th

    status_t test_Laplace16th_bc(int const echo=0) { // test a single iteration
        status_t stat(0);
        int8_t const BCs[] = {Isolated_Boundary, Periodic_Boundary};
        for (int8_t bz{0}; bz < 2; ++bz) {
        for (int8_t by{0}; by < 2; ++by) {
        for (int8_t bx{0}; bx < 2; ++bx) {
            int8_t const bc[] = {BCs[bx], BCs[by], BCs[bz]};
            stat += test_Laplace16th<double>(bc, echo);
            stat += test_Laplace16th<float> (bc, echo);
        }}} // bx by bz
        return stat;
    } // test_Laplace16th_bc

    status_t all_tests(int const echo) {
        auto const already_initialized = mpi_parallel::init();
        status_t stat(0);
        int n{0}; auto const t = int(control::get("parallel_poisson.select.test", -1.)); // -1:all
        if (t & (1 << n++)) stat += std::abs(test_parallel_grid(echo));
        if (t & (1 << n++)) stat += std::abs(test_solver<double>(echo)); // instantiation for both, double and float
        if (t & (1 << n++)) stat += std::abs(test_solver<float> (echo)); // compilation and convergence tests
        if (t & (1 << n++)) stat += std::abs(test_Laplace16th_bc(echo));
        if (!already_initialized) mpi_parallel::finalize();
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace parallel_poisson
