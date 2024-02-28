// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, ::snprintf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <algorithm> // std::swap<T>
#include <cmath> // std::sqrt
#include <type_traits> // std::is_same

#include "parallel_poisson.hxx"

#include "real_space.hxx" // ::grid_t
#include "data_view.hxx" // view2D<T>
#include "inline_math.hxx" // set, dot_product, pow2
#include "finite_difference.hxx" // ::stencil_t, ::apply
#include "constants.hxx" // ::pi
#include "boundary_condition.hxx" // Periodic_Boundary
#include "mpi_parallel.hxx" // MPI_Comm, MPI_COMM_WORLD, MPI_COMM_NULL
#include "load_balancer.hxx" // ::get, ::no_owner
#include "global_coordinates.hxx" // ::get
#include "boundary_condition.hxx" // Isolated_Boundary, Periodic_Boundary
#include "green_parallel.hxx" // ::exchange, ::RequestList_t
#include "print_tools.hxx" // printf_vector
#include "control.hxx" // ::get

#ifndef   NO_UNIT_TESTS
  #include "fourier_poisson.hxx" // ::solve
#endif // NO_UNIT_TESTS

namespace parallel_poisson {
  // solve the Poisson equation iteratively using the conjugate gradients method
  
  double constexpr m1over4pi = -.25/constants::pi; // -4*constants::pi is the electrostatics prefactor in Hartree atomic units

  template <typename real_t>
  double norm2(real_t const v[], size_t const n, MPI_Comm const comm=MPI_COMM_WORLD) {
      double s{0};
      for (size_t i{0}; i < n; ++i) { 
          s += pow2(v[i]);
      } // i 
      if (MPI_COMM_NULL != comm) mpi_parallel::sum(&s, 1, comm);
      return s;
  } // norm2

  template <typename real_t>
  double norm1(real_t const v[], size_t const n, MPI_Comm const comm=MPI_COMM_WORLD) {
      double s{0};
      for (size_t i{0}; i < n; ++i) {
          s += v[i];
      } // i
      if (MPI_COMM_NULL != comm) mpi_parallel::sum(&s, 1, comm);
      return s;
  } // norm1

  template <typename real_t>
  double scalar_product(real_t const v[], real_t const w[], size_t const n, MPI_Comm const comm=MPI_COMM_WORLD) {
      double dot{0};
      for (size_t i = 0; i < n; ++i) {
          dot += double(v[i])*double(w[i]); // conversion to double is different from dot_product define in inline_math.hxx
      } // i
      if (MPI_COMM_NULL != comm) mpi_parallel::sum(&dot, 1, comm);
      return dot;
  } // scalar_product

    class grid8x8x8_t {
    public:
        grid8x8x8_t() : n_local_blocks(0) { nb[0] = 0; nb[1] = 0; nb[2] = 0; } // default constructor
        grid8x8x8_t(uint32_t const ng[3], int8_t const bc[3], MPI_Comm const comm, int const echo=0) { // constructor
            if (echo > 3) std::printf("# %s(ng=[%d %d %d], bc=[%d %d %d])\n", __func__, ng[0], ng[1], ng[2], bc[0], bc[1], bc[2]);

            for (int d{0}; d < 3; ++d) {
                nb[d] = ng[d] >> 3; // divide by 8
                assert(nb[d]*8 == ng[d] && "grid numbers must be a positive multiple of 8");
                assert(nb[d] > 0 && "at least one block of 8 grid points needed");
            } // d

            // maybe the following 15 lines can be moved out of this constructor in the future, since load balancing should be done at program start
            view3D<uint16_t> owner_rank(nb[2],nb[1],nb[0], load_balancer::no_owner);
            int32_t const me = mpi_parallel::rank();
            double rank_center[4];
            { // scope: load balancing
                uint32_t const np = control::get("parallel_poisson.nprocs", double(mpi_parallel::size()));

                auto const load = load_balancer::get(np, me, nb, echo, rank_center, owner_rank.data());
                n_local_blocks = rank_center[3]; // the 4th component contains the number of items
                if (echo > 5) std::printf("# %s: load_balancer::get = %g, %d local blocks in rank#%i\n",
                                             __func__, load, n_local_blocks, me);

                // load_balancer has the mission to produce coherent domains with not too lengthy extents in space
                size_t const nall = size_t(nb[2])*size_t(nb[1])*size_t(nb[1]);
                if (MPI_COMM_NULL != comm) {
                    mpi_parallel::min(owner_rank.data(), nall, comm);
                    if (echo > 99) {
                        std::printf("# rank#%i owner_rank after  MPI_MIN ", me);
                        printf_vector(" %i", owner_rank.data(), nall);
                    } // echo
                } // not comm_null

            } // scope

            std::vector<int64_t> local_global_ids;
            { // scope: setup of star and remote_global_ids, determination of n_remote_blocks

                int32_t const iMAX = 2147483647;
                int32_t min_domain[] = {iMAX, iMAX, iMAX}, max_domain[] = {-1, -1, -1};
                double cnt_domain[] = {0, 0, 0};
                uint32_t nown{0};
                for (uint32_t iz = 0; iz < nb[2]; ++iz) {
                for (uint32_t iy = 0; iy < nb[1]; ++iy) {  // this triple loop does not scale well as it is the same work for all processes
                for (uint32_t ix = 0; ix < nb[0]; ++ix) {
                    if (me == owner_rank(iz,iy,ix)) {
                        ++nown;
                        int32_t const ixyz[] = {int32_t(ix), int32_t(iy), int32_t(iz)};
                        for (int d = 0; d < 3; ++d) {
                            min_domain[d] = std::min(min_domain[d], ixyz[d]);
                            max_domain[d] = std::max(max_domain[d], ixyz[d]);
                            cnt_domain[d] += ixyz[d];
                        } // d 
                    } // my
                }}} // iz iy ix
                if (nown != n_local_blocks) error("expected match between n_local_blocks=%d and count(owner_rank[]==me)=%d\n", n_local_blocks, nown);

                for (int d = 0; d < 3; ++d) {
                    cnt_domain[d] /= nown;
                    assert(max_domain[d] >= min_domain[d]);
                } // d
                if (echo > 7) std::printf("# rank#%i domain center %g %g %g rank center %g %g %g\n",
                        me, cnt_domain[0], cnt_domain[1], cnt_domain[2], rank_center[0], rank_center[1], rank_center[2]);

                int32_t constexpr HALO=1;
                int32_t const ndom[] = {max_domain[0] - min_domain[0] + 1 + 2*HALO,
                                        max_domain[1] - min_domain[1] + 1 + 2*HALO,
                                        max_domain[2] - min_domain[2] + 1 + 2*HALO};
                int32_t const ioff[] = {min_domain[0] - HALO, min_domain[1] - HALO, min_domain[2] - HALO};
                if (echo > 7) std::printf("# rank#%i domain(%d %d %d), halo= %d\n", me, ndom[0], ndom[1], ndom[2], HALO);

                uint8_t constexpr OUTSIDE=0, BORDER=1, INSIDE=2, VACUUM=7;
                view3D<uint8_t> domain(ndom[2],ndom[1],ndom[0], OUTSIDE);
                view3D<int32_t> domain_index(ndom[2],ndom[1],ndom[0], -1);

                local_global_ids.resize(n_local_blocks, -1);

                uint32_t ilb{0}; // index of the local block
                for (int32_t iz = HALO; iz < ndom[2] - HALO; ++iz) {
                for (int32_t iy = HALO; iy < ndom[1] - HALO; ++iy) {
                for (int32_t ix = HALO; ix < ndom[0] - HALO; ++ix) {
                    // translate domain indices {ix,iy,iz} into box indices
                    int32_t const ixyz[] = {ix + ioff[0], iy + ioff[1], iz + ioff[2]};
                    for (int d = 0; d < 3; ++d) { assert(ixyz[d] >= 0); assert(ixyz[d] < nb[d]); }
                    if (me == owner_rank(ixyz[2],ixyz[1],ixyz[0])) {
                        domain(iz,iy,ix) |= INSIDE;
                        domain_index(iz,iy,ix) = ilb; // domain index for inner elements

                        // mark lowest-order finite-difference sourrounding as border
                        domain(iz,iy,ix - 1) |= BORDER;
                        domain(iz,iy,ix + 1) |= BORDER;
                        domain(iz,iy - 1,ix) |= BORDER;
                        domain(iz,iy + 1,ix) |= BORDER;
                        domain(iz - 1,iy,ix) |= BORDER;
                        domain(iz + 1,iy,ix) |= BORDER;

                        local_global_ids[ilb] = global_coordinates::get(ixyz);

                        ++ilb; // count inside elements
                    } // me == owner
                }}} // iz iy ix
                assert(ilb <= (1ull << 31));
                if (ilb != n_local_blocks) error("expected match between n_local_blocks=%d and count(owner_rank[]==me)=%d\n", n_local_blocks, ilb);


                // loop over the domain again, this time including the halos
                uint32_t jrb{0}; // index for border-only blocks
                size_t st[4] = {0, 0, 0}; // statistics for display
                for (int32_t iz = 0; iz < ndom[2]; ++iz) {
                for (int32_t iy = 0; iy < ndom[1]; ++iy) {
                for (int32_t ix = 0; ix < ndom[0]; ++ix) {
                    auto const dom = domain(iz,iy,ix);
                    if (BORDER == dom) { // is border-only
                        ++jrb; // count remote elements
                        domain_index(iz,iy,ix) = ilb; // domain index for border elements
                        ++ilb; // continue counter for all elements
                    } // border
                    ++st[dom & 0x3]; // dom should be in [0, 3] anyway but better safe than sorry
                }}} // iz iy ix
                if (echo > 5) std::printf("# rank#%i has %ld outside, %ld border, %ld inside, %ld inside+border elements\n",
                                                  me, st[OUTSIDE], st[BORDER], st[INSIDE], st[INSIDE+BORDER]);
                uint32_t const n_remote_blocks = jrb;
                if (echo > 5) std::printf("# rank#%i has %d local and %d remote blocks, %d in total\n",
                                                  me, n_local_blocks, n_remote_blocks, n_local_blocks + n_remote_blocks);
                
                remote_global_ids.resize(n_remote_blocks, -1); // init remote element request lists
                star = view2D<int32_t>(n_local_blocks, 6, -1); // init finite-difference neighborhood lists

                size_t vacuum_assigned{0};
                // loop over domain again (3rd time), with halos
                for (int32_t iz = HALO; iz < ndom[2] - HALO; ++iz) {
                for (int32_t iy = HALO; iy < ndom[1] - HALO; ++iy) {
                for (int32_t ix = HALO; ix < ndom[0] - HALO; ++ix) {
                    // translate domain indices {ix,iy,iz} into box indices
                    int32_t const ixyz[] = {ix + ioff[0], iy + ioff[1], iz + ioff[2]};
                    for (int d = 0; d < 3; ++d) { assert(ixyz[d] >= 0); assert(ixyz[d] < nb[d]); } // should still hold...
                    if (domain(iz,iy,ix) & INSIDE) {
                        auto const id0 = domain_index(iz,iy,ix);
                        assert(id0 >= 0); assert(id0 < n_local_blocks);
                        for (int d = 0; d < 3; ++d) {
                            for (int i01 = 0; i01 < 2; ++i01) {
                                int32_t jxyz[] = {ixyz[0], ixyz[1], ixyz[2]}; // global coordinates
                                int32_t jdom[] = {ix, iy, iz}; // domain coordinates
                                jxyz[d] += (i01*2 - 1); // add or subtract 1
                                jdom[d] += (i01*2 - 1); // add or subtract 1
                                // if (echo > 7) std::printf("# rank#%i %d %d %d %c%c1\n", me, ixyz[0], ixyz[1], ixyz[2], 'x'+d, i01?'+':'-');
                                assert(domain(jdom[2],jdom[1],jdom[0]) & BORDER); // was set in 1st domain loop
                                auto const bc_shift = int(jxyz[d] < 0) - int(jxyz[d] >= int32_t(nb[d]));
                                if ((0 != bc_shift) && (Isolated_Boundary == bc[d])) {
                                    // id = -1; // not existing
                                    assert(BORDER == domain(jdom[2],jdom[1],jdom[0])); // must be border-only
                                    domain(jdom[2],jdom[1],jdom[0]) = VACUUM; // correct the mask
                                    ++vacuum_assigned;
                                } else {
                                    jxyz[d] += bc_shift*int32_t(nb[d]); // fold back into [0, nb)
                                    for (int d = 0; d < 3; ++d) { assert(jxyz[d] >= 0); assert(jxyz[d] < nb[d]); }
                                    auto const gid = global_coordinates::get(jxyz);
                                    auto const klb = domain_index(jdom[2],jdom[1],jdom[0]);
                                    star(id0,2*d + i01) = klb;
                                    assert(-1 != klb);
                                    if (klb >= n_local_blocks) {
                                        assert(BORDER == domain(jdom[2],jdom[1],jdom[0])); // must be a border-only element
                                        auto const irb = klb - n_local_blocks;
                                        assert(irb >= 0);
                                        assert((gid == remote_global_ids[irb]) || (-1 == remote_global_ids[irb])); // either unassigned(-1) or the same value
                                        remote_global_ids[irb] = gid;
                                    } // add to remote list
                                } // block exists
                            } // i01
                        } // d
                        if (echo > 19) std::printf("# rank#%i star[%i,:] = {%i %i %i %i %i %i}\n", me, id0,
                                        star(id0,0), star(id0,1), star(id0,2), star(id0,3), star(id0,4), star(id0,5));
                    } // is inside
                }}} // iz iy ix

                size_t vacuum_requested{0};
                for (auto id : remote_global_ids) {
                    vacuum_requested += (-1 == id);
                } // id
                if (vacuum_requested && echo > 3) std::printf("# rank#%i assigned %ld, request %ld vacuum cells\n", me, vacuum_assigned, vacuum_requested);
                assert(vacuum_assigned == vacuum_requested);

            } // scope

            if (echo > 8) {
                std::printf("# %s: requests={", __func__);
                for (auto rq : remote_global_ids) {
                    std::printf(" %lli", rq);
                } // rq
                std::printf(" }, %ld items\n", remote_global_ids.size());

                std::printf("# %s: offering={", __func__);
                for (auto of : local_global_ids) {
                    std::printf(" %lli", of);
                } // of
                std::printf(" }, %ld items\n", local_global_ids.size());
            } // echo

            requests = green_parallel::RequestList_t(remote_global_ids, local_global_ids, owner_rank.data(), nb, echo);

            if (echo > 8) {
                std::printf("# %s: RequestList.owner={", __func__);
                for (auto ow : requests.owner) {
                    std::printf(" %i", ow);
                } // ow
                std::printf(" }, %ld items\n", requests.owner.size());
            } // echo

        } // preferred constructor

        uint32_t n() const { return n_local_blocks; }
        uint32_t n_remote() const { return remote_global_ids.size(); }
        int32_t const* getStar() const { return (int32_t const*)star.data(); }
    public:
        green_parallel::RequestList_t requests;
    private:
        std::vector<int64_t> remote_global_ids; // may contain "-1"-entries
        view2D<int32_t> star; // local indices of 6 nearest finite-difference neighbors, star(n_local_blocks,6)
        uint32_t nb[3]; // box of blocks
        uint32_t n_local_blocks; // number of blocks owned by this MPI rank
    }; // grid8x8x8_t

    template <typename real_t>
    status_t data_exchange(
          real_t *v  // input and result array, data layout v[n_local_remote][8*8*8]
        , grid8x8x8_t const & g8 // descriptor
        , MPI_Comm const comm=MPI_COMM_WORLD // ToDo
        , int const echo=0 // log-level
    ) {
        auto const count = 8*8*8;
        auto const n_local = g8.n();
        auto const stat = green_parallel::exchange(v + count*n_local, v, g8.requests, count, echo);
        return stat;
    } // data_exchange


    template <typename real_t, typename double_t=double>
    status_t Laplace16th(
          real_t       *Av // result array, data layout Av[n_local_blocks][8*8*8]
        , real_t       *v  // input  array, data layout  v[n_local_remote][8*8*8], cannot be real_t const to call data_exchange
        , grid8x8x8_t const & g8 // descriptor
        , double const h2[3] // 1./grid_spacings^2
        , MPI_Comm const comm=MPI_COMM_WORLD // MPI communicator
        , int const echo=0 // log level
        , double const prefactor=1
    ) {
        if (echo > 9) std::printf("\n# %s start\n", __func__);

        auto const stat = data_exchange(v, g8, comm, echo);

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
        int const nlb = g8.n();
        auto const star = (int32_t const(*)[6])g8.getStar();
        for (uint32_t ilb = 0; ilb < nlb; ++ilb) { // loop over local blocks --> CUDA block-parallel
            size_t const i512 = ilb << 9; // block offset
            auto const nn = star[ilb]; // nearest-neighbor blocks of block ilb
            for (int iz = 0; iz < 8; ++iz) {
            for (int iy = 0; iy < 8; ++iy) { // loops over block elements --> CUDA thread-parallel
            for (int ix = 0; ix < 8; ++ix) {
                auto const ixyz = iz*64 + iy*8 + ix;
                auto const i0 = i512 + ixyz;
                double_t const av = cFD[0]*double_t(v[i0]);
                double_t ax{av}, ay{av}, az{av}; // accumulators
                // if (echo > 9) std::printf("# Av[%i][%3.3o] init as %g\n", ilb, ixyz, av);
                for (int ifd = 1; ifd <= 8; ++ifd) {
                    // as long as ifd is small enough, we take from the central block of x, otherwise from neighbors
                    auto const ixm = (ix >=    ifd) ? i0 - ifd : (nn[0] << 9) + ixyz + 8 - ifd;
                    auto const ixp = (ix + ifd < 8) ? i0 + ifd : (nn[1] << 9) + ixyz - 8 + ifd;
                    ax += cFD[ifd]*(double_t(v[ixm]) + double_t(v[ixp]));
                    auto const iym = (iy >=    ifd) ? i0 - 8*ifd : (nn[2] << 9) + ixyz + 64 - 8*ifd;
                    auto const iyp = (iy + ifd < 8) ? i0 + 8*ifd : (nn[3] << 9) + ixyz - 64 + 8*ifd;
                    ay += cFD[ifd]*(double_t(v[iym]) + double_t(v[iyp]));
                    auto const izm = (iz >=    ifd) ? i0 - 64*ifd : (nn[4] << 9) + ixyz + 512 - 64*ifd;
                    auto const izp = (iz + ifd < 8) ? i0 + 64*ifd : (nn[5] << 9) + ixyz - 512 + 64*ifd;
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
  status_t solve(real_t x[] // result to Laplace(x)/(-4*pi) == b
                , real_t const b[] // right hand side b
                , real_space::grid_t const & g // grid descriptor
                , char const method // use mixed precision as preconditioner
                , int const echo // =0 // log level
                , float const threshold // =3e-8 // convergence criterion
                , float *residual // =nullptr // residual that was reached
                , int const maxiter // =999 // maximum number of iterations 
                , int const miniter // =0   // minimum number of iterations
                , int restart // =4096 // number of iterations before restart, 1:steepest descent
                ) {

    // auto const comm = MPI_COMM_WORLD; // green_parallel::comm();
    // uint32_t const ng[] = {uint(g[0]), uint(g[1]), uint(g[2])};
    // grid8x8x8_t g8(ng, g.boundary_conditions(), comm, echo);

    size_t const nall = size_t(g[2])*size_t(g[1])*size_t(g[0]);
    // when parallel, we also need to define nloc, the number of grid elements in this local process
    // also, we want to group the grid into cubes of 8x8x8 grid points

    status_t ist(0);

    restart = ('s' == method) ? 1 : std::max(1, restart);

    if (std::is_same<real_t, double>::value) {
        view2D<float> xb(2, nall, 0.0); // get memory
        auto const x32 = xb[0], b32 = xb[1];
        set(b32, nall, b); // convert to float
        set(x32, nall, x); // convert to float
        if (echo > 5) std::printf("# %s solve in <float> precision first\n", __FILE__);
        ist += solve(x32, b32, g, method, echo, threshold, residual, maxiter, miniter, restart);
        if (echo > 5) std::printf("# %s switch back to <double> precision\n", __FILE__);
        set(x, nall, x32); // convert to double
    } // real_t == double

    // we use CG + order-16 FD

    int const nn_precond = 0; // 0:none, >0:range-1-stencil, <0:multi_grid (does not work properly yet)
    bool const use_precond = (0 != nn_precond);

    // find memory aligned nloc
    view2D<real_t> mem(4 + use_precond, nall, 0.0); // get memory
    auto const r=mem[0], p=mem[1], ax=mem[2], ap=mem[3], z=use_precond?mem[4]:r;    

    finite_difference::stencil_t<real_t> const Laplacian(g.h, 8, m1over4pi); // 8: use a 17-point stencil

    finite_difference::stencil_t<real_t> precond;
    if (nn_precond > 0) {
        precond = finite_difference::stencil_t<real_t>(g.h, nn_precond);
        auto const nn = precond.nearest_neighbors();
        double nrm{0};
        for (int d = 0; d < 3; ++d) {
            for (int i = 0; i < nn[d]; ++i) {
                precond.c2nd[d][i] = std::abs(precond.c2nd[d][i]);
                nrm += precond.c2nd[d][i] * (1 + (i > 0));
            } // i
        } // d
        nrm = 1./nrm;
        for (int d = 0; d < 3; ++d) {
            for (int i = 0; i < nn[d]; ++i) {
                precond.c2nd[d][i] *= nrm;
            } // i
        } // d
        if (echo > 6) std::printf("# %s use a diffusion preconditioner with %d %d %d neighbors\n", 
                                __FILE__, nn[0], nn[1], nn[2]);
    } // use_precond

    double const cell_volume = nall*g.dV();
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
    ist = finite_difference::apply(ax, x, g, Laplacian);
    if (ist) error("CG_solve: Laplacian failed with status %i", int(ist));

    // dump_to_file("cg_start", nall, x, nullptr, 1, 1, "x", echo);
    
    if (g.number_of_boundary_conditions(Periodic_Boundary) == 3) {
        double const bnorm = norm1(b, nall)/nall * g.dV(); // g.comm
        if (echo > 8) std::printf("# %s all boundary conditions are periodic but system is charged with %g electrons\n", __FILE__, bnorm);
    } // all_boundary_conditions_periodic

    // |r> = |b> - A|x> = |b> - |Ax>
    set(r, nall, b); add_product(r, nall, ax, real_t(-1));

    // res^2 = <r|r>
    double res2 = norm2(r, nall) * g.dV(); // g.comm
    double const res_start = std::sqrt(res2/cell_volume); // store staring residual
    if (echo > 8) std::printf("# %s start residual=%.1e\n", __FILE__, res_start);

    // |z> = |Pr> = P|r>
    if (use_precond) {
        ist = finite_difference::apply(z, r, g, precond);
        if (ist) error("CG_solve: Preconditioner failed with status %i", int(ist));
    } else assert(z == r);

    // rz_old = <r|z>
    double rz_old = scalar_product(r, z, nall) * g.dV(); // g.comm 

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
        ist = finite_difference::apply(ap, p, g, Laplacian);
        if (ist) error("CG_solve: Laplacian failed with status %i", int(ist));

        double const pAp = scalar_product(p, ap, nall) * g.dV(); // g.comm

        // alpha = rz_old / pAp
        double const alpha = (std::abs(pAp) < RZ_TINY) ? RZ_TINY : rz_old / pAp;

        // |x> = |x> + alpha |p>
        add_product(x, nall, p, real_t(alpha));

//       !============================================================
//       ! special treatment of completely periodic case
//       !============================================================
        if (g.number_of_boundary_conditions(Periodic_Boundary) == 3) {
            double const xnorm = norm1(x, nall)/nall; // g.comm
  //        xnorm = xnorm/real( g%ng_all(1)*g%ng_all(2)*g%ng_all(3) )
            // subtract the average potential
            for (size_t i{0}; i < nall; ++i) x[i] -= xnorm;
        } // 3 periodic BCs
//       !============================================================

        if (0 == it % restart) {
            // |Ax> = A|x>
            ist = finite_difference::apply(ax, x, g, Laplacian);
            if (ist) error("CG_solve: Laplacian failed with status %i", int(ist))
            // |r> = |b> - A|x> = |b> - |ax>
            set(r, nall, b); add_product(r, nall, ax, real_t(-1));
        } else {
            // |r> = |r> - alpha |ap>
            add_product(r, nall, ap, real_t(-alpha));
        } // restart?

        // res = <r|r>
        res2 = norm2(r, nall) * g.dV(); // g.comm

        // |z> = |Pr> = P|r>
        if (use_precond) {
            ist = finite_difference::apply(z, r, g, precond);
            if (ist) error("CG_solve: Preconditioner failed with status %i", int(ist));
        } else assert(z == r);

        // rz_new = <r|z>
        double const rz_new = scalar_product(r, z, nall) * g.dV(); // g.comm

        // beta = rz_new / rz_old
        double beta = rz_new / rz_old;
        if (rz_old < RS_TINY) {
            beta = 0;
            set(p, nall, z);
        } else {
            // |p> = |z> + beta |p>
            scale(p, nall, real_t(beta));
            add_product(p, nall, z, real_t(1));
        } // rz_old < tiny

        if (echo > 9) std::printf("# %s it=%i alfa=%g beta=%g\n", __FILE__, it, alpha, beta);
        if (echo > 7) std::printf("# %s it=%i res=%.2e E=%.15f\n", __FILE__, it, 
            std::sqrt(res2/cell_volume), scalar_product(x, b, nall)*g.dV());

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

    double const res = std::sqrt(res2/cell_volume);
    if (residual) *residual = res; // export

    // show the result
    if (echo > 2) std::printf("# %s %.2e -> %.2e e/Bohr^3%s in %d%s iterations\n", __FILE__,
        res_start, res, (res < threshold)?" converged":"", it, (it < maxiter)?"":" (maximum)");
    if (echo > 5) std::printf("# %s inner product <x|b> = %.15f\n", __FILE__, scalar_product(x, b, nall)*g.dV());

    return (res > threshold);
  } // solve

#ifdef    NO_UNIT_TESTS
  template // explicit template instantiation for double
  status_t solve(double*, double const*, real_space::grid_t const &, char, int, float, float*, int, int, int);

  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

  template <typename real_t>
  status_t test_solver(int const echo=9, int const ng_default=32) {
      int const ng[] = {ng_default, ng_default, ng_default};
      if (echo > 2) std::printf("\n# %s<%s> ng=[%d %d %d]\n", __func__, (8 == sizeof(real_t))?"double":"float", ng[0], ng[1], ng[2]);
      real_space::grid_t g(ng[0], ng[1], ng[2]); // grid spacing = 1.0
      g.set_boundary_conditions(1); // all boundary conditions periodic, ToDo: fails for isolated BCs
      view2D<real_t> xb(3, ng[2]*ng[1]*ng[0], 0.0); // get memory
      auto const x = xb[0], x_fft = xb[1], b = xb[2];
      double constexpr c1 = 1, a1=.125, c2 = -8 + 1.284139e-7, a2=.5; // parameters for two Gaussians, in total close to neutral
      double const cnt[] = {.5*ng[0], .5*ng[1], .5*ng[2]};
      {
        double integral{0};
        for (int iz = 0; iz < ng[2]; ++iz) {
        for (int iy = 0; iy < ng[1]; ++iy) {
        for (int ix = 0; ix < ng[0]; ++ix) {
            size_t const ixyz = (iz*ng[1] + iy)*ng[0] + ix;
            double const r2 = pow2(ix - cnt[0]) + pow2(iy - cnt[1]) + pow2(iz - cnt[2]);
            double const rho = c1*std::exp(-a1*r2) + c2*std::exp(-a2*r2);
            b[ixyz] = rho;
            integral += rho;
        }}} // ix iy iz
        if (echo > 3) std::printf("# %s integrated density %g\n", __FILE__, integral*g.dV());
      }

      float const threshold = (sizeof(real_t) > 4) ? 3e-8 : 5e-6;
      auto const method = control::get("parallel_poisson.test.method", "mix");
      int  const max_it = control::get("parallel_poisson.test.maxiter", 999.);
      float residual_reached{0};
      auto const stat = solve(x, b, g, *method, echo, threshold, &residual_reached, max_it);

      { // scope: create a reference solution by FFT (not MPI parallel)
          auto constexpr pi = constants::pi;
          double const mat[3][4] = {{2*pi/ng[0],0,0, 0},{0,2*pi/ng[1],0, 0}, {0,0,2*pi/ng[2], 0}};
          fourier_poisson::solve(x_fft, b, ng, mat);
      } // scope

      if (echo > 7) { // get a radial representation through Bessel transform
          std::printf("\n## r, V_fd, V_fft, rho (all in a.u.)\n"); // show all grid values (dots should not be connected by a line)
          for (int iz = 0; iz < ng[2]; ++iz) {
           for (int iy = 0; iy < ng[1]; ++iy) {
            for (int ix = 0; ix < ng[0]; ++ix) {
                size_t const ixyz = (iz*ng[1] + iy)*ng[0] + ix;
                double const r2 = pow2(ix - cnt[0]) + pow2(iy - cnt[1]) + pow2(iz - cnt[2]);
                std::printf("%g %g %g %g\n", std::sqrt(r2), x[ixyz], x_fft[ixyz], b[ixyz]); // point cloud
          }}} // ix iy iz
          // also plot the radial function of the right-hand-side
          std::printf("\n# r, rho\n");
          for (int ir{0}; ir <= 10*ng_default; ++ir) {
              auto const r = 0.1*ir, r2 = r*r;
              std::printf("%g %g\n", r, c1*std::exp(-a1*r2) + c2*std::exp(-a2*r2));
          } // ir
      } // echo
      return stat;
  } // test_solver

  status_t test_grid8(int const echo=0) {
      // test all combinations of isolated and periodic boundary conditions
      uint32_t const gm = control::get("parallel_poisson.grid.max", 9.); // and grids up to this number^3
      int8_t const BCs[] = {Isolated_Boundary, Periodic_Boundary};
      for (int8_t bz = 0; bz < 2; ++bz) {
      for (int8_t by = 0; by < 2; ++by) {
      for (int8_t bx = 0; bx < 2; ++bx) {
        int8_t const bc[] = {BCs[bx], BCs[by], BCs[bz]};
        if (echo > 7) std::printf("\n#\n");
        if (echo > 3) std::printf("# %s with boundary conditions [%d %d %d]\n", __func__, bc[0], bc[1], bc[2]);
        // test various combinations of grids
        for (uint32_t gz = 1; gz <= gm; ++gz) {
        for (uint32_t gy = 1; gy <= gm; ++gy) {
        for (uint32_t gx = 1; gx <= gm; ++gx) {
          uint32_t const ng[] = {8*gx, 8*gy, 8*gz};
          if (echo > 7) std::printf("\n#\n# %s with grid [%d %d %d]\n", __func__, gx, gy, gz);

          grid8x8x8_t g8(ng, bc, MPI_COMM_WORLD, echo/8); // reacts to +parallel_poisson.nprocs (fake MPI processes)

        }}} // gx gy gz
      }}} // bx by bz
      return 0;
  } // test_grid8

    template <typename real_t>
    status_t test_Laplace16th(int8_t const bc[3], int const echo=9) {
        if (echo > 4) std::printf("\n# %s<%s>(bc=[%d %d %d])\n", __func__, (8 == sizeof(real_t))?"double":"float", bc[0], bc[1], bc[2]);
        status_t stat(0);
        uint32_t const ng[] = {2*8, 2*8, 2*8};
        double const h2[] = {1, 1, 1}; // unity grid spacing
        grid8x8x8_t g8(ng, bc, MPI_COMM_WORLD, echo);
        auto const n = g8.n() + g8.n_remote();
        view3D<real_t> xAx(2, n, 512, real_t(0));
        auto const  x = (real_t (*)[512]) xAx(0,0);
        auto const Ax = (real_t (*)[512]) xAx(1,0);
        for (int i512 = 0; i512 < n*512; ++i512) { x[0][i512] = 1; }

        stat = Laplace16th(Ax[0], x[0], g8, h2, MPI_COMM_WORLD, echo);

        if (g8.n() > 0) {
            double diff{0};
            for (int i512 = 0; i512 < 512; ++i512) { diff += std::abs(Ax[0][i512]); }
            if (echo > 3) std::printf("# %s<%s>: diff= %g\n\n", __func__, (8 == sizeof(real_t))?"double":"float", diff/512);
        }
        // with double_t=float, we find diff=5.5e-7 per grid point,
        // with double_t=double         diff=6.3e-16 
        return stat;
    } // test_Laplace16th

    status_t test_Laplace16th_boundary_conditions(int const echo=0) {
        status_t stat(0);
        int8_t const BCs[] = {Isolated_Boundary, Periodic_Boundary};
        for (int8_t bz = 0; bz < 2; ++bz) {
        for (int8_t by = 0; by < 2; ++by) {
        for (int8_t bx = 0; bx < 2; ++bx) {
            int8_t const bc[] = {BCs[bx], BCs[by], BCs[bz]};
            stat += test_Laplace16th<double>(bc, echo);
            stat += test_Laplace16th<float> (bc, echo);
        }}} // bx by bz
        return stat;
    } // test_Laplace16th_boundary_conditions

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_grid8(echo);
      stat += test_solver<double>(echo); // instantiation for both, double and float
 //   stat += test_solver<float>(echo);  // compilation and convergence tests
      stat += test_Laplace16th_boundary_conditions(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace parallel_poisson
