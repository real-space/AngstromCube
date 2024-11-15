// This file is part of AngstromCube under MIT License
/*
 * load balancing for A43
 *    main purpose is the distribution of right-hand-side blocks
 *    to MPI processes and GPUs
 *
 */


#include <cstdio> // std::printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <cstdint> // size_t, uint32_t
#include <algorithm> // std::max, ::min, ::swap
#include <utility> // std::pair
#include <cmath> // std::ceil, ::pow, ::cbrt
#include <numeric> // std::iota
#ifndef   NO_UNIT_TESTS
    #include <cstdlib> // rand, RAND_MAX
    #include "progress_report.hxx" // ProgressReport
#endif // NO_UNIT_TESTS undefined

#include "load_balancer.hxx"

#include "inline_math.hxx" // set, pow2
#include "status.hxx" // status_t
#include "simple_stats.hxx" // ::Stats<>
#include "recorded_warnings.hxx" // warn
#include "constants.hxx" // ::pi
#include "print_tools.hxx" // printf_vector

#ifndef   NO_UNIT_TESTS
    #include "control.hxx" // ::get
    #ifdef    HAS_BITMAP_EXPORT
        #include "bitmap.hxx" // ::write_bmp_file
    #endif // HAS_BITMAP_EXPORT
#endif // NO_UNIT_TESTS undefined

namespace load_balancer {

#define   LOAD_BALANCER_DRAW_SVG
#ifdef    LOAD_BALANCER_DRAW_SVG
    static std::vector<double> draw2D; // global field
    int constexpr echo_record_planes = 10; // only record planes when echo level is high
#endif // LOAD_BALANCER_DRAW_SVG

    int constexpr X=0, Y=1, Z=2, W=3;

    template <typename real_t>
    double center_of_weight(
          double cow[4] // result: center of weight [0/1/2] and number of contributors [3]
        , size_t const nuna // number of unassigned work items
        , uint32_t const indirect[] // list of unassigned work items
        , real_t const (*const xyzw)[4] // positions of work items in space [0/1/2] and their weight [3]
    ) {
        set(cow, 4, 0.0); // initialize
        double w8sum{0};
        for (size_t iuna{0}; iuna < nuna; ++iuna) { // parallel, reduction(+:w8sum,cow)
            auto const iall = indirect[iuna];
            auto const *const xyz = xyzw[iall];
            double const w8 = xyz[W];
            w8sum += w8;
            // contributes if the weight is positive
            double const w8pos = double(w8 > 0);
            add_product(cow, 3, xyz, w8pos);
            cow[W] += w8pos;
        } // iall
        if (cow[W] > 0) scale(cow, 3, 1./cow[3]); // normalize
        return w8sum; // returns the sum of weights
    } // center_of_weight

    // idea for a stable load balancer:
    //    for np processes, uses celing(log_2(nprocs)) iterations
    //    in iteration #0, place the center at 0,0,0 and the diagonal opposite corner finding its longest extent
    //                     divide the work by a plane sweep into chunks proportional to ceiling(np/2) and floor(np/2)
    //    pass the new processor numbers to the next iteration ...
    //    last iteration: number of processors is 1 --> done or 2 --> as before

    // find the largest extent:
    // 1) find the center of weight (not accounting the load but only if there is a non-zero load)
    // 2) find the maximum distance^2 and its position
    // 3) assume that the largest extent is along the line through the center of weight and that position
    // --> a cuboid will always use the space diagonal

    template <typename real_t, typename real_w_t=double>
    double plane_balancer(
            int const nprocs // number of MPI processes to distribute the work items evenly to
        , int const rank   // my MPI rank
        , size_t const nall // number of all work items
        , real_t const (*const xyzw)[4] // xyzw[nall][4], positions [0/1/2] and weights [3] of the work items
        , real_w_t const w8s[]        // w8s[nall] weights separated
        , double const w8sum_all=1. // denominator of all weights
        , int const echo=0 // verbosity
        , double rank_center[4]=nullptr // export the rank center [0/1/2] and number of items [3]
        , uint16_t *const owner_rank=nullptr // export the rank of each task, [nall]
    ) {
        // complexity is order(N^2) as each processes loops over all tasks in the first iteration

        auto constexpr epsilon = 1e-6; // accuracy for sanity check

        bool constexpr UNASSIGNED = 1, ASSIGNED = 0;
        std::vector<bool> state(nall, UNASSIGNED);

        assert(nall <= (1ull << 32) && "using uint32_t indirection lists limits the total number to 2^32");
        std::vector<uint32_t> indirect(nall);
        std::iota(indirect.begin(), indirect.end(), uint32_t(0)); // initialize with {0, 1, ..., nall-1}

        size_t nuna{nall}; // number of unassigned blocks

        double load_now{w8sum_all};

        int np{nprocs}; // current number of processes among which we distribute
        int rank_offset{0}; // offset w.r.t. MPI ranks

        int tree_level{0};
        while (np > 1) {

            assert(rank_offset + np <= nprocs);
            // bisect the workload for np processors into (np + 1)/2 and np/2
            int const nhalf[] = {(np + 1) >> 1, np >> 1}; // larger part and smaller part

            // MPIrank ranges are 0:{off, ..., off+nhalf[0]-1} and 1:{off+nhalf[0], ..., off+np-1}
            int const i01 = (rank >= rank_offset + nhalf[0]);
            // i01 == 0: this rank is part of the first set of (np + 1)/2 processes
            // i01 == 1: this rank is part of the second set of (np)/2    processes
            if (echo > 19) std::printf("# rank#%i divides %d into %d and %d\n", rank, np, nhalf[i01], nhalf[1 - i01]);
            assert(nhalf[0] + nhalf[1] == np);

            // determine the direction of the largest extent

            // determine the center of weight
            double cow[4];
            auto const w8sum = center_of_weight(cow, nuna, indirect.data(), xyzw);

            // determine the largest distance^2 from the center of weight
            double maxdist2{-1}; int64_t imax{-1}; // block with the largest distance
            for (size_t iuna{0}; iuna < nuna; ++iuna) { // serial due to special reduction
                auto const iall = indirect[iuna];
                auto const *const xyz = xyzw[iall];
                auto const dist2 = pow2(xyz[X] - cow[X]) + pow2(xyz[Y] - cow[Y]) + pow2(xyz[Z] - cow[Z]);
                if (dist2 > maxdist2) { maxdist2 = dist2; imax = iall; }
            } // iuna

            if (imax < 0) {
                // this should only happens if the process is idle,
                //   i.e. there are no blocks assigned to this one
                load_now = 0;
                np = 1; // stop the while-loop
            } else {

                // determine a sorting direction
                add_product(cow, 3, xyzw[imax], -1.);
                auto const len2 = pow2(cow[X]) + pow2(cow[Y]) + pow2(cow[Z]);
                auto const norm = (len2 > 0) ? 1./std::sqrt(len2) : 0.0;
                double const vec[] = {cow[X]*norm, cow[Y]*norm, cow[Z]*norm};
                if (echo > 19) std::printf("# rank#%i sort along the [%g %g %g] direction\n", rank, vec[X], vec[Y], vec[Z]);

                using fui_t = std::pair<float,uint32_t>;
                std::vector<fui_t> v(nuna);

                for (size_t iuna{0}; iuna < nuna; ++iuna) { // parallel
                    auto const iall = indirect[iuna];
                    auto const *const xyz = xyzw[iall];
                    auto const f = xyz[X]*vec[X] // inner product
                                 + xyz[Y]*vec[Y]
                                 + xyz[Z]*vec[Z];
                    v[iuna].first  = f;
                    v[iuna].second = iall;
                } // iall

                auto const sort_lambda = [](fui_t i1, fui_t i2) { return i1.first < i2.first; };
                std::stable_sort(v.begin(), v.end(), sort_lambda);

                auto const target_load0 = nhalf[0]*w8sum; // relative target for load0 multiplied with np
                { // scope: distribute according to target loads
                    auto const state0 = i01 ? ASSIGNED : UNASSIGNED,
                               state1 = i01 ? UNASSIGNED : ASSIGNED;
                    double load0{0}, load1{0}; // relative load multiplied with np to avoid floating point errors
                    size_t isrt{0};
                    for (; load0 < target_load0; ++isrt) { // serial
                        auto const iall = v[isrt].second;
                        load0 += w8s[iall]*np;
                        state[iall] = state0;
                    } // while load1 < target_load1
                    auto const isrt_middle = isrt;
                    for(; isrt < nuna; ++isrt) { // parallel reduction(+:load1)
                        auto const iall = v[isrt].second;
                        load1 += w8s[iall]*np;
                        state[iall] = state1;
                    } // isrt
                    assert(std::abs((load0 + load1) - (w8sum*np)) < epsilon*(w8sum*np) && "Maybe failed due to accuracy issues");
                    load_now = (i01 ? load1 : load0)/np;

                    if (echo > 29) std::printf("# plane level=%d %g %g %g isrt=%lu %d|%d\n", tree_level, vec[X], vec[Y], vec[Z], isrt_middle, nhalf[0],nhalf[1]);
#ifdef    LOAD_BALANCER_DRAW_SVG
                    if (echo > echo_record_planes) { // show bisecting plane
                        // bisecting plane normal is the sorting vector vec, plane distance from the origin is ?
                        double pd{0}; int den{0};
                        if (isrt_middle < nuna) { pd += v[isrt_middle    ].first; ++den; } // distance of the point that is closest to the plane and belongs to load1
                        if (isrt_middle > 0)    { pd += v[isrt_middle - 1].first; ++den; } // distance of the point that is closest to the plane and belongs to load0
                        if (rank == rank_offset) { // only the "lower" half stores the separating plane
                            std::printf("# plane level=%d %g %g %g  dist= %g  isrt=%lu %d|%d\n", tree_level, vec[X], vec[Y], vec[Z], pd/den, isrt_middle, nhalf[0],nhalf[1]);
                            // store the 2D plane in a global variable to be drawn into an SVG later
                            auto const s = draw2D.size();
                            if (s > 0) {
                                assert(den > 0); // if 2==den we take the average between v[isrt_middle].first and v[isrt_middle-1].first
                                draw2D.resize(s + 4);
                                draw2D[s + 0] = vec[X];
                                draw2D[s + 1] = vec[Y];
                                draw2D[s + 2] = pd/den; // distance to origin
                                draw2D[s + 3] = tree_level;
                            } // s > 0
                        } // rank == rank_offset
                    } // echo
#endif // LOAD_BALANCER_DRAW_SVG
                } // scope

                if (echo > 19) std::printf("# rank#%i assign %g of %g (%.2f %%, target %.2f %%) to %d processes\n",
                                        rank, load_now, w8sum, load_now*100./w8sum, nhalf[i01]*100./np, nhalf[0]);

                // prepare for the next iteration
                rank_offset += i01*nhalf[0];
                np = nhalf[i01];
                // update the indirection list, determine which blocks are still unassigned
                size_t iunk{0}; // counter
                for (size_t iuna{0}; iuna < nuna; ++iuna) { // serial loop due to read-write access to array indirect
                    auto const iall = indirect[iuna]; // read from array indirect
                    if (UNASSIGNED == state[iall]) {
                        indirect[iunk] = iall; // write to array indirect at a lower address
                        ++iunk;
                    } // is_unassigned
                } // iall
                assert(iunk <= nuna);
                nuna = iunk; // set new number of unassigned blocks

            } // imax < 0

            ++tree_level;
        } // while np > 1

        if (nullptr != rank_center) {
            if (load_now > 0) {
                // compute the center of weight again, for display and export
                auto const w8sum = center_of_weight(rank_center, nuna, indirect.data(), xyzw);
                if (echo > 13) std::printf("# rank#%i assign %.3f %% center %g %g %g, %g items\n",
                    rank, w8sum*100/w8sum_all, rank_center[X], rank_center[Y], rank_center[Z], rank_center[W]);
            } else {
                set(rank_center, 4, 0.);
            } // load_now > 0
        } // rank_center

        if (echo > 9) std::printf("# rank#%i load target %.3f %%, assign %.3f %%\n",
                                    rank, 100./nprocs, load_now*100/w8sum_all);

        if (nullptr != owner_rank) {
            for (size_t iall{0}; iall < nall; ++iall) {
                if (UNASSIGNED == state[iall]) {
                    assert(no_owner == owner_rank[iall]);
                    owner_rank[iall] = rank;
                    assert(owner_rank[iall] == rank && "uint16_t too short for owner_ranks");
                } // unassigned
            } // iall
            // Beware: only the owned entries of owner_rank have been modified, so
            //         an MPI_MIN-Allreduce needs to be performed later. This is skipped here
            //         to maintain serial executability of this routine.
            if (echo > 99) { std::printf("# rank#%i owner_rank before MPI_MIN ", rank); printf_vector(" %i", owner_rank, nall); }
        } // owner_rank

        if (1) { // parallelized consistency check
            for (size_t iuna{0}; iuna < nuna; ++iuna) { // parallel
                auto const iall = indirect[iuna];
                auto const w8 = xyzw[iall][W];
                assert(w8 == w8s[iall] && "weights inconsistent");
            } // iuna
        } // consistency check

        return load_now;
    } // plane_balancer


    double get(
            uint32_t const comm_size // number of MPI processes in this communicator
        , int32_t  const comm_rank // rank of this MPI process
        , uint32_t const nb[3] // number of blocks in X/Y/Z direction
        , int const echo // =0, log level
        , double rank_center[4] // =nullptr, export the rank center [0/1/2] and number of items [3]
        , uint16_t *const owner_rank // =nullptr, export the owner rank of each task, [nb[Z]*nb[Y]*nb[X]]
    ) {
        // distribute a rectangular box of nb[X] x nb[Y] x nb[Z] with all weights 1

        assert(comm_rank >= 0);
        assert(comm_rank < comm_size);

        auto const nall = nb[X]*size_t(nb[Y])*size_t(nb[Z]);
        assert(nall > 0);
        assert(nall <= (size_t(1) << 32) && "uint32_t is not long enough!");

        double w8sum_all{0};
        auto const xyzw = new float[nall][4];
        std::vector<double> w8s(nall, 0);

        for (uint32_t iz{0}; iz < nb[Z]; ++iz) {
        for (uint32_t iy{0}; iy < nb[Y]; ++iy) {
        for (uint32_t ix{0}; ix < nb[X]; ++ix) {
            auto const iall = (size_t(iz)*nb[Y] + iy)*size_t(nb[X]) + ix;
//          assert(uint32_t(iall) == iall && "uint32_t is not long enough!");
            float const w8 = 1.f; // weight(ix,iy,iz); // WEIGHTS CAN BE INSERTED HERE
            w8s[iall]     = w8;
            w8sum_all    += w8;
            xyzw[iall][X] = ix + .5f;
            xyzw[iall][Y] = iy + .5f;
            xyzw[iall][Z] = iz + .5f;
            xyzw[iall][W] = w8;
        }}} // ix iy iz

        auto const load_now = plane_balancer(comm_size, comm_rank, nall, xyzw, w8s.data(), w8sum_all, echo,
                                            rank_center, owner_rank);
        delete[] xyzw;
        return load_now;
    } // get














#ifdef    NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    template <typename real_t>
    real_t analyze_load_imbalance(real_t const load[], int const nprocs, int const echo=1) {
        simple_stats::Stats<> st(0);
        for (int rank{0}; rank < nprocs; ++rank) {
            st.add(load[rank]);
            if (echo > 19) std::printf("# myrank=%i owns %g\n", rank, load[rank]);
        } // rank
        auto const mean = st.mean(), max = st.max();
        if (echo > 0) std::printf("# %ld processes own %g blocks, per process [%g, %.2f +/- %.2f, %g]\n",
                                        st.tim(), st.sum(), st.min(), mean, st.dev(), max);
        return (mean > 0) ? max/mean : -1;
    } // analyze_load_imbalance

    template <typename real_t>
    real_t distance_squared(real_t const a[], real_t const b[]) {
        return pow2(a[X] - b[X]) + pow2(a[Y] - b[Y]) + pow2(a[Z] - b[Z]);
    } // distance_squared

    double intersect(double xy[2], double const v1[3], double const v2[3], float const threshold=1e-7) {
        auto const d1 = v1[2], d2 = v2[2];
        // Given two lines with normal vectors v1 and v2 and distances to the origin d1 and d2, respectively.
        // Compute their intersection (if any)
        double const det2x2 = v1[0]*v2[1] - v1[1]*v2[0]; // 2d determinant
        xy[0] = 0; xy[1] = 0; // init result (should not be used if false is returned)
        // we look for a point (x,y) on both lines, i.e.
        //         v1 dot xy == d1
        // and simultaneousy
        //         v2 dot xy == d2
        // set result
        double const denom = 1./det2x2;
        if (denom == denom) {
            xy[0] = (v2[1]*d1 - v1[1]*d2)*denom;
            xy[1] = (v1[0]*d2 - v2[0]*d1)*denom;
//          std::printf("# line1 (%g,%g)*(x,y) == %g intersects with line2 (%g,%g)*(x,y) == %g at (%g,%g)\n", v1[0],v1[1], d1, v2[0],v2[1], d2, xy[0],xy[1]);
        } // denom is not NaN
        return std::abs(det2x2);
    } // intersect

    status_t test_plane_balancer(int const nprocs, int const n[3], int const echo=0) {
        status_t stat(0);

        if (nprocs < 1) return stat;

        auto const nall = (n[X])*size_t(n[Y])*size_t(n[Z]);
        assert(nall <= (size_t(1) << 32) && "uint32_t is not long enough!");

        double w8sum_all{0};
        int constexpr W = 3;
        auto const xyzw = new float[nall][4];
        std::vector<double> w8s(nall, 0);

        int const holes = control::get("load_balancer.test.holes", 0.);
        auto const hole_radius_squared = pow2(control::get("load_balancer.test.holes.radius", 8.));

        for (int iz{0}; iz < n[Z]; ++iz) {
        for (int iy{0}; iy < n[Y]; ++iy) {
        for (int ix{0}; ix < n[X]; ++ix) {
            auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
//          assert(uint32_t(iall) == iall && "uint32_t is not long enough!");
            float h{1};
            for (int ih{1 - holes}; ih < holes; ih += 2) {
                auto const x_hole = n[X]*ih/(2.*holes);
                auto const r2 = pow2(ix - .5*n[X] - x_hole) + pow2(iy - .5*n[Y]) + pow2(iz - .5*n[Z]);
                h *= (r2 > hole_radius_squared); // radius_squared
            } //
            float const w8 = 1.f*h; // weight(ix,iy,iz); // WEIGHTS CAN BE INSERTED HERE
            w8s[iall]     = w8;
            w8sum_all    += w8;
            xyzw[iall][W] = w8;
            xyzw[iall][X] = ix;
            xyzw[iall][Y] = iy;
            xyzw[iall][Z] = iz;
        }}} // ix iy iz
        double const longest_possible_distance = std::sqrt(pow2(n[X]) + pow2(n[Y]) + pow2(n[Z]));

        std::vector<double> load(nprocs, 0.0);
        auto const rank_center = new double[nprocs][4];
        bool constexpr compute_rank_centers = true;
        std::vector<uint16_t> owner_rank(nall, no_owner);

        if (echo > 0) std::printf("# %s: distribute %g blocks to %d processes\n\n", __func__, w8sum_all, nprocs);

#ifdef    LOAD_BALANCER_DRAW_SVG
        draw2D.resize(2); draw2D[0] = n[X]; draw2D[1] = n[Y]; // init
#endif // LOAD_BALANCER_DRAW_SVG

        int const echo_rank0 = control::get("load_balancer.test.echo.rank0", 0.); // increase the verbosity for rank0

        ProgressReport timer(__FILE__, __LINE__, 1.5, echo); // update the line every 1.5 seconds
        for (int rank{0}; rank < nprocs; ++rank) { // loop over all ranks serially
            load[rank] = plane_balancer(nprocs, rank, nall, xyzw, w8s.data(), w8sum_all, echo + (0 == rank)*echo_rank0,
                                        rank_center[rank], owner_rank.data());
            timer.report(rank, nprocs);
        } // rank

        analyze_load_imbalance(load.data(), nprocs, echo);

#ifdef    LOAD_BALANCER_DRAW_SVG
        int const nplanes = draw2D.size()/4;
        if (nplanes > 0) {
            // assume that the global array load_balancer::draw2D 
            // is an array of sets of 4 doubles which results from a depth-first traversal of the bisection tree
            assert(2 + 4*nplanes == draw2D.size() && "size of draw2D array must be 2+4*nplanes");
            auto const nx = int(draw2D[0]), ny = int(draw2D[1]);
            if (echo > 2) std::printf("\n# found %d planes for https://editsvgcode.com/\n", nplanes);
            auto const svg_filename = control::get("load_balancer.test.file", "plane_balancer.svg");
            auto const svg = std::fopen(svg_filename, "w");
            if (nullptr != svg) {
                int const comments = control::get("load_balancer.test.svg.comments", 1.);
                std::fprintf(svg, "<!-- SVG code generated by %s -->\n", __FILE__);
                std::fprintf(svg, "<svg viewBox=\"%d %d %d %d\" xmlns=\"http://www.w3.org/2000/svg\">\n", -10, -10, nx + 20, ny + 20);
                double const frame[4][4] = {{1,0,0,-1}, {0,1,0,-1}, {1,0,1.*nx,-1}, {0,1,1.*ny,-1}}; // frame has tree_level=-1
                if (comments > 2) {
                    // plot the frame first (additional to later)
                    std::fprintf(svg, "  <rect width=\"%d\" height=\"%d\" x=\"%g\" y=\"%g\" fill=\"none\" stroke=\"grey\" />\n", nx, ny, -.5, -.5);
                }
                assert(0 == draw2D[5] && "the 1st plane must be the origin");
                std::vector<int> ancestor(32, -1);
                int ip{0}; // plane index
                for (int rank{0}; rank < nprocs; ++rank) { // loop over all ranks serially

                    assert(rank_center[rank][0] < nx); assert(rank_center[rank][1] < ny);
                    assert(rank_center[rank][0] >= 0); assert(rank_center[rank][1] >= 0);
                    // replay the plane_balancer routine branching structure
                    int np{nprocs}, rank_offset{0}, tree_level{0};
                    while (np > 1) {
                        assert(rank_offset + np <= nprocs);
                        int const nhalf[] = {(np + 1) >> 1, np >> 1}; // larger part and smaller part
                        int const i01 = (rank >= rank_offset + nhalf[0]);
                        assert(nhalf[0] + nhalf[1] == np);

                        if (rank == rank_offset) {
                            std::printf("# line #%i level=%d      %d|%d\n", ip, tree_level, nhalf[0],nhalf[1]);
                            assert(ip < nplanes);

                            double const *const v1 = & draw2D.at(ip*4 + 2);

                            assert(tree_level == v1[3]); // ensure matching tree level from component #3
                            assert(tree_level >= 0 && tree_level < 32);
                            ancestor[tree_level] = ip;
                            if (comments > 3) {
                                std::fprintf(svg, "  <!-- I am plane #%i, level=%d, my ancestors are {", ip, tree_level);
                                for (int jp{0}; jp < tree_level; ++jp) {
                                    std::fprintf(svg, "%s%i", jp?",":"", ancestor[jp]);
                                } // jp
                                std::fprintf(svg, "} -->\n");
                            } // 0

                            double points[99][2];
                            int npoints{0}; int ipoint[99];
                            // determine who are my ancestors i.e. which lines are my parents and grandparents and so on.
                            // A line is not supposed to cross its ancestors
                            // The tree has been traversed depth-first
                            for (int jp{tree_level - 1}; jp >= -4; --jp) { // loops over ancestor lines (jp >= 0) and frame (jp in {-1, -2, -3, -4})
                                assert(ancestor[jp] >= 0);
                                double const *const v2 = (jp < 0) ? frame[jp + 4] : & draw2D.at(ancestor[jp]*4 + 2);
                                // compute intersection of the lines
                                auto const intersects = intersect(points[npoints], v1, v2);
                                if (intersects > 1e-12) {
                                    assert(npoints < 99);
                                    ipoint[npoints] = (jp < 0) ? jp : ancestor[jp];
                                    auto const x = points[npoints][0], y = points[npoints][1];
                                    double constexpr eps = 1e-9;
                                    // check if [x, y] are within the border rect [0...nx, 0...ny]
                                    if ((x > -eps) && (x < nx + eps) && (y > -eps) && (y < ny + eps)) {
                                        ++npoints; // accept the point
                                    } // inside rect
                                } // intersects
                            } // jp

                            if (npoints > 1) {
                                bool plot{true};
                                if (npoints > 2) {

                                    double cow0[] = {0, 0, 0, 0};
                                    double cow1[] = {0, 0, 0, 0};
                                    double den0{0}, den1{0};
                                    for (int r{0}; r < nhalf[0]; ++r) {
                                        add_product(cow0, 4, rank_center[r + rank_offset], 1.);
                                        den0 += 1;
                                    } // r
                                    if (den0 > 0) scale(cow0, 3, 1./den0);
                                    for (int r{0}; r < nhalf[1]; ++r) {
                                        add_product(cow1, 4, rank_center[r + rank_offset + nhalf[0]], 1.);
                                        den1 += 1;
                                    } // r
                                    assert(np == den0 + den1);
                                    if (den1 > 0) scale(cow1, 3, 1./den1);
                                    if (comments > 2) std::fprintf(svg, "  <!-- cow0: x=%g y=%g   cow1: x=%g y=%g-->\n", cow0[0], cow0[1], cow1[0], cow1[1]);

                                    double middle[] = {0, 0, 0, 0};
                                    add_product(middle, 3, cow0, 0.5);
                                    add_product(middle, 3, cow1, 0.5);
                                    if (comments > 4) std::fprintf(svg, "  <circle cx=\"%g\" cy=\"%g\" r=\".5\" stroke=\"blue\" />\n", middle[X], middle[Y]);

                                    // decide from which point to which other is the relevant section
                                    // project the middle point onto the plane, assume normalized normal vector (v1[0], v1[1])
                                    auto const plane_dist = v1[2];
                                    auto const project = v1[0]*middle[0] + v1[1]*middle[1] - plane_dist; // distance from the line in 2D
                                    double foot[] = {0, 0};
                                    set(foot, 2, middle);
                                    add_product(foot, 2, v1, -project);
                                    if (comments > 5) std::fprintf(svg, "  <circle cx=\"%g\" cy=\"%g\" r=\".5\" stroke=\"green\" />\n", foot[X], foot[Y]);

                                    // now find the pair of points that lie to the left and to the right of foot on the line (2D)
                                    // and have the minimum distance from each other. The inner product of (p_i - foot)*(p_j - foot) must be negative
                                    double dist2_min{9e299};
                                    int ji[2] = {-1, -1};
                                    for (int ik{1}; ik < npoints; ik++) {
                                        double const foot_i[2] = {points[ik][0] - foot[0], points[ik][1] - foot[1]};
                                        for (int jk{0}; jk < ik; ++jk) { // self-avoiding triangular loop
                                            double const foot_j[2] = {points[jk][0] - foot[0], points[jk][1] - foot[1]};
                                            double const inner = foot_i[0]*foot_j[0] + foot_i[1]*foot_j[1];
                                            double const dist2 = pow2(points[jk][0] - points[ik][0]) + pow2(points[jk][1] - points[ik][1]);
                                            if (comments > 6) std::fprintf(svg, "  <!-- line #%i has %d points, inner(%i,%i)= %g -->\n", ip, npoints, jk,ik, inner);
                                            if (inner < 0 && dist2 < dist2_min) {
                                                dist2_min = dist2;
                                                ji[0] = jk; ji[1] = ik;
                                            } // new extremum found
                                        } // jk
                                    } // ik
                                    if (comments > 3) std::fprintf(svg, "  <!-- line #%i has %d points, min= %g found at %i and %i -->\n",
                                                                                ip, npoints, std::sqrt(dist2_min), ji[0], ji[1]);
                                    if (ji[0] < ji[1]) {
                                        for (int k01{0}; k01 < 2; ++k01) { // loop must run forward!
                                            ipoint[k01]       = ipoint[ji[k01]];
                                            set(points[k01], 2, points[ji[k01]]);
                                        } // k01
                                    } else {
                                        warn("no matching pair found for plane#%i", ip);
                                        plot = false;
                                    }
                                } // npoints > 2
                                if (comments > 2) std::fprintf(svg, "  <!-- line #%i has %d points, take #%i and #%i -->\n", ip, npoints, ipoint[0], ipoint[1]);
                                if (plot) {
                                    std::fprintf(svg, "  <line x1=\"%g\" y1=\"%g\" x2=\"%g\" y2=\"%g\" stroke=\"black\" />\n",
                                                               points[0][0], points[0][1], points[1][0], points[1][1]);
                                } // plot
                            } else { // npoints > 1
                                if (echo > 0) std::printf("# strange case in SVG export: only %d points found\n", npoints);
                            } // npoints > 1

                            ++ip;
                        } // rank_offset

                        // replay the plane_balancer routine branching structure
                        rank_offset += i01*nhalf[0];
                        np = nhalf[i01];
                        ++tree_level;
                    } // while
                } // rank

                if (nplanes != ip) { warn("number of recorded planes %d but replay gave %d", nplanes, ip); }

                // plot the frame
                std::fprintf(svg, "  <rect width=\"%d\" height=\"%d\" x=\"%g\" y=\"%g\" fill=\"none\" stroke=\"grey\" />\n", nx, ny, -.5, -.5);

                if (comments > 0) {
                    std::fprintf(svg, "  <!-- show %d rank centers due to comments > 0, comments=%i -->\n", nprocs, comments);
                    for (int rank{0}; rank < nprocs; ++rank) {
                        auto const *const v = rank_center[rank];
                        std::fprintf(svg, "  <circle cx=\"%g\" cy=\"%g\" r=\"1\" fill=\"none\" stroke=\"red\" />\n", v[X], v[Y]);
                    } // rank
                } else { // comments > 0
                    std::fprintf(svg, "  <!-- do not show %d rank centers due to comments=%i -->\n", nprocs, comments);
                } // comments > 0
                std::fprintf(svg, "</svg>\n\n");
                std::fclose(svg);
                if (echo > 2) std::printf("# SVG file \'%s\' written\n\n", svg_filename);
            } // fopen successful
        } else { // nplanes > 0
            if (echo > 0) std::printf("\n# no planes found, increase to at least +verbosity=%d to perform plane recording\n", echo_record_planes + 1);
        } // nplanes > 0
#endif // LOAD_BALANCER_DRAW_SVG


        if (compute_rank_centers && nprocs > 1) {
            // analyze positions of rank centers
            double mindist2{9e300}; int ijmin[] = {-1, -1};
            double maxdist2{-1.0};  int ijmax[] = {-1, -1};
            simple_stats::Stats<> st2, st1;
            double const wbin = control::get("load_balancer.test.bin.width", 0.25), invbin = 1./wbin;
            int const nbin = 1 + int(longest_possible_distance/wbin);
            std::vector<uint32_t> hist(nbin, 0);
            int np{0}; // counter for the number of processes with a non-zero load
            for (int irank{0}; irank < nprocs; ++irank) {
                if (load[irank] > 0) {
                    ++np;
                    for (int jrank{0}; jrank < nprocs; ++jrank) { // self-avoiding triangular loop
                        if (load[jrank] > 0) {
                            auto const dist2 = distance_squared(rank_center[irank], rank_center[jrank]);
                            if (dist2 > 0 && dist2 < mindist2) { mindist2 = dist2; ijmin[0] = irank; ijmin[1] = jrank; }
                            if (dist2 > maxdist2) { maxdist2 = dist2; ijmax[0] = irank; ijmax[1] = jrank; }
                            auto const dist = std::sqrt(dist2);
                            if (echo > 15) std::printf("# distance-ij is %g\n", dist);
                            int const ibin = dist*invbin; // floor
                            ++hist[std::min(ibin, nbin - 1)];
                            st2.add(dist2);
                            st1.add(dist);
                        } // load
                    } // jrank
                } // load
            } // irank
            auto const maxdist = std::sqrt(maxdist2), mindist = std::sqrt(mindist2);
            if (echo > 1) std::printf("# shortest distance between centers is %g between rank#%i and #%i, longest is %g\n",
                                            mindist, ijmin[0], ijmin[1], maxdist);
            if (echo > 9) std::printf("# longest distance between centers is %g between rank#%i and #%i, shortest is %g\n",
                                            maxdist, ijmax[0], ijmax[1], mindist);
            if (echo > 12) {
                double const denom = 1./pow2(std::max(1, np));
                std::printf("## center-distance histogram, bin width %g\n", wbin);
                for (int ibin{0}; ibin < nbin; ++ibin) {
                    std::printf("%g %g\n", ibin*wbin, hist[ibin]*denom);
                } // ibin
                std::printf("\n\n");
            } // echo
            if (echo > 2) std::printf("# stats: distance [%g, %g +/- %g, %g]\n"
                                      "#        distance^2 [%g, %g +/- %g, %g]\n",
                                      st1.min(), st1.mean(), st1.dev(), st1.max(),
                                      st2.min(), st2.mean(), st2.dev(), st2.max());
        } // compute_rank_centers

        // check masks
        if (1) {
            int strange{0};
            for (int iz{0}; iz < n[Z]; ++iz) {
                for (int iy{0}; iy < n[Y]; ++iy) {
                for (int ix{0}; ix < n[X]; ++ix) {
                    auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
                    auto const owner = owner_rank[iall];
                    strange += (no_owner == owner); // under-assignement
                    if (no_owner == owner) warn("work item %d %d %d has not been assigned to any rank", ix,iy,iz);
                } // ix
                } // iy
            } // iz
            if (strange) warn("strange: %d under-assignments", strange);
        } // 1

        if (1 == n[Z] && echo > 5 && n[X] <= 300 && n[Y] <= 300) {
            std::printf("\n# visualize plane balancer %d x %d on %d processes:%s\n",
                n[Y], n[X], nprocs, (nprocs > 64) ? " (symbols are not unique!)" : "");
            int constexpr iz = 0;
            char const chars[65] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ<>";
            char constexpr mask_char = ' ';
            std::vector<char> line(n[X] + 1, ' '); line.at(n[X]) = '\0';
            for (int iy{0}; iy < n[Y]; ++iy) {
                for(int ix{0}; ix < n[X]; ++ix) {
                    auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
                    auto const owner = owner_rank[iall];
                    auto const c = (xyzw[iall][W] < 1) ? mask_char : chars[owner & 0x3f];
                    line[ix] = c;
                } // ix
                std::printf("# %s\n", line.data());
            } // iy
            std::printf("#\n\n");

#ifdef    HAS_BITMAP_EXPORT
            // create a color map
            int const n3 = std::ceil(std::cbrt(nprocs));      assert(n3*n3*n3 >= nprocs);
            std::vector<uint8_t> bmp_color(n3*n3*n3*4, 0);
            { // scope: generate >=nprocs distinct colors
                int iproc{0};
                for (int cx{1}; cx <= n3; ++cx) {
                    for (int cy{1}; cy <= n3; ++cy) {            // start from 1 to avoid black(0x000000)
                        for (int cz{1}; cz <= n3; ++cz) {
                            bmp_color[iproc*4    ] = (250*cz)/n3;
                            bmp_color[iproc*4 + 1] = (250*cy)/n3; // use 250 instead of 255 to avoid white(0xffffff)
                            bmp_color[iproc*4 + 2] = (250*cx)/n3;
                            ++iproc;
                }}} // cx cy cz
                assert(iproc >= nprocs && "We need enough distinct colors");
            } // scope

            // color the owned region
            std::vector<uint8_t> bmp_data(n[Y]*n[X]*4, 0);
            for (int iy{0}; iy < n[Y]; ++iy) {
                for(int ix{0}; ix < n[X]; ++ix) {
                    auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
                    auto const owner = owner_rank[iall];
                    for (int rgb{0}; rgb < 3; ++rgb) {
                        bmp_data[iall*4 + rgb] = bmp_color[owner*4 + rgb];
                        if (xyzw[iall][W] < 1) bmp_data[iall*4 + rgb] = 255; // white(0xffffff)
                    } // rgb
                } // ix
            } // iy

            // store as image
            char filename[64]; std::snprintf(filename, 64, "load_balancer-n%d-%dx%d", nprocs, n[X], n[Y]);
            if (echo > 5) std::printf("# try to create bitmap file %s with %d colors for %d domains\n", filename, n3*n3*n3, nprocs);
            stat += bitmap::write_bmp_file(filename, bmp_data.data(), n[Y], n[X], -1, 1.f);
#endif // HAS_BITMAP_EXPORT
        } // visualize

        if (echo > 5) {
            // visualize also 3D geometries
            double const p3d = 0.0625; // pseudo-3D: p3d controls how much we go to the right and down for each unit in z-direction
            double const camera[2][4] = {{1,0,p3d, 0}, {0,1,p3d, 0}}; // projection from 3D to 2D, only 3 components used
            std::printf("\n# visualize 3D plane balancer %d x %d x %d on %d processes\n", n[Z], n[Y], n[X], nprocs);
            auto const expl_max = control::get("load_balancer.test.cubes.factor", 2.);
            auto const expl_min = control::get("load_balancer.test.cubes.factor.min", expl_max);
            auto const expl_inc = control::get("load_balancer.test.cubes.factor.inc", 0.0625); // use 0.125 rather than 0.1 for rounding errors

            double cmin[] = {9e99, 9e99}, cmax[] = {-9e99, -9e99};
            for (int iz{0}; iz < n[Z]; ++iz) {
            for (int iy{0}; iy < n[Y]; ++iy) {
            for (int ix{0}; ix < n[X]; ++ix) {
                auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
                auto const rank = owner_rank[iall];
                auto const *const rc = rank_center[rank];
                auto const *const xyz = xyzw[iall];
                double const coords3D[] = {xyz[0] + expl_max*rc[0], xyz[1] + expl_max*rc[1], xyz[2] + expl_max*rc[2]};
                for (int d2{0}; d2 < 2; ++d2) {
                    auto const coord2D = dot_product(3, camera[d2], coords3D);
                    cmax[d2] = std::max(cmax[d2], coord2D);
                    cmin[d2] = std::min(cmin[d2], coord2D);
                } // d2
            }}} // ix iy iz
            if (echo > 5) std::printf("# with +load_balancer.test.cubes.factor=%g  new coords (x in [%g, %g], y in [%g, %g])\n",
                                                                                expl_max, cmin[0], cmax[0], cmin[1], cmax[1]);
            // double cow[] = {0, 0, 0, 0};
            // for (int rank{0}; rank < nprocs; ++rank) {
            //     add_product(cow, 4, rank_center[rank], 1.);
            // } // rank
            // scale(cow, 3, 1./nprocs);

            // create a color map (same as for the BMP export above)
            int const n3 = std::ceil(std::cbrt(nprocs));      assert(n3*n3*n3 >= nprocs);
            std::vector<char> svg_color(n3*n3*n3*8, '\0');
            std::vector<char> face_color(n3*n3*n3*16, '\0'); // 2 darker shades of svg_color
            { // scope: generate n3^3 >= nprocs distinct colors
                int iproc{0};
                for (int cx{1}; cx <= n3; ++cx) {
                    for (int cy{1}; cy <= n3; ++cy) {            // start from 1 to avoid black(0x000000)
                        for (int cz{1}; cz <= n3; ++cz) {
                            uint8_t const red   = (250*cz)/n3;
                            uint8_t const green = (250*cy)/n3; // use 250 instead of 255 to avoid white(0xffffff)
                            uint8_t const blue  = (250*cx)/n3;
                            std::snprintf(&svg_color[iproc*8], 8, "#%.2x%.2x%.2x", red, green, blue);
                            std::snprintf(&face_color[iproc*16    ], 8, "#%.2x%.2x%.2x", uint8_t(red*.875), uint8_t(green*.875), uint8_t(blue*.875));
                            std::snprintf(&face_color[iproc*16 + 8], 8, "#%.2x%.2x%.2x", uint8_t(red*.750), uint8_t(green*.750), uint8_t(blue*.750));
                            ++iproc;
                }}} // cx cy cz
                assert(iproc >= nprocs && "We need enough distinct colors");
            } // scope


            int n_expl_files{0};
            for (double expl{expl_min}; expl <= expl_max; expl += expl_inc) { auto const expl_now = expl;

                char svg_filename[96]; std::snprintf(svg_filename, 64, "plane_balancer-e%6.6d.svg", int(expl_now*1000));
                // std::snprintf(svg_filename, 64, "plane_balancer-n%d-%dx%dx%d-e%6.6d.svg", nprocs, n[X], n[Y], n[Z], int(expl_now*1000));
                auto const svg = std::fopen(svg_filename, "w");
                if (nullptr != svg) {
                    std::fprintf(svg, "<!-- SVG code generated by %s -->\n", __FILE__);
                    std::fprintf(svg, "<svg viewBox=\"%d %d %d %d\" xmlns=\"http://www.w3.org/2000/svg\">\n",
                                    int(cmin[X])-10, int(cmin[Y])-10, int(cmax[X]) + 20, int(cmax[Y]) + 20);

                    for (int iz{0}; iz < n[Z]; ++iz) {
                    for (int iy{0}; iy < n[Y]; ++iy) {
                    for (int ix{0}; ix < n[X]; ++ix) {
                        auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
                        auto const rank = owner_rank[iall];
                        auto const *const rc = rank_center[rank];
                        auto const *const xyz = xyzw[iall];
                        double const coords3D[] = {xyz[0] + expl_now*rc[0], xyz[1] + expl_now*rc[1], xyz[2] + expl_now*rc[2]};
#if       1
                        // plot a circle for each cube
                        double const coords2D[] = {dot_product(3, camera[0], coords3D),
                                                dot_product(3, camera[1], coords3D)};
                        std::fprintf(svg, "  <circle cx=\"%g\" cy=\"%g\" r=\".25\" stroke=\"%s\" fill=\"%s\" />\n", 
                                                coords2D[X], coords2D[Y], &svg_color[rank*8], &svg_color[rank*8]);
#else  // 1
                        // print each cube with 3 polygons -- needs some face ordering ... ToDo
                        int8_t constexpr icube[3][4] = {{7,6,4,5}, {7,5,1,3}, {7,3,2,6}};
                        double cube[8][2];
                        for (int i8{0}; i8 < 8; ++i8) {
                            int const bx = i8 & 0x1, by = (i8 >> 1) & 0x1, bz = (i8 >> 2) & 0x1;
                            double const c8[] = {coords3D[0] + bx, coords3D[1] + by, coords3D[2] + bz};
                            cube[i8][0] = dot_product(3, camera[0], c8);
                            cube[i8][1] = dot_product(3, camera[1], c8);
                        } // i8
                        char const *color = & svg_color[rank*8];
                        for (int face{0}; face < 3; ++face) {
                            auto const *const ic = icube[face];
                            std::fprintf(svg, "  <polygon points=\"");
                            for (int i5{0}; i5 < 5; ++i5) {
                                int const i4 = i5 & 0x3; // i4 goes {0,1,2,3,0}
                                std::fprintf(svg, "%g,%g ", cube[ic[i4]][0], cube[ic[i4]][1]);
                            } // i5
                            std::fprintf(svg, "\" style=\"fill:%s;stroke-width:0\" />\n", color);
                            color = & face_color[rank*16 + face*8];
                        } // face
#endif // 1
                    }}} // ix iy iz

                    std::fprintf(svg, "</svg>\n\n");
                    std::fclose(svg);
                    if (echo > 2) std::printf("# SVG file \'%s\' written for expl=%g\n", svg_filename, expl);
                    ++n_expl_files;
                } // fopen successful
            } // expl
            if (echo > 0) std::printf("# %d SVG files written for expl in [%g, %g] in steps of %g\n", n_expl_files, expl_min, expl_max, expl_inc);
        } // visualize 3D

        delete[] rank_center;
        delete[] xyzw;
        return stat;
    } // test_plane_balancer

    // example 5 processes -->
    //        rank#0      rank#1      rank#2      rank#3      rank#4
    //  take    3/5         3/5         3/5         2/5         2/5
    //  take    2/3         2/3         1/3         1/2         1/2
    //  take    1/2         1/2

    inline double random_between_0_and_1() {
        double constexpr rand_denom = 1./(RAND_MAX + 1.);
        return rand()*rand_denom;
    } // random_between_0_and_1


    status_t test_reference_point_cloud(int const nxyz[3], int const echo=9) {
        auto const maxdist_diagonal = std::sqrt(pow2(nxyz[X]) + pow2(nxyz[Y]) + pow2(nxyz[Z]));
        double const wbin = control::get("load_balancer.test.bin.width", 0.25), invbin = 1./wbin;
        int const nbin = int(maxdist_diagonal/wbin) + 1;
        if (echo > 0) std::printf("# reference point-cloud histogram for %d x %d x %d points, bin width %g, %d bins\n",
                                                                    nxyz[X], nxyz[Y], nxyz[Z],        wbin, nbin);
        std::vector<uint32_t> hist(nbin, 0);
        simple_stats::Stats<> st2, st1;
        for (int iz = 0; iz < nxyz[Z]; ++iz) {
        for (int iy = 0; iy < nxyz[Y]; ++iy) {
        for (int ix = 0; ix < nxyz[X]; ++ix) {
                for (int jz = 0; jz < nxyz[Z]; ++jz) {
                for (int jy = 0; jy < nxyz[Y]; ++jy) {
                for (int jx = 0; jx < nxyz[X]; ++jx) {
    //                    double const rnd[] = {0.5*random_between_0_and_1() - 0.25,
    //                                          0.5*random_between_0_and_1() - 0.25,
    //                                          0.5*random_between_0_and_1() - 0.25};
                        int constexpr rnd[] = {0, 0, 0};
                        double const dist2 = pow2(ix - jx + rnd[X])
                                            + pow2(iy - jy + rnd[Y])
                                            + pow2(iz - jz + rnd[Z]);
                        auto const dist = std::sqrt(dist2);
                        st2.add(dist2);
                        st1.add(dist);
                        if (echo > 15) std::printf("# distance-ij is %g\n", dist);
                        int const ibin = dist*invbin; // floor
                        ++hist[ibin];
                }}} // jx jy jz
        }}} // ix iy iz
        if (echo > 5) {
            double const by_n = 1./(1.*nxyz[X]*nxyz[Y]*nxyz[Z]);
            double const denom = pow2(by_n);
            std::printf("## point-distance histogram, bin width %g\n", wbin);
            for (int ibin = 0; ibin < nbin; ++ibin) {
                auto const radius = ibin*wbin;
                // number of points inside a sphere shell of from radius r - wbin to r is 4/3*pi*(r^3 - (r - wbin)^3)*rho
                // so in the limit of small wbin, this becomes 4*pi*r^2*wbin*rho
                auto const analytical = 4*constants::pi*pow2(radius)*wbin*by_n;
                std::printf("%g %g %g\n", radius, hist[ibin]*denom, analytical);
            } // ibin
            std::printf("\n\n");
        } // echo
        if (echo > 2) std::printf("# stats: distance [%g, %g +/- %g, %g]\n"
                                  "#        distance^2 [%g, %g +/- %g, %g]\n",
                                  st1.min(), st1.mean(), st1.dev(), st1.max(),
                                  st2.min(), st2.mean(), st2.dev(), st2.max());
        return 0;
    } // test_reference_point_cloud


    status_t all_tests(int const echo) {
        status_t stat(0);

        int const nprocs =      control::get("load_balancer.test.nprocs", 53.);
        auto const niso  =      control::get("load_balancer.test.n", 19.);
        int const nxyz[] = {int(control::get("load_balancer.test.nx", niso)), // number of blocks
                            int(control::get("load_balancer.test.ny", niso)),
                            int(control::get("load_balancer.test.nz", 1.))};
        if (echo > 0) std::printf("\n\n# %s start %d x %d x %d = %d with %d MPI processes\n",
                        __func__, nxyz[X], nxyz[Y], nxyz[Z], nxyz[X]*nxyz[Y]*nxyz[Z], nprocs);

        stat += test_plane_balancer(nprocs, nxyz, echo);
//      stat += test_reference_point_cloud(nxyz, echo);
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace load_balancer
