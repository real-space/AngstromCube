#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::sprintf
#include <cassert> // assert
#include <algorithm> // std::max
#include <cstdint> // int8_t, int16_t
#include <numeric> // std::iota

#include "simple_stats.hxx" // ::Stats
#include "status.hxx" // status_t
#include "control.hxx" // ::get
#include "recorded_warnings.hxx" // warn
#include "data_view.hxx" // view2D<T>
#include "mpi_parallel.hxx" // ::allreduce, ::sum, ::min, ::max
#include "global_coordinates.hxx" // ::get

namespace verify_benchmark {

    inline double show_symmetrized(
          view2D<double> const & density // assume to have 4^3 blocks of 4^3 grid points per cubic diamond unit cell
        , int64_t const gid[] // global block identifyer
        , size_t const nblocks // number of blocks on this MPI rank
        , uint8_t const index[][64]
        , size_t const n=40
        , uint8_t const indirection[]=nullptr
        , char const *const text="?"
        , int const echo=0 // log output level
    ) {
        assert(4*4*4 == density.stride());

        std::vector<simple_stats::Stats<double>> st(n);
        std::vector<int32_t> m444(4*4*4, 0);

        for (size_t iblock{0}; iblock < nblocks; ++iblock) {
            double const *const rho = density[iblock]; // a double[4*4*4] array
            uint32_t coords[3]; global_coordinates::get(coords, gid[iblock]);
            auto const i444 = (((coords[2] & 0x3)*4 + (coords[1] & 0x3))*4 + (coords[0] & 0x3)) & 63;
            ++m444[i444];
            auto const i16 = indirection ? indirection[i444] : i444;
            auto const *const idx = index[i16];
         // if (echo > 9) std::printf("# add block [%i %i %i], id= %lli with i16=%i\n", coords[0], coords[1], coords[2], gid[iblock], int(i16));
            for (int i64{0}; i64 < 64; ++i64) { // grid points in a 4x4x4 block
                auto const i40 = idx[i64];
                assert(i40 < n);
                st.at(i40).add(rho[i64]);
            } // i64
        } // iblock

        { // scope: make sure that there is a balanced number of blocks of each unit cell
            mpi_parallel::sum(m444.data(), m444.size()); // usues MPI_COMM_WORLD by default
            int32_t mini{2147483647}, maxi{-1};
            for (int i444{0}; i444 < 4*4*4; ++i444) {
                mini = std::min(mini, m444[i444]);
                maxi = std::max(maxi, m444[i444]);
            } // i444
            if (echo > 3) std::printf("# detected %d = %d complete diamond unit cells\n", mini, maxi);
            assert(mini == maxi); // all values must be the same
        } // scope

        for (int i40{0}; i40 < n; ++i40) {
            mpi_parallel::allreduce(st.at(i40)); // usues MPI_COMM_WORLD by default
        } // i40

        double maxdev{0}, totalsum{0}, totaltim{0};
        if (echo > 1) std::printf("\n## %s\n# i%lu avg avg+dev avg-dev min max dev tim\n", text, n);
        for (int i40{0}; i40 < n; ++i40) {
            auto const & s = st.at(i40);
            auto const avg = s.mean(), dev = s.dev();
            if (echo > 1) std::printf("%d  %g %g %g %g %g %g %lu\n", i40, avg, avg + dev, avg - dev, s.min(), s.max(), dev, s.tim());
            maxdev = std::max(maxdev, dev);
            totalsum += s.sum();
            totaltim += s.tim();
        } // i40
        if (echo > 1) std::printf("# relative deviation from symmetric is %g\n\n", maxdev*totaltim/std::abs(totalsum));

        return maxdev;
    } // show_symmetrized


    inline status_t verify_diamond(
          view2D<double> const & density // assume to have 4^3 blocks of 4^3 grid points per cubic diamond unit cell
        , int64_t const gid[] // global block identifyer
        , size_t const nblocks // number of blocks on this MPI rank
        , float const threshold=1e-5
        , int const echo=0 // log output level
    ) {

        // a 16x16x16 grid with 8 carbon atoms in periodic diamond symmetry has max 40 different field values
        // the following tables were generated by the local script test/understand_diamond_structure.cxx
        // assume that the coarse density has a diamond unit cell with 16x16x16 grid points, i.e. 4x4x4 blocks per unit cell
        uint8_t const index16[64] = {0,1,2,3,4,5,6,7,2,3,0,1,6,7,4,5,8,9,10,11,12,13,14,15,10,11,8,9,14,15,12,13,2,3,0,1,6,7,4,5,0,1,2,3,4,5,6,7,10,11,8,9,14,15,12,13,8,9,10,11,12,13,14,15};
        uint8_t const index40[16][64] = {
        {27,32,34,36,32,35,37,34,34,37,35,32,36,34,32,27,32,35,37,34,35,38,39,37,37,39,38,35,34,37,35,32,34,37,35,32,37,39,38,35,35,38,39,37,32,35,37,34,36,34,32,27,34,37,35,32,32,35,37,34,27,32,34,36},
        {33,28,23,18,31,25,19,12,26,20,14, 8,22,16,11, 5,31,25,19,12,30,24,17,10,29,21,15, 7,26,20,14, 8,26,20,14, 8,29,21,15, 7,30,24,17,10,31,25,19,12,22,16,11, 5,26,20,14, 8,31,25,19,12,33,28,23,18},
        {13, 9, 4, 1, 9, 6, 3, 4, 4, 3, 6, 9, 1, 4, 9,13, 9, 6, 3, 4, 6, 2, 0, 3, 3, 0, 2, 6, 4, 3, 6, 9, 4, 3, 6, 9, 3, 0, 2, 6, 6, 2, 0, 3, 9, 6, 3, 4, 1, 4, 9,13, 4, 3, 6, 9, 9, 6, 3, 4,13, 9, 4, 1},
        { 5,11,16,22, 8,14,20,26,12,19,25,31,18,23,28,33, 8,14,20,26, 7,15,21,29,10,17,24,30,12,19,25,31,12,19,25,31,10,17,24,30, 7,15,21,29, 8,14,20,26,18,23,28,33,12,19,25,31, 8,14,20,26, 5,11,16,22},
        {33,31,26,22,28,25,20,16,23,19,14,11,18,12, 8, 5,31,30,29,26,25,24,21,20,19,17,15,14,12,10, 7, 8,26,29,30,31,20,21,24,25,14,15,17,19, 8, 7,10,12,22,26,31,33,16,20,25,28,11,14,19,23, 5, 8,12,18},
        {18,12, 8, 5,12,10, 7, 8, 8, 7,10,12, 5, 8,12,18,23,19,14,11,19,17,15,14,14,15,17,19,11,14,19,23,28,25,20,16,25,24,21,20,20,21,24,25,16,20,25,28,33,31,26,22,31,30,29,26,26,29,30,31,22,26,31,33},
        { 5, 8,12,18,11,14,19,23,16,20,25,28,22,26,31,33, 8, 7,10,12,14,15,17,19,20,21,24,25,26,29,30,31,12,10, 7, 8,19,17,15,14,25,24,21,20,31,30,29,26,18,12, 8, 5,23,19,14,11,28,25,20,16,33,31,26,22},
        {22,26,31,33,26,29,30,31,31,30,29,26,33,31,26,22,16,20,25,28,20,21,24,25,25,24,21,20,28,25,20,16,11,14,19,23,14,15,17,19,19,17,15,14,23,19,14,11, 5, 8,12,18, 8, 7,10,12,12,10, 7, 8,18,12, 8, 5},
        {33,31,26,22,31,30,29,26,26,29,30,31,22,26,31,33,28,25,20,16,25,24,21,20,20,21,24,25,16,20,25,28,23,19,14,11,19,17,15,14,14,15,17,19,11,14,19,23,18,12, 8, 5,12,10, 7, 8, 8, 7,10,12, 5, 8,12,18},
        {18,12, 8, 5,23,19,14,11,28,25,20,16,33,31,26,22,12,10, 7, 8,19,17,15,14,25,24,21,20,31,30,29,26, 8, 7,10,12,14,15,17,19,20,21,24,25,26,29,30,31, 5, 8,12,18,11,14,19,23,16,20,25,28,22,26,31,33},
        { 5, 8,12,18, 8, 7,10,12,12,10, 7, 8,18,12, 8, 5,11,14,19,23,14,15,17,19,19,17,15,14,23,19,14,11,16,20,25,28,20,21,24,25,25,24,21,20,28,25,20,16,22,26,31,33,26,29,30,31,31,30,29,26,33,31,26,22},
        {22,26,31,33,16,20,25,28,11,14,19,23, 5, 8,12,18,26,29,30,31,20,21,24,25,14,15,17,19, 8, 7,10,12,31,30,29,26,25,24,21,20,19,17,15,14,12,10, 7, 8,33,31,26,22,28,25,20,16,23,19,14,11,18,12, 8, 5},
        {18,23,28,33,12,19,25,31, 8,14,20,26, 5,11,16,22,12,19,25,31,10,17,24,30, 7,15,21,29, 8,14,20,26, 8,14,20,26, 7,15,21,29,10,17,24,30,12,19,25,31, 5,11,16,22, 8,14,20,26,12,19,25,31,18,23,28,33},
        {36,34,32,27,34,37,35,32,32,35,37,34,27,32,34,36,34,37,35,32,37,39,38,35,35,38,39,37,32,35,37,34,32,35,37,34,35,38,39,37,37,39,38,35,34,37,35,32,27,32,34,36,32,35,37,34,34,37,35,32,36,34,32,27},
        {22,16,11, 5,26,20,14, 8,31,25,19,12,33,28,23,18,26,20,14, 8,29,21,15, 7,30,24,17,10,31,25,19,12,31,25,19,12,30,24,17,10,29,21,15, 7,26,20,14, 8,33,28,23,18,31,25,19,12,26,20,14, 8,22,16,11, 5},
        { 1, 4, 9,13, 4, 3, 6, 9, 9, 6, 3, 4,13, 9, 4, 1, 4, 3, 6, 9, 3, 0, 2, 6, 6, 2, 0, 3, 9, 6, 3, 4, 9, 6, 3, 4, 6, 2, 0, 3, 3, 0, 2, 6, 4, 3, 6, 9,13, 9, 4, 1, 9, 6, 3, 4, 4, 3, 6, 9, 1, 4, 9,13}};
        int16_t const multiplicity40[40] = {32,32,32,96,96,96,96,96,192,96,96,96,192,32,192,96,96,96,96,192,192,96,96,96,96,192,192,32,96,96,96,192,96,96,96,96,32,96,32,32};

        { // scope: perform some sanity checks on the number tables above
            std::vector<int> mult40(40, 0);
            for (int j16{0}; j16 < 16; ++j16) {
                int count{0};
                for (int i64{0}; i64 < 64; ++i64) {
                    count += (index16[i64] == j16);
                    ++mult40.at(index40[j16][i64]);
                } // i64
                assert(4 == count);
            } // j16
            size_t mult_sum{0};
            for (int i40{0}; i40 < 40; ++i40) {
                if (echo > 27) std::printf("# %s: multiplicity[%i]= %d found %d\n", __func__, i40, multiplicity40[i40], mult40[i40]);
                assert(multiplicity40[i40] == 4*mult40[i40]);
                mult_sum += size_t(mult40[i40]);
            } // m
            assert(16*16*16 == 4*mult_sum);
        } // scope
 
        auto const reldev = show_symmetrized(density, gid, nblocks, index40, 40, index16, "show value distribution for diamond symmetry", echo);

        auto const symmetry_broken = (reldev > threshold);
        if (echo > 1) std::printf("# diamond symmetry %s (reldev= %g threshold= %g)\n", symmetry_broken?"false":"true", reldev, threshold);
        return symmetry_broken ? 1 : 0;
    } // verify_diamond


    inline status_t verify(
          view2D<double> const & density // or local effective potential
        , int64_t const gid[]
        , size_t const nblocks
        , int const echo=0 // log output level
    ) {
        auto const verify_benchmark_diamond = "verify.benchmark.diamond";
        auto const dia = control::get(verify_benchmark_diamond, 0.);
        if (dia) {
            assert(density.stride() == 4*4*4);
            return verify_diamond(density, gid, nblocks, dia, echo);
        } else {
            if (echo > 9) std::printf("# %s: no verification found, activate e.g. +%s=1\n", __func__, verify_benchmark_diamond);
            return 0;
        }
    } // verify


    inline status_t all_tests(int const echo=0) {
        size_t const nblocks = 64;
        view2D<double> density(nblocks, 4*4*4, 1.);
        int64_t gid[nblocks]; std::iota(gid, gid + nblocks, 0);
        return verify(density, gid, nblocks, echo);
    } // all_tests

} // namespace verify_benchmark
