#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cmath> // std::exp
#include <vector> // std::vector<T>

#include "dyadic_plan.hxx" // dyadic_plan_t
#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_memory.hxx" // get_memory, free_memory, dim3, real_t_name
#include "green_sparse.hxx" // ::sparse_t<>
#include "inline_math.hxx" // pow2, pow3
#include "sho_tools.hxx" // ::nSHO, ::n2HO, ::n1HO
#include "constants.hxx" // ::sqrtpi
#include "green_parallel.hxx" // ::rank, ::size, ::dyadic_exchange
#include "recorded_warnings.hxx" // warn

#ifndef HAS_NO_CUDA
    #include <cuda/std/complex> // std::complex
    #define std__complex cuda::std::complex
#else
    #include <complex> // std::complex
    #define std__complex std::complex
#endif // HAS_NO_CUDA

namespace green_dyadic {

    int constexpr X = 0, Y = 1, Z = 2;

    template <typename real_t>
    float __host__ __device__
    Hermite_polynomials_1D(
          real_t (*const __restrict__ H1D)[3][4] // result H1D[nu][dir][i4]
        , float  (*const __restrict__ xi_squared)[4] // distance^2 xi_squared[dir][i4]
        , int    const ivec // thread index used for the vectorization (at least 3*4==12 threads must run)
        , int    const lmax
        , double const xyza[3+1] // atom image position in [0],[1],[2], sigma^{-1/2} in [3]
        , float  const xyzc[3]   // position of the target block (in units of 4*grid spacing)
        , double const hxyz[3+1] // grid spacings in [0],[1],[2], projection radius in [3]
    )
      // Evaluate the 1D Hermite-Gauss functions on a 3D block of 4^3 grid points in parallel
    {

        double const R2_projection = hxyz[3]*hxyz[3]; // projection radius can be controlled from outside

        if (lmax < 0) return R2_projection; // as float

        int constexpr l2b = 2;
        int constexpr n4 = 1 << l2b; // Block edge length is 2^l2b
        assert(4 == n4);
        if (ivec < 3*n4) { // use only the lowest 3*n4 threads

            // idea: we could use 60 of 64 threads to make M instances of the H1D-functions at slightly shifted positions
            // and then average over these M in order to filter out the highest frequency components.
            // this is along the line of the double grid technique with the refinement or can also be viewed
            // as an integration of a Haar-wavelet shaped basis function

            int const idir = ivec >> l2b; // divide by n4 --> 0:x, 1:y, 2:z-direction
            assert(idir >= 0); assert(idir < 3);
            int const i4   = ivec - n4*idir; // remainder --> 0, 1, 2 or 3 in the case l2b==2
            assert(i4 >= 0); assert(i4 < n4);

            double const grid_spacing  = hxyz[idir]; // load a float and convert to double
            double const cube_position = xyzc[idir]; // position of the target block, lower left front corner
            double const atom_position = xyza[idir]; // load atomic coordinate
            double const sigma_inverse = xyza[3]*xyza[3]; // sigma_inverse = inverse of the Gaussian width sigma

            double const xi = sigma_inverse*(grid_spacing*(cube_position*n4 + i4 + 0.5) - atom_position); // center of the grid cell // 6 flop
            double const xi2 = xi*xi; // 1 flop
            xi_squared[idir][i4] = xi2; // export the square of the distance in units of sigma
            double const H0 = (xi2 < R2_projection) ? std::exp(-0.5*xi2) : 0; // Gaussian envelope function // ? flop

            // Discuss:   We could export xi*xi as an array double (or maybe only float) xi2[3][4] and compare
            //              xi2[0][x] + xi2[1][y] + xi2[2][z] < R2_projection
            //            hence omitting the data transfer of BitMask entirely and save preparation code complexity


            H1D[0][idir][i4] = H0; // store Gaussian H_0

            double Hnp1, Hn{H0}, Hnm1{0};
//          std::printf("# %s idir=%d xi= %g H= %g", __func__, idir, xi, H0);
            for (int nu = 0; nu < lmax; ++nu) {
                double const nuhalf = 0.5 * nu; // nu/2. // 1 flop
                // use the two step recurrence relation to create the 1D Hermite polynomials
                // H1D[nu+1][idir][i4] = xi * H1D[nu][idir][i4] - nu/2. * H1D[nu-1][idir][i4];
                Hnp1 = xi * Hn - nuhalf * Hnm1; // see snippets/3Hermite.F90 for the derivation of this // 3 flop
                H1D[nu + 1][idir][i4] = Hnp1; // store
//              std::printf(" %g", Hnp1);
                // rotate registers for the next iteration
                Hnm1 = Hn; Hn = Hnp1; // ordering is important here
            } // nu
//          std::printf("\n");

            // Warning!
            // The Hermite polynomials are not normalized. To normalize w.r.t. the 2-norm,
            // we must divide by sqrt{sqrt{pi}*sigma*factorial(nu)} in an additional loop.
            // However, this is not necessary for their orthogonality

        } // ivec < 3*n4

        // here, start the tree reduction over the 5 elements for the grid-refinement

        return R2_projection; // as float
    } // Hermite_polynomials_1D


    int constexpr Lmax_default=7; // reduce this to lower the GPU register and shared memory usage

    template <typename real_t, int R1C2=2, int Noco=1, int Lmax=Lmax_default>
    void __global__ SHOprj( // launch SHOprj<real_t,R1C2,Noco> <<< {natoms, nrhs, 1}, {Noco*64, Noco, R1C2} >>> (...);
#ifdef HAS_NO_CUDA
          dim3 const & gridDim, dim3 const & blockDim,
#endif // HAS_NO_CUDA
          real_t         (*const __restrict__ Cpr)[R1C2][Noco   ][Noco*64] // result: projection coefficients, layout[natomcoeffs*nrhs][R1C2][Noco   ][Noco*64]
        , real_t   const (*const __restrict__ Psi)[R1C2][Noco*64][Noco*64] // input:  Green function,               layout[ncubes*nrhs][R1C2][Noco*64][Noco*64]
        , green_sparse::sparse_t<> const (*const __restrict__ sparse)
        , double   const (*const __restrict__ AtomPos)[3+1] // atomic positions [0],[1],[2], decay parameter [3]
        , int8_t   const (*const __restrict__ AtomLmax) // SHO basis size [iatom]
        , uint32_t const (*const __restrict__ AtomStarts) // prefix sum over nSHO(AtomLmax[:])
        , uint32_t const (*const __restrict__ irow_of_inzb) // row index of the Green function as a function of the non-zero index
        , float    const (*const __restrict__ CubePos)[3+1] // only [0],[1],[2] used, CubePos[irow][0:3]
        , double   const (*const __restrict__ hGrid) // grid spacings in [0],[1],[2], projection radius in [3]
    )
      // Compute the projection coefficients of wave/Green functions with atom-centered SHO-bases
    {
        assert(1       ==  gridDim.z);
        assert(Noco*64 == blockDim.x);
        assert(Noco    == blockDim.y);
        assert(R1C2    == blockDim.z);

        int const nrhs   = gridDim.y;

        __shared__ double hgrid[3+1]; // grid spacings

#ifndef HAS_NO_CUDA
        if (threadIdx.x < 4) hgrid[threadIdx.x] = hGrid[threadIdx.x]; // load grid spacing into shared memory
        int const irhs  = blockIdx.y;  // in [0, nrhs)
        int const iatom = blockIdx.x;  // in [0, natoms)
#else  // HAS_NO_CUDA
        set(hgrid, 4, hGrid);
        int const natoms = gridDim.x;
        for (int irhs  = 0; irhs < nrhs; ++irhs)
        for (int iatom = 0; iatom < natoms; ++iatom)
#endif // HAS_NO_CUDA
        { // block loops

        auto const bsr_of_iatom = sparse[irhs].rowStart();

#ifndef HAS_NO_CUDA
        if (bsr_of_iatom[iatom] >= bsr_of_iatom[iatom + 1]) return; // empty range in the sparse matrix, early return
#endif // HAS_NO_CUDA

        auto const inzb_of_bsr = sparse[irhs].colIndex();

        int const lmax = AtomLmax[iatom];
        if (lmax > Lmax) std::printf("# %s Error: lmax= %d but max. Lmax= %d, iatom=%d\n", __func__, lmax, Lmax, iatom);
        assert(lmax <= Lmax);

        __shared__ real_t H1D[Lmax + 1][3][4]; // non-normalized orthogonal 1-dimensional Hermite Gauss functions
        __shared__ float     xi_squared[3][4]; // distance along one Cartesian direction squared
        __shared__ double xyza[3+1]; // atom position
        __shared__ float  xyzc[3];   // cube position

#ifndef HAS_NO_CUDA
        // stage the atomic position of the atom treated in this CUDA block
        if (threadIdx.x < 4) xyza[threadIdx.x] = AtomPos[iatom][threadIdx.x]; // coalesced load
        int const reim = threadIdx.z; // in [0, R1C2)
        int const spin = threadIdx.y; // in [0, Noco)
        int const j    = threadIdx.x; // in [0, Noco*64)
#else // HAS_NO_CUDA
//      std::printf("# %s Here\n", __func__);
        set(xyza, 4, AtomPos[iatom]);
        for (int reim = 0; reim < R1C2; ++reim)
        for (int spin = 0; spin < Noco; ++spin)
        for (int j = 0; j < Noco*64; ++j)
#endif // HAS_NO_CUDA
        { // thread loops

        real_t czyx[sho_tools::nSHO(Lmax)]; // get nSHO accumulator registers
        for (int sho = 0; sho < sho_tools::nSHO(Lmax); ++sho) {
            czyx[sho] = 0; // init accumulator registers
        } // sho

        __syncthreads(); // sync necessary?

        // in this loop over bsr we accumulate atom projection coefficients over many cubes
        for (auto bsr = bsr_of_iatom[iatom]; bsr < bsr_of_iatom[iatom + 1]; ++bsr) {

            __syncthreads();

            auto const inzb = inzb_of_bsr[bsr]; // non-zero index of the Green function
            auto const irow = irow_of_inzb[inzb]; // row index of the Green function

            __shared__ float R2_proj;
#ifndef HAS_NO_CUDA
            if (threadIdx.x < 3) xyzc[threadIdx.x] = CubePos[irow][threadIdx.x]; // load into shared memory
            __syncthreads();

            // generate Hx, Hy, Hz up to lmax inside the 4^3 cube

            if (0 == threadIdx.y && 0 == threadIdx.z) // sufficient to be executed by the first 12 threads of a block
                R2_proj = Hermite_polynomials_1D(H1D, xi_squared, j, lmax, xyza, xyzc, hgrid);
            // how many times executed? nrhs * sparse[irhs].nNonzeros()

            __syncthreads();
#else  // HAS_NO_CUDA
            set(xyzc, 3, CubePos[irow]);

            for (int jj = 0; jj < 12; ++jj) {
                R2_proj = Hermite_polynomials_1D(H1D, xi_squared, jj, lmax, xyza, xyzc, hgrid);
            } // jj
#endif // HAS_NO_CUDA


            for (int z = 0; z < 4; ++z) { // loop over real-space grid in z-direction
                auto const d2z = xi_squared[Z][z] - R2_proj;
                if (d2z < 0) {

                    real_t byx[sho_tools::n2HO(Lmax)];
                    for (int iyx = 0; iyx < sho_tools::n2HO(lmax); ++iyx) { // loop over the combined index of ix and iy
                        byx[iyx] = 0; // init
                    } // iyx

                    for (int y = 0; y < 4; ++y) { // loop over real-space grid in y-direction
                        auto const d2yz = xi_squared[Y][y] + d2z;
                        if (d2yz < 0) {

                            real_t ax[Lmax + 1];
                            for (int ix = 0; ix <= lmax; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                ax[ix] = 0; // init
                            } // ix

                            // __unroll__
                            for (int x = 0; x < 4; ++x) { // loop over real-space grid in x-direction
                                auto const d2xyz = xi_squared[X][x] + d2yz;
                                if (d2xyz < 0) {
                                    int const xyz = (z*4 + y)*4 + x;
                                    real_t const ps = Psi[inzb][reim][spin*64 + xyz][j]; // load from global memory
                                    for (int ix = 0; ix <= lmax; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                        auto const Hx = H1D[ix][X][x]; // load Hx from shared memory
                                        ax[ix] += ps * Hx; // FMA: 2 flop * 4**3 * (L+1)
                                    } // ix
                                } // inside mask
                            } // x

                            for (int iy = 0; iy <= lmax; ++iy) { // loop over the 2nd Cartesian SHO quantum number
                                int const iyx0 = (iy*(2*lmax + 3 - iy)) >> 1;
                                auto const Hy = H1D[iy][Y][y]; // load Hy from shared memory
                                for (int ix = 0; ix <= lmax - iy; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                    byx[iyx0 + ix] += ax[ix] * Hy; // FMA: 2 flop * 4**2 * ((L+1)*(L+2))/2
                                } // ix
                            } // iy

                        } // inside mask
                    } // y

                    for (int sho = 0, iz = 0; iz <= lmax; ++iz) { // loop over the 3rd Cartesian SHO quantum number
                        auto const Hz = H1D[iz][Z][z]; // load Hz from shared memory
                        for (int iy = 0; iy <= lmax - iz; ++iy) { // loop over the 2nd Cartesian SHO quantum number
                            int const iyx0 = (iy*(2*lmax + 3 - iy)) >> 1;
                            for (int ix = 0; ix <= lmax - iz - iy; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                czyx[sho] += byx[iyx0 + ix] * Hz; // FMA: 2 flop * 4**1 * ((L+1)*(L+2)*(L+3))/6
                                ++sho;
                            } // ix
                        } // iy
                    } // iz

                } // inside mask
            } // z
            // now how many flop have been done?
            // FMA*( 4**3 * (L+1) + 4**2 * ((L+1)*(L+2))/2 + 4**1 * ((L+1)*(L+2)*(L+3))/6 )
            // so for L=5 this is 944 FMAs = 1888 flop

        } // bsr

        { // scope: store accumulators
            auto const a0 = AtomStarts[iatom];
            for (int sho = 0; sho < sho_tools::nSHO(lmax); ++sho) {
                Cpr[(a0 + sho)*nrhs + irhs][reim][spin][j] = czyx[sho]; // 0 flop
            } // sho
        } // scope

        }} // thread loops and block loops

    } // SHOprj

    template <typename real_t, int R1C2=2, int Noco=1>
    void __host__ SHOprj_driver(
          real_t         (*const __restrict__ Cpr)[R1C2][Noco   ][Noco*64] // result: projection coefficients
        , real_t   const (*const __restrict__ Psi)[R1C2][Noco*64][Noco*64] // input: Green function
        , double   const (*const __restrict__ AtomPos)[3+1] // atomic positions [0],[1],[2], decay parameter [3]
        , int8_t   const (*const __restrict__ AtomLmax) // SHO basis size [iatom]
        , uint32_t const (*const __restrict__ AtomStarts) // prefix sum over nSHO(AtomLmax[:])
        , uint32_t const natoms // number of atom images
        , green_sparse::sparse_t<> const (*const __restrict__ sparse)
        , uint32_t const (*const __restrict__ RowIndexCube) // row index of the Green function as a function of the non-zero index
        , float    const (*const __restrict__ CubePos)[3+1] // only [0],[1],[2] used
        , double   const (*const __restrict__ hGrid) // grid spacings in [0],[1],[2], projection radius in [3]
        , int      const nrhs // number of block columns in the Green function, max due to data type 2^16, due to GPU launch 2^16-1
        , int const echo=0
    ) {
        if (natoms*nrhs < 1) return;
        dim3 const gridDim(natoms, nrhs, 1), blockDim(Noco*64, Noco, R1C2);
        if (echo > 3) std::printf("# %s<%s,R1C2=%d,Noco=%d> <<< {natoms=%d, nrhs=%d, 1}, {%d, Noco=%d, R1C2=%d} >>>\n",
                            __func__, real_t_name<real_t>(), R1C2, Noco, natoms, nrhs, Noco*64, Noco, R1C2);
        SHOprj<real_t,R1C2,Noco> // launch <<< {natoms, nrhs, 1}, {Noco*64, Noco, R1C2} >>>
#ifndef HAS_NO_CUDA
              <<< gridDim, blockDim >>> (
#else  // HAS_NO_CUDA
                ( gridDim, blockDim,
#endif // HAS_NO_CUDA
               Cpr, Psi, sparse, AtomPos, AtomLmax, AtomStarts, RowIndexCube, CubePos, hGrid);
    } // SHOprj_driver


    // maybe this could be useful?
//     union mask64_t {
//         int64_t i;
//         uint16_t u[4];
//     }; // mask64_t

    template <typename real_t, int R1C2=2, int Noco=1, int Lmax=Lmax_default>
    void __global__ SHOadd( // launch SHOadd<real_t,R1C2,Noco> <<< {nnzb, 1, 1}, {Noco*64, Noco, R1C2} >>>
#ifdef HAS_NO_CUDA
          dim3 const & gridDim, dim3 const & blockDim,
#endif // HAS_NO_CUDA
          real_t         (*const __restrict__ Psi)[R1C2][Noco*64][Noco*64] // result: Green function to modify,       layout[nnzb][R1C2][Noco*64][Noco*64]
        , real_t   const (*const __restrict__ Cad)[R1C2][Noco   ][Noco*64] // input: addition coefficients, layout[natomcoeffs*nrhs][R1C2][Noco   ][Noco*64]
        , uint32_t const (*const __restrict__ bsr_of_inzb)
        , uint32_t const (*const __restrict__ iatom_of_bsr)
        , uint32_t const (*const __restrict__ irow_of_inzb) // row indices of the Green function
        , uint16_t const (*const __restrict__ irhs_of_inzb) // col indices of the Green function
        , double   const (*const __restrict__ AtomPos)[3+1] // atomic positions [0],[1],[2] and decay parameter [3]
        , int8_t   const (*const __restrict__ AtomLmax) // SHO basis size [iatom]
        , uint32_t const (*const __restrict__ AtomStarts) // prefix sum over nSHO(AtomLmax[:])
        , float    const (*const __restrict__ CubePos)[3+1] // only [0],[1],[2] used
        , double   const (*const __restrict__ hGrid) // grid spacings in [0],[1],[2], projection radius [3]
        , int      const nrhs // number of block columns in the Green function
    )
      // Add linear combinations of SHO-basis function to wave/Green functions
    {
        assert(1       ==  gridDim.y);
        assert(1       ==  gridDim.z);
        assert(Noco*64 == blockDim.x);
        assert(Noco    == blockDim.y);
        assert(R1C2    == blockDim.z);

        __shared__ double hgrid[3+1]; // grid spacings

#ifndef HAS_NO_CUDA
        if (threadIdx.x < 4) hgrid[threadIdx.x] = hGrid[threadIdx.x]; // load grid spacing into shared memory
        int const inzb = blockIdx.x;  // in [0, nnzb)
#else  // HAS_NO_CUDA
        set(hgrid, 4, hGrid);
        for (int inzb = 0; inzb < gridDim.x; ++inzb)
#endif // HAS_NO_CUDA
        { // block loop

#ifndef HAS_NO_CUDA
        if (bsr_of_inzb[inzb] >= bsr_of_inzb[inzb + 1]) return; // empty range in the sparse matrix, early return
#endif // HAS_NO_CUDA

        auto const irow = irow_of_inzb[inzb]; // Green function row index
        auto const irhs = irhs_of_inzb[inzb]; // Green function column index

        __shared__ real_t H1D[1 + Lmax][3][4]; // non-normalized orthogonal 1-dimensional Hermite Gauss functions
        __shared__ float     xi_squared[3][4]; // distance along one Cartesian direction squared
        __shared__ double xyza[3+1]; // atom position
        __shared__ float  xyzc[3+1]; // cube position

#ifndef HAS_NO_CUDA
        // stage the atomic position of the atom treated in this CUDA block
        if (threadIdx.x < 4) xyzc[threadIdx.x] = CubePos[irow][threadIdx.x]; // coalesced load
        int const reim = threadIdx.z; // in [0, R1C2)
        int const spin = threadIdx.y; // in [0, Noco)
        int const j    = threadIdx.x; // in [0, Noco*64)
#else  // HAS_NO_CUDA
        set(xyzc, 4, CubePos[irow]);
//      std::printf("# Here %s\n", __func__);
        for (int reim = 0; reim < R1C2; ++reim)
        for (int spin = 0; spin < Noco; ++spin)
        for (int j = 0; j < Noco*64; ++j)
#endif // HAS_NO_CUDA
        { // thread loops

        real_t czyx[4][4][4]; // get 64 accumulator registers --> 128 registers when real_t==double
        // __unroll__
        for (int z = 0; z < 4; ++z) {
            // __unroll__
            for (int y = 0; y < 4; ++y) {
                // __unroll__
                for (int x = 0; x < 4; ++x) {
                    czyx[z][y][x] = 0; // init accumulator registers
                } // x
            } // y
        } // z


        int64_t mask_all{0}; // mask accumulator that tells which grid points have to be updated

        // in this loop over bsr we accumulate atom projector functions over many atoms
        for (auto bsr = bsr_of_inzb[inzb]; bsr < bsr_of_inzb[inzb + 1]; ++bsr) {

            __syncthreads();
            auto const iatom = iatom_of_bsr[bsr];

            int const lmax = AtomLmax[iatom];
            auto const a0  = AtomStarts[iatom];
//          if (lmax > Lmax) std::printf("# %s Error: lmax=%d but max. Lmax=%d, iatom=%d, bsr=%d, inzb=%d\n", __func__, lmax, Lmax, iatom, bsr, inzb);
            assert(lmax <= Lmax);

            __shared__ float R2_proj;
#ifndef HAS_NO_CUDA
            if (threadIdx.x < 4) xyza[threadIdx.x] = AtomPos[iatom][threadIdx.x]; // coalesced load

            __syncthreads();

            // generate Hx, Hy, Hz up to Lmax inside the 4^3 cube
            if (0 == threadIdx.y && 0 == threadIdx.z) // sufficient to be executed by the first 12 threads of a block
                R2_proj = Hermite_polynomials_1D(H1D, xi_squared, j, lmax, xyza, xyzc, hgrid);

#else  // HAS_NO_CUDA
            set(xyza, 4, AtomPos[iatom]);

            for (int jj = 0; jj < 12; ++jj) {
                R2_proj = Hermite_polynomials_1D(H1D, xi_squared, jj, lmax, xyza, xyzc, hgrid);
            } // jj
#endif // HAS_NO_CUDA

            int sho{0};
            for (int iz = 0; iz <= lmax; ++iz) { // executes (L+1) times

                real_t byx[4][4]; // get 16 accumulator registers --> 32 registers if double
                // __unroll__
                for (int y = 0; y < 4; ++y) {
                    // __unroll__
                    for (int x = 0; x < 4; ++x) {
                        byx[y][x] = 0; // init
                    } // x
                } // y

                for (int iy = 0; iy <= lmax - iz; ++iy) { // executes (L+1)*(L+2)/2 times

                    real_t ax[4] = {0, 0, 0, 0}; // get 4 accumulator registers --> 8 32bit GPU registers if double

                    for (int ix = 0; ix <= lmax - iz - iy; ++ix) { // executes (L+1)*(L+2)*(L+3)/6 times
                        real_t const ca = Cad[(a0 + sho)*nrhs + irhs][reim][spin][j]; // load Noco*64 numbers from global memory
                        ++sho;
                        // __unroll__
                        for (int x = 0; x < 4; ++x) {
                            real_t const Hx = H1D[ix][X][x]; // load Hx from shared memory
                            ax[x] += Hx * ca; // FMA: 2 flop * 4 * (L+1)*(L+2)*(L+3)/6
                        } // x
                    } // ix

                    // __unroll__
                    for (int y = 0; y < 4; ++y) {
                        real_t const Hy = H1D[iy][Y][y]; // load Hy from shared memory
                        // __unroll__
                        for (int x = 0; x < 4; ++x) {
                            byx[y][x] += Hy * ax[x]; // FMA: 2 flop * 4**2 * (L+1)*(L+2)/2
                        } // x
                    } // y

                } // iy

                int64_t mask_atom{0}; // occupy 2 GPU registers
                // __unroll__
                for (int z = 0; z < 4; ++z) {
                    real_t const Hz = H1D[iz][Z][z]; // load Hz from shared memory
                    auto const d2z = xi_squared[Z][z] - R2_proj;
                    if (d2z < 0) {
                        // __unroll__
                        for (int y = 0; y < 4; ++y) {
                            auto const d2yz = xi_squared[Y][y] + d2z;
                            if (d2yz < 0) {
                                // __unroll__
                                for (int x = 0; x < 4; ++x) {
                                    auto const d2xyz = xi_squared[X][x] + d2yz;
                                    if (d2xyz < 0) {
                                        czyx[z][y][x] += Hz * byx[y][x]; // FMA: 2 flop * 4**3 * (L+1) are these real flops?
                                        int const xyz = (z*4 + y)*4 + x;
                                        mask_atom |= (int64_t(1) << xyz);
                                    } // d2xyz < 0
                                } // x
                            } // d2yz < 0
                        } // y
                    } // d2z < 0
                } // z
                mask_all |= mask_atom; // set mask bits for modified grid points to 1

            } // iz
            // now how many flop have been done?
            // FMA*( 4 * ((L+1)*(L+2)*(L+3))/6 + 4**2 * ((L+1)*(L+2))/2 + 4**3 * (L+1) )
            // so for Lmax=5 this is 944 FMAs = 1888 flop

        } // bsr

        if (mask_all) {
            for (int z = 0; z < 4; ++z) {
                for (int y = 0; y < 4; ++y) {
                    for (int x = 0; x < 4; ++x) {
                        int const xyz = (z*4 + y)*4 + x;
                        if ((mask_all >> xyz) & 0x1) { // probe the rightmost bit
                            Psi[inzb][reim][spin*64 + xyz][j] += czyx[z][y][x]; // up to 64 negligible flop
                            // possibility to use atomic add here
                        } // mask
                    } // x
                } // y
            } // z
        } // mask_all

        }} // thread loops and block loop

    } // SHOadd


    template <typename real_t, int R1C2=2, int Noco=1>
    void __host__ SHOadd_driver(
          real_t         (*const __restrict__ Psi)[R1C2][Noco*64][Noco*64] // result: Green functions to modify
        , real_t   const (*const __restrict__ Cad)[R1C2][Noco   ][Noco*64] // input: addition coefficients
        , double   const (*const __restrict__ AtomPos)[3+1] // atomic positions [0],[1],[2], decay parameter [3]
        , int8_t   const (*const __restrict__ AtomLmax) // SHO basis size [iatom]
        , uint32_t const (*const __restrict__ AtomStarts) // prefix sum over nSHO(AtomLmax[:])
        , uint32_t const (*const __restrict__ RowStartCubes) // rows==cubes
        , uint32_t const (*const __restrict__ ColIndexAtoms) // cols==atoms
        , uint32_t const (*const __restrict__ RowIndexCubes) // row index of the Green function [nnzb]
        , uint16_t const (*const __restrict__ ColIndexCubes) // col index of the Green function [nnzb]
        , float    const (*const __restrict__ CubePos)[3+1] // only [0],[1],[2] used
        , double   const (*const __restrict__ hGrid) // grid spacings in [0],[1],[2], projection radius in [3]
        , int      const nnzb // == number of all non-zero blocks in the Green function
        , int      const nrhs // == number of columns in the Green function
        , int const echo=0
    ) {
        if (nnzb < 1) return;
        if (nullptr == Psi) return;
        dim3 const gridDim(nnzb, 1, 1), blockDim(Noco*64, Noco, R1C2);
        if (echo > 3) std::printf("# %s<%s,R1C2=%d,Noco=%d> <<< {nnzb=%d, 1, 1}, {%d, Noco=%d, R1C2=%d} >>>\n",
                            __func__, real_t_name<real_t>(), R1C2, Noco, nnzb, Noco*64, Noco, R1C2);
        SHOadd<real_t,R1C2,Noco> // launch {nnzb, 1, 1}, {Noco*64, Noco, R1C2} >>>
#ifndef HAS_NO_CUDA
              <<< gridDim, blockDim >>> (
#else  // HAS_NO_CUDA
                ( gridDim, blockDim,
#endif // HAS_NO_CUDA
               Psi, Cad, RowStartCubes, ColIndexAtoms, RowIndexCubes, ColIndexCubes, AtomPos, AtomLmax, AtomStarts, CubePos, hGrid, nrhs);
    } // SHOadd_driver


    // For the Green function data layout X[nnzb][R1C2][Noco*64][Noco*64]
    // we need an additional level of parallelism
    // SHOprj: p_{ai}[64_i] * X[nnzb][R1C2][Noco*64_i][Noco*64_j] summed over 64_i -> c[nai*ncols][R1C2][Noco][Noco*64_j]
    // accumulating over Green function elements inzb that belong to the same icol, icol in [0, ncols).
    // SHOadd: Y[nnzb][R1C2][Noco*64_i][Noco*64_j] += p_{ai}[64_i] * d[nai*ncols][R1C2][Noco][Noco*64_j] summed over ai
    // so we need to introduce the information which ColIndex[inzb]

    // for SHOprj: RowStart needs another index [icol] which can be block-parallelised
    //             and ColIndex needs to encrypt the inzb index
    //             furthermore, we need to have an array to get the cube index from the inzb
    //             unless cube-position is replicated --> indirection latency vs memory capacity

    // what is the best data layout for c and d?
    // We need to multiply with atom-matrices
    // --> c/d[nai][ncols] instead of [ncols][nai]

    // how many instances are compiled?
    //  real_t=float/double --> 2x
    //  <R1C2=1, Noco=1>
    //  <R1C2=2, Noco=1>
    //  <R1C2=2, Noco=2>
    //  ---> 6 instances

    // The as-is state is suited for ncols=1, i.e. for only one block-column per Green function
    // The advantages of more than 1 block-column were given for the block-sparse Hamiltonian,
    // however, not for the ultra-sparse Hamiltonian.
    // with ncols=1:
    // For the Green function data layout X[ncubes][R1C2][Noco*64][Noco*64]
    // SHOprj: p_{ai}[64_i] * X[icube][R1C2][Noco*64_i][Noco*64_j] summed over 64_i -> c[nai][R1C2][Noco][Noco*64_j]
    // SHOadd: Y[icube][R1C2][Noco*64_i][Noco*64_j] += p_{ai}[64_i] * d[nai][R1C2][Noco][Noco*64_j] summed over ai

    // Memory estimate: double complex Noco=1 --> 2^16 Byte =  64 KiByte per block
    //                  double complex Noco=2 --> 2^18 Byte = 256 KiByte per block
    // GPU capacity 40 GiByte
    // tfQMRgpu needs 10 instances --> 4 GiByte = 2^32 Byte per instance
    // max number of blocks is 2^16 = 65,536 (Noco=1) or 2^14 = 16,384 (Noco=2)
    //
    // typical coarse grid spacing: 1 block edge = 1 Angstrom --> grid spacing 0.25 Angstrom == 0.47 Bohr
    // sphere with volume 65ki has radius 25.011 (Noco=1)
    // sphere with volume 16ki has radius 15.756 (Noco=2)
    // Lattice constant of gold: 4.078  Angstrom (fcc) --> 4 atoms per (4.078  Angstrom)^3 --> ~17 Angstrom^3 per atom, 1088 DoF per atom
    // 2^16 blocks --> 3755 gold atom centers inside the truncation cluster
    // 2^14 blocks -->  964 gold atom centers inside the truncation cluster
    // Lattice constant of iron: 2.8665 Angstrom (bcc) --> 2 atoms per (2.8665 Angstrom)^3 --> 11.7                      753 DoF

    // ==> try first with ncols=1 implementation because it is simpler

    template <typename real_t, int R1C2=2, int Noco=1, int n64=64>
    void __global__ SHOsum( // launch SHOsum<real_t,R1C2,Noco,n64> <<< {natoms, nrhs, 1}, {Noco*n64, Noco, 1} >>>
#ifdef HAS_NO_CUDA
          dim3 const & gridDim, dim3 const & blockDim,
#endif // HAS_NO_CUDA
          double         (*const __restrict__ aac)[R1C2][Noco][Noco*n64] // if (collect) result, atom coefficients else input
        , real_t         (*const __restrict__ aic)[R1C2][Noco][Noco*n64] // if (collect) input,  atom image coeffs else result
        , int8_t   const (*const __restrict__ AtomLmax) // SHO basis size [iatom]
        , uint32_t const (*const __restrict__ AtomStarts) // prefix sum over nSHO(AtomLmax[:])
        , uint32_t const (*const __restrict__ AtomImageStarts) // prefix sum over nSHO(AtomImageLmax[:])
        , double   const (*const __restrict__ AtomImagePhase)[4] // complex/magnetic phase [nAtomImages]
        , uint32_t const (*const __restrict__ bsr_of_iatom) // row starts
        , uint32_t const (*const __restrict__ iai_of_bsr)   // col index
        , bool     const collect=true // otherwise input and result are exchanged
    )
      // Collect the atomic projection coefficients from the atom images
      // or broadcast the atomic addition coefficients to the atom images
    {
        assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));
        auto const nrhs  = gridDim.y;
        assert(1        ==  gridDim.z);
        assert(Noco*n64 == blockDim.x);
        assert(Noco     == blockDim.y);
        assert(1        == blockDim.z);

        float const scale_imaginary = collect ? 1 : -1; // simple complex conjugate

#ifndef HAS_NO_CUDA
        auto const irhs  = blockIdx.y;
        auto const iatom = blockIdx.x;
#else  // HAS_NO_CUDA
        auto const natoms = gridDim.x;
        for (int irhs = 0; irhs < nrhs; ++irhs)
        for (int iatom = 0; iatom < natoms; ++iatom)
#endif // HAS_NO_CUDA
        { // block loops

        int const a0 = AtomStarts[iatom];
        int const lmax = AtomLmax[iatom];
        int const nSHO = sho_tools::nSHO(lmax);

#ifndef HAS_NO_CUDA
        auto const spin = threadIdx.y;
        auto const ivec = threadIdx.x;
#else  // HAS_NO_CUDA
        for (int spin = 0; spin < blockDim.y; ++spin)
        for (int ivec = 0; ivec < blockDim.x; ++ivec)
#endif // HAS_NO_CUDA
        { // thread loops

        int constexpr Re = 0, Im = R1C2 - 1; // indices for real and imaginary part (if complex)
        if (collect) {
            // clear atom result coefficients before accumulating
            for (int ai = 0; ai < nSHO; ++ai) {
                aac[(a0 + ai)*nrhs + irhs][Re][spin][ivec] = 0.0;
                aac[(a0 + ai)*nrhs + irhs][Im][spin][ivec] = 0.0; // 0 flop
            } // ai
        } // collect

        for (auto bsr = bsr_of_iatom[iatom]; bsr < bsr_of_iatom[iatom + 1]; ++bsr) {
            // accumulate projection coefficients for each contributing atom

            auto const iai = iai_of_bsr[bsr]; // index of atom image

            auto const phase = AtomImagePhase[iai];
            auto const ph = std__complex<double>(phase[0], phase[1]*scale_imaginary);
            auto const i0 = AtomImageStarts[iai];
            // assert(AtomImageLmax[iai] == lmax);

            for (int ai = 0; ai < nSHO; ++ai) {
                if (collect) {
                    // collect: aac += phase * aic
                    if (1 == R1C2) {
                        aac[(a0 + ai)*nrhs + irhs][Re][0][ivec] += phase[0] *
                        aic[(i0 + ai)*nrhs + irhs][Re][0][ivec]; // 2 flop
                    } else {
                        auto const c = ph * std__complex<double>(
                        aic[(i0 + ai)*nrhs + irhs][Re][spin][ivec],
                        aic[(i0 + ai)*nrhs + irhs][Im][spin][ivec]);
                        aac[(a0 + ai)*nrhs + irhs][Re][spin][ivec] += c.real();
                        aac[(a0 + ai)*nrhs + irhs][Im][spin][ivec] += c.imag(); // 8 flop
                        // if this is slow, TODO use __shared__ memory for aac
                    } // 1 == R1C2
                } else {
                    // broadcast: aic = phase * aac
                    if (1 == R1C2) {
                        aic[(i0 + ai)*nrhs + irhs][Re][0][ivec] = phase[0] *
                        aac[(a0 + ai)*nrhs + irhs][Re][0][ivec]; // 2 flop
                    } else {
                        auto const c = ph * std__complex<double>(
                        aac[(a0 + ai)*nrhs + irhs][Re][spin][ivec],
                        aac[(a0 + ai)*nrhs + irhs][Im][spin][ivec]);
                        aic[(i0 + ai)*nrhs + irhs][Re][spin][ivec] = c.real();
                        aic[(i0 + ai)*nrhs + irhs][Im][spin][ivec] = c.imag(); // 8 flop
                    } // 1 == R1C2
                } // collect
               // only spin diagonal operator implemented, TODO extend using phase[2] and phase[3]
            } // ai

        } // bsr
        }} // thread loops and block loops

    } // SHOsum

    template <typename real_t, int R1C2=2, int Noco=1, int n64=64>
    void __host__ SHOsum_driver(
          double         (*const __restrict__ aac)[R1C2][Noco][Noco*n64] // result, atom coefficients
        , real_t         (*const __restrict__ aic)[R1C2][Noco][Noco*n64] // input,  atom image coeffs
        , int8_t   const (*const __restrict__ AtomLmax) // SHO basis size [iatom]
        , uint32_t const (*const __restrict__ AtomStarts) // prefix sum over nSHO(AtomLmax[:])
        , uint32_t const (*const __restrict__ AtomImageStarts) // prefix sum over nSHO(AtomImageLmax[:])
        , double   const (*const __restrict__ AtomImagePhase)[4] // complex/magnetic phase [nAtomImages][4]
        , green_sparse::sparse_t<> const & sparse_SHOsum
        , uint32_t const nAtoms // number of atoms
        , int      const nrhs
        , bool     const collect=true // otherwise input and result are exchanged
        , int      const echo=0
    )
      // Collect the atomic projection coefficients from the atom images
      // or broadcast the atomic addition coefficients to the atom images
    {
        assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));

        dim3 const gridDim(nAtoms, nrhs, 1), blockDim(Noco*n64, Noco, 1);
        if (echo > 3) std::printf("# %s<%s,R1C2=%d,Noco=%d> <<< {nAtoms=%d, nrhs=%d, 1}, {%d, Noco=%d, 1} >>>\n",
                            __func__, real_t_name<real_t>(), R1C2, Noco,  nAtoms, nrhs,  Noco*n64, Noco);
        SHOsum<real_t,R1C2,Noco,n64> // launch SHOsum<real_t,R1C2,Noco,n64> <<< {nAtoms, nrhs, 1}, {Noco*n64, Noco, 1} >>>
#ifndef HAS_NO_CUDA
            <<< gridDim, blockDim >>> (
#else  // HAS_NO_CUDA
              ( gridDim, blockDim,
#endif // HAS_NO_CUDA
                aac, aic, AtomLmax, AtomStarts, AtomImageStarts, AtomImagePhase,
                sparse_SHOsum.rowStart(), sparse_SHOsum.colIndex(), collect);
    } // SHOsum_driver





    template <typename real_t, int R1C2=2, int Noco=1, int n64=64>
    void __global__ SHOmul( // launch SHOmul<real_t,R1C2,Noco,n64> <<< {natoms, nrhs, 1}, {Noco*n64, Noco, 1} >>>
#ifdef HAS_NO_CUDA
          dim3 const & gridDim, dim3 const & blockDim,
#endif // HAS_NO_CUDA
          real_t         (*const __restrict__ aac)[R1C2][Noco][Noco*n64] // result,  atom addition coefficients
        , real_t   const (*const __restrict__ apc)[R1C2][Noco][Noco*n64] // input, atom projection coefficients
        , double   const (*const *const __restrict__ AtomMatrices) // matrices
        , int8_t   const (*const __restrict__ AtomLmax) // SHO basis size [iatom]
        , uint32_t const (*const __restrict__ AtomStarts) // prefix sum over nSHO(AtomLmax[:])
    )
      // Multiply the atom-specific matrices onto the projection coefficients
    {
        assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));
        int const nrhs   = gridDim.y;
        assert(1       ==  gridDim.z);
        assert(Noco*n64== blockDim.x);
        assert(Noco    == blockDim.y);
        assert(1       == blockDim.z);

#ifndef HAS_NO_CUDA
        int const irhs  = blockIdx.y;
        int const iatom = blockIdx.x;
#else  // HAS_NO_CUDA
        int const natoms = gridDim.x;
        for (int irhs = 0; irhs < nrhs; ++irhs)
        for (int iatom = 0; iatom < natoms; ++iatom)
#endif // HAS_NO_CUDA
        { // block loops

        int const a0 = AtomStarts[iatom];
        int const lmax = AtomLmax[iatom];
        int const nSHO = sho_tools::nSHO(lmax);
        auto const AtomMat = AtomMatrices[iatom];

#ifndef HAS_NO_CUDA
        int const spin = threadIdx.y;
        int const ivec = threadIdx.x;
#else  // HAS_NO_CUDA
        for (int spin = 0; spin < blockDim.y; ++spin)
        for (int ivec = 0; ivec < blockDim.x; ++ivec)
#endif // HAS_NO_CUDA
        { // thread loops

            for (int ai = 0; ai < nSHO; ++ai) {
                int constexpr Real = 0;
                if (1 == R1C2) { // is real
                    // version for R1C2==1,Noco==1 is more readable
                    double cad{0};
                    for (int aj = 0; aj < nSHO; ++aj) {
                        double const cpr = apc[(a0 + aj)*nrhs + irhs][Real][0][ivec]; // load projection coefficient
                        double const am = AtomMat[ai*nSHO + aj];                      // load matrix element
                        cad += am * cpr; // 2 flop
                    } // aj
                    aac[(a0 + ai)*nrhs + irhs][Real][0][ivec] = cad;                  // store addition coefficient

                } else {
                    // matrix layout: AtomMat[Noco*Noco*R1C2*nSHO*nSHO]
                    assert(2 == R1C2); // is complex
                    int constexpr Imag = R1C2 - 1; // index for imaginary part
                    std__complex<double> cad(0.0, 0.0);
                    for (int spjn = 0; spjn < Noco; ++spjn) {
                        for (int aj = 0; aj < nSHO; ++aj) {
                            // load projection coefficient
                            std__complex<double> const cpr(     // load projection coefficient
                                       apc[(a0 + aj)*nrhs + irhs][Real][spjn][ivec],
                                       apc[(a0 + aj)*nrhs + irhs][Imag][spjn][ivec]);
                            // load matrix element
                            std__complex<double> const am(      // load matrix element
                                       AtomMat[(((spin*Noco + spjn)*R1C2 + Real)*nSHO + ai)*nSHO + aj],
                                       AtomMat[(((spin*Noco + spjn)*R1C2 + Imag)*nSHO + ai)*nSHO + aj]);
                            cad += am * cpr; // 8 flop
                            // if (0 == ivec && 0 == irhs) {
                            //     std::printf("# %s atom=%i cpr=(%g, %g),\tam[%2i,%2i]=(%g, %g)\n",
                            //               __func__, iatom, cpr.real(), cpr.imag(), ai, aj, am.real(), am.imag());
                            // }
                        } // aj
                    } // spjn
                    aac[(a0 + ai)*nrhs + irhs][Real][spin][ivec] = cad.real(); // store addition coefficient
                    aac[(a0 + ai)*nrhs + irhs][Imag][spin][ivec] = cad.imag(); // store addition coefficient

                } // 1 == R1C2
            } // ai

        }} // thread loops and block loops

    } // SHOmul

    template <typename real_t, int R1C2=2, int Noco=1, int n64=64>
    void __host__ SHOmul_driver(
          real_t         (*const __restrict__ aac)[R1C2][Noco][Noco*n64] // result,  atom addition coefficients
        , real_t   const (*const __restrict__ apc)[R1C2][Noco][Noco*n64] // input, atom projection coefficients
        , double   const (*const *const __restrict__ AtomMatrices)
        , int8_t   const (*const __restrict__ AtomLmax)
        , uint32_t const (*const __restrict__ AtomStarts)
        , int      const nAtoms // number of atoms
        , int      const nrhs // number of block columns in the Green function
        , int const echo=0 // log level
    ) {
        assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));
        if (nAtoms*nrhs < 1) return;

        dim3 const gridDim(nAtoms, nrhs, 1), blockDim(Noco*n64, Noco, 1);
        if (echo > 3) std::printf("# %s<%s,R1C2=%d,Noco=%d> <<< {nAtoms=%d, nrhs=%d, 1}, {%d, Noco=%d, 1} >>>\n",
                           __func__, real_t_name<real_t>(), R1C2, Noco,  nAtoms, nrhs,  Noco*n64, Noco);
        SHOmul<real_t,R1C2,Noco,n64> // launch <<< {nAtoms, nrhs, 1}, {Noco*n64, Noco, 1} >>>
#ifndef HAS_NO_CUDA
            <<< gridDim, blockDim >>> (
#else //  HAS_NO_CUDA
           (    gridDim, blockDim,
#endif // HAS_NO_CUDA
            aac, apc, AtomMatrices, AtomLmax, AtomStarts);
    } // SHOmul_driver



    template <typename real_t, int R1C2=2, int Noco=1>
    size_t __host__ multiply(
          real_t         (*const __restrict__ Ppsi)[R1C2][Noco*64][Noco*64] // result,  modified Green function blocks [nnzb][R1C2][Noco*64][Noco*64]
        , real_t         (*const __restrict__  Cpr)[R1C2][Noco]   [Noco*64] // projection coefficients     [natomcoeffs*nrhs][R1C2][Noco   ][Noco*64]
        , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input, unmodified Green function blocks [nnzb][R1C2][Noco*64][Noco*64]
        , dyadic_plan_t const & p
        , uint32_t const (*const __restrict__ RowIndexCubes) // Green functions rowIndex[nnzb]
        , uint16_t const (*const __restrict__ ColIndexCubes) // Green functions colIndex[nnzb]
        , float    const (*const __restrict__ CubePos)[3+1] // cube positions +alignment
        , uint32_t const nnzb // number of blocks of the Green function
        , int const echo=0 // log-level
        , double (*const __restrict__  Cpr_export)[R1C2][Noco][Noco*64]=nullptr // optional result, reduced projection coefficients [ncoeffs*nrhs]
    ) {
        assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));
        if (p.nAtomImages*p.nrhs < 1) return 0; // empyt GPU kernels may not run

        assert(p.AtomImageStarts);
        size_t const natomcoeffs = p.AtomImageStarts[p.nAtomImages];
        if (echo > 6) std::printf("# %s<%s,R1C2=%d,Noco=%d> nAtoms=%d nAtomImages=%d nrhs=%d ncoeffs=%ld\n",
                  __func__, real_t_name<real_t>(), R1C2, Noco, p.nAtoms, p.nAtomImages, p.nrhs, natomcoeffs);

        SHOprj_driver<real_t,R1C2,Noco>(Cpr, psi, p.AtomImagePos, p.AtomImageLmax, p.AtomImageStarts, p.nAtomImages,
                                                p.sparse_SHOprj, RowIndexCubes, CubePos, p.grid_spacing, p.nrhs, echo);

        real_t (*Cad)[R1C2][Noco][Noco*64]{nullptr};

        size_t const ncoeffs = p.AtomStarts[p.nAtoms];
        if (p.nAtomImages > p.nAtoms) {
            // we have to distinguish between atoms and their periodic images, this case: more images than atoms

            auto cprj = Cpr_export ? Cpr_export : get_memory<double[R1C2][Noco][Noco*64]>(ncoeffs*p.nrhs, echo, "cprj");

            SHOsum_driver<real_t,R1C2,Noco>(cprj, Cpr, p.AtomLmax, p.AtomStarts, p.AtomImageStarts, p.AtomImagePhase,
                                            p.sparse_SHOsum, p.nAtoms, p.nrhs, true, echo); // true:collect
            if (Cpr_export) { return 0; }

            auto cadd = get_memory<double[R1C2][Noco][Noco*64]>(ncoeffs*p.nrhs, echo, "cadd");

            SHOmul_driver<double,R1C2,Noco>(cadd, cprj, p.AtomMatrices, p.AtomLmax, p.AtomStarts, p.nAtoms, p.nrhs, echo);

            free_memory(cprj);

            Cad = get_memory<real_t[R1C2][Noco][Noco*64]>(natomcoeffs*p.nrhs, echo, "Cad"); // rectangular storage of atoms x right-hand-sides

            SHOsum_driver<real_t,R1C2,Noco>(cadd, Cad, p.AtomLmax, p.AtomStarts, p.AtomImageStarts, p.AtomImagePhase,
                                            p.sparse_SHOsum, p.nAtoms, p.nrhs, false, echo); // false:broadcast

            free_memory(cadd);

        } else {
            assert(p.nAtomImages == p.nAtoms); // there is at most one relevant image of each atom so we can ignore Bloch phases
            assert(natomcoeffs == ncoeffs);
            if (Cpr_export) { set(Cpr_export[0][0][0], ncoeffs*p.nrhs*R1C2*Noco*Noco*64, Cpr[0][0][0]); return 0; }

            Cad = get_memory<real_t[R1C2][Noco][Noco*64]>(natomcoeffs*p.nrhs, echo, "Cad"); // rectangular storage of atoms x right-hand-sides

            SHOmul_driver<real_t,R1C2,Noco>(Cad, Cpr, p.AtomMatrices, p.AtomLmax, p.AtomStarts, p.nAtoms, p.nrhs, echo);

        } // more images than atoms

        SHOadd_driver<real_t,R1C2,Noco>(Ppsi, Cad, p.AtomImagePos, p.AtomImageLmax, p.AtomImageStarts,
                                        p.sparse_SHOadd.rowStart(), p.sparse_SHOadd.colIndex(),
                                        RowIndexCubes, ColIndexCubes, CubePos, p.grid_spacing, nnzb, p.nrhs, echo);
        free_memory(Cad);

        return p.get_flop_count(R1C2, Noco, echo);
    } // multiply (dyadic operations)


    inline std::vector<double> __host__ sho_normalization(int const lmax, double const sigma=1) {
        int const n1ho = sho_tools::n1HO(lmax);
        std::vector<double> v1(n1ho, constants::sqrtpi*sigma);
        {
            double fac{1};
            for (int nu = 0; nu <= lmax; ++nu) {
                // fac == factorial(nu) / 2^nu
                v1[nu] *= fac;
                // now v1[nu] == sqrt(pi) * sigma * factorial(nu) / 2^nu
                fac *= 0.5*(nu + 1); // update fac for the next iteration
            } // nu
        }

        std::vector<double> vec(sho_tools::nSHO(lmax), 1.);
        {
            int sho{0};
            for (int iz = 0; iz <= lmax; ++iz) {
                for (int iy = 0; iy <= lmax - iz; ++iy) {
                    for (int ix = 0; ix <= lmax - iz - iy; ++ix) {
                        vec[sho] = v1[ix] * v1[iy] * v1[iz];
                        ++sho;
            }}} // ix iy iz
            assert(vec.size() == sho);
        }
        return vec;
    } // sho_normalization


    template <int R1C2=2, int Noco=1>
    std::vector<std::vector<double>> __host__ SHOprj_right( // result: density matrices
          double   const (*const __restrict__ pGreen)[R1C2][Noco][Noco*64] // input: projected Green function <p|G, layout[natomcoeffs*nrhs][R1C2][Noco][Noco*64]
        , double   const (*const __restrict__ AtomImagePos)[3+1] // atomic positions [0],[1],[2], decay parameter [3]
        , int8_t   const (*const __restrict__ AtomImageLmax) // SHO basis size [nAtomImages]
        , uint32_t const (*const __restrict__ AtomImageStarts) // prefix sum over nSHO(AtomImageLmax[:])
        , uint32_t const (*const __restrict__ AtomImageIndex) // [nAtomImages]
        , double   const (*const __restrict__ AtomImagePhase)[4] // [nAtomImages][4]
        , uint32_t const nAtomImages
        , int8_t   const (*const __restrict__ AtomLmax) // SHO basis size [nAtoms]
        , uint32_t const nAtoms
        , float    const (*const __restrict__ colCubePos)[3+1] // only [0],[1],[2] used, CubePos[irhs][0:3] of RHS cubes
        , double   const (*const __restrict__ hGrid) // grid spacings in [0],[1],[2], projection radius in [3]
        , uint32_t const nrhs // number or right hand sides
        , int      const echo=0 // log-level
        , double   const factor=1.0
    )
        // Other than SHOprj, this version of Spherical Harmonic Oscillator projects onto the Green functions right index
    {
        // SimpleTimer timer(__FILE__, __LINE__, __func__, echo);
        assert(1 == Noco || 2 == Noco);

        if (echo > 0) std::printf("# %s for %d atoms, %d atom images\n", __func__, nAtoms, nAtomImages);
        std::vector<std::vector<double>> pGp(nAtoms); // result: projection coefficients, layout[nAtoms][nSHO*nSHO*Noco*Noco]
        for (uint32_t iatom = 0; iatom < nAtoms; ++iatom) {
            int const lmax = AtomLmax[iatom];
            int const nSHO = sho_tools::nSHO(lmax);
            pGp[iatom].resize(nSHO*nSHO*Noco*Noco, 0.0); // init
        } // iai

        typedef std::complex<double> complex_t;
        complex_t constexpr zero = complex_t(0, 0);

        for (uint32_t iai = 0; iai < nAtomImages; ++iai) {

        int const lmax = AtomImageLmax[iai];
        int constexpr Lmax = 7;
        if (lmax > Lmax) std::printf("# %s Error: lmax= %d but max. Lmax= %d, iai=%d, iatom=%d\n", __func__, lmax, Lmax, iai, AtomImageIndex[iai]);
        assert(lmax <= Lmax);

        auto const a0 = AtomImageStarts[iai];
        int const nSHO = sho_tools::nSHO(lmax);
        int const n2HO = sho_tools::n2HO(lmax);
        int const n1HO = sho_tools::n1HO(lmax);

        std::vector<complex_t> piGp(nSHO*nSHO*Noco*Noco, zero);

        double H1D[1 + Lmax][3][4]; // non-normalized orthogonal 1-dimensional Hermite Gauss functions
        float     xi_squared[3][4]; // distance along one Cartesian direction squared

        // ToDo: would it be good to move these 3 loops inside?
        for (int spin = 0; spin < Noco; ++spin) {
        for (int spjn = 0; spjn < Noco; ++spjn) {
        for (int isho = 0; isho < nSHO; ++isho) {

        std::vector<complex_t> czyx(nSHO, zero);

        // in this loop over RHS we accumulate atom projection coefficients over many cubes, however, we also need MPI_SUM later
        for (uint32_t irhs = 0; irhs < nrhs; ++irhs) {

            float R2_projection;
            for (int i12 = 0; i12 < 3*4; ++i12) { // 3 directions times 4 grid points in a block edge
                R2_projection = Hermite_polynomials_1D(H1D, xi_squared, i12, lmax, AtomImagePos[iai], colCubePos[irhs], hGrid);
            } // i12

            if (echo > 16) std::printf("# %s  iai=%d, iatom=%d, lmax=%d, irhs=%d, R^2=%g, H0= %g %g %g\n",
                __func__, iai, AtomImageIndex[iai], lmax, irhs, R2_projection, H1D[0][X][0], H1D[0][Y][0], H1D[0][Z][0]);

            for (int z4 = 0; z4 < 4; ++z4) { // loop over real-space grid in z-direction
                auto const d2z = xi_squared[Z][z4] - R2_projection;
                if (d2z < 0) {
                    std::vector<complex_t> byx(n2HO, zero);
                    for (int y4 = 0; y4 < 4; ++y4) { // loop over real-space grid in y-direction
                        auto const d2yz = xi_squared[Y][y4] + d2z;
                        if (d2yz < 0) {
                            std::vector<complex_t> ax(n1HO, zero);
                            for (int x4 = 0; x4 < 4; ++x4) { // loop over real-space grid in x-direction
                                auto const d2xyz = xi_squared[X][x4] + d2yz;
                                if (d2xyz < 0) {
                                    int const xyz = (z4*4 + y4)*4 + x4;
                                    int constexpr RealPart = 0, ImagPart= R1C2 - 1;
                                    auto const pG = complex_t(pGreen[(a0 + isho)*nrhs + irhs][RealPart][spin][spjn*64 + xyz],  // load Re{<p|G}
                                                              pGreen[(a0 + isho)*nrhs + irhs][ImagPart][spin][spjn*64 + xyz]); // load Im{<p|G}
                                    for (int ix = 0; ix <= lmax; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                        ax[ix] += pG * H1D[ix][X][x4];
                                    } // ix
                                } // inside mask x
                            } // x4
                            for (int iy = 0; iy <= lmax; ++iy) { // loop over the 2nd Cartesian SHO quantum number
                                int const iyx0 = (iy*(2*lmax + 3 - iy)) >> 1;
                                auto const Hy = H1D[iy][Y][y4]; // load Hy
                                for (int ix = 0; ix <= lmax - iy; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                    byx[iyx0 + ix] += ax[ix] * Hy;
                                } // ix
                            } // iy
                        } // inside mask y
                    } // y4
                    for (int jsho = 0, iz = 0; iz <= lmax; ++iz) { // loop over the 3rd Cartesian SHO quantum number
                        auto const Hz = H1D[iz][Z][z4]; // load Hz
                        for (int iy = 0; iy <= lmax - iz; ++iy) { // loop over the 2nd Cartesian SHO quantum number
                            int const iyx0 = (iy*(2*lmax + 3 - iy)) >> 1;
                            for (int ix = 0; ix <= lmax - iz - iy; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                czyx[jsho] += byx[iyx0 + ix] * Hz; // form <p|G|p>
                                ++jsho;
                            } // ix
                        } // iy
                    } // iz
                } // inside mask z
            } // z4

        } // irhs

        for (int jsho = 0; jsho < nSHO; ++jsho) {
            piGp[((isho*nSHO + jsho)*Noco + spin)*Noco + spjn] = czyx[jsho]; // TODO: inefficient memory access to piGp, maybe move {isho, spin, spjn} loops to inside
        } // jsho

        } // isho
        } // spjn
        } // spin

            auto const iatom = AtomImageIndex[iai];
            assert(lmax                == AtomLmax[iatom]   && "Inconsistency between AtomImageLmax[:] and AtomLmax[AtomImageIndex[:]]");
            assert(nSHO*nSHO*Noco*Noco == pGp[iatom].size() && "Inconsistency between AtomImageLmax[:] and AtomLmax[AtomImageIndex[:]]");
            auto const phase = AtomImagePhase[iai];
            auto const ph = complex_t(phase[0], phase[1]); // Block phase factor, ToDo: needs a conjugation here?, how to deal with spin?
            for (int ij = 0; ij < nSHO*nSHO*Noco*Noco; ++ij) {
                complex_t const c = piGp[ij] * ph; // complex multiplication
                pGp[iatom][ij] += c.imag()*factor; // add only the imaginary part, accumulate over images of the same atom
            } // ij

        } // iai

        return pGp;
    } // SHOprj_right



    template <typename real_t, int R1C2=2, int Noco=1>
    std::vector<std::vector<double>> __host__ get_projection_coefficients( // result: atomic density matrices
          real_t   const (*const __restrict__ Green)[R1C2][Noco*64][Noco*64] // input, unmodified Green function blocks [nnzb][R1C2][Noco*64][Noco*64]
        , dyadic_plan_t const & p
        , uint32_t const (*const __restrict__ RowIndexCubes) // Green functions rowIndex[nnzb]
        , float    const (*const __restrict__ rowCubePos)[3+1] // cube positions +alignment of row indices
        , float    const (*const __restrict__ colCubePos)[3+1] // cube positions +alignment of col indices
        , int const echo=0 // log-level
    )
        // The typical projection operation of green_dyadic is SHOprj performing <p|G, however,
        // SHOprj_right contracts the localized Gauss-Hermite functions with the right Green function index
        // which is necessary to find atomic density matrices as <p_i|G|p_j>.
    {
        if (echo > 0) std::printf("# %s\n", __func__);
        size_t const ncoeffs = p.AtomStarts[p.nAtoms];
        auto Cpr_export = get_memory<double[R1C2][Noco][Noco*64]>(ncoeffs*p.nrhs, echo, "Cpr_export");
        size_t const natomcoeffs = p.AtomImageStarts[p.nAtomImages];
        auto Cpr = get_memory<real_t[R1C2][Noco][Noco*64]>(natomcoeffs*p.nrhs, echo, "Cpr");
        auto const nflop = multiply<real_t,R1C2,Noco>(nullptr, Cpr, Green, p, RowIndexCubes, nullptr, rowCubePos, 0, echo, Cpr_export);
        free_memory(Cpr);
        assert(0 == nflop);
        if (echo > 3) std::printf("# %s from projected Green function\n", __func__);
        assert(p.nAtomImages == p.nAtoms && "reduction over Bloch phases not yet implemented!");
        auto const result = SHOprj_right(Cpr_export, p.AtomImagePos, p.AtomImageLmax, p.AtomImageStarts, p.AtomImageIndex,
                                p.AtomImagePhase, p.nAtomImages, p.AtomLmax, p.nAtoms, colCubePos, p.grid_spacing, p.nrhs, echo);
        free_memory(Cpr_export);
        warn("missing is an MPI_Allreduce onto each atom matrix for %d local atoms", int(p.nAtoms));
        return result;
    } // get_projection_coefficients

    status_t all_tests(int const echo=0); // declaration only

} // namespace green_dyadic
