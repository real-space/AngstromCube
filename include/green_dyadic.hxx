#pragma once

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_memory.hxx" // get_memory, free_memory, dim3, real_t_name
#include "green_sparse.hxx" // ::sparse_t<>
#include "inline_math.hxx" // pow2, pow3
#include "sho_tools.hxx" // ::nSHO, ::n2HO, ::n1HO

#ifndef HAS_NO_CUDA
    #include <cuda/std/complex> // std::complex
    #define std__complex cuda::std::complex
#else
    #include <complex> // std::complex
    #define std__complex std::complex
#endif // HAS_NO_CUDA

namespace green_dyadic {


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

        double const R2_projection = pow2(double(hxyz[3])); // projection radius can be controlled from outside

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

            double const xi = sigma_inverse*(grid_spacing*(cube_position*n4 + i4 + 0.5) - atom_position); // center of the grid cell
            double const xi2 = xi*xi;
            xi_squared[idir][i4] = xi2; // export the square of the distance in units of sigma
            double const H0 = (xi2 < R2_projection) ? std::exp(-0.5*xi2) : 0; // Gaussian envelope function
            
            // Discuss:   We could export xi*xi as an array double (or maybe only float) xi2[3][4] and compare
            //              xi2[0][x] + xi2[1][y] + xi2[2][z] < R2_projection
            //            hence omitting the data transfer of BitMask entirely and save preparation code complexity


            H1D[0][idir][i4] = H0; // store Gaussian H_0

            double Hnp1, Hn{H0}, Hnm1{0};
//          std::printf("# %s idir=%d xi= %g H= %g", __func__, idir, xi, H0);
            for (int nu = 0; nu < lmax; ++nu) {
                double const nuhalf = 0.5 * nu; // nu/2.
                // use the two step recurrence relation to create the 1D Hermite polynomials
                // H1D[nu+1][idir][i4] = xi * H1D[nu][idir][i4] - nu/2. * H1D[nu-1][idir][i4];
                Hnp1 = xi * Hn - nuhalf * Hnm1; // see snippets/3Hermite.F90 for the derivation of this
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

#ifndef HAS_NO_CUDA
            if (threadIdx.x < 3) xyzc[threadIdx.x] = CubePos[irow][threadIdx.x]; // load into shared memory
            __syncthreads();

            // generate Hx, Hy, Hz up to lmax inside the 4^3 cube
            auto const R2_proj = Hermite_polynomials_1D(H1D, xi_squared, j, lmax, xyza, xyzc, hgrid); // sufficient to be executed by (threadIdx.y == 0, threadIdx.z == 0)

            __syncthreads();
#else  // HAS_NO_CUDA
            set(xyzc, 3, CubePos[irow]);

            float R2_proj;
            for (int jj = 0; jj < 12; ++jj) {
                R2_proj = Hermite_polynomials_1D(H1D, xi_squared, jj, lmax, xyza, xyzc, hgrid);
            } // jj
#endif // HAS_NO_CUDA


            for (int z = 0; z < 4; ++z) { // loop over real-space grid in z-direction
                auto const d2z = xi_squared[2][z] - R2_proj;
                if (d2z < 0) {

                    real_t byx[sho_tools::n2HO(Lmax)];
                    for (int iyx = 0; iyx < sho_tools::n2HO(lmax); ++iyx) { // loop over the combined index of ix and iy
                        byx[iyx] = 0; // init
                    } // iyx

                    for (int y = 0; y < 4; ++y) { // loop over real-space grid in y-direction
                        auto const d2yz = xi_squared[1][y] + d2z;
                        if (d2yz < 0) {

                            real_t ax[Lmax + 1];
                            for (int ix = 0; ix <= lmax; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                ax[ix] = 0; // init
                            } // ix

                            // __unroll__
                            for (int x = 0; x < 4; ++x) { // loop over real-space grid in x-direction
                                auto const d2xyz = xi_squared[0][x] + d2yz;
                                if (d2xyz < 0) {
                                    int const xyz = (z*4 + y)*4 + x;
                                    real_t const ps = Psi[inzb][reim][spin*64 + xyz][j]; // load from global memory
                                    for (int ix = 0; ix <= lmax; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                        auto const Hx = H1D[ix][0][x]; // load Hx from shared memory
                                        ax[ix] += ps * Hx; // FMA: 2 flop * 4**3 * (L+1)
                                    } // ix
                                } // inside mask
                            } // x

                            for (int iy = 0; iy <= lmax; ++iy) { // loop over the 2nd Cartesian SHO quantum number
                                int const iyx0 = (iy*(2*lmax + 3 - iy)) >> 1;
                                for (int ix = 0; ix <= lmax - iy; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                    auto const Hy = H1D[iy][1][y]; // load Hy from shared memory
                                    byx[iyx0 + ix] += ax[ix] * Hy; // FMA: 2 flop * 4**2 * ((L+1)*(L+2))/2
                                } // ix
                            } // iy

                        } // inside mask
                    } // y

                    for (int sho = 0, iz = 0; iz <= lmax; ++iz) { // loop over the 3rd Cartesian SHO quantum number
                        auto const Hz = H1D[iz][2][z]; // load Hz from shared memory
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
            // FMA*( 4**3 * ((L+1)) + 4**2 * ((L+1)*(L+2))/2 + 4**1 * ((L+1)*(L+2)*(L+3))/6 )
            // so for Lmax=5 this is 944 FMAs = 1888 flop

        } // bsr

        { // scope: store accumulators
            auto const a0 = AtomStarts[iatom];
            // ToDo: use compressed storage to allow for lmax[iatom] without waisting memory for Cpr
            for (int sho = 0; sho < sho_tools::nSHO(lmax); ++sho) {
                Cpr[(a0 + sho)*nrhs + irhs][reim][spin][j] = czyx[sho];
            } // sho
        } // scope

        }} // thread loops and block loops
        
    } // SHOprj

    template <typename real_t, int R1C2=2, int Noco=1>
    size_t __host__ SHOprj_driver(
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

        dim3 const gridDim(natoms, nrhs, 1), blockDim(Noco*64, Noco, R1C2);
        if (echo > 3) std::printf("# %s<%s,R1C2=%d,Noco=%d> <<< {natoms=%d, nrhs=%d, 1}, {%d, %d, %d} >>>\n",
                            __func__, real_t_name<real_t>(), R1C2, Noco, natoms, nrhs, Noco*64, Noco, R1C2);
        SHOprj<real_t,R1C2,Noco> // launch <<< {natoms, nrhs, 1}, {Noco*64, Noco, R1C2} >>>
#ifndef HAS_NO_CUDA
              <<< gridDim, blockDim >>> (
#else  // HAS_NO_CUDA
                ( gridDim, blockDim,
#endif // HAS_NO_CUDA
               Cpr, Psi, sparse, AtomPos, AtomLmax, AtomStarts, RowIndexCube, CubePos, hGrid);

        return 0; // TODO estimate the total number of floating point operations
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
            if (lmax > Lmax) std::printf("# %s Error: lmax=%d but max. Lmax=%d, iatom=%d, bsr=%d, inzb=%d\n", __func__, lmax, Lmax, iatom, bsr, inzb);
            assert(lmax <= Lmax);

#ifndef HAS_NO_CUDA
            if (threadIdx.x < 4) xyza[threadIdx.x] = AtomPos[iatom][threadIdx.x]; // coalesced load

            __syncthreads();

            // generate Hx, Hy, Hz up to Lmax inside the 4^3 cube
            auto const R2_proj = Hermite_polynomials_1D(H1D, xi_squared, j, lmax, xyza, xyzc, hgrid); // sufficient to be executed by (threadIdx.y == 0, threadIdx.z == 0)

#else  // HAS_NO_CUDA
            set(xyza, 4, AtomPos[iatom]);

            float R2_proj;
            for (int jj = 0; jj < 12; ++jj) {
                R2_proj = Hermite_polynomials_1D(H1D, xi_squared, jj, lmax, xyza, xyzc, hgrid);
            } // jj
#endif // HAS_NO_CUDA
              
            int sho{0};
            for (int iz = 0; iz <= lmax; ++iz) { // executes (Lmax + 1) times

                real_t byx[4][4]; // get 16 accumulator registers --> 32 registers if double
                // __unroll__
                for (int y = 0; y < 4; ++y) {
                    // __unroll__
                    for (int x = 0; x < 4; ++x) {
                        byx[y][x] = 0; // init
                    } // x
                } // y

                for (int iy = 0; iy <= lmax - iz; ++iy) { // executes (Lmax + 1)*(Lmax + 2)/2 times

                    real_t ax[4] = {0, 0, 0, 0}; // get 4 accumulator registers --> 8 32bit GPU registers if double

                    for (int ix = 0; ix <= lmax - iz - iy; ++ix) { // executes (Lmax + 1)*(Lmax + 2)*(Lmax + 3)/6 times
                        real_t const ca = Cad[(a0 + sho)*nrhs + irhs][reim][spin][j]; // load Noco*64 numbers from global memory
                        ++sho;
                        // __unroll__
                        for (int x = 0; x < 4; ++x) {
                            real_t const Hx = H1D[ix][0][x]; // load Hx from shared memory
                            ax[x] += Hx * ca; // FMA: 2 flop * 4 * (Lmax + 1)*(Lmax + 2)*(Lmax + 3)/6
                        } // x
                    } // ix

                    // __unroll__
                    for (int y = 0; y < 4; ++y) {
                        real_t const Hy = H1D[iy][1][y]; // load Hy from shared memory
                        // __unroll__
                        for (int x = 0; x < 4; ++x) {
                            byx[y][x] += Hy * ax[x]; // FMA: 2 flop * 4**2 * (Lmax + 1)*(Lmax + 2)/2
                        } // x
                    } // y

                } // iy

                int64_t mask_atom{0}; // occupy 2 GPU registers
                // __unroll__
                for (int z = 0; z < 4; ++z) {
                    real_t const Hz = H1D[iz][2][z]; // load Hz from shared memory
                    auto const d2z = xi_squared[2][z] - R2_proj;
                    if (d2z < 0) {
                        // __unroll__
                        for (int y = 0; y < 4; ++y) {
                            auto const d2yz = xi_squared[1][y] + d2z;
                            if (d2yz < 0) {
                                // __unroll__
                                for (int x = 0; x < 4; ++x) {
                                    auto const d2xyz = xi_squared[0][x] + d2yz;
                                    if (d2xyz < 0) {
                                        czyx[z][y][x] += Hz * byx[y][x]; // FMA: 2 flop * 4**3 * (Lmax + 1) are these real flops?
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
    size_t __host__ SHOadd_driver(
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

        dim3 const gridDim(nnzb, 1, 1), blockDim(Noco*64, Noco, R1C2);
        if (echo > 3) std::printf("# %s<%s,R1C2=%d,Noco=%d> <<< {nnzb=%d, 1, 1}, {%d, %d, %d} >>>\n",
                            __func__, real_t_name<real_t>(), R1C2, Noco, nnzb, Noco*64, Noco, R1C2);
        SHOadd<real_t,R1C2,Noco> // launch {nnzb, 1, 1}, {Noco*64, Noco, R1C2} >>>
#ifndef HAS_NO_CUDA
              <<< gridDim, blockDim >>> (
#else  // HAS_NO_CUDA
                ( gridDim, blockDim,
#endif // HAS_NO_CUDA
               Psi, Cad, RowStartCubes, ColIndexAtoms, RowIndexCubes, ColIndexCubes, AtomPos, AtomLmax, AtomStarts, CubePos, hGrid, nrhs);

        return 0; // TODO estimate the total number of floating point operations
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
          double         (*const __restrict__ aac)[R1C2][Noco][Noco*n64] // result, atom coefficients
        , real_t         (*const __restrict__ aic)[R1C2][Noco][Noco*n64] // input,  atom image coeffs
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
        int const nrhs   = gridDim.y;
        assert(1       ==  gridDim.z);
        assert(Noco*n64== blockDim.x);
        assert(Noco    == blockDim.y);
        assert(1       == blockDim.z);

        float const scale_imaginary = collect ? 1 : -1; // simple complex conjugate
        
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

#ifndef HAS_NO_CUDA
        int const spin = threadIdx.y;
        int const ivec = threadIdx.x;
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
                aac[(a0 + ai)*nrhs + irhs][Im][spin][ivec] = 0.0;
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
                        aic[(i0 + ai)*nrhs + irhs][Re][0][ivec];
                    } else {
                        auto const c = ph * std__complex<double>(
                        aic[(i0 + ai)*nrhs + irhs][Re][spin][ivec],
                        aic[(i0 + ai)*nrhs + irhs][Im][spin][ivec]);
                        aac[(a0 + ai)*nrhs + irhs][Re][spin][ivec] += c.real();
                        aac[(a0 + ai)*nrhs + irhs][Im][spin][ivec] += c.imag();
                        // if this is slow, TODO use __shared__ memory for aac
                    } // 1 == R1C2
                } else {
                    // broadcast: aic = phase * aac
                    if (1 == R1C2) {
                        aic[(i0 + ai)*nrhs + irhs][Re][0][ivec] = phase[0] *
                        aac[(a0 + ai)*nrhs + irhs][Re][0][ivec];
                    } else {
                        auto const c = ph * std__complex<double>(
                        aac[(a0 + ai)*nrhs + irhs][Re][spin][ivec],
                        aac[(a0 + ai)*nrhs + irhs][Im][spin][ivec]);
                        aic[(i0 + ai)*nrhs + irhs][Re][spin][ivec] = c.real();
                        aic[(i0 + ai)*nrhs + irhs][Im][spin][ivec] = c.imag();
                    } // 1 == R1C2
                } // collect
               // only spin diagonal operator implemented, TODO extend using phase[2] and phase[3]
            } // ai

        } // bsr
        }} // thread loops and block loops

    } // SHOsum

    template <typename real_t, int R1C2=2, int Noco=1, int n64=64>
    size_t __host__ SHOsum_driver(
          double         (*const __restrict__ aac)[R1C2][Noco][Noco*n64] // result, atom coefficients
        , real_t         (*const __restrict__ aic)[R1C2][Noco][Noco*n64] // input,  atom image coeffs
        , int8_t   const (*const __restrict__ AtomLmax) // SHO basis size [iatom]
        , uint32_t const (*const __restrict__ AtomStarts) // prefix sum over nSHO(AtomLmax[:])
        , uint32_t const (*const __restrict__ AtomImageStarts) // prefix sum over nSHO(AtomImageLmax[:])
        , double   const (*const __restrict__ AtomImagePhase)[4] // complex/magnetic phase [nAtomImages]
        , green_sparse::sparse_t<> const & sparse_SHOsum
        , uint32_t const nAtoms
        , int      const nrhs
        , bool     const collect=true // otherwise input and result are exchanged
        , int      const echo=0
    )
      // Collect the atomic projection coefficients from the atom images
      // or broadcast the atomic addition coefficients to the atom images
    {
        assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));

        dim3 const gridDim(nAtoms, nrhs, 1), blockDim(Noco*n64, Noco, 1);
        if (echo > 3) std::printf("# %s<%s,R1C2=%d,Noco=%d,%d> <<< {nAtoms=%d, nrhs=%d, 1}, {%d, %d, 1} >>>\n",
                            __func__, real_t_name<real_t>(), R1C2, Noco, n64,  nAtoms, nrhs,  Noco*n64, Noco);
        SHOsum<real_t,R1C2,Noco,n64> // launch SHOsum<real_t,R1C2,Noco,n64> <<< {nAtoms, nrhs, 1}, {Noco*n64, Noco, 1} >>>
#ifndef HAS_NO_CUDA
            <<< gridDim, blockDim >>> (
#else  // HAS_NO_CUDA
              ( gridDim, blockDim,
#endif // HAS_NO_CUDA
                aac, aic, AtomLmax, AtomStarts, AtomImageStarts, AtomImagePhase,
                sparse_SHOsum.rowStart(), sparse_SHOsum.colIndex(), collect);

        return 0; // TODO estimate the total number of floating point operations
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
                int constexpr Re = 0;
                if (1 == R1C2) { // is real
                    // version for Noco==1, R1C2==1, more readable
                    double cad{0};
                    for (int aj = 0; aj < nSHO; ++aj) {
                        auto const cpr = double(apc[(a0 + aj)*nrhs + irhs][Re][0][ivec]); // load projection coefficient
                        auto const am = AtomMat[ai*nSHO + aj];                            // load matrix element
                        cad += am * cpr;
                    } // aj
                    aac[(a0 + ai)*nrhs + irhs][Re][0][ivec] = cad;                        // store addition coefficient

                } else {
                    // matrix layout: AtomMat[Noco*Noco*R1C2*nSHO*nSHO]
                    assert(2 == R1C2); // is complex
                    int constexpr Im = R1C2 - 1; // index for imaginary part
                    std__complex<double> cad = 0;
                    for (int spjn = 0; spjn < Noco; ++spjn) {
                        for (int aj = 0; aj < nSHO; ++aj) {
                            // load projection coefficient
                            auto const cpr = std__complex<double>(
                                       apc[(a0 + aj)*nrhs + irhs][Re][spjn][ivec],
                                       apc[(a0 + aj)*nrhs + irhs][Im][spjn][ivec]);
                            // load matrix element
                            auto const am = std__complex<double>(
                                       AtomMat[(((spin*Noco + spjn)*R1C2 + 0)*nSHO + ai)*nSHO + aj],
                                       AtomMat[(((spin*Noco + spjn)*R1C2 + 1)*nSHO + ai)*nSHO + aj]);
                            cad += am * cpr;
                        } // aj
                    } // spjn
                    // store addition coefficient
                    aac[(a0 + ai)*nrhs + irhs][Re][spin][ivec] = cad.real();
                    aac[(a0 + ai)*nrhs + irhs][Im][spin][ivec] = cad.imag();

                } // 1 == R1C2
            } // ai

        }} // thread loops and block loops

    } // SHOmul

    template <typename real_t, int R1C2=2, int Noco=1, int n64=64>
    size_t __host__ SHOmul_driver(
          real_t         (*const __restrict__ aac)[R1C2][Noco][Noco*n64] // result,  atom addition coefficients
        , real_t   const (*const __restrict__ apc)[R1C2][Noco][Noco*n64] // input, atom projection coefficients
        , double   const (*const *const __restrict__ AtomMatrices)
        , int8_t   const (*const __restrict__ AtomLmax)
        , uint32_t const (*const __restrict__ AtomStarts)
        , int      const natoms
        , int      const nrhs // number of block columns in the Green function
        , int const echo=0 // log level
    ) {
        assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));

        dim3 const gridDim(natoms, nrhs, 1), blockDim(Noco*n64, Noco, 1);
        if (echo >-1) std::printf("# %s<%s,R1C2=%d,Noco=%d,%d> <<< {natoms=%d, nrhs=%d, 1}, {%d, %d, 1} >>>\n",
                           __func__, real_t_name<real_t>(), R1C2, Noco, n64,  natoms, nrhs,  Noco*n64, Noco);
        SHOmul<real_t,R1C2,Noco,n64> // launch <<< {natoms, nrhs, 1}, {Noco*n64, Noco, 1} >>>
#ifndef HAS_NO_CUDA
            <<< gridDim, blockDim >>> (
#else //  HAS_NO_CUDA
           (    gridDim, blockDim,
#endif // HAS_NO_CUDA
            aac, apc, AtomMatrices, AtomLmax, AtomStarts);

        return 0; // TODO estimate the total number of floating point operations
    } // SHOmul_driver




            
  class dyadic_plan_t {
  public: // members

      uint32_t* AtomStarts     = nullptr; // [nAtoms + 1]
      int8_t*   AtomLmax       = nullptr; // [nAtoms]
      double**  AtomMatrices   = nullptr; // [nAtoms][2*nc^2] atomic matrices, nc number of SHO coefficients of this atom
      uint32_t nAtoms = 0;

      uint32_t* AtomImageStarts     = nullptr; // [nAtomImages + 1]
      double  (*AtomImagePos)[3+1]  = nullptr; // [nAtomImages]
      int8_t*   AtomImageLmax       = nullptr; // [nAtomImages]
      double  (*AtomImagePhase)[4]  = nullptr; // [nAtomImages]
      int8_t  (*AtomImageShift)[4]  = nullptr; // [nAtomImages]
      uint32_t nAtomImages  = 0;

      double *grid_spacing   = nullptr; // [3+1] hx,hy,hz,rcut/sigma

      int32_t nrhs = 0;
      green_sparse::sparse_t<> * sparse_SHOprj = nullptr; // [nrhs]
      green_sparse::sparse_t<>   sparse_SHOadd;
      green_sparse::sparse_t<>   sparse_SHOsum;

  public: // methods

      dyadic_plan_t() {} // default constructor

      dyadic_plan_t(int const echo) {
          if (echo > 0) std::printf("# construct %s\n", __func__);
      } // constructor

      ~dyadic_plan_t() {
          std::printf("# destruct %s\n", __func__);
          free_memory(AtomImageStarts);
          free_memory(AtomStarts);
          if (AtomMatrices) for (uint32_t ia = 0; ia < nAtoms; ++ia) free_memory(AtomMatrices[ia]);
          free_memory(AtomMatrices);
          free_memory(grid_spacing);
          free_memory(AtomImagePos);
          free_memory(AtomImageLmax);
          free_memory(AtomImagePhase);
          free_memory(AtomImageShift);
          free_memory(AtomLmax);
          if (sparse_SHOprj) for (int32_t irhs = 0; irhs < nrhs; ++irhs) sparse_SHOprj[irhs].~sparse_t<>();
          free_memory(sparse_SHOprj);
      } // constructor
      
      status_t consistency_check() const {
          status_t stat(0);
          assert(nAtomImages >= nAtoms);
          auto const rowStart = sparse_SHOsum.rowStart();
          auto const colIndex = sparse_SHOsum.colIndex();
          stat += (sparse_SHOsum.nRows() != nAtoms);
          if (!rowStart) return stat;
          for (uint32_t ia = 0; ia < nAtoms; ++ia) {
              for (auto bsr = rowStart[ia]; bsr < rowStart[ia + 1]; ++bsr) {
                  auto const iai = colIndex[bsr];
                  stat += (AtomImageLmax[iai] != AtomLmax[ia]);
              } // bsr
          } // ia
          return stat; // number of errors
      } // consistency_check

  }; // dyadic_plan_t

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
    ) {
        assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));
        size_t nops{0};
        
        assert(p.AtomImageStarts);
        auto const natomcoeffs = p.AtomImageStarts[p.nAtomImages];
        if (echo > 6) std::printf("# %s<%s,R1C2=%d,Noco=%d> nAtoms=%d nAtomImages=%d nrhs=%d ncoeffs=%d\n",
                  __func__, real_t_name<real_t>(), R1C2, Noco, p.nAtoms, p.nAtomImages, p.nrhs, natomcoeffs);

        nops += SHOprj_driver<real_t,R1C2,Noco>(Cpr, psi, p.AtomImagePos, p.AtomImageLmax, p.AtomImageStarts, p.nAtomImages,
                                                p.sparse_SHOprj, RowIndexCubes, CubePos, p.grid_spacing, p.nrhs, echo);

        auto Cad = get_memory<real_t[R1C2][Noco][Noco*64]>(natomcoeffs*p.nrhs, echo, "Cad"); // rectangular storage of atoms x right-hand-sides

        auto const ncoeffs = p.AtomStarts[p.nAtoms];
        if (p.nAtomImages > p.nAtoms) {
            // we have to distinguish between atoms and their periodic images

            auto cprj = get_memory<double[R1C2][Noco][Noco*64]>(ncoeffs*p.nrhs, echo, "cprj");

            SHOsum_driver<real_t,R1C2,Noco>(cprj, Cpr, p.AtomLmax, p.AtomStarts, p.AtomImageStarts, p.AtomImagePhase,
                                            p.sparse_SHOsum, p.nAtoms, p.nrhs, true, echo); // true:collect

            auto cadd = get_memory<double[R1C2][Noco][Noco*64]>(ncoeffs*p.nrhs, echo, "cadd");
            
            nops += SHOmul_driver<double,R1C2,Noco>(cadd, cprj, p.AtomMatrices, p.AtomLmax, p.AtomStarts, p.nAtoms, p.nrhs, echo);

            free_memory(cprj);

            SHOsum_driver<real_t,R1C2,Noco>(cadd, Cad, p.AtomLmax, p.AtomStarts, p.AtomImageStarts, p.AtomImagePhase,
                                            p.sparse_SHOsum, p.nAtoms, p.nrhs, false, echo); // false:broadcast

            free_memory(cadd);

        } else {
            assert(p.nAtomImages == p.nAtoms); // there is at most one relevant image of each atom so we can ignore phases
            assert(natomcoeffs == ncoeffs);

            nops += SHOmul_driver<real_t,R1C2,Noco>(Cad, Cpr, p.AtomMatrices, p.AtomLmax, p.AtomStarts, p.nAtoms, p.nrhs, echo);

        } // more images than atoms

        nops += SHOadd_driver<real_t,R1C2,Noco>(Ppsi, Cad, p.AtomImagePos, p.AtomImageLmax, p.AtomImageStarts, 
                                                p.sparse_SHOadd.rowStart(), p.sparse_SHOadd.colIndex(),
                                                RowIndexCubes, ColIndexCubes, CubePos, p.grid_spacing, nnzb, p.nrhs, echo);
        free_memory(Cad);

        return nops; // total number of floating point operations performed
    } // multiply (dyadic operations)


    std::vector<double> sho_normalization(int8_t const lmax, double const sigma=1) {
        int const n1ho = sho_tools::n1HO(lmax);
        std::vector<double> v1(n1ho, constants::sqrtpi * sigma);
        double fac{1};
        for (int nu = 0; nu <= lmax; ++nu) {
            v1[nu] *= fac;
            fac *= 0.5*(nu + 1);
        } // nu

        std::vector<double> vec(sho_tools::nSHO(lmax), 1.);
        int sho{0};
        for (int iz = 0; iz <= lmax; ++iz) {
            for (int iy = 0; iy <= lmax - iz; ++iy) {
                for (int ix = 0; ix <= lmax - iz - iy; ++ix) {
                    vec[sho] = v1[ix] * v1[iy] * v1[iz];
                    ++sho;
        }}}
        assert(vec.size() == sho);
        return vec;
    } // sho_normalization
            

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS


  inline status_t test_Hermite_polynomials_1D(int const echo=0, double const sigma=16) {
      int constexpr Lmax = 7;
      double H1D[Lmax + 1][3][4]; // 1D-Hermite polynomials times Gaussian decay function for 3 directions and grid points
      float  h1D[Lmax + 1][3][4]; // H1D in single precision
      float     xi_squared[3][4]; // (x/sigma)^2
      float  const xyzc[] = {0, 0, 0}; // cube position vector [in units of 4 grid points]
      double const hxyz[] = {1, 1, 1,   6.2832}; // grid spacing vector [in Bohr] and truncation_radius/sigma
      double       xyza[] = {0, 0, 0,   1./std::sqrt(sigma)}; // atom position vector [in Bohr] and spread^{-1/2}
      if (echo > 10) std::printf("\n## %s xi, H1D[0:%d]:", __func__, Lmax);
      double maxdev{0};
      for (int it = 0; it < 9; ++it) {
          for (int i3 = 0; i3 < 3; ++i3) {
              xyza[i3] = -4*(it*3 + i3);
              for (int i4 = 0; i4 < 4; ++i4) {
                  int const i = i3*4 + i4; // == threadIdx.x
                  Hermite_polynomials_1D(h1D, xi_squared, i, Lmax, xyza, xyzc, hxyz); // float
                  Hermite_polynomials_1D(H1D, xi_squared, i, Lmax, xyza, xyzc, hxyz); // double
//                if (echo > 10) std::printf("# %s idir=%d i4=%d ivec=%i xi^2=%g H1D= ", __func__, i3, i4, i, xi_squared[i3][i4]);
                  if (echo > 10) std::printf("\n%g  ", std::sqrt(xi_squared[i3][i4]));
                  for (int l = 0; l <= Lmax; ++l) {
                      if (echo > 10) std::printf(" %g", H1D[l][i3][i4]);
                      auto const dev = h1D[l][i3][i4] - H1D[l][i3][i4];
                      if (H1D[l][i3][i4] != 0.0) {
                          auto const reldev = std::abs(dev/H1D[l][i3][i4]);
                          maxdev = std::max(maxdev, reldev);
                      }
                  } // l
              } // i4
          } // i3
      } // it translations
      if (echo > 10) std::printf("\n\n\n");
      if (echo > 3) std::printf("# %s largest relative deviation is %.1e\n", __func__, maxdev);
      return 0;
  } // test_Hermite_polynomials_1D




    template <typename real_t, int R1C2=2, int Noco=1>
    size_t __host__ multiply(
          real_t         (*const __restrict__ Ppsi)[R1C2][Noco*64][Noco*64] // result,  modified Green function blocks [nnzb][R1C2][Noco*64][Noco*64]
        , real_t         (*const __restrict__  Cpr)[R1C2][Noco]   [Noco*64] // projection coefficients     [natomcoeffs*nrhs][R1C2][Noco   ][Noco*64]
        , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input, unmodified Green function blocks [nnzb][R1C2][Noco*64][Noco*64]
        , double   const (*const __restrict__ AtomPos)[3+1] // atomic positions and spread parameter
        , int8_t   const (*const __restrict__ AtomLmax) // SHO basis size [nAtomImages]
        , uint32_t const (*const __restrict__ AtomStarts) // prefetch sum over nSHO(AtomLmax[:])
        , uint32_t const natoms
        , green_sparse::sparse_t<> const (*const __restrict__ sparse_SHOprj)
        , double   const (*const          *const __restrict__ AtomMatrices)
        , green_sparse::sparse_t<> const &                    sparse_SHOadd
        , uint32_t const (*const __restrict__ RowIndexCubes) // Green functions rowIndex[nnzb]
        , uint16_t const (*const __restrict__ ColIndexCubes) // Green functions colIndex[nnzb]
        , float    const (*const __restrict__ CubePos)[3+1] // cube positions +alignment
        , double   const (*const __restrict__ hGrid) // grid spacings + projection radius
        , uint32_t const nnzb // number of blocks of the Green function
        , uint32_t const nrhs // number of block columns of the Green function
        , int const echo=0 // log-level
    ) {
        assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));
        size_t nops{0};

        auto const natomcoeffs = AtomStarts[natoms];
        if (echo > 6) std::printf("# %s<%s> R1C2=%d Noco=%d natoms=%d nrhs=%d ncoeffs=%d\n", __func__, real_t_name<real_t>(), R1C2, Noco, natoms, nrhs, natomcoeffs);

        nops += SHOprj_driver<real_t,R1C2,Noco>(Cpr, psi, AtomPos, AtomLmax, AtomStarts, natoms, sparse_SHOprj, RowIndexCubes, CubePos, hGrid, nrhs, echo);
        
        // ToDo: in the future, nAtomImages and natoms can be different. Then, we need an extra step of accumulating the projection coefficients
        //       with the corresponding phases and, after SHOmul, distributing the addition coefficients to the images with the inverse phases.        

        auto Cad = get_memory<real_t[R1C2][Noco][Noco*64]>(natomcoeffs*nrhs, echo, "Cad"); // rectangular storage of atoms x right-hand-sides

        // very rigid version, R1C2 and Noco could be function arguments instead of template parameters
        nops += SHOmul_driver<real_t,R1C2,Noco,64>(Cad, Cpr, AtomMatrices, AtomLmax, AtomStarts, natoms, nrhs, echo);

        nops += SHOadd_driver<real_t,R1C2,Noco>(Ppsi, Cad, AtomPos, AtomLmax, AtomStarts, sparse_SHOadd.rowStart(), sparse_SHOadd.colIndex(), RowIndexCubes, ColIndexCubes, CubePos, hGrid, nnzb, nrhs, echo);
        
        free_memory(Cad);
        
        return nops; // total number of floating point operations performed
    } // multiply (deprecated)



  template <typename real_t, int R1C2=2, int Noco=1>
  inline status_t test_SHOprj_and_SHOadd(int const echo=0, double const sigma=1, int8_t const lmax=5) {
      // check if drivers compile and the normalization of the lowest (up to 64) SHO functions
      auto const rc = control::get("green_dyadic.test.rc", 7.);
      int  const nb = control::get("green_dyadic.test.nb", 11.), natoms = 1, nrhs = 1, nnzb = pow3(nb);
      int  const nsho = sho_tools::nSHO(lmax);
      auto psi = get_memory<real_t[R1C2][Noco*64][Noco*64]>(nnzb, echo, "psi");
      set(psi[0][0][0], nnzb*R1C2*pow2(Noco*64ull), real_t(0)); // clear
      auto apc = get_memory<real_t[R1C2][Noco   ][Noco*64]>(natoms*nsho*nrhs, echo, "apc");
      set(apc[0][0][0], natoms*nsho*nrhs*R1C2*pow2(Noco)*64, real_t(0)); // clear

      auto sparse_SHOprj = get_memory<green_sparse::sparse_t<>>(nrhs, echo, "sparse_SHOprj");
      {   
          std::vector<std::vector<uint32_t>> SHO_prj(natoms);
          SHO_prj[0].resize(nnzb); for (int inzb = 0; inzb < nnzb; ++inzb) SHO_prj[0][inzb] = inzb;
          sparse_SHOprj[0] = green_sparse::sparse_t<>(SHO_prj, false, __func__, echo - 9);
      }
      green_sparse::sparse_t<> sparse_SHOadd;
      {
          std::vector<std::vector<uint32_t>> SHO_add(nnzb);
          for (int inzb = 0; inzb < nnzb; ++inzb) SHO_add[inzb].resize(1, 0);
          sparse_SHOadd    = green_sparse::sparse_t<>(SHO_add, false, __func__, echo - 9);
      }

      auto ColIndexCubes = get_memory<uint16_t>(nnzb, echo, "ColIndexCubes");     set(ColIndexCubes, nnzb, uint16_t(0));
      auto RowIndexCubes = get_memory<uint32_t>(nnzb, echo, "RowIndexCubes");     for (int inzb = 0; inzb < nnzb; ++inzb) RowIndexCubes[inzb] = inzb;
      auto hGrid         = get_memory<double>(3+1, echo, "hGrid");                set(hGrid, 3, 0.25); hGrid[3] = rc;
      auto AtomPos       = get_memory<double[3+1]>(natoms, echo, "AtomPos");      set(AtomPos[0], 3, hGrid, 0.5*4*nb);  AtomPos[0][3] = 1./std::sqrt(sigma);
      auto AtomLmax      = get_memory<int8_t>(natoms, echo, "AtomLmax");          set(AtomLmax, natoms, lmax);
      auto AtomStarts    = get_memory<uint32_t>(natoms + 1, echo, "AtomStarts");  for(int ia = 0; ia <= natoms; ++ia) AtomStarts[ia] = ia*nsho;
      auto CubePos       = get_memory<float[3+1]>(nnzb, echo, "CubePos");
      for (int iz = 0; iz < nb; ++iz) {
      for (int iy = 0; iy < nb; ++iy) {
      for (int ix = 0; ix < nb; ++ix) {
          int const xyz0[] = {ix, iy, iz, 0};
          set(CubePos[(iz*nb + iy)*nb + ix], 4, xyz0);
      }}} // ix iy iz

      // see if these drivers compile and can be executed without segfaults
      if (echo > 11) std::printf("# here %s:%d\n", __func__, __LINE__);

      auto const dV = hGrid[0]*hGrid[1]*hGrid[2];
      auto const sho_norm = sho_normalization(lmax, sigma);
      for (int isho = 0; isho < std::min(nsho, 64); ++isho) apc[isho*nrhs][0][0][isho] = dV/sho_norm[isho]; // set "unit matrix" but normalized

      SHOadd_driver<real_t,R1C2,Noco>(psi, apc, AtomPos, AtomLmax, AtomStarts, sparse_SHOadd.rowStart(), sparse_SHOadd.colIndex(), RowIndexCubes, ColIndexCubes, CubePos, hGrid, nnzb, nrhs, echo);
      if (0) {
          size_t nz{0};
          for (int i = 0; i < nnzb*64; ++i) {
              auto const value = psi[i >> 6][0][i & 63][0];
              nz += (0 == value);
              if (echo > 19) std::printf(" %g", value);
          } // isho
          if (echo > 3) std::printf("\n# %ld non-zeros of %d\n", nz, nnzb*64);
      } // echo

      SHOprj_driver<real_t,R1C2,Noco>(apc, psi, AtomPos, AtomLmax, AtomStarts, natoms, sparse_SHOprj, RowIndexCubes, CubePos, hGrid, nrhs, echo);

      double maxdev[2] = {0, 0}; // {off-diagonal, diagonal}
      { // scope
          auto const msho = std::min(nsho, 64); // show projection coefficients
          if (echo > 5) std::printf("# %d projection coefficients ", msho);
          for (int isho = 0; isho < msho; ++isho) {
              if (echo > 9) std::printf("\n# projection coefficients[%2d]: ", isho);
              for (int jsho = 0; jsho < msho; ++jsho) {
                  double const value = apc[isho*nrhs][0][0][jsho];
                  int const diag = (isho == jsho);
                  if (echo > 9) std::printf(" %g", value);
                  maxdev[diag] = std::max(maxdev[diag], std::abs(value - diag));
              } // isho
              if (echo > 9) std::printf(" diagonal=");
              if (echo > 5) std::printf(" %.3f", apc[isho*nrhs][0][0][isho]);
          } // isho
          if (echo > 5) std::printf("\n");
      } // scope
      if (echo > 2) std::printf("# %s<real%ld> orthogonality error %.1e, normalization error %.1e\n", __func__, sizeof(real_t), maxdev[0], maxdev[1]);

      if (0) {
          // also test the deprecated interface 'multiply'
          auto AtomMatrices = get_memory<double*>(natoms, echo, "AtomMatrices");
          for (int ia = 0; ia < natoms; ++ia) {
              if (echo > 9) std::printf("# %s atom image #%d has lmax= %d and %d coefficients starting at %d\n", __func__, ia, AtomLmax[ia], nsho, AtomStarts[ia]);
              AtomMatrices[ia] = get_memory<double>(pow2(Noco)*2*pow2(nsho), echo, "AtomMatrix[ia]");
              set(AtomMatrices[ia], pow2(Noco)*2*pow2(nsho), 0.0);
          } // ia
          multiply<real_t,R1C2,Noco>(psi, apc, psi, AtomPos, AtomLmax, AtomStarts, natoms, sparse_SHOprj, AtomMatrices, sparse_SHOadd, RowIndexCubes, ColIndexCubes, CubePos, hGrid, nnzb, nrhs, echo);
          for (int ia = 0; ia < natoms; ++ia) free_memory(AtomMatrices[ia]);
          free_memory(AtomMatrices);
      } // 0

      free_memory(ColIndexCubes);
      free_memory(RowIndexCubes);
      free_memory(sparse_SHOprj);
      free_memory(CubePos);
      free_memory(AtomStarts);
      free_memory(AtomLmax);
      free_memory(AtomPos);
      free_memory(apc);
      free_memory(psi);
      return 0;
  } // test_SHOprj_and_SHOadd

  inline status_t test_SHOprj_and_SHOadd(int const echo=0, double const sigma=1) {
      status_t stat(0);
      stat += test_SHOprj_and_SHOadd<float ,1,1>(echo, sigma); // real
//       stat += test_SHOprj_and_SHOadd<float ,2,1>(echo, sigma); // complex
//       stat += test_SHOprj_and_SHOadd<float ,2,2>(echo, sigma); // non-collinear
      stat += test_SHOprj_and_SHOadd<double,1,1>(echo, sigma); // real
//       stat += test_SHOprj_and_SHOadd<double,2,1>(echo, sigma); // complex
//       stat += test_SHOprj_and_SHOadd<double,2,2>(echo, sigma); // non-collinear
      return stat;
  } // test_SHOprj_and_SHOadd

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
//    stat += test_Hermite_polynomials_1D(echo);
      stat += test_SHOprj_and_SHOadd(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_dyadic

// 
// General thoughts:
//      
//    We assume that a Green function G(r,r') comes in as psi.
//    The spatial arguments \vec r of the Green function are
//    sampled on an equidistant Cartesian real-space grids.
//    Always 4x4x4 grid points are grouped to a "cube".
//
//    The original codes benchmarked in the PASC19 paper by Baumeister & Tsukamoto
//    had the following data layout.
//
//        projection/addition coefficients layout[natoms*nSHO*nrhs][Noco*64]
//        wave functions                   layout[ncubes* 64 *nrhs][Noco*64]
//
//    However, we will have to reformulate to
//        projection/addition coefficients layout[natomcoeffs*nrhs*R1C2*Noco][Noco*64]
//        Green function                   layout[nnzb       ][R1C2][Noco*64][Noco*64]
//
//    where natomcoeffs == natoms*nSHO(numax) only if all atoms have the same numax
//    and   nnzb == ncubes*nrhs only if the Green function is dense.
//
//    For considerations of geometry(*), the coefficient matrix will probably stay dense (natomcoeffs \times nrhs)
//    Then, we can still  launch SHOprj <<< {nrhs, natoms, 1}, {Noco*64, R1C2*Noco, 1} >>> providing
//                                   bsr in RowStartAtoms[irhs][iatom +0 .. +1], icube = irhs_of_inzb[irhs][bsr] --> one sparse_t per irhs
//    But we will need to launch SHOadd <<< {nnzb,    1,   1}, {Noco*64, R1C2*Noco, 1} >>> providing irhs = colIndex[inzb], 
//                                  irow = rowIndex[inzb] and iatom = ColIndexAtoms[irow][bsr], ColIndexAtoms as before
//
//    (*) The geometry consideration assumes the case that there is a compact cluster of
//        right hand side (source) cubes assigned to one MPI-process and isolated boundary conditions
//        for the cell. If the radius of target cubes is large compared to the cluster extent
//        the number of zero blocks inside the coefficient matrix is small compared to its size
//        since also each atomic projection region typically spans over several cubes.
//
