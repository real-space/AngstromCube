#pragma once

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cmath> // std::sqrt
#include <algorithm> // std::max
#include <utility> // std::swap //, std::move
#include <vector> // std::vector<T>

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_memory.hxx" // get_memory, free_memory
#include "data_view.hxx" // view3D<T>
#include "inline_math.hxx" // set
#include "simple_stats.hxx" // ::Stats<>
#include "sho_tools.hxx" // ::nSHO, ::n2HO, ::n1HO

namespace green_projection {
  
#define NO_CUDA
#define __global__
#define ____restrict____
#define __device__
#define __shared__
#define __unroll__
#define __host__

  
    template <typename real_t>
    inline float __host__ __device__
    Hermite_polynomials_1D(
          real_t (*const __restrict__ H1D)[3][4] // result
          float  (*const __restrict__ xi_squared)[4] // distance^2
        , int    const ivec // thread index used for the vectorization
        , int    const Lmax
        , double const xyza[4] // atomic position in [0],[1] and [2], sigma^{-1/2} in [3]
        , float  const xyzc[3] // position of the target block
        , float  const hxyz[4] // grid spacings in [0],[1] and [2], projection radius in [3]
    ) {

        double const R2_projection = pow2(double(hxyz[3])); // Rprj can be controlled from outside

        if (Lmax < 0) return R2_projection;

        int constexpr l2b = 2;
        int constexpr n4 = 1 << l2b; // Block edge length is 2^l2b
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
            double const cube_position = xyzc[idir]; // position of the target block
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

            double Hnp1, Hn = H0, Hnm1{0};
            for (int nu = 0; nu < Lmax; ++nu) {
                double const nuhalf = 0.5 * nu; // nu/2.
                // use the two step recurrence relation to create the 1D Hermite polynomials
                // H1D[nu+1][idir][i4] = xi * H1D[nu][idir][i4] - nu/2. * H1D[nu-1][idir][i4];
                Hnp1 = xi * Hn - nuhalf * Hnm1; // see snippets/3Hermite.F90 for the derivation of this
                H1D[nu + 1][idir][i4] = Hnp1; // store
                // rotate registers for the next iteration
                Hnm1 = Hn; Hn = Hnp1; // ordering is important here
            } // nu

            // Warning!
            // The Hermite polynomials are not normalized. To normalize w.r.t. the 2-norm,
            // we must divide by sqrt{sqrt{pi}*sigma*factorial(nu)} in an additional loop.
            // But this is not necessary for their orthogonality

        } // ivec < 3*n4

        // here, start the tree reduction over the 5 elements for the grid-refinement

        return R2_projection;
    } // Hermite_polynomials_1D

    
    int constexpr Lmax_default=7; // reduce this to lower the GPU register usage
    
    template <typename real_t, int nvec=64, int Lmax=Lmax_default>
    void __global__ SHOprj( // launch SHOprj<real_t,nvec> <<< {nRHSs/nvec, natoms, 1}, {nvec} >>> (...);
          real_t        (*const __restrict__ Cpr)[nvec] // result: projection coefficients
        , real_t  const (*const __restrict__ Psi)[nvec] // input:  wave functions
        , double  const (*const __restrict__ AtomPos)[4] // atomic positions [0],[1],[2] and decay parameter [3]
        , int     const (*const __restrict__ RowStart) // rows==atoms
        , int     const (*const __restrict__ ColIndex) // cols==cubes
//      , int64_t const (*const __restrict__ BitMask) // can be omitted
        , float   const (*const __restrict__ CubePos)[4] // only [0],[1],[2] used
        , float   const (*const __restrict__ hGrid) // grid spacings in [0],[1] and [2], projection radius in [3]
        , int nrhs=0
        , int const natoms=0
    ) {

#ifdef NO_CUDA
        for (int block_y = 0; block_y < natoms; ++block_y) {
        for (int block_x = 0; block_x < nrhs; ++block_x) {
        for (int thread_x = 0; thread_x < nvec; ++thread_x) {
        assert(false); // does not work without CUDA due to usage of shared memory elements
#else
        assert(nvec == blockDim.x);
        nrhs = gridDim.x;
        { int const block_y = blockIdx.y;
        { int const block_x = blockIdx.x;
        { int const thread_x = threadIdx.x;
#endif

        int const iatom = block_y; // in [0, natoms)
        if (RowStart[iatom] >= RowStart[iatom + 1]) return; // empty range in the sparse matrix, early return

        int const lmax = Lmax; // ToDo: make this atom-dependent in the future
        assert(lmax <= Lmax);

        int const J = block_x;  // in [0, nRHSs/nvec) ceiling
        int const j = thread_x; // in [0, nvec)

        __shared__ real_t H1D[Lmax + 1][3][4]; // non-normalized orthogonal 1-dimensional Hermite Gauss functions
        __shared__ float  xi_squared[3][4]; // distance along one Cartesian direction squared
        __shared__ double xyza[4]; // atom position
        __shared__ float  xyzc[4]; // cube position
        __shared__ float  hgrid[4];

        // stage the atomic position of the atom treated in this CUDA block
        if (thread_x < 4) {
            xyza[thread_x] = AtomPos[iatom][thread_x]; // coalesced load
            hgrid[thread_x] = hGrid[thread_x]; // load grid spacing into shared memory
        } // lowest 4 threads

        real_t czyx[sho_tools::nSHO(Lmax)]; // get nSHO accumulator registers
        for (int sho = 0; sho < sho_tools::nSHO(Lmax); ++sho) {
            czyx[sho] = 0; // init accumulator registers
        } // sho

        __syncthreads(); // sync necessary?

        // in this loop over bsr we accumulate atom projection coefficients over many cubes
        for (int bsr = RowStart[iatom]; bsr < RowStart[iatom + 1]; ++bsr) {

            __syncthreads();

            int const icube = ColIndex[bsr];

            if (thread_x < 4) {
                xyzc[thread_x] = CubePos[icube][thread_x]; // load into shared memory
            } // lowest 4 threads

            __syncthreads();

            // generate Hx, Hy, Hz up to lmax inside the 4^3 cube
            auto const R2_proj = Hermite_polynomials_1D(H1D, xi_squared, j, lmax, xyza, xyzc, hgrid);

//          int64_t const mask_atom = BitMask[bsr]; // needs 2 GPU registers (32bit each)
 
            __syncthreads();

            for (int z = 0; z < 4; ++z) { // loop over real-space grid in z-direction
//              int const mask_z = (mask_atom >> (16*z)) & 0xffff; // allows to work on a 32bit register using 16bit
//              if (0 != mask_z) {
                auto const d2z = xi_squared[2][z] - R2_proj;
                if (d2z < 0) {

                    real_t byx[sho_tools::n2HO(Lmax)];
                    __unroll__
                    for (int iyx = 0; iyx < sho_tools::n2HO(Lmax); ++iyx) { // loop over the combined index of ix and iy
                        byx[iyx] = 0; // init
                    } // iyx

                    for (int y = 0; y < 4; ++y) { // loop over real-space grid in y-direction
//                      int const mask_y = (mask_z >> (4*y)) & 0xf; // allows to work on a 32bit register using 4bit
//                      if (0 != mask_y) {
                        auto const d2yz = xi_squared[1][y] + d2z;
                        if (d2yz < 0) {

                            real_t a[Lmax + 1];
                            __unroll__
                            for (int ix = 0; ix <= Lmax; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                ax[ix] = 0; // init
                            } // ix

                            __unroll__
                            for (int x = 0; x < 4; ++x) { // loop over real-space grid in x-direction
//                              int const mask_x = (mask_y >> (1*x)) & 0x1; // uses only 1 bit
//                              if (0 != mask_x) {
                                auto const d2xyz = xi_squared[0][x] + d2yz;
                                if (d2xyz < 0) {
                                    int const xyz = (z*4 + y)*4 + x;
                                    real_t const ps = Psi[(icube*64 + xyz)*nrhs + J][j]; // load from global memory
                                    for (int ix = 0; ix <= lmax; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                        auto const Hx = H1D[ix][0][x]; // load Hx from shared memory
                                        ax[ix] += ps * Hx; // FMA: 2 flop * 4**3 * (L+1)
                                    } // ix
                                } // 0 != mask_x
                            } // x

                            for (int iy = 0; iy <= lmax; ++iy) { // loop over the 2nd Cartesian SHO quantum number
                                int const iyx0 = (iy*(2*lmax + 3 - iy)) >> 1;
                                for (int ix = 0; ix <= lmax - iy; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                    auto const Hy = H1D[iy][1][y]; // load Hy from shared memory
                                    byx[iyx0 + ix] += ax[ix] * Hy; // FMA: 2 flop * 4**2 * ((L+1)*(L+2))/2
                                } // ix
                            } // iy

                        } // 0 != mask_y
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

                } // 0 != mask_z
            } // z
            // now how many flop have been done?
            // FMA*( 4**3 * ((L+1)) + 4**2 * ((L+1)*(L+2))/2 + 4**1 * ((L+1)*(L+2)*(L+3))/6 )
            // so for Lmax=5 this is 944 FMAs = 1888 flop

        } // bsr

        { // store accumulators
            int const nSHO = sho_tools::nSHO(Lmax);
            // ToDo: use compressed storage to allow for lmax[iatom] without waisting memory for Cpr
            for (int sho = 0; sho < sho_tools::nSHO(lmax); ++sho) {
                Cpr[(iatom*nSHO + sho)*nrhs + J][j] = czyx[sho];
            } // sho
        } // store

        }}} // close 3 loops

    } // SHOprj
    
    
    
    
    // maybe this could be useful?
    union mask64_t {
        int64_t i;
        uint16_t u[4];
    }; // mask64_t

    
    template <typename real_t, int nvec=64, int Lmax=Lmax_default>
    void __global__ SHOadd( // launch SHOadd<real_t,nvec> <<< {nRHSs/nvec, ncubes, 1}, {nvec} >>> (...);
          real_t        (*const __restrict__ Psi)[nvec] // result: modified wave functions
        , real_t  const (*const __restrict__ Cad)[nvec] // input: addition coefficients
        , double  const (*const __restrict__ AtomPos)[4] // atomic positions [0],[1],[2] and decay parameter [3]
        , int     const (*const __restrict__ RowStart) // rows==cubes
        , int     const (*const __restrict__ ColIndex) // cols==atoms
        , float   const (*const __restrict__ CubePos)[4] // only [0],[1],[2] used
        , float   const (*const __restrict__ hGrid) // grid spacings in [0],[1] and [2], projection radius in [3]
        , int nrhs=0
        , int const ncubes=0
    ) {

#ifdef NO_CUDA
        for (int block_y = 0; block_y < ncubes; ++block_y) {
        for (int block_x = 0; block_x < nrhs; ++block_x) {
        for (int thread_x = 0; thread_x < nvec; ++thread_x) {
        assert(false); // does not work without CUDA due to usage of shared memory elements
#else
        assert(nvec == blockDim.x);
        nrhs = gridDim.x;
        { int const block_y = blockIdx.y;
        { int const block_x = blockIdx.x;
        { int const thread_x = threadIdx.x;
#endif

        int const icube = block_y; // in [0, ncubes)
        if (RowStart[icube] >= RowStart[icube + 1]) return; // empty range in the sparse matrix, early return

        int const lmax = Lmax; // ToDo: make this atom-dependent in the future
        assert(lmax <= Lmax);

        int const J = block_x;  // in [0, nrhs)
        int const j = thread_x; // in [0, nvec)

        __shared__ real_t H1D[1 + Lmax][3][4]; // non-normalized orthogonal 1-dimensional Hermite Gauss functions
        __shared__ float  xi_squared[3][4]; // distance along one Cartesian direction squared
        __shared__ double xyza[4]; // atom position
        __shared__ float  xyzc[4]; // cube position
        __shared__ float  hgrid[4]; // grid spacings

        // stage the atomic position of the atom treated in this CUDA block
        if (thread_x < 4) {
            xyzc[thread_x] = CubePos[icube][thread_x]; // coalesced load
            hgrid[thread_x] = hGrid[thread_x]; // load grid spacing into shared memory
        } // lowest 4 threads

        real_t czyx[4][4][4]; // get 64 accumulator registers --> 128 registers if double
        __unroll__
        for (int z = 0; z < 4; ++z) {
            __unroll__
            for (int y = 0; y < 4; ++y) {
                __unroll__
                for (int x = 0; x < 4; ++x) {
                    czyx[z][y][x] = 0; // init accumulator registers
                } // x
            } // y
        } // z

        int64_t mask_all{0}; // mask accumulator that tells which grid points have to be updated

        // in this loop over bsr we accumulate atom projector functions over many atoms
        for (int bsr = RowStart[icube]; bsr < RowStart[icube + 1]; ++bsr) {

            __syncthreads();

            int const iatom = ColIndex[bsr];

            if (thread_x < 4) {
                xyza[thread_x] = AtomPos[iatom][thread_x]; // load, convert to float
            } // lowest 4 threads

            __syncthreads();

            // generate Hx, Hy, Hz up to Lmax inside the 4^3 cube
            auto const R2_proj = Hermite_polynomials_1D(H1D, xi_squared, j, xyza, xyzc, hgrid);

            int sho{0};
            int const nSHO = sho_tools::nSHO(Lmax);
            for (int iz = 0; iz <= lmax; ++iz) { // executes (Lmax + 1) times

                real_t byx[4][4]; // get 16 accumulator registers --> 32 registers if double
                __unroll__
                for (int y = 0; y < 4; ++y) {
                    __unroll__
                    for (int x = 0; x < 4; ++x) {
                        byx[y][x] = 0; // init
                    } // x
                } // y

                for (int iy = 0; iy <= lmax - iz; ++iy) { // executes (Lmax + 1)*(Lmax + 2)/2 times

                    real_t ax[4] = {0, 0, 0, 0}; // get 4 accumulator registers --> 8 32bit GPU registers if double

                    for (int ix = 0; ix <= lmax - iz - iy; ++ix) { // executes (Lmax + 1)*(Lmax + 2)*(Lmax + 3)/6 times
                        real_t const ca = Cad[(iatom*nSHO + sho)*nrhs + J][j]; // load nvec numbers from global memory
                        ++sho;
                        __unroll__
                        for (int x = 0; x < 4; ++x) {
                            real_t const Hx = H1D[ix][0][x]; // load Hx from shared memory
                            ax[x] += Hx * ca; // FMA: 2 flop * 4 * (Lmax + 1)*(Lmax + 2)*(Lmax + 3)/6
                        } // x
                    } // ix

                    __unroll__
                    for (int y = 0; y < 4; ++y) {
                        real_t const Hy = H1D[iy][1][y]; // load Hy from shared memory
                        __unroll__
                        for (int x = 0; x < 4; ++x) {
                            byx[y][x] += Hy * ax[x]; // FMA: 2 flop * 4**2 * (Lmax + 1)*(Lmax + 2)/2
                        } // x
                    } // y

                } // iy

                int64_t mask_atom{0}; // occupy 2 GPU registers
                __unroll__
                for (int z = 0; z < 4; ++z) {
                    real_t const Hz = H1D[iz][2][z]; // load Hz from shared memory
                    auto const d2z = xi_squared[2][z] - R2_proj;
                    if (d2z < 0) {
                        __unroll__
                        for (int y = 0; y < 4; ++y) {
                            auto const d2yz = xi_squared[1][y] + d2z;
                            if (d2yz < 0) {
                                __unroll__
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
            if (lmax == Lmax) assert(nSHO == sho);
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
                            Psi[(icube*64 + xyz)*nrhs + J][j] += czyx[z][y][x]; // up to 64 negligible flop
                        } // mask
                    } // x
                } // y
            } // z
        } // mask_all

        }}} // close 3 loops
        
    } // SHOadd

#undef __global__


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

    template <typename real_t, int R1C2=2, int Noco=1>
    size_t multiply(
          real_t         (*const ____restrict____ apc)[R1C2][Noco][Noco*64] // result
        , real_t   const (*const ____restrict____ psi)[R1C2][Noco*64][Noco*64] // input
        , double   const hgrid=1 // grid spacing, ToDo make it an array of X,Y,Z
    ) {
      
        return 0ul; // total number of floating point operations performed
	} // multiply (projection operations)

#undef NO_CUDA

  
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_simple_projection(int const echo=0) {
      if (echo > 0) std::printf("# %s: no test included!\n", __func__);
      return STATUS_TEST_NOT_INCLUDED;
  } // test_simple_projection

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_simple_projection(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_projection
