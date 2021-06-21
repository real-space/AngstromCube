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
    Hermite_1D_polys(
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
    } // Hermite_1D_polys

    
    
    template <typename real_t, int nvec=64, int Lmax=7>
    void __global__ SHOprj( // launch SHOprj<real_t,nvec> <<< {nRHSs/nvec, natoms, 1}, {nvec} >>> (...);
          real_t        (*const __restrict__ Cpr)[nvec] // result: projection coefficients
        , real_t  const (*const __restrict__ Psi)[nvec] // input:  wave functions
        , double  const (*const __restrict__ AtomPos)[4] // atomic positions [0],[1],[2] and decay parameter [3]
        , int     const (*const __restrict__ RowStart) // rows==atoms
        , int     const (*const __restrict__ ColIndex) // cols==cubes
        , int64_t const (*const __restrict__ BitMask) // can be omitted in the future
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
        nrhs = gridDim.z;
        { int const block_y = blockIdx.y;
        { int const block_x = blockIdx.x;
        { int const thread_x = threadIdx.x;
#endif

        int const iatom = block_y;  // in [0, natoms)
        if (RowStart[iatom] >= RowStart[iatom + 1]) return; // empty range in the sparse matrix, early return

        int const lmax = Lmax; // ToDo: make this atom-dependent in the future
        assert(lmax <= Lmax);

        int const J = block_x;  // in [0, nRHSs/nvec) ceiling
        int const j = thread_x; // in [0, nvec)

        __shared__ real_t H1D[Lmax + 1][3][4]; // non-normalized orthogonal 1-dimensional Hermite Gauss functions
        __shared__ float  xi_squared[3][4];
        __shared__ double xyza[4]; // atom position
        __shared__ float  xyzc[4]; // cube position
        __shared__ float  hgrid[4];

        // stage the atomic position of the atom treated in this CUDA block
        if (thread_x < 4) {
            xyza[thread_x] = AtomPos[iatom][thread_x]; // coalesced load
            hgrid[thread_x] = hGrid[thread_x]; // load grid spacing into shared memory
        } // lowest 4 threads

        real_t c[sho_tools::nSHO(Lmax)]; // get nSHO accumulator registers
        for (int sho = 0; sho < sho_tools::nSHO(Lmax); ++sho) {
            c[sho] = 0; // init accumulator registers
        } // sho

        __syncthreads(); // sync necessary?

        for (int bsr = RowStart[iatom]; bsr < RowStart[iatom + 1]; ++bsr) {

            __syncthreads();

            int const icube = ColIndex[bsr];

            if (thread_x < 4) {
                xyzc[thread_x] = CubePos[icube][thread_x]; // load into shared memory
            } // lowest 4 threads

            __syncthreads();

            // generate Hx, Hy, Hz up to lmax inside the 4^3 cube
            auto const R2_proj = Hermite_1D_polys(H1D, xi_squared, j, lmax, xyza, xyzc, hgrid);

            int64_t const mask_atom = BitMask[bsr]; // needs 2 GPU registers (32bit)

            __syncthreads();

            for (int z = 0; z < 4; ++z) { // loop over real-space grid in z-direction
                int const mask_z = (mask_atom >> (16*z)) & 0xffff; // allows to work on a 32bit register using 16bit
                // auto const d2z = xi_squared[2][z] - R2_proj;
                if (0 != mask_z) { // alternative: if (d2z < 0)

                    real_t b[sho_tools::n2HO(Lmax)];
                    __unroll__
                    for (int iyx = 0; iyx < sho_tools::n2HO(Lmax); ++iyx) { // loop over the combined index of ix and iy
                        b[iyx] = 0; // init
                    } // iyx

                    for (int y = 0; y < 4; ++y) { // loop over real-space grid in y-direction
                        int const mask_y = (mask_z >> (4*y)) & 0xf; // allows to work on a 32bit register using 4bit
                        // auto const d2yz = xi_squared[1][y] + d2z;
                        if (0 != mask_y) { // alterative: if (d2yz < 0)

                            real_t a[Lmax + 1];
                            __unroll__
                            for (int ix = 0; ix <= Lmax; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                a[ix] = 0; // init
                            } // ix

                            __unroll__
                            for (int x = 0; x < 4; ++x) { // loop over real-space grid in x-direction
                                int const mask_x = (mask_y >> (1*x)) & 0x1; // uses only 1 bit
                                // auto const d2xyz = xi_squared[0][x] + d2yz;
                                if (0 != mask_x) { // alterative: if (d2xyz < 0)
                                    int const xyz = (z*4 + y)*4 + x;
                                    real_t const ps = Psi[(icube*64 + xyz)*nrhs + J][j]; // load from global memory
                                    for (int ix = 0; ix <= lmax; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                        auto const Hx = H1D[ix][0][x]; // load Hx from shared memory
                                        a[ix] += ps * Hx; // FMA: 2 flop * 4**3 * (L+1)
                                    } // ix
                                } // 0 != mask_x
                            } // x

                            for (int iy = 0; iy <= lmax; ++iy) { // loop over the 2nd Cartesian SHO quantum number
                                int const iyx0 = (iy*(2*lmax + 3 - iy)) >> 1;
                                for (int ix = 0; ix <= lmax - iy; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                    auto const Hy = H1D[iy][1][y]; // load Hy from shared memory
                                    b[iyx0 + ix] += a[ix] * Hy; // FMA: 2 flop * 4**2 * ((L+1)*(L+2))/2
                                } // ix
                            } // iy

                        } // 0 != mask_y
                    } // y

                    for (int sho = 0, iz = 0; iz <= lmax; ++iz) { // loop over the 3rd Cartesian SHO quantum number
                        auto const Hz = H1D[iz][2][z]; // load Hz from shared memory
                        for (int iy = 0; iy <= lmax - iz; ++iy) { // loop over the 2nd Cartesian SHO quantum number
                            int const iyx0 = (iy*(2*lmax + 3 - iy)) >> 1;
                            for (int ix = 0; ix <= lmax - iz - iy; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                c[sho] += b[iyx0 + ix] * Hz; // FMA: 2 flop * 4**1 * ((L+1)*(L+2)*(L+3))/6
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
                Cpr[(iatom*nSHO + sho)*nrhs + J][j] = c[sho];
            } // sho
        } // store

        }}} // close 3 loops

    } // SHOprj
    
#undef __global__

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
