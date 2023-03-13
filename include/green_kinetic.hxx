#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // uint32_t, int32_t, uint16_t, int16_t
#include <cassert> // assert

#include "status.hxx" // status_t
#include "green_memory.hxx" // get_memory, free_memory --> used in multiply, ToDo: move into planning phase
#include "green_sparse.hxx" // ::sparse_t<T>
#include "data_view.hxx" // view3D<T>

  // ToDo: move to green_utils.hxx or similar
  template <typename uint_t, typename int_t> inline
  size_t index3D(uint_t const n[3], int_t const i[3]) {
      // usual 3D indexing
      return size_t(i[2]*n[1] + i[1])*n[0] + i[0];
  } // index3D


namespace green_kinetic {

    int constexpr nhalo = 4; // a maximum of 4 blocks (i.e. 16 grid points) is the range of the FD stencil.

    int32_t constexpr BLOCK_EXISTS = 1;
    int32_t constexpr BLOCK_IS_ZERO = 0;
    int32_t constexpr BLOCK_NEEDS_PHASE = -1;
   

    template <typename real_t, int R1C2=2, int Noco=1> // Stride is determined by the lattice dimension along which we derive
    void __global__ Laplace8th( // GPU kernel, must be launched with <<< {Nrows, 16, 1}, {Noco*64, Noco, R1C2} >>>
#ifdef    HAS_NO_CUDA
          dim3 const & gridDim, dim3 const & blockDim,
#endif // HAS_NO_CUDA
          real_t        (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // intent(inout)
        , real_t  const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // intent(in)
        , int32_t const (*const *const __restrict__ index_list) // index list that brings the blocks in order,
                                                                // list must contain at least 4+1+4 elements
        , double const prefactor
        , int const Stride // 4^0, 4^1 or 4^2
        , double const phase[2][2]
    ) {
        // prepare finite-difference coefficients
        double const norm = prefactor/5040.; // {-14350, 8064, -1008, 128, -9}/5040. --> 8th order
        real_t const  c0 =     -14350*norm,
                      c1 =       8064*norm,
                      c2 =      -1008*norm,
                      c3 =        128*norm,
                      c4 =         -9*norm; //-> NStep = 4
        // the NStep=4 8th order stencil has a ratio of c0/c4 = 2^10.64 so half precision (11 significant bits) is the limit

        // the following list gives the indices of blocks that belong to the same right hand side
        // and are neighbors in the direction of derivation

        assert(16 == gridDim.y);
        assert(1  == gridDim.z);
        assert(Noco*64 == blockDim.x);
        assert(Noco    == blockDim.y);
        assert(R1C2    == blockDim.z);

#ifdef    HAS_NO_CUDA
        dim3 blockIdx(0,0,0);
        for (blockIdx.y = 0; blockIdx.y < gridDim.y; ++blockIdx.y)
        for (blockIdx.x = 0; blockIdx.x < gridDim.x; ++blockIdx.x)
#endif // HAS_NO_CUDA
        { // block loops

        auto const *const list = index_list[blockIdx.x]; // abbreviate pointer

        int const i16 = blockIdx.y;
        int const i64 = (16==Stride)? i16 : ( (4==Stride)? (16*(i16 >> 2) + (i16 & 0x3)) : (4*i16) );

#ifdef    HAS_NO_CUDA
        dim3 threadIdx(0,0,0);
        for (threadIdx.z = 0; threadIdx.z < blockDim.z; ++threadIdx.z)
        for (threadIdx.y = 0; threadIdx.y < blockDim.y; ++threadIdx.y)
        for (threadIdx.x = 0; threadIdx.x < blockDim.x; ++threadIdx.x)
#endif // HAS_NO_CUDA
        { // thread loops

#define INDICES(i4) [threadIdx.z][threadIdx.y*64 + i64 + Stride*i4][threadIdx.x]

        real_t w0{0}, w1{0}, w2{0}, w3{0}; // 4 registers for one block

        int ilist{nhalo - 1}; // counter for index_list, use only the last entry of the halo

        int ii = list[ilist++]; // the 1st block can be 0 (non-existent) or <0 for a periodic image
        // =========================================================================================
        // === periodic boundary conditions ========================================================
        if (BLOCK_IS_ZERO == ii) {
            // block does not exist (isolated/vacuum boundary condition), registers w0...w3 are initialized 0
        } else { // is periodic
            if (ii > 0) std::printf("# Error: ii= %d ilist= %i\n", ii, ilist);
            assert(ii <= BLOCK_NEEDS_PHASE && "list[3] must be either 0 (isolated BC) or negative (periodic BC)");
            int const jj = BLOCK_NEEDS_PHASE*(ii + BLOCK_EXISTS); // index of the periodic image of a block
            assert(jj >= 0); // must be a valid index to dereference psi[]
            assert(phase && "a phase must be given for complex BCs");
            real_t const ph_Re = phase[0][0]; // real part of the left complex phase factor
            // inital block load
            w0 = ph_Re * psi[jj]INDICES(0);
            w1 = ph_Re * psi[jj]INDICES(1);
            w2 = ph_Re * psi[jj]INDICES(2);
            w3 = ph_Re * psi[jj]INDICES(3);
            if (2 == R1C2) {
                real_t const ph_Im = phase[0][1] * (1 - 2*threadIdx.z); // imaginary part of the left complex phase factor
#define INDICES_Im(i4) [R1C2 - 1 - threadIdx.z][threadIdx.y*64 + i64 + Stride*i4][threadIdx.x]
                // inital load of imaginary part
                w0 -= ph_Im * psi[jj]INDICES_Im(0);
                w1 -= ph_Im * psi[jj]INDICES_Im(1);
                w2 -= ph_Im * psi[jj]INDICES_Im(2);
                w3 -= ph_Im * psi[jj]INDICES_Im(3);
            } // is complex
        } // is periodic
        // === periodic boundary conditions ========================================================
        // =========================================================================================

        real_t w4, w5, w6, w7, wn; // 4 + 1 registers, wn is the register that always receives the most recently loaded value

        ii = list[ilist++] - BLOCK_EXISTS; assert(ii >= 0); // this block must be a regular index since it is the 1st block for which we store a result
        // initially load one block in advance
        w4 = psi[ii]INDICES(0);
        w5 = psi[ii]INDICES(1);
        w6 = psi[ii]INDICES(2);
        w7 = psi[ii]INDICES(3);

        // main loop
        assert(5 == ilist);
        while (ii >= 0) {
            int const i0 = ii; // set central index
            ii = list[ilist++] - BLOCK_EXISTS; // get next index
//          if (0 == threadIdx.x && 0 == blockIdx.x) std::printf("# loop: ii= %d ilist= %i\n", ii, ilist);
            bool const load = (ii >= 0);
            // use a rotating register file, see figs/rotating_register_file.fig or *.pdf

            // FD9POINT = load?, compute, store, update rotating register file
#define FD9POINT(i4,  M4, M3, M2, M1, W0, P1, P2, P3, P4) \
            P4 = load ? psi[ii]INDICES(i4) : 0; \
            Tpsi[i0]INDICES(i4) += c0*W0 + c1*M1 + c1*P1 + c2*M2 + c2*P2 + c3*M3 + c3*P3 + c4*M4 + c4*P4; \
            M4 = P4;

            if (0 == (ilist & 0x1)) { // even
                FD9POINT(0,  w0, w1, w2, w3, w4, w5, w6, w7, wn)
                FD9POINT(1,  w1, w2, w3, w4, w5, w6, w7, w0, wn)
                FD9POINT(2,  w2, w3, w4, w5, w6, w7, w0, w1, wn)
                FD9POINT(3,  w3, w4, w5, w6, w7, w0, w1, w2, wn)
            } else {                // odd
                FD9POINT(0,  w4, w5, w6, w7, w0, w1, w2, w3, wn)
                FD9POINT(1,  w5, w6, w7, w0, w1, w2, w3, w4, wn)
                FD9POINT(2,  w6, w7, w0, w1, w2, w3, w4, w5, wn)
                FD9POINT(3,  w7, w0, w1, w2, w3, w4, w5, w6, wn)
            } // ilist is even or odd
#undef  FD9POINT
        } // while loop

        // =========================================================================================
        // === periodic boundary conditions ========================================================
        // correct for the tail part if periodic
        ii = list[ilist - 1]; // recover the list entry which stopped the while-loop
        if (BLOCK_IS_ZERO == ii) {
            // block does not exist (isolated/vacuum boundary condition), no further action necessary
        } else {
            assert(ii < 0 && "last list item must be either 0 (isolated BC) or negative (periodic BC)");
            int const jj = BLOCK_NEEDS_PHASE*(ii + BLOCK_EXISTS); // index of the periodic image of a block
            assert(jj >= 0); // must be a valid index to dereference psi[]
            int const i0 = list[ilist - 2] - BLOCK_EXISTS; // recover the last central index
            assert(i0 >= 0); // must be a valid index to dereference Tpsi[]
            assert(phase && "no (right) phase given"); // should have failed already above
            real_t const ph_Re = phase[1][0]; // real part of the right complex phase factor
            // final block load
            w0 = ph_Re * psi[jj]INDICES(0);
            w1 = ph_Re * psi[jj]INDICES(1);
            w2 = ph_Re * psi[jj]INDICES(2);
            w3 = ph_Re * psi[jj]INDICES(3);
            if (2 == R1C2) {
                real_t const ph_Im = phase[1][1] * (1 - 2*threadIdx.z); // imaginary part of the right complex phase factor
                // final load of imaginary part
                w0 -= ph_Im * psi[jj]INDICES_Im(0);
                w1 -= ph_Im * psi[jj]INDICES_Im(1);
                w2 -= ph_Im * psi[jj]INDICES_Im(2);
                w3 -= ph_Im * psi[jj]INDICES_Im(3);
#undef  INDICES_Im
            } // is complex
            // compute correction and add
            Tpsi[i0]INDICES(0) += c4*w0;
            Tpsi[i0]INDICES(1) += c3*w0 + c4*w1;
            Tpsi[i0]INDICES(2) += c2*w0 + c3*w1 + c4*w2;
            Tpsi[i0]INDICES(3) += c1*w0 + c2*w1 + c3*w2 + c4*w3;
        } // is periodic
        // === periodic boundary conditions ========================================================
        // =========================================================================================

#undef  INDICES
        }} // thread and block loops

//      std::printf("# %s(Stride=%d)\n", __func__, Stride);

    } // Laplace8th



    template <typename real_t, int R1C2=2, int Noco=1>
    void __global__ Laplace16th( // GPU kernel, must be launched with <<< {Nrows, 16, 1}, {Noco*64, Noco, R1C2} >>>
#ifdef    HAS_NO_CUDA
          dim3 const & gridDim, dim3 const & blockDim,
#endif // HAS_NO_CUDA
          real_t        (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // intent(inout)
        , real_t  const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // intent(in)
        , int32_t const (*const *const __restrict__ index_list) // index list that brings the blocks in order,
                                                                // list must contain at least 4+1+4 elements
        , double const prefactor
        , int const Stride // Stride is determined by the lattice dimension along which we derive: 1, 4 or 4^2
        , double const phase[2][2]=nullptr // WARNING this routine cannot treat periodic boundary conditions
    ) {
        assert(nullptr == phase);          // WARNING this routine cannot treat periodic boundary conditions
        // prepare finite-difference coefficients
        //            c_0        c_1         c_2       c_3        c_4      c_5       c_6     c_7     c_8
        // FD16th = [-924708642, 538137600, -94174080, 22830080, -5350800, 1053696, -156800, 15360, -735] / 302702400
        double const norm = prefactor/302702400.;
        real_t const  c0 = -924708642*norm,
                      c1 =  538137600*norm,
                      c2 =  -94174080*norm,
                      c3 =   22830080*norm,
                      c4 =   -5350800*norm,
                      c5 =    1053696*norm,
                      c6 =    -156800*norm,
                      c7 =      15360*norm,
                      c8 =       -735*norm;

        assert(16 == gridDim.y);
        assert(1  == gridDim.z);
        assert(Noco*64 == blockDim.x);
        assert(Noco    == blockDim.y);
        assert(R1C2    == blockDim.z);

#ifdef    HAS_NO_CUDA
        dim3 blockIdx(0,0,0);
        for (blockIdx.y = 0; blockIdx.y < gridDim.y; ++blockIdx.y)
        for (blockIdx.x = 0; blockIdx.x < gridDim.x; ++blockIdx.x)
#endif // HAS_NO_CUDA
        { // block loops

        auto const *const list = index_list[blockIdx.x]; // abbreviate pointer

        int const i16 = blockIdx.y; // in [0, 16)
        int const i64 = (16==Stride)? i16 : ( (4==Stride)? (16*(i16 >> 2) + (i16 & 0x3)) : (4*i16) );

#ifdef    HAS_NO_CUDA
        dim3 threadIdx(0,0,0);
        for (threadIdx.z = 0; threadIdx.z < blockDim.z; ++threadIdx.z)
        for (threadIdx.y = 0; threadIdx.y < blockDim.y; ++threadIdx.y)
        for (threadIdx.x = 0; threadIdx.x < blockDim.x; ++threadIdx.x)
#endif // HAS_NO_CUDA
        { // thread loops

#define INDICES(i4) [threadIdx.z][threadIdx.y*64 + i64 + Stride*i4][threadIdx.x]

        real_t w0{0}, w1{0}, w2{0}, w3{0}, w4{0}, w5{0}, w6{0}, w7{0}, // initialize two non-existing blocks (isolated boundary condition)
               w8, w9, wa, wb, wc, wd, we, wf, wn; // 8 + 8 + 1 registers

        int ilist{0}; // counter for index_list
        for (int i4 = 0; i4 < nhalo; ++i4) {
            assert(0 == list[ilist++] && "Laplace16th can only deal with isolated boundary conditions");
        } // i4

        // initially load two blocks in advance
        int i0 = list[ilist++] - BLOCK_EXISTS; // load index for 1st non-zero block
        assert(i0 >= 0); // 1st central block must exist
        w8 = psi[i0]INDICES(0); // inital load
        w9 = psi[i0]INDICES(1); // inital load
        wa = psi[i0]INDICES(2); // inital load
        wb = psi[i0]INDICES(3); // inital load

        int i1 = list[ilist++] - BLOCK_EXISTS; // load index for the 2nd block
        if (i1 >= 0) {
            wc = psi[i1]INDICES(0); // inital load
            wd = psi[i1]INDICES(1); // inital load
            we = psi[i1]INDICES(2); // inital load
            wf = psi[i1]INDICES(3); // inital load
        } else {
            wc = 0; wd = 0; we = 0; wf = 0; // second block is already non-existing
        } // i1 valid

        // main loop
        int i2 = list[ilist++] - BLOCK_EXISTS; // get next index

        assert(ilist == 7);

        while (i0 >= 0) {
            bool const load = (i2 >= 0); // load or not?
            // use a rotating register file, compare version of Laplace8th above

            // FD17POINT = load?, compute, store, update rotating register file
#define FD17POINT(i4,  M8, M7, M6, M5, M4, M3, M2, M1, W0, P1, P2, P3, P4, P5, P6, P7, P8) \
            P8 = load ? psi[i2]INDICES(i4) : 0; \
            Tpsi[i0]INDICES(i4) += c0*W0 + c1*M1 + c1*P1 + c2*M2 + c2*P2 + c3*M3 + c3*P3 + c4*M4 + c4*P4 \
                                         + c5*M5 + c5*P5 + c6*M6 + c6*P6 + c7*M7 + c7*P7 + c8*M8 + c8*P8; \
            M8 = P8;

            int const mod4 = ilist & 0x3; // binary modulo 4
            // ToDo: use a switch statement or nested if branches
            if        (0x3 == mod4) { // 3
                FD17POINT(0,  w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, wa, wb, wc, wd, we, wf, wn)
                FD17POINT(1,  w1, w2, w3, w4, w5, w6, w7, w8, w9, wa, wb, wc, wd, we, wf, w0, wn)
                FD17POINT(2,  w2, w3, w4, w5, w6, w7, w8, w9, wa, wb, wc, wd, we, wf, w0, w1, wn)
                FD17POINT(3,  w3, w4, w5, w6, w7, w8, w9, wa, wb, wc, wd, we, wf, w0, w1, w2, wn)
            } else if (0x0 == mod4) { // 0
                FD17POINT(0,  w4, w5, w6, w7, w8, w9, wa, wb, wc, wd, we, wf, w0, w1, w2, w3, wn)
                FD17POINT(1,  w5, w6, w7, w8, w9, wa, wb, wc, wd, we, wf, w0, w1, w2, w3, w4, wn)
                FD17POINT(2,  w6, w7, w8, w9, wa, wb, wc, wd, we, wf, w0, w1, w2, w3, w4, w5, wn)
                FD17POINT(3,  w7, w8, w9, wa, wb, wc, wd, we, wf, w0, w1, w2, w3, w4, w5, w6, wn)
            } else if (0x1 == mod4) { // 1
                FD17POINT(0,  w8, w9, wa, wb, wc, wd, we, wf, w0, w1, w2, w3, w4, w5, w6, w7, wn)
                FD17POINT(1,  w9, wa, wb, wc, wd, we, wf, w0, w1, w2, w3, w4, w5, w6, w7, w8, wn)
                FD17POINT(2,  wa, wb, wc, wd, we, wf, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, wn)
                FD17POINT(3,  wb, wc, wd, we, wf, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, wa, wn)
            } else                   { // 2
                FD17POINT(0,  wc, wd, we, wf, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, wa, wb, wn)
                FD17POINT(1,  wd, we, wf, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, wa, wb, wc, wn)
                FD17POINT(2,  we, wf, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, wa, wb, wc, wd, wn)
                FD17POINT(3,  wf, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, wa, wb, wc, wd, we, wn)
            } // ilist mod(4)
#undef  FD17POINT
            i0 = i1; i1 = i2; i2 = list[ilist++] - BLOCK_EXISTS; // rotate and get next index
        } // while loop

#undef  INDICES
        }} // thread and block loops

    } // Laplace16th



//  NVIDIA V100 GPU:
//   Maximum number of threads per multiprocessor:  2048
//   Maximum number of threads per block:           1024
//   Max dimension size of a thread block (x,y,z): (1024, 1024, 64)
//   Max dimension size of a grid size    (x,y,z): (2147483647, 65535, 65535)
//   Maximum memory pitch:                          2147483647 bytes



    template <typename real_t, int R1C2=2, int Noco=1>
    int Laplace_driver(
          real_t        (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // intent(inout)
        , real_t  const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // intent(in)
        , int32_t const (*const *const __restrict__ index_list) // index list that brings the blocks in order,
                                            // list must contain at least one element and is finalized with -1
        , double const prefactor
        , uint32_t const num
        , int const Stride // Stride is determined by the lattice dimension along which we derive: 1, 4 or 4^2
        , double const phase[2][2]=nullptr
        , int const FD_range=4 // 4 or 8 are implemented
    ) {
        if (num < 1 || FD_range < 1) return 0;
        assert(1 == Stride || 4 == Stride || 16 == Stride);
        auto const kernel_ptr = (8 == FD_range) ? Laplace16th<real_t,R1C2,Noco> : Laplace8th<real_t,R1C2,Noco>;
        if (8 == FD_range) phase = nullptr; //    Laplace16th cannot handle Bloch phases
        dim3 const gridDim(num, 16, 1), blockDim(Noco*64, Noco, R1C2);
        kernel_ptr // GPU kernel, must be launched with <<< {16, Nrows, 1}, {Noco*64, Noco, R1C2} >>>
#ifdef    HAS_NO_CUDA
                  (    gridDim, blockDim,
#else  // HAS_NO_CUDA
                   <<< gridDim, blockDim >>> (
#endif // HAS_NO_CUDA
                   Tpsi, psi, index_list, prefactor, Stride, phase);
        return (8 == FD_range) ? 17 : 9; // returns the number of stencil coefficients
    } // Laplace_driver


    void __host__ set_phase(
          double phase[3][2][2]
        , double const phase_angles[3]=nullptr
        , int const echo=0
    ); // declaration only

    // template <typename real_t, int R1C2=2, int Noco=1>
    // size_t multiply(
    //       real_t         (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // result
    //     , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input
    //     , uint32_t const num[3] // number of sticks in X,Y,Z direction
    //     , int32_t  const *const *const __restrict__ x_list //
    //     , int32_t  const *const *const __restrict__ y_list //
    //     , int32_t  const *const *const __restrict__ z_list //
    //     , double   const hgrid[3] // grid spacing in X,Y,Z
    //     , int      const FD_range[3] // finite-difference stencil range (in grid points)
    //     , double   const phase[3][2][2] // complex Bloch phase factors
    //     , size_t   const nnzb=1 // total number of non-zero blocks (to get the operations count correct)
    //     , int      const echo=0
    // ) {
    //     int nFD[3]; // numbers of stencil coefficients
    //     size_t nops{0};
    //     for (int dd = 0; dd < 3; ++dd) { // derivative direction, serial due to update of Tpsi
    //         auto const f = -0.5/pow2(hgrid[dd]); // prefactor of the kinetic energy in Hartree atomic units
    //         int  const stride = 1 << (2*dd);
    //         auto const list = (2 == dd) ? z_list : ((1 == dd) ? y_list : x_list);
    //         nFD[dd] = Laplace_driver<real_t,R1C2,Noco>(Tpsi, psi, list, f, num[dd], stride, phase[dd], FD_range[dd]);
    //         nops += nnzb*nFD[dd]*R1C2*pow2(Noco*64ul)*2ul; // total number of floating point operations performed
    //     } // dd
    //     char const fF = (8 == sizeof(real_t)) ? 'F' : 'f'; // Mflop:float, MFlop:double
    //     if (echo > 7) std::printf("# green_kinetic::%s nFD= %d %d %d, numbers= %d %d %d, %.3f M%clop\n",
    //                           __func__, nFD[0], nFD[1], nFD[2], num[0], num[1], num[2], nops*1e-6, fF);
    //     return nops;
    // } // multiply (kinetic energy operator)


    status_t finite_difference_plan(
            green_sparse::sparse_t<int32_t> & sparse // result: derivation list of block indices
          , int16_t & FD_range // side result: recommended finite-difference range
          , int const dd // direction of derivative
          , bool const boundary_is_periodic
          , uint32_t const num_target_coords[3]
          , uint32_t const RowStart[]
          , uint16_t const ColIndex[]
          , view3D<int32_t> const & iRow_of_coords // (Z,Y,X) look-up table: row index of the Green function as a function of internal 3D coordinates, -1:non-existent
          , std::vector<bool> const sparsity_pattern[] // memory saving bit-arrays sparsity_pattern[irhs][idx3]
          , unsigned const nrhs=1 // number of right hand sides
          , int const echo=0 // log level
    ); // declaration only




    class kinetic_plan_t {
    public:

        kinetic_plan_t() {} // default constructor
#if 0
        kinetic_plan_t(
            green_sparse::sparse_t<int32_t> & sparse // result
          , int16_t & FD_range // side result
          , int const dd // direction of derivative, 0:X, 1:Y, 2:Z
          , bool const boundary_is_periodic
          , uint32_t const num_target_coords[3]
          , uint32_t const RowStart[]
          , uint16_t const ColIndex[]
          , view3D<int32_t> const & iRow_of_coords // (Z,Y,X) look-up table: row index of the Green function as a function of internal 3D coordinates, -1:non-existent
          , std::vector<bool> const sparsity_pattern[] // memory saving bit-arrays sparsity_pattern[irhs][idx3]
          , unsigned const nrhs=1 // number of right hand sides
          , double const grid_spacing=1 // grid spacing in derivative direction
          , int const echo=0 // log level
      ) 
        : derivative_direction(dd)
      {
          auto const stat = finite_difference_plan(sparse, FD_range, // results
                    dd, boundary_is_periodic, num_target_coords, RowStart, ColIndex,
                    iRow_of_coords, sparsity_pattern, nrhs, echo);
          if (0 == stat) {
              prefactor = -0.5/(grid_spacing*grid_spacing);
              lists = get_memory<int32_t const *>(sparse.nRows(), echo, "lists[dd]");
          } else {
              if (echo > 2) std::printf("# failed to set up finite_difference_plan for %c-direction, status= %i\n", 'x' + dd, int(stat));
          }
      } // constructor
#endif // does not work because the copy assignment operator of sparse is deleted

      void set(
            int const dd
          , double const grid_spacing=1
          , size_t const nnzbX=1
          , int const echo=0 // log level
      ) {
            derivative_direction = dd;
            nnzb = nnzbX;
            prefactor = -0.5/(grid_spacing*grid_spacing);
            char lists_[8] = "list[?]"; lists_[5] = 'x' + dd;
            lists = get_memory<int32_t const *>(sparse.nRows(), echo, lists_);
            auto const rowStart = sparse.rowStart();
            auto const colIndex = sparse.colIndex();
            for (uint32_t il = 0; il < sparse.nRows(); ++il) {
                lists[il] = &colIndex[rowStart[il]];
            } // il
      } // set

      ~kinetic_plan_t() {
          if (lists) {
              std::printf("# free list for dd=%i\n", derivative_direction);
              free_memory(lists);
          }
      } // destructor

    public:
        // members
        green_sparse::sparse_t<int32_t> sparse;
        double prefactor = 0; // = -0.5/h^2,  h:grid spacing in X,Y,Z // prefactor of the kinetic energy in Hartree atomic units
        size_t nnzb; // total number of non-zero blocks (to get the operations count correct)
        int16_t FD_range = 8; // 4 or 8
        int32_t const ** lists = nullptr; // in device memory
        int derivative_direction = -1; // derivative direction

    public:

        template <typename real_t, int R1C2=2, int Noco=1>
        size_t multiply(
              real_t         (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // result
            , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input
            , double   const phase[2][2]=nullptr // complex Bloch phase factors
            , int      const echo=0
        ) const { // members of the kinetic_plan_t are not changed
            int  const stride = 1 << (2*derivative_direction); // 4^dd: X:1, Y:4, Z:16
            auto const nFD = Laplace_driver<real_t,R1C2,Noco>(Tpsi, psi, lists, prefactor, sparse.nRows(), stride, phase, FD_range);
            size_t const nops = nnzb*nFD*R1C2*pow2(Noco*64ul)*2ul;
            if (echo > 7) {
                char const fF = (8 == sizeof(real_t)) ? 'F' : 'f'; // Mflop:float, MFlop:double
                std::printf("# green_kinetic::%s dd=\'%c\', nFD= %d, number= %d, %.3f M%clop\n",
                    __func__, 'x' + derivative_direction, nFD, sparse.nRows(), nops*1e-6, fF);
            } // echo
            return nops; // return the total number of floating point operations performed
        } // multiply (kinetic energy operator)

    }; // class kinetic_plan_t



#if 0
    template <typename real_t, int R1C2=2, int Noco=1>
    size_t multiply( // ToDo: delete this and see if it is still in use
          real_t         (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // result
        , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input
        , green_sparse::sparse_t<int32_t> const kinetic_plan[3]
        , double const hgrid[3] // grid spacing in X,Y,Z
        , double const phase[3][2][2] // complex Bloch phase factors [direction][forward:backward][real:imag]
        , int16_t const FD_range[3] // finite-difference stencil range (in grid points)
        , size_t const nnzb=1 // total number of non-zero blocks (to get the operations count correct)
        , int const echo=0
    ) {
        int const nFD[] = {FD_range[0], FD_range[1], FD_range[2]};
        uint32_t num[3];
        int32_t const ** lists[3];
        for (int dd = 0; dd < 3; ++dd) {
            // convert fd_plan to arrays of GPU pointers
            num[dd] = kinetic_plan[dd].nRows();
            lists[dd] = get_memory<int32_t const *>(num[dd], echo, "num[dd]");
            auto const rowStart = kinetic_plan[dd].rowStart();
            auto const colIndex = kinetic_plan[dd].colIndex();
            for (uint32_t il = 0; il < num[dd]; ++il) {
                lists[dd][il] = &colIndex[rowStart[il]];
            } // il
        } // dd
        if (echo > 13) std::printf("# green_kinetic::%s FD_range= %d %d %d, numbers= %d %d %d\n",
                   __func__, FD_range[0], FD_range[1], FD_range[2], num[0], num[1], num[2]);

        auto const nops = multiply<real_t,R1C2,Noco>(Tpsi, psi, num, lists[0], lists[1], lists[2], hgrid, nFD, phase, nnzb, echo);

        for (int dd = 0; dd < 3; ++dd) {
            free_memory(lists[dd]);
        } // dd
        return nops;
    } // multiply (kinetic energy operator)
#endif

    status_t all_tests(int const echo=0); // declaration only


} // namespace green_kinetic
