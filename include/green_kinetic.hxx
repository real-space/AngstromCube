#pragma once

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cmath> // std::sqrt
#include <algorithm> // std::max
#include <utility> // std::swap //, ::move
#include <vector> // std::vector<T>

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_memory.hxx" // get_memory, free_memory
#include "green_sparse.hxx" // ::sparse_t<>
#include "data_view.hxx" // view3D<T>
#include "inline_math.hxx" // set
#include "simple_stats.hxx" // ::Stats<>
#include "print_tools.hxx" // printf_vector(format, ptr, number [, ...])
#include "constants.hxx" // ::pi

  // ToDo: move to green_utils.hxx or similar
  template <typename uint_t, typename int_t> inline
  size_t index3D(uint_t const n[3], int_t const i[3]) {
      // usual 3D indexing
      return size_t(i[2]*n[1] + i[1])*n[0] + i[0];
  } // index3D


namespace green_kinetic {

    int constexpr nhalo = 4; // a maximum of 4 bocks (i.e. 16 grid points) is the range of the FD stencil.
  
    int32_t get_inz(
            uint32_t const idx[3]
          , view3D<int32_t> const & iRow_of_coords // (Z,Y,X) look-up table: row index of the Green function as a function of internal 3D coordinates, -1:non-existent
          , uint32_t const RowStart[]
          , uint16_t const ColIndex[]
          , unsigned const irhs
          , char const dd='?' // derivative direction
    ) {
          int constexpr X=0, Y=1, Z=2;
                          auto const iRow = iRow_of_coords(idx[Z], idx[Y], idx[X]);
                          assert(iRow >= 0 && "sparsity_pattern[irhs][idx3] does not match iRow_of_coords[iz][iy][ix]");

                          int32_t inz_found{-1};
                          for (auto inz = RowStart[iRow]; inz < RowStart[iRow + 1]; ++inz) {
                              if (ColIndex[inz] == irhs) {
                                  inz_found = inz; // store where it was found
                                  inz = RowStart[iRow + 1]; // stop search loop
                              } // found
                          } // search
                          return inz_found;
    } // get_inz
  
  
    inline status_t finite_difference_plan(
            green_sparse::sparse_t<int32_t> & sparse // result
          , int const dd // direction of derivative
          , bool const boundary_is_periodic
          , uint32_t const num_target_coords[3]
          , uint32_t const RowStart[]
          , uint16_t const ColIndex[]
          , view3D<int32_t> const & iRow_of_coords // (Z,Y,X) look-up table: row index of the Green function as a function of internal 3D coordinates, -1:non-existent
          , std::vector<bool> const sparsity_pattern[] // memory saving bit-arrays sparsity_pattern[irhs][idx3]
          , unsigned const nrhs=1 // number of right hand sides
          , int const echo=0
      )
        // Preparation of Finite-Difference index lists
        // 2D example: non-zero index -1 means non-existent
        // 
        //                            --> x-direction
        //        0  1  2  3  4
        //     5  6  7  8  9 10 11        |
        //    12 13 14 15 16 17 18 19     |
        //    20 21 22 23 24 25 26 27     v
        //    28 29 30 31 32 33 34        y-direction
        //       35 36 37 38 39
        //
        //
        //  6 x-lists:
        //    list[0] == { 0  0  0  0  1  2  3  4  5  0  0  0  0}
        //    list[1] == { 0  0  0  0  6  7  8  9 10 11 12  0  0  0  0}
        //    list[2] == { 0  0  0  0 13 14 15 16 17 18 19 20  0  0  0  0}
        //    list[3] == { 0  0  0  0 21 22 23 24 25 26 27 28  0  0  0  0}
        //    list[4] == { 0  0  0  0 29 30 31 32 33 34 35  0  0  0  0}
        //    list[5] == { 0  0  0  0 36 37 38 39 40  0  0  0  0}
        //
        //  8 y-lists:
        //    list[0] == { 0  0  0  0  6 13 21 29  0  0  0  0}
        //    list[1] == { 0  0  0  0  1  7 14 22 30 36  0  0  0  0}
        //    list[2] == { 0  0  0  0  2  8 15 23 31 37  0  0  0  0}
        //    list[3] == { 0  0  0  0  3  9 16 24 32 38  0  0  0  0}
        //    list[4] == { 0  0  0  0  4 10 17 25 33 39  0  0  0  0}
        //    list[5] == { 0  0  0  0  5 11 18 26 34 40  0  0  0  0}
        //    list[6] == { 0  0  0  0 12 19 27 35  0  0  0  0}
        //    list[7] == { 0  0  0  0 20 28  0  0  0  0}
        //
        // Preparation of Finite-Difference index lists
        // 2D example with a periodic x-direction
        // 
        //                            --> x-direction
        //        0  1  2  3  4
        //        5  6  7  8  9
        //
        //  2 x-lists:
        //    list[0] == {-2 -3 -4 -5   1  2  3  4  5  -1 -2 -3 -4}
        //    list[1] == {-7 -8 -9 -10  6  7  8  9 10  -6 -7 -8 -9}
        //
      {
          int constexpr X=0, Y=1, Z=2;
          // prepare the finite-difference sequence lists
          char const direction = 'x' + dd;
          assert(X == dd || Y == dd || Z == dd); 
          assert(num_target_coords[X] > 0); assert(num_target_coords[Y] > 0); assert(num_target_coords[Z] > 0);
          int num[3];
          set(num, 3, num_target_coords);
          uint32_t const num_dd = num[dd];
          num[dd] = 1; // replace number of target blocks in derivative direction
          if (echo > 0) std::printf("# FD lists in %c-direction %d %d %d\n", direction, num[X], num[Y], num[Z]);
          simple_stats::Stats<> length_stats;
          size_t const max_lists = nrhs*size_t(num[Z])*size_t(num[Y])*size_t(num[X]);
          std::vector<std::vector<int32_t>> list(max_lists);
          size_t ilist{0};
          for (unsigned irhs = 0; irhs < nrhs; ++irhs) {
//                if (echo > 0) std::printf("# FD list for RHS #%i\n", irhs);
              auto const & sparsity_rhs = sparsity_pattern[irhs];
              for (uint32_t iz = 0; iz < num[Z]; ++iz) { //  
              for (uint32_t iy = 0; iy < num[Y]; ++iy) { //   one of these 3 loops has range == 1
              for (uint32_t ix = 0; ix < num[X]; ++ix) { // 
                  uint32_t idx[3] = {ix, iy, iz}; // non-const
                  assert(ilist < max_lists);
                  assert(0 == list[ilist].size()); // make sure the list is empty at start
                  int32_t last_id{-1};
                  for (uint32_t id = 0; id < num_dd; ++id) { // loop over direction to derive
                      idx[dd] = id; // replace index in the derivate direction
//                    if (echo > 0) std::printf("# FD list for RHS #%i test coordinates %i %i %i\n", irhs, idx[X], idx[Y], idx[Z]);
                      auto const idx3 = index3D(num_target_coords, idx); // flat index
                      if (sparsity_rhs[idx3]) {
                          auto const inz_found = get_inz(idx, iRow_of_coords, RowStart, ColIndex, irhs, 'x' + dd);
                          assert(inz_found >= 0 && "sparsity_pattern[irhs][idx3] inconsistent with iRow_of_coords[iz][iy][ix]");

                          if (0 == list[ilist].size()) {
                              list[ilist].reserve(nhalo + num_dd + nhalo); // makes push_back operation faster
                              list[ilist].resize(nhalo, 0); // prepend {0, 0, 0, 0}
                              // =====================================================================================================
                              // =====================================================================================================
                              if (boundary_is_periodic) {
                                  if (0 == id) { // only when we are at the leftmost boundary
                                      for(int ihalo = 0; ihalo < nhalo; ++ihalo) {
                                          uint32_t jdx[] = {idx[X], idx[Y], idx[Z]};
                                          jdx[dd] = (nhalo*num_target_coords[dd] - nhalo + ihalo) % num_target_coords[dd]; // periodic wrap around
                                          auto const jdx3 = index3D(num_target_coords, jdx); // flat index
                                          if (sparsity_rhs[jdx3]) {
                                              auto const jnz_found = get_inz(jdx, iRow_of_coords, RowStart, ColIndex, irhs, 'X' + dd);
                                              if (jnz_found >= 0) {
                                                  list[ilist][0 + ihalo] = -(jnz_found + 1); // negative index marks that a (left) Bloch phase needs to be multiplied
                                              } // block exists
                                          } // block exists
                                      } // ihalo
                                  } // leftmost boundary
                              } // boundary_is_periodic
                              // =====================================================================================================
                              // =====================================================================================================
                          } // list was empty

                          list[ilist].push_back(inz_found + 1); // +1 so that non-exiting elements are marked by 0
                          last_id = id;
                      } // sparsity pattern
                  } // id
                  int const list_length = list[ilist].size();
                  if (list_length > nhalo) {
                      length_stats.add(list_length - nhalo);
//                    if (echo > 0) std::printf("# FD list of length %d for the %c-direction %i %i %i\n", list_length, direction, idx[X], idx[Y], idx[Z]);
                      // add nhalo end-of-sequence markers
                      for (int ihalo = 0; ihalo < nhalo; ++ihalo) {
                          list[ilist].push_back(0); // append {0, 0, 0, 0} to mark the end of the derivative sequence
                      } // ihalo
                      // =====================================================================================================
                      // =====================================================================================================
                      if (boundary_is_periodic) {
                          if (num_target_coords[dd] - 1 == last_id) { // only when the last entry was at the rightmost boundary
                                      for(int ihalo = 0; ihalo < nhalo; ++ihalo) {
                                          uint32_t jdx[] = {idx[X], idx[Y], idx[Z]};
                                          jdx[dd] = ihalo % num_target_coords[dd]; // periodic wrap around
                                          auto const jdx3 = index3D(num_target_coords, jdx); // flat index
                                          if (sparsity_rhs[jdx3]) {
                                              auto const jnz_found = get_inz(jdx, iRow_of_coords, RowStart, ColIndex, irhs, 'X' + dd);
                                              if (jnz_found >= 0) {
                                                  list[ilist][list_length + ihalo] = -(jnz_found + 1); // negative index marks that a (right) Bloch phase needs to be multiplied
                                              } // block exists
                                          } // block exists
                                      } // ihalo
                          } // rightmost boundary
                      } // boundary_is_periodic
                      // =====================================================================================================
                      // =====================================================================================================
                      assert(list_length + 4 == list[ilist].size());

                      ++ilist; // create next list index
                  } // list_length > 0
              }}} // ixyz
          } // irhs

          // store the number of lists
          uint32_t const n_lists = ilist; assert(n_lists == ilist && "too many lists, max. 2^32-1");
          if (echo > 0) std::printf("# %d FD lists for the %c-direction (%.2f %%), length %.3f +/- %.3f, min %g max %g\n",
                                n_lists, direction, n_lists/(max_lists*.01),
                                length_stats.mean(), length_stats.dev(), length_stats.min(), length_stats.max());
          list.resize(n_lists);

          sparse = green_sparse::sparse_t<int32_t>(list, false, "finite_difference_list", echo);
          return 0;
      } // finite_difference_plan
  

    template <typename real_t, int R1C2=2, int Noco=1> // Stride is determined by the lattice dimension along which we derive
    void __global__ Laplace8th( // GPU kernel, must be launched with <<< {Nrows, 16, 1}, {Noco*64, Noco, R1C2} >>>
#ifdef HAS_NO_CUDA
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

#ifdef HAS_NO_CUDA
        dim3 blockIdx(0,0,0);
        for (blockIdx.y = 0; blockIdx.y < gridDim.y; ++blockIdx.y)
        for (blockIdx.x = 0; blockIdx.x < gridDim.x; ++blockIdx.x)
#endif // HAS_NO_CUDA
        { // block loops

        auto const *const list = index_list[blockIdx.x]; // abbreviate pointer

        int const i16 = blockIdx.y;
        int const i64 = (16==Stride)? i16 : ( (4==Stride)? (16*(i16 >> 2) + (i16 & 0x3)) : (4*i16) );

#ifdef HAS_NO_CUDA
        dim3 threadIdx(0,0,0);
        for (threadIdx.z = 0; threadIdx.z < blockDim.z; ++threadIdx.z)
        for (threadIdx.y = 0; threadIdx.y < blockDim.y; ++threadIdx.y)
        for (threadIdx.x = 0; threadIdx.x < blockDim.x; ++threadIdx.x)
#endif // HAS_NO_CUDA
        { // thread loops

#define INDICES(i4) [threadIdx.z][threadIdx.y*64 + i64 + Stride*i4][threadIdx.x]

        real_t w0{0}, w1{0}, w2{0}, w3{0}; // 4 registers for one block

        int ilist{3}; // counter for index_list

        int ii = list[ilist++]; // the 1st block can be 0 (non-existent) or <0 for a periodic image
        // =========================================================================================
        // === periodic boundary conditions ========================================================
        if (ii) {
            if (ii > 0) std::printf("# Error: ii= %d ilist= %i\n", ii, ilist);
            assert(ii < 0 && "list[3] must be either 0 (isolated BC) or negative (periodic BC)");
            int const jj = -ii - 1; // index of the periodic image of a block
            assert(jj >= 0);
            assert(phase && "no phase given");
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
        } else { // is periodic
            // block does not exist (isolated/vacuum boundary condition), registers w0...w3 are initialized 0
        } // is periodic
        // === periodic boundary conditions ========================================================
        // =========================================================================================

        real_t w4, w5, w6, w7, wn; // 4 + 1 registers, wn is the register that always receives the most recently loaded value

        ii = list[ilist++] - 1; assert(ii >= 0); // this block must be a regular index since it is the 1st block for which we store a result
        // initially load one block in advance
        w4 = psi[ii]INDICES(0);
        w5 = psi[ii]INDICES(1);
        w6 = psi[ii]INDICES(2);
        w7 = psi[ii]INDICES(3);

        // main loop
        assert(5 == ilist);
        while (ii >= 0) {
            int const i0 = ii; // set central index 
            ii = list[ilist++] - 1; // get next index
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
        if (ii) {
            assert(ii < 0 && "last list item must be either 0 (isolated BC) or negative (periodic BC)");
            int const jj = -ii - 1; // index of the periodic image of a block
            assert(jj >= 0);
            int const i0 = list[ilist - 2] - 1; // recover the last central index
            assert(i0 >= 0);
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
#ifdef HAS_NO_CUDA
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
        // prepare finite-difference coefficients
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

#ifdef HAS_NO_CUDA
        dim3 blockIdx(0,0,0);
        for (blockIdx.y = 0; blockIdx.y < gridDim.y; ++blockIdx.y)
        for (blockIdx.x = 0; blockIdx.x < gridDim.x; ++blockIdx.x)
#endif // HAS_NO_CUDA
        { // block loops

        auto const *const list = index_list[blockIdx.x]; // abbreviate pointer

        int const i16 = blockIdx.y;
        int const i64 = (16==Stride)? i16 : ( (4==Stride)? (16*(i16 >> 2) + (i16 & 0x3)) : (4*i16) );

#ifdef HAS_NO_CUDA
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
        for (int i4 = 0; i4 < 4; ++i4) {
            assert(0 == list[ilist++] && "Laplace16th can only deal with isolated boundary conditions");
        } // i4

        // initially load two blocks in advance
        int i0 = list[ilist++] - 1; // load index for 1st non-zero block
        assert(i0 >= 0); // 1st central block must exist
        w8 = psi[i0]INDICES(0); // inital load
        w9 = psi[i0]INDICES(1); // inital load
        wa = psi[i0]INDICES(2); // inital load
        wb = psi[i0]INDICES(3); // inital load

        int i1 = list[ilist++] - 1; // load index for the 2nd block
        if (i1 >= 0) {
            wc = psi[i1]INDICES(0); // inital load
            wd = psi[i1]INDICES(1); // inital load
            we = psi[i1]INDICES(2); // inital load
            wf = psi[i1]INDICES(3); // inital load
        } else {
            wc = 0; wd = 0; we = 0; wf = 0; // second block is already non-existing
        } // i1 valid

        // main loop
        int i2 = list[ilist++] - 1; // get next index

        assert(ilist == 7);

        while (i0 >= 0) {
            bool const load = (i2 >= 0); // load or not?
            // use a rotating register file, see above version of Laplace8th

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
            i0 = i1; i1 = i2; i2 = list[ilist++] - 1; // rotate and get next index
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
        , int const nFD=4
    ) {
        assert(1 == Stride || 4 == Stride || 16 == Stride);
        auto const kernel_ptr = (8 == nFD) ? Laplace16th<real_t,R1C2,Noco> : Laplace8th<real_t,R1C2,Noco>;
        dim3 const gridDim(num, 16, 1), blockDim(Noco*64, Noco, R1C2);
        kernel_ptr // GPU kernel, must be launched with <<< {16, Nrows, 1}, {Noco*64, Noco, R1C2} >>>
#ifdef HAS_NO_CUDA
                  (    gridDim, blockDim,
#else  // HAS_NO_CUDA
                   <<< gridDim, blockDim >>> (
#endif // HAS_NO_CUDA                     
                   Tpsi, psi, index_list, prefactor, Stride, phase);
        return (8 == nFD) ? 8 : 4;
    } // Laplace_driver


    template <typename real_t=double>
    double __host__ set_phase(real_t phase[2][2], double const phase_angle=0, char const direction='?', int const echo=0) {
        // phase_angle comes in units of 2*pi
        // phase_angle ==     0.00 --> phase = 1
        // phase_angle == +/- 0.25 --> phase = +/-i
        // phase_angle == +/- 0.50 --> phase = -1
        // phase_angle == +/- 0.75 --> phase = -/+i
        assert(std::abs(phase_angle) < 0.503); // 181 degrees are accepted for test purposes

        double dev{0};
        
        auto const arg = 2*constants::pi*phase_angle;
        auto cs = std::cos(arg), sn = std::sin(arg);

//         auto const ph = std::pow(std::complex<double>(0, 1), 4*phase_angle);
//         auto cs = ph.real(), sn = ph.imag();

        // purify the values for multiples of 15 degrees
        int const phase24 = phase_angle*24;
        char const *corrected = "";
        if (phase_angle*24 == phase24) { // exactly on the grid of 15 degree steps
            int const phase4 = (phase24 + 120) /  6;
            int const phase2 = (phase24 + 120) % 12;
            double const sh = std::sqrt(0.5), s3 = std::sqrt(3);
            double const abs_values6[7] = {0.0, (s3 - 1)*sh*.5, .5, sh, s3*.5, (s3 + 1)*sh*.5, 1.0}; // == sin(i*15degrees)

            int64_t const sin_values = 0x0123456543210; // look-up table with 12 entries (4bit wide each)
            int const i6_sin = (sin_values >> (4*phase2)) & 0x7;
            int const i6_cos = 6 - i6_sin;
//          std::printf("angle= %g degrees = %d * 15 degrees, ph4=%d ph2=%d i6_cos=%i i6_sin=%i\n", phase_angle*360, phase24, phase4, phase2, i6_cos, i6_sin);

            int const sgn_sin = 1 - ( phase4      & 0x2), // {1, 1, -1, -1}[phase4 % 4]
                      sgn_cos = 1 - ((phase4 + 1) & 0x2); // {1, -1, -1, 1}[phase4 % 4]
            sn = sgn_sin * abs_values6[i6_sin];
            cs = sgn_cos * abs_values6[i6_cos];
            corrected = " corrected";
        }
// #ifdef  DEBUG
        if ('\0' != *corrected) {
            if (std::abs(cs - std::cos(arg)) > 2.3e-16) error("cosine for phase_angle= %g degrees deviates after purification, cos= %g expected %g", phase_angle*360, cs, std::cos(arg));
            if (std::abs(sn - std::sin(arg)) > 2.8e-16) error(  "sine for phase_angle= %g degrees deviates after purification, sin= %g expected %g", phase_angle*360, sn, std::sin(arg));
            dev = std::max(std::abs(cs - std::cos(arg)), std::abs(sn - std::sin(arg)));
        }
// #endif // DEBUG
        
        int constexpr Re=0, Im=1, Left=0, Right=1;
        phase[Left][Re] = cs;
        phase[Left][Im] = sn;
        if (echo > 9) std::printf("# %s angle in %c-direction is %g degree, (%g, %g)%s\n",
                              __func__, direction, phase_angle*360,  cs, sn, corrected);
        phase[Right][Re] =  cs;
        phase[Right][Im] = -sn; // right phase factor is conjugate of the left

        return dev;
    } // set_phase

    template <typename real_t=double>
    void __host__ set_phase(real_t phase[3][2][2], double const phase_angle[3]=nullptr, int const echo=0) {
        for (int dd = 0; dd < 3; ++dd) {
            set_phase(phase[dd], phase_angle ? phase_angle[dd] : 0, 'x' + dd, echo);
        } // dd
    } // set_phase
    
    
    template <typename real_t, int R1C2=2, int Noco=1>
    size_t multiply(
          real_t         (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // result
        , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // 
        , uint32_t const num[3] // number of sticks in X,Y,Z direction
        , int32_t  const *const *const __restrict__ x_list // 
        , int32_t  const *const *const __restrict__ y_list // 
        , int32_t  const *const *const __restrict__ z_list // 
        , double   const hgrid[3] // grid spacing in X,Y,Z
        , int      const FD_range[3] // finite-difference stencil range (in grid points)
        , double   const phase[3][2][2] // complex phase factors
        , size_t   const nnzb=1 // total number of non-zero blocks (to get the operations count correct)
        , int      const echo=0
    ) {
        int nFD[] = {FD_range[0], FD_range[1], FD_range[2]};
        size_t nops{0};
        for (int dd = 0; dd < 3; ++dd) { // derivative direction, serial due to update of Tpsi
            double const f = -0.5/(hgrid[dd]*hgrid[dd]); // prefactor of the kinetic energy in Hartree atomic units
            int const Stride = 1 << (2*dd);
            auto const list = (2 == dd) ? z_list : ((1 == dd) ? y_list : x_list);
            nFD[dd] = Laplace_driver<real_t,R1C2,Noco>(Tpsi, psi, list, f, num[dd], Stride, phase[dd], nFD[dd]);
            nops += nnzb*(2*nFD[dd] + 1)*R1C2*(Noco*64ul)*(Noco*64ul)*2ul; // total number of floating point operations performed
        } // dd
        char const fF = (8 == sizeof(real_t)) ? 'F' : 'f';
        if (echo > 1) std::printf("# green_kinetic::%s nFD= %d %d %d, %.3f M%clop\n", __func__, nFD[0], nFD[1], nFD[2], nops*1e-6, fF);
        return nops;
    } // multiply (kinetic energy operator)


    template <typename real_t, int R1C2=2, int Noco=1>
    size_t multiply(
          real_t         (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // result
        , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // 
        , green_sparse::sparse_t<int32_t> const kinetic_plan[3]
        , double   const hgrid[3] // grid spacing in X,Y,Z
        , int const FD_range=4 // finite-difference stencil range (in grid points)
        , size_t const nnzb=1 // total number of non-zero blocks (to get the operations count correct)
        , int const echo=0
    ) {
        int const nFD[] = {FD_range, FD_range, FD_range};
        auto phase = get_memory<double[2][2]>(3, echo, "phase"); // --> TODO move into argument list
        set_phase(phase, nullptr, echo); // neutral (Gamma-point) phase factors
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
        if (echo > 3) std::printf("# green_kinetic::%s nFD= %d, number= %d %d %d\n", __func__, FD_range, num[0], num[1], num[2]);

        auto const nops = multiply<real_t,R1C2,Noco>(Tpsi, psi, num, lists[0], lists[1], lists[2], hgrid, nFD, phase, nnzb, echo);
        
        for (int dd = 0; dd < 3; ++dd) {
            free_memory(lists[dd]);
        } // dd
        free_memory(phase);
        return nops;
   } // multiply (kinetic energy operator)
                  

                  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <typename real_t, int R1C2=2, int Noco=1>
  inline status_t test_finite_difference(int const echo=0) {
      size_t const nnzb = 15;
      auto Tpsi = get_memory<real_t[R1C2][Noco*64][Noco*64]>(nnzb, echo, "Tpsi");
      auto  psi = get_memory<real_t[R1C2][Noco*64][Noco*64]>(nnzb, echo,  "psi");
      for (size_t i = 0; i < nnzb*R1C2*pow2(Noco*64); ++i) psi[0][0][0][i] = 0; // clear
      auto indx = get_memory<int32_t>(4 + nnzb + 4);
      set(indx, 4 + nnzb + 4, 0); for (size_t i = 0; i < nnzb; ++i) indx[4 + i] = i + 1;
      uint32_t const num[] = {1, 1, 1}; // number of lists
      int const nFD4[] = {4, 4, 4}, nFD8[] = {8, 8, 8};
      double const hgrid[] = {1, 1, 1};
      // merely check if memory accesses are allowed
      multiply<real_t,R1C2,Noco>(Tpsi, psi, num, &indx, &indx, &indx, hgrid, nFD4, nullptr, nnzb, echo);
      multiply<real_t,R1C2,Noco>(Tpsi, psi, num, &indx, &indx, &indx, hgrid, nFD8, nullptr, nnzb, echo);


      // test also the periodic case (nFD4 only)
      for (int i = 0; i < 4; ++i) { indx[i] = 4 - (nnzb + i + 1); indx[4 + nnzb + i] = -(i + 1); }
      if (echo > 7) { std::printf("# periodic indices: "); printf_vector(" %d", indx, 4 + nnzb + 4); }
      double const phase_angle[] = {(2 == R1C2) ? -1./3. : -.5, (2 == R1C2) ? 0.125 : 0, 0.5}; // in units of 2*pi
      // if we use real numbers only, phase_angle must be a multiple of 0.5

      double const wave_vector[] = {-2*constants::pi*phase_angle[0]/(nnzb*4), 0, 0};
      // set up a periodic wave in the lowest rhs
      int const Re = 0, Im = R1C2 - 1;
      for (size_t inzb = 0; inzb < nnzb; ++inzb) {
          int const iz = 0, iy = 0;
          for (int ix = 0; ix < 4; ++ix) {
              auto const arg = wave_vector[0]*(inzb*4 + ix)*hgrid[0];
              psi[inzb][Re][(iz*4 + iy)*4 + ix][0] = std::cos(arg);
              psi[inzb][Im][(iz*4 + iy)*4 + ix][0] = std::sin(arg);
          } // ix
      } // inzb
      auto phase = get_memory<double[2][2]>(3, echo, "phase"); // --> TODO move into argument list
      set_phase(phase, phase_angle, 2*echo);
      uint32_t const num100[] = {1, 0, 0}; // only derive in x-direction
      multiply<real_t,R1C2,Noco>(Tpsi, psi, num100, &indx, &indx, &indx, hgrid, nFD4, phase, nnzb, echo);
      free_memory(phase);

      double deviation{0};
      { // scope: check
          auto const t = 0.5*pow2(wave_vector[0]);
          double dev2{0}, deva{0}, norm2{0};
          if (echo > 9) std::printf("\n## x  Tpsi_Re Tpsi_Im  psi_Re, psi_Im:"
                      " (wave vector k= %g, k^2/2= %g)\n", wave_vector[0], t);
          for (size_t inzb = 0; inzb < nnzb; ++inzb) {
              int const iz = 0, iy = 0;
              for (int ix = 0; ix < 4; ++ix) {
                  int const i64 = (iz*4 + iy)*4 + ix;
                  if (echo > 9) {
                      double const x = (inzb*4 + ix)*hgrid[0];
                      std::printf("%g  %g %g  %g %g\n", x,
                           Tpsi[inzb][Re][i64][0],  Tpsi[inzb][Im][i64][0]*(2 == R1C2),
                          t*psi[inzb][Re][i64][0], t*psi[inzb][Im][i64][0]*(2 == R1C2));
                  } // echo
                  for(int reim = 0; reim < R1C2; ++reim) {
                      auto const dev = Tpsi[inzb][reim][i64][0] - t*psi[inzb][reim][i64][0]; // deviation
                      deva  += std::abs(dev); // measure the absolute deviation
                      dev2  += pow2(dev); // measure the square deviation
                      norm2 += pow2(psi[inzb][reim][i64][0]);
                  } // reim
              } // ix
          } // inzb
          deviation = std::sqrt(dev2);
          if (echo > 5) std::printf("# deviation of a periodic Bloch wave with wave vector k= %g is %.1e (abs)"
                            " and %.1e (rms) at norm^2 %g\n\n", wave_vector[0], deva, deviation, norm2);
      } // scope

      free_memory(Tpsi);
      free_memory(psi);
      free_memory(indx);
      return (deviation > ((8 == sizeof(real_t)) ? 2e-15 : 5e-6));
  } // test_finite_difference

  inline status_t test_set_phase(int const echo=0) {
      double phase[2][2], maxdev{0};
//    for (int iangle = -181; iangle <= 181; ++iangle) {
      for (int iangle = -180; iangle <= 180; iangle += 5) {
          auto const dev = set_phase(phase, iangle/360., 't', echo);
          maxdev = std::max(maxdev, dev);
      } // iangle
      return (maxdev > 3e-16);
  } // test_set_phase

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_set_phase(echo);
      stat += test_finite_difference<float ,1,1>(echo);
      stat += test_finite_difference<float ,2,1>(echo);
      stat += test_finite_difference<float ,2,2>(echo);
      stat += test_finite_difference<double,1,1>(echo);
      stat += test_finite_difference<double,2,1>(echo);
      stat += test_finite_difference<double,2,2>(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_kinetic
