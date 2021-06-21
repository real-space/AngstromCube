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

  // ToDo: move to green_utils.hxx or similar
  template <typename uint_t, typename int_t> inline
  size_t index3D(uint_t const n[3], int_t const i[3]) {
      // usual 3D indexing
      return size_t(i[2]*n[1] + i[1])*n[0] + i[0];
  } // index3D


namespace green_kinetic {

  class finite_difference_plan_t {
  private:
      uint32_t *prefix; // in managed memory
      int32_t *fd_list; // in managed memory
      uint32_t n_lists;

  public:
      finite_difference_plan_t() : prefix(nullptr), fd_list(nullptr), n_lists(0) {}

      finite_difference_plan_t(
            int const dd // direction of derivative
          , uint16_t const num_target_coords[3]
          , uint32_t const RowStart[]
          , uint16_t const ColIndex[]
          , view3D<int32_t> const & iRow_of_coords // (Z,Y,X)
          , std::vector<bool> const sparsity_pattern[]
          , unsigned const nRHSs=1
          , int const echo=0
      ) {
          int constexpr X=0, Y=1, Z=2;
          // prepare the finite-difference sequence lists
          char const direction = 'x' + dd;
          assert(X == dd || Y == dd || Z == dd); 
          int num[3];
          set(num, 3, num_target_coords);
          int const num_dd = num[dd];
          num[dd] = 1; // replace number of target blocks in derivative direction
          if (echo > 0) std::printf("# FD lists for the %c-direction %d %d %d\n", direction, num[X], num[Y], num[Z]);
          simple_stats::Stats<> length_stats;
          std::vector<std::vector<int32_t>> list;
          size_t const max_lists = nRHSs*size_t(num[Z])*size_t(num[Y])*size_t(num[X]);
          list.resize(max_lists);
          size_t ilist{0};
          for (int iRHS = 0; iRHS < nRHSs; ++iRHS) {
//                if (echo > 0) std::printf("# FD list for RHS #%i\n", iRHS);
              auto const & sparsity_RHS = sparsity_pattern[iRHS];
              for (int iz = 0; iz < num[Z]; ++iz) { //  
              for (int iy = 0; iy < num[Y]; ++iy) { //   only 2 of these 3 loops have a range > 1
              for (int ix = 0; ix < num[X]; ++ix) { // 
                  int idx[3] = {ix, iy, iz};
                  for (int id = 0; id < num_dd; ++id) { // loop over direction to derive
                      idx[dd] = id; // replace index in the derivate direction
//                           if (echo > 0) std::printf("# FD list for RHS #%i test coordinates %i %i %i\n",
//                                                   iRHS, idx[X], idx[Y], idx[Z]);
                      auto const idx3 = index3D(num_target_coords, idx);
                      if (sparsity_RHS[idx3]) {
                          auto const iRow = iRow_of_coords(idx[Z], idx[Y], idx[X]);
                          assert(iRow >= 0);

                          int32_t inz_found{-1};
                          for (auto inz = RowStart[iRow]; inz < RowStart[iRow + 1]; ++inz) {
                              if (ColIndex[inz] == iRHS) {
                                  inz_found = inz; // store where it was found
                                  inz = RowStart[iRow + 1]; // stop search loop
                              } // found
                          } // search
                          assert(inz_found >= 0); // fails at inconsistency between sparsity_pattern and the BSR tables

                          assert(ilist < max_lists);
                          list[ilist].push_back(inz_found);
                      } // sparsity pattern
                  } // id
                  int const list_length = list[ilist].size();
                  if (list_length > 0) {
                      length_stats.add(list_length);
//                           if (echo > 0) std::printf("# FD list of length %d for the %c-direction %i %i %i\n",
//                                                   list_length, direction, idx[X], idx[Y], idx[Z]);
                      // add 4 end-of-sequence markers
                      for (int i = 0; i < 4; ++i) {
                          list[ilist].push_back(-1);
                      } // i
                      
                      ++ilist; // create a new list index
                  } // list_length > 0
              }}} // ixyz
          } // iRHS
          n_lists = ilist; assert(n_lists == ilist && "too many lists, max. 2^32-1");
          if (echo > 0) std::printf("# %d FD lists for the %c-direction (%.2f %%), length %.3f +/- %.3f, min %g max %g\n",
                                n_lists, direction, n_lists/(max_lists*.01),
                                length_stats.avg(), length_stats.var(), length_stats.min(), length_stats.max());

          // store index starts in managed memory
          prefix = get_memory<uint32_t>(n_lists + 1); // create in GPU memory
          prefix[0] = 0;
          for (int ilist = 0; ilist < n_lists; ++ilist) {
              int const n = list[ilist].size();
              prefix[ilist + 1] = prefix[ilist] + n;
          } // ilist
          size_t const ntotal = prefix[n_lists];
          if (echo > 0) std::printf("# FD lists for the %c-direction require %ld uint32_t, i.e. %.3f kByte\n",
                                  direction, ntotal, ntotal*sizeof(uint32_t)*1e-3);
          // store indices in managed memory
          fd_list = get_memory<int32_t>(ntotal); // create in GPU memory
          { // scope: copy indices into managed memory
              size_t ntotal_check{0};
              for (int ilist = 0; ilist < n_lists; ++ilist) {
                  int const n = list[ilist].size();
                  assert(ntotal_check == prefix[ilist]); // sanity
                  ntotal_check += n;
                  set(&fd_list[prefix[ilist]], n, list[ilist].data()); // copy into GPU memory
              } // ilist
              assert(ntotal == ntotal_check); // sanity
          } // scope
 
      } // constructor

      finite_difference_plan_t& operator= (finite_difference_plan_t && rhs) {
          // std::printf("# finite_difference_plan_t& operator= (finite_difference_plan_t && rhs);\n");
          std::swap(fd_list, rhs.fd_list);
          std::swap(prefix , rhs.prefix);
          std::swap(n_lists, rhs.n_lists);
          return *this;
      } // move assignment

      finite_difference_plan_t(finite_difference_plan_t && rhs) = delete;
//       {   std::printf("# finite_difference_plan_t(finite_difference_plan_t && rhs);\n");
//           *this = std::move(rhs);
//       } // move constructor

      finite_difference_plan_t(finite_difference_plan_t const & rhs) = delete; // copy constructor

      finite_difference_plan_t& operator= (finite_difference_plan_t const & rhs) = delete; // move assignment

      ~finite_difference_plan_t() {
          std::printf("# destruct %s, pointers= %p and %p\n", __func__, (void*)fd_list, (void*)prefix); std::fflush(stdout);
          free_memory(fd_list);
          free_memory(prefix);
      } // destructor

      uint32_t size() const { return n_lists; } // number of lists
      int32_t const * list(uint32_t const i) const { return fd_list + prefix[i]; }

  }; // class finite_difference_plan_t

  
#define NO_CUDA
#define __global__
#define __restrict__

	template <typename real_t, int Stride, int R1C2=2, int Noco=1> // Stride is determined by the lattice dimension along which we derive
	void __global__ Laplace8th( // GPU kernel, must be launched with <<< <<< {16, Nrows, 1}, {Noco*64, Noco, R1C2} >>>
          real_t        (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // intent(inout)
		, real_t  const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // intent(in)
        , int32_t const (*const *const __restrict__ index_list) // index list that brings the blocks in order,
                                          // list must contain at least one element and is finalized with -1
		, double const prefactor=1
		, int const Nrows=1 // only needed when running whithout CUDA
    ) {
        // prepare finite-difference coefficients
		double const norm = prefactor/5040.; // {-14350, 8064, -1008, 128, -9}/5040. --> 8th order
		real_t const  c0 = -14350*norm,
                      c1 =   8064*norm,
                      c2 =  -1008*norm,
                      c3 =    128*norm,
                      c4 =     -9*norm; //-> NStep = 4
		// the NStep=4 8th order stencil has a ratio of c0/c4 = 2^10.64 so half precision (11 significant bits) is the limit

		// the following list gives the indices of blocks that belong to the same right hand side 
		// and are neighbors in the direction of derivation

#ifdef NO_CUDA
        for (int block_y = 0; block_y < Nrows; ++block_y) {
        for (int block_x = 0; block_x < 16; ++block_x) {
        for (int thread_z = 0; thread_z < R1C2; ++thread_z) {
        for (int thread_y = 0; thread_y < Noco; ++thread_y) {
        for (int thread_x = 0; thread_x < Noco*64; ++thread_x) {
#else
        { int const block_y = blockIdx.y;
        { int const block_x = blockIdx.x;
        { int const thread_z = threadIdx.z;
        { int const thread_y = threadIdx.y;
        { int const thread_x = threadIdx.x;
#endif

		auto const *const list = index_list[block_y]; // abbreviate pointer

		int const i16 = block_x;
		int const i64 = (16==Stride)? i16 : ( (4==Stride)? (16*(i16 >> 2) + (i16 & 0x3)) : (4*i16) );
#define INDICES(i4) [thread_z][thread_y*64 + i64 + Stride*i4][thread_x]

		real_t w0{0}, w1{0}, w2{0}, w3{0}, // initialize one non-existing block
		       w4, w5, w6, w7, wn; // 4 + 4 + 1 registers, wn is the register that always receives the most recently loaded value
        int ilist{0}; // counter for index_list
        int ii, ip; // elements of index_list

		// initially load one block in advance
            ii = list[ilist++]; // at least on element must be nonzero (can be issued earlier)
        w4 = psi[ii]INDICES(0); // inital load
        w5 = psi[ii]INDICES(1); // inital load
        w6 = psi[ii]INDICES(2); // inital load
        w7 = psi[ii]INDICES(3); // inital load

        // main loop
        // so far ilist == 1
		while (ii > -1) {
            ip = ii; // set previous index 
            ii = list[ilist++]; // get next index
            bool const load = (ii > -1);
            // use a rotating register file, see figs/rotating_register_file.fig or *.pdf

            // FD9POINT = load?, compute, store, update rotating register file
#define FD9POINT(i4,  M4, M3, M2, M1, W0, P1, P2, P3, P4) \
            P4 = load ? psi[ii]INDICES(i4) : 0; \
            Tpsi[ip]INDICES(i4) += c0*W0 + c1*M1 + c1*P1 + c2*M2 + c2*P2 + c3*M3 + c3*P3 + c4*M4 + c4*P4; \
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

#undef  INDICES

        }}}}} // close 5 loops

    } // Laplace8th
    
    
    
	template <typename real_t, int Stride, int R1C2=2, int Noco=1> // Stride is determined by the lattice dimension along which we derive: 1, 4 or 16
	void __global__ Laplace16th( // GPU kernel, must be launched with <<< {16, Nrows, 1}, {Noco*64, Noco, R1C2} >>>
          real_t        (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // intent(inout)
		, real_t  const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // intent(in)
        , int32_t const (*const *const __restrict__ index_list) // index list that brings the blocks in order,
                                            // list must contain at least one element and is finalized with -1
		, double const prefactor=1 // alternatively, we could precalculate and load the 9 coefficients c0...c8
		, int const Nrows=1 // only needed when running whithout CUDA
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

#ifdef NO_CUDA
        for (int block_y = 0; block_y < Nrows; ++block_y) {
        for (int block_x = 0; block_x < 16; ++block_x) {
        for (int thread_z = 0; thread_z < R1C2; ++thread_z) {
        for (int thread_y = 0; thread_y < Noco; ++thread_y) {
        for (int thread_x = 0; thread_x < Noco*64; ++thread_x) {
#else
        { int const block_y = blockIdx.y;
        { int const block_x = blockIdx.x;
        { int const thread_z = threadIdx.z;
        { int const thread_y = threadIdx.y;
        { int const thread_x = threadIdx.x;
#endif

		auto const *const list = index_list[block_y]; // abbreviate pointer

		int const i16 = block_x;
		int const i64 = (16==Stride)? i16 : ( (4==Stride)? (16*(i16 >> 2) + (i16 & 0x3)) : (4*i16) );
#define INDICES(i4) [thread_z][thread_y*64 + i64 + Stride*i4][thread_x]

		real_t w0{0}, w1{0}, w2{0}, w3{0}, w4{0}, w5{0}, w6{0}, w7{0}, // initialize two non-existing blocks
               w8, w9, wa, wb, wc, wd, we, wf, wn; // 8 + 8 + 1 registers
        int ilist{0}; // counter for index_list
        int ii, ip, pp; // elements of index_list

		// initially load two blocks in advance
            ii = list[ilist++]; // load index for 1st non-zero block
        w8 = psi[ii]INDICES(0); // inital load
        w9 = psi[ii]INDICES(1); // inital load
        wa = psi[ii]INDICES(2); // inital load
        wb = psi[ii]INDICES(3); // inital load

            ip = ii; // set previous index
            ii = list[ilist++]; // load index for the second block
        if (ii > -1) {
            wc = psi[ii]INDICES(0); // inital load
            wd = psi[ii]INDICES(1); // inital load
            we = psi[ii]INDICES(2); // inital load
            wf = psi[ii]INDICES(3); // inital load
        } else {
           wc = 0; wd = 0; we = 0; wf = 0; // second block is already non-existing
        } // ii valid

        // main loop
            pp = ip; // set previous previous index 
            ip = ii; // set previous index 
            ii = list[ilist++]; // get next index
        // now ilist == 3
        while (ip > -1) {
            bool const load = (ii > -1); // load or not?
            // use a rotating register file, see above version of Laplace8th
            
            // FD17POINT = load?, compute, store, update rotating register file
#define FD17POINT(i4,  M8, M7, M6, M5, M4, M3, M2, M1, W0, P1, P2, P3, P4, P5, P6, P7, P8) \
            P8 = load ? psi[ii]INDICES(i4) : 0; \
            Tpsi[pp]INDICES(i4) += c0*W0 + c1*M1 + c1*P1 + c2*M2 + c2*P2 + c3*M3 + c3*P3 + c4*M4 + c4*P4 \
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
            pp = ip; // set previous previous index 
            ip = ii; // set previous index 
            ii = list[ilist++]; // get next index
        } // while loop

#undef  INDICES

        }}}}} // close 5 loops

    } // Laplace16th

    
#undef __global__

    template <typename real_t, int R1C2=2, int Noco=1>
    size_t multiply(
          real_t         (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // result
        , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // 
        , uint32_t const num[] // number of sticks in X,Y,Z direction
        , int32_t  const *const *const __restrict__ x_list // 
        , int32_t  const *const *const __restrict__ y_list // 
        , int32_t  const *const *const __restrict__ z_list // 
        , double   const hgrid=1 // grid spacing, ToDo make it an array of X,Y,Z
        , int nFD=4 // finite-difference stencil range (in grid points)
        , size_t const nnzb=1 // total number of non-zero blocks
    ) {

        double const f = -0.5/(hgrid*hgrid); // prefactor of the kinetic energy in Hartree atomic units
        // Laplace launch params <<< {n4*n4, Nrow, 1}, {Noco*64, R1C2, Noco} >>>
        if (8 == nFD) {
#ifdef NO_CUDA
            Laplace16th<real_t, 1,R1C2,Noco> (Tpsi, psi, x_list, f, num[0]); // x-direction
            Laplace16th<real_t, 4,R1C2,Noco> (Tpsi, psi, y_list, f, num[1]); // y-direction
            Laplace16th<real_t,16,R1C2,Noco> (Tpsi, psi, z_list, f, num[2]); // z-direction
        } else {
            Laplace8th <real_t, 1,R1C2,Noco> (Tpsi, psi, x_list, f, num[0]); // x-direction
            Laplace8th <real_t, 4,R1C2,Noco> (Tpsi, psi, y_list, f, num[1]); // y-direction
            Laplace8th <real_t,16,R1C2,Noco> (Tpsi, psi, z_list, f, num[2]); // z-direction
#else
            Laplace16th<real_t, 1,R1C2,Noco> <<< {16, num[0], 1}, {Noco*64, Noco, R1C2} >>> (Tpsi, psi, x_list, f);
            Laplace16th<real_t, 4,R1C2,Noco> <<< {16, num[1], 1}, {Noco*64, Noco, R1C2} >>> (Tpsi, psi, y_list, f);
            Laplace16th<real_t,16,R1C2,Noco> <<< {16, num[2], 1}, {Noco*64, Noco, R1C2} >>> (Tpsi, psi, z_list, f);
        } else {
            Laplace8th <real_t, 1,R1C2,Noco> <<< {16, num[0], 1}, {Noco*64, Noco, R1C2} >>> (Tpsi, psi, x_list, f);
            Laplace8th <real_t, 4,R1C2,Noco> <<< {16, num[1], 1}, {Noco*64, Noco, R1C2} >>> (Tpsi, psi, y_list, f);
            Laplace8th <real_t,16,R1C2,Noco> <<< {16, num[2], 1}, {Noco*64, Noco, R1C2} >>> (Tpsi, psi, z_list, f);
#endif
            nFD = 4;
        } // nFD

		return nnzb*3ul*(nFD+1+nFD)*R1C2*(Noco*64ul)*(Noco*64ul)*2ul; // total number of floating point operations performed
	} // multiply (kinetic energy operator)

#undef NO_CUDA

  
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_finite_difference(int const echo=0) {
      if (echo > 0) std::printf("# %s: no test included!\n", __func__);
      return STATUS_TEST_NOT_INCLUDED;
  } // test_finite_difference

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_finite_difference(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_kinetic
