#pragma once

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <complex> // std::complex

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_memory.hxx" // dim3, get_memory, free_memory
#include "inline_math.hxx" // pow2, set
#include "mpi_parallel.hxx" // MPI_Comm, MPI_UINT16, MPI_COMM_WORLD
#include "global_coordinates.hxx" // ::get
#include "recorded_warnings.hxx" // error
#include "print_tools.hxx" // printf_vector

namespace green_potential {

    template <typename real_t, int R1C2=2, int Noco=1>
    void __global__ Potential( // GPU kernel, must be launched with <<< {64, any, 1}, {Noco*64, Noco, R1C2} >>>
#ifdef HAS_NO_CUDA
          dim3 const & gridDim, dim3 const & blockDim,
#endif // HAS_NO_CUDA
          real_t        (*const __restrict__ Vpsi)[R1C2][Noco*64][Noco*64] // result
        , real_t  const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input Green function
        , double  const (*const *const __restrict__ Vloc)[64] // local potential, Vloc[Noco*Noco][iloc][4*4*4]
        , int32_t const (*const __restrict__ iloc_of_inzb) // translation from inzb to iloc, [inzb]
        , int16_t const (*const __restrict__ shift)[3+1] // 3D block shift vector (target minus source), 4th component unused, [inzb][0:2]
        , double  const (*const __restrict__ hxyz) // grid spacing in X,Y,Z direction
        , int     const nnzb // number of all blocks to be treated
        , float   const Vconf // prefactor for the confinement potential
        , float   const rcut2 // cutoff radius^2 for the confinement potential, negative for no confinement
        , real_t  const E_real // real      part of the energy parameter, this could be subtracted from Vloc beforehand
        , real_t  const E_imag // imaginary part of the energy parameter
    ) {
        assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));

        assert(64      ==  gridDim.x);
        assert(1       ==  gridDim.z);
        assert(Noco*64 == blockDim.x);
        assert(Noco    == blockDim.y);
        assert(R1C2    == blockDim.z);

        bool const imaginary = ((2 == R1C2) && (0 != E_imag));

#ifndef HAS_NO_CUDA
        int const inz0 = blockIdx.y;  // start of the grid-stride loop on y-blocks
        int const i64  = blockIdx.x;  // target grid point index == inner column dimension, in [0, 64)
        int const reim = threadIdx.z; // real or imaginary part of the Green function
        int const spin = threadIdx.y; // non-collinear spin index
        int const j64  = threadIdx.x; // source grid point index == right hand side vectorization, in [0, Noco*64)
#else  // HAS_NO_CUDA
        assert(1 == gridDim.y && "CPU kernel Potential needs increment 1 for grid stride loop");
        int constexpr inz0 = 0;
        for (int i64 = 0; i64 < 64; ++i64)
        for (int reim = 0; reim < R1C2; ++reim)
        for (int spin = 0; spin < Noco; ++spin)
        for (int j64 = 0; j64 < Noco*64; ++j64)
#endif // HAS_NO_CUDA
        { // threads loops and block loop

#define CONFINEMENT_POTENTIAL
#ifdef  CONFINEMENT_POTENTIAL
        // generate the position difference of grid points (target minus source)
        int const x = ( i64       & 0x3) - ( j64       & 0x3);
        int const y = ((i64 >> 2) & 0x3) - ((j64 >> 2) & 0x3);
        int const z = ((i64 >> 4) & 0x3) - ((j64 >> 4) & 0x3);
        // due to the masking, we ignore the Noco-spin index inside j64
//      std::printf("# block %d thread %d has in-block shifts %d %d %d \n", i64, j64, x, y, z); // debug
#endif // CONFINEMENT_POTENTIAL

        assert(R1C2 >= Noco); // static_assert?

        auto const V_imag = E_imag * real_t(1 - 2*reim);

        for (int inzb = inz0; inzb < nnzb; inzb += gridDim.y) { // grid-stride loop on y-blocks

            real_t Vconfine = 0; // non-const
#ifdef  CONFINEMENT_POTENTIAL
            if (rcut2 >= 0.f) {
                int constexpr n4 = 4;
                auto const s = shift[inzb]; // shift vectors between target minus source cube
                double const d2 = pow2((s[0]*n4 + x)*hxyz[0])
                                + pow2((s[1]*n4 + y)*hxyz[1])
                                + pow2((s[2]*n4 + z)*hxyz[2]);
                auto const d2out = real_t(d2 - rcut2);
                Vconfine = (d2out > 0) ? Vconf*pow2(d2out) : 0; // quartic confinement potential, V ~ d^4
//              std::printf("%.1f %.3f  %d %d %d  %d\n", d2, Vconfine, s[0], s[1], s[2], s[3]); // s[3] == source block index
            } // rcut^2 >= 0
#endif // CONFINEMENT_POTENTIAL

            auto const iloc = iloc_of_inzb[inzb]; // target index for the local potential, can be -1 for non-existing
            real_t const Vloc_diag = (iloc < 0) ? 0 : Vloc[spin][iloc][i64];

            // gather all real-valued and spin-diagonal contributions
            real_t const Vtot = Vloc_diag + Vconfine - E_real; // diagonal part of the potential

            auto vpsi = Vtot * psi[inzb][reim][spin*64 + i64][j64]; // non-const, potential is diagonal in real-space

            if (imaginary) {
                // V-E has an imaginary part V_Im = -E_imag
                // then explicitly:
                //    Vpsi_Re = V_Re * psi_Re - V_Im * psi_Im = V_Re * psi_Re + E_imag * psi_Im (reim=0)
                //    Vpsi_Im = V_Re * psi_Im + V_im * psi_Re = V_Re * psi_Im - E_imag * psi_Re (reim=0)
                vpsi += V_imag * psi[inzb][1 - reim][spin*64 + i64][j64];
            } // imaginary

            if (2 == Noco && iloc >= 0) { // the other spin component is (1 - spin)
                /*                                                                   */
                /*  this code would be correct if a noco potential had real values,  */
                /*  however, it has 4 components (V_0, V_x, V_y, V_z)                */
                /*                                                                   */
                /*           /  V_0 + V_z   V_x-i*V_y \    / V_dndn   V_x-i*V_y \    */
                /*  V_noco = |                        | =  |                    |    */
                /*           \  V_x+i*V_y   V_0 - V_z /    \ V_x+i*V_y   V_upup /    */
                /*                                                                   */
                /*  however, to avoid memory accesses we store                       */
                /*  these four combinations in Vloc[0..3][icube][i64]                */
                /*  Vdndn = V_1 + V_z, Vupup = V_1 - V_z, V_x and V_y                */
                /*                                                                   */
                /*  Explicitly:                                                      */
                /*    Vpsi_dn_Re = V_dndn psi_dn_Re + V_x psi_up_Re + V_y psi_up_Im  */
                /*    Vpsi_dn_Im = V_dndn psi_dn_Im + V_x psi_up_Im - V_y psi_up_Re  */
                /*                                                                   */
                /*    Vpsi_up_Re = V_upup psi_up_Re + V_x psi_dn_Re - V_y psi_dn_Im  */
                /*    Vpsi_up_Im = V_upup psi_up_Im + V_x psi_dn_Im + V_y psi_dn_Re  */
                /*                                                                   */
                /*                                                                   */
                real_t const cs = (1 - 2*(reim ^ spin)); // complex sign is -1 if (reim != spin)

                vpsi += Vloc[2][iloc][i64] * psi[inzb][    reim][(1 - spin)*64 + i64][j64];    // V_x
                vpsi += Vloc[3][iloc][i64] * psi[inzb][1 - reim][(1 - spin)*64 + i64][j64]*cs; // V_y
            } // non-collinear

            Vpsi[inzb][reim][spin*64 + i64][j64] = vpsi; // store

        } // inzb

        } // thread loops and block loop

    } // Potential


    template <typename real_t, int R1C2=2, int Noco=1>
    size_t multiply(
          real_t         (*const __restrict__ Vpsi)[R1C2][Noco*64][Noco*64] // result
        , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input
        , double   const (*const *const __restrict__ Vloc)[64] // local potential, Vloc[Noco*Noco][iloc][4*4*4]
        , int32_t  const (*const __restrict__ vloc_index) // iloc_of_inzb[nnzb]
        , int16_t  const (*const __restrict__ shift)[3+1] // 3D block shift vector (target minus source), 4th component unused
        , double   const (*const __restrict__ hxyz) // grid spacing in X,Y,Z direction
        , int      const nnzb // number of all blocks to be treated
        , std::complex<double> const E_param=0 // energy parameter, the real part could be subtracted from Vloc beforehand
        , float    const Vconf=0  // prefactor for the confinement potential
        , float    const rcut2=-1 // cutoff radius^2 for the confinement potential, -1: no confinement
        , int const echo=0
    ) {

        if (echo > 11) {
            std::printf("# %s<%s,R1C2=%d,Noco=%d> Vpsi=%p, psi=%p, Vloc=%p, vloc_index=%p, shift=%p, hxyz=%p, nnzb=%d, Vconf=%g, rcut2=%.f, E=(%g, %g)\n",
                           __func__, (4 == sizeof(real_t))?"float":"double", R1C2, Noco, (void*)Vpsi, (void*)psi,
                           (void*)Vloc, (void*)vloc_index, (void*)shift, (void*)hxyz, nnzb, Vconf, rcut2, E_param.real(), E_param.imag());
        } // echo

        Potential<real_t,R1C2,Noco>
#ifndef HAS_NO_CUDA
            <<< dim3(64, 7, 1), dim3(Noco*64, Noco, R1C2) >>> ( // 7=any, maybe find a function for a good choice
#else  // HAS_NO_CUDA
              ( dim3(64, 1, 1), dim3(Noco*64, Noco, R1C2),
#endif // HAS_NO_CUDA
            Vpsi, psi, Vloc, vloc_index, shift, hxyz, nnzb, Vconf, rcut2, E_param.real(), E_param.imag());

        return 0ul; // total number of floating point operations performed
    } // multiply potential

  inline char const * spin_name(int const Noco, int const spin) {
      if (1 == Noco && 0 == spin) return "";
      if (2 == Noco) {
          switch (spin) {
            case 0: return " V_down";
            case 1: return " V_up";
            case 2: return " V_x";
            case 3: return " V_y";
          } // spin
      } // 2 == Noco
      return " ???";
  } // spin_name

  // For MPI parallel calculations the potential values need to be exchanged
  template <typename rank_int_t=uint16_t>
  status_t exchange(
        double    (*const Veff[4])[64]  // output effective potentials, data layout Veff[Noco^2][nrows][64]
      , std::vector<int64_t> const & requests  // indices requested by this MPI process, [nrows]
      , double const (*const Vinp)[64]  //  input effective potentials, data layout Vinp[ncols*Noco^2 ][64]
      , std::vector<int64_t> const & offerings // indices offered by this MPI process, [ncols]
      , rank_int_t const owner_rank[] // where to find it, [nb[Z]*nb[Y]*nb[X]]
      , uint32_t const nb[3] // global bounding box, maybe group this with owner_rank into a translator object global_index --> owner_rank
      , int const Noco=1 // 1:no spin, 2: (non-collinear) spin
      , bool const debug=true
      , MPI_Comm const comm=MPI_COMM_WORLD // MPI communicator handle
      , int const echo=0 // log-level
  ) {
      if (echo > 0) std::printf("# MPI exchange of potential, MPI one-sided communication, Noco=%d\n", Noco);
      assert(1 == Noco || 2 == Noco);

      assert(Veff && "may not be called with a nullptr for output");
      for (int spin = 0; spin < Noco*Noco; ++spin) {
          assert(Veff[spin] && "may not be called with a nullptr for spin output");
      } // spin
      assert(Vinp && "may not be called with a nullptr for input");

      status_t status(0);

      auto const me = mpi_parallel::rank(comm);
      auto const ncols = offerings.size(); // number of offerings
      auto const nrows =  requests.size(); // number of requests
      size_t stats[] = {0, 0};

#ifndef   HAS_NO_MPI

      int constexpr X=0, Y=1, Z=2;
      auto const nall = nb[Z]*size_t(nb[Y])*size_t(nb[X]);

      std::vector<uint16_t> local_index(nall, 0);
      std::vector<uint16_t> local_check(nall*int(debug), 0);

      assert(ncols <= (1ul << 16) && "each process can hold max 2^16 Green function columns");
      for (size_t col = 0; col < ncols; ++col) {
          auto const global_id = offerings[col];

          // translate global_id into an box index iall
          uint32_t xyz[3]; global_coordinates::get(xyz, global_id);
          for (int d = 0; d < 3; ++d) {
              assert(xyz[d] < nb[d]);
          } // d
          auto const iall = (xyz[Z]*size_t(nb[Y]) + xyz[Y])*nb[X] + xyz[X];
          assert(iall < nall);

          assert(me == owner_rank[iall] && "offerings must be owned");

          local_index[iall] = col;
          if (debug) ++local_check[iall];
      } // col

      if (debug) {
          if (echo > 7) { std::printf("# local_check before "); printf_vector(" %i", local_check); }
          for (size_t iall = 0; iall < nall; ++iall) {
              assert(local_check[iall] <= 1 && "duplicates found");
          } // iall
          status += MPI_Allreduce(MPI_IN_PLACE, local_check.data(), nall, MPI_UINT16, MPI_SUM, comm);

          if (echo > 7) { std::printf("# local_check after  "); printf_vector(" %i", local_check); }
          for (size_t iall = 0; iall < nall; ++iall) {
              assert(1 == local_check[iall] && "not all covered");
          } // iall
          local_check.resize(0); // DEBUG
      } // debug

      // get a global list of which local index is where
      status += MPI_Allreduce(MPI_IN_PLACE, local_index.data(), nall, MPI_UINT16, MPI_MAX, comm);

      // set up a memory window to read from
      MPI_Win window;
      size_t const disp_unit = 64*sizeof(double); // in Bytes
      size_t const win_size  = ncols*Noco*Noco*disp_unit; // in Bytes
      status += MPI_Win_create((void*)Vinp, win_size, disp_unit, MPI_INFO_NULL, comm, &window);

      // synchronize processes
      int const assertions = 0; // use bitwise or, e.g. MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE | MPI_MODE_NOSUCCEED;
      status += MPI_Win_fence(assertions, window);

#endif // HAS_NO_MPI

      // look up who owns the requested cube
      for (size_t row = 0; row < nrows; ++row) {
          auto const global_id = requests[row];

#ifndef   HAS_NO_MPI

          // translate global_id into an owner rank
          uint32_t xyz[3]; global_coordinates::get(xyz, global_id);
          for (int d = 0; d < 3; ++d) {
              assert(xyz[d] < nb[d]);
          } // d
          auto const iall = (xyz[Z]*size_t(nb[Y]) + xyz[Y])*nb[X] + xyz[X];
          assert(iall < nall);
          auto const owner = owner_rank[iall];
          auto const iloc = local_index[iall];

#else  // HAS_NO_MPI

          // version without MPI
          auto const owner = me;
          // search global_id in offerings
          int64_t iloc{-1};
          for (size_t col = 0; col < ncols && iloc < 0; ++col) {
              if (global_id == offerings[col]) iloc = col;
          } // col
          assert(iloc > -1 && "index not found in local offerings");

#endif // HAS_NO_MPI

          for (int spin = 0; spin < Noco*Noco; ++spin) {
              auto const target_disp = iloc*Noco*Noco + spin;
              if (me == owner) {
                  if (echo > 5 + 5*spin) std::printf("# exchange: rank #%i local copy of potential at cube id o%llo%s\n", me, global_id, spin_name(Noco, spin));
                  set(Veff[spin][row], 64, Vinp[target_disp]);
                  ++stats[0];
              } else { // me == owner
                  ++stats[1];
#ifndef   HAS_NO_MPI
                  if (echo > 7 + 5*spin) std::printf("# exchange: rank #%i get potential at cube id o%llo from rank #%i element %i%s\n", me, global_id, owner, iloc, spin_name(Noco, spin));
                  status += MPI_Get(Veff[spin][row], 64, MPI_DOUBLE, owner, target_disp, 64, MPI_DOUBLE, window);
                  /*
                   *  Criticism for MPI_Get-solution:
                   *    each time called it pushes a request into the remote process and receives an MPI_Put with the data from that process.
                   *    After the construction of the Green function, the structure of requested potential elements does not change any more.
                   *    Latencies would be lower if we separated it into two phases:
                   *    a) communicate where to push, e.g. with MPI_Alltoall and MPI_Alltoallv before the loop over energy parameters.
                   *    b) MPI_Put the data
                   */
#else  // HAS_NO_MPI
                  error("Without MPI all potential elements must reside in the same process, me=%i, owner=%i", me, owner);
#endif // HAS_NO_MPI
              } // me == owner
          } // spin
      } // row

#ifndef   HAS_NO_MPI
      // synchronize processes
      status += MPI_Win_fence(assertions, window);
      // MPI_Win_free is only needed if we use MPI_Win_alloc above
#endif // HAS_NO_MPI

      if (echo > 7) std::printf("# rank #%i copied %.3f k and pulled %.3f k elements\n", me, stats[0]*.001, stats[1]*.001);
      assert(stats[0] + stats[1] == Noco*Noco*nrows);

      return status;
  } // exchange


#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <typename real_t, int R1C2=2, int Noco=1>
  inline status_t test_multiply(int const echo=0) {
      auto psi   = get_memory<real_t[R1C2][Noco*64][Noco*64]>(1);
      auto Vloc  = get_memory<double[64]>(1*Noco*Noco);
      auto vloc_index = get_memory<int32_t>(1);         vloc_index[0] = 0;
      auto shift = get_memory<int16_t[3+1]>(1);         set(shift[0], 3+1, int16_t(0));
      auto hGrid = get_memory<double>(3+1);             set(hGrid, 3+1, 1.);
      int const nnzb = 1;

      multiply<real_t,R1C2,Noco>(psi, psi, &Vloc, vloc_index, shift, hGrid, nnzb);

      free_memory(hGrid);
      free_memory(shift);
      free_memory(vloc_index);
      free_memory(Vloc);
      free_memory(psi);
      return 0;
  } // test_multiply

  inline status_t test_multiply(int const echo=0) {
      status_t stat(0);
      stat += test_multiply<float ,1,1>(echo);
      stat += test_multiply<float ,2,1>(echo);
      stat += test_multiply<float ,2,2>(echo);
      stat += test_multiply<double,1,1>(echo);
      stat += test_multiply<double,2,1>(echo);
      stat += test_multiply<double,2,2>(echo);
      return stat;
  } // test_multiply

  template <int Noco=1>
  inline status_t test_exchange(int echo=0) {
      status_t stat(0);
      uint32_t const nb[] = {2, 2, 1};
      std::vector<int64_t> requests = {0, 3, 2, 1};
      auto const nrows = requests.size();
      double (*Veff[Noco*Noco])[64];
      for (int spin = 0; spin < Noco*Noco; ++spin) Veff[spin] = (double(*)[64])malloc(nrows*64*sizeof(double));
      auto const Vinp = new double[nrows][64];
      int const me = mpi_parallel::rank(),
                np = mpi_parallel::size(); assert(np > 0);
      std::vector<int64_t> offerings(0);
      std::vector<uint16_t> owner_rank(nrows, 0);
      for(int row{0}; row < nrows; ++row) {
          int const rank = row % np; // block-cyclic distribution
          owner_rank[row] = rank;
          if (me == rank) offerings.push_back(row);
      } // row
      stat += exchange(Veff, requests, Vinp, offerings, owner_rank.data(), nb, Noco, true, MPI_COMM_WORLD, echo);
      for (int spin = 0; spin < Noco*Noco; ++spin) free(Veff[spin]);
      delete[] Vinp;
      return stat;
  } // test_exchange

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
//       stat += test_multiply(echo);
      stat += test_exchange<1>(echo);
      stat += test_exchange<2>(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_potential
