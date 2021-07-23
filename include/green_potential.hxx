#pragma once

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_memory.hxx" // get_memory, free_memory
#include "sho_tools.hxx" // ::nSHO, ::n2HO, ::n1HO

namespace green_potential {
  
#define NO_CUDA
#define __global__
#define ____restrict____
#define __device__
#define __shared__
#define __unroll__
#define __host__

    
    template <typename real_t, int R1C2=2, int Noco=1>
    void __global__ Potential( // GPU kernel, must be launched with <<< {64, any, 1}, {Noco*64, Noco, R1C2} >>>
          real_t        (*const restrict Vpsi)[R1C2][Noco*64][Noco*64] // result
        , real_t  const (*const restrict  psi)[R1C2][Noco*64][Noco*64] //
        , real_t  const (*const restrict Vloc)[Noco*Noco][64] // local potential
        , size_t  const nnzb // number of all blocks to be treated
        , int     const (*const restrict CubeIndex) // translation from inzb to itrg, with uint16_t, we can cover a radius of 25 blocks
        , int16_t const (*const restrict shift)[4] // 3D block shift vector (target minus source), 4th component unused
        , float   const rcut2 // cutoff radius^2 for the confinement potential
        , real_t  const E_real // real      part of the energy parameter, this could be subtracted from Vloc beforehand
        , real_t  const E_imag // imaginary part of the energy parameter
    ) {

        int const reim = threadIdx.z; // real or imaginary part of the Green function
        int const spin = threadIdx.y; // non-collinear spin index
        int const j64 = threadIdx.x; // source == right hand side vectorization, in [0, Noco*64)
        int const i64 = blockIdx.x;  // target == inner column dimension, in [0, 64)
        assert(1       == gridDim.z);
        assert(64      == gridDim.x);
        assert(Noco*64 == blockDim.x);
        assert(Noco    == blockDim.y);
        assert(R1C2    == blockDim.z);

#define CONFINEMENT_POTENTIAL
#ifdef  CONFINEMENT_POTENTIAL
        // generate the position difference of grid points (target minus source)
        int const x = ( i64       & 0x3) - ( j64 >>    & 0x3);
        int const y = ((i64 >> 2) & 0x3) - ((j64 >> 2) & 0x3);
        int const z = ((i64 >> 4) & 0x3) - ((j64 >> 4) & 0x3);
        // due to the masking, we ignore the Noco-spin index inside j64
//      printf("# block %d thread %d has in-block shifts %d %d %d \n", i64, j64, x, y, z); // debug
#endif

        assert(R1C2 >= Noco); // static_assert?

        bool const imaginary = ((2 == R1C2) && (0 != E_imag));
        auto const V_imag = E_imag * real_t(1 - 2*reim);

        for (int inzb = blockIdx.y; inzb < nnzb; inzb += gridDim.y) { // grid-stride loop on y.blocks

#ifdef  CONFINEMENT_POTENTIAL
            auto const s = shift[inzb]; // shift vectors between target minus source cube
            double const d2 = pow2(x + n4*s[0]) + pow2(y + n4*s[1]) + pow2(z + n4*s[2]); // metric missing here
            auto const d2out = real_t(d2 - rcut2);
            real_t const Vconfine = (d2out > 0)? pow2(d2out) : 0; // confinement potential
//          printf("%.1f %.3f  %d %d %d  %d \n", d2, Vconfine, s[0], s[1], s[2], s[3]); // s[3] == source block index
#else
            real_t const Vconfine = 0; // no confinement potential
#endif

            int const icube = CubeIndex[inzb]; // target index for the local potential
            real_t const Vloc_diag = Vloc[icube][spin][i64]; // ToDo: activate and find the bug where the out-of-bounds happens

            // gather all real-valued and spin-diagonal contributions
            real_t const Vtot = Vloc_diag + Vconfine - E_real; // diagonal part of the potential

            auto vpsi = Vtot * psi[inzb][reim][spin*64 + i64][j64]; // potential is diagonal in real-space

            if (imaginary) {
                // V-E has an imaginary part V_Im = -E_imag
                // then explicitly:
                //    Vpsi_Re = V_Re * psi_Re - V_Im * psi_Im
                //    Vpsi_Im = V_Re * psi_Im + V_im * psi_Re
                vpsi += V_imag * psi[inzb][1 - reim][spin*64 + i64][j64]; // ToDo: check the sign again!
            } // imaginary

            if (2 == Noco) { // the other spin component is (1 - spin)
                /*                                                                   */
                /*  this code would be correct if a noco potential had real values,  */
                /*  however, it has 4 components (V_1, V_x, V_y, V_z)                */
                /*                                                                   */
                /*           /  V_1 + V_z   V_x-i*V_y \    / V_dndn   V_x-i*V_y \    */
                /*  V_noco = |                        | =  |                    |    */
                /*           \  V_x+i*V_y   V_1 - V_z /    \ V_x+i*V_y   V_upup /    */
                /*                                                                   */
                /*  however, to avoid memory accesses we store                       */
                /*  these four combinations in Vloc[0..3][i64]                       */
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
                vpsi += Vloc[icube][2][i64] * psi[inzb][reim    ][(1 - spin)*64 + i64][j64]; // V_x
                vpsi += Vloc[icube][3][i64] * psi[inzb][1 - reim][(1 - spin)*64 + i64][j64]  // V_y
                            * (1 - 2*(reim ^ spin)); // complex sign is -1 if (reim != spin)
            } // non-collinear

            Vpsi[inzb][reim][spin*64 + i64][j64] = vpsi; // store
        } // inzb

    } // Potential

    template <typename real_t, int R1C2=2, int Noco=1>
    size_t multiply(
          real_t         (*const ____restrict____ Vpsi)[R1C2][Noco*64][Noco*64] // result
        , real_t   const (*const ____restrict____  psi)[R1C2][Noco*64][Noco*64] // input
        , double   const hgrid=1 // grid spacing, ToDo make it an array of X,Y,Z
    ) {
      
        return 0ul; // total number of floating point operations performed
    } // multiply potential

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

} // namespace green_potential
