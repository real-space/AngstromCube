#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int16_t, int8_t
#include <cassert> // assert
#include <complex> // std::complex

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_cuda.hxx" // __global__
#include "inline_math.hxx" // pow2


//// control here, if we use a confinement potential at all?
#define   CONFINEMENT_POTENTIAL


namespace green_potential {

    template <typename real_t, int R1C2=2, int Noco=1>
    void __global__ Potential( // GPU kernel, must be launched with <<< {64, any, 1}, {Noco*64, Noco, R1C2} >>>
#ifdef    HAS_NO_CUDA
          dim3 const & gridDim, dim3 const & blockDim,
#endif // HAS_NO_CUDA
          real_t        (*const __restrict__ Vpsi)[R1C2][Noco*64][Noco*64] // result
        , real_t  const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input Green function
        , double  const (*const *const __restrict__ Vloc)[64] // local potential, Vloc[Noco*Noco][iloc][4*4*4]
        , int32_t const (*const __restrict__ iloc_of_inzb) // translation from inzb to iloc, [inzb]
        , int16_t const (*const __restrict__ target_minus_source)[3+1] // 3D cube shift vector, 4th component unused, [inzb][0:2]
        // , float   const (*const __restrict__ rowCubePos)[3+1] // target cube coords
        // , float   const (*const __restrict__ colCubePos)[3+1] // source cube coords
        , double  const (*const __restrict__ hxyz) // grid spacing in X,Y,Z direction
        , int     const nnzb // number of all cubes to be treated
        , float   const Vconf // prefactor for the confinement potential
        , float   const rconf2 // confinement radius^2 for the confinement potential, negative for no confinement
        , float   const rcut2  // hard cutoff radius^2 for masking --> this modifies the input vector!
        , double  const E_real // real      part of the energy parameter, this could be subtracted from Vloc beforehand
        , double  const E_imag // imaginary part of the energy parameter
    ) {
        assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));

        assert(64      ==  gridDim.x);
        assert(1       ==  gridDim.z);
        assert(Noco*64 == blockDim.x);
        assert(Noco    == blockDim.y);
        assert(R1C2    == blockDim.z);

        bool const imaginary = ((2 == R1C2) && (0 != E_imag));

#ifndef   HAS_NO_CUDA
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

#ifdef    CONFINEMENT_POTENTIAL
        // generate the position difference of grid points (target minus source)
        int const dx = ( i64       & 0x3) - ( j64       & 0x3); // difference in x-direction in [-3, 3]
        int const dy = ((i64 >> 2) & 0x3) - ((j64 >> 2) & 0x3); // difference in y-direction in [-3, 3]
        int const dz = ((i64 >> 4) & 0x3) - ((j64 >> 4) & 0x3); // difference in z-direction in [-3, 3]
        // due to the masking, we ignore the Noco-spin index inside j64
//      std::printf("# block %d thread %d has in-block shifts %d %d %d \n", i64, j64, x, y, z); // debug
#endif // CONFINEMENT_POTENTIAL

        assert(R1C2 >= Noco); // static_assert?

        double const V_imag = E_imag * double(1 - 2*reim);

        for (int inzb = inz0; inzb < nnzb; inzb += gridDim.y) { // grid-stride loop on y-blocks

#ifdef    CONFINEMENT_POTENTIAL
            real_t f{1}, Vconfine{0};
            if (rcut2 >= 0.f) {
                int constexpr n4 = 4;
                auto const s = target_minus_source[inzb]; // shift vectors between target minus source cube, ToDo: we could generate it from rowCubePos - colCubePos (which needs more bandwidth if nCols==1)
                auto const r2 = pow2((s[0]*n4 + dx)*hxyz[0])
                              + pow2((s[1]*n4 + dy)*hxyz[1])
                              + pow2((s[2]*n4 + dz)*hxyz[2]);
                if (r2 > rcut2) {
                    f = 0; // masking (avoid too large numeric values)
                } else if (r2 > rconf2) {
                    // confinement potential: Vconf*(r^2 - rconf^2)^2 for r > rconf i.e. r^2 > rconf^2
                    auto const r2out = real_t(r2 - rconf2);
                    Vconfine = Vconf*pow2(r2out); // quartic confinement potential, V ~ r^4 for large r
//                  std::printf("%g %g\n", r2, Vconfine); // plot
                }
//              std::printf("%.1f %.3f  %d %d %d  %d\n", r2, Vconfine, s[0], s[1], s[2], s[3]); // s[3] == source cube index
            } // rcut^2 >= 0
#else  // CONFINEMENT_POTENTIAL
            real_t constexpr f = 1, Vconfine = 0;
#endif // CONFINEMENT_POTENTIAL

            auto const iloc = iloc_of_inzb[inzb]; // target index for the local potential, can be -1 for non-existing
            double const Vloc_diag = (iloc < 0) ? 0 : Vloc[spin][iloc][i64];

            // gather all real-valued and spin-diagonal contributions
            double const Vtot = Vloc_diag + Vconfine - E_real; // diagonal part of the potential

            auto vpsi = Vtot * psi[inzb][reim][spin*64 + i64][j64]; // non-const, potential is diagonal in real-space

            if (imaginary) {
                // V-E has an imaginary part V_Im = -E_imag
                // then explicitly:
                //    Vpsi_Re = V_Re * psi_Re - V_Im * psi_Im = V_Re * psi_Re + E_imag * psi_Im (reim=0)
                //    Vpsi_Im = V_Re * psi_Im + V_im * psi_Re = V_Re * psi_Im - E_imag * psi_Re (reim=1)
                vpsi += V_imag * psi[inzb][1 - reim][spin*64 + i64][j64];
            } // imaginary

            if (2 == Noco) {
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
                if (iloc >= 0) { // the other spin component is (1 - spin)
                    real_t const cs = (1 - 2*(reim ^ spin)); // complex sign is -1 if (reim != spin)

                    vpsi += Vloc[2][iloc][i64] * psi[inzb][    reim][(1 - spin)*64 + i64][j64];    // V_x
                    vpsi += Vloc[3][iloc][i64] * psi[inzb][1 - reim][(1 - spin)*64 + i64][j64]*cs; // V_y
                } // iloc > 0
            } // non-collinear

            Vpsi[inzb][reim][spin*64 + i64][j64] = f*vpsi; // store

        } // inzb

        } // thread loops and block loop

    } // Potential


    template <typename real_t, int R1C2=2, int Noco=1>
    size_t multiply(
          real_t         (*const __restrict__ Vpsi)[R1C2][Noco*64][Noco*64] // result
        , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input
        , double   const (*const *const __restrict__ Vloc)[64] // local potential, Vloc[Noco*Noco][iloc][4*4*4]
        , int32_t  const (*const __restrict__ vloc_index) // iloc_of_inzb[nnzb]
        , int16_t  const (*const __restrict__ target_minus_source)[3+1] // 3D cube shift vector (target minus source), 4th component unused
        , double   const (*const __restrict__ hxyz) // grid spacing in X,Y,Z direction
        , uint32_t const nnzb // number of all cubes to be treated
        , std::complex<double> const E_param=0 // energy parameter, the real part could be subtracted from Vloc beforehand
        , float    const Vconf=0  // prefactor for the confinement potential
        , float    const rconf2=-1
        , float    const rcut2=9e37 // cutoff radius^2 for the confinement potential, -1: no confinement
        , int const echo=0
    ) {

        if (echo > 11) {
            std::printf("# %s<%s,R1C2=%d,Noco=%d> Vpsi=%p, psi=%p, Vloc=%p, vloc_index=%p, target_minus_source=%p, hxyz=%p, nnzb=%d, Vconf=%g, rcut2=%.f, E=(%g, %g)\n",
                           __func__, (4 == sizeof(real_t))?"float":"double", R1C2, Noco, (void*)Vpsi, (void*)psi,
                           (void*)Vloc, (void*)vloc_index, (void*)target_minus_source, (void*)hxyz, nnzb, Vconf, rcut2, E_param.real(), E_param.imag());
        } // echo

        Potential<real_t,R1C2,Noco>
#ifndef   HAS_NO_CUDA
            <<< dim3(64, 7, 1), dim3(Noco*64, Noco, R1C2) >>> ( // 7=any, maybe find a function for a good choice
#else  // HAS_NO_CUDA
              ( dim3(64, 1, 1), dim3(Noco*64, Noco, R1C2),
#endif // HAS_NO_CUDA
            Vpsi, psi, Vloc, vloc_index, target_minus_source, hxyz, nnzb, Vconf, rconf2, rcut2, E_param.real(), E_param.imag());

        return (  1ul
               +  2ul*(0 != E_param.imag())
               +  4ul*(2 == Noco)
#ifdef    CONFINEMENT_POTENTIAL
               + 10ul*(rcut2 >= 0.f)
#endif // CONFINEMENT_POTENTIAL
               )*nnzb*pow2(64ul*Noco)*R1C2; // total number of floating point operations performed
    } // multiply potential


    template <typename real_t, int R1C2=2, int Noco=1>
    void __global__ Mask( // GPU kernel, must be launched with <<< {64, any, 1}, {Noco*64, Noco, R1C2} >>>
#ifdef    HAS_NO_CUDA
          dim3 const & gridDim, dim3 const & blockDim,
#endif // HAS_NO_CUDA
          real_t        (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input Green function (not real_t const any more due to rcut2 masking)
        , int16_t const (*const __restrict__ target_minus_source)[3+1] // 3D cube shift vector, 4th component unused, [inzb][0:2]
        , double  const (*const __restrict__ hxyz) // grid spacing in X,Y,Z direction
        , int     const nnzb // number of all cubes to be treated
        , float   const rconf2 // confinement radius^2 for the confinement potential, negative for no confinement
        , float   const rcut2  // hard cutoff radius^2 for masking --> this modifies the input vector!
    ) {
        assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));

        assert(64      ==  gridDim.x);
        assert(1       ==  gridDim.z);
        assert(Noco*64 == blockDim.x);
        assert(Noco    == blockDim.y);
        assert(R1C2    == blockDim.z);

        double const denom = (rcut2 > rconf2) ? 1./(rcut2 - rconf2) : 0;

#ifndef   HAS_NO_CUDA
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

        // generate the position difference of grid points (target minus source)
        int const dx = ( i64       & 0x3) - ( j64       & 0x3); // difference in x-direction in [-3, 3]
        int const dy = ((i64 >> 2) & 0x3) - ((j64 >> 2) & 0x3); // difference in y-direction in [-3, 3]
        int const dz = ((i64 >> 4) & 0x3) - ((j64 >> 4) & 0x3); // difference in z-direction in [-3, 3]

        assert(R1C2 >= Noco); // static_assert?

        for (int inzb = inz0; inzb < nnzb; inzb += gridDim.y) { // grid-stride loop on y-blocks

            { // scope
                int constexpr n4 = 4;
                auto const s = target_minus_source[inzb]; // shift vectors between target minus source cube
                double const r2 = pow2((s[0]*n4 + dx)*hxyz[0])
                                + pow2((s[1]*n4 + dy)*hxyz[1])
                                + pow2((s[2]*n4 + dz)*hxyz[2]);
                // masking
                if (r2 > rcut2) {
                    psi[inzb][reim][spin*64 + i64][j64] = 0; // enforce truncation of the Green function
                } else if (r2 > rconf2) {
                    // masking function must go smoothly from 1 to 0 as r2 goes from rconf2 to rcut2:
                    // transform onto the interval x in [0, 1] with x = (r2-rconf2)/(rcut2 - rconf2)
                    auto const x = (r2 - rconf2)*denom;
                    // ansatz:
                    //      f(x) = (x - 1)^2 * (2x + 1) == 2*x^3 - 3*x^2 + 1
                    //                            f'(x) == 6*x^2 - 6*x
                    // conditions:
                    //      f (0) == 1
                    //      f'(0) == 0
                    //      f (1) == 0
                    //      f'(1) == 0
                    real_t const f = x*x*(2*x - 3) + 1;
                    // scale the input vector with this factor f in [0, 1]
                    psi[inzb][reim][spin*64 + i64][j64] *= f;
                }
            } // scope

        } // inzb

        } // thread loops and block loop

    } // Mask


    template <typename real_t, int R1C2=2, int Noco=1>
    size_t multiply_mask( // mask
          real_t         (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input and output
        , int16_t  const (*const __restrict__ target_minus_source)[3+1] // 3D cube shift vector (target minus source), 4th component unused
        , double   const (*const __restrict__ hxyz) // grid spacing in X,Y,Z direction
        , uint32_t const nnzb // number of all cubes to be treated
        , float    const rconf2=-1
        , float    const rcut2=9e37 // cutoff radius^2 for the confinement potential, -1: no confinement
        , int const echo=0
    ) {

        if (echo > 11) {
            std::printf("# %s<%s,R1C2=%d,Noco=%d> psi=%p, target_minus_source=%p, hxyz=%p, nnzb=%d, rconf^2=%.f, rcut^2=%.f\n",
                           __func__, (4 == sizeof(real_t))?"float":"double", R1C2, Noco, (void*)psi,
                           (void*)target_minus_source, (void*)hxyz, nnzb, rconf2, rcut2);
        } // echo

#ifdef    CONFINEMENT_POTENTIAL
        if (rcut2 >= 0.f) {
            Mask<real_t,R1C2,Noco>
#ifndef   HAS_NO_CUDA
                <<< dim3(64, 7, 1), dim3(Noco*64, Noco, R1C2) >>> ( // 7=any, maybe find a function for a good choice
#else  // HAS_NO_CUDA
                  ( dim3(64, 1, 1), dim3(Noco*64, Noco, R1C2),
#endif // HAS_NO_CUDA
                psi, target_minus_source, hxyz, nnzb, rconf2, rcut2);
        }
#endif // CONFINEMENT_POTENTIAL

        return 0; // total number of floating point operations performed
    } // multiply_mask




    status_t all_tests(int const echo=0); // declaration only

} // namespace green_potential
