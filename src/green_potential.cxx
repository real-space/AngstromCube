// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <complex> // std::complex

#include "green_potential.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_memory.hxx" // dim3, get_memory, free_memory
#include "inline_math.hxx" // pow2, set
#include "green_parallel.hxx" // ::potential_exchange
#include "global_coordinates.hxx" // ::get
#include "recorded_warnings.hxx" // error
#include "print_tools.hxx" // printf_vector

namespace green_potential {

#ifdef    NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    template <typename real_t, int R1C2=2, int Noco=1>
    status_t test_multiply(
          double   const (*const *const __restrict__ Vloc)[64] // local potential, Vloc[Noco*Noco][iloc][4*4*4]
        , int32_t  const (*const __restrict__ vloc_index) // iloc_of_inzb[nnzb]
        , int16_t  const (*const __restrict__ shift)[3+1] // 3D block shift vector (target minus source), 4th component unused
        , double   const (*const __restrict__ hxyz) // grid spacing in X,Y,Z direction
        , uint32_t const nnzb=1
        , int const echo=0
    ) {
        auto psi   = get_memory<real_t[R1C2][Noco*64][Noco*64]>(nnzb, echo, "psi");
        auto Vpsi  = get_memory<real_t[R1C2][Noco*64][Noco*64]>(nnzb, echo, "Vpsi");

        if (echo > 5) std::printf("# %s<real_t=%s,R1C2=%d,Noco=%d>(Vpsi=%p, psi=%p, ...)\n",
              __func__, (8 == sizeof(real_t)) ? "double" : "float", R1C2, Noco, (void*)Vpsi, (void*)psi);
        multiply<real_t,R1C2,Noco>(Vpsi, psi, Vloc, vloc_index, shift, hxyz, nnzb);

        free_memory(Vpsi);
        free_memory(psi);
        return 0;
    } // test_multiply

    status_t test_multiply(int const echo=0, int const Noco=2) {
        status_t stat(0);
        uint32_t const nnzb = 1;
        auto Vloc = get_memory<double(*)[64]>(Noco*Noco, echo, "Vloc");
        for (int mag = 0; mag < Noco*Noco; ++mag) Vloc[mag] = get_memory<double[64]>(1, echo, "Vloc[mag]");
        auto vloc_index = get_memory<int32_t>(1, echo, "vloc_index"); vloc_index[0] = 0;
        auto shift = get_memory<int16_t[3+1]>(1, echo, "shift");      set(shift[0], 3+1, int16_t(0));
        auto hxyz = get_memory<double>(3+1, echo, "hxyz");            set(hxyz, 3+1, 1.);

        stat += test_multiply<float ,1,1>(Vloc, vloc_index, shift, hxyz, nnzb, echo);
        stat += test_multiply<float ,2,1>(Vloc, vloc_index, shift, hxyz, nnzb, echo);
        stat += test_multiply<float ,2,2>(Vloc, vloc_index, shift, hxyz, nnzb, echo);
        stat += test_multiply<double,1,1>(Vloc, vloc_index, shift, hxyz, nnzb, echo);
        stat += test_multiply<double,2,1>(Vloc, vloc_index, shift, hxyz, nnzb, echo);
        stat += test_multiply<double,2,2>(Vloc, vloc_index, shift, hxyz, nnzb, echo);

        free_memory(hxyz);
        free_memory(shift);
        free_memory(vloc_index);
        for (int mag = 0; mag < Noco*Noco; ++mag) free_memory(Vloc[mag]);
        free_memory(Vloc);
        return stat;
    } // test_multiply

    status_t all_tests(int const echo) {
        status_t stat(0);
        stat += test_multiply(echo);
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_potential
