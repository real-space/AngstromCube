// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cmath> // std::sqrt
#include <algorithm> // std::max
#include <utility> // std::swap
#include <vector> // std::vector<T>

#include "green_kinetic.hxx" // index3D, ::nhalo, ::Laplace_driver

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_memory.hxx" // get_memory, free_memory
#include "green_sparse.hxx" // ::sparse_t<>
#include "kinetic_plan.hxx" // ::set_phase
#include "data_view.hxx" // view3D<T>
#include "inline_math.hxx" // set
#include "simple_stats.hxx" // ::Stats<>
#include "print_tools.hxx" // printf_vector(format, ptr, number [, ...])
#include "constants.hxx" // ::pi
#include "control.hxx" // ::get

namespace green_kinetic {

#ifdef    NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    template <typename real_t, int R1C2=2, int Noco=1>
    size_t test_multiply(
          real_t         (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // result
        , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input
        , uint32_t const num[3] // number of sticks in X,Y,Z direction
        , int32_t  const *const *const __restrict__ list // for the test below we use the same index list for x, y and z
        , double   const hgrid[3] // grid spacing in X,Y,Z
        , int      const FD_range[3] // finite-difference stencil range (in grid points)
        , double   const phase[3][2][2]=nullptr // complex Bloch phase factors
        , size_t   const nnzb=1 // total number of non-zero blocks (to get the operations count correct)
        , int      const echo=0
    ) {
        int nFD[3]; // numbers of non-zero stencil coefficients
        size_t nops{0};
        for (int dd = 0; dd < 3; ++dd) { // derivative direction, serial due to update of Tpsi
            auto const f = -0.5/pow2(hgrid[dd]); // prefactor of the kinetic energy in Hartree atomic units
            int  const stride = 1 << (2*dd); // stride is 4^dd = 1('x'), 4('y') or 16('z')
            nFD[dd] = green_kinetic::Laplace_driver<real_t,R1C2,Noco>(Tpsi, psi, list, f, num[dd], stride, phase ? phase[dd] : nullptr, FD_range[dd]);
            cuCheck( cudaPeekAtLastError() ); // If this fails, probably the Laplace-kernel required too many registers. Try to compile with optimization
            cuCheck( cudaDeviceSynchronize() );
            nops += nnzb*nFD[dd]*R1C2*pow2(Noco*64ul)*2ul; // total number of floating point operations performed
        } // dd
        char const fF = (8 == sizeof(real_t)) ? 'F' : 'f'; // Mflop:float, MFlop:double
        if (echo > 7) std::printf("# green_kinetic::%s nFD= %d %d %d, numbers= %d %d %d, %.3f M%clop\n",
                              __func__, nFD[0], nFD[1], nFD[2], num[0], num[1], num[2], nops*1e-6, fF);
        return nops;
    } // test_multiply (this interface is only used for tests)


    template <typename real_t, int R1C2=2, int Noco=1>
    status_t test_finite_difference(int const echo=0, int const dd=0) {
        uint32_t const nnzb = std::max(2., control::get("green_kinetic.test.nnzb", 35.)); // as test we use a 1D chain of nnzb blocks, TODO: segfault for nnzb=1
        if (echo > 4) std::printf("\n# %s<%s,R1C2=%d,Noco=%d>(derivative_direction=\'%c\'), nnzb= %d\n",
                                        __func__, real_t_name<real_t>(), R1C2, Noco, 'x' + dd, nnzb);
        assert(nnzb > 0); // later we will divide by nnzb

        status_t stat(0);

        // allocate GPU memory in the sparse Green function 4x4x4 data layout
        auto Tpsi = get_memory<real_t[R1C2][Noco*64][Noco*64]>(nnzb, echo, "Tpsi");
        auto  psi = get_memory<real_t[R1C2][Noco*64][Noco*64]>(nnzb, echo,  "psi");
        set(Tpsi[0][0][0], nnzb*R1C2*pow2(Noco*64ul), real_t(0)); // clear
        set( psi[0][0][0], nnzb*R1C2*pow2(Noco*64ul), real_t(0)); // clear

        // create an index list with halo buffers on each side
        auto indx = get_memory<int32_t>(nhalo + nnzb + nhalo, echo, "indx");
        set(indx, nhalo + nnzb + nhalo, CUBE_IS_ZERO); // init all indices as non-existing
        for (size_t i = 0; i < nnzb; ++i) {
            indx[nhalo + i] = i + CUBE_EXISTS; // set the inner indices as existing
        } // i
        uint32_t const num111[] = {1, 1, 1}; // number of lists, derive in all 3 directions
        int const FD_range4[] = {4, 4, 4}, FD_range8[] = {8, 8, 8};
        double const hgrid[] = {1, 1, 1};
        // merely check if memory accesses are allowed
        auto index_list = get_memory<int32_t const*>(1, echo, "index_list");
        index_list[0] = indx;
        cudaDeviceSynchronize();
        test_multiply<real_t,R1C2,Noco>(Tpsi, psi, num111, index_list, hgrid, FD_range8, nullptr, nnzb, echo);
        cudaDeviceSynchronize();
        test_multiply<real_t,R1C2,Noco>(Tpsi, psi, num111, index_list, hgrid, FD_range4, nullptr, nnzb, echo);
        cudaDeviceSynchronize();

        size_t const n1D_all = nnzb*4; // total number of sampling points
        std::vector<double> reference_function[R1C2], reference_result[R1C2];
        for (int reim = 0; reim < R1C2; ++reim) {
            reference_function[reim].resize(n1D_all, 0.0);
            reference_result[reim].resize(n1D_all, 0.0);
        } // reim


        auto phase = get_memory<double[2][2]>(3, echo, "phase");

//      for (int j64 = 0; j64 < 64; ++j64) {
        { int constexpr j64 = 0; // perform all derivatives only in the vector component zero       
        int constexpr Re = 0, Im = R1C2 - 1;

        auto const is_double = int(8 == sizeof(real_t));
        float const threshold[][4] = {{1.6e-5f, 1.5e-5f, 4.5e-3f, 0.f},  // thresholds for float
                                      {5.5e-9f, 4.0e-14f, 8e-12f, 0.f}}; // thresholds for double

        for (int itest = 0; itest < 3; ++itest) { // 3 different tests 0:(FD=4,iso), 1:(FD=8,iso), 2:(FD=4,peri)
//      for (int itest = 2; itest < 3; ++itest) { // only 2:(FD=4,peri)
//      for (int itest = 1; itest < 2; ++itest) { // only 1:(FD=8,iso)

            auto const periodic = (itest > 1);
            auto const FD_range = (itest & 1) ? FD_range8 : FD_range4;
            auto const border   = (itest + 1)*(itest < 2);
                // itest == 0:(FD=4,iso)    we should not compare the first and the last block,         border = 1
                // itest == 1:(FD=8,iso)    we should not compare the first 2 and the last 2 blocks,    border = 2
                // itest == 2:(FD=4,peri)   we can compare all n1D_all entries,                         border = 0

            double kvec{1};
            char const* what = nullptr;
            if (periodic) {

                what = "periodic Bloch wave"; // test the periodic case, only FD_range=4 implemented
                for (int i = 0; i < nhalo; ++i) {
                    indx[nhalo + nnzb + i] = CUBE_NEEDS_PHASE*(i + CUBE_EXISTS); // set the upper halo
                    indx[i] = CUBE_NEEDS_PHASE*(nnzb - nhalo + i + CUBE_EXISTS); // set the lower halo
                } // i
                if (echo > 7) { std::printf("# periodic indices: "); printf_vector(" %d", indx, nhalo + nnzb + nhalo); }
                double const phase_angle[] = {(2 == R1C2) ? .25 : .5, (2 == R1C2) ? -1/3. : -.5, .5}; // in units of 2*pi
                // if we use real numbers only, phase_angle must be half integer, i.e. 0.0:Gamma or 0.5:X-point
                kinetic_plan::set_phase(phase, phase_angle, echo/2);

                kvec = 2*constants::pi*(0 + phase_angle[dd])/(n1D_all*hgrid[dd]); // wave vector in units of inverse grid points
                auto const Ek = 0.5*pow2(kvec); // kinetic energy in Hartree units
                if (echo > 9) std::printf("# %s: wave vector k= %g sqRy, k^2/2= %g Ha\n", __func__, kvec, Ek);
                for (size_t i1D = 0; i1D < n1D_all; ++i1D) {
                    auto const arg = kvec*(i1D + .5)*hgrid[dd];
                    reference_function[Im][i1D] = std::sin(-arg);               // plane wave
                    reference_function[Re][i1D] = std::cos( arg);               // plane wave
                    reference_result[Im][i1D] = Ek*reference_function[Im][i1D]; // plane wave
                    reference_result[Re][i1D] = Ek*reference_function[Re][i1D]; // plane wave
                } // ix

            } else {  // periodic

                what = "localized Gauss-Hermite function"; // test with a localized function set, FD_range=8
                for (int i = 0; i < nhalo; ++i) {
                    indx[nhalo + nnzb + i] = CUBE_IS_ZERO; // set upper halo
                    indx[               i] = CUBE_IS_ZERO; // set lower halo
                } // i
                kvec = 14./(n1D_all*hgrid[dd]); // inverse sigma (to be displayed as k= below)
                auto const Ek = 0.5*pow2(kvec); // kinetic energy in Hartree units
                if (echo > 9) std::printf("# %s: sigma = %g Bohr\n", __func__, 1./kvec);
                for (size_t i1D = 0; i1D < n1D_all; ++i1D) {
                    auto const x = kvec*(i1D - .5*(n1D_all - 1))*hgrid[dd];
                    auto const Gaussian = std::exp(-.5*x*x);
                    reference_function[Im][i1D] = Gaussian*x;               // Gauss-Hermite-function #1
                    reference_function[Re][i1D] = Gaussian;                 // Gauss-Hermite-function #0
                    reference_result[Im][i1D] = Ek*Gaussian*(3. - x*x)*x;   // Gauss-Hermite-function #3
                    reference_result[Re][i1D] = Ek*Gaussian*(1. - x*x);     // Gauss-Hermite-function #2
                } // ix

            } // periodic

            int const stride = (1ul << (2*dd));

            // set up a periodic wave in the lowest rhs j64=0
            for (size_t inzb = 0; inzb < nnzb; ++inzb) {
                for (int i4 = 0; i4 < 4; ++i4) {
                    auto const i64 = i4*stride; // two block indices are zero
                    auto const i1D = inzb*4 + i4;
                    psi[inzb][Im][i64][j64] = reference_function[Im][i1D];
                    psi[inzb][Re][i64][j64] = reference_function[Re][i1D];
                } // i4
            } // inzb
            // Caveat: Treats R1C2==2 without errors, but does not test the upper components if Noco==2

            uint32_t num100[] = {0, 0, 0}; num100[dd] = 1; // only derive in dd-direction
            set(Tpsi[0][0][0], nnzb*R1C2*pow2(Noco*64ul), real_t(0)); // clear
            cudaDeviceSynchronize();
            test_multiply<real_t,R1C2,Noco>(Tpsi, psi, num100, index_list, hgrid, FD_range, periodic ? phase : nullptr, nnzb, echo);
            cudaDeviceSynchronize();

            int failed{0};
            { // scope: compare result of multiply with reference_result
                double dev2{0}, deva{0}, norm2{0};
                if (echo > 9) std::printf("\n## x  Tpsi_Re Tpsi_Im,  reference Tpsi_Re Tpsi_Im,  psi_Re psi_Im:\n"); // plot header
                for (size_t inzb = border; inzb < nnzb - border; ++inzb) {
                    for (int i4 = 0; i4 < 4; ++i4) {
                        auto const i64 = i4*stride; // two block indices are zero
                        auto const i1D = inzb*4 + i4;
                        if (echo > 9) { // plot
                            std::printf("%g  %.15e %g  %.15e %g  %.15e %g\n", i1D*hgrid[0],
                                Tpsi[inzb][Re][i64][j64],   Tpsi[inzb][Im][i64][j64]*Im,
                                reference_result[Re][i1D], reference_result[Im][i1D]*Im,
                                psi[inzb][Re][i64][j64],     psi[inzb][Re][i64][j64]*Im);
                        } // echo
                        for (int reim = 0; reim < R1C2; ++reim) {
                            auto const dev = Tpsi[inzb][reim][i64][j64] - reference_result[reim][i1D]; // deviation
                            deva  += std::abs(dev); // measure the absolute deviation
                            dev2  += pow2(dev); // measure the square deviation
                            norm2 += pow2(reference_result[reim][i1D]);
                        } // reim
                    } // i4
                } // inzb
                auto const dev = std::sqrt(dev2/norm2);
                deva          /= std::sqrt(norm2);
                auto const thres = threshold[is_double][itest];
                failed = (dev > thres);
                if (echo > 5) std::printf("# <%s>deviation of a %s in %c-direction with k= %g is %.1e (abs) and %.1e (rms), norm^2= %g, FD=%d, rms_threshold=%.1e %s\n",
                    is_double?"double":"float", what, 'x'+dd, kvec, deva, dev, norm2, FD_range[0], thres, failed?"FAILED":"ok");
            } // scope
            stat += failed;

        } // itest
        } // j64

        free_memory(phase);
        free_memory(Tpsi);
        free_memory(psi);
        free_memory(indx);
        free_memory(index_list);
        return stat;
    } // test_finite_difference

  status_t test_set_phase(int const echo=0) {
      double phase[2][2], maxdev{0};
      for (int iangle = -180; iangle <= 180; iangle += 5) {
          auto const dev = kinetic_plan::set_phase(phase, iangle/360., 't', echo);
          maxdev = std::max(maxdev, dev);
      } // iangle
      return (maxdev > 3e-16);
  } // test_set_phase

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_set_phase(echo);
      for (int dd = 0; dd < 3; ++dd) {
          stat += test_finite_difference<float ,1,1>(echo, dd);
          stat += test_finite_difference<float ,2,1>(echo, dd);
          stat += test_finite_difference<float ,2,2>(echo, dd);
          stat += test_finite_difference<double,1,1>(echo, dd);
          stat += test_finite_difference<double,2,1>(echo, dd);
          stat += test_finite_difference<double,2,2>(echo, dd);
      } // dd
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_kinetic
