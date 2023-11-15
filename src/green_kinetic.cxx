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
#include "data_view.hxx" // view3D<T>
#include "inline_math.hxx" // set
#include "simple_stats.hxx" // ::Stats<>
#include "print_tools.hxx" // printf_vector(format, ptr, number [, ...])
#include "constants.hxx" // ::pi

namespace green_kinetic {

    int32_t get_inz( // search a column index in a BSR structure
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

        for (auto inz = RowStart[iRow]; inz < RowStart[iRow + 1]; ++inz) {
            if (ColIndex[inz] == irhs) return inz; // ToDo: bisection search would be smarter...
        } // search
        return -1; // not found
    } // get_inz


    status_t finite_difference_plan(
            green_sparse::sparse_t<int32_t> & sparse // result
          , int16_t & FD_range // side result
          , int const dd // direction of derivative
          , bool const boundary_is_periodic
          , uint32_t const num_target_coords[3]
          , uint32_t const RowStart[]
          , uint16_t const ColIndex[]
          , view3D<int32_t> const & iRow_of_coords // (Z,Y,X) look-up table: row index of the Green function as a function of internal 3D coordinates, -1:non-existent
          , std::vector<bool> const sparsity_pattern[] // memory saving bit-arrays sparsity_pattern[irhs][idx3]
          , unsigned const nrhs // =1 // number of right hand sides
          , int const echo // =0 // log level
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
        //    list[0] ==    { 0  0  0  0  1  2  3  4  5  0  0  0  0}
        //    list[1] == { 0  0  0  0  6  7  8  9 10 11 12  0  0  0  0}
        //    list[2] == { 0  0  0  0 13 14 15 16 17 18 19 20  0  0  0  0}
        //    list[3] == { 0  0  0  0 21 22 23 24 25 26 27 28  0  0  0  0}
        //    list[4] == { 0  0  0  0 29 30 31 32 33 34 35  0  0  0  0}
        //    list[5] ==    { 0  0  0  0 36 37 38 39 40  0  0  0  0}
        //
        //  8 y-lists:
        //    list[0] ==    { 0  0  0  0  6 13 21 29  0  0  0  0}
        //    list[1] == { 0  0  0  0  1  7 14 22 30 36  0  0  0  0}
        //    list[2] == { 0  0  0  0  2  8 15 23 31 37  0  0  0  0}
        //    list[3] == { 0  0  0  0  3  9 16 24 32 38  0  0  0  0}
        //    list[4] == { 0  0  0  0  4 10 17 25 33 39  0  0  0  0}
        //    list[5] == { 0  0  0  0  5 11 18 26 34 40  0  0  0  0}
        //    list[6] ==    { 0  0  0  0 12 19 27 35  0  0  0  0}
        //    list[7] ==       { 0  0  0  0 20 28  0  0  0  0}
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
                              list[ilist].resize(nhalo, CUBE_IS_ZERO); // prepend {0, 0, 0, 0}
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
                                                  list[ilist][0 + ihalo] = CUBE_NEEDS_PHASE*(jnz_found + CUBE_EXISTS);
                                              } // block exists
                                          } // block exists
                                      } // ihalo
                                  } // leftmost boundary
                              } // boundary_is_periodic
                              // =====================================================================================================
                              // =====================================================================================================
                          } // list was empty

                          list[ilist].push_back(inz_found + CUBE_EXISTS);
                          last_id = id;
                      } // sparsity pattern
                  } // id
                  int const list_length = list[ilist].size();
                  if (list_length > nhalo) {
                      length_stats.add(list_length - nhalo);
//                    if (echo > 0) std::printf("# FD list of length %d for the %c-direction %i %i %i\n", list_length, direction, idx[X], idx[Y], idx[Z]);
                      // add nhalo end-of-sequence markers
                      for (int ihalo = 0; ihalo < nhalo; ++ihalo) {
                          list[ilist].push_back(CUBE_IS_ZERO); // append {0, 0, 0, 0} to mark the end of the derivative sequence
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
                                                  list[ilist][list_length + ihalo] = CUBE_NEEDS_PHASE*(jnz_found + CUBE_EXISTS);
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
                                n_lists, direction, n_lists/std::max(max_lists*.01, .01),
                                length_stats.mean(), length_stats.dev(), length_stats.min(), length_stats.max());
          list.resize(n_lists);

          sparse = green_sparse::sparse_t<int32_t>(list, false, "finite_difference_list", echo);

          if (boundary_is_periodic && FD_range > 4) {
              warn("boundary is periodic in %c-direction, reduce finite-difference range to 4", 'x' + dd);
              FD_range = 4;
          }
          if (FD_range < 1) {
              warn("kinetic energy switched off in %c-direction", 'x' + dd);
          }
          return 0;
      } // finite_difference_plan

    double __host__ set_phase(double phase[2][2], double const phase_angle=0, char const direction='?', int const echo=0) {
        // phase_angle in units of 2*pi
        // phase_angle ==     0.00 --> phase = 1
        // phase_angle == +/- 0.25 --> phase = +/-i
        // phase_angle == +/- 0.50 --> phase = -1
        // phase_angle == +/- 0.75 --> phase = -/+i
        assert(std::abs(phase_angle) < 0.503); // 181 degrees are accepted for test purposes


        double const arg = 2*constants::pi*phase_angle;
        double cs = std::cos(arg), sn = std::sin(arg);

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
        double dev{0};
// #ifdef    DEBUG
        if ('\0' != *corrected) {
            if (std::abs(cs - std::cos(arg)) > 2.8e-16) error("cosine for phase_angle= %g degrees deviates after purification, cos= %g expected %g", phase_angle*360, cs, std::cos(arg));
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

    void __host__ set_phase(
          double phase[3][2][2]
        , double const phase_angles[3] // =nullptr
        , int const echo // =0
    ) {
        for (int dd = 0; dd < 3; ++dd) {
            set_phase(phase[dd], phase_angles ? phase_angles[dd] : 0, 'x' + dd, echo);
        } // dd
    } // set_phase

#ifdef    NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    template <typename real_t, int R1C2=2, int Noco=1>
    size_t multiply(
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
    } // multiply (this interface is only used for tests)


    template <typename real_t, int R1C2=2, int Noco=1>
    status_t test_finite_difference(int const echo=0, int const dd=0) {
        int const nnzb = 35; // as test we use a 1D chain of nnzb blocks
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
        multiply<real_t,R1C2,Noco>(Tpsi, psi, num111, index_list, hgrid, FD_range8, nullptr, nnzb, echo);
        cudaDeviceSynchronize();
        multiply<real_t,R1C2,Noco>(Tpsi, psi, num111, index_list, hgrid, FD_range4, nullptr, nnzb, echo);
        cudaDeviceSynchronize();

        size_t const n1D_all = nnzb*4; // total number of sampling points
        std::vector<double> reference_function[R1C2], reference_result[R1C2];
        for (int reim = 0; reim < R1C2; ++reim) {
            reference_function[reim].resize(n1D_all, 0.0);
            reference_result[reim].resize(n1D_all, 0.0);
        } // reim

        int constexpr Re = 0, Im = R1C2 - 1, j64 = 0;

        auto phase = get_memory<double[2][2]>(3, echo, "phase");

        auto const is_double = (8 == sizeof(real_t));
        float const threshold[][4] = {{1.6e-5, 1.5e-5, 4.5e-3, 0}, // thresholds for float
                                      {5.5e-9, 4.0e-14, 8e-12, 0}}; // thresholds for double
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
                set_phase(phase, phase_angle, echo/2);

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
                kvec = 14./(n1D_all*hgrid[dd]); // inverse sigma
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
            // Caveat: Treats R1C2==2 correctly, but does not test the upper components if Noco==2

            uint32_t num100[] = {0, 0, 0}; num100[dd] = 1; // only derive in dd-direction
            set(Tpsi[0][0][0], nnzb*R1C2*pow2(Noco*64ul), real_t(0)); // clear
            cudaDeviceSynchronize();
            multiply<real_t,R1C2,Noco>(Tpsi, psi, num100, index_list, hgrid, FD_range, periodic ? phase : nullptr, nnzb, echo);
            cudaDeviceSynchronize();

            int failed{0};
            { // scope: compare result of multiply with reference_result
                double dev2{0}, deva{0}, norm2{0};
                if (echo > 9) std::printf("\n## x  Tpsi_Re Tpsi_Im  psi_Re, psi_Im:\n"); // plot header
                for (size_t inzb = border; inzb < nnzb - border; ++inzb) {
                    for (int i4 = 0; i4 < 4; ++i4) {
                        auto const i64 = i4*stride; // two block indices are zero
                        auto const i1D = inzb*4 + i4;
                        if (echo > 9) { // plot
                            std::printf("%g  %.15e %g  %.15e %g\n", i1D*hgrid[0],
                                Tpsi[inzb][Re][i64][j64],  Tpsi[inzb][Im][i64][j64]*Im,
                                reference_result[Re][i1D], reference_result[Im][i1D]*Im);
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
                failed = (dev > threshold[is_double][itest]);
                if (echo > 5) std::printf("# deviation of a %s with k= %g is %.1e (abs) and %.1e (rms), norm^2= %g, FD=%d, threshold=%.1e %s\n",
                                             what, kvec, deva, dev, norm2, FD_range[0], threshold[is_double][itest], failed?"FAILED":"ok");
            } // scope
            stat += failed;

        } // itest

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
          auto const dev = set_phase(phase, iangle/360., 't', echo);
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
