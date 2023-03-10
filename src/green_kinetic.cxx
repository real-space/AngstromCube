// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cmath> // std::sqrt
#include <algorithm> // std::max
#include <utility> // std::swap //, ::move
#include <vector> // std::vector<T>

#include "green_kinetic.hxx" // index3D, ::nhalo, multiply

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_memory.hxx" // get_memory, free_memory
#include "green_sparse.hxx" // ::sparse_t<>
#include "data_view.hxx" // view3D<T>
#include "inline_math.hxx" // set
#include "simple_stats.hxx" // ::Stats<>
#include "print_tools.hxx" // printf_vector(format, ptr, number [, ...])
#include "constants.hxx" // ::pi

namespace green_kinetic {

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
                                n_lists, direction, n_lists/std::max(max_lists*.01, .01),
                                length_stats.mean(), length_stats.dev(), length_stats.min(), length_stats.max());
          list.resize(n_lists);

          sparse = green_sparse::sparse_t<int32_t>(list, false, "finite_difference_list", echo);

          if (boundary_is_periodic && FD_range > 4) {
              warn("boundary is perodic in %c-direction, reduce range to 4", 'x' + dd);
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

    void __host__ set_phase(
          double phase[3][2][2]
        , double const phase_angles[3] // =nullptr
        , int const echo // =0
    ) {
        for (int dd = 0; dd < 3; ++dd) {
            set_phase(phase[dd], phase_angles ? phase_angles[dd] : 0, 'x' + dd, echo);
        } // dd
    } // set_phase

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <typename real_t, int R1C2=2, int Noco=1>
  status_t test_finite_difference(int const echo=0) {
      size_t const nnzb = 15;
      if (echo > 4) std::printf("\n# %s<%s,R1C2=%d,Noco=%d> start with nnzb= %ld\n", __func__, real_t_name<real_t>(), R1C2, Noco, nnzb);
      auto Tpsi = get_memory<real_t[R1C2][Noco*64][Noco*64]>(nnzb, echo, "Tpsi");
      auto  psi = get_memory<real_t[R1C2][Noco*64][Noco*64]>(nnzb, echo,  "psi");
      set(Tpsi[0][0][0], nnzb*R1C2*pow2(Noco*64ul), real_t(0)); // clear
      set( psi[0][0][0], nnzb*R1C2*pow2(Noco*64ul), real_t(0)); // clear
      auto indx = get_memory<int32_t>(4 + nnzb + 4);
      set(indx, 4 + nnzb + 4, 0); for (size_t i = 0; i < nnzb; ++i) indx[4 + i] = i + 1;
      uint32_t const num[] = {1, 1, 1}; // number of lists, derive in all 3 directions
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

      double const k_x = -2*constants::pi*phase_angle[0]/(nnzb*4); // wave vector
      // set up a periodic wave in the lowest rhs
      int constexpr Re = 0, Im = R1C2 - 1;
      for (size_t inzb = 0; inzb < nnzb; ++inzb) {
          int constexpr iz = 0, iy = 0;
          for (int ix = 0; ix < 4; ++ix) {
              auto const arg = k_x*(inzb*4 + ix)*hgrid[0];
              psi[inzb][Im][(iz*4 + iy)*4 + ix][0] = std::sin(arg);
              psi[inzb][Re][(iz*4 + iy)*4 + ix][0] = std::cos(arg);
          } // ix
      } // inzb
      auto phase = get_memory<double[2][2]>(3, echo, "phase"); // --> TODO move into argument list
      set_phase(phase, phase_angle, 2*echo);
      uint32_t const num100[] = {1, 0, 0}; // only derive in x-direction
      set(Tpsi[0][0][0], nnzb*R1C2*pow2(Noco*64ul), real_t(0)); // clear
      multiply<real_t,R1C2,Noco>(Tpsi, psi, num100, &indx, &indx, &indx, hgrid, nFD4, phase, nnzb, echo);
      free_memory(phase);

      double deviation{0};
      { // scope: check
          auto const Ek = 0.5*pow2(k_x); // kinetic energy
          double dev2{0}, deva{0}, norm2{0};
          if (echo > 9) std::printf("\n## x  Tpsi_Re Tpsi_Im  psi_Re, psi_Im: (wave vector k= %g sqRy, k^2/2= %g Ha)\n", k_x, Ek);
          for (size_t inzb = 0; inzb < nnzb; ++inzb) {
              int const iz = 0, iy = 0;
              for (int ix = 0; ix < 4; ++ix) {
                  int const i64 = (iz*4 + iy)*4 + ix;
                  if (echo > 9) {
                      std::printf("%g  %g %g  %g %g\n", (inzb*4 + ix)*hgrid[0],
                            Tpsi[inzb][Re][i64][0],   Tpsi[inzb][Im][i64][0]*Im,
                          Ek*psi[inzb][Re][i64][0], Ek*psi[inzb][Im][i64][0]*Im);
                  } // echo
                  for (int reim = 0; reim < R1C2; ++reim) {
                      auto const dev = Tpsi[inzb][reim][i64][0] - Ek*psi[inzb][reim][i64][0]; // deviation
                      deva  += std::abs(dev); // measure the absolute deviation
                      dev2  += pow2(dev); // measure the square deviation
                      norm2 += pow2(psi[inzb][reim][i64][0]);
                  } // reim
              } // ix
          } // inzb
          deviation = std::sqrt(dev2);
          if (echo > 5) std::printf("# deviation of a periodic Bloch wave with wave vector k= %g is %.1e (abs)"
                                  " and %.1e (rms) at norm^2 %g\n\n", k_x, deva, deviation, norm2);
      } // scope

      free_memory(Tpsi);
      free_memory(psi);
      free_memory(indx);
      return (deviation > ((8 == sizeof(real_t)) ? 2e-15 : 5e-6));
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
