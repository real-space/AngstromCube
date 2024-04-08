// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // uint32_t, int32_t, uint16_t, int16_t
#include <cassert> // assert

#include "kinetic_plan.hxx"

#include "green_memory.hxx" // get_memory, free_memory --> used in multiply, ToDo: move into planning phase
#include "green_sparse.hxx" // ::sparse_t<T>
#include "data_view.hxx" // view3D<T>
#include "constants.hxx" // ::pi
#include "recorded_warnings.hxx" // warn
#include "print_tools.hxx" // printf_vector

namespace kinetic_plan {


    int32_t get_inz( // search a column index in a BSR structure
          uint32_t const idx[3] // {ix,iy,iz}-coordinates into iRow_of_coords
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
        , uint32_t const periodicity // periodicity in derivative direction, 0 if not periodic
        , std::vector<int32_t> const target_axes[3] // contain the global indices w.r.t. the box and -1 if outside, not used except for their .size()
        , uint32_t const RowStart[]
        , uint16_t const ColIndex[]
        , view3D<int32_t> const & iRow_of_coords // (Z,Y,X) look-up table: row index of the Green function as a function of internal 3D coordinates, -1:non-existent
        , std::vector<std::vector<bool>> const & sparsity_pattern // memory saving bit-arrays sparsity_pattern[nrhs][idx3]
        , int const echo // =0 // log level
    ) {
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

        int constexpr X=0, Y=1, Z=2;
        // prepare the finite-difference sequence lists
        char const direction = 'x' + dd;
        assert(X == dd || Y == dd || Z == dd);
        uint32_t const num_target_coords[] = {uint32_t(target_axes[X].size()),
                                              uint32_t(target_axes[Y].size()),
                                              uint32_t(target_axes[Z].size())};
        assert(num_target_coords[X] > 0); assert(num_target_coords[Y] > 0); assert(num_target_coords[Z] > 0);
   //   auto const & target_axis = target_axes[dd]; // in derivative direction (not used)

        auto const number_all_target_coords = (num_target_coords[Z])*size_t(num_target_coords[Y])*size_t(num_target_coords[X]);
        int num[3];
        set(num, 3, num_target_coords);
        uint32_t const num_dd = num[dd];
        num[dd] = 1; // replace number of target blocks in derivative direction
        if (echo > 0) std::printf("# FD lists in %c-direction %d %d %d\n", direction, num[X], num[Y], num[Z]);
        simple_stats::Stats<> length_stats;
        uint32_t const nrhs = sparsity_pattern.size();
        auto const max_lists = nrhs*size_t(num[Z])*size_t(num[Y])*size_t(num[X]);
        std::vector<std::vector<int32_t>> list(max_lists);

        size_t ilist{0}, boundary_block_missing{0};
        for (uint32_t irhs = 0; irhs < nrhs; ++irhs) {
//          if (echo > 0) std::printf("# FD list for RHS #%i\n", irhs);
            auto const & sparsity_rhs = sparsity_pattern[irhs];
            assert(number_all_target_coords == sparsity_rhs.size());
            for (uint32_t iz = 0; iz < num[Z]; ++iz) { //
            for (uint32_t iy = 0; iy < num[Y]; ++iy) { //   one of these 3 loops has range == 1
            for (uint32_t ix = 0; ix < num[X]; ++ix) { //
                uint32_t idx[] = {ix, iy, iz}; // non-const
                assert(ilist < max_lists);
                assert(0 == list[ilist].size()); // make sure the list is empty at start
                int32_t last_id{-1}, first_id{-1}; // index in derivative direction of the first and last non-zero entry
                for (uint32_t id{0}; id < num_dd; ++id) { // loop over direction to derive
                    idx[dd] = id; // replace index in the derivate direction
//                  if (echo > 0) std::printf("# FD list for RHS #%i test coordinates %i %i %i\n", irhs, idx[X], idx[Y], idx[Z]);
                    auto const idx3 = index3D(num_target_coords, idx); // flat index
                    if (sparsity_rhs[idx3]) {
                        auto const inz_found = get_inz(idx, iRow_of_coords, RowStart, ColIndex, irhs, 'x' + dd);
                        assert(inz_found >= 0 && "sparsity_pattern[irhs][idx3] inconsistent with iRow_of_coords[iz][iy][ix]");

                        if (0 == list[ilist].size()) {
                            list[ilist].reserve(nhalo + num_dd + nhalo); // makes push_back operation faster
                            list[ilist].resize(nhalo, CUBE_IS_ZERO); // prepend {0, 0, 0, 0}
                            first_id = id; // first nonzero id
                        } // list was empty

                        list[ilist].push_back(inz_found + CUBE_EXISTS);
                        last_id = id; // last nonzero id
                    } // sparsity pattern
                } // id

                auto const list_length = list[ilist].size();
                if (list_length > nhalo) {
                    length_stats.add(list_length - nhalo);
//                  if (echo > 0) std::printf("# FD list of length %d for RHS#%i in %c-direction %i %i %i\n", list_length - nhalo, irhs, direction, idx[X], idx[Y], idx[Z]);
                    // add nhalo end-of-sequence markers
                    for (int ihalo = 0; ihalo < nhalo; ++ihalo) {
                        list[ilist].push_back(CUBE_IS_ZERO); // append {0, 0, 0, 0} to mark the end of the derivative sequence
                    } // ihalo

                    if (periodicity) {
//                      if (echo > 0) std::printf("# FD list of length %d for RHS#%i in %c-direction first=%i last=%i\n", list_length - nhalo, irhs, direction, first_id, last_id);
                        assert(first_id >= 0); assert(last_id >= 0);
                        if(last_id + 1 - first_id != periodicity) error("last_id= %i first_id= %i periodicity= %d RHS#%i xyz={%i %i %i} dd=%c",
                                                                         last_id,    first_id,    periodicity,   irhs,        ix,iy,iz, dd+'x');
                        assert(last_id + 1 - first_id == periodicity); // if the direction is periodic, exactly periodicity adjacent blocks must be nonzero
                        for (int i01{0}; i01 < 2; ++i01) { // i01==0:start and i01==1:end
                            for (int ihalo = 0; ihalo < nhalo; ++ihalo) {
                                int32_t const jdx_dd = i01 ? (last_id - periodicity + ihalo + 1):
                                                            (first_id + periodicity + ihalo - nhalo);
                                if (jdx_dd >= 0 && jdx_dd < num_target_coords[dd]) {
                                    uint32_t jdx[] = {idx[X], idx[Y], idx[Z]}; // non-const
                                    jdx[dd] = jdx_dd;
                                    assert(jdx_dd == jdx[dd] && "safe assign");

                                    auto const jdx3 = index3D(num_target_coords, jdx); // flat index
                                    if (sparsity_rhs[jdx3]) {
                                        auto const jnz_found = get_inz(jdx, iRow_of_coords, RowStart, ColIndex, irhs, 'X' + dd);
                                        if (jnz_found >= 0) {
                                            list[ilist][ihalo + i01*list_length] = CUBE_NEEDS_PHASE*(jnz_found + CUBE_EXISTS);
                                        } // block exists
                                    } else {
                                        boundary_block_missing += (ihalo == (1 - i01)*(nhalo - 1)); // only warn if the one block close to the central regions is not found
                                    } // block should exists
                                } // halo index inside target box
                            } // ihalo
                        } // i01
                    } // periodicity
                    assert(list_length + nhalo == list[ilist].size());

                    ++ilist; // create next list index
                } // list_length > 0

            }}} // ixyz
        } // irhs

        if (boundary_block_missing) {
            warn("in %ld situations halo blocks were not found", boundary_block_missing);
        } // boundary_block_missing

        // store the number of lists
        uint32_t const n_lists = ilist; assert(n_lists == ilist && "too many lists, max. 2^32-1");
        if (echo > 0) std::printf("# %d FD lists for the %c-direction (%.2f %%), lengths in [%g, %g +/- %g, %g]\n",
                            n_lists, direction, n_lists/std::max(max_lists*.01, .01),
                            length_stats.min(), length_stats.mean(), length_stats.dev(), length_stats.max());
        list.resize(n_lists);

        sparse = green_sparse::sparse_t<int32_t>(list, false, "finite_difference_list", echo); // this assignment uses the move assignment

        if (periodicity && FD_range > 4) {
            warn("boundary is periodic in %c-direction, reduce finite-difference range to 4", 'x' + dd);
            FD_range = 4;
        }
        if (FD_range < 1) {
            warn("kinetic energy switched off in %c-direction", 'x' + dd);
        }
        return 0;
    } // finite_difference_plan



    double set_phase(
          double phase[2][2]
        , double const phase_angle // =0
        , char const direction // ='?'
        , int const echo // =0
    ) {
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
            if (std::abs(cs - std::cos(arg)) > 3.4e-16) error("cosine for phase_angle= %g degrees deviates after purification, cos= %g expected %g diff %.1e", phase_angle*360, cs, std::cos(arg), cs - std::cos(arg));
            if (std::abs(sn - std::sin(arg)) > 3.4e-16) error(  "sine for phase_angle= %g degrees deviates after purification, sin= %g expected %g diff %.1e", phase_angle*360, sn, std::sin(arg), sn - std::sin(arg));
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

    void set_phase(
          double phase[3][2][2]
        , double const phase_angles[3] // =nullptr
        , int const echo // =0
    ) {
        for (int dd = 0; dd < 3; ++dd) {
            set_phase(phase[dd], phase_angles ? phase_angles[dd] : 0, 'x' + dd, echo);
        } // dd
    } // set_phase

} // namespace kinetic_plan


    kinetic_plan_t::kinetic_plan_t(
            int16_t & FD_range // side result
          , int const dd // direction of derivative, 0:X, 1:Y, 2:Z
          , uint32_t const periodicity
          , std::vector<int32_t> const target_axes[3]
          , uint32_t const RowStart[]
          , uint16_t const ColIndex[]
          , view3D<int32_t> const & iRow_of_coords // (Z,Y,X) look-up table: row index of the Green function as a function of internal 3D coordinates, -1:non-existent
          , std::vector<std::vector<bool>> const & sparsity_pattern // memory saving bit-arrays sparsity_pattern[nrhs][idx3]
          , int const echo // =0 // log level
      )
        : derivative_direction_(dd)
      {
          auto const stat = kinetic_plan::finite_difference_plan(sparse_, FD_range, // results
                    dd, periodicity, target_axes, RowStart, ColIndex,
                    iRow_of_coords, sparsity_pattern, echo);
          if (stat) warn("failed to set up finite_difference_plan for %c-direction, status= %i", 'x'+dd, int(stat));
      } // constructor


      void kinetic_plan_t::set(
            int const dd
          , double const grid_spacing // =1
          , size_t const nnzbX // =1
          , int const echo // =0 // log level
      ) {
            derivative_direction_ = dd;
            nnzb_ = nnzbX;
            prefactor_ = -0.5/(grid_spacing*grid_spacing);
            char _lists[8] = "list[?]"; _lists[5] = 'x' + dd;
            if (lists_) free_memory(lists_);
            lists_ = get_memory<int32_t const *>(sparse_.nRows(), echo, _lists);
            auto const rowStart = sparse_.rowStart();
            auto const colIndex = sparse_.colIndex();
            for (uint32_t il = 0; il < sparse_.nRows(); ++il) {
                lists_[il] = &colIndex[rowStart[il]];
            } // il
      } // set


      kinetic_plan_t::~kinetic_plan_t() {
          if (lists_) {
              std::printf("# free list for %c-direction\n", 'x'+derivative_direction_);
              free_memory(lists_);
          }
      } // destructor
