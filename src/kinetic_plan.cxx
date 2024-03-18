// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // uint32_t, int32_t, uint16_t, int16_t
#include <cassert> // assert

#include "kinetic_plan.hxx"

#include "green_memory.hxx" // get_memory, free_memory --> used in multiply, ToDo: move into planning phase
#include "green_sparse.hxx" // ::sparse_t<T>
#include "data_view.hxx" // view3D<T>

namespace kinetic_plan {


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

    int constexpr nhalo = kinetic_plan_t::nhalo;

    int32_t constexpr CUBE_EXISTS = 1;
    int32_t constexpr CUBE_IS_ZERO = 0;
    int32_t constexpr CUBE_NEEDS_PHASE = -1;


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
        if (echo > 0) std::printf("# %d FD lists for the %c-direction (%.2f %%), lengths in [%g, %g +/- %g, %g]\n",
                            n_lists, direction, n_lists/std::max(max_lists*.01, .01),
                            length_stats.min(), length_stats.mean(), length_stats.dev(), length_stats.max());
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

} // namespace kinetic_plan

#if 0
    kinetic_plan_t::kinetic_plan_t(
            int16_t & FD_range // side result
          , int const dd // direction of derivative, 0:X, 1:Y, 2:Z
          , bool const boundary_is_periodic
          , uint32_t const num_target_coords[3]
          , uint32_t const RowStart[]
          , uint16_t const ColIndex[]
          , view3D<int32_t> const & iRow_of_coords // (Z,Y,X) look-up table: row index of the Green function as a function of internal 3D coordinates, -1:non-existent
          , std::vector<bool> const sparsity_pattern[] // memory saving bit-arrays sparsity_pattern[irhs][idx3]
          , unsigned const nrhs // =1 // number of right hand sides
          , double const grid_spacing // =1 // grid spacing in derivative direction
          , int const echo // =0 // log level
      )
        : derivative_direction_(dd)
      {
          auto const stat = finite_difference_plan(sparse_, FD_range, // results
                    dd, boundary_is_periodic, num_target_coords, RowStart, ColIndex,
                    iRow_of_coords, sparsity_pattern, nrhs, echo);
          if (0 == stat) {
              prefactor_ = -0.5/(grid_spacing*grid_spacing);
              lists_ = get_memory<int32_t const *>(sparse_.nRows(), echo, "lists[dd]");
          } else {
              if (echo > 2) std::printf("# failed to set up finite_difference_plan for %c-direction, status= %i\n", 'x' + dd, int(stat));
          }
      } // constructor
#endif // does not work because the copy assignment operator of sparse is deleted

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


        // template <typename real_t, int R1C2=2, int Noco=1>
        // size_t kinetic_plan_t::multiply(
        //       real_t         (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // result
        //     , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input
        //     , double   const phase[2][2]=nullptr // complex Bloch phase factors
        //     , int      const echo=0
        // ) const { // members of the kinetic_plan_t are not changed
        //     int  const stride = 1 << (2*derivative_direction); // 4^dd: X:1, Y:4, Z:16
        //     auto const nFD = Laplace_driver<real_t,R1C2,Noco>(Tpsi, psi, lists, prefactor, sparse.nRows(), stride, phase, FD_range);
        //     size_t const nops = nnzb*nFD*R1C2*pow2(Noco*64ul)*2ul;
        //     if (echo > 7) {
        //         char const fF = (8 == sizeof(real_t)) ? 'F' : 'f'; // Mflop:float, MFlop:double
        //         std::printf("# green_kinetic::%s dd=\'%c\', nFD= %d, number= %d, %.3f M%clop\n",
        //             __func__, 'x' + derivative_direction, nFD, sparse.nRows(), nops*1e-6, fF);
        //     } // echo
        //     return nops; // return the total number of floating point operations performed
        // } // multiply (kinetic energy operator)

