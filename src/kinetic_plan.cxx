// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // uint32_t, int32_t, uint16_t, int16_t
#include <cassert> // assert

#include "kinetic_plan.hxx"

#include "green_memory.hxx" // get_memory, free_memory --> used in multiply, ToDo: move into planning phase
#include "green_sparse.hxx" // ::sparse_t<T>
#include "data_view.hxx" // view3D<T>

#if 0
    kinetic_plan_t::kinetic_plan_t(
            green_sparse::sparse_t<int32_t> & sparse // result
          , int16_t & FD_range // side result
          , int const dd // direction of derivative, 0:X, 1:Y, 2:Z
          , bool const boundary_is_periodic
          , uint32_t const num_target_coords[3]
          , uint32_t const RowStart[]
          , uint16_t const ColIndex[]
          , view3D<int32_t> const & iRow_of_coords // (Z,Y,X) look-up table: row index of the Green function as a function of internal 3D coordinates, -1:non-existent
          , std::vector<bool> const sparsity_pattern[] // memory saving bit-arrays sparsity_pattern[irhs][idx3]
          , unsigned const nrhs=1 // number of right hand sides
          , double const grid_spacing=1 // grid spacing in derivative direction
          , int const echo=0 // log level
      )
        : derivative_direction(dd)
      {
          auto const stat = finite_difference_plan(sparse, FD_range, // results
                    dd, boundary_is_periodic, num_target_coords, RowStart, ColIndex,
                    iRow_of_coords, sparsity_pattern, nrhs, echo);
          if (0 == stat) {
              prefactor = -0.5/(grid_spacing*grid_spacing);
              lists = get_memory<int32_t const *>(sparse.nRows(), echo, "lists[dd]");
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
            derivative_direction = dd;
            nnzb = nnzbX;
            prefactor = -0.5/(grid_spacing*grid_spacing);
            char lists_[8] = "list[?]"; lists_[5] = 'x' + dd;
            lists = get_memory<int32_t const *>(sparse.nRows(), echo, lists_);
            auto const rowStart = sparse.rowStart();
            auto const colIndex = sparse.colIndex();
            for (uint32_t il = 0; il < sparse.nRows(); ++il) {
                lists[il] = &colIndex[rowStart[il]];
            } // il
      } // set

      kinetic_plan_t::~kinetic_plan_t() {
          if (lists) {
              std::printf("# free list for %c-direction\n", 'x'+derivative_direction);
              free_memory(lists);
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

