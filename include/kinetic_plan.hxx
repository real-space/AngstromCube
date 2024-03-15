#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // uint32_t, int32_t, uint16_t, int16_t
#include <cassert> // assert

#include "green_sparse.hxx" // ::sparse_t<T>
// #include "data_view.hxx" // view3D<T>

  // ToDo: move to green_utils.hxx or similar
  template <typename uint_t, typename int_t> inline
  size_t index3D(uint_t const n[3], int_t const i[3]) {
      // usual 3D indexing
      return size_t(i[2]*n[1] + i[1])*n[0] + i[0];
  } // index3D


    class kinetic_plan_t {
    public:

        kinetic_plan_t() {} // default constructor
#if 0
        kinetic_plan_t(
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
      );
#endif // does not work because the copy assignment operator of sparse is deleted
      ~kinetic_plan_t(); // destructor

      void set(
            int const dd
          , double const grid_spacing=1
          , size_t const nnzbX=1
          , int const echo=0 // log level
      );

    public: // ToDo: check which members could be private
        // members
        green_sparse::sparse_t<int32_t> sparse;
        double prefactor = 0; // = -0.5/h^2,  h:grid spacing in X,Y,Z // prefactor of the kinetic energy in Hartree atomic units
        size_t nnzb; // total number of non-zero blocks (to get the operations count correct)
        int16_t FD_range = 8; // 4 or 8
        int32_t const ** lists = nullptr; // in device memory
        int derivative_direction = -1; // derivative direction

    public:

        // template <typename real_t, int R1C2=2, int Noco=1>
        // size_t multiply(
        //       real_t         (*const __restrict__ Tpsi)[R1C2][Noco*64][Noco*64] // result
        //     , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input
        //     , double   const phase[2][2]=nullptr // complex Bloch phase factors
        //     , int      const echo=0
        // ) const;

    }; // class kinetic_plan_t

