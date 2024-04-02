#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // uint32_t, int32_t, uint16_t, int16_t
#include <cassert> // assert

#include "green_sparse.hxx" // ::sparse_t<T>
#include "data_view.hxx" // view3D<T>

    template <typename uint_t, typename int_t> inline
    size_t index3D(uint_t const n[3], int_t const i[3]) {
        return size_t(i[2]*n[1] + i[1])*n[0] + i[0]; // usual 3D indexing
    } // index3D

namespace kinetic_plan {

    double set_phase(
          double phase[2][2]
        , double const phase_angle=0
        , char const direction='?'
        , int const echo=0
    ); // declaration only

    void set_phase(
          double phase[3][2][2]
        , double const phase_angles[3]=nullptr
        , int const echo=0
    ); // declaration only

} // namespace kinetic_plan


    class kinetic_plan_t {
    public:

        static int constexpr nhalo = 4; // a maximum of 4 blocks (i.e. 16 grid points) is the range of the FD stencil.

        kinetic_plan_t() {} // default constructor

        kinetic_plan_t(
            int16_t & FD_range // side result
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
      ); // constructor

      ~kinetic_plan_t(); // destructor

      void set(
            int const dd
          , double const grid_spacing=1
          , size_t const nnzbX=1
          , int const echo=0 // log level
      ); // declaration only

    public: // ToDo: check which members could be private
        // members
        green_sparse::sparse_t<int32_t> sparse_;
        double prefactor_ = 0; // = -0.5/h^2,  h:grid spacing in X,Y,Z // prefactor of the kinetic energy in Hartree atomic units
        size_t nnzb_; // total number of non-zero blocks (to get the operations count correct)
        int16_t FD_range_ = 8; // 4 or 8
        int32_t const ** lists_ = nullptr; // in device memory
        int derivative_direction_ = -1; // derivative direction

    }; // class kinetic_plan_t

