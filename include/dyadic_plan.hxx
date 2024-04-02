#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // FILE
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <vector> // std::vector<T>
#include <utility> // std::swap

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_sparse.hxx" // ::sparse_t<>
#include "inline_math.hxx" // pow2, pow3
#include "data_view.hxx" // view2D<T>

class dyadic_plan_t {
public: // members

    uint32_t* AtomStarts          = nullptr; // [nAtoms + 1]
    int8_t*   AtomLmax            = nullptr; // [nAtoms]
    double**  AtomMatrices        = nullptr; // [nAtoms][2*nc^2] atomic matrices in GPU memory, nc: number of SHO coefficients of this atom
    uint32_t nAtoms               = 0;
    // std::vector<double> AtomSigma;

    uint32_t* AtomImageIndex      = nullptr; // [nAtomImages]
    uint32_t* AtomImageStarts     = nullptr; // [nAtomImages + 1]
    double  (*AtomImagePos)[3+1]  = nullptr; // [nAtomImages][3+1]
    int8_t*   AtomImageLmax       = nullptr; // [nAtomImages]
    double  (*AtomImagePhase)[4]  = nullptr; // [nAtomImages][4]
    int8_t  (*AtomImageShift)[4]  = nullptr; // [nAtomImages][3+1]
    uint32_t nAtomImages          = 0;

    double* grid_spacing          = nullptr; // [3+1] hx,hy,hz,rcut/sigma

    int32_t nrhs                  = 0;
    green_sparse::sparse_t<>* sparse_SHOprj = nullptr; // [nrhs]
    green_sparse::sparse_t<>  sparse_SHOadd;
    green_sparse::sparse_t<>  sparse_SHOsum;

    std::vector<int64_t> global_atom_ids; // [nAtoms]
    std::vector<int32_t> global_atom_index;
    std::vector<int32_t> original_atom_index;

    view2D<double> AtomMatrices_; // dim1=nAtoms, stride=MPI_MAX(2*nc[ia]^2)

    size_t  flop_count_SHOgen = 0,
            flop_count_SHOsum = 0,
            flop_count_SHOmul = 0,
            flop_count_SHOadd = 0;

public: // methods

    dyadic_plan_t(int const echo=0); // default constructor

    dyadic_plan_t( // constructor
        double const cell[3]
      , int8_t const boundary_condition[3]
      , double const grid_spacing[3]
      , std::vector<double> const & xyzZinso // [natoms*8]
      , uint32_t const nRowsGreen
      , uint32_t const nrhs
      , uint32_t const *const rowStartGreen
      , uint16_t const *const colIndexGreen
      , int16_t const (*internal_target_coords)[3+1]
      , int32_t const global_internal_offset[3]
      , double const r_block_circumscribing_sphere
      , double const max_distance_from_center
      , double const r_trunc
      , int const echo=0 // verbosity
      , int const Noco=1 // 1:collinear spins, 2:Non-collinear
      , FILE* const svg=nullptr // for exporting Scalabe Vector Graphics
    ); // declaration only

    ~dyadic_plan_t(); // destructor

    dyadic_plan_t(dyadic_plan_t const &) = delete; // { std::cout << "A(A&)\n"; } // copy constructor
    dyadic_plan_t(dyadic_plan_t &&) = delete; // { std::cout << "A(A&&)\n"; } // move constructor
    dyadic_plan_t & operator=(dyadic_plan_t const &) = delete; // { std::cout << "A=(A&)\n"; return *this; } // copy assignment

//  dyadic_plan_t & operator=(dyadic_plan_t &&) = delete; // { std::cout << "A=(A&&)\n"; return *this; } // move assignment
    dyadic_plan_t & operator=(dyadic_plan_t && rhs); // custom move assignment operator

    status_t consistency_check() const; // declaration only

    int flop_count_SHOprj_SHOadd(int const L) {
        return 2*(4*4*4 * (L+1) + 4*4 * (((L+1)*(L+2))/2) + 4 * (((L+1)*(L+2)*(L+3))/6));
    } // flop_count_SHOprj_SHOadd

    int flop_count_Hermite_Gauss(int const L) {
        int constexpr flop_exp = 0; // TODO find out how many flop are needed for std::exp
        return 7 + flop_exp + 4*L;
    } // flop_count_Hermite_Gauss

    void update_flop_counts(int const echo=0); // declaration only

    size_t get_flop_count(int const R1C2, int const Noco, int const echo=0) const {
        size_t nops{0};
//        nops += 0*flop_count_SHOgen; // Hermite Gauss functions
//        nops += 0*flop_count_SHOprj*R1C2*pow2(Noco); // projection
//        nops += 0*flop_count_SHOsum*pow2(R1C2)*pow2(Noco); // collect
        nops +=   flop_count_SHOmul*pow2(R1C2)*pow3(Noco); // small matrix multiplication
        nops += 2*flop_count_SHOsum*pow2(R1C2)*pow2(Noco); // broadcast
        nops += 2*flop_count_SHOgen; // Hermite Gauss functions
        nops += 2*flop_count_SHOadd*R1C2*pow2(Noco); // addition
        return nops;
    } // get_flop_count

}; // dyadic_plan_t
