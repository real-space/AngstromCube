#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <vector> // std::vector<T>

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_sparse.hxx" // ::sparse_t<>

class dyadic_plan_t {
public: // members

    uint32_t* AtomStarts          = nullptr; // [nAtoms + 1]
    int8_t*   AtomLmax            = nullptr; // [nAtoms]
    double**  AtomMatrices        = nullptr; // [nAtoms][2*nc^2] atomic matrices, nc: number of SHO coefficients of this atom
    uint32_t nAtoms               = 0;
    // std::vector<double> AtomSigma;

    uint32_t* AtomImageIndex      = nullptr; // [nAtomImages]
    uint32_t* AtomImageStarts     = nullptr; // [nAtomImages + 1]
    double  (*AtomImagePos)[3+1]  = nullptr; // [nAtomImages][4]
    int8_t*   AtomImageLmax       = nullptr; // [nAtomImages]
    double  (*AtomImagePhase)[4]  = nullptr; // [nAtomImages][4]
    int8_t  (*AtomImageShift)[4]  = nullptr; // [nAtomImages][4]
    uint32_t nAtomImages          = 0;

    double* grid_spacing          = nullptr; // [3+1] hx,hy,hz,rcut/sigma

    int32_t nrhs                  = 0;
    green_sparse::sparse_t<>* sparse_SHOprj = nullptr; // [nrhs]
    green_sparse::sparse_t<>  sparse_SHOadd;
    green_sparse::sparse_t<>  sparse_SHOsum;

    std::vector<int64_t> global_atom_ids; // [nAtoms]
    std::vector<int32_t> global_atom_index;
    std::vector<int32_t> original_atom_index;

    size_t  flop_count_SHOgen = 0,
            flop_count_SHOsum = 0,
            flop_count_SHOmul = 0,
            flop_count_SHOadd = 0;

public: // methods

    dyadic_plan_t(int const echo=0); // constructor
    ~dyadic_plan_t(); // destructor

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
