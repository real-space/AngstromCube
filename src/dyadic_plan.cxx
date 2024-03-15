// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert

#include "dyadic_plan.hxx"

#include "status.hxx" // status_t
#include "green_memory.hxx" // get_memory, free_memory, dim3, real_t_name
#include "green_sparse.hxx" // ::sparse_t<>
#include "inline_math.hxx" // pow2, pow3
#include "sho_tools.hxx" // ::nSHO, ::n2HO, ::n1HO

    dyadic_plan_t::dyadic_plan_t(int const echo) { // constructor
        if (echo > 0) std::printf("# construct %s\n", __func__);
        // please see construct_dyadic_plan in green_function.hxx for the construction of the dyadic_plan_t
    } // empty and default constructor

    dyadic_plan_t::~dyadic_plan_t() { // destructor
#ifdef DEBUG
        std::printf("# destruct %s\n", __func__);
#endif // DEBUG
        free_memory(AtomImageIndex);
        free_memory(AtomImageStarts);
        free_memory(AtomStarts);
        if (AtomMatrices) for (uint32_t ia = 0; ia < nAtoms; ++ia) free_memory(AtomMatrices[ia]);
        free_memory(AtomMatrices);
        free_memory(grid_spacing);
        free_memory(AtomImagePos);
        free_memory(AtomImageLmax);
        free_memory(AtomImagePhase);
        free_memory(AtomImageShift);
        free_memory(AtomLmax);
        if (sparse_SHOprj) for (int32_t irhs = 0; irhs < nrhs; ++irhs) sparse_SHOprj[irhs].~sparse_t<>();
        free_memory(sparse_SHOprj);
    } // constructor

    status_t dyadic_plan_t::consistency_check() const {
        status_t stat(0);
        assert(nAtomImages >= nAtoms);
        auto const rowStart = sparse_SHOsum.rowStart();
        auto const colIndex = sparse_SHOsum.colIndex();
        stat += (sparse_SHOsum.nRows() != nAtoms);
        if (!rowStart) return stat;
        for (uint32_t ia = 0; ia < nAtoms; ++ia) {
            for (auto bsr = rowStart[ia]; bsr < rowStart[ia + 1]; ++bsr) {
                auto const iai = colIndex[bsr];
                stat += (AtomImageLmax[iai] != AtomLmax[ia]);
            } // bsr
        } // ia
        return stat; // number of errors
    } // consistency_check


    void dyadic_plan_t::update_flop_counts(int const echo) {

        flop_count_SHOgen = 0;
        flop_count_SHOadd = 0;
        {
            auto const iai_of_bsr = sparse_SHOadd.colIndex();
            auto const nnz        = sparse_SHOadd.nNonzeros();
            for (size_t bsr = 0; bsr < nnz; ++bsr) {
                int const lmax = AtomImageLmax[iai_of_bsr[bsr]];
                flop_count_SHOadd += 64 * flop_count_SHOprj_SHOadd(lmax); // 64 threads per block (real version)
                flop_count_SHOgen += 12 * flop_count_Hermite_Gauss(lmax); // 12 threads per block eval the Hermite polynomials
                // this flop count does not account for masked blocks
            } // bsr
        }
//        flop_count_SHOprj = flop_count_SHOadd; // symmetric, both missing factor R1C2 Noco^2

        flop_count_SHOmul = 0;
        for (int ia = 0; ia < nAtoms; ++ia) {
            int const lmax = AtomImageLmax[ia];
            flop_count_SHOmul += pow2(sho_tools::nSHO(lmax));
        } // ia
        flop_count_SHOmul *= 2*nrhs*64; // missing factor R1C2^2 Noco^3

        flop_count_SHOsum = 0;
        if (nAtomImages > nAtoms) {
            for (int iai = 0; iai < nAtomImages; ++iai) {
                int const lmax = AtomImageLmax[iai];
                flop_count_SHOsum += sho_tools::nSHO(lmax);
            } // iai
            flop_count_SHOsum *= 2*nrhs*64; // missing factor R1C2^2 Noco^2
        } else { /* no need to run SHOsum */ }

    } // update_flop_counts

