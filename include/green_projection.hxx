#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cstdio> // std::printf
#include <cmath> // std::exp
#include <vector> // std::vector<T>
#include <complex> // std::complex<T>

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "sho_tools.hxx" // ::nSHO, ::n2HO, ::n1HO
#include "green_parallel.hxx" // ::rank, ::size, ::dyadic_exchange

#include "control.hxx" // ::get

namespace green_projection {

    inline float Hermite_polynomials_1D( // compare green_dyadic::Hermite_polynomials_1D
          double (*const __restrict__ H1D)[3][4] // result H1D[nu][dir][i4]
        , float  (*const __restrict__ xi_squared)[4] // distance^2 xi_squared[dir][i4]
        , int    const lmax
        , double const xyza[3+1] // atom image position in [0],[1],[2], sigma^{-1/2} in [3]
        , float  const xyzc[3]   // position of the target block (in units of 4*grid spacing)
        , double const hxyz[3+1] // grid spacings in [0],[1],[2], projection radius in [3]
    ) {
        double const R2_projection = pow2(double(hxyz[3])); // projection radius can be controlled from outside

        if (lmax < 0) return R2_projection; // as float

        double const sigma_inverse = xyza[3]*xyza[3]; // sigma_inverse = inverse of the Gaussian width sigma
        for (int idir = 0; idir < 3; ++idir) {
            double const grid_spacing  = hxyz[idir]; // load a float and convert to double
            double const cube_position = xyzc[idir]; // position of the target block, lower left front corner
            double const atom_position = xyza[idir]; // load atomic coordinate
            for (int i4 = 0; i4 < 4; ++i4) {

                double const xi = sigma_inverse*(grid_spacing*(cube_position*4 + i4 + 0.5) - atom_position); // center of the grid cell // 6 flop
                double const xi2 = xi*xi; // 1 flop
                xi_squared[idir][i4] = xi2; // export the square of the distance in units of sigma
                double const H0 = (xi2 < R2_projection) ? std::exp(-0.5*xi2) : 0; // Gaussian envelope function // ? flop

                H1D[0][idir][i4] = H0; // store Gaussian H_0

                double Hnp1, Hn{H0}, Hnm1{0};
    //          std::printf("# %s idir=%d xi= %g H= %g", __func__, idir, xi, H0);
                for (int nu = 0; nu < lmax; ++nu) {
                    double const nuhalf = 0.5 * nu; // nu/2. // 1 flop
                    // use the two step recurrence relation to create the 1D Hermite polynomials
                    // H1D[nu+1][idir][i4] = xi * H1D[nu][idir][i4] - nu/2. * H1D[nu-1][idir][i4];
                    Hnp1 = xi * Hn - nuhalf * Hnm1; // see snippets/3Hermite.F90 for the derivation of this // 3 flop
                    H1D[nu + 1][idir][i4] = Hnp1; // store
    //              std::printf(" %g", Hnp1);
                    // rotate registers for the next iteration
                    Hnm1 = Hn; Hn = Hnp1; // ordering is important here
                } // nu
    //          std::printf("\n");

                // Warning!
                // The Hermite polynomials are not normalized. To normalize w.r.t. the 2-norm,
                // we must divide by sqrt{sqrt{pi}*sigma*factorial(nu)} in an additional loop.
                // However, this is not necessary for their orthogonality

            } // i4
        } // idir

        return R2_projection;
    } // Hermite_polynomials_1D

    template <int R1C2=2, int Noco=1>
    std::vector<std::vector<double>> SHOprj( // compare green_dyadic::SHOprj
          double   const (*const __restrict__ pGreen)[R1C2][Noco][Noco*64] // input: projected Green function, layout[natomcoeffs*nrhs][R1C2][Noco][Noco*64]
        , double   const (*const __restrict__ AtomImagePos)[3+1] // atomic positions [0],[1],[2], decay parameter [3]
        , int8_t   const (*const __restrict__ AtomImageLmax) // SHO basis size [iai]
        , uint32_t const (*const __restrict__ AtomImageStarts) // prefix sum over nSHO(AtomLmax[:])
        , uint32_t const (*const __restrict__ AtomImageIndex) // [nAtomImages]
        , double   const (*const __restrict__ AtomImagePhase)[4] // [nAtomImages][4]
        , uint32_t const nAtomImages
        , int8_t   const (*const __restrict__ AtomLmax)
        , uint32_t const nAtoms
        , float    const (*const __restrict__ colCubePos)[3+1] // only [0],[1],[2] used, CubePos[irhs][0:3] of RHS cubes
        , double   const (*const __restrict__ hGrid) // grid spacings in [0],[1],[2], projection radius in [3]
        , uint32_t const nrhs
        , int      const echo=0 // log-level
    ) {
        if (echo > 0) std::printf("# %s for %d atoms\n", __func__, nAtoms);
        std::vector<std::vector<double>> pGp(nAtoms); // result: projection coefficients, layout[nAtoms][nSHO^2*Noco*Noco]
        for (int iatom = 0; iatom < nAtoms; ++iatom) {
            int const lmax = AtomLmax[iatom];
            int const nSHO = sho_tools::nSHO(lmax);
            pGp[iatom].resize(nSHO*nSHO*Noco*Noco, 0.0); // init
        } // iatom

        typedef std::complex<double> complex_t;
        complex_t constexpr zero = complex_t(0, 0);

        for (int iai = 0; iai < nAtomImages; ++iai) {

        int const lmax = AtomLmax[iai];
        int constexpr Lmax = 7;
        if (lmax > Lmax) std::printf("# %s Error: lmax= %d but max. Lmax= %d, iai=%d\n", __func__, lmax, Lmax, iai);
        assert(lmax <= Lmax);

        auto const a0 = AtomImageStarts[iai];
        int const nSHO = sho_tools::nSHO(lmax);
        int const n2HO = sho_tools::n2HO(lmax);

        std::vector<complex_t> piGp(nSHO*nSHO*Noco*Noco, zero);

        double H1D[1 + Lmax][3][4]; // non-normalized orthogonal 1-dimensional Hermite Gauss functions
        float     xi_squared[3][4]; // distance along one Cartesian direction squared

        for (int spin = 0; spin < Noco; ++spin) {
        for (int spjn = 0; spjn < Noco; ++spjn) {
        for (int jsho = 0; jsho < nSHO; ++jsho) {

        std::vector<complex_t> czyx(nSHO, zero);

        // in this loop over RHS we accumulate atom projection coefficients over many cubes, however, we also need to MPI-accumulate later
        for (uint32_t irhs = 0; irhs < nrhs; ++irhs) {

            auto const R2_proj = Hermite_polynomials_1D(H1D, xi_squared, lmax, AtomImagePos[iai], colCubePos[irhs], hGrid);

            for (int z4 = 0; z4 < 4; ++z4) { // loop over real-space grid in z-direction
                auto const d2z = xi_squared[2][z4] - R2_proj;
                if (d2z < 0) {
                    std::vector<complex_t> byx(n2HO, zero);
                    for (int y4 = 0; y4 < 4; ++y4) { // loop over real-space grid in y-direction
                        auto const d2yz = xi_squared[1][y4] + d2z;
                        if (d2yz < 0) {
                            std::vector<complex_t> ax(lmax + 1, zero);
                            for (int x4 = 0; x4 < 4; ++x4) { // loop over real-space grid in x-direction
                                auto const d2xyz = xi_squared[0][x4] + d2yz;
                                if (d2xyz < 0) {
                                    int const xyz = (z4*4 + y4)*4 + x4;
                                    int constexpr RealPart = 0, ImagPart= R1C2 - 1;
                                    auto const pG = complex_t(pGreen[(a0 + jsho)*nrhs + irhs][RealPart][spin][spjn*64 + xyz],  // load Re{<p|G>}
                                                              pGreen[(a0 + jsho)*nrhs + irhs][ImagPart][spin][spjn*64 + xyz]); // load Im{<p|G>}
                                    for (int ix = 0; ix <= lmax; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                        ax[ix] += pG * H1D[ix][0][x4];
                                    } // ix
                                } // inside mask
                            } // x4
                            for (int iy = 0; iy <= lmax; ++iy) { // loop over the 2nd Cartesian SHO quantum number
                                int const iyx0 = (iy*(2*lmax + 3 - iy)) >> 1;
                                auto const Hy = H1D[iy][1][y4]; // load Hy
                                for (int ix = 0; ix <= lmax - iy; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                    byx[iyx0 + ix] += ax[ix] * Hy;
                                } // ix
                            } // iy
                        } // inside mask
                    } // y4
                    for (int isho = 0, iz = 0; iz <= lmax; ++iz) { // loop over the 3rd Cartesian SHO quantum number
                        auto const Hz = H1D[iz][2][z4]; // load Hz
                        for (int iy = 0; iy <= lmax - iz; ++iy) { // loop over the 2nd Cartesian SHO quantum number
                            int const iyx0 = (iy*(2*lmax + 3 - iy)) >> 1;
                            for (int ix = 0; ix <= lmax - iz - iy; ++ix) { // loop over the 1st Cartesian SHO quantum number
                                czyx[isho] += byx[iyx0 + ix] * Hz; // form <p|G|p>
                                ++isho;
                            } // ix
                        } // iy
                    } // iz
                } // inside mask
            } // z4

        } // irhs

        for (int isho = 0; isho < nSHO; ++isho) {
            piGp[((isho*nSHO + jsho)*Noco + spin)*Noco + spjn] = czyx[isho];
        } // isho

        } // jsho
        } // spjn
        } // spin

            auto const iatom = AtomImageIndex[iai];
            assert(nSHO*nSHO*Noco*Noco == pGp[iatom].size() && "Inconsistency between AtomLmax[:] and AtomImageLmax[AtomImageIndex[:]]");
            auto const phase = AtomImagePhase[iai];
            auto const ph = complex_t(phase[0], phase[1]); // Block phase factor, ToDo: needs a conjugation here?
            for (int ij = 0; ij < nSHO*nSHO*Noco*Noco; ++ij) {
                complex_t const c = piGp[ij] * ph; // complex multiplication
                pGp[iatom][ij] += c.imag(); // add only the imaginary part
            } // ij

        } // iai

        return pGp;
    } // SHOprj

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_projection
