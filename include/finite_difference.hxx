#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // uint32_t, int16_t
#include <vector> // std::vector<T>
#ifdef    GENERAL_CELL
    #include <cmath> // std::pow(std::complex<T>, std::complex<T>)
#endif // GENERAL_CELL

#include "real_space.hxx" // ::grid_t
#include "boundary_condition.hxx" // *_Boundary
#include "inline_math.hxx" // intpow, pow2
#include "status.hxx" // status_t
#include "uniform_laplacian.hxx" // ::get

namespace finite_difference {

    int constexpr nnArraySize = 16;

    template <typename real_t> // real_t may be float or double
    class stencil_t {
    public:
        real_t c2nd[3][nnArraySize]; // coefficients for the 2nd derivative
    private:
        int8_t _nn[3]; // number of FD neighbors

        void _constructor(double const grid_spacing[3], int const nneighbors[3], double const scale_factor); // declaration only

    public:

        stencil_t(double const grid_spacing[3], int const nneighbors[3], double const scale_factor=1) {
            _constructor(grid_spacing, nneighbors, scale_factor);
        } // preferred constructor

        stencil_t(double const grid_spacing[3], int const nn=4, double const scale_factor=1) {
            int const nns[3] = {nn, nn, nn};
            _constructor(grid_spacing, nns, scale_factor);
        } // isotropic nn constructor

        stencil_t(double const h=1, int const nn=4, double const scale_factor=1) {
            double const hgs[3] = {h, h, h};
            int const nns[3] = {nn, nn, nn};
            _constructor(hgs, nns, scale_factor);
        } // isotropic constructor, default constructor

        double clear_diagonal_elements() { // modifies the coefficients c2nd[][]
            double diag{0};
            for (int d = 0; d < 3; ++d) {
                diag += c2nd[d][0];
                c2nd[d][0] = 0; // clear diagonal elements
            } // d
            return diag;
        } // clear_diagonal_elements

        void scale_coefficients(double const f[3]) {
            for (int d = 0; d < 3; ++d) {
                for (int i = 0; i < nnArraySize; ++i) {
                    c2nd[d][i] *= f[d];
                } // i
            } // d
        } // scale_coefficients

        void scale_coefficients(double const f) { double const f3[] = {f, f, f}; scale_coefficients(f3); }
        int8_t const * nearest_neighbors() const { return _nn; }
        int nearest_neighbors(int const d) const { assert(d >= 0); assert(d < 3); return _nn[d]; }
        real_t const * Laplace_coefficients(int const d) const { assert(d >= 0); assert(d < 3); return c2nd[d]; }

    }; // class stencil_t

    template <typename real_t>
    double polar(std::complex<real_t> const c) { return std::atan2(c.imag(), c.real())*(180./constants::pi); }
    template <typename real_t> double polar(real_t const r) { return (r < 0)*180; }

    template <typename complex_out_t // result is stored in this precision
             ,typename complex_in_t // input comes in this precision
             ,typename real_fd_t> // computations are executed in this precision
    status_t apply(
          complex_out_t out[]
        , complex_in_t const in[]
        , real_space::grid_t const & g
        , stencil_t<real_fd_t> const & fd
        , double const factor=1
        , complex_in_t const boundary_phase[3][2]=nullptr
    ) {

        int constexpr n16 = nnArraySize; // max number of finite difference neighbors, typically 16
        typedef int16_t list_integ_t;
        std::vector<list_integ_t> list[3]; // can be of type int16_t
        std::vector<complex_in_t> phas[3];
        for (int d = 0; d < 3; ++d) {
            int const n = g[d];
            // check that n is smaller than the upper limit of int
            assert(0 < n); assert(n <= std::numeric_limits<list_integ_t>::max());
            int const bc = g.boundary_condition(d);
            int const nf = fd.nearest_neighbors(d);
            if (nf > n) error("finite-difference range (%d) in %c-direction is larger than grid (%d grid points)", nf, 'x' + d, n);
            assert(nf <= n);
            assert(nf <= n16);
            int const nh = n16 + n + n16; // number including largest halos
            list[d] = std::vector<list_integ_t>(nh, -1); // get memory, init as -1:non-existing
            phas[d] = std::vector<complex_in_t>(nh, 0); // get memory, init zero

            // core region
            for (int j = 0; j < n; ++j) {
                list[d][n16 + j] = j; // identity
                phas[d][n16 + j] = 1; // neutral
            } // j

            complex_in_t const phase_low = boundary_phase ? boundary_phase[d][0] : 1;
            complex_in_t const phase_upp = boundary_phase ? boundary_phase[d][1] : 1;

            // lower boundary
            if (Periodic_Boundary == bc || Shifted_Boundary == bc) { // periodic BC
                for (int j = -nf; j < 0; ++j) {
                    list[d][n16 + j] = (n + j) % n; // wrap around
                    phas[d][n16 + j] = phase_low; // incorrect if nf > n
                } // j
            } else if (Mirrored_Boundary == bc) { // mirror BC
                for (int j = -nf; j < 0; ++j) {
                    list[d][n16 + j] = - 1 - j; // mirror at -1 | 0
                    phas[d][n16 + j] = phase_low; // incorrect if nf > n
                } // j
            } // else open BC, list[:] = -1, phas[:] = 0

            // upper boundary
            if (Periodic_Boundary == bc || Shifted_Boundary == bc) { // periodic BC
                for (int j = 0; j < nf; ++j) {
                    list[d][n16 + n + j] = (n + j) % n; // wrap around
                    phas[d][n16 + n + j] = phase_upp; // incorrect if nf > n
                } // j
            } else if (Mirrored_Boundary == bc) { // mirror BC
                for (int j = 0; j < nf; ++j) {
                    list[d][n16 + n + j] = n - 1 - j; // mirror at n-1 | n
                    phas[d][n16 + n + j] = phase_upp; // incorrect if nf > n
                } // j
            } // else open BC, list[:] = -1, phas[:] = 0

            if (0) { // DEBUG: show indirection list and phase factors
                std::printf("# indirection list for %c-direction ", 'x'+d);
                for (int j = -nf; j < n + nf; ++j) {
                    if (0 == j || n == j) std::printf(" |");
                    std::printf(" %i", list[d][n16 + j]);
                } // j
                std::printf("\n");
                std::printf("# phase factor list for %c-direction ", 'x'+d);
                for (int j = -nf; j < n + nf; ++j) {
                    if (0 == j || n == j) std::printf(" |");
                    std::printf("  %g %g", std::real(phas[d][n16 + j]), std::imag(phas[d][n16 + j]));
                } // j
                std::printf("\n");
            } // show indirection list

        } // spatial direction d

        auto const nx = g('x'), ny = g('y'), nz = g('z');
#ifdef    GENERAL_CELL
        complex_in_t const phase_xy_low = std::pow(boundary_phase ? boundary_phase[0][0] : 1, complex_in_t(g.shift_yx/double(nx))),
                           phase_xy_upp = std::pow(boundary_phase ? boundary_phase[0][1] : 1, complex_in_t(g.shift_yx/double(nx)));
        if (0) { std::printf("# %s: low= %g %g = %g degrees, upp= %g %g = %g degrees\n", __func__,
                                std::real(phase_xy_low), std::imag(phase_xy_low), polar(phase_xy_low),
                                std::real(phase_xy_upp), std::imag(phase_xy_upp), polar(phase_xy_upp)); }
        assert(0 <= g.shift_yx); assert(g.shift_yx < nx);

        complex_in_t const phase_xz_low = std::pow(boundary_phase ? boundary_phase[0][0] : 1, complex_in_t(g.shift_zx/double(nx))),
                           phase_xz_upp = std::pow(boundary_phase ? boundary_phase[0][1] : 1, complex_in_t(g.shift_zx/double(nx)));
        complex_in_t const phase_yz_low = std::pow(boundary_phase ? boundary_phase[1][0] : 1, complex_in_t(g.shift_zy/double(ny))),
                           phase_yz_upp = std::pow(boundary_phase ? boundary_phase[1][1] : 1, complex_in_t(g.shift_zy/double(ny)));
#endif // GENERAL_CELL

        real_fd_t const scale_factor = factor;
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {

                    complex_out_t t(0); // init result

                    for (int d = 0; d < 3; ++d) {
                        int const nf = fd.nearest_neighbors(d);
                        int zyx[3] = {x, y, z};
                        int const i_center = zyx[d];
                        for (int jmi = -nf; jmi <= nf; ++jmi) {
                            int const j = i_center + jmi;
                            int const index = list[d][n16 + j]; // indirection to capture boundary
                            if (index >= 0) {
                                zyx[d] = index;
                                auto phase = phas[d][n16 + j]; // phase factor when crossing a periodic boundary
#ifdef    GENERAL_CELL
                                // ToDo: put this part into a section which is taken out in regular versions by a boolean template parameter

                                // allow shift-rectangular cells from lower triangular cell matrices
                                if (1 == d) { // derive in y-direction

                                    if (j >= ny) {
                                        zyx[0] = (x - g.shift_yx + nx) % nx; // modify the x-coordinate of the source
                                        phase *= phase_xy_upp;
                                        if (zyx[0] > x) { phase *= phas[0][n16 - 1]; }
                                    } else
                                    if (j < 0) {
                                        zyx[0] = (x + g.shift_yx + nx) % nx; // modify the x-coordinate of the source
                                        phase *= phase_xy_low;
                                        if (zyx[0] < x) { phase *= phas[0][n16 + nx]; } 
                                    } else {
                                        zyx[0] = x;
                                    }

                                } else // 'y'
                                if (2 == d) { // derive in z-direction

                                    if (j >= nz) {
                                        zyx[0] = (x - g.shift_zx + nx) % nx; // modify the x-coordinate of the source
                                        phase *= phase_xz_upp;
                                        if (zyx[0] > x) { phase *= phas[0][n16 - 1]; }
                                    } else
                                    if (j < 0) {
                                        zyx[0] = (x + g.shift_zx + nx) % nx; // modify the x-coordinate of the source
                                        phase *= phase_xz_low;
                                        if (zyx[0] < x) { phase *= phas[0][n16 + nx]; } 
                                    } else {
                                        zyx[0] = x;
                                    }

                                    if (j >= nz) {
                                        zyx[1] = (y - g.shift_zy + ny) % ny; // modify the y-coordinate of the source
                                        phase *= phase_yz_upp;
                                        if (zyx[1] > y) { phase *= phas[1][n16 - 1]; }
                                    } else
                                    if (j < 0) {
                                        zyx[1] = (y + g.shift_zy + ny) % ny; // modify the y-coordinate of the source
                                        phase *= phase_yz_low;
                                        if (zyx[1] < y) { phase *= phas[1][n16 + ny]; } 
                                    } else {
                                        zyx[1] = y;
                                    }

                                } // 'z'
                                assert(0 <= zyx[0]); assert(zyx[0] < nx);
                                assert(0 <= zyx[1]); assert(zyx[1] < ny);
                                assert(0 <= zyx[2]); assert(zyx[2] < nz);
#endif // GENERAL_CELL
                                int const jzyx = (zyx[2]*ny + zyx[1])*nx + zyx[0]; // source index
                                auto const coeff = fd.c2nd[d][std::abs(jmi)];
                                auto const contrib = phase * in[jzyx];
                                //   if (29 == x && 0 == y) { std::printf("# FD add %6.1f + %6.1f [%9.3f] = %6.1f \t from [%i %i %i]\n",
                                //                 polar(phase), polar(in[jzyx]), coeff, polar(contrib), zyx[2], zyx[1], zyx[0]); }
                                t += contrib*coeff;
                            } // index exists
                        } // jmi
                    } // d direction of the derivative

                    int const izyx = (z*ny + y)*nx + x;
                    out[izyx] = t * scale_factor; // store

                } // x
            } // y
        } // z

        return 0; // success
    } // apply


    status_t all_tests(int const echo=0); // declaration only

} // namespace finite_difference
