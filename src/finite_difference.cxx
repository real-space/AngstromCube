// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <vector> // std::vector<T>

#include "finite_difference.hxx"

#include "uniform_laplacian.hxx" // ::get
#include "status.hxx" // status_t
#include "real_space.hxx" // ::grid_t
#include "boundary_condition.hxx" // *_Boundary
#include "recorded_warnings.hxx" // warn
#include "inline_math.hxx" // pow2
#include "control.hxx" // ::get

namespace finite_difference {

    template <typename real_t> // real_t may be float or double
    void stencil_t<real_t>::_constructor(double const grid_spacing[3], int const nneighbors[3], double const scale_factor) {
        for (int d = 0; d < 3; ++d) {
            for (int i = 0; i < nnArraySize; ++i) { c2nd[d][i] = 0; } // clear
            _nn[d] = uniform_laplacian::get(c2nd[d], nneighbors[d], grid_spacing[d], 'x' + d);
            if (_nn[d] < nneighbors[d]) {
                warn("In stencil_t requested nn=%i but use nn=%i for %c-direction", nneighbors[d], _nn[d], 'x' + d);
            }
        } // d spatial direction
        scale_coefficients(scale_factor);
    } // _constructor






#ifdef    NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    status_t test_create_and_destroy(int const echo=9) {
        auto const f = new stencil_t<float>();
        f->~stencil_t();
        stencil_t<double> d;
        return 0;
    } // test_create_and_destroy

    template <typename real_t>
    status_t test_Laplacian(int const echo=3) {
        status_t stat(0);
        double const h[3] = {1, 1, 1}; // unit grid spacings
        for (int dir = 0; dir < 3; ++dir) {
            int nn[3] = {0,0,0}; nn[dir] = 12; // switch FD off for the two perpendicular directions
            stencil_t<real_t> Laplacian(h, nn);
            int dims[] = {1,1,1}; dims[dir] = 127 + dir;
            real_space::grid_t g(dims);
            g.set_boundary_conditions(Periodic_Boundary);
            double const k = (1 + dir)*2*constants::pi/g[dir]; // wave vector of a single plane wave
            std::vector<real_t> values(g.all()), result(g.all());
            for (size_t i = 0; i < g.all(); ++i) values[i] = std::cos(k*i); // fill with some non-zero values
            stat += finite_difference::apply(result.data(), values.data(), g, Laplacian);
            if (echo > 5) std::printf("\n# in, result, ref values:\n");
            double dev{0};
            for (size_t i = 0; i < g.all(); ++i) {
                auto const ref = -k*k*values[i]; // analytic solution to the Laplacian operator applied to a plane wave
                if (echo > 5) std::printf("%ld %g %g %g\n", i, values[i], result[i], ref);
                // compare in the middle range result and ref values
                dev += std::abs(result[i] - ref);
            } // i
            if (echo > 2) std::printf("# %s %c-direction: dev = %g\n", __func__, 'x'+dir, dev);
        } // direction
        return stat;
    } // test_Laplacian

    template <typename real_t>
    status_t test_Bloch_wave(int const echo=3) {
        status_t stat(0);
        double const h[3] = {1, 1, 1}; // unit grid spacings
        std::complex<real_t> boundary_phase[3][2] = {{-1,-1}, {-1,-1}, {-1,-1}};
        double maxdev{0};
        for (int dir = 0; dir < 3; ++dir) {
            int nn[3] = {0,0,0}; nn[dir] = 12; // switch FD off for the two perpendicular directions
            stencil_t<real_t> Laplacian(h, nn);
            int dims[] = {1,1,1}; dims[dir] = 127 + dir;
            real_space::grid_t g(dims);
            g.set_boundary_conditions(Periodic_Boundary);
            std::vector<std::complex<real_t>> values(g.all()), result(g.all());
            for (int iphase = 0; iphase <= 180; iphase += 20) {
                double const k = (1 + dir + iphase/360.)*2*constants::pi/g[dir]; // wave vector of a single plane wave
                double const arg = iphase*constants::pi/180.;
                boundary_phase[dir][0] = std::complex<real_t>(std::cos(arg), -std::sin(arg));
                boundary_phase[dir][1] = real_t(1)/boundary_phase[dir][0];
                for (size_t i = 0; i < g.all(); ++i) {
                    values[i] = std::complex<real_t>(std::cos(k*i), std::sin(k*i));
                } // i
                stat += finite_difference::apply(result.data(), values.data(), g, Laplacian, 1, boundary_phase);
                double dev{0};
                for (size_t i = 0; i < g.all(); ++i) {
                    auto const ref = -k*k*values[i]; // analytic solution to the Laplacian operator applied to a plane wave
                    // compare in the middle range result and ref values
                    dev += std::abs(result[i] - ref);
                } // i
                if (echo > 3) std::printf("# %s %c-direction: dev = %g\n", __func__, 'x'+dir, dev);
                maxdev = std::max(maxdev, std::abs(dev));
            } // iphase
        } // direction
        if (echo > 1) std::printf("\n# %s largest deviation is %.1e\n", __func__, maxdev);
        return stat;
    } // test_Bloch_wave

#ifdef    GENERAL_CELL
    template <typename real_t>
    status_t test_general_cell(int const echo=3) {
        // ToDo: not covered: more than one shift non-zero at the same time
        status_t stat(0);
        // create a shifted Cartesian unit cell with angles 60 degree
        int const dims[] = {15, 13, 12};
        double const h[3] = {1./dims[0], std::sqrt(3/4.)/dims[1], std::sqrt(2/3.)/dims[2]}; // grid spacings
        int const nn[3] = {12, 12, 12}; // FD-order
        stencil_t<real_t> Laplacian(h, nn);
        real_space::grid_t g(dims);
        std::vector<std::complex<real_t>> values(g.all()), result(g.all());
        std::complex<real_t> boundary_phase[3][2];
        g.set_boundary_conditions(Periodic_Boundary, Shifted_Boundary, Shifted_Boundary);
        if (echo > 1) std::printf("\n# %s: start\n", __func__);
        char const shift_name[8][8] = {"noshift", "xy", "xz", "xz+yz", "yz", "xy+yz", "xy+xz", "xyxzyz"};
                               // is == 000        001   010    011     100    101      110      111
        int const ntests = control::get("finite_difference.general_cell.test", 4.); // 0:nothing, 1:unshifted, 4:minimum, 7:mixed directions, 8:all
        int const incdeg = control::get("finite_difference.general_cell.increment", 15.); // 1 degree:very fine grained, 90 degrees: only Cartesian directions
        int8_t const is_of_itest[] = {0, 1,2,4, 3,5,6, 7};
        double maxdevall{0};
        for (int itest{0}; itest < std::min(8, ntests); ++itest) {
            int const is = is_of_itest[itest];
            int const shift_dims[] = {((is >> 0) & 0x1)*(dims[0] - 1),      // max shift along x-direction on crossing a y-boundary
                                      ((is >> 1) & 0x1)*(dims[0] - 1),      // max shift along x-direction on crossing a z-boundary
                                      ((is >> 2) & 0x1)*(dims[1] - 1)};     // max shift along y-direction on crossing a z-boundary
            if (echo > 3) { std::printf("# %s: test xy_shift in [0, %d], xz_shift in [0, %d], yz_shift in [0, %d], name= %s\n",
                                __func__, shift_dims[0], shift_dims[1], shift_dims[2], shift_name[is]); }
            double maxdev_is{0};
            for (int idirection{0}; idirection <= 90; idirection += incdeg) { // angle
                // prepare plane wave vector (this could be moved outside the shift loops to save time)
                double const k = 1.6; // sqRy
                auto constexpr arc = constants::pi/180.;
                auto const k_cos = k*std::cos(idirection*arc), k_sin = k*std::sin(idirection*arc);
                double const kv_dir[] = {k_cos, k_sin, .1, k_cos, k_sin, .1, k_cos, k_sin, .1, k_cos};
                double const *const kv = kv_dir + is; // only reference *kv as kv[3]
                // prepare wave values
                for (int iz{0}; iz < g('z'); ++iz) {
                    for (int iy{0}; iy < g('y'); ++iy) {
                        for (int ix{0}; ix < g('x'); ++ix) {
                            auto const arg = kv[0]*ix*h[0] + kv[1]*iy*h[1] + kv[2]*iz*h[2]; // prepare an arbitrary plane wave
                            values[(iz*g('y') + iy)*g('x') + ix] = std::complex<real_t>(std::cos(arg), std::sin(arg));
                        } // ix
                    } // iy
                } // iz
                // prepare matching boundary phases
                for (int dir{0}; dir < 3; ++dir) {
                    auto const arg = kv[dir]*g[dir]*h[dir];
                    boundary_phase[dir][1] = std::complex<real_t>(std::cos(arg), std::sin(arg));
                    boundary_phase[dir][0] = real_t(1)/boundary_phase[dir][1];
                    // if (echo > 11) { std::printf("# %s: %d deg, %c-phase= %g %g\n", __func__, idirection, 'x'+dir, boundary_phase[dir][1].real(), boundary_phase[dir][1].imag()); }
                } // dir

                double maxdev{0};
                for (int xy_shift{0}; xy_shift <= shift_dims[0]; ++xy_shift) {          // shift along x-direction on crossing a y-boundary
                for (int xz_shift{0}; xz_shift <= shift_dims[1]; ++xz_shift) {          // shift along x-direction on crossing a z-boundary
                for (int yz_shift{0}; yz_shift <= shift_dims[2]; ++yz_shift) {          // shift along y-direction on crossing a z-boundary
                    if (echo > 13) { std::printf("# %s: test %s-shifts\n", __func__, shift_name[is]); }
                    double const cell_shape[3][4] = {{h[0]*dims[0],             0,             0, 0},
                                                     {h[0]*xy_shift, h[1]*dims[1],             0, 0},  // lower triangular matrix
                                                     {h[0]*xz_shift, h[1]*yz_shift, h[2]*dims[2], 0}};
                    g.set_cell_shape(cell_shape, echo/4);

                    // apply
                    stat += finite_difference::apply(result.data(), values.data(), g, Laplacian, 1, boundary_phase);

                    // compare to anaytical Laplacian
                    auto const k2 = pow2(kv[0]) + pow2(kv[1]) + pow2(kv[2]);
                    double dev{0};
                    for (size_t i{0}; i < g.all(); ++i) {
                        auto const val = values[i], res = result[i], ref = -k2*val; // reference is the analytic solution to the Laplacian operator applied to a plane wave
                        dev += std::abs(res - ref);
                    } // i
                    if (echo > 9) { std::printf("# %s: direction=%4d degrees xy= %d/%d, xz= %d/%d, yz= %d/%d, dev= %.2e\n", __func__,
                                                idirection, xy_shift, dims[0], xz_shift, dims[0], yz_shift, dims[1], dev/g.all()); }
                    maxdev = std::max(maxdev, std::abs(dev/g.all()));
                }}} // *_shift

                stat += (maxdev > 1e-12);
                if (echo > 4) std::printf("# %s: direction=%4d degrees, dev= %g\n", __func__, idirection, maxdev);
                maxdev_is = std::max(maxdev_is, maxdev);
            } // idirection

            if (echo > 0) { std::printf("# %s: largest deviation tests \'%s\' is %.1e\n", __func__, shift_name[is], maxdev_is); }
            maxdevall = std::max(maxdevall, maxdev_is);
        } // itest --> is
        if (echo > 0) { std::printf("# %s: largest deviation of all tests is %.1e\n", __func__, maxdevall); }
        return stat;
    } // test_general_cell
#endif // GENERAL_CELL

    status_t all_tests(int const echo) {
        status_t stat(0);
#ifdef    GENERAL_CELL
        stat += test_general_cell<double>(echo);
#endif // GENERAL_CELL
        stat += test_create_and_destroy(echo);
        stat += test_Laplacian<double>(echo);
        stat += test_Bloch_wave<double>(echo);
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace finite_difference
