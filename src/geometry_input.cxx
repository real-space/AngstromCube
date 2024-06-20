// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, ::snprintf
#include <fstream> // std::ifstream
#include <string> // std::string, ::getline
#include <sstream> // std::istringstream
#include <cstdint> // int32_t

#include "geometry_input.hxx"

#include "real_space.hxx" // ::grid_t
#include "chemical_symbol.hxx" // ::decode
#include "control.hxx" // ::get, ::set, ::echo_set_without_warning
#include "mpi_parallel.hxx" // ::comm, ::rank, ::broadcast
#include "data_view.hxx" // view2D<T>
#include "unit_system.hxx" // ::length_unit, ::energy_unit

namespace geometry_input {

    double constexpr Bohr2Angstrom = 0.52917724924;
    double constexpr Angstrom2Bohr = 1./Bohr2Angstrom;

    status_t read_xyz_file(
          view2D<double> & xyzZ
        , int32_t & n_atoms
        , double cell[3][4]
        , int8_t bc[3] // =nullptr
        , char const *const filename // ="atoms.xyz"
        , int const echo // =5 log-level
    ) {
        auto const comm = mpi_parallel::comm(); // default_communicator
        auto const me = mpi_parallel::rank(comm);
        int return_status{0};
        if (0 == me) { // MPI master task

            std::ifstream infile(filename, std::ifstream::in);
            if (infile.fail()) { error("unable to open file '%s' for reading coordinates", filename); }

            int64_t natoms{0}, linenumber{2};
            infile >> natoms; // read the number of atoms
            if (natoms < 0) error("%s:1 indicates a negative number of atoms, %ld", filename, natoms);
            std::string line;
            std::getline(infile, line); // read line #1 again
            std::getline(infile, line); // read line #2
            if (echo > 2) std::printf("# %s: expect %lld atoms, cell parameters= \"%s\"\n", filename, natoms, line.c_str());
            { // scope: parse line #2
                std::istringstream iss(line);
                std::string word;
                if ('%' == line[0]) {
                    double L[3][4] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}};
                    if (echo > 0) std::printf("# 1st char in 2nd line is %%, read periodic unit cell in file \'%s\': %s\n", filename, line.c_str());
                    iss >> word >> L[0][0] >> L[0][1] >> L[0][2] >> L[1][0] >> L[1][1] >> L[1][2] >> L[2][0] >> L[2][1] >> L[2][2];
                    if (nullptr != cell) set(cell[0], 12, L[0], Angstrom2Bohr);
                    if (nullptr != bc) set(bc, 3, Periodic_Boundary);
                } else {
                    double L[3] = {0,0,0};
                    std::string B[3];
                    iss >> word >> L[0] >> L[1] >> L[2] >> B[0] >> B[1] >> B[2]; // Cartesian mode
                    if (nullptr != cell) set(cell[0], 12, 0.0); // clear
                    for (int d{0}; d < 3; ++d) {
                        if (nullptr != cell) {
                            cell[d][d] = L[d] * Angstrom2Bohr;
                            assert(cell[d][d] > 0);
                        }
                        if (nullptr != bc) bc[d] = boundary_condition::fromString(B[d].c_str(), echo, 'x' + d);
                    } // d
                } // cell == Basis
            } // scope

            xyzZ = view2D<double>(natoms, 4, 0.0);
            int64_t ia{0}; // global index of atoms
            size_t ncommented_lines{0}, ncommented_atoms{0}, nignored_lines{0}, nempty_lines{0}, nignored_atoms{0}, ninvalid_atoms{0};
            while (std::getline(infile, line)) {
                ++linenumber;
                std::istringstream iss(line);
                std::string Symbol;
                double pos[3]; // position vector
                if (iss >> Symbol >> pos[0] >> pos[1] >> pos[2]) {
                    char const *const Sy = Symbol.c_str();
                    char const S = Sy[0];
                    assert('\0' != S); // Symbol should never be an empty string
                    if ('#' != S) {
                        // add an atom
                        char const y = (Sy[1]) ? Sy[1] : ' ';
                        if (echo > 7) std::printf("# %c%c  %16.9f %16.9f %16.9f\n", S,y, pos[0], pos[1], pos[2]);
                        int const iZ = chemical_symbol::decode(S, y);
                        if (iZ >= 0) {
                            if (ia < natoms) {
                                // add atom to list
                                set(xyzZ[ia], 3, pos, Angstrom2Bohr);
                                xyzZ(ia,3) = iZ;
                            } else {
                                ++nignored_atoms;
                                if (echo > 3) std::printf("# %s:%lli contains \"%s\" as atom #%lli but only %lld atoms expected!\n",
                                                            filename, linenumber, line.c_str(), ia, natoms);
                            }
                            ++ia;
                        } else {
                            ++ninvalid_atoms;
                            if (echo > 3) std::printf("# %s:%lli contains invalid atom entry \"%s\"!\n",
                                                        filename, linenumber, line.c_str());
                        }
                    } else {
                        ++ncommented_atoms; // this is a comment line that could be interpreted as atom if it did not start from '#'
                        if (echo > 4) std::printf("# %s:%lli \t commented atom=\"%s\"\n", filename, linenumber, line.c_str());
                    }
                } else { // iss
                    if ('#' == Symbol[0]) {
                        ++ncommented_lines; // all exected atoms have been found, maybe this is a comment
                        if (echo > 1) std::printf("# %s:%lli \t comment=\"%s\"\n", filename, linenumber, line.c_str());
                    } else if ('\0' == Symbol[0]) {
                        ++nempty_lines; // ignore empty lines silently
                    } else {
                        ++nignored_lines;
                        if (echo > 5) std::printf("# %s:%lli \t ignored line=\"%s\"\n", filename, linenumber, line.c_str());
                    }
                } // iss
            } // parse file line by line

            // show irregularites
            if (echo > 0 && ncommented_lines > 0) std::printf("# %s: found %ld commented lines\n", filename, ncommented_lines);
            if (echo > 0 && ncommented_atoms > 0) std::printf("# %s: found %ld commented atoms\n", filename, ncommented_atoms);
            if (echo > 0 && ninvalid_atoms   > 0) std::printf("# %s: found %ld invalid atoms\n",   filename, ninvalid_atoms);
            if (echo > 0 && nignored_lines   > 0) std::printf("# %s: ignored %ld lines\n",         filename, nignored_lines);
            if (echo > 0 && nignored_atoms   > 0) std::printf("# %s: ignored %ld valid atoms\n",   filename, nignored_atoms);
            if (echo > 5 && nempty_lines     > 0) std::printf("# %s: found %ld blank lines\n",     filename, nempty_lines);

            n_atoms = std::min(ia, natoms); // export the number of atoms
            return_status = ia - natoms;

        } // end MPI master task

        mpi_parallel::broadcast(&n_atoms, comm);
        mpi_parallel::broadcast(bc, comm, 3);
        mpi_parallel::broadcast(cell[0], comm, 3*4);
        mpi_parallel::broadcast(&return_status, comm);

        if (0 != me) {
            if (echo > 0) std::printf("# rank#%i received n_atoms= %d, bc= %d %d %d, status= %i from master\n",
                                                me,         n_atoms,     bc[0], bc[1], bc[2],  return_status);
            xyzZ = view2D<double>(n_atoms, 4, 0.0);
        } // not MPI master

        // send the atomic positions and atomic numbers from master to all others
        mpi_parallel::broadcast(xyzZ.data(), comm, n_atoms*4);

        return return_status; // returns 0 if the expected number of atoms has been found
    } // read_xyz_file



    inline int even(int const any, unsigned const n_even=2) { return ((any - 1 + n_even)/n_even)*n_even; }
    inline int n_grid_points(double const suggest, unsigned const n_even=2) { return even(int(std::ceil(suggest)), n_even); }

    status_t init_geometry_and_grid(
            real_space::grid_t & g // output grid descriptor
          , view2D<double> & xyzZ // output atom coordinates and core charges Z
          , int32_t & natoms // output number of atoms found
          , unsigned const n_even // =2
          , int const echo // =0 log-level
    ) {
        status_t stat(0);
        assert(n_even > 0);

        natoms = 0; // number of atoms
        int8_t bc[3]; // boundary conditions
        double cell[3][4] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}}; // general cell parameters
        auto const geo_file = control::get("geometry.file", "atoms.xyz");
        if (echo > 3) std::printf("# +geometry.file=%s\n", geo_file);
        stat += read_xyz_file(xyzZ, natoms, cell, bc, geo_file, echo);

        auto const grid_spacing_unit_name = control::get("grid.spacing.unit", "Bohr");
        char const *_lu;
        auto const lu = unit_system::length_unit(grid_spacing_unit_name, &_lu);
        assert(lu > 0);
        auto const in_lu = 1./lu;
        int32_t ng_input[] = {0, 0, 0};
        double  default_hg_used{0};
        int32_t isotropic_ng_used{0};
        auto const keyword_ng = "grid.points";
        auto const keyword_hg = "grid.spacing";

        auto const ng_isotropic = control::get(keyword_ng, 0.); // 0 is not a usable default value, --> try to use grid spacings

        { // scope: determine grid spacings and number of grid points

            // precedence:
            //    highest:  grid.points.x, .y, .z
            //           :  grid.points
            //           :  grid.spacing.x, .y, .z
            //     lowest:  grid.spacing             default value = 0.125 Angstrom

            auto const default_grid_spacing = 0.23621577; // == 0.125 Angstrom
            auto const hg_isotropic = std::abs(control::get(keyword_hg, -default_grid_spacing*lu));

            int default_hg_times{0}, default_hg_xyz{0};
            int isotropic_ng_times{0}, isotropic_ng_xyz{0};
            int32_t ng[3] = {0, 0, 0};
            for (int d{0}; d < 3; ++d) { // directions x, y, z
                char keyword[32];
                std::snprintf(keyword, 32, "%s.%c", keyword_ng, 'x'+d);
                ng[d] = int32_t(control::get(keyword, ng_isotropic)); // "grid.points.x", ".y", ".z"
                ng_input[d] = ng[d];
                if (ng[d] < 1) {
                    // try to initialize the grid via "grid.spacing"
                    std::snprintf(keyword, 32, "%s.%c", keyword_hg, 'x'+d);
                    double const hg_lu = control::get(keyword, hg_isotropic); // "grid.spacing.x", ".y", ".z"
                    bool const is_default_hg = (hg_lu == hg_isotropic);
                    double const hg = hg_lu*in_lu;
                    if (echo > 8) std::printf("# input grid spacing in %c-direction is %g %s = %g %s%s\n",
                                'x'+d, hg_lu, _lu, hg*Ang, _Ang, is_default_hg?" (default)":"");
                    default_hg_times +=  int(is_default_hg);
                    default_hg_xyz   += (int(is_default_hg) << d);
                    if (hg <= 0) error("grid spacings must be positive, found %g %s in %c-direction", hg*Ang, _Ang, 'x'+d);
                    ng[d] = n_grid_points(std::abs(cell[d][d])/hg, n_even);
                    if (ng[d] < 1) error("no grid points with grid spacings %g %s in %c-direction", hg*Ang, _Ang, 'x'+d);
                    ng[d] = even(ng[d], n_even); // if odd, increment to nearest higher even number
                    if (is_default_hg) default_hg_used = std::abs(cell[d][d]/ng[d]);
                } else {
                    ng[d] = even(ng[d], n_even); // if odd, increment to nearest higher even number
                }
                if (ng_isotropic == ng[d]) {
                    ++isotropic_ng_times;
                    isotropic_ng_xyz += (1 << d);
                    isotropic_ng_used = ng_isotropic;
                } // isotropic_ng
                if (echo > 8) std::printf("# use %d grid points in %c-direction\n", ng[d], 'x'+d);
            } // d

            char const which_ones[8][16] = {"", ", x", ", y", ", x and y", ", z", ", x and z", ", y and z", ", x, y, and z"};
            if (default_hg_times > 0) {
                if (echo > 6) std::printf("# default grid.spacing %g %s used for %d directions%s\n",
                    default_hg_used*Ang, _Ang, default_hg_times, which_ones[default_hg_xyz & 7]);
            } // default_hg_times
            if (isotropic_ng_times > 0) {
                if (echo > 6) std::printf("# isotropic grid.points=%d used for %d directions%s\n",
                    isotropic_ng_used, isotropic_ng_times, which_ones[isotropic_ng_xyz & 7]);
            } // isotropic_ng_times

            g = real_space::grid_t(ng[0], ng[1], ng[2]);
        } // scope

        if (echo > 1) std::printf("# use  %d x %d x %d  grid points\n", g[0], g[1], g[2]);
        g.set_boundary_conditions(bc[0], bc[1], bc[2]);
        g.set_cell_shape(cell, echo);
        assert(g.is_shifted() && "self_consistency module can only treat lower triangular cell matrices!");

        double const max_grid_spacing = std::max(std::max(std::max(1e-9, g.h[0]), g.h[1]), g.h[2]);
        if (echo > 1) std::printf("# use  %g %g %g  %s  dense grid spacing, corresponds to %.1f Ry\n",
              g.h[0]*Ang, g.h[1]*Ang, g.h[2]*Ang, _Ang, pow2(constants::pi/max_grid_spacing));
        for (int d{0}; d < 3; ++d) {
            assert(std::abs(g.h[d]*g.inv_h[d] - 1) < 4e-16 && "grid spacing and its inverse do not match");
        } // d
        if (g.is_Cartesian()) {
            if (echo > 1) std::printf("# cell is  %.9f %.9f %.9f  %s\n", cell[0][0]*Ang, cell[1][1]*Ang, cell[2][2]*Ang, _Ang);
            for (int d{0}; d < 3; ++d) {
                if (std::abs(g.h[d]*g[d] - cell[d][d]) >= 1e-6) {
                    warn("grid in %c-direction seems inconsistent, %d * %g differs from %g %s",
                            'x'+d, g[d], g.h[d]*Ang, cell[d][d]*Ang, _Ang);
                    ++stat;
                } // deviates
            } // d
        } else if (echo > 1) {
            for (int d{0}; d < 3; ++d) {
                std::printf("# cell is  %15.9f %15.9f %15.9f  %s\n", g.cell[d][0]*Ang, g.cell[d][1]*Ang, g.cell[d][2]*Ang, _Ang);
            } // d
        } // is_Cartesian

        { // scope: correct the global variables, suppress warnings
            for (int d{0}; d < 3; ++d) { // directions x, y, z
                char keyword[32];
                if (ng_input[d] != g[d]) { // correct only if not matching (correcting removes the data origin visible with +control.show=-6)
                    std::snprintf(keyword, 32, "%s.%c", keyword_ng, 'x'+d);
                    control::set(keyword, g[d], control::echo_set_without_warning);
                } // correct only if not matching
                // always correct the grid spacing to the exact floating point value used
                std::snprintf(keyword, 32, "%s.%c", keyword_hg, 'x'+d);
                control::set(keyword, g.h[d]*lu, control::echo_set_without_warning);
            } // d
            double const show_isotropic_ng = (isotropic_ng_used > 0) ? isotropic_ng_used : std::cbrt(g[2]*1.*g[1]*1.*g[0]);
            if (ng_isotropic != show_isotropic_ng) {
                control::set(keyword_ng, show_isotropic_ng, control::echo_set_without_warning);
            } // correct only if not matching
            double const show_default_hg = (default_hg_used > 0) ? default_hg_used : std::cbrt(g.dV());
            control::set(keyword_hg, show_default_hg*lu,              control::echo_set_without_warning);
        } // scope

        return stat;
    } // init_geometry_and_grid


    double get_temperature(int const echo, double const def) { // def=1e-3
        auto const unit = control::get("electronic.temperature.unit", "Ha");
        char const *_eu;
        auto const eu = unit_system::energy_unit(unit, &_eu);
        auto const temp = control::get("electronic.temperature", def*eu)/eu;
        if (echo > 0) std::printf("# electronic.temperature= %g %s == %g %s\n", temp*eu, _eu, temp*eV, _eV);
        return temp;
    } // get_temperature

#ifdef    NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    status_t test_init(int const echo=3) {
        real_space::grid_t g;
        view2D<double> xyzZ;
        int32_t natoms;
        return init_geometry_and_grid(g, xyzZ, natoms, 8, echo);
    } // test_init

    status_t all_tests(int const echo) {
        status_t stat(0);
        stat += test_init(echo);
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace geometry_input
