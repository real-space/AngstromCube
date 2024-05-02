#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::sprintf
#include <cmath> // std::exp, ::sqrt
#include <cassert> // assert
#include <algorithm> // std::max

#include "radial_grid.h" // radial_grid_t
#include "quantum_numbers.h" // ell_QN_t

#include "status.hxx" // status_t

namespace atom_core {

    status_t solve(
          double const Z // number of protons
        , int const echo=0 // log output level
        , char const config='a' // a:auto or custom
        , radial_grid_t const *rg=nullptr // optional: use this radial grid
        , double *export_Zeff=nullptr // optional result: -r*V_eff(r)
    ); // declaration only

    status_t scf_atom(
          radial_grid_t const & g // radial grid descriptor
        , double const Z // number of protons
        , int const echo=0 // log output level
        , double const occupations[][2]=nullptr // occupation numbers by nl_index and spin
        , double *export_Zeff=nullptr
    ); // declaration only

    status_t read_Zeff_from_file(
          double Zeff[] // -r*V_eff(r), -r*V_eff(r=0) should be ~= Z
        , radial_grid_t const & g // radial grid descriptor
        , double const Z // number of protons
        , char const *basename="pot/Zeff" // beginning of the filename
        , double const factor=1 // optional factor, e.g. -1 if the output is r*V(r)
        , int const echo=0 // log output level
        , char const *prefix="" // logging prefix
    ); // declaration only

    status_t store_Zeff_to_file(
          double const Zeff[] // -r*V_eff(r), -r*V_eff(r=0) should be ~= Z
        , double const r[] // radial grid support points
        , int const nr // number of radial grid points
        , double const Z // number of protons
        , char const *basename="pot/Zeff" // beginning of the filename
        , double const factor=1 // optional factor, e.g. -1 if the input is r*V(r)
        , int const echo=0 // log output level
        , char const *prefix="" // ="" // logging prefix
    ); // declaration only

    inline void get_Zeff_file_name(
          char *filename // result
        , char const *basename // filename before the dot
        , float const Z // number of protons
        , size_t const nchars=128
    ) {
        std::snprintf(filename, nchars, "%s.%03g", basename, Z);
    } // get_Zeff_file_name

    void rad_pot(
          double rV[] // result: r*V(r)
        , radial_grid_t const & g // radial grid descriptor
        , double const rho4pi[] // 4*pi*rho(r)
        , double const Z=0 // number of protons
        , double *energies=nullptr // result: energy contribution break down
    ); // declaration only

    inline double guess_energy(double const Z, int const enn) {
        auto const Zn2 = (Z*Z)/double(enn*enn);
        return -.5*Zn2 *  // Hydrogen-like energies in the Hartree unit system
              (.783517 + 2.5791E-5*Zn2) * // fit for the correct 1s energy
              std::exp(-.01*(enn - 1)*Z);
    } // guess_energy

    inline int nl_index(int const enn, int const ell) {
        assert(ell >= 0 && "angular momentum quantum number");
        assert(enn > ell && "atomic quantum numbers");
        return (enn*(enn - 1))/2 + ell;
    } // nl_index

    double neutral_atom_total_energy(double const Z); // declaration only

    status_t all_tests(int const echo=0); // declaration only

} // namespace atom_core
