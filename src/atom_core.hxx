#pragma once

#include <cstdio> // std::sprintf
#include <cmath> // std::exp, std::sqrt
#include <cassert> // assert
#include <algorithm> // std::max

#include "radial_grid.h" // radial_grid_t
#include "quantum_numbers.h" // ell_QN_t

#include "status.hxx" // status_t

namespace atom_core {

  status_t scf_atom(
      radial_grid_t const &g, // radial grid descriptor
      double const Z, // atomic number
      int const echo=0, // log output level
      double const occupations[][2]=nullptr); // occupation numbers by nl_index

  status_t read_Zeff_from_file(double Zeff[], radial_grid_t const &g, double const Z,
                    char const basename[]="pot/Zeff", double const factor=1, int const echo=0, char const *prefix="");
  
  status_t store_Zeff_to_file(double const Zeff[], double const r[], int const nr, double const Z,
                    char const basename[]="pot/Zeff", double const factor=1, int const echo=9);

  inline void get_Zeff_file_name(char *filename, char const *basename, float const Z) {
      std::sprintf(filename, "%s.%03g", basename, Z); }

  double initial_density(double r2rho[], radial_grid_t const &g, double const Z, double const charged=0);

  void rad_pot(double rV[], radial_grid_t const &g, double const rho4pi[], double const Z=0, double *energies=nullptr);

  inline double guess_energy(double const Z, int const enn) {
      return -.5*(Z/enn)*(Z/enn) *  // Hydrogen-like energies in the Hartree unit system
            (.783517 + 2.5791E-5*(Z/enn)*(Z/enn)) * // fit for the correct 1s energy
            std::exp(-.01*(enn - 1)*Z); // guess energy
  } // guess_energy

  inline int nl_index(int const enn, int const ell) { 
      assert(ell >= 0); assert(enn > ell); // atomic eigenstates
      return (enn*(enn - 1))/2 + ell;
  } // nl_index

  inline char ellchar(ell_QN_t const ell) {
      char constexpr special_ellchars[] = "spd";
      if (ell < 0) return '?';
      if (ell < 3) return special_ellchars[ell];
      return 99 + ell; // "fghijk ..." // ToDo: actually 'i' does not belong in here
  } // ellchar
  
  inline double neutral_atom_total_energy_LDA(double const Z) {
      // fitting LDA total energies/Z^2 for Z=10..120
      double const a0 = 0.18094;
      double const a1 = 0.383205;
      double const a2 = -0.0109251;
      double const a3 = 4.75216e-05;
// or fitting LDA total energies/Z^2.5 for Z=10..120
// 	a0 = 0.178902
// 	a1 = 0.384094
// 	a2 = -0.0110257
// 	a3 = 4.78395e-05
      return -std::max(0., Z)*Z*(a0 + a1*std::sqrt(std::max(0., Z)) + a2*Z + a3*Z*Z);
  } // neutral_atom_total_energy_LDA

  inline double neutral_atom_total_energy(double const Z) {
      if (Z <= 0) return 0;
      double const E_LDA[16] = {0,-0.445893560,-2.834410555,-7.333797749,
        -14.448933773,-24.350492007,-37.440386817,-54.053760337,
        -74.524727413,-99.186432180, -128.371547297,-161.650489998,
        -199.451741467,-241.760703146,-288.815708267,-340.781209719};
      if (Z <= 0) return 0;
      if (Z > 16) return neutral_atom_total_energy_LDA(Z);
      int const iZ = int(Z);
      assert(iZ >= 0 && iZ <= 16);
      // compute weights for a cubic Lagrange polynomial through E[iZ-1], E[iZ] and E[iZ+1]
      double const xm1 = iZ - 1, x_0 = iZ, xp1 = iZ + 1;
      double const wm1 = ((Z - x_0)*(Z - xp1))/((xm1 - x_0)*(xm1 - xp1));
      double const w_0 = ((Z - xp1)*(Z - xm1))/((x_0 - xp1)*(x_0 - xm1));
      double const wp1 = ((Z - xm1)*(Z - x_0))/((xp1 - xm1)*(xp1 - x_0));
      double const ym1 = (iZ >  0)? E_LDA[iZ - 1] : 0;
      double const y_0 = (iZ < 16)? E_LDA[iZ]     : neutral_atom_total_energy_LDA(x_0);
      double const yp1 = (iZ < 15)? E_LDA[iZ + 1] : neutral_atom_total_energy_LDA(xp1);
      return wm1*ym1 + w_0*y_0 + wp1*yp1;
  } // neutral_atom_total_energy

  status_t all_tests(int const echo=0); // declaration only

} // namespace atom_core
