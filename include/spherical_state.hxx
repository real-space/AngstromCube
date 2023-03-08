#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cmath> // std::sqrt
#include <algorithm> // std::max

#include "energy_level.hxx" // energy_level_t, TRU_ONLY
#include "display_units.h" // eV, _eV, Ang, _Ang
#include "radial_grid.hxx" // radial_grid_t

  typedef struct energy_level_t<TRU_ONLY> spherical_state_t; // Pseudo=TRU_ONLY=1: spherical states, e.g. core states, only a TRU wave is stored

  int constexpr core=0, semicore=1, valence=2, csv_undefined=3; // ToDo: as enum to distinguish the different energy level classes

  inline char const * csv_name(int const csv) { // 0:core, 1:semicore, 2:valence, else:?
      return (core     == csv) ? "core" : (
             (valence  == csv) ? "valence" : (
             (semicore == csv) ? "semicore" : "?" ) );
  } // csv_name

  inline void show_state(
        char const *label // log-prefix
      , char const *csv_class // classification as {core, semicore, valence, ?}
      , char const *tag // name of the state
      , double const occ // occupation number
      , double const energy // energy eigenvalue or energy parameter
      , char const final='\n' // line end
  ) {
      std::printf("# %s %-9s%-4s%6.1f E=%16.6f %s%c", label, csv_class, tag, occ, energy*eV, _eV, final);
  } // show_state

  inline double show_state_analysis( // returns the charge outside the sphere
        int const echo // log-level
      , char const *label // log-prefix
      , radial_grid_t const & rg // radial grid descriptor, should be rg[TRU]
      , double const wave[] // true radial wave function (Mind: not scaled by r)
      , char const *tag // name of the state
      , double const occ // occupation number
      , double const energy // energy eigenvalue or energy parameter
      , char const *csv_class // classification as {core, semicore, valence, ?}
      , int const ir_cut=0 // radial grid index of the augmentation radius
  ) {
      // display stat information about a radial wave functions

      double q{0}, qr{0}, qr2{0}, qrm1{0}, qout{0};
      for (int ir = 0; ir < rg.n; ++ir) {
          double const rho_wf = wave[ir]*wave[ir];
          double const dV = rg.r2dr[ir];
          double const r = rg.r[ir];
          double const r_inv_dV = rg.rdr[ir]; // (1/r)*r^2 dr = r*dr
          q    += rho_wf*dV; // charge
          qr   += rho_wf*r*dV; // for <r>
          qr2  += rho_wf*r*r*dV; // for variance
          qrm1 += rho_wf*r_inv_dV; // Coulomb integral without -Z
          qout += rho_wf*dV*(ir >= ir_cut); // masked
      } // ir
      double const qinv = (q > 0) ? 1./q : 0;
      double const charge_outside = qout*qinv;

      if (echo > 0) {
//        std::printf("# %s %-9s%-4s%6.1f E=%16.6f %s ", label, csv_class, tag, occ, energy*eV, _eV);
          show_state(label, csv_class, tag, occ, energy, ' ');
          if (echo > 4) { // detailed information
              auto const rms = std::sqrt(std::max(0., qr2*qinv));
              std::printf(" <r>=%g rms=%g %s <r^-1>=%g %s\t", qr*qinv*Ang, rms*Ang, _Ang, qrm1*qinv*eV, _eV);
          } // echo
          std::printf((charge_outside > 1e-3) ? "out=%6.2f %%\n"
                                              : "out= %.0e %%\n", 100*charge_outside);
      } // echo

      return charge_outside; // fraction of charge outside the augmentation radius
  } // show_state_analysis

  // ToDo: Would it make sense to define this function?
  // template <int Pseudo> show_state_analysis(echo, label, rg, energy_level<Pseudo> const & state);

