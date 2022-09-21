#pragma once

#include <cassert> // assert
#include <cstdio> // std::printf, ::snprintf
#include <vector> // std::vector
#include <algorithm> // std::min, ::max

#include "status.hxx" // status_t
#include "control.hxx" // ::get
#include "atom_core.hxx" // ::guess_energy, ::nl_index
#include "quantum_numbers.h" // enn_QN_t
#include "radial_grid.hxx" // radial_grid_t
#include "sigma_config.hxx" // ::element_t
#include "spherical_state.hxx" // core, semicore, valence, csv_undefined, csv_name,
//       "spherical_state.hxx" // show_state_analysis, show_state
#include "chemical_symbol.h" // element_symbols
#include "recorded_warnings.hxx" // warn
#include "radial_eigensolver.hxx" // ::shooting_method
#include "inline_math.hxx" // set
#include "sho_tools.hxx" // ::nn_max

namespace element_config {

  inline char const* element_symbol(int const Z) {
      int const z128 = Z & 127; // == Z % 128
      return &(element_symbols[2*z128]); // this and the following char belong to this element
  } // element_symbol

  char const ellchar[] = "spdfgh+";

  inline sigma_config::element_t const & get(
        double const Z_core // nuclear charge
      , double const ionization
      , radial_grid_t const & rg // TRU radial grid descriptor
      , double const rV[]   // TRU radial potential r*V(r)
      , char const *const Sy="X"
      , double const core_state_localization=-1
      , int const echo=3 // log-level
      , int const SRA=1 // 1: Scalar Relativistic Approximation
  ) {
      if (echo > 0) std::printf("# %s element_config for Z= %g\n", Sy, Z_core);

      auto & e = *(new sigma_config::element_t); // is this memory ever released?

      e.Z = Z_core;

      char Sy_config[32];
      std::snprintf(Sy_config, 31, "element_%s.rcut", Sy);
      auto const rcut_default = control::get("element_config.rcut", 2.0);
      e.rcut = control::get(Sy_config, rcut_default); // augmentation radius in Bohr

      std::snprintf(Sy_config, 31, "element_%s.sigma", Sy);
      auto const sigma_default = control::get("element_config.sigma", 0.5);
      e.sigma = control::get(Sy_config, sigma_default); // spread for projectors in Bohr

      std::snprintf(Sy_config, 31, "element_%s.numax", Sy);
      auto const numax_default = control::get("element_config.numax", 3.);
      e.numax = int(control::get(Sy_config, numax_default));

      for (int ell = 0; ell < 8; ++ell) {
          e.nn[ell] = std::max(0, sho_tools::nn_max(e.numax, ell));
      } // ell

      double hole_charge{0}, hole_charge_used{0};
      int inl_hole{-1}, ell_hole{-1}; // invalid
      std::snprintf(Sy_config, 31, "element_%s.hole.enn", Sy);
      int const enn_hole = control::get(Sy_config, 0.);
      if (enn_hole > 0) {
          if (enn_hole < 9) {
              std::snprintf(Sy_config, 31, "element_%s.hole.ell", Sy);
              ell_hole = control::get(Sy_config, -1.);
              if (ell_hole >= 0 && ell_hole < enn_hole) {
                  inl_hole = atom_core::nl_index(enn_hole, ell_hole);
                  std::snprintf(Sy_config, 31, "element_%s.hole.charge", Sy);
                  hole_charge = std::min(std::max(0.0, control::get(Sy_config, 1.)), 2.*(2*ell_hole + 1));
              } else warn("%s=%d is out of range [0, %d]", Sy_config, ell_hole, enn_hole - 1);
          } else warn("%s=%d is too large", Sy_config, enn_hole);
      }

      set(e.method, 16, '\0'); // clear
      auto const method_default = control::get("element_config.method", "sinc");
      std::snprintf(Sy_config, 31, "element_%s.method", Sy);
      std::snprintf(e.method, 15, "%s", control::get(Sy_config, method_default));

      set(e.occ[0], 32*2, 0.0); // clear occupation numbers

      double core_valence_separation{0}, core_semicore_separation{0}, semi_valence_separation{0};
      if (core_state_localization > 0) {
          if (echo > 0) std::printf("# %s use core state localization criterion with %g %%\n", Sy, core_state_localization*100);
      } else {
          core_valence_separation  = control::get("element_config.core.valence", -2.0); // in Hartree always
          core_semicore_separation = control::get("element_config.core.semicore",    core_valence_separation);
          semi_valence_separation  = control::get("element_config.semicore.valence", core_valence_separation);
          if (core_semicore_separation > semi_valence_separation) {
              warn("%s element_config.core.semicore=%g %s may not be higher than ~.semicore.valence=%g %s, correct for it",
                    Sy, core_semicore_separation*eV, _eV, semi_valence_separation*eV, _eV);
              core_semicore_separation = semi_valence_separation; // correct -. no semicore states possible
          } // warning
      } // criterion

      int const ir_cut = radial_grid::find_grid_index(rg, e.rcut);
      if (echo > 6) std::printf("# %s cutoff radius %g %s at grid point %d of max. %d\n",
                                    Sy, e.rcut*Ang, _Ang, ir_cut, rg.n);

      std::vector<int8_t> as_valence(40, -1);
      std::vector<enn_QN_t> enn_core_ell(8, 0); // energy quantum number of the highest occupied core level
      std::vector<double> wave(rg.n, 0.0); // get memory for the true radial wave functions
      std::vector<double> r2rho(rg.n, 0.0);
      double csv_charge[3] = {0, 0, 0};

      { // scope:
          int highest_occupied_state_index{-1};
          int ics{0}; // init counter for core states

          // loop through the entire configuration of elements up to Z=120
          double n_electrons{Z_core + ionization}; // init number of electrons to be distributed
          for (int nq_aux = 0; nq_aux < 8; ++nq_aux) {    // auxiliary quantum number, allows Z up to 120
              enn_QN_t enn = (nq_aux + 1)/2;              // principal quantum number n
              for (int ell = nq_aux/2; ell >= 0; --ell) { // orbital angular momentum l
                  ++enn; // update principal quantum number
//                for (int jj = 2*ell; jj >= 2*ell; jj -= 2) // total angular momentum j
                  {   int const jj = 2*ell;
                      double const max_occ = 2*(jj + 1); // largest state occupation

                      char tag[4]; std::snprintf(tag, 3, "%d%c", enn, ellchar[ell]);
                      set(r2rho.data(), rg.n, 0.0); // clear

                      double E{atom_core::guess_energy(Z_core, enn)}; // init with a guess
                      // solve the eigenvalue problem with a spherical potential
                      radial_eigensolver::shooting_method(SRA, rg, rV, enn, ell, E, wave.data(), r2rho.data());

                      int const inl = atom_core::nl_index(enn, ell);

                      double const hole = hole_charge*(inl_hole == inl);
                      double const occ_no_hole = std::min(std::max(0., n_electrons),        max_occ);
                      double const occ         = std::min(std::max(0., occ_no_hole - hole), max_occ);
                      double const real_hole_charge = occ_no_hole - occ;
                      if (real_hole_charge > 0) {
                          hole_charge_used = real_hole_charge;
                          if (echo > 1) std::printf("# %s %s has an occupation hole of element_%s.hole.charge=%g electrons\n", 
                                                       Sy, tag,                                Sy, real_hole_charge);
                          assert(enn_hole == enn && ell_hole == ell);
                      } // core hole active

                      int csv{csv_undefined};
                      auto const charge_outside = show_state_analysis(echo, Sy, rg, wave.data(), tag, occ, E, "?", ir_cut);
                      if (core_state_localization > 0) {
                          // criterion based on the charge outside the sphere
                          if (charge_outside > core_state_localization) {
                              csv = valence; // mark as valence state
                          } else {
                              csv = core; // mark as core state
                          } // stay in the core
                      } else { // criterion
                          // energy criterions
                          if (E > semi_valence_separation) {
                              csv = valence; // mark as valence state
                          } else if (E > core_semicore_separation) {
                              csv = semicore; // mark as semicore state
                          } else { // move to the semicore band
                              csv = core; // mark as core state
                          } // stay in the core
                      } // criterion
                      assert(csv_undefined != csv);
                      if (echo > 15) std::printf("# as_%s[nl_index(enn=%d, ell=%d) = %d] = %d\n", csv_name(csv), enn, ell, inl, ics);

                      if (valence == csv) as_valence[inl] = ics; // mark as good for the valence band, store the core state index

                      if (occ > 0) {
                          highest_occupied_state_index = ics; // store the index of the highest occupied core state
                          if (echo > 5) show_state(Sy, csv_name(csv), tag, occ, E);
                          if (as_valence[inl] < 0) {
                              enn_core_ell[ell] = std::max(enn, enn_core_ell[ell]); // find the largest enn-quantum number of the occupied core states
                          } // not as valence
                          csv_charge[csv] += occ;

                          double const has_norm = dot_product(rg.n, r2rho.data(), rg.dr);
                          if (has_norm <= 0) {
                              warn("%s %i%c-state cannot be normalized! integral= %g electrons", Sy, enn, ellchar[ell], has_norm);
                          } // cannot be normalized

                          e.occ[inl][0] = e.occ[inl][1] = 0.5*occ; // split occupation number between dn-spin and up-spin
                          e.csv[inl] = csv;
                      } // occupied

                      n_electrons -= occ; // subtract as many electrons as have been assigned to this state
                      ++ics;
                  } // jj
              } // ell
          } // nq_aux
          int const nstates = highest_occupied_state_index + 1; // correct the number of core states to those occupied
          if (echo > 0) std::printf("# %s found %d spherical states\n", Sy, nstates);

          if (n_electrons > 0) warn("%s after distributing %g, %g electrons remain",
                                     Sy, Z_core - ionization, n_electrons);

          auto const total_n_electrons = csv_charge[core] + csv_charge[semicore] + csv_charge[valence];
          if (echo > 2) std::printf("# %s initial occupation with %g electrons: %g core, %g semicore and %g valence electrons\n", 
                                  Sy, total_n_electrons, csv_charge[core], csv_charge[semicore], csv_charge[valence]);

          if (inl_hole >= 0) {
              auto const diff = hole_charge_used - hole_charge;
              if (std::abs(diff) > 5e-16) {
                  warn("hole.charge=%g requested in %s-%d%c state but used %g electrons (difference %.1e)",
                        hole_charge, Sy, enn_hole,ellchar[ell_hole], hole_charge_used, diff);
              } // deviation
              if (echo > 0) std::printf("# %s occupation hole of %g electrons in the %d%c state\n",
                                           Sy, hole_charge_used, enn_hole,ellchar[ell_hole]);
          } // warning when deviates

      } // scope
//    set(e.ncmx, 4, enn_core_ell.data());

      return e;
  } // get

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_show_all_elements(int const echo=2) {
      if (echo > 0) std::printf("# %s:  %s\n", __FILE__, __func__);
      float iocc[32];

      for (int spin = 1; spin >= 0; --spin) { // 1:with and 0:without spin-orbit coupling
          auto const echo_occ = 9; // how high must the log-level be to display also the occupation numbers?
          if (echo > 2 + 3*spin) { // how high must the log-level be to display also the spin-orbit table?
              if (spin > 0) std::printf("# including spin-orbit\n");
              std::printf("\n#  nl    j   occ    index\n");
              for (int inl = 0; inl < 32; ++inl) iocc[inl] = 0; // clear array
              int iZ{0}; // atomic number
              int ishell{0}; // shell index
              for (int nq_aux = 0; nq_aux < 8; ++nq_aux) {    // auxiliary quantum number
                  int enn{(nq_aux + 1)/2};                    // principal quantum number n
                  for (int ell = nq_aux/2; ell >= 0; --ell) { // orbital angular momentum l
                      ++enn;
                      for (int jj = 2*ell + spin; jj >= std::max(0, 2*ell - spin); jj -= 2) { // quantum number j = l + s, jj=2*j
                          std::printf("%4d%c%6.1f%4d%9d    %c", enn, ellchar[ell], jj*.5f,
                                      (2 - spin)*(jj + 1), ishell, (echo > echo_occ)?'\n':' ');
                          for (int mm = -jj; mm <= jj; mm += 2) { // angular momentum z-component emm, mm=2*emm
                              for (int s = 0; s <= 1 - spin; ++s) { // when 0==spin, this loop fills each state with 2 electrons
                                  iocc[ishell] += 1;
                                  ++iZ; // next element
                                  auto const El = element_symbol(iZ);
                                  if (echo > echo_occ) {
                                      std::printf("# %c%c occupation", El[0], El[1]);
                                      for (int inl = 0; inl <= ishell; ++inl) {
                                          std::printf(" %g", iocc[inl]);
                                      } // inl
                                      std::printf("\n");
                                  } else {
                                      std::printf(" %c%c", El[0], El[1]);
                                  } // echo_occ
                              } // s
                          } // mm
                          ++ishell; // next shell
                          std::printf("\n");
                      } // jj
                  } // ell
              } // nq_aux
              double sum_occ{0};
              for (int inl = 0; inl < ishell; ++inl) {
                  sum_occ += iocc[inl];
              } // inl
              std::printf("# sum(occ) = %.3f\n", sum_occ);
              assert(120 == sum_occ); // sanity check
              assert(20 + 12*spin == ishell); // sanity check
          } // echo > echos
      } // spin
      if (echo > 0) std::printf("# table of elements according to Aco Z. Muradjan\n\n");
      return 0;
  } // test_show_all_elements

//  core electron configurations for predefined cores:
//
//  rare gases        [frozen d]        [frozen f]
//  __ '    '   0
//  He '1   '   2
//  Ne '22  '  10
//  Ar '33  '  18     3d '333 '  28
//  Kr '443 '  36     4d '444 '  46
//  Xe '554 '  54     5d '5554'  78     4f '5544'  68
//  Rn '6654'  86     6d '6665' 110     6d '6655' 100
//  uo '7765' 118
//
//   1s
//   2s 2p
//   3s 3p 3d
//   4s 4p 4d 4f
//   5s 5p 5d 5f
//   6s 6p 6d
//   7s 7p
//   8s

  inline status_t all_tests(int const echo=0) {
      return test_show_all_elements(echo);
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace element_config
