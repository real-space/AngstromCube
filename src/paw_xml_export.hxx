#pragma once

#include <cstdio> // printf, std::fprintf, std::FILE, std::fopen, std::fclose, std::snprintf

#include "chemical_symbol.hxx" // ::get
#include "radial_grid.h" // radial_grid_t
#include "energy_level.hxx" // partial_wave_t, TRU, SMT, TRU_AND_SMT

namespace paw_xml_export {

  template <class PartialWave>
  int write_to_file(
        double const Z // number of protons in the core
      , radial_grid_t const rg[TRU_AND_SMT] // TRU and SMT
      , std::vector<PartialWave> const & valence_states
      , char const valence_states_active[]
      , view3D<double> const & kinetic_energy_differences // (ts,nln,nln)
      , double const n_electrons[3] // core, semicore, valence
      , view2D<double> const spherical_density[TRU_AND_SMT] // spherical_density[ts](csv,ir), ToDo: extend to have KIN for meta-GGAs
      , view2D<double> const & projector_functions // ==r[ir]*Projectors(nln,ir)
      , double const r_cut=2
      , double const sigma_cmp=.5
      , double const zero_potential[]=nullptr
      , int const echo=0 // log-level
      , char const *filename=nullptr // filename, nullptr: use a default name
      , char const *pathname="."
  ) { 
      double constexpr Y00 = .28209479177387817; // == 1/sqrt(4*pi)
      double const E[8] = {0, 0, 0, 0,   0, 0, 0, 0};
      char const ts_label[TRU_AND_SMT][8] = {"ae", "pseudo"};
      char const ts_tag[TRU_AND_SMT][4] = {"ae", "ps"};

      char Sy[4], file_name_buffer[512];
      int const iZ = chemical_symbol::get(Sy, Z);
      if (nullptr == filename) {
          std::snprintf(file_name_buffer, 511, "%s/%s.xml", pathname, Sy);
          filename = file_name_buffer;
      } // generate a default file name
      if (echo > 0) printf("# %s Sy=%s Z=%g iZ=%d filename=%s\n", __func__, Sy, Z, iZ, filename);
      char const *const config_string=Sy; // ToDo: insert the config string here

      std::FILE *const f = std::fopen(filename, "w");
      if (nullptr == f) {
          if (echo > 1) printf("# %s Error opening file %s!\n", __func__, filename);
          return __LINE__;
      } // failed to open

      // XML file header
      std::fprintf(f, "<?xml version=\"%.1f\"?>\n", 1.0);
      std::fprintf(f, "<paw_setup version=\"%.1f\">\n", 0.6);
      std::fprintf(f, "  <!-- Z=%g %s setup for the Projector Augmented Wave method. -->\n", Z, Sy);
      std::fprintf(f, "  <!-- Units: Hartree and Bohr radii.                      -->\n");
      

      std::fprintf(f, "  <atom symbol=\"%s\" Z=\"%g\" core=\"%g\" semicore=\"%g\" valence=\"%g\"/>\n",
                                        Sy,  Z,  n_electrons[0], n_electrons[1], n_electrons[2]);
      std::fprintf(f, "  <xc_functional type=\"LDA\" name=\"PZ81\"/>\n");
      std::fprintf(f, "  <generator type=\"scalar-relativistic\" name=\"Koelling-Harmon\">\n");
      std::fprintf(f, "     %s\n  </generator>\n", config_string);
      std::fprintf(f, "  <ae_energy kinetic=\"%.6f\" xc=\"%.6f\"\n             "
                          "electrostatic=\"%.6f\" total=\"%.6f\"/>\n", E[3], E[2], E[1], E[0]);
      std::fprintf(f, "  <core_energy kinetic=\"%.6f\"/>\n", E[4]);
      
      std::fprintf(f, "  <valence_states>\n");
      for(size_t iln = 0; iln < valence_states.size(); ++iln) {
          if (valence_states_active[iln]) {
              auto const & vs = valence_states[iln];
              std::fprintf(f, "    <state");
              std::fprintf(f, " n=\"%d\"", vs.enn); // ToDo: do not show principal quantum number for excited states
              std::fprintf(f, " l=\"%d\"", vs.ell);
              if (vs.occupation > 1e-24) std::fprintf(f, " f=\"%g\"", vs.occupation);
              std::fprintf(f, " rc=\"%.3f\" e=\"%9.6f\" id=\"%s-%s\"/>\n", r_cut, vs.energy, Sy,vs.tag);
          } // active
      } // iln
      std::fprintf(f, "  </valence_states>\n");
      
      for(int ts = TRU; ts <= SMT; ++ts) {
          double const prefactor = rg[ts].rmax/(std::exp(rg[ts].anisotropy * (rg[ts].n - 1)) - 1.);
          std::fprintf(f, "  <radial_grid eq=\"r=f*(exp(d*i)-1)\""
              " d=\"%g\" f=\"%g\" n=\"%d\" istart=\"%d\" iend=\"%d\" id=\"g_%s\"/>\n",
              rg[ts].anisotropy, prefactor, rg[ts].n - 1, 1, rg[ts].n - 1, ts_tag[ts]);
      } // ts
      double const rc_cmp = sigma_cmp*std::sqrt(2.); // ToDo: check formula
      std::fprintf(f, "  <shape_function type=\"gauss\" rc=\"%.12e\"/>\n", rc_cmp);

      auto constexpr m = 1; // 0.005; // DEBUG, scale the number of radial grid entries down to say one line

      // core electron density
      if (n_electrons[0] > 0) {
          int constexpr csv = 0; // 0:core
          for(int rk = 0; rk < 1; ++rk) { // rk==0 density, rk==1 kinetic energy density, ToDo
              auto const rk_tag = rk ? "_kinetic_energy" : "";
              for(int ts = TRU; ts <= SMT; ++ts) {
                  std::fprintf(f, "  <%s_core%s_density grid=\"g_%s\">\n    ", ts_label[ts], rk_tag, ts_tag[ts]);
                  for(int ir = 1; ir < rg[ts].n*m; ++ir) {
                      std::fprintf(f, " %.12e", spherical_density[ts](csv,ir)*Y00);
                  } // ir
                  std::fprintf(f, "\n  </%s_core%s_density>\n", ts_label[ts], rk_tag);
              } // ts
          } // rk
      } // core electrons present

      // spherical valence density
      if (n_electrons[2] > 0) {
          int constexpr csv = 2; // 2:valence
          {   auto const ts = SMT; // (smooth only)
              std::fprintf(f, "  <%s_valence_density grid=\"g_%s\">\n    ", ts_label[ts], ts_tag[ts]);
              for(int ir = 1; ir < rg[ts].n*m; ++ir) {
                  std::fprintf(f, " %.12e", spherical_density[ts](csv,ir)*Y00);
              } // ir
              std::fprintf(f, "\n  </%s_valence_density>\n", ts_label[ts]);
          } // ts
      } // valence electrons presents

      if (zero_potential != nullptr) { // scope: zero potential
          auto const ts = SMT;
          std::fprintf(f, "  <zero_potential grid=\"g_%s\">\n    ", ts_tag[ts]);
          for(int ir = 1; ir < rg[ts].n*m; ++ir) {
              std::fprintf(f, " %.12e", zero_potential[ir]);
          } // ir
          std::fprintf(f, "\n  </zero_potential>\n");
      } // scope

      for(size_t iln = 0; iln < valence_states.size(); ++iln) {
          if (valence_states_active[iln]) {
              auto const & vs = valence_states[iln];
              for(int ts = TRU; ts <= SMT; ++ts) {
                  std::fprintf(f, "  <%s_partial_wave state=\"%s-%s\" grid=\"g_%s\">\n    ", ts_label[ts], Sy,vs.tag, ts_tag[ts]);
                  if (vs.wave[ts] == nullptr) return __LINE__; // error
                  for(int ir = 1; ir < rg[ts].n*m; ++ir) {
                      std::fprintf(f, " %.12e", vs.wave[ts][ir]);
                  } // ir
                  std::fprintf(f, "\n  </%s_partial_wave>\n", ts_label[ts]);
              } // ts
              {   auto const ts = SMT;
                  std::fprintf(f, "  <projector_function state=\"%s-%s\" grid=\"g_%s\">\n    ", Sy,vs.tag, ts_tag[ts]);
                  for(int ir = 1; ir < rg[ts].n*m; ++ir) {
                      std::fprintf(f, " %.12e", projector_functions(iln,ir)*rg[ts].rinv[ir]);
                  } // ir
                  std::fprintf(f, "\n  </projector_function>\n");
              } // ts
          } // active
      } // iln


      { // scope: kinetic energy differences
          std::fprintf(f, "  <kinetic_energy_differences>\n");
          for(size_t iln = 0; iln < valence_states.size(); ++iln) {
              if (valence_states_active[iln]) {
                  std::fprintf(f, "    ");
                  for(size_t jln = 0; jln < valence_states.size(); ++jln) {
                      if (valence_states_active[jln]) {
                          std::fprintf(f, " %.12e", kinetic_energy_differences(TRU,iln,jln)
                                                  - kinetic_energy_differences(SMT,iln,jln));
                      } // active
                  } // jln
                  std::fprintf(f, " \n");
              } // active
          } // iln
          std::fprintf(f, "  </kinetic_energy_differences>\n");
      } // scope
      
      std::fprintf(f, "  <!-- exact_exchange_X_matrix not included -->\n");
      std::fprintf(f, "  <exact_exchange core-core=\"0\"/>\n");

      // XML file tail
      std::fprintf(f, "</paw_setup>\n");
      std::fclose(f);
      
      if (echo > 3) printf("# file %s written\n", filename);
      return 0; // 0:success
  } // write_to_file

} // namespace paw_xml_export
