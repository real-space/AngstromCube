#pragma once

#include <cstdio> // std::printf, ::fprintf, ::fopen
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <algorithm> // std::max, ::min
#include <time.h> // std::time, ::localtime, ::strftime

#include "status.hxx" // status_t

// #include "inline_math.hxx" // set, pow2
#include "recorded_warnings.hxx" // warn
// #include "constants.hxx" // ::pi
// #include "print_tools.hxx" // printf_vector
#include "control.hxx" // ::get

namespace export_upf {

  inline char c4(int const i) { return (i & 0x3) ? ' ' : '\n'; }

  inline status_t write_to_file(
      // , double const Z=3.0
      // , int const n_proj=4
      // , int const n_mesh=20
      // , int const n_wfc=2
      // , int const l_max=1
      double const Z // number of protons in the core
    , radial_grid_t const rg[TRU_AND_SMT] // TRU and SMT radial grid
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
    , double const core_kinetic_energy=0
    , double const kinetic_energy=0
    , double const xc_energy=0
    , double const es_energy=0
    , double const total_energy=0
    , char const *custom_configuration_string=nullptr
    , char const *xc_functional="LDA"
    , char const *pathname="." // does not work
    , char const *filename=nullptr // nullptr: use a default name
  ) {
      status_t stat(0);

      int const l_max = 1;
      int const n_wfc = 2;
      int const n_mesh = rg[SMT].n;
      int const n_proj =

      auto const *const filename = control::get("export_upf.filename", "X.upf");
      if (nullptr == filename) error("needs a filename", 0);
      if ('\0' == *filename) error("needs a valid filename", 0);

      if (echo > 1) std::printf("# try to open file %s for writing.\n", filename);
      auto *const f = std::fopen(filename, "w");
      if (nullptr == f) {
          if (echo > 0) std::printf("# %s Error opening file %s for writing!\n", __func__, filename);
          return __LINE__;
      } // failed to open

      std::printf("# generate UPF according to pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format\n");
      std::fprintf(f, "<UPF version=\"2.0.1\">\n");
          std::fprintf(f, "<PP_INFO>\n");
              std::fprintf(f, "<PP_INPUTFILE>\n");
                  auto const config_string = "";
                  std::fprintf(f, "%s\n", config_string);
              std::fprintf(f, "</PP_INPUTFILE>\n");
          std::fprintf(f, "</PP_INFO>\n");
          std::fprintf(f, "<!--                               -->\n");
          std::fprintf(f, "<!-- END OF HUMAN READABLE SECTION -->\n");
          std::fprintf(f, "<!--                               -->\n");
          std::fprintf(f, "<PP_HEADER");
              std::fprintf(f, "\ngenerated=\"libliveatom %s\"", control::get("git.key", ""));
              std::fprintf(f, "\nauthor=\"Paul Baumeister\"");
              auto const now = std::time(0);
              char date[80]; std::strftime(date, sizeof(date), "%Y-%m-%d", std::localtime(&now));
              std::fprintf(f, "\ndate=\"%s\"", date);
              std::fprintf(f, "\ncomment=\"%s\"", "");
              std::fprintf(f, "\nelement=\"%s\"", "X");
              std::fprintf(f, "\npseudo_type=\"NC\""); // should be USPP and PAW in the future
              std::fprintf(f, "\nrelativistic=\"scalar\"");
              std::fprintf(f, "\nis_ultrasoft=\"F\"");
              std::fprintf(f, "\nis_paw=\"F\"");
              std::fprintf(f, "\nis_coulomb=\"F\"");
              std::fprintf(f, "\nhas_so=\"F\""); // SO:spin-orbit
              std::fprintf(f, "\nhas_wfc=\"F\"");
              std::fprintf(f, "\nhas_gipaw=\"F\""); // gauge-invariant PAW
              std::fprintf(f, "\ncore_correction=\"F\""); // non-linear core correction
              std::fprintf(f, "\nfunctional=\"%s\"", xc_functional); // density functional used
              std::fprintf(f, "\nz_valence=\"%.3f\"", 3.0);
              std::fprintf(f, "\ntotal_psenergy=\"%g\"", 0.0);
              std::fprintf(f, "\nrho_cutoff=\"%g\"", 0.0);
              std::fprintf(f, "\nl_max=\"%d\"", l_max);
              std::fprintf(f, "\nl_local=\"-1\"");
              std::fprintf(f, "\nmesh_size=\"%d\"", n_mesh);
              std::fprintf(f, "\nnumber_of_wfc=\"%d\"", n_wfc);
              std::fprintf(f, "\nnumber_of_proj=\"%d\"", n_proj);
          std::fprintf(f, "/>\n"); // end of PP_HEADER
          std::fprintf(f, "<PP_MESH>\n");
              std::fprintf(f, "<PP_R type=\"real\" size=\"%d\" columns=\"4\">", n_mesh);
              for (int ir = 0; ir < n_mesh; ++ir) {
                  std::fprintf(f, "%c%.12g", c4(ir), 0.0); // r
              } // ir
              std::fprintf(f, "\n</PP_R>\n");
              std::fprintf(f, "<PP_RAB type=\"real\" size=\"%d\" columns=\"4\">", n_mesh);
              for (int ir = 0; ir < n_mesh; ++ir) {
                  std::fprintf(f, "%c%.12g", c4(ir), 0.0); // dr
              } // ir
              std::fprintf(f, "\n</PP_RAB>\n");
          std::fprintf(f, "</PP_MESH>\n");

          std::fprintf(f, "<PP_LOCAL type=\"real\" size=\"%d\" columns=\"4\">", n_mesh);
          for (int ir = 0; ir < n_mesh; ++ir) {
              std::fprintf(f, "%c%.12g", c4(ir), 0.0); // dr
          } // ir
          std::fprintf(f, "\n</PP_LOCAL>\n");

          std::fprintf(f, "<PP_NONLOCAL>\n");
          for (int ibeta = 0; ibeta < n_proj; ++ibeta) {
              std::fprintf(f, "<PP_BETA.%d type=\"real\" size=\"%d\" columns=\"4\"\n", ibeta+1, n_mesh);
              std::fprintf(f, "index=\"%i\" angular_momentum=\"%d\"\n", ibeta+1, 0);
              std::fprintf(f, "cutoff_radius_index=\"%i\" cutoff_radius=\"%g\">",
                              radial_grid::find_grid_index(rg[SMT], r_cut), r_cut);
              for (int ir = 0; ir < n_mesh; ++ir) {
                  std::fprintf(f, "%c%.12g", c4(ir), 0.0); // dr
              } // ir
              std::fprintf(f, "\n</PP_BETA.%d>\n", ibeta+1);
          } // ibeta
              std::fprintf(f, "<PP_DIJ type=\"real\" size=\"%d\" columns=\"4\">", n_proj*n_proj);
              for (int ij = 0; ij < n_proj*n_proj; ++ij) {
                  std::fprintf(f, "%c%.12g", c4(ij), 0.0);
              } // ij
              std::fprintf(f, "\n</PP_DIJ>\n");
          std::fprintf(f, "</PP_NONLOCAL>\n");

          if (echo > 5) std::printf("# not implemented: PP_PSWFC, PP_RHOATOM\n");

      std::fprintf(f, "</UPF>\n");
      std::fclose(f);
      if (echo > 1) std::printf("# file %s written.\n", filename);
      return 0;
  } // write_to_file

//  inline status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }

} // namespace export_upf
