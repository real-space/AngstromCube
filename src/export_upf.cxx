#include <cstdio> // std::printf, ::fprintf, ::fopen
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <cstdint> // size_t, uint32_t
#include <algorithm> // std::max, ::min, ::swap
#include <time.h> // localtime

#include "export_upf.hxx"

#include "inline_math.hxx" // set, pow2
#include "recorded_warnings.hxx" // warn
#include "constants.hxx" // ::pi
#include "print_tools.hxx" // printf_vector
#include "control.hxx" // ::get

namespace export_upf {

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline char c4(int const i) { return (i & 0x3) ? ' ' : '\n'; }

  status_t test1(int const echo=0
      , double const Z=3.0
      , int const n_proj=4
      , int const n_mesh=20
      , int const n_wfc=2
      , int const l_max=1
  ) {
      status_t stat(0);

      auto const *const filename = control::get("export_upf.filename", "X.upf");
      if (nullptr == filename) error("needs a filename", 0);
      if ('\0' == *filename) error("needs a valid filename", 0);

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
              std::fprintf(f, "\nfunctional=\"%s\"", "PZ"); // density functional used
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
              std::fprintf(f, "cutoff_radius_index=\"%i\" cutoff_radius=\"%g\">", n_mesh/2, 1.0);
              for (int ir = 0; ir < n_mesh; ++ir) {
                  std::fprintf(f, "%c%.12g", c4(ir), 0.0); // dr
              } // ir
              std::fprintf(f, "\n</PP_BETA.%d>\n", ibeta+1);
          } // ibeta
              std::fprintf(f, "<PP_DIJ type=\"real\" size=\"%d\" columns=\"4\">", n_proj*nproj);
              for (int ij = 0; ij < n_proj*n_proj; ++ij) {
                  std::fprintf(f, "%c%.12g", c4(ij), 0.0);
              } // ij
              std::fprintf(f, "</PP_DIJ>\n");
          std::fprintf(f, "</PP_NONLOCAL>\n");

      std::fprintf(f, "</UPF>\n");
      std::fclose(f);
      if (echo > 3) std::printf("# file %s written\n", filename);
      return 0;
  } // test1

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test1(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace export_upf
