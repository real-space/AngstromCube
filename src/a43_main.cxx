// This file is part of AngstromCube under MIT License

#ifndef NO_UNIT_TESTS
  #include "exchange_correlation.hxx" // ::all_tests
  #include "spherical_harmonics.hxx" // ::all_tests
  #include "conjugate_gradients.hxx" // ::all_tests
  #include "potential_generator.hxx" // ::all_tests
  #include "hermite_polynomial.hxx" // ::all_tests
  #include "global_coordinates.hxx" // ::all_tests
  #include "radial_eigensolver.hxx" // ::all_tests
  #include "boundary_condition.hxx" // ::all_tests
  #include "fermi_distribution.hxx" // ::all_tests
  #include "recorded_warnings.hxx" // ::all_tests
  #include "finite_difference.hxx" // ::all_tests
  #include "radial_integrator.hxx" // ::all_tests
  #include "geometry_analysis.hxx" // ::all_tests
  #include "density_generator.hxx" // ::all_tests
  #include "fourier_transform.hxx" // ::all_tests
  #include "iterative_poisson.hxx" // ::all_tests
  #include "self_consistency.hxx" // ::all_tests
  #include "radial_potential.hxx" // ::all_tests
  #include "bessel_transform.hxx" // ::all_tests
  #include "parallel_domains.hxx" // ::all_tests
  #include "structure_solver.hxx" // ::all_tests
  #include "green_projection.hxx" // ::all_tests
  #include "scattering_test.hxx" // ::all_tests
  #include "davidson_solver.hxx" // ::all_tests
  #include "progress_report.hxx" // ::all_tests
  #include "chemical_symbol.hxx" // ::all_tests
  #include "linear_operator.hxx" // ::all_tests
  #include "sho_hamiltonian.hxx" // ::all_tests
  #include "fourier_poisson.hxx" // ::all_tests
  #include "solid_harmonics.hxx" // ::all_tests
  #include "bisection_tools.hxx" // ::all_tests
  #include "green_parallel.hxx" // ::all_tests
  #include "poisson_solver.hxx" // ::all_tests
  #include "brillouin_zone.hxx" // ::all_tests
  #include "sho_projection.hxx" // ::all_tests
  #include "shift_boundary.hxx" // ::all_tests
  #include "linear_algebra.hxx" // ::all_tests
  #include "grid_operators.hxx" // ::all_tests
  #include "dense_operator.hxx" // ::all_tests
  #include "element_config.hxx" // ::all_tests
  #include "symmetry_group.hxx" // ::all_tests
  #include "complex_tools.hxx" // ::all_tests
  #include "vector_layout.hxx" // ::all_tests
  #include "sho_potential.hxx" // ::all_tests
  #include "pawxml_import.hxx" // ::all_tests
  #include "load_balancer.hxx" // ::all_tests
  #include "simple_stats.hxx" // ::all_tests
  #include "mpi_parallel.hxx" // ::all_tests
  #include "angular_grid.hxx" // ::all_tests
  #include "pseudo_tools.hxx" // ::all_tests
  #include "simple_timer.hxx" // ::all_tests
  #include "sigma_config.hxx" // ::all_tests
  #include "dense_solver.hxx" // ::all_tests
  #include "json_reading.hxx" // ::all_tests
  #include "green_input.hxx" // ::all_tests
  #include "xml_reading.hxx" // ::all_tests
  #include "unit_system.hxx" // ::all_tests
  #include "simple_math.hxx" // ::all_tests
  #include "sho_overlap.hxx" // ::all_tests
  #include "radial_grid.hxx" // ::all_tests
  #include "single_atom.hxx" // ::all_tests
  #include "inline_math.hxx" // ::all_tests
  #include "sho_unitary.hxx" // ::all_tests
  #include "plane_wave.hxx" // ::all_tests
  #include "atom_image.hxx" // ::all_tests
  #include "real_space.hxx" // ::all_tests
  #include "multi_grid.hxx" // ::all_tests
  #include "sho_radial.hxx" // ::all_tests
  #include "sho_tools.hxx" // ::all_tests
  #include "sho_basis.hxx" // ::all_tests
  #include "atom_core.hxx" // ::all_tests
  #include "data_view.hxx" // ::all_tests
  #include "control.hxx" // ::all_tests

  // "green_*.hxx" headers requiring CUDA are included in green_tests.cu/green_tests.cxx
  #include "green_tests.hxx" // ::add_tests
//   #include "green_experiments.hxx" // ::all_tests
//   #include "green_potential.hxx" // ::all_tests
//   #include "green_function.hxx" // ::all_tests
//   #include "green_kinetic.hxx" // ::all_tests
//   #include "green_sparse.hxx" // ::all_tests
//   #include "green_dyadic.hxx" // ::all_tests
//   #include "green_action.hxx" // ::all_tests

#endif // not NO_UNIT_TESTS

#include <cstdlib> // std::abs, ::abort

#include "recorded_warnings.hxx" // warn, ::show_warnings, ::clear_warnings
#include "simple_timer.hxx" // SimpleTimer
#include "unit_system.hxx" // ::set_output_units
#include "control.hxx" // ::command_line_interface, ::get

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

#include <cstdio> // std::printf
#include <cassert> // assert
#include <vector> // std::vector
#include <string> // std::string
#include <tuple> // std::tuple<...>, ::make_tuple, ::get

  status_t run_unit_tests(char const *const module=nullptr, int const echo=0) {

#ifdef  NO_UNIT_TESTS
      error("version was compiled with -D NO_UNIT_TESTS but try to test \'%s\'", module);
      return STATUS_TEST_NOT_INCLUDED;
#else // NO_UNIT_TESTS

      SimpleTimer unit_test_timer(__func__, 0, module, 0); // timer over all tests

      std::string const input_name(module ? module : "");
      bool const show = ('?' == input_name[0]);
      bool const all  = ( 0  == input_name[0]) || show;
      bool const chapters = all && (!show) && (echo > 0); // show chapter separators
      if (echo > 0) {
          if (show) { std::printf("\n# show available module tests:\n"); } else
          if (all)  { std::printf("\n# run all tests!\n\n"); }
          else      { std::printf("\n# run unit tests for module '%s'\n\n", input_name.c_str()); }
      } // echo

      double constexpr ChapterMarker = -1.; // passed in the time entry
      std::vector<std::tuple<char const*, double, status_t>> results;
      { // testing scope

#include "add_module_test.h" // macro definition of add_module_test(MODULE_NAME)

#define   start_a_chapter(chapter_name) {                                           \
              if (chapters) std::printf("\n\n\n\n#\n# %s modules\n#\n\n\n\n", chapter_name); \
              if (all) results.push_back(std::make_tuple(chapter_name, ChapterMarker, 0)); \
          } // start_a_chapter

          start_a_chapter("general"); // *****************************************
          add_module_test(control);
          add_module_test(recorded_warnings);
          add_module_test(simple_stats);
          add_module_test(simple_timer);
          add_module_test(progress_report);
          add_module_test(json_reading);
          add_module_test(xml_reading);
          add_module_test(data_view);

          start_a_chapter("math"); // *****************************************
          add_module_test(inline_math);
          add_module_test(simple_math);
          add_module_test(complex_tools);
          add_module_test(linear_algebra);
          add_module_test(bisection_tools);
          add_module_test(solid_harmonics);
          add_module_test(finite_difference);
          add_module_test(hermite_polynomial);
          add_module_test(spherical_harmonics);
          add_module_test(dense_operator);
          add_module_test(angular_grid);

          start_a_chapter("inactive"); // *****************************************
          add_module_test(vector_layout);
          add_module_test(linear_operator);

          start_a_chapter("input"); // *****************************************
          add_module_test(chemical_symbol);
          add_module_test(boundary_condition);
          add_module_test(geometry_analysis);
          add_module_test(shift_boundary);
          add_module_test(unit_system);
          add_module_test(real_space);

          start_a_chapter("parallelization"); // *****************************************
          add_module_test(mpi_parallel);
          add_module_test(parallel_domains);
          add_module_test(global_coordinates);
          add_module_test(load_balancer);

          start_a_chapter("electrostatics"); // *****************************************
          add_module_test(multi_grid);
          add_module_test(fourier_transform);
          add_module_test(fourier_poisson);
          add_module_test(radial_potential);
          add_module_test(iterative_poisson);
          add_module_test(poisson_solver);

          start_a_chapter("radial"); // *****************************************
          add_module_test(radial_grid);
          add_module_test(radial_integrator);
          add_module_test(radial_eigensolver);
          add_module_test(atom_core);
          add_module_test(sigma_config);
          add_module_test(element_config);
          add_module_test(pseudo_tools);
          add_module_test(bessel_transform);
          add_module_test(pawxml_import);
          add_module_test(scattering_test);
          add_module_test(single_atom);

          start_a_chapter("SHO-specific"); // *****************************************
          add_module_test(sho_tools);
          add_module_test(sho_unitary);
          add_module_test(sho_overlap);
          add_module_test(sho_projection);
          add_module_test(sho_potential);
          add_module_test(sho_radial);
          add_module_test(sho_basis);
          add_module_test(sho_hamiltonian);

          start_a_chapter("Hamiltonian"); // *****************************************
          add_module_test(atom_image);
          add_module_test(plane_wave);
          add_module_test(grid_operators);
          add_module_test(green_input);
          add_module_test(green_parallel);
          add_module_test(green_projection);
//  The 7 green_* modules below are grouped for faster development
//           add_module_test(green_sparse);
//           add_module_test(green_function);
//           add_module_test(green_kinetic);
//           add_module_test(green_potential);
//           add_module_test(green_dyadic);
//           add_module_test(green_action);
//           add_module_test(green_experiments);
//  The 7 green_* modules above are grouped for faster development
          green_tests::add_tests(results, input_name, show, all, echo);

          start_a_chapter("eigensolver"); // *****************************************
          add_module_test(conjugate_gradients);
          add_module_test(davidson_solver);
          add_module_test(dense_solver);
          add_module_test(structure_solver);

          start_a_chapter("DFT-specific"); // *****************************************
          add_module_test(symmetry_group);
          add_module_test(brillouin_zone);
          add_module_test(fermi_distribution);
          add_module_test(exchange_correlation);
          add_module_test(potential_generator);
          add_module_test(density_generator);
          add_module_test(self_consistency);

          start_a_chapter("================"); // last separator
#undef    start_a_chapter
#undef    add_module_test

      } // testing scope

      status_t status(0);
      if (results.size() < 1) { // nothing has been tested
          error("test for '%s' not found, use -t '?' to see available modules!", module);
      } else {
          if (echo > 0) std::printf("\n\n");
          int const show_timings = control::get("timings.show", 0.);
          int nonzero_status{0};
          int nmodules{0}; // count the module without counting chapters
          for (auto result : results) {
              auto const name = std::get<0>(result);
              auto const time = std::get<1>(result);
              auto const stat = std::get<2>(result);
              if (ChapterMarker == time) {
                  if (echo > 1) std::printf("#  ==== %s modules ====\n", name); // this is a chapter marker
              } else {
                  if (echo > 1) {
                      if (show) {
                          std::printf("#    module= %s\n", name);
                      } else {
                          std::printf("#    module= %-24s status= %i", name, int(stat));
                          if (0 != stat)    std::printf(" FAILED");
                          if (show_timings) std::printf(" \ttime=%9.3f seconds", time);
                          std::printf("\n");
                      }
                  } // echo
                  status += std::abs(int(stat)); // we must take |stat| so that positive and negative statuses do not cancel out
                  nonzero_status += (0 != int(stat));
                  ++nmodules;
              } // chapter marker
          } // result
          if (show) {
              if (echo > 0) std::printf("\n# %d modules can be tested\n", nmodules);
              warn("display mode only, none of %d modules has been tested", nmodules);
          } else { // show
              if (nmodules > 1 && echo > 0) {
                  std::printf("\n#%3d modules have been tested,  total status= %d", nmodules, int(status));
                  if (show_timings) std::printf(" \t %13.3f seconds", unit_test_timer.stop()); // total time
                  std::printf("\n\n");
              } // show total status if many modules have been tested
              if (status > 0) warn("Tests for %d module%s failed!", nonzero_status, (1 == nonzero_status)?"":"s");
          } // show
      } // something has been tested
      return status;
#endif // NO_UNIT_TESTS
  } // run_unit_tests

  int show_help(char const *executable) {
      std::printf("Usage %s [OPTION]\n"
        "   --help           [-h]\tThis help message\n"
        "   --version            \tShow version number\n"
#ifndef  NO_UNIT_TESTS
        "   --test <module>  [-t]\tRun module unit test\n"
#endif // NO_UNIT_TESTS
        "   --verbose        [-V]\tIncrement verbosity level\n"
        "   +<name>=<value>      \tModify variable environment\n"
        "\n", executable);
      return 0;
  } // show_help

  int show_version(char const *executable="#", int const echo=0) {
#ifdef    _GIT_KEY
      // stringify the value of a macro, two expansion levels needed
      #define macro2string(a) stringify(a)
      #define stringify(b) #b
      auto const git_key = macro2string(_GIT_KEY);
      #undef  stringify
      #undef  macro2string
      control::set("git.key", git_key); // store in the global variable environment
      if (echo > 0) std::printf("# %s git checkout %s\n\n", executable, git_key);
#endif // _GIT_KEY
      return 0;
  } // show_version

  int main(int const argc, char const *argv[]) {
      if (argc < 2) {
          std::printf("%s: no arguments passed!\n", (argc < 1) ? __FILE__ : argv[0]);
          return -1;
      } // no argument passed to executable
      status_t stat(0);
      char const *test_unit = ""; // the name of the unit to be tested
      int run_tests{0};
      int verbosity{3}; // set default verbosity low
      control::set("executable.name", argv[0]);
      for (int iarg = 1; iarg < argc; ++iarg) {
          assert(nullptr != argv[iarg]);
          char const ci0 = *argv[iarg]; // char #0 of command line argument #iarg
          if ('-' == ci0) {

              // options (short or long)
              char const ci1 = *(argv[iarg] + 1); // char #1 of command line argument #iarg
              char const IgnoreCase = 'a' - 'A'; // use with | to convert upper case chars into lower case chars
              if ('-' == ci1) {

                  // long options with "--"
                  std::string option(argv[iarg] + 2); // + 2 to remove "--" in front
                  if ("help" == option) {
                      return show_help(argv[0]);
                  } else
                  if ("version" == option) {
                      return show_version(argv[0], 1);
                  } else
                  if ("verbose" == option) {
                      verbosity = 6; // set verbosity high
                  } else
                  if ("test" == option) {
                      ++run_tests; if (iarg + 1 < argc) test_unit = argv[iarg + 1];
                  } else {
                      ++stat; warn("ignored unknown command line option --%s", option.c_str());
                  } // option

              } else { // ci1

                  // short options with "-"
                  if ('h' == (ci1 | IgnoreCase)) {
                      return show_help(argv[0]);
                  } else
                  if ('v' == (ci1 | IgnoreCase)) {
                      verbosity += 1 + 3*('V' == ci1); // increment by 'V':4, 'v':1
                  } else
                  if ('t' == (ci1 | IgnoreCase)) {
                      ++run_tests; if (iarg + 1 < argc) test_unit = argv[iarg + 1];
                  } else {
                      ++stat; warn("ignored unknown command line option -%c", ci1);
                  } // ci1

              } // ci1

          } else // ci0
          if ('+' == ci0) {
              stat += control::command_line_interface(argv[iarg] + 1, iarg); // start after the '+' char
          } else
          if (argv[iarg] != test_unit) {
              ++stat; warn("ignored command line argument \'%s\'", argv[iarg]);
          } // ci0

      } // iarg
      //
      if (verbosity > 0) {
          std::printf("\n#");
          for (int iarg = 0; iarg < argc; ++iarg) {
              std::printf(" %s", argv[iarg]); // repeat all command line arguments for completeness of the log file
          } // iarg
          std::printf("\n");
      } // verbosity
      //
      // in addition to command_line_interface, we can modify the control environment by a file
      stat += control::read_control_file(control::get("control.file", ""), verbosity);
      //
      int const echo = control::get("verbosity", double(verbosity)); // verbosity may have been defined in the control file
      //
      stat += show_version(argv[0], echo);
      //
      if (echo > 0) std::printf("\n# verbosity=%d\n", echo);
      //
      stat += unit_system::set(control::get("output.length.unit", "Bohr"),
                               control::get("output.energy.unit", "Ha"), echo);
      // run
      if (run_tests) {
          stat += run_unit_tests(test_unit, echo);
      } // run_tests

      // finalize
      {   int const control_show = control::get("control.show", 0.); // 0:show none, 1:show used, 2:show unused, 4:show defaults
          if (echo > 3) std::printf("\n# control.show=%d     0:none 1:used 2:unused 4:defaults\n", control_show);
          if (control_show && echo > 0) {
              stat += control::show_variables(control_show);
          }
      } // show all variable names defined in the control environment

      if (echo > 0) recorded_warnings::show_warnings(3);
      recorded_warnings::clear_warnings(1);
      return int(stat);
  } // main
