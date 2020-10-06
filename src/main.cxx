#include <cstdio> // printf
#include <cassert> // assert
#include <cstdlib> // std::abs
#include <vector> // std::vector
#include <string> // std::string
#include <utility> // std::pair

#include "recorded_warnings.hxx" // warn, ::show_warnings, ::clear_warnings
#include "simple_timer.hxx" // SimpleTimer
#include "unit_system.hxx" // ::set_input_units, ::set_output_units
#include "control.hxx" // ::command_line_interface

#include "status.hxx" // status_t

#ifndef NO_UNIT_TESTS
#include "recorded_warnings.hxx" // ::all_tests
#include "finite_difference.hxx" // ::all_tests
#include "hermite_polynomial.hxx" // ::all_tests
#include "spherical_harmonics.hxx" // ::all_tests
#include "conjugate_gradients.hxx" // ::all_tests
#include "potential_generator.hxx" // ::all_tests
#include "fermi_distribution.hxx" // ::all_tests
#include "radial_eigensolver.hxx" // ::all_tests
#include "boundary_condition.hxx" // ::all_tests
#include "radial_integrator.hxx" // ::all_tests
#include "geometry_analysis.hxx" // ::all_tests
#include "density_generator.hxx" // ::all_tests
#include "fourier_transform.hxx" // ::all_tests
#include "iterative_poisson.hxx" // ::all_tests
// #include "sigma_autoconfig.hxx" // ::all_tests
#include "radial_potential.hxx" // ::all_tests
#include "bessel_transform.hxx" // ::all_tests
#include "parallel_domains.hxx" // ::all_tests
#include "scattering_test.hxx" // ::all_tests
#include "davidson_solver.hxx" // ::all_tests
#include "chemical_symbol.hxx" // ::all_tests
#include "linear_operator.hxx" // ::all_tests
#include "sho_hamiltonian.hxx" // ::all_tests
#include "fourier_poisson.hxx" // ::all_tests
#include "solid_harmonics.hxx" // ::all_tests
#include "bisection_tools.hxx" // ::all_tests
#include "pw_hamiltonian.hxx" // ::all_tests
#include "sho_projection.hxx" // ::all_tests
#include "shift_boundary.hxx" // ::all_tests
#include "linear_algebra.hxx" // ::all_tests
#include "grid_operators.hxx" // ::all_tests
#include "dense_operator.hxx" // ::all_tests
#include "element_config.hxx" // ::all_tests
#include "complex_tools.hxx" // ::all_tests
#include "vector_layout.hxx" // ::all_tests
#include "sho_potential.hxx" // ::all_tests
// #include "isolated_atom.hxx" // ::all_tests
#include "angular_grid.hxx" // ::all_tests
// #include "pseudo_tools.hxx" // ::all_tests
#include "inline_tools.hxx" // ::all_tests
#include "simple_timer.hxx" // ::all_tests
#include "sigma_config.hxx" // ::all_tests
#include "unit_system.hxx" // ::all_tests
#include "simple_math.hxx" // ::all_tests
#include "sho_overlap.hxx" // ::all_tests
#include "radial_grid.hxx" // ::all_tests
#include "single_atom.hxx" // ::all_tests
#include "inline_math.hxx" // ::all_tests
#include "sho_unitary.hxx" // ::all_tests
#include "atom_image.hxx" // ::all_tests
#include "real_space.hxx" // ::all_tests
#include "multi_grid.hxx" // ::all_tests
#include "sho_radial.hxx" // ::all_tests
#include "sho_tools.hxx" // ::all_tests
#include "atom_core.hxx" // ::all_tests
#include "data_view.hxx" // ::all_tests
#include "control.hxx" // ::all_tests
#endif


#ifndef _Output_Units_Fixed
      #include "display_units.h" // extern definitions
      // global variables
      double eV  = 1; char const *_eV  = ""; // dynamic energy unit
      double Ang = 1; char const *_Ang = ""; // dynamic length unit
      // end global variables
#endif

  status_t run_unit_tests(char const *module=nullptr, int const echo=0) {
      status_t status(0);
#ifdef  NO_UNIT_TESTS
      error("version was compiled with -D NO_UNIT_TESTS");
      status = -1;
#else
      bool const all = (nullptr == module);
      auto const m = std::string(all ? "" : module);
      if (all) { if (echo > 0) printf("\n# run all tests!\n\n"); } 
      else     { if (echo > 0) printf("# run unit tests for module '%s'\n\n", m.c_str()); }

      std::vector<std::pair<std::string,status_t>> run;
      { // testing scope
#define   module_test(NAME, FUN) \
          if (all || (0 == std::string(NAME).compare(m))) { \
              SimpleTimer timer("module test for", 0, NAME, 0); \
              run.push_back(make_pair(std::string(NAME), FUN(echo))); \
          }
//            if (echo > -1) printf("\n# Module test for %s\n\n", NAME);
          module_test("recorded_warnings.",     recorded_warnings::all_tests);
          module_test("finite_difference.",     finite_difference::all_tests);
          module_test("hermite_polynomial.",   hermite_polynomial::all_tests);
          module_test("spherical_harmonics.", spherical_harmonics::all_tests);
          module_test("conjugate_gradients.", conjugate_gradients::all_tests);
          module_test("potential_generator.", potential_generator::all_tests);
          module_test("radial_eigensolver.",   radial_eigensolver::all_tests);
          module_test("boundary_condition.",   boundary_condition::all_tests);
          module_test("fermi_distribution.",   fermi_distribution::all_tests);
          module_test("radial_integrator.",     radial_integrator::all_tests);
          module_test("geometry_analysis.",     geometry_analysis::all_tests);
          module_test("density_generator.",     density_generator::all_tests);
          module_test("fourier_transform.",     fourier_transform::all_tests);
          module_test("iterative_poisson.",     iterative_poisson::all_tests);
          module_test("radial_potential.",       radial_potential::all_tests);
//           module_test("sigma_autoconfig.",       sigma_autoconfig::all_tests);
          module_test("bessel_transform.",       bessel_transform::all_tests);
          module_test("parallel_domains.",       parallel_domains::all_tests);
          module_test("scattering_test.",         scattering_test::all_tests);
          module_test("davidson_solver.",         davidson_solver::all_tests);
          module_test("chemical_symbol.",         chemical_symbol::all_tests);
          module_test("linear_operator.",         linear_operator::all_tests);
          module_test("sho_hamiltonian.",         sho_hamiltonian::all_tests);
          module_test("fourier_poisson.",         fourier_poisson::all_tests);
          module_test("solid_harmonics.",         solid_harmonics::all_tests);
          module_test("bisection_tools.",         bisection_tools::all_tests);
          module_test("pw_hamiltonian.",           pw_hamiltonian::all_tests);
          module_test("sho_projection.",           sho_projection::all_tests);
          module_test("shift_boundary.",           shift_boundary::all_tests);
          module_test("linear_algebra.",           linear_algebra::all_tests);
          module_test("grid_operators.",           grid_operators::all_tests);
          module_test("dense_operator.",           dense_operator::all_tests);
          module_test("element_config.",           element_config::all_tests);
          module_test("complex_tools.",             complex_tools::all_tests);
          module_test("vector_layout.",             vector_layout::all_tests);
          module_test("sho_potential.",             sho_potential::all_tests);
//           module_test("isolated_atom.",             isolated_atom::all_tests);
          module_test("angular_grid.",               angular_grid::all_tests);
          module_test("inline_tools.",               inline_tools::all_tests);
//           module_test("pseudo_tools.",               pseudo_tools::all_tests);
          module_test("simple_timer.",               simple_timer::all_tests);
          module_test("sigma_config.",               sigma_config::all_tests);
          module_test("unit_system.",                 unit_system::all_tests);
          module_test("simple_math.",                 simple_math::all_tests);
          module_test("sho_overlap.",                 sho_overlap::all_tests);
          module_test("radial_grid.",                 radial_grid::all_tests);
          module_test("single_atom.",                 single_atom::all_tests);
          module_test("inline_math.",                 inline_math::all_tests);
          module_test("sho_unitary.",                 sho_unitary::all_tests);
          module_test("atom_image.",                   atom_image::all_tests);
          module_test("real_space.",                   real_space::all_tests);
          module_test("multi_grid.",                   multi_grid::all_tests);
          module_test("sho_radial.",                   sho_radial::all_tests);
          module_test("sho_tools.",                     sho_tools::all_tests);
          module_test("atom_core.",                     atom_core::all_tests);
          module_test("data_view.",                     data_view::all_tests);
          module_test("control.",                         control::all_tests);
#undef    module_test
      } // testing scope

      if (run.size() < 1) { // nothing has been tested
          if (echo > 0) printf("# ERROR: test for '%s' not found!\n", module);
          status = -1;
      } else {
          if (echo > 0) printf("\n\n#%3ld modules have been tested:\n", run.size());
          for(auto r : run) {
              auto const stat = r.second;
              if (echo > 0) printf("#    module= %-24s status= %i\n", r.first.c_str(), int(stat));
              status += std::abs(int(stat));
          } // r
          if (echo > 0) {
              printf("\n#%3ld modules have been tested,  total status= %d\n\n", run.size(), int(status));
              if (status > 0) printf("# Warning! At least one module test failed!\n");
          } // echo
      } // something has been tested
#endif
      return status;
  } // run_unit_tests

  int show_help(char const *executable) {
      printf("Usage %s [OPTION]\n"
        "   --help           [-h]\tThis help message\n"
        "   --version            \tShow version number\n"
#ifndef  NO_UNIT_TESTS
        "   --test <module.> [-t]\tTest module\n"
#endif
        "   --verbose        [-v]\tIncrement verbosity level\n"
        "   +<name>=<value>      \tModify variable environment\n"
        "\n", executable);
      return 0;
  } // show_help

  int show_version(char const *executable="#") {
#ifdef _GIT_KEY
      // stringify the value of a macro, two expansion levels needed
      #define macro2string(a) stringify(a)
      #define stringify(b) #b
      printf("%s git checkout " macro2string(_GIT_KEY) "\n\n", executable);
      #undef  stringify
      #undef  macro2string
#endif
      return 0;
  } // show_version
  
  int main(int const argc, char const *argv[]) {
      status_t stat(0);
      char const *test_unit = nullptr; // the name of the unit to be tested
      int run_tests{0};
      int verbosity{3}; // set default verbosity low
      if (argc < 2) {
          printf("%s: no arguments passed!\n", (argc < 1)?__FILE__:argv[0]); 
          return -1;
      } // no argument passed to executable
      for(int iarg = 1; iarg < argc; ++iarg) {
          assert(nullptr != argv[iarg]);
          char const ci0 = *argv[iarg]; // char #0 of command line argument #1
          if ('-' == ci0) {

              // options (short or long)
              char const ci1 = *(argv[iarg] + 1); // char #1 of command line argument #1
              char const IgnoreCase = 32; // use with | to convert upper case chars into lower case chars
              if ('-' == ci1) {

                  // long options
                  std::string option(argv[iarg] + 2); // remove two '-' in front
                  if ("help" == option) {
                      return show_help(argv[0]);
                  } else 
                  if ("version" == option) {
                      return show_version(argv[0]);
                  } else 
                  if ("verbose" == option) {
                      verbosity = 6; // set high
                  } else
                  if ("test" == option) {
                      ++run_tests; if (iarg + 1 < argc) test_unit = argv[iarg + 1];
                  } else {
                      ++stat; warn("# ignored unknown command line option --%s", option.c_str());
                  } // option

              } else { // ci1

                  // short options
                  if ('h' == (ci1 | IgnoreCase)) {
                      return show_help(argv[0]);
                  } else
                  if ('v' == (ci1 | IgnoreCase)) {
                      ++verbosity; verbosity += 3*('V' == ci1); // increment by 'V':4, 'v':1
                  } else
                  if ('t' == (ci1 | IgnoreCase)) {
                      ++run_tests; if (iarg + 1 < argc) test_unit = argv[iarg + 1];
                  } else {
                      ++stat; warn("# ignored unknown command line option -%c", ci1);
                  } // ci1

              } // ci1

          } else // ci0
          if ('+' == ci0) {
              stat += control::command_line_interface(argv[iarg] + 1); // start after the '+' char
          } else
          if (argv[iarg] != test_unit) {
              ++stat; warn("# ignored command line argument \'%s\'", argv[iarg]);
          } // ci0

      } // iarg
      int echo{verbosity}; // define verbosity for repeating arguments and control file entries
      if (echo > 0) {
          printf("\n#");
          for(int iarg = 0; iarg < argc; ++iarg) {
              printf(" %s", argv[iarg]); // repeat the command line arguments
          }   printf("\n");
      } // echo
      //
      // in addition to command_line_interface, we can modify the control environment by a file
      stat += control::read_control_file(control::get("control.file", ""), echo);
      //
      echo = int(control::get("verbosity", double(verbosity))); // redefine verbosity here
      //
      if (echo > 0) {
          printf("\n");
          show_version();
          printf("\n# verbosity = %d\n", echo);
      } // echo
      stat += unit_system::set_output_units(
                  control::get("output.energy.unit", "Ha"),
                  control::get("output.length.unit", "Bohr"));
      if (run_tests) stat += run_unit_tests(test_unit, echo);
      if (echo > 0) recorded_warnings::show_warnings(3);
      recorded_warnings::clear_warnings(1);
      return int(stat);
  } // main
