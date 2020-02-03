#include <cstdio> // printf
#include <cassert> // assert
#include <cstdlib> // std::abs
#include <vector> // std::vector
#include <string> // std::string
#include <utility> // std::pair

#include "recorded_warnings.hxx" // warn, ::show_warnings, ::clear_warnings
#include "simple_timer.hxx" // SimpleTimer
#include "control.hxx" // ::cli

typedef int status_t;

#include "recorded_warnings.hxx" // ::all_tests
#include "finite_difference.hxx" // ::all_tests
#include "hermite_polynomial.hxx" // ::all_tests
#include "spherical_harmonics.hxx" // ::all_tests
#include "conjugate_gradients.hxx" // ::all_tests
#include "radial_eigensolver.hxx" // ::all_tests
#include "boundary_condition.hxx" // ::all_tests
#include "radial_integrator.hxx" // ::all_tests
#include "geometry_analysis.hxx" // ::all_tests
#include "radial_potential.hxx" // ::all_tests
#include "bessel_transform.hxx" // ::all_tests
#include "scattering_test.hxx" // ::all_tests
#include "real_space_grid.hxx" // ::all_tests
#include "davidson_solver.hxx" // ::all_tests
#include "chemical_symbol.hxx" // ::all_tests
#include "linear_operator.hxx" // ::all_tests
#include "fourier_poisson.hxx" // ::all_tests
#include "spherical_atoms.hxx" // ::all_tests
#include "solid_harmonics.hxx" // ::all_tests
#include "sho_projection.hxx" // ::all_tests
#include "shift_boundary.hxx" // ::all_tests
#include "linear_algebra.hxx" // ::all_tests
#include "grid_operators.hxx" // ::all_tests
#include "element_config.hxx" // ::all_tests
#include "vector_layout.hxx" // ::all_tests
#include "sho_potential.hxx" // ::all_tests
#include "angular_grid.hxx" // ::all_tests
#include "simple_timer.hxx" // ::all_tests
#include "simple_math.hxx" // ::all_tests
#include "sho_overlap.hxx" // ::all_tests
#include "radial_grid.hxx" // ::all_tests
#include "single_atom.hxx" // ::all_tests
#include "inline_math.hxx" // ::all_tests
#include "sho_unitary.hxx" // ::all_tests
#include "atom_image.hxx" // ::all_tests
#include "sho_radial.hxx" // ::all_tests
#include "sho_tools.hxx" // ::all_tests
#include "atom_core.hxx" // ::all_tests
#include "data_view.hxx" // ::all_tests
#include "control.hxx" // ::all_tests

  status_t run_unit_tests(char const *module, int const echo=0) 
  {
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
          module_test("radial_eigensolver.",   radial_eigensolver::all_tests);
          module_test("boundary_condition.",   boundary_condition::all_tests);
          module_test("radial_integrator.",     radial_integrator::all_tests);
          module_test("geometry_analysis.",     geometry_analysis::all_tests);
          module_test("radial_potential.",       radial_potential::all_tests);
          module_test("bessel_transform.",       bessel_transform::all_tests);
          module_test("scattering_test.",         scattering_test::all_tests);
          module_test("real_space_grid.",         real_space_grid::all_tests);
          module_test("davidson_solver.",         davidson_solver::all_tests);
          module_test("chemical_symbol.",         chemical_symbol::all_tests);
          module_test("linear_operator.",         linear_operator::all_tests);
          module_test("fourier_poisson.",         fourier_poisson::all_tests);
          module_test("spherical_atoms.",         spherical_atoms::all_tests);
          module_test("solid_harmonics.",         solid_harmonics::all_tests);
          module_test("sho_projection.",           sho_projection::all_tests);
          module_test("shift_boundary.",           shift_boundary::all_tests);
          module_test("linear_algebra.",           linear_algebra::all_tests);
          module_test("grid_operators.",           grid_operators::all_tests);
          module_test("element_config.",           element_config::all_tests);
          module_test("vector_layout.",             vector_layout::all_tests);
          module_test("sho_potential.",             sho_potential::all_tests);
          module_test("angular_grid.",               angular_grid::all_tests);
          module_test("simple_timer.",               simple_timer::all_tests);
          module_test("simple_math.",                 simple_math::all_tests);
          module_test("sho_overlap.",                 sho_overlap::all_tests);
          module_test("radial_grid.",                 radial_grid::all_tests);
          module_test("single_atom.",                 single_atom::all_tests);
          module_test("inline_math.",                 inline_math::all_tests);
          module_test("sho_unitary.",                 sho_unitary::all_tests);
          module_test("atom_image.",                   atom_image::all_tests);
          module_test("sho_radial.",                   sho_radial::all_tests);
          module_test("sho_tools.",                     sho_tools::all_tests);
          module_test("atom_core.",                     atom_core::all_tests);
          module_test("data_view.",                     data_view::all_tests);
          module_test("control.",                         control::all_tests);
#undef    module_test
      } // testing scope

      int status = 0;
      if (run.size() < 1) { // nothing has been tested
          if (echo > 0) printf("# ERROR: test for '%s' not found!\n", module);
          status = -1;
      } else {
          if (echo > 0) printf("\n\n#%3ld modules have been tested:\n", run.size());
          for(auto r : run) {
              auto const stat = r.second;
              if (echo > 0) printf("#    module= %-24s status= %i\n", r.first.c_str(), stat);
              status += std::abs(stat);
          } // r
          if (echo > 0) {
              printf("\n#%3ld modules have been tested,  total status= %d\n\n", run.size(), status);
              if (status > 0) printf("# Warning! At least one module test failed!\n");
          } // echo
      } // something has been tested
      return status;
  } // run_unit_tests


  int main(int const argc, char const *argv[]) {
      status_t stat = 0;
      char const *test_unit = nullptr;
      bool run_tests = false;
      if (argc < 2) { 
          printf("%s: no arguments passed!\n", (argc < 1)?__FILE__:argv[0]); 
          return -1;
      } // no argument passed to executable
      for(int iarg = 1; iarg < argc; ++iarg) {
          char const ci0 = *argv[iarg]; // char #0 of command line argument #1
          if ('-' == ci0) {
              char const ci1 = *(argv[iarg] + 1); // char #1 of command line argument #1
              char const IgnoreCase = 32; // use with | to convert upper case chars into lower case chars
              if ('h' == (ci1 | IgnoreCase)) {
                  printf("Usage %s [OPTION]\n", argv[0]);
                  printf("   -h, -H      \tThis help message\n"
                         "   -t <module> \tTest module\n"
                         "   +<var>=<val>\tModify variable environment\n"
                         "\n");
                  return 0;
              } else if ('t' == (ci1 | IgnoreCase)) {
                  run_tests = true;
                  if (argc > iarg + 1) test_unit = argv[iarg + 1]; // the name of the unit to be tested
              } else {
                  warn("# ignored unknown command line option %c%c", ci0, ci1);
                  ++stat; // error
              } // help or test
          } // '-'
          else
          if ('+' == ci0) {
              stat += control::cli(argv[iarg] + 1); // start after the '+' char
          } // '+'
          else 
          if (argv[iarg] != test_unit) {
              warn("# ignored command line argument %s", argv[iarg]);
          }
      } // iarg
      int const echo = control::get("verbosity", 3.); // define default verbosity here
      if (echo > 0) { 
          printf("\n#");
          for(int iarg = 0; iarg < argc; ++iarg) {
              printf(" %s", argv[iarg]); // repeat the command line arguments
          }   printf("\n\n");
      } // echo
      if (run_tests) stat += run_unit_tests(test_unit, echo);
      if (echo > 0) recorded_warnings::show_warnings(3);
      recorded_warnings::clear_warnings(1);
      return stat;
  } // main
