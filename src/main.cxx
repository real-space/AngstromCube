#include <cstdio> // printf
#include <vector> // std::vector
#include <string> // std::string
#include <cassert> // assert

#include "element_configuration.hxx" // all_tests
#include "hermite_polynomials.hxx" // all_tests
#include "radial_eigensolver.hxx" // all_tests
#include "radial_integrator.hxx" // all_tests
#include "finite_difference.hxx" // all_tests
#include "radial_potential.hxx" // all_tests
#include "bessel_transform.hxx" // all_tests
#include "real_space_grid.hxx" // all_tests
#include "fourier_poisson.hxx" // all_tests
#include "sho_projection.hxx" // all_tests
#include "angular_grid.hxx" // all_tests
#include "radial_grid.hxx" // all_tests
#include "single_atom.hxx" // all_tests
#include "inline_math.hxx" // all_tests
#include "sho_unitary.hxx" // all_tests
#include "sho_radial.hxx" // all_tests
#include "sho_tools.hxx" // all_tests
#include "atom_core.hxx" // all_tests
#include "overlap.hxx" // all_tests

  int run_unit_tests(char const *module) 
  {
      bool const all = (nullptr == module);
      auto const m = std::string(all?"":module);
      if (all) { printf("# run all tests!\n"); } 
      else     { printf("# run unit tests for module '%s'\n\n", m.c_str()); }
      
      std::vector<std::pair<std::string,int>> run;
      { // testing scope
#define   module_test(NAME, FUN) if (all || (0 == std::string(NAME).compare(m))) \
                           { run.push_back(make_pair(std::string(NAME), FUN())); }
          module_test("element_configuration.", element_configuration::all_tests);
          module_test("hermite_polynomials.",     hermite_polynomials::all_tests);
          module_test("radial_eigensolver.",       radial_eigensolver::all_tests);
          module_test("finite_difference.",         finite_difference::all_tests);
          module_test("radial_integrator.",         radial_integrator::all_tests);
          module_test("radial_potential.",           radial_potential::all_tests);
          module_test("bessel_transform.",           bessel_transform::all_tests);
          module_test("real_space_grid.",             real_space_grid::all_tests);
          module_test("fourier_poisson.",             fourier_poisson::all_tests);
          module_test("sho_projection.",               sho_projection::all_tests);
          module_test("angular_grid.",                   angular_grid::all_tests);
          module_test("radial_grid.",                     radial_grid::all_tests);
          module_test("single_atom.",                     single_atom::all_tests);
          module_test("inline_math.",                     inline_math::all_tests);
          module_test("sho_unitary.",                     sho_unitary::all_tests);
          module_test("sho_radial.",                       sho_radial::all_tests);
          module_test("sho_tools.",                         sho_tools::all_tests);
          module_test("atom_core.",                         atom_core::all_tests);
          module_test("overlap.",                             overlap::all_tests);
#undef    module_test
      } // testing scope

      int status = 0;
      if (run.size() < 1) { // nothing has been tested
          printf("# ERROR: test for '%s' not found!\n", module);
          status = -1;
      } else {
          printf("\n# %ld modules have been tested:\n", run.size());
          for(auto r : run) {
              auto const stat = r.second;
              printf("#   module = '%s' \t status = %d\n", r.first.c_str(), stat);
              status += abs(stat);
          } // r
          printf("# total status = %d\n\n", status);
          if (status > 0) printf("# Warning! At least one module test failed!\n");
      } // something has been tested
      return status;
  } // run_unit_tests


  int main(int const argc, char const * argv[]) 
  {
      if (argc < 2) { printf("%s: no arguments passed!\n", (argc < 1)?__FILE__:argv[0]); return -1; }
//    printf("%s: argument #1 is %s\n", argv[0], argv[1]);
      char const c10 = *argv[1]; // char #0 of command line argument #1
      if ('-' == c10) {
          char const c11 = *(argv[1] + 1); // char #1 of command line argument #1
          char const IgnoreCase = 32; // use with | to convert upper case chars into lower case chars
          if ('h' == (c11 | IgnoreCase)) {
              printf("Usage %s [OPTION]\n", argv[0]);
              printf(" -h, -H     \tThis Help message\n"
                      " -t <module> \tTest module\n");
          } else if ('t' == (c11 | IgnoreCase)) {
              return run_unit_tests((argc > 2)?argv[2]:nullptr);
          } // help or test
      } // option expected
      return 0;
  } // main
