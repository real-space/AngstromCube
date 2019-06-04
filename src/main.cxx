#include <cstdio> // printf
#include <vector> // std::vector
#include <string> // std::string
#include <cassert> // assert

#include "element_configuration.hxx" // all_tests
#include "radial_eigensolver.hxx" // all_tests
#include "radial_integrator.hxx" // all_tests
#include "radial_potential.hxx" // all_tests
#include "angular_grid.hxx" // all_tests
#include "radial_grid.hxx" // all_tests
#include "single_atom.hxx" // all_tests
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
          module_test("radial_eigensolver.",       radial_eigensolver::all_tests);
          module_test("radial_integrator.",         radial_integrator::all_tests);
          module_test("radial_potential.",           radial_potential::all_tests);
          module_test("angular_grid.",                   angular_grid::all_tests);
          module_test("radial_grid.",                     radial_grid::all_tests);
          module_test("single_atom.",                     single_atom::all_tests);
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
