#include <cstdio> // printf
#include <vector> // std::vector
#include <string> // std::string
#include <cassert> // assert

typedef int status_t;

#include "recorded_warnings.hxx" // all_tests
#include "hermite_polynomial.hxx" // all_tests
#include "radial_eigensolver.hxx" // all_tests
#include "boundary_condition.hxx" // all_tests
#include "radial_integrator.hxx" // all_tests
#include "finite_difference.hxx" // all_tests
#include "geometry_analysis.hxx" // all_tests
#include "radial_potential.hxx" // all_tests
#include "bessel_transform.hxx" // all_tests
#include "real_space_grid.hxx" // all_tests
#include "sho_hamiltonian.hxx" // all_tests
#include "chemical_symbol.hxx" // all_tests
#include "fourier_poisson.hxx" // all_tests
#include "spherical_atoms.hxx" // all_tests
#include "solid_harmonics.hxx" // all_tests
#include "sho_projection.hxx" // all_tests
#include "grid_operators.hxx" // all_tests
#include "element_config.hxx" // all_tests
#include "angular_grid.hxx" // all_tests
#include "radial_grid.hxx" // all_tests
#include "single_atom.hxx" // all_tests
#include "inline_math.hxx" // all_tests
#include "sho_unitary.hxx" // all_tests
#include "sho_radial.hxx" // all_tests
#include "sho_tools.hxx" // all_tests
#include "atom_core.hxx" // all_tests
#include "overlap.hxx" // all_tests
#include "spherical_harmonics.hxx" // no test implemented

#include "recorded_warnings.hxx" // show_warnings, clear_warnings

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
          module_test("recorded_warnings.",   recorded_warnings::all_tests);
          module_test("hermite_polynomial.", hermite_polynomial::all_tests);
          module_test("radial_eigensolver.", radial_eigensolver::all_tests);
          module_test("boundary_condition.", boundary_condition::all_tests);
          module_test("finite_difference.",   finite_difference::all_tests);
          module_test("geometry_analysis.",   geometry_analysis::all_tests);
          module_test("radial_integrator.",   radial_integrator::all_tests);
          module_test("radial_potential.",     radial_potential::all_tests);
          module_test("bessel_transform.",     bessel_transform::all_tests);
          module_test("real_space_grid.",       real_space_grid::all_tests);
          module_test("sho_hamiltonian.",       sho_hamiltonian::all_tests);
          module_test("chemical_symbol.",       chemical_symbol::all_tests);
          module_test("fourier_poisson.",       fourier_poisson::all_tests);
          module_test("spherical_atoms.",       spherical_atoms::all_tests);
          module_test("solid_harmonics.",       solid_harmonics::all_tests);
          module_test("sho_projection.",         sho_projection::all_tests);
          module_test("grid_operators.",         grid_operators::all_tests);
          module_test("element_config.",         element_config::all_tests);
          module_test("angular_grid.",             angular_grid::all_tests);
          module_test("radial_grid.",               radial_grid::all_tests);
          module_test("single_atom.",               single_atom::all_tests);
          module_test("inline_math.",               inline_math::all_tests);
          module_test("sho_unitary.",               sho_unitary::all_tests);
          module_test("sho_radial.",                 sho_radial::all_tests);
          module_test("sho_tools.",                   sho_tools::all_tests);
          module_test("atom_core.",                   atom_core::all_tests);
          module_test("overlap.",                       overlap::all_tests);
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
      int stat = 0;
      if (argc < 2) { printf("%s: no arguments passed!\n", (argc < 1)?__FILE__:argv[0]); return -1; }
//    printf("%s: argument #1 is %s\n", argv[0], argv[1]);
      char const c10 = *argv[1]; // char #0 of command line argument #1
      if ('-' == c10) {
          char const c11 = *(argv[1] + 1); // char #1 of command line argument #1
          char const IgnoreCase = 32; // use with | to convert upper case chars into lower case chars
          if ('h' == (c11 | IgnoreCase)) {
              printf("Usage %s [OPTION]\n", argv[0]);
              printf(" -h, -H      \tThis Help message\n"
                     " -t <module> \tTest module\n");
          } else if ('t' == (c11 | IgnoreCase)) {
              stat = run_unit_tests((argc > 2)?argv[2]:nullptr);
          } // help or test
      } // option expected
      recorded_warnings::show_warnings(3);
      recorded_warnings::clear_warnings(1);
      return stat;
  } // main
