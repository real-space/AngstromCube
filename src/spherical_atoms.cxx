#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::floor

#include "spherical_atoms.hxx"

#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set
#include "constants.hxx" // pi
#include "real_space_grid.hxx" // grid_t
#include "chemical_symbol.h" // element_symbols
// #include "single_atom.hxx"

// #define FULL_DEBUG
// #define DEBUG

namespace spherical_atoms {
  // this module makes a DFT calculation based on atoms
  // that live in a spherical potential which is found
  // by projecting the 3D potential.
  // their wave functions do not hybridize but they 
  // feel the effect of the density of neighboring atoms
  
  status_t init(int const echo=9) {
      int const dims[] = {96, 80, 80};
      real_space_grid::grid_t<double,1> g(dims);
      g.set_grid_spacing(0.25);
      set(g.values, g.all(), 0.0);
      int const na = 2;
      double const xyzZ[na][4] = {{-2,0,0, 13}, {2,0,0, 15}}; // Al-P
      printf("# %s List of Atoms: (coordinates in %s)\n", __func__,_Ang);
      for(int ia = 0; ia < na; ++ia) {
          int const iZ = (int)std::round(xyzZ[ia][3]);
          printf("%c%c  %16.9f%16.9f%16.9f\n", element_symbols[2*iZ], element_symbols[2*iZ + 1], 
                  xyzZ[ia][0]*Ang, xyzZ[ia][1]*Ang, xyzZ[ia][2]*Ang);
      } // ia
      
      // compute the self-consistent solution of a single_atom, all states in the core
      // get the spherical core_density and bring it to the 3D grid
      // get the ell=0 compensator charge and add it to the 3D grid
      // envoke exchange_correlation and fourier_poisson
      // add XC and electrostatic potential and zero_potential contributions
      // project the total effective potential to each center using bessel_transforms
      // feed back spherical potential into single_atom
      
      return 0;
  } // init

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo=9) {
      return init();
  } // test_create_and_destroy

  status_t all_tests() {
    auto status = 0;
    status += test_create_and_destroy();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace spherical_atoms
