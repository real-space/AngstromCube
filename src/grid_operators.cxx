#include <cstdio> // printf
#include <vector> // std::vector<T>
#include <cassert> // assert
 
#include "grid_operators.hxx"

#include "real_space_grid.hxx" // ::grid_t
#include "finite_difference.hxx" // ::finite_difference_t, ::Laplacian
#include "atom_image.hxx" // ::atom_image_t, ::sho_atom_t
#include "inline_math.hxx" // set
#include "sho_projection.hxx" // ::sho_project, ::sho_add
#include "sho_tools.hxx" // ::nSHO

// #include "display_units.h" // eV, _eV, Ang, _Ang

// #define FULL_DEBUG
#define DEBUG

#ifdef  DEBUG
    #include "debug_output.hxx" // dump_to_file
#endif

#ifdef FULL_DEBUG
    #define full_debug(print) print 
#else
    #define full_debug(print)
#endif

#ifdef DEBUG
    #define debug(print) print 
#else
    #define debug(print)
#endif


namespace grid_operators {
  // setup of the real-space grid-based Hamiltonian and overlap operator
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t basic_test(int const echo=9) {
      status_t stat = 0;
      int constexpr D0 = 2; // vectorization
      int const dims[] = {12, 13, 14}, nn[] = {8, 8, 8};
      real_space_grid::grid_t<D0> g(dims);
      std::vector<double> psi(2*g.all(), 1.0);
      std::vector<double> potential(dims[2]*dims[1]*dims[0], 0.5);
      std::vector<atom_image::sho_atom_t> a(1);
      std::vector<atom_image::atom_image_t> ai(1);
      a[0]  = atom_image::sho_atom_t(3, 0.5, 999); // numax=3, sigma=0.5, atom_id=999
      ai[0] = atom_image::atom_image_t(dims[0]*g.h[0]/2, dims[1]*g.h[1]/2, dims[2]*g.h[2]/2, 999, 0);
      finite_difference::finite_difference_t<double> kinetic(g.h, nn);
      kinetic.scale_coefficients(-0.5);
      stat += grid_Hamiltonian(psi.data(), &psi[g.all()], g, a, ai, kinetic, potential.data());
      stat += grid_Overlapping(psi.data(), &psi[g.all()], g, a, ai);
      return stat;
  } // basic_test

  status_t class_test(int const echo=9) {
      status_t stat = 0;
      int constexpr D0 = 1; // vectorization
      int const dims[] = {12, 13, 14};
      int const all = dims[0]*dims[1]*dims[2];
      std::vector<double> psi(all, 1.0), Hpsi(all);
      grid_operator_t<double, double, D0> op(dims);
      stat += op.Hamiltonian(Hpsi.data(), psi.data(), echo);
      stat += op.Overlapping(psi.data(), Hpsi.data(), echo);
      stat += op.Conditioner(Hpsi.data(), psi.data(), echo);
      return stat;
  } // class_test
  
  status_t all_tests(int const echo) {
    auto status = 0;
    status += class_test(echo);
    status += basic_test(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace grid_operators
