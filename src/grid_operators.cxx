#include <cstdio> // printf
#include <vector> // std::vector<T>
#include <cassert> // assert
 
#include "grid_operators.hxx"

#include "real_space.hxx" // ::grid_t
#include "finite_difference.hxx" // ::stencil_t, ::derive
#include "atom_image.hxx" // ::atom_image_t, ::sho_atom_t
#include "sho_projection.hxx" // ::sho_prefactor

namespace grid_operators {
  // setup of the real-space grid-based Hamiltonian and overlap operator
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t basic_test(int const echo=9) {
      status_t stat(0);
      int const dims[] = {12, 13, 14};
      real_space::grid_t g(dims);
      std::vector<double> psi(2*g.all(), 1.0);
      std::vector<double> potential(dims[2]*dims[1]*dims[0], 0.5);
      std::vector<atom_image::sho_atom_t> a(1);
      a[0] = atom_image::sho_atom_t(0.5, 3, 999); // sigma=0.5, numax=3, atom_id=999, no position given!
      double const apos[] = {0,0,0};
      a[0].set_image_positions(apos);
      finite_difference::stencil_t<double> kinetic(g.h, 8, -0.5);
//       // decativated: interfaces deprecated, use op.Hamiltonian in the future!
//       stat += grid_Hamiltonian(psi.data(), &psi[g.all()], g, a, kinetic, potential.data());
//       stat += grid_Overlapping(psi.data(), &psi[g.all()], g, a);
      return stat;
  } // basic_test

  
  status_t projector_normalization_test(int const echo=9) {
      status_t stat(0);
      real_space::grid_t g(32, 32, 32); // grid spacings = {1,1,1} by default
      int const numax = 3;
      double const sigma = 2.0;
      std::vector<atom_image::sho_atom_t> a;
      //                          x   y   z   Z   id   numax  sigma
      double const xyzZinso[] = {16, 16, 16, 7.5, 747, numax, sigma};
      stat += list_of_atoms(a, xyzZinso, 1, 7, g, echo);
      int const ncoeff = sho_tools::nSHO(numax);
      std::vector<double> psi(g.all(), 0.0), a_vec(ncoeff, 0.0), p_vec(ncoeff);
      auto const a_coeff = a_vec.data();
      auto       p_coeff = p_vec.data();
      grid_operator_t<double> op(g, a);
      double const f2 = pow2(sho_projection::sho_prefactor(0, 0, 0, sigma));
      double dev{0};
      for(int ic = 0; ic < ncoeff; ++ic) {
          set(psi.data(), g.all(), 0.0);
          a_coeff[ic] = 1;
          stat += op.get_start_waves(psi.data(), &a_coeff, echo);
          stat += op.get_atom_coeffs(&p_coeff, psi.data(), echo);
          if (echo > 5) printf("# %s%3i  ", __func__, ic);
          for(int jc = 0; jc < ncoeff; ++jc) {
              if (ic == jc) { 
                  if (echo > 5) printf(" %7.3f", p_coeff[jc]*f2);
              } else {
                  dev = std::max(dev, std::abs(p_coeff[jc]*f2));
              }
          } // jc
          if (echo > 5) printf("\n");
          a_coeff[ic] = 0;
      } // ic
      if (echo > 2) printf("# %s: largest deviation is %.1e\n", __func__, dev);
      return stat + (dev > 3e-14);
  } // projector_normalization_test
  
  status_t class_test(int const echo=9) {
      status_t stat(0);
      int const dims[] = {12, 13, 14};
      int const all = dims[0]*dims[1]*dims[2];
      std::vector<double> psi(all, 1.0), Hpsi(all);
      grid_operator_t<double> op(dims);
      stat += op.Hamiltonian(Hpsi.data(), psi.data(), echo);
      stat += op.Overlapping(psi.data(), Hpsi.data(), echo);
      stat += op.Conditioner(Hpsi.data(), psi.data(), echo);
      return stat;
  } // class_test

  status_t class_with_atoms_test(int const echo=9) {
      status_t stat(0);
      real_space::grid_t g(36, 25, 24);
      std::vector<atom_image::sho_atom_t> a;
      //                          x   y   z    Z     id nu sigma dummy
      double const xyzZinso[] = {.1, .2, -4,  13.0,  767, 3, 1.5, 9e9,
                                -.1, .2,  3,  15.1,  757, 4, 1.7, 8e8};
      stat += list_of_atoms(a, xyzZinso, 2, 8, g, echo);
      std::vector<double> psi(g.all(), 1.0), Hpsi(g.all());
      grid_operator_t<double> op(g, a);
      stat += op.Hamiltonian(Hpsi.data(), psi.data(), echo);
      stat += op.Overlapping(psi.data(), Hpsi.data(), echo);
      stat += op.Conditioner(Hpsi.data(), psi.data(), echo);
      return stat;
  } // class_with_atoms_test
  
  status_t all_tests(int const echo) {
    status_t status(0);
    status += class_test(echo);
    status += basic_test(echo);
    status += class_with_atoms_test(echo);
    status += projector_normalization_test(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace grid_operators
