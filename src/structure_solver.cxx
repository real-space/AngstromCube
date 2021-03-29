#include <cstdio> // std::printf, std::snprintf
#include <vector> // std::vector<T>

#include "structure_solver.hxx" // RealSpaceKohnSham

#include "real_space.hxx" // ::grid_t
#include "fermi_distribution.hxx" // ::FermiLevel_t, ::Fermi_level
#include "data_list.hxx" // data_list<T>

namespace structure_solver {

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo) {
      status_t stat(0);
//       real_space::grid_t g(7, 8, 9);
//       RealSpaceKohnSham KS(g, list_of_atoms, 1, echo);
      return stat;
  } // test_create_and_destroy

  status_t test_free_electrons(int const echo) {
      status_t stat(0);
//       real_space::grid_t g(8, 8, 8);
//       std::vector<atom_image::sho_atom_t> list_of_atoms(1);
//       list_of_atoms[0] = atom_image::sho_atom_t(1., 1, -1, nullptr, 0);
//       double const pos[3] = {0, 0, 0};
//       list_of_atoms[0].set_image_positions(pos);
//       RealSpaceKohnSham KS(g, list_of_atoms, 1, echo);
//       std::vector<int32_t> n_atom_rho(1, 4);
//       data_list<double> atom_mat(n_atom_rho);
//       fermi_distribution::FermiLevel_t Fermi;
//       std::vector<double> Vtot(g.all(), 0.0);
//       KS.solve(Fermi, g, Vtot.data(), n_atom_rho, atom_mat);
      return stat;
  } // test_free_electrons

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_create_and_destroy(echo);
      stat += test_free_electrons(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace structure_solver
