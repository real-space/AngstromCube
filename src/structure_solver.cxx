// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, ::snprintf
#include <vector> // std::vector<T>

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

#include "structure_solver.hxx" // RealSpaceKohnSham

#include "real_space.hxx" // ::grid_t
#include "fermi_distribution.hxx" // ::FermiLevel_t, ::Fermi_level
#include "data_list.hxx" // data_list<T>
#include "data_view.hxx" // view2D<T>

namespace structure_solver {

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo) {
      status_t stat(0);
      if (echo > 5) std::printf("\n# %s start\n\n", __func__);
      real_space::grid_t g(7, 8, 9);
      int const natoms = 0;
      view2D<double> xyzZinso(natoms, 8, 0.0);
      RealSpaceKohnSham KS(g, xyzZinso, natoms, 0, echo);
      if (echo > 5) std::printf("\n# %s done\n\n", __func__);
      return stat;
  } // test_create_and_destroy

  status_t test_free_electrons(int const echo) {
      status_t stat(0);
      if (echo > 5) std::printf("\n# %s start\n\n", __func__);
      real_space::grid_t g(16, 16, 16); // so far boundary condition is isolated
      g.set_boundary_conditions(1, 1, 1);
      int const natoms = 0;
      view2D<double> xyzZinso(natoms, 8, 0.0);
      assert(control::get("bands.extra", 1.) > 0); // overwrite the default for "bands.extra" which is 0 in the RealSpaceKohnSham constructor.
      RealSpaceKohnSham KS(g, xyzZinso, natoms, 1, echo);
      double charges[] = {0, 0, 0, 0};
      view2D<double> valence_rho(2, g.all(), 0.0);
      std::vector<int32_t> n_atom_rho(natoms, 1); // list of pow2(nSHO(numax[ia]))
      data_list<double> atom_rho[2];
      atom_rho[0] = data_list<double>(n_atom_rho); // density
      atom_rho[1] = data_list<double>(n_atom_rho); // response
      std::vector<int32_t> n_atom_mat(natoms, 2); // list of 2*pow2(nSHO(numax[ia]))
      data_list<double> atom_mat(n_atom_mat, 0.0);
      fermi_distribution::FermiLevel_t Fermi(control::get("valence.electrons", 0.5));
      std::vector<double> Vtot(g.all(), 0.0); // empty potential
      stat += KS.solve(valence_rho, atom_rho, charges, Fermi, g, Vtot.data(), natoms, atom_mat, 'e', -1, echo);
      if (echo > 5) std::printf("\n# %s done\n\n", __func__);
      return stat;
  } // test_free_electrons

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_free_electrons(echo);
      stat += test_create_and_destroy(echo); // this needs to run seconds as it asks for control::get("bands.extra", 0.);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace structure_solver
