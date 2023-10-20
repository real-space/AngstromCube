// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <cstdlib> // std::atoi, ::atof
#include <cstdint> // int8_t
#include <string> // std::string, ::string:npos

#include "green_input.hxx"

#include "status.hxx" // status_t

#ifdef  HAS_RAPIDXML
  // git clone https://github.com/dwd/rapidxml or https://github.com/discord/rapidxml
  #include "rapidxml.hpp" // ::xml_document<>
  #include "rapidxml_utils.hpp" // rapidxml::file<>

  #include "xml_reading.hxx" // ::find_attribute, ::find_child
#endif // HAS_RAPIDXML

#include "sho_tools.hxx" // ::nSHO
#include "xml_reading.hxx" // ::read_sequence
#include "inline_math.hxx" // set
#include "control.hxx" // ::get
#include "json_reading.hxx" // ::load_Hamiltonian

namespace green_input {

  status_t load_Hamiltonian(
        uint32_t ng[3] // numbers of grid points
      , int8_t bc[3] // boundary conditions
      , double hg[3] // grid spacings
      , std::vector<double> & Veff
      , int & natoms
      , std::vector<double> & xyzZinso
      , std::vector<std::vector<double>> & atom_mat
      , char const *const filename // ="Hmt.json", input
      , int const echo // =0, log-level
  ) {
      assert(nullptr != filename);
      if (std::string::npos != std::string(filename).find(".json")) {
          if (echo > 0) std::printf("# filename \"%s\" looks like .json formatted\n", filename);
          auto const json_stat = json_reading::load_Hamiltonian(ng, bc, hg, Veff, natoms, xyzZinso, atom_mat, filename, echo);
          if (0 == json_stat) {
              return 0;
          } else {
              if (echo > 0) std::printf("# failed to read file \"%s\" in .json format, try .xml\n", filename);
          }
      } else {
          if (echo > 0) std::printf("# filename \"%s\" looks like .xml formatted\n", filename);
      } // filename contains ".json"

#ifndef HAS_RAPIDXML
      warn("Unable to load_Hamiltonian when compiled without -D HAS_RAPIDXML", 0);
      return STATUS_TEST_NOT_INCLUDED;
#else  // HAS_RAPIDXML

      rapidxml::xml_document<> doc;
      try {
          rapidxml::file<> infile(filename);
          try {
              doc.parse<0>(infile.data());
          } catch (...) {
              warn("failed to parse \"%s\"", filename);
              return -2; // error
          } // try + catch
      } catch (...) {
          warn("failed to open \"%s\"", filename);
          return -1; // error
      } // catch

      set(ng, 3, 0u);
      set(hg, 3, 1.);
      Veff.resize(0);
      natoms = 0;
      xyzZinso.resize(0);
      atom_mat.resize(0);

      auto const grid_Hamiltonian = doc.first_node("grid_Hamiltonian");
      if (grid_Hamiltonian) {
          auto const sho_atoms = xml_reading::find_child(grid_Hamiltonian, "sho_atoms", echo);
          if (sho_atoms) {
              auto const number = xml_reading::find_attribute(sho_atoms, "number", "0", echo);
              if (echo > 5) std::printf("# expect %s sho_atoms\n", number);
              natoms = std::atoi(number);
              xyzZinso.resize(natoms*8);
              atom_mat.resize(natoms);
              int ia{0};
              for (auto atom = sho_atoms->first_node(); atom; atom = atom->next_sibling()) {
                  auto const gid = xml_reading::find_attribute(atom, "gid", "-1");
                  if (echo > 5) std::printf("# <%s gid=%s>\n", atom->name(), gid);
                  xyzZinso[ia*8 + 4] = std::atoi(gid);

                  double pos[3] = {0, 0, 0};
                  auto const position = xml_reading::find_child(atom, "position", echo);
                  for (int d = 0; d < 3; ++d) {
                      char const axyz[] = {char('x' + d), '\0'}; // "x", "y", "z"
                      auto const value = xml_reading::find_attribute(position, axyz);
                      if (*value != '\0') {
                          pos[d] = std::atof(value);
                          if (echo > 5) std::printf("# %s= %.15g\n", axyz, pos[d]);
                      } else warn("no attribute '%c' found in <atom><position> in file \"%s\"", *axyz, filename);
                      xyzZinso[ia*8 + d] = pos[d];
                  } // d

                  auto const projectors = xml_reading::find_child(atom, "projectors", echo);
                  int numax{-1};
                  {
                      auto const value = xml_reading::find_attribute(projectors, "numax", "-1");
                      if (*value != '\0') {
                          numax = std::atoi(value);
                          if (echo > 5) std::printf("# numax= %d\n", numax);
                      } else warn("no attribute 'numax' found in <projectors> in file \"%s\"", filename);
                  }
                  double sigma{-1};
                  {
                      auto const value = xml_reading::find_attribute(projectors, "sigma", "");
                      if (*value != '\0') {
                          sigma = std::atof(value);
                          if (echo > 5) std::printf("# sigma= %g\n", sigma);
                      } else warn("no attribute 'sigma' found in <projectors> in file \"%s\"", filename);
                  }
                  xyzZinso[ia*8 + 5] = numax;
                  xyzZinso[ia*8 + 6] = sigma;
                  int const nSHO = sho_tools::nSHO(numax);
                  atom_mat[ia].resize(2*nSHO*nSHO);
                  for (int h0s1 = 0; h0s1 < 2; ++h0s1) {
                      auto const matrix_name = h0s1 ? "overlap" : "hamiltonian";
                      auto const matrix = xml_reading::find_child(atom, matrix_name, echo);
                      if (matrix) {
                          if (echo > 22) std::printf("# %s.values= %s\n", matrix_name, matrix->value());
                          auto const v = xml_reading::read_sequence<double>(matrix->value(), echo, nSHO*nSHO);
                          if (echo > 5) std::printf("# %s matrix has %ld values, expect %d x %d = %d\n",
                              matrix_name, v.size(), nSHO, nSHO, nSHO*nSHO);
                          assert(v.size() == nSHO*nSHO);
                          double maxdev{0}; // measure deviation from a symmetric matrix
                          for (int i = 0; i < nSHO; ++i) {
                              for (int j = 0; j < i; ++j) {
                                  maxdev = std::max(maxdev, std::abs(v[i*nSHO + j] - v[j*nSHO + i]));
                              } // j
                          } // i
                          if (echo > 3) std::printf("# %s matrix of atom #%s has a max deviation of %.1e from symmetric\n",
                                                      matrix_name, gid, maxdev);
                          if (maxdev > 1e-6) warn("%s matrix of atom #%s has a max deviation of %.1e from symmetric",
                                                      matrix_name, gid, maxdev);
                          set(atom_mat[ia].data() + h0s1*nSHO*nSHO, nSHO*nSHO, v.data()); // copy
                      } else warn("atom with global_id=%s has no %s matrix in file \"%s\"", gid, matrix_name, filename);
                  } // h0s1
                  ++ia; // count up the number of atoms
              } // atom
              assert(natoms == ia && "Inconsistency in file"); // sanity
          } else warn("no <sho_atoms> found in grid_Hamiltonian in file \"%s\"", filename);

          auto const spacing = xml_reading::find_child(grid_Hamiltonian, "spacing", echo);
          for (int d = 0; d < 3; ++d) {
              char axyz[] = {0, 0}; axyz[0] = 'x' + d; // "x", "y", "z"
              auto const value = xml_reading::find_attribute(spacing, axyz);
              if (*value != '\0') {
                  hg[d] = std::atof(value);
                  if (echo > 3) std::printf("# h%s = %.15g\n", axyz, hg[d]);
              } // value != ""
          } // d

          auto const boundary = xml_reading::find_child(grid_Hamiltonian, "boundary", echo);
          for (int d = 0; d < 3; ++d) {
              char axyz[] = {0, 0}; axyz[0] = 'x' + d; // "x", "y", "z"
              auto const value = xml_reading::find_attribute(boundary, axyz);
              if (*value != '\0') {
                  bc[d] = std::atoi(value);
                  if (echo > 3) std::printf("# BC%s = %d\n", axyz, bc[d]);
              } // value != ""
          } // d

          auto const potential = xml_reading::find_child(grid_Hamiltonian, "potential", echo);
          // if (echo > 0) std::printf("# potential pointer %p\n", (void*)potential);
          // ToDo: debug: segfaults on Mac M1 if not about 650 chars are inside the content of potential
          if (potential) {
              for (int d = 0; d < 3; ++d) {
                  char axyz[] = {'n', 0, 0}; axyz[1] = 'x' + d; // "nx", "ny", "nz"
                  auto const value = xml_reading::find_attribute(potential, axyz);
                  if (*value != '\0') {
                      ng[d] = std::atoi(value);
                      if (echo > 5) std::printf("# %s = %d\n", axyz, ng[d]);
                  } // value != ""
              } // d
              auto const ngall = (ng[0])*size_t(ng[1])*size_t(ng[2]);
              if (echo > 33) std::printf("# potential.values= %s\n", potential->value());
              Veff = xml_reading::read_sequence<double>(potential->value(), echo, ngall);
              if (echo > 2) std::printf("# potential has %ld values, expect %d x %d x %d = %ld\n",
                                                        Veff.size(), ng[0], ng[1], ng[2], ngall);
              if (Veff.size() != ngall) {
                  if (echo > 0) std::printf("# expected %d*%d*%d = %ld potential values but found %ld\n",
                                                        ng[2], ng[1], ng[0], ngall, Veff.size());
                  auto const empty_keyword = "green_input.empty.potential";
                  auto const empty = control::get(empty_keyword, 1.); // 0:error, 1:warn, 2:okay
                  if (empty > 0) { // ok
                      Veff.resize(ngall, 0.0); // very useful for testing
                      if (empty < 2) warn("in file \"%s\" %d*%d*%d = %ld potential values set to zero due to +%s=%g",
                                           filename, ng[2], ng[1], ng[0], ngall, empty_keyword, empty);
                  } else {
                      error("expected %d*%d*%d = %ld potential values but found %ld, try +%s=1 to override",
                              ng[2], ng[1], ng[0], ngall, Veff.size(), empty_keyword);
                      return -1; // failure
                  }
              } // Veff.size != ngall
          } else warn("grid_Hamiltonian has no potential in file \"%s\"", filename);

      } else warn("no grid_Hamiltonian found in file \"%s\"", filename);

#endif // HAS_RAPIDXML
      return 0; // success
  } // load_Hamiltonian

  status_t test_loading(int const echo=0) {
      uint32_t ng[3]; // numbers of grid points
      int8_t bc[3]; // boundary conditions
      double hg[3]; // grid spacings
      std::vector<double> Veff;
      int natoms;
      std::vector<double> xyzZinso;
      std::vector<std::vector<double>> atom_mat;
      auto const filename = control::get("hamiltonian.file", "Hmt.json");
      return load_Hamiltonian(ng, bc, hg, Veff, natoms, xyzZinso, atom_mat, filename, echo);
  } // test_loading

  status_t all_tests(int echo) {
      status_t stat(0);
      stat += test_loading(echo);
      return stat;
  } // all_tests

} // namespace green_input
