#pragma once

#ifdef HAS_RAPIDXML
//   #include <string>
  #include <cstdlib> // std::atof, std::strtod
//   #include <fstream>
//   #include <streambuf>
  #include <cstring> // std::strcmp

  // from https://sourceforge.net/projects/rapidxml/files/latest/download
  #include "tools/rapidxml/rapidxml.hpp" // ::xml_document<>
  #include "tools/rapidxml/rapidxml_utils.hpp" // ::file<>
#endif

#include <cerrno> // errno, ERANGE


#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "recorded_warnings.hxx" // warn

#include "green_function.hxx" // ::construct_Green_function

namespace xml_reading {

#ifdef HAS_RAPIDXML
    char const empty_string[] = "";
  
    inline char const * find_attribute(
          rapidxml::xml_node<> const *node
        , char const *const name
        , char const *const default_value=""
        , int const echo=0
    ) { 
        if (nullptr == node) return empty_string;
        for (auto attr = node->first_attribute(); attr; attr = attr->next_attribute()) {
            if (0 == std::strcmp(name, attr->name())) {
                return attr->value();
            } // found
        } // attr
        return default_value;
    } // find_attribute 
  
    inline rapidxml::xml_node<> const * find_child(
          rapidxml::xml_node<> const *node
        , char const *const name
        , int const echo=0
    ) { 
        if (nullptr == node) return nullptr;
        for (auto child = node->first_node(); child; child = child->next_sibling()) {
            if (0 == std::strcmp(name, child->name())) {
                return child;
            } // found
        } // attr
        return nullptr;
    } // find_attribute 
#endif

    template <typename real_t>
    std::vector<real_t> read_sequence(
          char const *sequence
        , int const echo=0
        , size_t const reserve=0
    ) {
        char *end;
        char const *seq{sequence};
        std::vector<real_t> v;
        v.reserve(reserve);
        for (double f = std::strtod(seq, &end); seq != end; f = std::strtod(seq, &end)) {
            seq = end;
            if (errno == ERANGE){
                warn("range error, got %g", f);
                errno = 0;
            } else {
                v.push_back(real_t(f));
            }
        } // f
        return v;
    } // read_sequence

  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_xml_reader(int const echo=0) {
      status_t stat(0);
#ifdef HAS_RAPIDXML
      rapidxml::file<> infile("Hmt.xml");

      rapidxml::xml_document<> doc;
      doc.parse<0>(infile.data());

      double hg[3] = {1, 1, 1};
      int ng[3] = {0, 0, 0};
      std::vector<double> Veff;
      std::vector<double> xyzZinso;
      std::vector<std::vector<double>> atom_mat;
      int natoms{0};

      auto const grid_Hamiltonian = doc.first_node("grid_Hamiltonian");
      if (grid_Hamiltonian) {
          auto const sho_atoms = find_child(grid_Hamiltonian, "sho_atoms", echo);
          if (sho_atoms) {
              auto const number = find_attribute(sho_atoms, "number", "0", echo);
              if (echo > 5) printf("# found number=%s\n", number);
              natoms = std::atoi(number);
              xyzZinso.resize(natoms*8);
              atom_mat.resize(natoms);
              int ia{0};
              for (auto atom = sho_atoms->first_node(); atom; atom = atom->next_sibling()) {
                  auto const gid = find_attribute(atom, "gid", "-1");
                  if (echo > 5) printf("# <%s gid=%s>\n", atom->name(), gid);
                  xyzZinso[ia*8 + 4] = std::atoi(gid);

                  double pos[3] = {0, 0, 0};
                  auto const position = find_child(atom, "position", echo);
                  for(int d = 0; d < 3; ++d) {
                      char axyz[] = {0, 0}; axyz[0] = 'x' + d; // "x", "y", "z"
                      auto const value = find_attribute(position, axyz);
                      if (*value != '\0') {
                          pos[d] = std::atof(value);
                          if (echo > 5) printf("# %s = %.15g\n", axyz, pos[d]);
                      } // value != ""
                      xyzZinso[ia*8 + d] = pos[d];
                  } // d

                  auto const projectors = find_child(atom, "projectors", echo);
                  int numax{-1};
                  {
                      auto const value = find_attribute(projectors, "numax", "-1");
                      if (*value != '\0') {
                          numax = std::atoi(value);
                          if (echo > 5) printf("# numax= %d\n", numax);
                      } // value != ""
                  }
                  double sigma{-1};
                  {
                      auto const value = find_attribute(projectors, "sigma", "-1");
                      if (*value != '\0') {
                          sigma = std::atof(value);
                          if (echo > 5) printf("# sigma= %g\n", sigma);
                      } // value != ""
                  }
                  xyzZinso[ia*8 + 5] = numax;
                  xyzZinso[ia*8 + 6] = sigma;
                  int const nSHO = sho_tools::nSHO(numax);
                  atom_mat[ia].resize(2*nSHO*nSHO);
                  for(int h0s1 = 0; h0s1 < 2; ++h0s1) {
                      auto const matrix_name = h0s1 ? "overlap" : "hamiltonian";
                      auto const matrix = find_child(atom, matrix_name, echo);
                      if (matrix) {
                          if (echo > 22) printf("# %s.values= %s\n", matrix_name, matrix->value());
                          auto const v = read_sequence<double>(matrix->value(), echo, nSHO*nSHO);
                          if (echo > 5) printf("# %s matrix has %d values, expect %d x %d = %d\n",
                              matrix_name, v.size(), nSHO, nSHO, nSHO*nSHO);
                          assert(v.size() == nSHO*nSHO);
                          set(atom_mat[ia].data() + h0s1*nSHO*nSHO, nSHO*nSHO, v.data()); // copy
                      } else warn("atom with global_id=%s has no %s matrix!", gid, matrix_name);
                  } // h0s1
                  ++ia; // count up the number of atoms
              } // atom
              assert(natoms == ia); // sanity
          } else warn("no <sho_atoms> found in grid_Hamiltonian");

          auto const spacing = find_child(grid_Hamiltonian, "spacing", echo);
          for(int d = 0; d < 3; ++d) {
              char axyz[] = {0, 0}; axyz[0] = 'x' + d; // "x", "y", "z"
              auto const value = find_attribute(spacing, axyz);
              if (*value != '\0') {
                  hg[d] = std::atof(value);
                  if (echo > 5) printf("# h%s = %.15g\n", axyz, hg[d]);
              } // value != ""
          } // d

          auto const potential = find_child(grid_Hamiltonian, "potential", echo);
          if (potential) {
              for(int d = 0; d < 3; ++d) {
                  char axyz[] = {'n', 0, 0}; axyz[1] = 'x' + d; // "nx", "ny", "nz"
                  auto const value = find_attribute(potential, axyz);
                  if (*value != '\0') {
                      ng[d] = std::atoi(value);
                      if (echo > 5) printf("# %s = %d\n", axyz, ng[d]);
                  } // value != ""
              } // d
              if (echo > 33) printf("# potential.values= %s\n", potential->value());
              Veff = read_sequence<double>(potential->value(), echo, ng[2]*ng[1]*ng[0]);
              if (echo > 5) printf("# potential has %d values, expect %d x %d x %d = %d\n",
                  Veff.size(), ng[0], ng[1], ng[2], ng[2]*ng[1]*ng[0]);
              assert(Veff.size() == ng[2]*ng[1]*ng[0]);
          } else warn("grid_Hamiltonian has no potential!");

      } else warn("no grid_Hamiltonian found");

      return green_function::construct_Green_function(ng, hg, Veff, xyzZinso, atom_mat, echo);
#else
      warn("Unable to check usage of rapidxml when compiled without -D HAS_RAPIDXML");
#endif
      return stat;
  } // test_xml_reader

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_xml_reader(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace xml_reading
