#pragma once

#ifdef HAS_RAPIDXML
  #include <string>
  #include <fstream>
  #include <streambuf>
  #include <cstring> // std::strcmp
  #include "tools/rapidxml-1.13/rapidxml.hpp" // ::xml_document<>
  #include "tools/rapidxml-1.13/rapidxml_utils.hpp" // ::file<>
#endif



#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "recorded_warnings.hxx" // warn

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

  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_xml_reader(int const echo=0) {
      status_t stat(0);
#ifdef HAS_RAPIDXML
      rapidxml::file<> infile("Hmt.xml");

      rapidxml::xml_document<> doc;
      doc.parse<0>(infile.data());

      auto const grid_Hamiltonian = doc.first_node("grid_Hamiltonian");
      auto const sho_atoms = find_child(grid_Hamiltonian, "sho_atoms", echo);
      if (sho_atoms) {
//        auto const number = find_attribute(sho_atoms, "number", "0", echo);
//        printf("# found number=%s\n", number);
          for (auto atom = sho_atoms->first_node(); atom; atom = atom->next_sibling()) {
              auto const gid = find_attribute(atom, "gid", "-1");
              printf("# <%s gid=%s>\n", atom->name(), gid);
              
              auto const position = find_child(atom, "position", echo);
              for(int d = 0; d < 3; ++d) {
                  char axyz[] = {0, 0}; axyz[0] = 'x' + d; // "x", "y", "z"
                  auto const value = find_attribute(position, axyz);
                  if (*value != '\0') {
                      double const pxyz = std::atof(value);
                      printf("# %s = %.15g\n", axyz, pxyz);
                  } // value != ""
              } // d
              
              auto const projectors = find_child(atom, "projectors", echo);
              int numax{-1};
              double sigma{-1};
              {
                  auto const value = find_attribute(projectors, "numax", "-1");
                  if (*value != '\0') {
                      numax = std::atoi(value);
                      printf("# numax= %d\n", numax);
                  } // value != ""
              }
              {
                  auto const value = find_attribute(projectors, "sigma", "-1");
                  if (*value != '\0') {
                      sigma = std::atof(value);
                      printf("# sigma= %g\n", sigma);
                  } // value != ""
              }

              int const nSHO = sho_tools::nSHO(numax);
              for(int h0s1 = 0; h0s1 < 2; ++h0s1) {
                  auto const matrix_name = h0s1 ? "overlap" : "hamiltonian";
                  auto const matrix = find_child(atom, matrix_name, echo);
                  if (nullptr == matrix) {
                      warn("atom with global_id=%s has no %s matrix!", gid, matrix_name);
                  } else {
//                       printf("# %s.values= %s\n", matrix_name, matrix->value());
                  }
              } // h0s1

          } // atom
      } else warn("no <sho_atoms> found in grid_Hamiltonian");

      { 
          {
              auto const spacing = find_child(grid_Hamiltonian, "spacing", echo);
              for(int d = 0; d < 3; ++d) {
                  char axyz[] = {0, 0}; axyz[0] = 'x' + d; // "x", "y", "z"
                  auto const value = find_attribute(spacing, axyz);
                  if (*value != '\0') {
                      double const hxyz = std::atof(value);
                      printf("# h%s = %.15g\n", axyz, hxyz);
                  } // value != ""
              } // d
          }
          {
              auto const potential = find_child(grid_Hamiltonian, "potential", echo);
              for(int d = 0; d < 3; ++d) {
                  char axyz[] = {'n', 0, 0}; axyz[1] = 'x' + d; // "nx", "ny", "nz"
                  auto const value = find_attribute(potential, axyz);
                  if (*value != '\0') {
                      int const nxyz = std::atoi(value);
                      printf("# %s = %d\n", axyz, nxyz);
                  } // value != ""
              } // d
              if (nullptr == potential) {
                  warn("grid_Hamiltonian has no potential!");
              } else {
//                   printf("# potential.values= %s\n", potential->value());
              }              
          }
      }
              
      

            
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
