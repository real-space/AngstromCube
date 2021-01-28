#pragma once

#include <cstdio> // std::printf, std::snprintf, std::fprintf, std::fopen, std::close, std::remove
#include <vector> // std::vector<T>
#include <ctime> // std::strftime, ::tm, ::gmtime

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

// #include "green_function.hxx" // ::construct_Green_function

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
    } // find_child 
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

    inline status_t temporary_file_extension(char str[], size_t const count) {
        auto const now = std::time(nullptr);
        auto const bytes_written = std::strftime(str, count, "%Y%b%d_%H%M%S", std::gmtime(&now));
        return (bytes_written < 1); // report failure
    } // temporary_file_extension
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_xml_reader(int const echo=0) {
#ifdef HAS_RAPIDXML
      status_t stat(0);
      
// #define USE_TEMPORARY_FILE
#ifdef  USE_TEMPORARY_FILE

      // create a temporary xml file before trying to read it
      char filename[96], temporary[96];
      stat += temporary_file_extension(temporary, 95);
      std::snprintf(filename, 95, "test_%s.xml", temporary);
      if (echo > 0) printf("# %s: use test file \"%s\"\n", __func__, filename);

      std::FILE *const f = std::fopen(filename, "w");
      if (nullptr == f) {
          if (echo > 0) printf("# %s Error opening file %s for writing!\n", __func__, filename);
          return __LINE__;
      } // failed to open

      #define print2file(...) std::fprintf(f, __VA_ARGS__)
#else

      size_t const buffer = 2048;
      char mutable_string[buffer]; // internal file
      size_t nchars{0};

      #define print2file(...) nchars += std::snprintf(mutable_string + nchars, buffer - nchars, __VA_ARGS__);
#endif // USE_TEMPORARY_FILE

      print2file("<?xml version=\"%.1f\"?>\n", 1.0);
      print2file("<grid_Hamiltonian version=\"%.1f\">\n", 0.);
      print2file("  <!-- Units: Hartree and Bohr radii. -->\n");
      print2file("  <sho_atoms number=\"%d\">\n", 1);
      for(int ia = 0; ia < 1; ++ia) {
          print2file("    <atom gid=\"%i\">\n", ia);
          print2file("      <position x=\"%.6f\" y=\"%.6f\" z=\"%.6f\"/>\n", 0., 0., 0.);
          print2file("      <projectors type=\"sho\" numax=\"%d\" sigma=\"%.3f\"/>\n", 1, 1);
          int const nSHO = 4;
          for(int h0s1 = 0; h0s1 < 2; ++h0s1) {
              auto const tag = h0s1 ? "overlap" : "hamiltonian";
              print2file("      <%s>", tag);
              for(int i = 0; i < nSHO; ++i) {
                  print2file("\n        ");
                  for(int j = 0; j < nSHO; ++j) {
                      print2file(" %.1f", i + 0.1*j);
                  } // j
              } // i
              print2file("\n      </%s>\n", tag);
          } // h0s1
          print2file("    </atom>\n");
      } // ia
      print2file("  </sho_atoms>\n");

      print2file("  <spacing x=\"%g\" y=\"%g\" z=\"%g\"/>\n", .25, .25, .25);
      print2file("  <potential nx=\"%d\" ny=\"%d\" nz=\"%d\">", 2, 2, 2);
      int const nzyx = 2*2*2;
      for(int izyx = 0; izyx < nzyx; ++izyx) {
          if (0 == (izyx & 3)) print2file("\n    ");
          print2file(" %.6f", izyx*100./nzyx);
      } // ia
      print2file("\n  </potential>\n");
      print2file("</grid_Hamiltonian>\n");

#undef  print2file
#ifdef  USE_TEMPORARY_FILE

      std::fclose(f);
      if (echo > 3) printf("# file %s written\n", filename);
      
      // read the file in again
      rapidxml::file<> infile(filename);
      auto const mutable_string = infile.data();
      auto const nchars         = infile.size();

      // delete the temporary file after reading
      stat += std::remove(filename);

#endif // USE_TEMPORARY_FILE

      rapidxml::xml_document<> doc;
      // start to parse the test file
      doc.parse<0>(mutable_string);

      if (echo > 29) { // show where rapidxml replaced '<' and '=' by '\0' 
          printf("# %s: string after parsing (%ld chars):\n", __func__, nchars);
          for(size_t ic = 0; ic < nchars; ++ic) {
              char const c = mutable_string[ic];
              printf("%c", ('\0' == c)?' ':c); // print string characters 1-by-1 to replace '\0'-characters
          } // ic
          printf("# %s: end string\n", __func__);
      } // echo

      // data to structures to be filled:
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
                  int const nSHO = ((1 + numax)*(2 + numax)*(3 + numax))/6; // nSHO
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
                          for(int ij = 0; ij < nSHO*nSHO; ++ij) {
                              atom_mat[ia][h0s1*nSHO*nSHO + ij] = v[ij]; // copy
                          } // ij
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

      return stat;
#else
      warn("Unable to check usage of rapidxml when compiled without -D HAS_RAPIDXML");
      return STATUS_TEST_NOT_INCLUDED;
#endif
  } // test_xml_reader

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_xml_reader(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace xml_reading
