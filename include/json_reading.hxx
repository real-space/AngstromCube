#pragma once

#include <cstdio> // std::printf
#include <vector> // std::vector<T>

#ifdef  HAS_RAPIDJSON
  #include <cassert> // assert
  #include <algorithm> // std::min
  #include <fstream> // std::ifstream
  #include <sstream> // std::ostringstream
  // git clone https://github.com/Tencent/rapidjson.git
  #include "rapidjson/include/rapidjson/document.h" // rapidjson::Document
#endif // HAS_RAPIDJSON

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "recorded_warnings.hxx" // warn
#include "sho_tools.hxx" // ::nSHO

namespace json_reading {

  template <class C>
  std::vector<double> read_json_array(C const & object, char const *const name="?", int const echo=0) {
      std::vector<double> vector(0);
      if (object.IsArray()) {
          auto const array = object.GetArray();
          auto const array_size = array.Size();
          vector.resize(array_size);
          for (unsigned i = 0; i < array_size; ++i) {
              vector[i] = array[i].GetDouble();
              if (echo > 5) std::printf("# %s %s[%i] = %g\n", __func__, name, i, vector[i]);
          } // i
          if (echo > 2) std::printf("# %s %s has %ld entries\n", __func__, name, size_t(array_size));
      } else {
          if (echo > 0) std::printf("# %s object(\"%s\") is not an array!", __func__, name);
      }
      return vector;
  } // read_json_array

  template <class C>
  std::vector<std::vector<double>> read_json_matrix(C const & object, char const *const name="?", int const echo=0) {
      std::vector<std::vector<double>> matrix(0);
      if (object.IsArray()) {
          auto const array = object.GetArray();
          auto const array_size = array.Size();
          matrix.resize(array_size);
          char row_name[32];
          for (unsigned i = 0; i < array_size; ++i) {
              std::snprintf(row_name, 32, "%s[%i]", name, i);
              matrix[i] = read_json_array(array[i], row_name, echo/2);
              if (echo > 5) std::printf("# %s %s has %ld entries\n", __func__, row_name, matrix[i].size());
          } // i
          if (echo > 2) std::printf("# %s %s has %ld entries\n", __func__, name, size_t(array_size));
      } else {
          if (echo > 0) std::printf("# %s object(\"%s\") is not an array!", __func__, name);
      }
      return matrix;
  } // read_json_matrix


  inline status_t load_Hamiltonian(
        uint32_t grid[3] // numbers of grid points
      , int8_t boundary[3] // boundary conditions
      , double spacing[3] // grid spacings
      , std::vector<double> & potential
      , int & natoms
      , std::vector<double> & xyzZinso
      , std::vector<std::vector<double>> & atom_mat
      , char const *const filename="Hmt.xml" // input
      , int const echo=0 // log-level
  ) {
#ifdef  HAS_RAPIDJSON
      std::ifstream const file_stream(filename); // taking file as inputstream
      if (file_stream) {
          if (echo > 7) std::printf("# %s: opening \"%s\"\n", __func__, filename);
      } else {
          if (echo > 0) std::printf("# %s: failed to open \"%s\"\n", __func__, filename);
          return -1;
      } // file_stream
      std::ostringstream ss;
      ss << file_stream.rdbuf(); // reading data
      auto const infile = ss.str();

      rapidjson::Document doc;
      if (doc.Parse(infile.data()).HasParseError()) {
          if (echo > 0) std::printf("# %s: rapidjson.Parse failed, try in situ parsing!\n", __func__);
          std::vector<char> json(infile.size() + 1, '\0');
          std::strncpy(json.data(), infile.data(), json.size());
          if (doc.ParseInsitu(json.data()).HasParseError()) {
              if (echo > 0) std::printf("# %s: rapidjson.ParseInsitu failed for \"%s\"!\n", __func__, filename);
              return -1;
          } // ParseInsitu
      } // Parse
      assert(doc.IsObject());

      if (doc.HasMember("comment")) {
          assert(doc["comment"].IsString());
          if (echo > 0) std::printf("# %s: comment = \"%s\"\n", filename, doc["comment"].GetString());
      } // has comment

      if (doc.HasMember("spacing")) {
          auto const array = read_json_array(doc["spacing"], "spacing", echo);
          for (unsigned i = 0; i < array.size(); i++) {
              if (i < 3) spacing[i] = array[i];
          } // i
      } // has spacing
      if (echo > 0) std::printf("# %s: grid spacing %g %g %g Bohr\n", filename, spacing[0], spacing[1], spacing[2]);

      if (doc.HasMember("boundary")) {
          auto const array = read_json_array(doc["boundary"], "boundary", echo);
          for (unsigned i = 0; i < array.size(); i++) {
              if (i < 3) boundary[i] = int(array[i]);
          } // i
      } // has boundary
      if (echo > 0) std::printf("# %s: cell boundary %d %d %d\n", filename, boundary[0], boundary[1], boundary[2]);

      if (doc.HasMember("potential")) {
          assert(doc["potential"].IsObject());
          auto const pot = doc["potential"].GetObject();
          if (pot.HasMember("grid")) {
              auto const array = read_json_array(pot["grid"], "grid", echo);
              for (unsigned i = 0; i < array.size(); ++i) {
                  if (i < 3) grid[i] = int(array[i]);
              } // i
              if (echo > 0) std::printf("# %s: potential grid = %d %d %d\n", filename, grid[0], grid[1], grid[2]);
          } // has grid
          if (pot.HasMember("values")) {
              potential = read_json_array(pot["values"], "potential values", echo/2); // no echo here, too long
          } // has values
      } // has potential
      {
          auto const nzyx = grid[2]*size_t(grid[1])*size_t(grid[0]);
          if (potential.size() != nzyx) {
              warn("in %s found %ld potential values, expect %ld, set missing to zero\n", filename, potential.size(), nzyx);
              potential.resize(nzyx, 0.0);
          } // sizes differ
      }

      double const cell[] = {grid[0]*spacing[0], grid[1]*spacing[1], grid[2]*spacing[2]};
      double const rel[] = {1/cell[0], 1/cell[1], 1/cell[2]}; // for relative positions

      natoms = 0;
      if (doc.HasMember("sho_atoms")) {
          assert(doc["sho_atoms"].IsObject());
          auto const sho_atoms = doc["sho_atoms"].GetObject();
          if (sho_atoms.HasMember("number")) {
              natoms = sho_atoms["number"].GetInt();
          }
          if (sho_atoms.HasMember("atoms")) {
              assert(sho_atoms["atoms"].IsArray());
              auto const atoms = sho_atoms["atoms"].GetArray();
              assert(natoms == atoms.Size());
              xyzZinso.resize(8*natoms, 0.0);
              atom_mat.resize(natoms);
              std::vector<int32_t> atom_id(natoms, -1);
              std::vector<int8_t> numax(natoms, -1);
              std::vector<double> sigma(natoms, 1.);
              for (unsigned ia = 0; ia < natoms; ++ia) {
                  assert(atoms[ia].IsObject());
                  auto const atom = atoms[ia].GetObject();
                  int nsho{0};
                  if (atom.HasMember("atom_id")) {
                      if (atom["atom_id"].IsInt()) {
                          atom_id[ia] = atom["atom_id"].GetInt();
                          xyzZinso[ia*8 + 4] = atom_id[ia];
                          if (echo > 7) std::printf("# atom_id = %i\n", atom_id[ia]);
                      }
                  } // has atom_id
                  if (atom.HasMember("projectors")) {
                      assert(atom["projectors"].IsObject());
                      auto const prj = atom["projectors"].GetObject();
                      if (prj.HasMember("type")) {
                          assert(prj["type"].IsString());
                          auto const type = prj["type"].GetString();
                          if (echo > 7) std::printf("# atom #%i projector type = \"%s\"\n", atom_id[ia], type);
                      } // has type
                      if (prj.HasMember("numax")) {
                          assert(prj["numax"].IsInt());
                          numax[ia] = prj["numax"].GetInt();
                          xyzZinso[ia*8 + 5] = numax[ia];
                          nsho = sho_tools::nSHO(numax[ia]);
                          atom_mat[ia].resize(2*nsho*nsho, 0.0);
                      } // has numax
                      if (prj.HasMember("sigma")) {
                          assert(prj["sigma"].IsFloat());
                          sigma[ia] = prj["sigma"].GetDouble();
                          xyzZinso[ia*8 + 6] = sigma[ia];
                      } // has sigma
                  } // has projectors
                  for (int h0s1 = 0; h0s1 < 2; ++h0s1) {
                      auto const *const hs = h0s1 ? "overlap" : "hamiltonian";
                      if (atom.HasMember(hs)) {
                          auto const values = read_json_matrix(atom[hs], hs, echo/8);
                          double maxdev{0}; // measure deviation from a symmetric matrix
                          for (int i = 0; i < nsho; ++i) {
                              assert(nsho == values[i].size());
                              for (int j = 0; j < nsho; ++j) {
                                  maxdev = std::max(maxdev, std::abs(values[i][j] - values[j][i]));
                                  atom_mat[ia][(h0s1*nsho + i)*nsho + j] = values[i][j];
                              } // j
                          } // i
                          if (echo > 3) std::printf("# %s matrix of atom #%i has a max deviation of %.1e from symmetric\n",
                                                      hs, atom_id[ia], maxdev);
                          if (maxdev > 1e-6) warn("%s matrix of atom #%i has a max deviation of %.1e from symmetric",
                                                      hs, atom_id[ia], maxdev);
                      } // has matrix
                  } // {hamiltonian, overlap}
                  if (atom.HasMember("position")) {
                      auto const pos = read_json_array(atom["position"], "atom position", echo/4);
                      assert(3 == pos.size());
                      if (echo > 3) std::printf("# atom #%i numax= %i sigma= %.3f Bohr at position %g %g %g\n",
                              atom_id[ia], numax[ia], sigma[ia], pos[0]*rel[0], pos[1]*rel[1], pos[2]*rel[2]);
                      for (int d = 0; d < 3; ++d) { xyzZinso[ia*8 + d] = pos[d]; }
                      xyzZinso[ia*8 + 3] = 0; // number of protons in the core is unknown
                      xyzZinso[ia*8 + 7] = 0; // spare entry
                  } // has position

              } // ia
          } // has atoms
      } // has sho_atoms
      if (echo > 0) std::printf("# %s: found %d SHO atoms\n", filename, natoms);
      if (echo > 0) std::printf("# %s: cell size is %g %g %g Bohr\n", filename, cell[0], cell[1], cell[2]);

      return 0;
#else  // HAS_RAPIDJSON
      warn("Unable to check usage of rapidjson when compiled without -D HAS_RAPIDJSON", 0);
      return -1;
#endif // HAS_RAPIDJSON
  } // load_Hamiltonian

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_json_reader(int const echo=0) {
      uint32_t ng[3];
      int8_t bc[3];
      double hg[3];
      std::vector<double> Veff;
      int natoms{0};
      std::vector<double> xyzZinso;
      std::vector<std::vector<double>> atom_mat;
      return load_Hamiltonian(ng, bc, hg, Veff, natoms, xyzZinso, atom_mat, "Hmt.json", echo);
  } // test_json_reader

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_json_reader(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace json_reading
