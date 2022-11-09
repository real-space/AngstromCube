#pragma once

#include <cstdio> // std::printf
#ifdef  HAS_RAPIDJSON
  #include <algorithm> // std::min
  #include <cassert> // assert
  #include <fstream> // std::ifstream
  #include <sstream> // std::ostringstream
  // git clone https://github.com/Tencent/rapidjson.git
  #include "rapidjson/include/rapidjson/document.h" // rapidjson::Document
#endif // HAS_RAPIDJSON

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "recorded_warnings.hxx" // warn

namespace json_reading {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <class C>
  std::vector<double> read_json_array(C const & object, char const *const name="?", int const echo=0) {
      std::vector<double> vector(0);
      if (object.IsArray()) {
          auto const array = object.GetArray();
          auto const array_size = array.Size();
          vector.resize(array_size);
          for (unsigned i = 0; i < array_size; ++i) {
              vector[i] = array[i].GetDouble();
              if (echo > 0) std::printf("# %s %s[%i] = %g\n", __func__, name, i, vector[i]);
          } // i
      } else {
          if (echo > 0) std::printf("# %s object(\"%s\") is not an array!", __func__, name);
      }
      return vector;
  } // read_json_array


  inline status_t test_json_reader(int const echo=0) {
      status_t stat(0);
#ifdef  HAS_RAPIDJSON

      auto const filename = "Hmt.json";
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

      double spacing[] = {1, 1, 1};
      if (doc.HasMember("spacing")) {
          auto const array = read_json_array(doc["spacing"], "spacing", echo);
          for (unsigned i = 0; i < array.size(); i++) {
              if (i < 3) spacing[i] = array[i];
          } // i
      } // has spacing
      if (echo > 0) std::printf("# %s: grid spacing %g %g %g Bohr\n", filename, spacing[0], spacing[1], spacing[2]);

      int boundary[] = {0, 0, 0};
      if (doc.HasMember("boundary")) {
          auto const array = read_json_array(doc["boundary"], "boundary", echo);
          for (unsigned i = 0; i < array.size(); i++) {
              if (i < 3) boundary[i] = int(array[i]);
          } // i
      } // has boundary
      if (echo > 0) std::printf("# %s: cell boundary %d %d %d\n", filename, boundary[0], boundary[1], boundary[2]);

      int grid[] = {4, 4, 4};
      std::vector<double> potential(0);
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
              potential = read_json_array(pot["values"], "potential values", echo*0); // no echo here, too long
          } // has values
      } // has potential
      {
          auto const nzyx = grid[2]*size_t(grid[1])*size_t(grid[0]);
          if (potential.size() != nzyx) {
              warn("in %s found %ld potential values, expect %ld, set missing to zero\n", filename, potential.size(), nzyx);
              potential.resize(nzyx, 0.0);
          } // sizes differ
      }

      unsigned natoms{0};
      if (doc.HasMember("sho_atoms")) {
          assert(doc["sho_atoms"].IsObject());
          auto const sho_atoms = doc["sho_atoms"].GetObject();
          if (sho_atoms.HasMember("number")) {
              natoms = sho_atoms["number"].GetInt();
          }
          if (sho_atoms.HasMember("atoms")) {
              assert(sho_atoms["atoms"].IsArray());
              auto const atoms_object = sho_atoms["atoms"].GetObject();
              auto const atoms = atoms_object.GetArray();
              assert(natoms == atoms.Size());
              for (unsigned ia = 0; ia < natoms; ++ia) {
                  auto const atom = atoms[ia].GetObject();
                  assert(atom.IsObject());
              } // ia
          } // has atoms
      } // has sho_atoms
      if (echo > 0) std::printf("# %s: found %d SHO atoms\n", filename, natoms);

#else  // HAS_RAPIDJSON
      warn("Unable to check usage of rapidjson when compiled without -D HAS_RAPIDJSON", 0);
#endif // HAS_RAPIDJSON
      return stat;
  } // test_json_reader

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_json_reader(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace json_reading
