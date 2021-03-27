#pragma once

#ifdef HAS_RAPIDJSON
  #include <cassert> // assert
//   #include <cstdlib> // std::atof, std::strtod
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

  inline status_t test_json_reader(int const echo=0) {
      status_t stat(0);
#ifdef HAS_RAPIDJSON

      std::string infile;
      std::ifstream file_stream("Hmt.json"); // taking file as inputstream
      if (file_stream) {
          std::ostringstream ss;
          ss << file_stream.rdbuf(); // reading data
          infile = ss.str();
      } // file_stream

#if 0
      if (echo > 0) {
          char head[512];
          std::strncpy(head, infile.data(), std::min(infile.size(), size_t(511)));
          head[511] = '\0';
          std::printf("# %s: after reading json file:\n# %s\n", __func__, head);
      } // echo
#endif // NEVER

      rapidjson::Document doc;
      if (doc.Parse(infile.data()).HasParseError()) {
          if (echo > 0) std::printf("# %s: rapidjson.Parse failed, try in situ parsing!\n", __func__);
          std::vector<char> json(infile.size() + 1);
          std::strncpy(json.data(), infile.data(), json.size());
          if (doc.ParseInsitu(json.data()).HasParseError()) {
              if (echo > 0) std::printf("# %s: rapidjson.ParseInsitu failed!\n", __func__);
              return -1;
          } // ParseInsitu
      } // Parse
      assert(doc.IsObject()); // fails although documentation says this should not

      assert(doc["comment"].IsString());
      assert(doc["spacing"].IsArray());
      assert(doc["potential"].IsObject());
      assert(doc["sho_atoms"].IsObject());
      assert(doc.HasMember("potential"));
      assert(doc.HasMember("sho_atoms"));

#else
      warn("Unable to check usage of rapidjson when compiled without -D HAS_RAPIDJSON", 0);
#endif
      return stat;
  } // test_json_reader

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_json_reader(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace json_reading
