#pragma once


#ifdef HAS_RAPIDJSON
  #include <cassert> // assert
//   #include <cstdlib> // std::atof, std::strtod
  #include <fstream> // std::ifstream
  #include <sstream> // std::ostringstream

  // git clone https://github.com/Tencent/rapidjson.git
  #include "tools/rapidjson/include/rapidjson/document.h" // rapidjson::Document
#endif

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

      rapidjson::Document doc;
      doc.Parse(infile.data());
      assert(doc.IsObject());

      assert(doc["comment"].IsString());
      assert(doc["spacing"].IsArray());
      assert(doc["sho_atoms"].IsObject());
      assert(doc["potential"].IsObject());

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
