#include <cstdio> // std::printf, std::snprintf
#include <cassert> // assert
#include <string> // std::string
#include <cstdlib> // std::atof
#include <map> // std::map<T,T2>
#include <cstring> // std::strchr, std::strncpy
#include <cmath> // std::sqrt
#include <fstream> // std::fstream

#include "control.hxx" // default_echo_level
// declarations: set, get, command_line_interface, read_control_file, all_tests

#include "recorded_warnings.hxx" // warn

namespace control {

  inline char* find_equal_sign(char const * string) { return (char*)std::strchr(string, '='); }

  int constexpr MaxNameLength = 64; // max variable name length

  // hidden function: _manage_variables
  char const* _manage_variables(char const *name, char op='?', char const *value=nullptr, int const echo=0) {

      static std::map<std::string, std::string> _map; // hidden archive

      assert(nullptr == std::strchr(name, '=')); // make sure that there is no '=' sign in the name

      std::string const varname(name);
      if ('?' == op) { // get
          // query only
          auto const & oldvalue = _map[varname];
          if (echo > 7) std::printf("# control found \"%s\" = \"%s\"\n", name, oldvalue.c_str());
          assert((nullptr == value) && "value pointer must be null");
          return oldvalue.c_str();
      } else
      if ('=' == op) { // set
          assert((nullptr != value) && "value pointer may not be null");
          auto const newvalue = std::string(value);
          bool constexpr warn_redefines = true;
          if (warn_redefines) {
              auto const & oldvalue = _map[varname];
              auto const old = oldvalue.c_str();
              bool const redefined = (old && ('\0' != *old));
              if (echo > 7) {
                  std::printf("# control sets \"%s\"", name);
                  if (redefined) std::printf(" from \"%s\"", old);
                  std::printf(" to \"%s\"\n", value);
              } // echo
              if (redefined) {
                  warn("variable \"%s\" was redefined from \"%s\" to \"%s\"", name, old, value);
              } // redefined
          } // warn_redefines
          _map[varname] = newvalue;
          return value;
      } else
      if ('!' == op) { // list all defined variables
          if (echo > 0) {
              std::printf("# control has the following variables defined:\n");
              for (auto const & pair : _map) {
                  double const numeric = std::atof(pair.second.c_str());
                  char buffer[32]; std::snprintf(buffer, 31, "%.16e", numeric);
                  if (pair.second == buffer) {
                      std::printf("# %s=%g\n", pair.first.c_str(), numeric);
                  } else {
                      std::printf("# %s=%s\n", pair.first.c_str(), pair.second.c_str());
                  }
              } // pair
              std::printf("#\n\n");
          } // echo
      } else {
          error("allowed operations are '?','=','!' but no other, found op='%c'", op);
      }
      return nullptr;

  } // _manage_variables

  char const* set(char const *name, char const *value, int const echo) {
      if (echo > 5) std::printf("# control::set(\"%s\", \"%s\")\n", name, value);
      assert(nullptr != value);
      return _manage_variables(name, '=', value, echo);
  } // set<string>

  status_t command_line_interface(char const *statement, int const echo) {
      auto const equal = find_equal_sign(statement);
      if (nullptr != equal) {
          auto const equal_char = equal - statement;
          assert('=' == statement[equal_char]);
          char name[MaxNameLength]; // get a mutable string
          std::strncpy(name, statement, std::min(MaxNameLength, int(equal_char))); // copy the statement up to '='
          name[equal_char] = '\0'; // delete the '=' sign to mark the name
          char const *value = equal + 1; // everything after the '=' marker
          if (echo > 7) std::printf("# control::set(statement=\"%s\") found name=\"%s\", value=\"%s\"\n", statement, name, value);
          set(name, value, echo); // value comes after '='
          return 0;
      } else {
          warn("# ignored statement \"%s\", maybe missing \'=\'", statement);
          return 1; // error, no '=' sign given
      }
  } // command_line_interface

  char const* get(char const *name, char const *default_value) {
      int const echo = default_echo_level;
      auto const value = _manage_variables(name, '?', nullptr, echo);
      if (nullptr != value) {
          if ('\0' != *value) {
              if (echo > 5) std::printf("# control::get(\"%s\", default=\"%s\") = \"%s\"\n", name, default_value, value);
              return value;
          }
      }
      if (echo > 5) std::printf("# control::get(\"%s\") defaults to \"%s\"\n", name, default_value);
      bool constexpr set_on_get = true;
      if (set_on_get) _manage_variables(name, '=', default_value, echo);
      return default_value;
  } // get<string>

  char const* set(char const *name, double const value, int const echo) {
      char buffer[32]; std::snprintf(buffer, 31, "%.16e", value);
      return set(name, buffer, echo);
  } // set<double>

  double get(char const *name, double const default_value) {
      char buffer[32]; std::snprintf(buffer, 31, "%.16e", default_value);
      return std::atof(get(name, buffer));
  } // get<double>


  std::string left_trim(std::string const & s)  {
      std::string const WhiteSpaceChars = " \n\r\t\f\v";
      size_t const start = s.find_first_not_of(WhiteSpaceChars);
      return (start == std::string::npos) ? "" : s.substr(start);
  } // left_trim

//   std::string trim(std::string const & s)  {
//       size_t const end = s.find_last_not_of(WhiteSpaceChars);
//       return (end == std::string::npos) ? "" : s.substr(0, end + 1);
//   } // right_trim

  status_t read_control_file(char const *filename, int const echo) {
      status_t stat(0);
      char const CommentChar = '#'; // commented lines in control files
      char const EchoComment = '!'; // comments that should appear in the log

      if (nullptr == filename) {
          if (echo > 1) std::printf("# no control file passed\n");
          return stat; // 0
      }
      
      if ('\0' == *filename) {
          if (echo > 1) std::printf("# no control file given\n");
          return stat; // 0
      }

      std::ifstream infile(filename, std::ifstream::in);
      if (infile.fail()) {
          warn("Unable to open file '%s' for reading controls", filename);
          return -1;
      } // failed

      if (echo > 1) std::printf("\n# reading '%s' ...\n\n", filename);
      int linenumber{0}, ncomments{0}, nempty{0};
      std::string line;
      while (std::getline(infile, line)) {
          ++linenumber;
          if (echo > 18) std::printf("# %s:%d\t  %s\n", filename, linenumber, line.c_str());
          auto const tlin = left_trim(line);
          if (CommentChar == tlin[0]) {
              ++ncomments;
              if (echo > 9) std::printf("# %s:%d\t comment: %s\n",
                                           filename, linenumber, tlin.c_str());
              if (EchoComment == tlin[1]) {
                  if (echo > 0) std::printf("%s\n", tlin.c_str());
              }
          } else if ("" == tlin) {
              ++nempty;
              if (echo > 11) std::printf("# %s:%d\t is empty\n", filename, linenumber);
          } else {
              if (echo > 8) std::printf("# %s:%d\t  %s\n", filename, linenumber, tlin.c_str());
              auto const line_stat = command_line_interface(tlin.c_str());
              if (line_stat) {
                  warn("failure parsing %s:%d \'%s\'", filename, linenumber, line.c_str());
              } else {
                 if (echo > 0) std::printf("# %s\n", tlin.c_str()); // show the valid commands
              } 
              stat += line_stat;
          }
      } // parse file line by line

      if (echo > 3) std::printf("# %s found %d comments in file '%s', status=%i\n\n",
                                   __func__, ncomments, filename, int(stat)); 

      return stat;
  } // read_control_file

  status_t show_variables(int const echo) {
      _manage_variables("", '!', nullptr, echo);
      return 0;
  } // show_variables

  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_control(int const echo=9) {
      if (echo > 1) std::printf("\n# %s: %s\n", __FILE__, __func__);
      status_t stat(0);
      // _manage_variables("a", "5"); // manual set
      // auto const five = _manage_variables("a"); // manual get
      // if (echo > 1) std::printf("# a = %s\n", five);
      // auto const b_defined = _manage_variables("b"); // get
      // if (echo > 1) std::printf("# b = %s (%p)\n", b_defined, b_defined);

      set("a", "5");
      stat += command_line_interface("a=6"); // launches warning about redefining
      auto const a = get("a", "defaultA");
      if (echo > 1) std::printf("# a = %s\n", a);

      auto const b = get("b", "defaultB");
      if (echo > 1) std::printf("# b = %s\n", b);

      auto const c = get("c", "3.14");
      if (echo > 1) std::printf("# c = %s\n", c);

      auto const c_double = get("c", 3.1415);
      if (echo > 1) std::printf("# c<double> = %g\n", c_double);

      return stat;
  } // test_control

  status_t test_precision(int const echo=3, int const nmax=100) {
      // check if there are rounding errors arising from the 
      //    ASCII representation of double precision numbers
      if (echo > 2) std::printf("\n# %s: %s\n", __FILE__, __func__);
      status_t stat(0);
      double d{.2};
      for (int i = 0; i < nmax; ++i) {
          set("d", d, echo); // warning about redefining "d" in the second iteration
          double const g = get("d", 1.);
          stat += (g != d);
          d *= std::sqrt(33/32.); // some irrational number close to 1
      } // i
      if (echo > 1) std::printf("# %s: for %i of %i cases double precision numbers are not retrieved\n", __func__, stat, nmax);
      return stat;
  } // test_precision

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_control(echo);
//       stat += test_precision(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace control
