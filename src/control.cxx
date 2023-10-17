// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, ::snprintf
#include <cassert> // assert
#include <string> // std::string
#include <cstdlib> // std::atof
#include <map> // std::map<T1,T2>
#include <tuple> // std::tuple<T1,...,Tn>, ::get
#include <cstring> // std::strchr, ::strncpy
#include <cmath> // std::sqrt
#include <fstream> // std::fstream

#include "control.hxx" // ::default_echo_level, ::set, ::get, ::command_line_interface, ::read_control_file, ::all_tests

#include "recorded_warnings.hxx" // warn

namespace control {

  int constexpr MaxNameLength = 64; // max variable name length

  // convert numeric values to strings without precision loss
  inline void double2string(char buffer[32], double const number) {
      std::snprintf(buffer, 32, "%.16e", number);
  } // double2string

  // convert a string to a scalar double value
  inline double string2double(char const *const string) {
      return std::atof(string);
  } // string2double


  int32_t constexpr default_value_tag = 2e9; // pass this as linenumber to _environment

  // hidden function:
  //    _environment(echo, name, value) --> set
  //    _environment(echo, name, value, linenumber) --> set_to_default
  //    _environment(echo, name) --> get
  //    _environment(echo) --> show_variables
  char const* _environment(
        int const echo
      , char const *const name=nullptr
      , char const *const value=nullptr
      , int const linenumber=0
  ) {
      bool constexpr warn_about_redefinitons = true;

      static std::map<std::string, std::tuple<std::string,uint32_t,int32_t>> _map; // hidden archive
      // the hidden archive is sorted by the name string and contains the value string,
      // an access counter and a reference to the line number to track where the value was set

      if (nullptr != name) {
          assert(nullptr == std::strchr(name, '=')); // make sure that there is no '=' sign in the name

          std::string const varname(name);
          auto & tuple = _map[varname];
          if (nullptr != value) {

              // set
              if (warn_about_redefinitons) {
                  auto const oldvalue = std::get<0>(tuple).c_str();
                  assert(nullptr != oldvalue);
                  bool const redefined = ('\0' != *oldvalue);
                  if (echo > 7) {
                      std::printf("# control sets \"%s\"", name);
                      if (redefined) std::printf(" from \"%s\"", oldvalue);
                      std::printf(" to \"%s\"\n", value);
                  } // echo
                  if (redefined) {
                      warn("variable \"%s\" was redefined from \"%s\" to \"%s\"", name, oldvalue, value);
                  } // redefined
              } // warn_about_redefinitons
              std::get<0>(tuple) = value;
              std::get<1>(tuple) = (default_value_tag == linenumber); // counter how many times this variable was evaluated: init as 1 for defaults, 0 otherwise
              std::get<2>(tuple) = linenumber; // store line number in input file
                                               // or (if negative) command line argument number
              return value;

          } else { // value

              // get
              auto const oldvalue = std::get<0>(tuple).c_str();
              ++std::get<1>(tuple); // increment reading counter
              if (echo > 7) std::printf("# control found \"%s\" = \"%s\"\n", name, oldvalue);
              return oldvalue;

          } // value

      } else { // name

          // show_variables
          if (echo) {
              int const show = echo;
              std::printf("\n# control.show=%d (0:none, 1:minimal, 2:unused, 4:defaults, negative for details)\n", show);
              bool const show_unused  = std::abs(show) & 0x2; // all or only accessed ones
              bool const show_default = std::abs(show) & 0x4; // list also variables that are at their default value
              bool const show_details = (show < 0); // show access count and is_default
              std::printf("# control has the following variables defined:\n#\n");
              int listed{0};
              for (auto const & pair : _map) {
                  auto const times_used = std::get<1>(pair.second); // how many times was this value used?
                  if (show_unused || times_used > 0) {
                      auto const line = std::get<2>(pair.second);
                      bool const is_default = (default_value_tag == line);
                      if (show_default || !is_default) {
                          if (show_details) {
                              std::printf("# used %dx, %s %3d\t", times_used,
                                  is_default?"def ":((line > 0)?"argv":(line?"line":"set ")),
                                  is_default?0:std::abs(line));
                          } // show_details
                          auto const & string = std::get<0>(pair.second);
                          double const numeric = string2double(string.c_str());
                          char buffer[32]; double2string(buffer, numeric);
                          if (string == buffer) {
                              // can be parsed as double, print with %g format
                              std::printf("# %s=%g\n", pair.first.c_str(), numeric);
                          } else {
                              std::printf("# %s=%s\n", pair.first.c_str(), string.c_str());
                          }
                          ++listed;
                      } // variable is at its default value
                  } // variable has been used or we want to show all variables defined
              } // pair
              std::printf("#\n# %d variables listed for control.show=%d\n", listed, show);
          } // show
          return nullptr;

      } // name

  } // _environment

  // define a pair of strings (name,value) for the variable environment
  void set(char const *const name, char const *const value, int const echo) {
      if (echo > 5) std::printf("# control::set(\"%s\", \"%s\")\n", name, value);
      assert(nullptr != name  && "control::set(name, value) needs a valid string as name!");
      assert(nullptr != value && "control::set(name, value) needs a valid string as value!");
      _environment(echo, name, value); // set
  } // set<string>

  // look up a name in the variable environment, return default string if not defined
  char const* get(char const *const name, char const *const default_value) {
      assert(nullptr != name && "control::get(name, default_value) needs a valid string as name!");
      int const echo = control::default_echo_level;
      auto const value = _environment(echo, name); // get
      if (nullptr != value && '\0' != *value) {
          if (echo > 5) std::printf("# control::get(\"%s\", default=\"%s\") = \"%s\"\n", name, default_value, value);
          return value;
      } else {
          if (echo > 5) std::printf("# control::get(\"%s\") defaults to \"%s\"\n", name, default_value);
          return _environment(echo, name, default_value, default_value_tag); // set_to_default
      }
  } // get<string>

  // print a list of all variables. +control.show decides by (1:minimal, 2:unused, 4:defaults)
  status_t show_variables(int const echo) {
      _environment(echo); // show_variables
      return 0;
  } // show_variables

  status_t command_line_interface(char const *const statement, int const iarg, int const echo) {
      auto const equal = (char const*)std::strchr(statement, '='); // find the position of '=' in statement
      if (nullptr == equal) {
          warn("ignored statement \"%s\", maybe missing \'=\'", statement);
          return 1; // error, no '=' sign given
      } // no '='-sign in statement
      auto const equal_char = equal - statement;
      assert('=' == statement[equal_char]);
      char name[MaxNameLength]; // get a mutable string
      std::strncpy(name, statement, std::min(MaxNameLength, int(equal_char))); // copy the statement up to '='
      name[equal_char] = '\0'; // delete the '=' sign to mark the name
      char const *const value = equal + 1; // everything after the '=' marker
      if (echo > 7) std::printf("# control::set(statement=\"%s\") found name=\"%s\", value=\"%s\"\n", statement, name, value);
      _environment(echo, name, value, iarg); // set
      return 0;
  } // command_line_interface

  // define a (name,value) pair, value is converted from double to a string
  void set(char const *const name, double const value, int const echo) {
      /* This version of set is only used in the tests below */
      char buffer[32]; double2string(buffer, value);
      return set(name, buffer, echo);
  } // set<double>

  // look up a name in the variable environment, return default double if not defined
  double get(char const *const name, double const default_value) {
      char buffer[32]; double2string(buffer, default_value);
      return string2double(get(name, buffer));
  } // get<double>

  // specify entire vectors of doubles, e.g. name.x, name.y, ...
  double get(
        double vec[] // result array
      , char const *const name // base keyword
      , char const *const xyz // ="xyz" // single characters to append to the base keyword after a '.'
      , double const default_value // =0 // default value
  ) {
      double const def_val = get(name, default_value);
      char name_x[96];
      for (int d = 0; '\0' != xyz[d]; ++d) {
          std::snprintf(name_x, 96, "%s.%c", name, xyz[d]);
          vec[d] = get(name_x, def_val);
      } // d
      return def_val;
  } // get

  // remove whitespace character at the start of a string
  std::string left_trim(std::string const & s)  {
      std::string const WhiteSpaceChars = " \n\r\t\f\v";
      size_t const start = s.find_first_not_of(WhiteSpaceChars);
      return (start == std::string::npos) ? "" : s.substr(start);
  } // left_trim

//   // remove whitespace character at the end of a string
//   std::string trim(std::string const & s)  {
//       size_t const end = s.find_last_not_of(WhiteSpaceChars);
//       return (end == std::string::npos) ? "" : s.substr(0, end + 1);
//   } // right_trim

  // read definitions of variables from an input file
  status_t read_control_file(char const *const filename, int const echo) {
      status_t stat(0);
      char const CommentChar = '#'; // commented lines in control files
      char const EchoComment = '!'; // comments that should appear in the log

      assert(nullptr != filename);
      if ('\0' == *filename) {
          if (echo > 1) std::printf("# no control file given\n");
          return stat; // 0
      } // filename != ""

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
              auto const line_stat = command_line_interface(tlin.c_str(), -linenumber);
              if (line_stat) {
                  warn("failure parsing %s:%d \'%s\'", filename, linenumber, line.c_str());
              } else {
                 if (echo > 0) std::printf("# %s\n", tlin.c_str()); // show the valid commands
              }
              stat += line_stat;
          }
      } // parse file line by line
      if (echo > 7 && nempty > 0) std::printf("# %s has %d empty lines\n", filename, nempty);
      if (echo > 1) std::printf("\n");
      if (echo > 3) std::printf("# %s found %d comments in file '%s', status=%i\n\n",
                                   __func__, ncomments, filename, int(stat));
      return stat;
  } // read_control_file









#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_control(int const echo=9) {
      if (echo > 1) std::printf("\n# %s: %s\n", __FILE__, __func__);
      status_t stat(0);

      set("a", "5", echo); // use the string-setter routine
      if (echo > 1) std::printf("# a = %s\n", get("a", ""));

      set("b", 6., echo); // use the double-setter routine
      if (echo > 1) std::printf("# b = %s\n", get("b", ""));

      // or use the command line interface
      stat += command_line_interface("a=6", 99, echo); // launches warning about redefining "a"

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

  status_t test_precision(int const echo=3, int const nmax=106) {
      // check if there are rounding errors arising from the
      //    ASCII representation of double precision numbers
      if (echo > 2) std::printf("\n# %s: %s\n", __FILE__, __func__);
      status_t stat(0);
      double d{0.2}; // choose 1/5 which cannot be represented exactly in binary
      for (int i = 0; i < nmax; ++i) {
          set("d", d, echo/2); // warning about redefining "d" in the second iteration
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
      stat += test_precision(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace control
