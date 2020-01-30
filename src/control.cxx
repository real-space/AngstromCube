#include <cstdio> // printf, sprintf
#include <cassert> // assert
#include <string>
#include <cstdlib> // std::atof
#include <map> // std::map<T,T2>
#include <cstring> // strchr, strncpy
#include <cmath> // std::sqrt

#include "control.hxx"

#include "recorded_warnings.hxx" // warn

// #define FULL_DEBUG
// #define DEBUG

namespace control {

  // public functions: set and get

  inline char* find_equal_sign(char const * string) { return (char*)strchr(string, '='); }

  int constexpr MaxNameLength = 64;
  
  // hidden function: _manage_variables
  char const* _manage_variables(char const *name, char const *value=nullptr, int const echo=0) {

      static std::map<std::string, std::string> _map; // hidden archive

      assert(nullptr == strchr(name, '=')); // make sure that there is no '=' sign in the name

      auto const varname = std::string(name);
      if (nullptr == value) {
          // query
          auto const & oldvalue = _map[varname];
          if (echo > 7) printf("# control found \"%s\" = \"%s\"\n", name, oldvalue.c_str());
          return oldvalue.c_str();
      } else {
          // set
          auto const newvalue = std::string(value);
          if (echo > 7) {
              auto const & oldvalue = _map[varname];
              auto const old = oldvalue.c_str();
              printf("# control sets \"%s\" ", name);
              if (old) { if (0 != *old) printf("from \"%s\" ", old); }
              printf("to \"%s\"\n", value);
          } // echo
          _map[varname] = newvalue;
          return value;
      }

  } // _manage_variables

  char const* set(char const *name, char const *value, int const echo) {
      if (echo > 5) printf("# control::set(\"%s\", \"%s\")\n", name, value);
      return _manage_variables(name, value, echo);
  } // set

  status_t cli(char const *statement, int const echo) {
      auto const equal = find_equal_sign(statement);
      if (nullptr != equal) {
          auto const equal_char = equal - statement;
          assert('=' == statement[equal_char]);
          char name[MaxNameLength]; // get a mutable string
          std::strncpy(name, statement, std::min(MaxNameLength, (int)equal_char)); // copy the statement up to '='
          name[equal_char] = 0; // delete the '=' sign to mark the name
          char const *value = equal + 1; // everything after the '=' marker
          if (echo > 7) printf("# control::set(statement=\"%s\") found name=\"%s\", value=\"%s\"\n", statement, name, value);
          set(name, value, echo); // value comes after '='
          return 0;
      } else {
          warn("# ignored statement \"%s\", maybe missing \'=\'", statement);
          return 1; // error, no '=' sign given
      }
  } // cli
  
  char const* get(char const *name, char const *default_value, int const echo) {
      auto const value = _manage_variables(name, nullptr, echo);
      if (nullptr != value) {
          if (0 != value[0]) {
              if (echo > 5) printf("# control::get(\"%s\", default=\"%s\") = \"%s\"\n", name, default_value, value);
              return value;
          }
      }
      if (echo > 5) printf("# control::get(\"%s\") defaults to \"%s\"\n", name, default_value);
      return default_value;
  } // get

  char const* set(char const *name, double const value, int const echo) {
      char buffer[32]; std::sprintf(buffer, "%.16e", value);
      return set(name, buffer, echo);
  } // set<double>

  double get(char const *name, double const default_value, int const echo) {
      char buffer[32]; std::sprintf(buffer, "%.16e", default_value);
      return std::atof(get(name, buffer, echo));
  } // get<double>

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_control(int const echo=9) {
    if (echo > 1) printf("\n# %s %s\n", __FILE__, __func__);
    status_t stat = 0;
    // manage_variables("a", "5"); // set
    // auto const five = manage_variables("a"); // get
    // if (echo > 1) printf("# a = %s\n", five);
    // auto const defined = manage_variables("b"); // get
    // if (echo > 1) printf("# b = %s (%p)\n", defined, defined);

    set("a", "5");
    cli("a=6");
    auto const a = get("a", "defaultA");
    if (echo > 1) printf("# a = %s\n", a);

    auto const b = get("b", "defaultB");
    if (echo > 1) printf("# b = %s\n", b);

    auto const c = get("c", "3.14");
    if (echo > 1) printf("# c = %s\n", c);

    auto const c_double = get("c", 3.1415);
    if (echo > 1) printf("# c<double> = %g\n", c_double);

    return stat;
  } // test_control


  status_t test_precision(int const echo=3) {
    // check if there are rounding errors arising from the 
    //    ASCII representation of double precision numbers
    if (echo > 2) printf("\n# %s %s\n", __FILE__, __func__);
    status_t stat = 0;
    auto const nmax = 100;
    double d = .2;
    for(int i = 0; i < nmax; ++i) {
        set("d", d, echo);
        double const g = get("d", 1., echo);
        stat += (g != d);
        d *= std::sqrt(33/32.); // some irrational number close to 1
    } // i
    if (echo > 1) printf("# %s: for %i of %i cases double precision numbers are not retrieved\n", __func__, stat, nmax);
    return stat;
  } // test_precision
  
  status_t all_tests(int const echo) {
    auto status = 0;
    status += test_control(echo);
    status += test_precision(echo);
    return status;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace control
