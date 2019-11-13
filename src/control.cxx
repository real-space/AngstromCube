#include <cstdio> // printf
#include <cassert> // assert
#include <string>
#include <cstdlib> // std::atof
#include <map> // std::map<T,T2>

#include "control.hxx"

// #define FULL_DEBUG
// #define DEBUG

namespace control {

  // public functions: set and get

  // hidden function: _manage_variables
  char const* _manage_variables(char const *name, char const *value=nullptr) {

      static std::map<std::string, std::string> map_; // hidden archive

      int constexpr echo = 9;

      assert(nullptr == strchr(name, '=')); // make sure that there is no '=' sign in the name

      auto const varname = std::string(name);
      if (nullptr == value) {
          // query
          auto const & oldvalue = map_[varname];
          if (echo > 7) printf("# control found \"%s\" = \"%s\"\n", name, oldvalue.c_str());
          return oldvalue.c_str();
      } else {
          // set
          auto const newvalue = std::string(value);
          if (echo > 7) {
              auto const & oldvalue = map_[varname];
              auto const old = oldvalue.c_str();
              printf("# control sets \"%s\" ", name);
              if (old) { if (0 != *old) printf("from \"%s\" ", old); }
              printf("to \"%s\"\n", value);
          } // echo
          map_[varname] = newvalue;
          return value;
      }

  } // _manage_variables


  char const* set(char const *name, char const *value, int const echo) {
      if (echo > 5) printf("# control::set(\"%s\", \"%s\")\n", name, value);
      return _manage_variables(name, value);
  } // set

  char const* get(char const *name, char const *default_value, int const echo) {
      auto const value = _manage_variables(name);
      if (0 == value[0]) {
          if (echo > 5) printf("# control::get(\"%s\") defaults to \"%s\"\n", name, default_value);
          return default_value;
      } else {
          if (echo > 5) printf("# control::get(\"%s\", default=\"%s\") = \"%s\"\n", name, default_value, value);
          return value;
      }
  } // get

  char const* set(char const *name, double const value, int const echo) {
      char buffer[32]; std::sprintf(buffer, "%30.20e", value);
      return set(name, buffer, echo);
  } // set<double>

  double get(char const *name, double const default_value, int const echo) {
      char buffer[32]; std::sprintf(buffer, "%30.20e", default_value);
      return std::atof(get(name, buffer, echo));
  } // get<double>

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_control(int const echo=9) {
    status_t stat = 0;
    // manage_variables("a", "5"); // set
    // auto const five = manage_variables("a"); // get
    // if (echo > 1) printf("# a = %s\n", five);
    // auto const defined = manage_variables("b"); // get
    // if (echo > 1) printf("# b = %s (%p)\n", defined, defined);

    set("a", "5");
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

  status_t all_tests() {
    auto status = 0;
    status += test_control();
    return status;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace control
