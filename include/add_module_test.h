// This header serves to #define add_module_test needed in main.cxx, green.cxx and green_tests.cu

// requirements: (the following include statements and varaiable definitions must
//                appear earlier in the file where 'add_module_test.h' is included)
//    #include <cstdio> // std::printf
//    #include "simple_timer.hxx" // SimpleTimer
//    std::vector<std::tuple<char const*, double, status_t>> results;
//    std::string const input_name;
//    int const echo;
//    bool const all;
//    bool const show;

#define   add_module_test(MODULE_NAME) {                                            \
              auto const module_name = #MODULE_NAME;                                \
              if (all || (input_name == module_name)) {                             \
                  SimpleTimer timer(module_name, 0, "", 0);                         \
                  if (echo > 3) std::printf("\n\n\n# ============= Module test"     \
                     " for %s ==================\n\n", module_name);                \
                  auto const stat = show ? 0 : MODULE_NAME::all_tests(echo);        \
                  double const time = timer.stop();                                 \
                  results.push_back(std::make_tuple(module_name, time, stat));      \
              }                                                                     \
          } // add_module_test

// Please remember to #undef add_module_test after usage
