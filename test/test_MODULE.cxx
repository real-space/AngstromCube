// g++ -std=c++11 -O0 -g -pedantic -Wall -Wno-format-security -Wno-format test_MODULE.cxx && ./a.out
#include <cstdio> // std::printf
#include <cassert> // assert
#include "status.hxx" // status_t
#include "control.hxx" // ::command_line_interface, ::get
// #include "recorded_warnings.hxx" // warn, error
#include "MODULE.hxx" // ::all_tests
int main(int const argc, char const *argv[]) {
    status_t stat(0);
    for (int iarg = 1; iarg < argc; ++iarg) {
        assert(nullptr != argv[iarg]);
        if ('+' == *argv[iarg]) { // char #0 of command line argument #iarg
            stat += control::command_line_interface(argv[iarg] + 1, iarg); // start after the '+' char
        } // '+'
    } // iarg
    stat += MODULE::all_tests(int(control::get("verbosity", 5.)));
    if (0 != int(stat)) std::printf("\n# %s: all_tests = %i\n", __FILE__, int(stat));
    return int(stat);
} // main
