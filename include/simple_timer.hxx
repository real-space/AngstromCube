#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <chrono> // std::chrono::high_resolution_clock
#include <string> // std::string

#include "status.hxx" // status_t

inline char const * strip_path(char const *const str=nullptr, char const search='/') {
    // search a string for e.g. '/' in path names, return pointer to the following chars
    if (nullptr == str) return "";
    char const *pt = str;
    for (char const *c = str; *c != '\0'; ++c) {
        if (*c == search) { pt = c + 1; }
    } // c
    return pt;
} // strip_path


class SimpleTimer {
    // This timer object prints the time elapsed between construction and destruction 
    // into stdout immediately when destructed. See test_usage() below.
    // When stop() is called, the timer returns the elapsed time in seconds as double.
private: // members
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::string file;
    std::string func;
    int line;
    int echo;

public:

    SimpleTimer(char const *sourcefile, int const sourceline=0, char const *function=nullptr, int const echo=2)
    : file(strip_path(sourcefile)), func(function), line(sourceline), echo(echo) {
        start_time = std::chrono::high_resolution_clock::now(); // start
    } // constructor

    double stop(int const stop_echo=0) const { 
        auto const stop_time = std::chrono::high_resolution_clock::now(); // stop
        auto const musec = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time).count();
        if (stop_echo > 0) std::printf("# timer started at %s:%d %s took %.5f sec\n", file.c_str(), line, func.c_str(), 1e-6*musec);
        return 1e-6*musec;
    } // stop

    ~SimpleTimer() { stop(echo); } // destructor

}; // class SimpleTimer



namespace simple_timer {

#ifdef    NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS
  status_t all_tests(int const echo=0); // declaration only
#endif // NO_UNIT_TESTS

} // namespace simple_timer
