#pragma once

#include "status.hxx" // status_t

#include <cstdio> // fprintf, std::sprintf
#include <utility> // std::forward

#define warn(...) std::sprintf(recorded_warnings::_new_warning(__FILE__, __LINE__, __func__), __VA_ARGS__); 

#define error(...) { \
    recorded_warnings::_print_error_message(stdout, __FILE__, __LINE__, __VA_ARGS__ ); \
    recorded_warnings::_print_error_message(stderr, __FILE__, __LINE__, __VA_ARGS__ ); \
    exit(__LINE__); }

namespace recorded_warnings {

  char* _new_warning(char const *file, int const line, char const *func); // hidden, please use the macro above

  template <class... Args>
  void _print_error_message(FILE* os, char const *srcfile, int const srcline, Args &&... args) {
        fprintf(os, "\n\n# Error in %s:%i Message:\n#   ", srcfile, srcline);
        fprintf(os, std::forward<Args>(args)...);
        fprintf(os, "\n\n");
        fflush(os);
  } // _print_error_message

  status_t show_warnings(int const echo=1);

  status_t clear_warnings(int const echo=1);

  status_t all_tests(int const echo=0);

} // namespace recorded_warnings
