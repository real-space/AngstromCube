#pragma once

#include "status.hxx" // status_t

#include <cstdio> // std::fprintf, std::sprintf, stdout, stderr, std::fflush
#include <utility> // std::forward



namespace recorded_warnings {

// #define warn(...) std::sprintf(recorded_warnings::_new_warning(__FILE__, __LINE__, __func__), __VA_ARGS__);
  char* _new_warning(char const *file, int const line, char const *func); // hidden, please use the macro above

#define error(...) { \
    recorded_warnings::_print_error_message(stdout, __FILE__, __LINE__, __VA_ARGS__ ); \
    recorded_warnings::_print_error_message(stderr, __FILE__, __LINE__, __VA_ARGS__ ); \
    exit(__LINE__); }
    
  template <class... Args>
  void _print_error_message(FILE* os, char const *srcfile, int const srcline, Args &&... args) {
        std::fprintf(os, "\n\n# Error in %s:%i  Message:\n#   ", srcfile, srcline);
        std::fprintf(os, std::forward<Args>(args)...);
        std::fprintf(os, "\n\n");
        std::fflush(os);
  } // _print_error_message
  
#define warn(...) recorded_warnings::_print_warning_message(__FILE__, __LINE__, __func__, __VA_ARGS__);
  template <class... Args>
  int _print_warning_message(char const *srcfile, int const srcline, char const* func, Args &&... args) {
      auto const str = _new_warning(__FILE__, __LINE__, __func__);
      // generate the warning message
      int const nchars = std::sprintf(str, std::forward<Args>(args)...);
      return nchars;
  } // _print_warning_message
  
  status_t show_warnings(int const echo=1);

  status_t clear_warnings(int const echo=1);

  status_t all_tests(int const echo=0);

} // namespace recorded_warnings
