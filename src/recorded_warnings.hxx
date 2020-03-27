#pragma once

#include "status.hxx" // status_t

#include <cstdio> // printf, std::fprintf, std::sprintf, stdout, stderr, std::fflush
#include <utility> // std::forward, std::pair<T1,T1>

namespace recorded_warnings {

// #define warn(...) std::sprintf(recorded_warnings::_new_warning(__FILE__, __LINE__, __func__), __VA_ARGS__);
//   char* _new_warning(char const *file, int const line, char const *func); // hidden, please use the macro above
  
   std::pair<char*,int> _new_warning(char const *file, int const line, char const *func); // hidden, please use the macro above

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

  inline char const * after_last_slash(char const *path_and_file, char const slash='/') {
      auto const has_slash = std::strrchr(path_and_file, slash);
      return has_slash ? (has_slash + 1) : path_and_file;
  } // after_last_slash

  template <class... Args>
  int _print_warning_message(char const *srcfile, unsigned const srcline, char const* func, Args &&... args) {
      auto const str_int = _new_warning(srcfile, srcline, func);
      char*   message = str_int.first;
      int const flags = str_int.second;
      // generate the warning message
      int const nchars = std::sprintf(message, std::forward<Args>(args)...);

      // check decisions if the message should be printed
      if (flags & 1) { // warning to stdout
            printf("# Warning: %s\n", message);
            if(flags & 4) printf("# This warning will not be shown again!\n");
      } // message to stdout
      if (flags & 2) {
          std::fprintf(stderr, "%s:%d warn(\"%s\")\n", after_last_slash(srcfile), srcline, message);
      } // message to stderr
      if (flags & 1) printf("\n"); // give more optical weight to the warning lines

      return nchars;
  } // _print_warning_message
  
  status_t show_warnings(int const echo=1);

  status_t clear_warnings(int const echo=1);

  status_t all_tests(int const echo=0);

} // namespace recorded_warnings
