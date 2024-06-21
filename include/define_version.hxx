#pragma once

#include <cstdio> // std::printf
#include "control.hxx" // ::set

    inline void define_version(char const *const keyword, int const echo=0) {
#ifdef    _GIT_KEY
      // stringify the value of a macro, two expansion levels needed
      #define macro2string(a) stringify(a)
      #define stringify(b) #b
      auto const git_key = macro2string(_GIT_KEY);
      #undef  stringify
      #undef  macro2string
      control::set(keyword, git_key); // store in the global variable environment
      if (echo > 0) std::printf("# %s: git checkout %s\n\n", keyword, git_key);
#endif // _GIT_KEY
    } // define_git_key
