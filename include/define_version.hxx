#pragma once

    inline char const * define_version() {
#ifdef    _GIT_KEY
        // stringify the value of a macro, two expansion levels needed
        #define macro2string(a) stringify(a)
        #define stringify(b) #b
        auto const git_key = macro2string(_GIT_KEY);
        #undef  stringify
        #undef  macro2string
        return git_key;
#else  // _GIT_KEY
        return "";
#endif // _GIT_KEY
    } // define_git_key
