#ifdef    _GIT_KEY
        // stringify the value of a macro, two expansion levels needed
        #define macro2string(a) stringify(a)
        #define stringify(b) #b
        auto const version_key = macro2string(_GIT_KEY);
        #undef  stringify
        #undef  macro2string
#else  // _GIT_KEY
        auto const version_key = "";
#endif // _GIT_KEY
