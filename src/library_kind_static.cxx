// This file is part of AngstromCube under MIT License

#include <cstdint> // int32_t

// #include "control.hxx" // ::set, ::echo_set_without_warning
#include "define_version.hxx" // define_version

extern "C" {
    #include "single_atom.h"

    void live_atom_is_a_dynamic_library_(int32_t *is_dynamic) {
        static bool set_version{true};
        if (set_version) {
// #ifdef    _GIT_KEY
//             // stringify the value of a macro, two expansion levels needed
//             #define macro2string(a) stringify(a)
//             #define stringify(b) #b
//             auto const git_key = macro2string(_GIT_KEY);
//             #undef  stringify
//             #undef  macro2string
//             control::set("git.key.atom", git_key); // store in the global variable environment
// #endif // _GIT_KEY
            define_version("version.latom", 0);
            set_version = false; // only once
        } // set_version
        if (nullptr != is_dynamic) *is_dynamic = 0; // static
    } // live_atom_is_a_dynamic_library_
} // extern "C"
