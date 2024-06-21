// This file is part of AngstromCube under MIT License

#include <cstdint> // int32_t

#include "define_version.hxx" // define_version
#include "control.hxx" // ::set

extern "C" {
    #include "single_atom.h"

    void live_atom_is_a_dynamic_library_(int32_t *is_dynamic) {

        if (nullptr != is_dynamic) *is_dynamic = 0; // static

        static bool set_version{true};
        if (set_version) {
            auto const version_key = define_version();
            control::set("version.latom", version_key);
            set_version = false; // only once
        } // set_version

    } // live_atom_is_a_dynamic_library_

} // extern "C"
