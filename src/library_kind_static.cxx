// This file is part of AngstromCube under MIT License

#include <cstdint> // int32_t

extern "C" {
    #include "single_atom.h"

    void live_atom_is_a_dynamic_library_(int32_t *status) {
        if (nullptr != status) *status = 0; // is static
    } // live_atom_is_a_dynamic_library_
} // extern "C"


