#pragma once
// This file is part of AngstromCube under MIT License

#ifdef    DEBUG
    #include <cstdio> // std::printf
#endif // DEBUG

#include "status.hxx" // status_t

namespace green_memory {

    void* malloc(size_t const size_in_Bytes, char const *const name=""); // declaration only
    void free(void* ptr, char const *const name=""); // declaration only

    status_t all_tests(int const echo=0); // declaration only

} // namespace green_memory


    template <typename T>
    T* get_memory(size_t const size=1, int const echo=0, char const *const name="") {

        size_t const total = size*sizeof(T);
#ifdef    DEBUG
        if (echo > 0) {
            std::printf("# managed memory: %lu x %.3f kByte = \t", size, sizeof(T)*1e-3);
            if (total > 1e9) { std::printf("%.9f GByte\n", total*1e-9); } else 
            if (total > 1e6) { std::printf("%.6f MByte\n", total*1e-6); } else 
                             { std::printf("%.3f kByte\n", total*1e-3); }
        } // echo
#endif // DEBUG

        auto const ptr = green_memory::malloc(total, name);

#ifdef    DEBUGGPU
        std::printf("# get_memory \t%lu x %.3f kByte = \t%.3f kByte, %s at %p\n", size, sizeof(T)*1e-3, total*1e-3, name, (void*)ptr);
#endif // DEBUGGPU

        return (T*)ptr;
    } // get_memory


    template <typename T>
    void _free_memory(T* & ptr, char const *const name="") {
//      std::printf("# free_memory %s at %p\n", name, (void*)ptr);
        if (nullptr != ptr) {
#ifdef    DEBUGGPU
            std::printf("# free_memory %s at %p\n", name, (void*)ptr);
#endif // DEBUGGPU
            green_memory::free((void*)ptr, name);
        } // d
        ptr = nullptr;
    } // _free_memory

#define free_memory(PTR) _free_memory(PTR, #PTR)

    template <typename real_t=float>
    inline char const* real_t_name() { return (8 == sizeof(real_t)) ? "double" : "float"; }
