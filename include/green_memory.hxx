#pragma once
// This file is part of AngstromCube under MIT License

// #doc
// This offers an allocator for managed GPU memory.
//

#include <cstddef> // size_t
#ifdef    DEBUG
    #include <cstdio> // std::printf
#endif // DEBUG

#include "status.hxx" // status_t

namespace green_memory {

    // wrapper to allocate new managed GPU memory (or CPU memory if compiled without CUDA)
    void* malloc(size_t const size_in_Bytes, char const *const name=""); // declaration only

    // wrapper to free the managed GPU memory
    void free(void* ptr, char const *const name=""); // declaration only

    // show how much memory is in use now
    size_t total_memory_now(); // declaration only

    // show the maximum of memory that has been used since the program start
    size_t high_water_mark(); // declaration only

    // self-tests
    status_t all_tests(int const echo=0); // declaration only

} // namespace green_memory

    // useful wrapper around malloc to allocate arrays, e.g. get_memory<float[3][4]>
    template <typename T>
    T* get_memory(size_t const size=1, int const echo=0, char const *const name="") {

        size_t const total = size*sizeof(T);
#ifdef    DEBUG
        if (echo > 7) {
            double const f = (total > 1e9) ? 1e-9 : ((total > 1e6) ? 1e-6 : 1e-3);
            char  const _f = (total > 1e9) ? 'G'  : ((total > 1e6) ? 'M'  : 'k' );
            std::printf("# managed memory: %.3f k x %.3f kByte = \t%g %cByte\n", size*1e-3, sizeof(T)*1e-3, total*f, _f);
        } // echo
#endif // DEBUG

        auto const ptr = green_memory::malloc(total, name);

#ifdef    DEBUGGPU
        std::printf("# get_memory \t%lu x %.3f kByte = \t%.3f kByte, %s at %p\n", size, sizeof(T)*1e-3, total*1e-3, name, (void*)ptr);
#endif // DEBUGGPU

        return (T*)ptr;
    } // get_memory

    // corresponding wrapper for get_memory<T>, please use the macro free_memory (without leading underscore) to get a debug string
    template <typename T>
    void _free_memory(T* & ptr, char const *const name="") {
#ifdef    DEBUG
//      std::printf("# free_memory \'%s\' at %p\n", name, (void*)ptr);
#endif // DEBUG
        if (nullptr != ptr) {
#ifdef    DEBUGGPU
            std::printf("# free_memory \'%s\' at %p\n", name, (void*)ptr);
#endif // DEBUGGPU
            green_memory::free((void*)ptr, name);
        } // d
        ptr = nullptr;
    } // _free_memory

#define free_memory(PTR) _free_memory(PTR, #PTR)

    // display name according to floating.point precision
    template <typename real_t=float>
    inline char const* real_t_name() { return (8 == sizeof(real_t)) ? "double" : "float"; }
