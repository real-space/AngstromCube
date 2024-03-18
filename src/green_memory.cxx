// This file is part of AngstromCube under MIT License

#include <cstdlib> // std::size_t
#include <cstdio> // std::printf

#include "green_memory.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_cuda.hxx" // cuCheck, cudaMallocManaged, cudaFree


namespace green_memory {

    void* malloc(std::size_t const size_in_Bytes, char const *const name) {
#ifndef HAS_NO_CUDA
        char* ptr{nullptr};
        cuCheck(cudaMallocManaged(&ptr, size_in_Bytes));
#else  // HAS_NO_CUDA
        auto const *const ptr = new char[size_in_Bytes];
#endif // HAS_NO_CUDA

#ifdef    DEBUGGPU
        std::printf("# green_memory::malloc %.3f kByte, %s at %p\n", size_in_Bytes*1e-3, name, (void*)ptr);
#endif // DEBUGGPU

        return (void*)ptr;
    } // malloc

    void free(void* ptr, char const *const name) {
        if (ptr) {
#ifdef    DEBUGGPU
          std::printf("# green_memory::free %s at %p\n", name, ptr);
#endif // DEBUGGPU

#ifndef HAS_NO_CUDA
          cuCheck(cudaFree((void*)ptr));
#else  // HAS_NO_CUDA
          delete[] (char*)ptr;
#endif // HAS_NO_CUDA
        } else {
            std::printf("# green_memory::free %s got nullptr\n", name);
        }
    } // free


#ifdef  NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

    status_t test_green_memory(int const echo=0) {
        status_t stat(0);
        auto const nrand = (std::rand()*365)/RAND_MAX;
        auto const mem0 = green_memory::malloc(nrand, "mem0");
        green_memory::free(mem0, "mem0");
        return stat;
    } // test_green_memory

    status_t all_tests(int const echo) {
        status_t stat(0);
        stat += test_green_memory(echo);
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_memory
