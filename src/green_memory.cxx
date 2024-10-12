// This file is part of AngstromCube under MIT License

#define   HAS_MEMORY_COUNTER

#include <cstdlib> // std::size_t
#include <cstdio> // std::printf
#ifdef    HAS_MEMORY_COUNTER
    #include <map> // std::map
    #include <algorithm> // std::max 
#endif // HAS_MEMORY_COUNTER

#include "green_memory.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_cuda.hxx" // cuCheck, cudaMallocManaged, cudaFree


namespace green_memory {

#ifdef    HAS_MEMORY_COUNTER
    size_t _memory_counter{0}; // global variable
    size_t _memory_maximum{0}; // global variable
    std::map<void*,size_t> _memory_map;
#endif // HAS_MEMORY_COUNTER

    void* malloc(std::size_t const size_in_Bytes, char const *const name) {
#ifndef   HAS_NO_CUDA
        void* ptr{nullptr};
        cuCheck(cudaMallocManaged(&ptr, size_in_Bytes));
#else  // HAS_NO_CUDA
        void *const ptr = new char[size_in_Bytes];
#endif // HAS_NO_CUDA

#ifdef    HAS_MEMORY_COUNTER
        _memory_counter += size_in_Bytes;
        _memory_maximum = std::max(_memory_maximum, _memory_counter);
        _memory_map[ptr] = size_in_Bytes;
#endif // HAS_MEMORY_COUNTER

#ifdef    DEBUGGPU
        std::printf("# green_memory::malloc %.3f kByte, \'%s\' at %p\n", size_in_Bytes*1e-3, name, ptr);
#endif // DEBUGGPU

        return ptr;
    } // malloc

    void free(void* ptr, char const *const name) {
        if (ptr) {
#ifdef    DEBUGGPU
          std::printf("# green_memory::free \'%s\' at %p\n", name, ptr);
#endif // DEBUGGPU

#ifdef    HAS_MEMORY_COUNTER
        auto const it = _memory_map.find(ptr);
        if (it != _memory_map.end()) {
            auto const size_in_Bytes = it->second;
            _memory_counter -= size_in_Bytes;
            _memory_map.erase(it);
        } // found
#endif // HAS_MEMORY_COUNTER

#ifndef   HAS_NO_CUDA
          cuCheck(cudaFree((void*)ptr));
#else  // HAS_NO_CUDA
          delete[] (char*)ptr;
#endif // HAS_NO_CUDA
        } else {
            std::printf("# green_memory::free \'%s\' got nullptr\n", name);
        }
    } // free


    size_t total_memory_now() {
#ifdef    HAS_MEMORY_COUNTER
        return _memory_counter;
#endif // HAS_MEMORY_COUNTER
        return 0;
    } // total_memory_now

    size_t high_water_mark() {
#ifdef    HAS_MEMORY_COUNTER
        return _memory_maximum;
#endif // HAS_MEMORY_COUNTER
        return 0;
    } // high_water_mark


#ifdef    NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

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
