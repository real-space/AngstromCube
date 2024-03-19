#pragma once
// This file is part of AngstromCube under MIT License

#ifndef   HAS_NO_CUDA

  #include <cuda.h> // dim3, cudaStream_t, __syncthreads, cuda*

  __host__ inline
  void __cudaSafeCall(cudaError_t err, char const *file, int const line, char const *call) { // Actual check function
      if (cudaSuccess != err) {
          auto const msg = cudaGetErrorString(err);
          std::fprintf(stdout, "[ERROR] CUDA call to %s at %s:%d\n%s\n", call, file, line, msg); std::fflush(stdout);
          std::fprintf(stderr, "[ERROR] CUDA call to %s at %s:%d\n%s\n", call, file, line, msg); std::fflush(stderr);
          std::exit(line);
      }
  } // __cudaSafeCall
  #define cuCheck(err) __cudaSafeCall((err), __FILE__, __LINE__, #err) // Syntactic sugar to enhance output

#else  // HAS_NO_CUDA

  // replace CUDA specifics
  #define __global__
  #define __restrict__
  #define __device__
  #define __shared__
  #define __unroll__
  #define __host__

  struct dim3 {
      int x, y, z;
      dim3(int xx, int yy=1, int zz=1) : x(xx), y(yy), z(zz) {}
  }; // dim3

#ifndef   HAS_TFQMRGPU
    inline void __syncthreads(void) {} // dummy
    typedef int cudaError_t;
    inline cudaError_t cudaDeviceSynchronize(void) { return 0; } // dummy
    inline cudaError_t cudaPeekAtLastError(void) { return 0; } // dummy
#else  // HAS_TFQMRGPU
    #define gpuStream_t cudaStream_t
    #include "tfqmrgpu_cudaStubs.hxx" // cuda... (dummies)
#endif // HAS_TFQMRGPU

  #define cuCheck(err) ;

#endif // HAS_NO_CUDA
