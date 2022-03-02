#pragma once

#ifdef HAS_NO_CUDA
  #define __global__
  #define __restrict__
  #define __device__
  #define __shared__
  #define __unroll__
  #define __host__
#endif // HAS_NO_CUDA

#ifdef HAS_NO_CUDA
  struct dim3 {
      int x, y, z;
      dim3(int xx, int yy=1, int zz=1) : x(xx), y(yy), z(zz) {}
  }; // dim3

  inline void __syncthreads(void) {} // dummy
#endif // HAS_NO_CUDA      

#define DEBUG

  template <typename T>
  T* get_memory(size_t const size=1, int const echo=0, char const *const name="") {

#ifdef DEBUG
      if (echo > 0) {
          size_t const total = size*sizeof(T);
          std::printf("# managed memory: %lu x %.3f kByte = \t", size, sizeof(T)*1e-3);
          if (total > 1e9) { std::printf("%.9f GByte\n", total*1e-9); } else 
          if (total > 1e6) { std::printf("%.6f MByte\n", total*1e-6); } else 
                           { std::printf("%.3f kByte\n", total*1e-3); }
      } // echo
#endif // DEBUG

      T* d{nullptr};
#ifdef  HAS_CUDA
      CCheck(cudaMallocManaged(&d, size*sizeof(T)));
#else  // HAS_CUDA
      d = new T[size];
#endif // HAS_CUDA

#if 1 // def DEBUG
      std::printf("# get_memory \t%lu x %.3f kByte = \t%.3f kByte, %s at %p\n", size, sizeof(T)*1e-3, size*sizeof(T)*1e-3, name, (void*)d);
#endif // DEBUG

      return d;
  } // get_memory


  template <typename T>
  void _free_memory(T* &d, char const *const name="") {
      if (nullptr != d) {
#if 1 // def DEBUG
          std::printf("# free_memory %s at %p\n", name, (void*)d);
#endif // DEBUG

#ifdef  HAS_CUDA
          CCheck(cudaFree(d));
#else  // HAS_CUDA
          delete[] d;
#endif // HAS_CUDA
      } // d
      d = nullptr;
  } // free_memory

#define free_memory(PTR) _free_memory(PTR, #PTR)
  
  template <typename real_t=float>
  inline char const* real_t_name() { return (8 == sizeof(real_t)) ? "double" : "float"; } 

  //
  // Memory layout for Green function and atomic projection coefficients
  //
  //  version G(*)[R1C2][Noco*64][Noco*64]
  //          a(*)[R1C2][Noco   ][Noco*64]
  //
  //            kinetic:    <<< {16, Nrows, 1}, {Noco*64, Noco, R1C2} >>>
  //            add:        <<< {nrhs, ncubes, 1}, {Noco*64, 1, 1} >>>
  //            prj:        <<< {nrhs, natoms, 1}, {Noco*64, 1, 1} >>>
  //            potential:  <<< {64, any, 1}, {Noco*64, Noco, R1C2} >>>
  //  --> if 2 == R1C2, tfQMRgpu with LM=Noco*64
  //
