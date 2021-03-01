#pragma once

  template <typename T>
  T* get_memory(size_t const size=1) {
#ifdef DEBUG
      double const GByte = 1e-9; 
      char const *const _GByte = "GByte";
      std::printf("# managed memory: %lu x %.3f kByte = \t%.6f %s\n", size, 1e-3*sizeof(T), size*sizeof(T)*GByte, _GByte);
#endif
      T* d{nullptr};
#ifdef HAS_CUDA
      CCheck(cudaMallocManaged(&d, size*sizeof(T)));
#else
      d = new T[size];
#endif
      return d;
  } // get_memory

  template <typename T>
  void free_memory(T* &d) {
      if (nullptr != d) {
#ifdef HAS_CUDA
          CCheck(cudaFree(d));
#else
          delete[] d;
#endif
      } // d
      d = nullptr;
  } // free_memory
