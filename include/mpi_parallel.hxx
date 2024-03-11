#pragma once
// This file is part of AngstromCube under MIT License

#ifndef   HAS_NO_MPI

  #include <mpi.h> // MPI_*

#else  // HAS_NO_MPI

  // define MPI stubs

  #include <cassert> // assert
  #include <cstdint> // int64_t
  #include <cstring> // std::memcpy

  #include "recorded_warnings.hxx" // warn, error

  typedef int64_t MPI_Comm;
  MPI_Comm constexpr MPI_COMM_WORLD = -1, MPI_COMM_NULL = 0;
  int constexpr MPI_SUCCESS = 0;
  typedef int64_t MPI_Request;
  MPI_Request constexpr MPI_REQUEST_NULL = 0;

  typedef char MPI_Op;
  MPI_Op constexpr MPI_SUM = '+', MPI_MAX = 'M', MPI_MIN = 'm', MPI_OP_NULL = 0;//, MPI_PROD = '*';

  typedef int MPI_Datatype;
    MPI_Datatype constexpr MPI_DOUBLE   = 53,
                           MPI_FLOAT    = 24,
                           MPI_INT8_T   =  7,
                           MPI_UINT8_T  =  8,
                           MPI_INT16_T  = 15,
                           MPI_UINT16_T = 16,
                           MPI_INT32_T  = 31,
                           MPI_UINT32_T = 32,
                           MPI_INT64_T  = 63,
                           MPI_UINT64_T = 64;
  void const * const MPI_IN_PLACE = nullptr;

  inline size_t const size_of(MPI_Datatype const datatype) {
     switch (datatype) {
       case 0:                                      return 1; // 1 Byte
       case MPI_INT8_T : case MPI_UINT8_T :         return 1; // 1 Byte
       case MPI_INT16_T: case MPI_UINT16_T:         return 2; // 2 Byte
       case MPI_INT32_T: case MPI_UINT32_T:         return 4; // 4 Byte
       case MPI_INT64_T: case MPI_UINT64_T:         return 8; // 8 Byte
       case MPI_FLOAT:                  return sizeof(float); // 4 Byte
       case MPI_DOUBLE:                return sizeof(double); // 8 Byte
     }
     warn("unknown MPI_Datatype %d", int(datatype));
     return 0;
  } // size_of

  // This is a replacement if you do not have MPI installed
  #define ok   return MPI_SUCCESS
  inline int MPI_Init(int *argc, char ***argv) { ok; }
  inline int MPI_Finalize(void) { ok; }
  inline int MPI_Comm_rank(MPI_Comm comm, int *rank) { assert(rank); *rank = 0; ok; }
  inline int MPI_Comm_size(MPI_Comm comm, int *size) { assert(size); *size = 1; ok; }
  inline int MPI_Allreduce(void const *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
      if (sendbuf) { std::memcpy(recvbuf, sendbuf, count*size_of(datatype)); } ok; }
  inline int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) { ok; }
  inline int MPI_Barrier(MPI_Comm comm) { ok; }
  inline double MPI_Wtime(void) { return 0; } // ToDo
  // add more MPI_ replacement functions here ...
  #undef  ok

#endif // HAS_NO_MPI

#include <cstdio> // std::printf
#include "simple_stats.hxx" // ::Stats<double>
#include "status.hxx" // status_t

namespace mpi_parallel {

  #define MPI_Check(MPI_Function) \
    mpi_parallel::__check_MPI_call((MPI_Function), __FILE__, __LINE__, #MPI_Function)

  inline int __check_MPI_call( // --> always envoke via the macro MPI_Check
        int const MPI_status
      , char const *const file
      , unsigned const line
      , char const* const name
  ) {
#ifdef FULLDEBUG
      std::printf("# calling %s in %s:%d returned status= %i\n", name, file, line, MPI_status);
#endif
      if (MPI_SUCCESS != MPI_status) {
          std::printf("\n# in %s:%d failed with status= %i calling %s\n", file, line, MPI_status, name);
          // add here how to react to 
      } // failure
      return MPI_status; // for futher use
  } // __check_MPI_call --> always envoke via the macro MPI_Check



  // MPI utilities

  inline int init(int argc=0, char **argv=nullptr) { // forward the arguments of main
      static bool already{false};
      if (already) return 1; // has already been initialized
      already = true;
      return MPI_Check(MPI_Init(&argc, &argv));
  } // init

  inline MPI_Comm comm() { return MPI_COMM_WORLD; }

  inline unsigned size(MPI_Comm const comm=MPI_COMM_WORLD) {
      int size{0};
      MPI_Check(MPI_Comm_size(comm, &size));
      assert( size > 0 );
      return size;
  } // size

  inline int rank(MPI_Comm const comm=MPI_COMM_WORLD, unsigned const check_size=0) { // check_size=0: do not check
      int rank{-1};
      MPI_Check(MPI_Comm_rank(comm, &rank));
      assert( rank >= 0 );
      if (check_size > 0) assert( rank < check_size );
      return rank;
  } // rank

  template <typename T> MPI_Datatype get(T t=0);
  template <> inline MPI_Datatype get<char>    (char     t) { return MPI_UINT8_T;  }
  template <> inline MPI_Datatype get<int8_t>  (int8_t   t) { return MPI_INT8_T;   }
  template <> inline MPI_Datatype get<uint8_t> (uint8_t  t) { return MPI_UINT8_T;  }
  template <> inline MPI_Datatype get<int16_t> (int16_t  t) { return MPI_INT16_T;  }
  template <> inline MPI_Datatype get<uint16_t>(uint16_t t) { return MPI_UINT16_T; }
  template <> inline MPI_Datatype get<int32_t> (int32_t  t) { return MPI_INT32_T;  }
  template <> inline MPI_Datatype get<uint32_t>(uint32_t t) { return MPI_UINT32_T; }
  template <> inline MPI_Datatype get<int64_t> (int64_t  t) { return MPI_INT64_T;  }
  template <> inline MPI_Datatype get<size_t>  (size_t   t) { return MPI_UINT64_T; }
  template <> inline MPI_Datatype get<float>   (float    t) { return MPI_FLOAT;    }
  template <> inline MPI_Datatype get<double>  (double   t) { return MPI_DOUBLE;   }

  template <typename T>
  inline int allreduce(T *recv, MPI_Op const op=MPI_SUM, MPI_Comm const comm=MPI_COMM_WORLD, size_t const count=1, T const *send=nullptr) {
      return MPI_Allreduce(send ? send : MPI_IN_PLACE, recv, count, get<T>(), op, comm);
  } // allreduce

  template <typename T>
  inline int broadcast(T *buffer, MPI_Comm const comm=MPI_COMM_WORLD, size_t const count=1, int const root=0) {
      return MPI_Bcast(buffer, count, get<T>(), root, comm);
  } // broadcast

  template <typename T>
  inline int max(T *recv, size_t const count=1, MPI_Comm const comm=MPI_COMM_WORLD) {
      return allreduce(recv, MPI_MAX, comm, count);
  } // max

  template <typename T>
  inline int min(T *recv, size_t const count=1, MPI_Comm const comm=MPI_COMM_WORLD) {
      return allreduce(recv, MPI_MIN, comm, count);
  } // min

  template <typename T>
  inline int sum(T *recv, size_t const count=1, MPI_Comm const comm=MPI_COMM_WORLD) {
      return allreduce(recv, MPI_SUM, comm, count);
  } // sum

  template <typename T>
  inline T max(T const in, MPI_Comm const comm=MPI_COMM_WORLD) {
      T out{in}; allreduce(&out, MPI_MAX, comm, 1); return out;
  } // max (scalars)

  template <typename T>
  inline T min(T const in, MPI_Comm const comm=MPI_COMM_WORLD) {
      T out{in}; allreduce(&out, MPI_MIN, comm, 1); return out;
  } // min (scalars)

  template <typename T>
  inline T sum(T const in, MPI_Comm const comm=MPI_COMM_WORLD) {
      T out{in}; allreduce(&out, MPI_SUM, comm, 1); return out;
  } // sum (scalars)

  inline int barrier(MPI_Comm const comm=MPI_COMM_WORLD) { 
      return MPI_Check(MPI_Barrier(comm));
  } // barrier

  inline int finalize(void) {
      return MPI_Check(MPI_Finalize());
  } // finalize


  inline int allreduce(simple_stats::Stats<double> & stats, MPI_Comm const comm=MPI_COMM_WORLD) {
      double v[8];
      stats.get(v);
      auto const status_sum = sum(v, 5, comm); // MPI_SUM on {v[0], v[1], v[2], v[3], v[4]}
      auto const status_max = max(v + 6, 2, comm); // MPI_MAX on {v[6], v[7]}
      stats.set(v);
      return status_sum + status_max;
  } // allreduce


  status_t all_tests(int const echo=0); // declaration only

} // namespace mpi_parallel
