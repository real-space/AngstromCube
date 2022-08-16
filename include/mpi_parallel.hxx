#pragma once

/*
 * Where does A43 need MPI parallelization?
 * 
 * - Assume that eigenvector mode runs on a single CPU only
 * 
 * There are two levels of work: potential generation and KS-equation solving
 * - potential_generator. needs to deal with internal boundary conditions
 *   and works on a small domain of the entire supercell.
 *   For load balancing of systems with large vacuum regions that we want to
 *   sparsify, we could even think of mirco-domain, e.g. 8x8x8 grid points
 *   and let an MPI process deal with several of those.
 *   
 * - The Kohn-Sham equation can be solved in order(N) using
 *   a Green function or density matrix method.
 *   In case of the Green function, MPI communication is only needed
 *   to communicate the potential V(\vec r) and the atomic matrices
 *   (assuming that the atom centers and sigmas did not move/change).
 * 
 * A first implementation of the commmunication could be
 *    MPI_Alltoallv
 * where a different message length per process pair needs to be agreed on.
 * In order to scale, this commmunication pattern could be replaced
 * by MPI-one-sided communication. However, we still need to bookkeep,
 * where the data resides.
 * 
 */

#ifndef HAS_NO_MPI

  #include <mpi.h> // MPI_*

#else // HAS_NO_MPI

  // define MPI stubs

  #include <cassert> // assert
  #include <cstdint> // int64_t
  #include <cstring> // std::memcpy

  #include "recorded_warnings.hxx" // warn, error

  typedef int64_t MPI_Comm;
  MPI_Comm constexpr MPI_COMM_WORLD = -1;
  int constexpr MPI_SUCCESS = 0;

  typedef char MPI_Op;
  MPI_Op constexpr MPI_SUM = '+', MPI_MAX = 'M', MPI_MIN = 'm', MPI_OP_NULL = 0;//, MPI_PROD = '*';

  typedef int MPI_Datatype;
  inline size_t const size_of(MPI_Datatype const datatype) {
     switch (datatype) {
       case 0: return 1;
     }
     warn("unknown MPI_Datatype %d", datatype);
     return 0;
  } // size_of

  // This is a replacement if you do not have MPI installed
  #define ok   return MPI_SUCCESS
  int MPI_Init(int *argc, char ***argv) { ok; }
  int MPI_Finalize(void) { ok; }
  int MPI_Comm_rank(MPI_Comm comm, int *rank) { assert(rank); *rank = 0; ok; }
  int MPI_Comm_size(MPI_Comm comm, int *size) { assert(size); *size = 1; ok; }
  int MPI_Allreduce(void const *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
      std::memcpy(recvbuf, sendbuf, count*size_of(datatype)); ok; }
  int MPI_Barrier(MPI_Comm comm) { ok; }
  // add more MPI_ replacement functions here ...
  #undef  ok



#endif // HAS_NO_MPI

#include <cstdio> // std::printf

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

  int init(int argc=0, char **argv=nullptr) { // forward the arguments of main
      static bool already{false};
      if (already) return 1; // has already been initialized
      already = true;
      return MPI_Check(MPI_Init(&argc, &argv));
  } // init

  MPI_Comm comm() { return MPI_COMM_WORLD; }

  unsigned size(MPI_Comm const comm=MPI_COMM_WORLD) {
      int size{0};
      MPI_Check(MPI_Comm_size(comm, &size));
      assert( size > 0 );
      return size;
  } // size

  int rank(MPI_Comm const comm=MPI_COMM_WORLD, unsigned const check_size=0) { // check_size=0: do not check
      int rank{-1};
      MPI_Check(MPI_Comm_rank(comm, &rank));
      assert( rank >= 0 );
      if (check_size > 0) assert( rank < check_size );
      return rank;
  } // rank

  template <typename T> MPI_Datatype get();
//   template <> MPI_Datatype get<int8_t>  { return MPI_INTEGER1; }
//   template <> MPI_Datatype get<int16_t> { return MPI_INTEGER2; }
//   template <> MPI_Datatype get<int32_t> { return MPI_INTEGER4; }
//   template <> MPI_Datatype get<int64_t> { return MPI_INTEGER8; }

  template <typename T>
  int allreduce(T *recv, MPI_Op const op=MPI_SUM, MPI_Comm const comm=MPI_COMM_WORLD, int const count=1, T const *send=nullptr) {
      if (!send) send = recv; // need a deep copy here?
      return MPI_Allreduce(send, recv, count, get<T>(), op, comm);
  } // allreduce

  int barrier(MPI_Comm const comm=MPI_COMM_WORLD) { 
      return MPI_Check(MPI_Barrier(comm));
  } // barrier

  int finalize(void) {
      return MPI_Check(MPI_Finalize());
  } // finalize

} // namespace mpi_parallel


#include "status.hxx" // status_t

namespace mpi_parallel {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_simple(int const echo=0) {
      status_t stat(0);
//    int argc{0}; char **argv{nullptr};
//    MPI_Check(MPI_Init(&argc, &argv));
      MPI_Comm const comm{MPI_COMM_WORLD};
      int size{0};  MPI_Check(MPI_Comm_size(comm, &size));
      int rank{-1}; MPI_Check(MPI_Comm_rank(comm, &rank));
      MPI_Check(MPI_Barrier(comm));
      if (echo > 0) std::printf("# %s: MPI Comm_size= %d  rank= %i\n", __FILE__, size, rank);
//    MPI_Check(MPI_Finalize());

      stat += (size < 1); // error if less than 1 process
      stat += (rank < 0); // error if negative rank
      stat += (rank >= size); // error if rank out of range
      return stat;
  } // test_simple

  inline status_t test_wrappers(int const echo=0) {
      status_t stat(0);
//    stat += mpi_parallel::init();
      auto const comm = mpi_parallel::comm();
      auto const size = mpi_parallel::size(comm);
      auto const rank = mpi_parallel::rank(comm, size);
      stat += mpi_parallel::barrier(comm);
      if (echo > 0) std::printf("# %s: MPI Comm_size= %d  rank= %i\n", __FILE__, size, rank);
//    stat += mpi_parallel::finalize();

      stat += (size < 1); // error if less than 1 process
      stat += (rank < 0); // error if negative rank
      stat += (rank >= size); // error if rank out of range
      return stat;
  } // test_wrappers

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      bool const already_initialized = mpi_parallel::init();
      stat += test_wrappers(echo);
      stat += test_simple(echo);
      if (!already_initialized) mpi_parallel::finalize();
      return stat;
  } // all_tests

#endif

} // namespace mpi_parallel
