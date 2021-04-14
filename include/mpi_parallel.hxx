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

  #include <cassert> // assert
  #include <cstdint> // int64_t

  typedef int64_t MPI_Comm;
  MPI_Comm constexpr MPI_COMM_WORLD = -1;
  int constexpr MPI_SUCCESS = 0;

  // This is a replacement if you do not have MPI installed
  #define ok   return MPI_SUCCESS
  int MPI_Init(int *argc, char ***argv) {
      static bool already{false};
      int const return_value = int(already);
      already = true;
      return return_value;
  } // MPI_Init
  int MPI_Finalize(void) { ok; };
  int MPI_Comm_rank(MPI_Comm comm, int *rank) { assert(rank); *rank = 0; ok; }
  int MPI_Comm_size(MPI_Comm comm, int *size) { assert(size); *size = 1; ok; }
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

  int initialize(int const argc=0, char *argv[]=nullptr) { // forward the arguments of main
      int argument_count{argv ? argc : 0};
      char **argument_values{(argc > 0) ? argv : nullptr};
      return MPI_Check(MPI_Init(&argument_count, &argument_values));
  } // initialize

  MPI_Comm comm() { return MPI_COMM_WORLD; }

  unsigned size(MPI_Comm comm=MPI_COMM_WORLD) {
      int size{0};
      MPI_Check(MPI_Comm_size(comm, &size));
      assert( size > 0 );
      return size;
  } // size

  int rank(MPI_Comm comm=MPI_COMM_WORLD, unsigned const size=0) { // size=0: do not check
      int rank{-1};
      MPI_Check(MPI_Comm_rank(comm, &rank));
      assert( rank >= 0 );
      if (size > 0) assert( rank < size );
      return rank;
  } // rank

  int barrier(MPI_Comm comm=MPI_COMM_WORLD) { 
      return MPI_Check(MPI_Barrier(comm));
  } // barrier
  
  int finalize() {
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
      int argc{0};
      char ** argv{nullptr};
      MPI_Check(MPI_Init(&argc, &argv));
      MPI_Comm const comm{MPI_COMM_WORLD};
      int size{0};  MPI_Check(MPI_Comm_size(comm, &size));
      int rank{-1}; MPI_Check(MPI_Comm_rank(comm, &rank));
      MPI_Check(MPI_Barrier(comm));
      if (echo > 0) std::printf("# %s: MPI Comm_size= %d  rank= %i\n", __FILE__, size, rank);
      MPI_Check(MPI_Finalize());

      stat += (size < 1); // error if less than 1 process
      stat += (rank < 0); // error if negative rank
      stat += (rank >= size); // error if rank out of range
      return stat;
  } // test_simple

  inline status_t test_wrappers(int const echo=0) {
      status_t stat(0);
      stat += mpi_parallel::initialize();
      auto const comm = mpi_parallel::comm();
      auto const size = mpi_parallel::size(comm);
      auto const rank = mpi_parallel::rank(comm, size);
      stat += mpi_parallel::barrier(comm);
      if (echo > 0) std::printf("# %s: MPI Comm_size= %d  rank= %i\n", __FILE__, size, rank);
      stat += mpi_parallel::finalize();

      stat += (size < 1); // error if less than 1 process
      stat += (rank < 0); // error if negative rank
      stat += (rank >= size); // error if rank out of range
      return stat;
  } // test_wrappers

  inline status_t all_tests(int const echo=0) { 
      status_t stat(0);
      stat += test_simple(echo);
      stat += test_wrappers(echo);
      return stat;
  } // all_tests

#endif

} // namespace mpi_parallel
