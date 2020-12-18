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
  #define MPI_(FunctionName)  MPI_##FunctionName

#else // HAS_NO_MPI

  #include <cassert> // assert
  #include <cstdint> // int64_t

  typedef int64_t MPI_Comm;
  MPI_Comm constexpr MPI_COMM_WORLD = -1;
  int constexpr MPI_SUCCESS = 0;

  #define MPI_(FunctionName)  mpi_parallel::FunctionName
namespace mpi_parallel {
  // This is a replacement if you do not have MPI installed
#define ok   return MPI_SUCCESS
  int Init(int *argc, char ***argv) { ok; };
  int Finalize(void) { ok; };
  int Comm_rank(MPI_Comm comm, int *rank) { *rank = 0; ok; }
  int Comm_size(MPI_Comm comm, int *size) { *size = 1; ok; }
#undef  ok

} // namespace mpi_parallel

#endif // HAS_NO_MPI

#include <cstdio> // printf
#include "status.hxx" // status_t

namespace mpi_parallel {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS
  inline status_t simple_test(int const echo=0) {
      status_t stat(0);
      int argc{0};
      char ** argv{nullptr};
      stat += MPI_(Init)(&argc, &argv);
      MPI_Comm const comm{MPI_COMM_WORLD};
      int size{0}, rank{-1};
      stat += MPI_(Comm_size)(comm, &size);
      stat += MPI_(Comm_rank)(comm, &rank);
      printf("# %s: MPI Comm_size= %d  rank= %i\n", __FILE__, size, rank);
      stat += MPI_(Finalize)();
      return stat;
  } // simple_test

  inline status_t all_tests(int const echo=0) { 
      status_t stat(0);
      stat += simple_test(echo);
      return stat;
  } // all_tests
#endif

} // namespace mpi_parallel
