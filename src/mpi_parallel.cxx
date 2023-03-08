// This file is part of AngstromCube under MIT License
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

  #include <cassert> // assert
  #include <cstdint> // int64_t
  #include <cstring> // std::memcpy

  #include "mpi_parallel.hxx"

  #include "recorded_warnings.hxx" // warn, error


namespace mpi_parallel {

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_simple(int const echo=0) {
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

  status_t test_wrappers(int const echo=0) {
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

  template <typename T=uint16_t>
  status_t test_constants(int const echo=0) {
      // check that get<T>() produces the right constant, in particular MPI_UINT16 is correct for uint16_t
      T a[] = {0, 1, 2, 3};
      auto const np = mpi_parallel::size();
      if (np > (1u << 14)) {
          if (echo > 1) std::printf("# %s %s cannot test with more than 2^14 processes, found %d\n", __FILE__, __func__, np);
          return 1;
      }
      MPI_Allreduce(MPI_IN_PLACE, a, 4, get<T>(), MPI_SUM, MPI_COMM_WORLD);
      return (0 != a[0]) + (np != a[1]) + (2*np != a[2]) + (3*np != a[3]);
  } // test_constants

  status_t all_tests(int const echo) {
      status_t stat(0);
      bool const already_initialized = mpi_parallel::init();
      stat += test_wrappers(echo);
      stat += test_simple(echo);
      stat += test_constants(echo);
      stat += test_constants<double>(echo);
      stat += test_constants<size_t>(echo);
      if (!already_initialized) mpi_parallel::finalize();
      return stat;
  } // all_tests

#endif

} // namespace mpi_parallel
