#pragma once

/*
 * Where does A43 need load balancing?
 * 
 * 
 */

#include <cassert> // assert
#include <cstdint> // int64_t, size_t
#include <cstdio> // std::printf
#include <algorithm> // std::max
#include <vector> // std::vector<T>
#include <cmath> // std::floor, ::ceil

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
// #include "mpi_parallel.hxx" // ...
#include "simple_stats.hxx" // ::Stats<>
#include "recorded_warnings.hxx" // error

namespace load_balancer {

  template <typename number_t, int sgn=1>
  int largest_of_3(number_t const n[3]) {
      if (n[0]*sgn >= n[1]*sgn) { // answer is 0 or 2
          return (n[0]*sgn >= n[2]*sgn) ? 0 : 2;
      } else {               // answer is 1 or 2
          return (n[1]*sgn >= n[2]*sgn) ? 1 : 2;
      }
  } // largest_of_3

  status_t bisection_balancer(
        int nd[3] // result and input for the number of blocks to be distributed
      , int & level
      , int const myrank=0
      , int const nprocs=1
      , int const echo=0 // log-level
  ) {
      int jproc{myrank};
      int np{nprocs};
      level = 0;
      while (nd[0]*size_t(nd[1])*size_t(nd[2]) > 1 && np > 1) {
          // divide the number of processes into potentially unequal halfs
                                                    assert(np > 1);
          int const right_np = np/2;                assert(right_np > 0);
          int const left_np = np - right_np;        assert(left_np > 0);
                                                    assert(left_np + right_np == np);
          // divide along direction d
          int const d = largest_of_3(nd);           assert(nd[d] > 1);
          int const right_nd = std::max(1, (nd[d] * right_np)/np); // flooring
                                                    assert(right_nd > 0);
          int const left_nd = nd[d] - right_nd;     assert(left_nd > 0);
                                                    assert(left_nd + right_nd == nd[d]);
          int const is_right = (jproc >= left_np);
          // apply division
          jproc -= is_right * left_np;              assert(jproc >= 0);
          nd[d] = is_right ? right_nd : left_nd;    assert(nd[d] > 0);
          np    = is_right ? right_np : left_np;    assert(jproc < np);
          ++level;
      } // while
      if (echo > 99) std::printf("# myrank=%i grid %d %d %d, level= %d\n", myrank, nd[0], nd[1], nd[2], level);
      
      return np - 1;
  } // bisection_balancer
  
} // namespace load_balancer


#ifndef  NO_UNIT_TESTS
    #include "control.hxx" // ::get
#endif // NO_UNIT_TESTS

namespace load_balancer {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_bisection_balancer(int const echo=0) {
      status_t stat(0);
      int constexpr X=0, Y=1, Z=2;

      // iterative process for balacing the computational load onto MPI processes.
      // The load can be the number and costs of compute tasks, memory consumption, etc.

      int const n[] = {int(control::get("load_balancer.test.nx", 18.)), // number of blocks
                       int(control::get("load_balancer.test.ny", 19.)),
                       int(control::get("load_balancer.test.nz", 20.))};
      int const nprocs = control::get("load_balancer.test.nprocs", 53.);
      if (echo > 0) std::printf("\n\n# %s start %d x %d x %d with %d MPI processes\n", __func__, n[X], n[Y], n[Z], nprocs);

      if (echo > 0) std::printf("# division of %d x %d x %d = %d blocks onto %d processes, expected %.2f blocks per process\n",
                                  n[X], n[Y], n[Z], n[X]*n[Y]*n[Z], nprocs, n[X]*double(n[Y])*double(n[Z])/nprocs);

      simple_stats::Stats<> level_stats(0), grid_stats[4]; // statistics over MPI processes

for (int myrank = 0; myrank < nprocs; ++myrank) {

      // recursively take the longest dimension,
      // divide it by "half" and distribute the remaining processes as equal as possible

      int level{0};
      int nd[] = {n[X], n[Y], n[Z]};
      stat += bisection_balancer(nd, level, myrank, nprocs, echo);
      for (int d = 0; d < 3; ++d) grid_stats[d].add(nd[d]);
      grid_stats[3].add(nd[0]*size_t(nd[1])*size_t(nd[2])); // volume
      level_stats.add(level);

} // myrank

      if (echo > 0) std::printf("# division into [%g, %.2f +/- %.2f, %g] levels\n",
                          level_stats.min(), level_stats.mean(), level_stats.dev(), level_stats.max());
      if (echo > 0) std::printf("# 2^level = [%g, %g, %g]\n",
                          std::pow(2., level_stats.min()), std::pow(2., level_stats.mean()), std::pow(2., level_stats.max()));
      for (int d = 0; d < 3; ++d) {
          if (echo > 5) std::printf("# grid in %c-direction has [%g, %.2f +/- %.2f, %g] blocks\n", 'x' + d, 
                  grid_stats[d].min(), grid_stats[d].mean(), grid_stats[d].dev(), grid_stats[d].max());
      } // d
      if (echo > 4) std::printf("# grid volumes in [%g, %.2f +/- %.2f, %g] blocks\n",
                grid_stats[3].min(), grid_stats[3].mean(), grid_stats[3].dev(), grid_stats[3].max());
      auto const load_balance = grid_stats[3].mean()/grid_stats[3].max();
      auto const load_imbalance = grid_stats[3].max()/grid_stats[3].min();
      if (echo > 0) std::printf("# load balance %.2f %%, imbalance %.2f\n", load_balance*100, load_imbalance);

      return stat;
  } // test_bisection_balancer

  inline status_t all_tests(int const echo=0) { 
      status_t stat(0);
      stat += test_bisection_balancer(echo);
      return stat;
  } // all_tests

#endif

} // namespace load_balancer
