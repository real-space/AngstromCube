
/*
 * Where does A43 need load balancing?
 * 
 * 
 */

#include <cassert> // assert
#include <cstdint> // int64_t, size_t, uint8_t
#include <cstdio> // std::printf
#include <algorithm> // std::max
#include <vector> // std::vector<T>
#include <cmath> // std::floor, ::ceil

#include "load_balancer.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
// #include "mpi_parallel.hxx" // ...
#include "simple_stats.hxx" // ::Stats<>
#include "data_view.hxx" // view4D<T>
#include "inline_math.hxx" // pow2
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

namespace load_balancer {

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

  
  inline status_t test_diffusion_balancer(int const echo=0) {
      status_t stat(0);
      int constexpr X=0, Y=1, Z=2;

      // iterative process for balacing the computational load onto MPI processes.
      // The load can be the number and costs of compute tasks, memory consumption, etc.

      int const n[] = {int(control::get("load_balancer.test.nx", 18.)), // number of blocks
                       int(control::get("load_balancer.test.ny", 19.)),
                       int(control::get("load_balancer.test.nz", 20.))};
      int const nprocs = control::get("load_balancer.test.nprocs", 53.);
      if (echo > 0) std::printf("\n\n# %s start %d x %d x %d with %d MPI processes\n", __func__, n[X], n[Y], n[Z], nprocs);

      assert(nprocs > 0);
      auto const blockvolume = n[X]*size_t(n[Y])*size_t(n[Z]);
      auto const blocksperproc = blockvolume/double(nprocs);
      auto const cube_root = std::cbrt(blocksperproc);
      if (echo > 0) std::printf("# division of %d x %d x %d = %ld blocks onto %d processes, expected %.2f^3 = %.2f blocks per process\n",
                                  n[X], n[Y], n[Z], blockvolume, nprocs, cube_root, blocksperproc);
      assert(blockvolume > 0);

      // now approximatelty factorize nprocs
      auto const inv_cbrt = 1./cube_root;
      int npd[3];
      for (int d = 0; d < 3; ++d) {
          auto const npdf = n[d]*inv_cbrt;
          npd[d] = std::ceil(npdf);
          if (echo > 2) std::printf("# in %c-direction [%d, %.2f, %d] blocks per process\n", 'x' + d, npd[d] - 1, npdf, npd[d]);
          assert(npd[d] > 0);
      } // d
      int const nprocs_more = npd[X] * npd[Y] * npd[Z];
      assert(nprocs_more > 0);
      float const by_npd[] = {n[X]/float(npd[X]), n[Y]/float(npd[Y]), n[Z]/float(npd[Z])};

//    float const decay2 = pow2(cube_root); // set decay such that the distribution is clearly nonzero on the positions of direct neighbors

      if (echo > 0) std::printf("# division into %d x %d x %d = %d processes is closest to %d processes\n",
                                        npd[X], npd[Y], npd[Z], nprocs_more,             nprocs);

      // skip positions that are too much
      view4D<float> pop(nprocs, n[Z], n[Y], n[X], 0.f);

for (int myrank = 0; myrank < nprocs; ++myrank) {
  
      int const iproc_more = (myrank*nprocs_more)/nprocs;     assert(iproc_more < nprocs_more);
      int const coords[] = {iproc_more % npd[X], (iproc_more / npd[X]) % npd[Y], iproc_more / (npd[X]*npd[Y])};

      if (echo > 9) std::printf("# myrank=%i pseudorank=%i of %d, coords %d %d %d\n", myrank, iproc_more, nprocs_more, coords[X], coords[Y], coords[Z]);
      // place centroid on the map
      float const position[] = {(coords[X] + .5f)*by_npd[X],
                                (coords[Y] + .5f)*by_npd[Y],
                                (coords[Z] + .5f)*by_npd[Z]};

      auto pop_myrank = pop[myrank];

      for (int iz = 0; iz < n[Z]; ++iz) {
      for (int iy = 0; iy < n[Y]; ++iy) {
      for (int ix = 0; ix < n[X]; ++ix) {
          auto const dist2 = pow2(ix - position[X])
                           + pow2(iy - position[Y])
                           + pow2(iz - position[Z]);
//        pop_myrank(iz,iy,ix) = 250.f/pow2(dist2 + decay2);
          pop_myrank(iz,iy,ix) = -dist2;
      }}} // ix iy iz

} // myrank

      view3D<int32_t> owner_rank(n[Z], n[Y], n[X], -1);
      std::vector<uint32_t> n_own(nprocs, 0);
      for (int iz = 0; iz < n[Z]; ++iz) {
      for (int iy = 0; iy < n[Y]; ++iy) {
      for (int ix = 0; ix < n[X]; ++ix) {
          float largest{-3e33f};
          int l_rank{-1};
          for (int myrank = 0; myrank < nprocs; ++myrank) {
              auto const p = pop(myrank,iz,iy,ix);
              if (p > largest) { largest = p; l_rank = myrank; }
          } // myrank
          assert(l_rank >= 0);
          owner_rank(iz,iy,ix) = l_rank;
          ++n_own[l_rank];
      }}} // ix iy iz

      simple_stats::Stats<> st(0);
      for (int myrank = 0; myrank < nprocs; ++myrank) {
          st.add(n_own[myrank]);
          if (echo > 9) std::printf("# myrank=%i owns %d\n", myrank, n_own[myrank]);
      } // myrank
      if (echo > 0) std::printf("# %ld processes own %g blocks, per process [%g, %.2f +/- %.2f, %g]\n",
                                    st.tim(), st.sum(), st.min(), st.mean(), st.dev(), st.max());

      auto const load_balance = st.mean()/st.max(), load_imbalance = st.max()/st.min();
      if (echo > 0) std::printf("# load balance %.2f %%, imbalance %.2f\n", load_balance*100, load_imbalance);

      return stat;
  } // test_diffusion_balancer


  status_t all_tests(int const echo) { 
      status_t stat(0);
//    stat += test_diffusion_balancer(echo);
      stat += test_bisection_balancer(echo);
      return stat;
  } // all_tests

} // namespace load_balancer

#endif // NO_UNIT_TESTS
