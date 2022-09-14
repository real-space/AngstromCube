#pragma once

#include <cstdio> // std::printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <numeric> // std::iota
#include <cstdint> // uint32_t

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "data_view.hxx" // view2D<T>
#include "inline_math.hxx" // set, pow2, scale


namespace load_balancer {

  int constexpr X=0, Y=1, Z=2, W=3;

  template <typename real_t>
  inline double center_of_weight(
        double cow[4] // result: center of weight [0/1/2] and number of contributors [3]
      , size_t const nuna // number of unassigned work items
      , uint32_t const indirect[] // list of unassigned work items
      , view2D<real_t> const & xyzw // positions of work items in space [0/1/2] and their weight [3]
  ) {
      set(cow, 4, 0.0); // initialize
      double w8sum{0};
      for (size_t iuna = 0; iuna < nuna; ++iuna) { // parallel, reduction(+:w8sum,cow)
          auto const iall = indirect[iuna];
          auto const *const xyz = xyzw[iall];
          double const w8 = xyz[W];
          w8sum += w8;
          // contributes if the weight is positive
          double const w8pos = double(w8 > 0);
          add_product(cow, 3, xyz, w8pos);
          cow[3] += w8pos;
      } // iall
      if (cow[3] > 0) scale(cow, 3, 1./cow[3]); // normalize
      return w8sum; // returns the sum of weights
  } // center_of_weight

  // idea for a stable load balancer:
  //    for np processes, uses celing(log_2(nprocs)) iterations
  //    in iteration #0, place the center at 0,0,0 and the diagonal opposite corner finding its longest extent
  //                     divide the work by a plane sweep into chunks proportional to ceiling(np/2) and floor(np/2)
  //    pass the new processor numbers to the next iteration ...
  //    last iteration: number of processors is 1 --> done or 2 --> as before

  // find the largest extent:
  // 1) find the center of weight (not accounting the load but only if there is a non-zero load)
  // 2) find the maximum distance^2 and its position
  // 3) assume that the largest extent is along the line through the center of weight and that position
  // --> a cuboid will always use the space diagonal

  template <typename real_t, typename real_w_t=double>
  inline double plane_balancer(
        int const nprocs // number of MPI processes to distribute the work items evenly to
      , int const rank   // my MPI rank
      , size_t const nall // number of all work items
      , view2D<real_t> const & xyzw // xyzw[nall][4], positions [0/1/2] and weights [3] of the work items
      , real_w_t const w8s[]        // w8s[nall] weights separated
      , double const w8sum_all=1. // denominator of all weights
      , int const echo=0 // verbosity
      , double rank_center[4]=nullptr // export the rank center [0/1/2] and number of items [3]
      , int32_t *owner_rank=nullptr // export the rank of each task
  ) {
      // complexity is order(N^2) as each processes loops over all tasks in the first iteration

      auto constexpr epsilon = 1e-6; // accuracy for sanity check

      bool constexpr UNASSIGNED = 1, ASSIGNED = 0;
      std::vector<bool> state(nall, UNASSIGNED);

      std::vector<uint32_t> indirect(nall, 0);
      std::iota(indirect.begin(), indirect.end(), 0); // initialize with 0,1,....,nall-1

      size_t nuna{nall}; // number of unassigned blocks

      double load_now{w8sum_all};

      int np{nprocs}; // current number of processes among which we distribute
      int rank_offset{0}; // offset w.r.t. MPI ranks

      while (np > 1) {

          assert(rank_offset + np <= nprocs);
          // bisect the workload for np processors into (np + 1)/2 and np/2
          int const nhalf[] = {(np + 1) >> 1, np >> 1};
          // MPIrank ranges are {off ... off+nhalf[0]-1} and {off+nhalf[0] ...off+np-1}
          int const i01 = (rank >= rank_offset + nhalf[0]);
          // i01 == 0: this rank is part of the first (np + 1)/2 processes
          // i01 == 1: this rank is part of the second np/2 processes
          if (echo > 19) std::printf("# rank#%i divides %d into %d and %d\n", rank, np, nhalf[i01], nhalf[1 - i01]);
          assert(nhalf[0] + nhalf[1] == np);

          // determine the direction of the largest extent

          // determine the center of weight
          double cow[4];
          auto const w8sum = center_of_weight(cow, nuna, indirect.data(), xyzw);

          // determine the largest distance^2 from the center
          double maxdist2{-1}; int64_t imax{-1}; // block with the largest distance
          for (size_t iuna = 0; iuna < nuna; ++iuna) { // serial due to special reduction
              auto const iall = indirect[iuna];
              auto const *const xyz = xyzw[iall];
              auto const dist2 = pow2(xyz[X] - cow[X])
                               + pow2(xyz[Y] - cow[Y])
                               + pow2(xyz[Z] - cow[Z]);
              if (dist2 > maxdist2) { maxdist2 = dist2; imax = iall; }
          } // iuna

          if (imax < 0) {
              // this should only happens if the process is idle,
              //   i.e. there are no blocks assigned to this one
              load_now = 0;
              np = 1; // stop the while-loop
          } else {

              // determine a sorting direction
              add_product(cow, 3, xyzw[imax], -1.);
              auto const len2 = pow2(cow[X]) + pow2(cow[Y]) + pow2(cow[Z]);
              auto const norm = (len2 > 0) ? 1./std::sqrt(len2) : 0.0;
              double const vec[] = {cow[X]*norm, cow[Y]*norm, cow[Z]*norm};
              if (echo > 19) std::printf("# rank#%i sort along the [%g %g %g] direction\n", rank, vec[X], vec[Y], vec[Z]);

              using fui_t = std::pair<float,uint32_t>;
              std::vector<fui_t> v(nuna);

              for (size_t iuna = 0; iuna < nuna; ++iuna) { // parallel
                  auto const iall = indirect[iuna];
                  auto const *const xyz = xyzw[iall];
                  auto const f = xyz[X]*vec[X] + xyz[Y]*vec[Y] + xyz[Z]*vec[Z]; // inner product
                  v[iuna] = std::make_pair(float(f), uint32_t(iall));
              } // iall

              auto lambda = [](fui_t i1, fui_t i2) { return i1.first < i2.first; };
              std::stable_sort(v.begin(), v.end(), lambda);

              auto const by_np = 1./np, target_load0 = nhalf[0]*by_np*w8sum;
              {
                  auto const state0 = i01 ? ASSIGNED : UNASSIGNED,
                             state1 = i01 ? UNASSIGNED : ASSIGNED;
                  double load0{0}, load1{0};
                  size_t isrt;
                  for (isrt = 0; load0 < target_load0; ++isrt) { // serial
                      auto const iall = v[isrt].second;
                      load0 += w8s[iall];
                      state[iall] = state0;
                  }
                  for(; isrt < nuna; ++isrt) { // parallel, reduction(+:load1)
                      auto const iall = v[isrt].second;
                      load1 += w8s[iall];
                      state[iall] = state1;
                  } // i
                  assert(std::abs(load0 + load1 - w8sum) < epsilon*w8sum && "Maybe failed due to accuracy issues");
                  load_now = i01 ? load1 : load0;
              }

              if (echo > 19) std::printf("# rank#%i assign %g of %g (%.2f %%, target %.2f %%) to %d processes\n",
                                            rank, load_now, w8sum, load_now*100./w8sum, nhalf[i01]*by_np*100, nhalf[0]);

              // prepare for the next iteration
              rank_offset += i01*nhalf[0];
              np = nhalf[i01];
              // update the indirection list, determine which blocks are still unassigned
              size_t iunk{0}; // counter
              for (size_t iuna = 0; iuna < nuna; ++iuna) { // serial loop due to read-write access to array indirect
                  auto const iall = indirect[iuna]; // read from array indirect
                  if (UNASSIGNED == state[iall]) {
                      indirect[iunk] = iall; // write to array indirect at a lower address
                      ++iunk;
                  } // is_unassigned
              } // iall
              assert(iunk <= nuna);
              nuna = iunk; // set new number of unassigned blocks

          } // imax < 0

      } // while np > 1

      if (rank_center && load_now > 0) {
          // compute the center of weight again, for display and export
          double cow[4];
          auto const w8sum = center_of_weight(cow, nuna, indirect.data(), xyzw);
          if (echo > 13) std::printf("# rank#%i assign %.3f %% center %g %g %g, %g items\n",
                                        rank, w8sum*100/w8sum_all, cow[X], cow[Y], cow[Z], cow[W]);
          set(rank_center, 4, cow); // export
      } // rank_center

      if (echo > 9) std::printf("# rank#%i assign %.3f %%, target %.3f %%\n\n",
                                   rank, load_now*100/w8sum_all, 100./nprocs);

      if (owner_rank) { // export mask
          for (size_t iall = 0; iall < nall; ++iall) {
              if (UNASSIGNED == state[iall]) owner_rank[iall] = rank;
          } // iall
      } // owner_rank

      if (1) { // parallelized consistency check
          for (size_t iuna = 0; iuna < nuna; ++iuna) { // parallel
              auto const iall = indirect[iuna];
              auto const w8 = xyzw(iall,W);
              assert(w8 == w8s[iall] && "weights inconsistent");
          } // iuna
      } // consistency check

      return load_now;
  } // plane_balancer


  template <typename int_t>
  inline double get(
        uint32_t const comm_size
      , int32_t  const comm_rank
      , int_t const nb[3]
      , int const echo=0
      , double rank_center[4]=nullptr // export the rank center [0/1/2] and number of items [3]
      , int32_t *owner_rank=nullptr // export the rank of each task
  ) {
      // distribute a rectangular box of nb[X] x nb[Y] x nb[Z] with all weights 1

      assert(0 <= comm_rank && comm_rank < comm_size);
      auto const nall = size_t(nb[X])*size_t(nb[Y])*size_t(nb[Z]);
      assert(nall > 0);
      assert(nall <= (size_t(1) << 32) && "uint32_t is not long enough!");

      double w8sum_all{0};
      int constexpr W = 3;
      view2D<float> xyzw(nall, 4, 0.f);
      std::vector<double> w8s(nall, 0);

      for (int iz = 0; iz < nb[Z]; ++iz) {
      for (int iy = 0; iy < nb[Y]; ++iy) {
      for (int ix = 0; ix < nb[X]; ++ix) {
          auto const iall = size_t(iz*nb[Y] + iy)*nb[X] + ix;
//        assert(uint32_t(iall) == iall && "uint32_t is not long enough!");
          auto const w8 = 1.f; // weight(ix,iy,iz); // WEIGHTS CAN BE INSERTED HERE
          w8s[iall]    = w8;
          w8sum_all   += w8;
          xyzw(iall,W) = w8;
          xyzw(iall,X) = ix;
          xyzw(iall,Y) = iy;
          xyzw(iall,Z) = iz;
      }}} // ix iy iz

      return plane_balancer(comm_size, comm_rank, nall, xyzw, w8s.data(), w8sum_all, echo
                                          , rank_center, owner_rank);
  } // get


#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS
  status_t all_tests(int const echo=0); // declaration only
#endif // NO_UNIT_TESTS

} // namespace load_balancer
