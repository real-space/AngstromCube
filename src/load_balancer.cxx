// This file is part of AngstromCube under MIT License
/*
 * load balancing for A43
 *    main purpose is the distribution of right-hand-side blocks
 *    to MPI processes and GPUs
 *
 */


#include <cstdio> // std::printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <cstdint> // size_t, uint32_t
#include <algorithm> // std::max, ::min, ::swap
#include <utility> // std::pair
#include <cmath> // std::ceil, ::pow, ::cbrt
#include <numeric> // std::iota
#ifndef   NO_UNIT_TESTS
  #include <cstdlib> // rand, RAND_MAX
  #include "progress_report.hxx" // ProgressReport
#endif // NO_UNIT_TESTS

#include "load_balancer.hxx"

#include "inline_math.hxx" // set, pow2
#include "status.hxx" // status_t
#include "simple_stats.hxx" // ::Stats<>
#include "recorded_warnings.hxx" // warn
#include "constants.hxx" // ::pi
#include "print_tools.hxx" // printf_vector

#ifndef   NO_UNIT_TESTS
  #include "control.hxx" // ::get
  #ifdef    HAS_BITMAP_EXPORT
    #include "bitmap.hxx" // ::write_bmp_file
  #endif // HAS_BITMAP_EXPORT
#endif // NO_UNIT_TESTS

namespace load_balancer {

// #define   LOAD_BALANCER_DRAW_SVG
#ifdef    LOAD_BALANCER_DRAW_SVG
  static std::vector<double> draw2D; // global field
#endif // LOAD_BALANCER_DRAW_SVG

  int constexpr X=0, Y=1, Z=2, W=3;

  template <typename real_t>
  double center_of_weight(
        double cow[4] // result: center of weight [0/1/2] and number of contributors [3]
      , size_t const nuna // number of unassigned work items
      , uint32_t const indirect[] // list of unassigned work items
      , real_t const (*const xyzw)[4] // positions of work items in space [0/1/2] and their weight [3]
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
          cow[W] += w8pos;
      } // iall
      if (cow[W] > 0) scale(cow, 3, 1./cow[3]); // normalize
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
  double plane_balancer(
        int const nprocs // number of MPI processes to distribute the work items evenly to
      , int const rank   // my MPI rank
      , size_t const nall // number of all work items
      , real_t const (*const xyzw)[4] // xyzw[nall][4], positions [0/1/2] and weights [3] of the work items
      , real_w_t const w8s[]        // w8s[nall] weights separated
      , double const w8sum_all=1. // denominator of all weights
      , int const echo=0 // verbosity
      , double rank_center[4]=nullptr // export the rank center [0/1/2] and number of items [3]
      , uint16_t *const owner_rank=nullptr // export the rank of each task, [nall]
  ) {
      // complexity is order(N^2) as each processes loops over all tasks in the first iteration

      auto constexpr epsilon = 1e-6; // accuracy for sanity check

      bool constexpr UNASSIGNED = 1, ASSIGNED = 0;
      std::vector<bool> state(nall, UNASSIGNED);

      assert(nall <= (1ull << 32) && "using uint32_t indirection lists limits the total number to 2^32");
      std::vector<uint32_t> indirect(nall);
      std::iota(indirect.begin(), indirect.end(), uint32_t(0)); // initialize with {0, 1, ..., nall-1}

      size_t nuna{nall}; // number of unassigned blocks

      double load_now{w8sum_all};

      int np{nprocs}; // current number of processes among which we distribute
      int rank_offset{0}; // offset w.r.t. MPI ranks

      int tree_level{0};
      while (np > 1) {

          assert(rank_offset + np <= nprocs);
          // bisect the workload for np processors into (np + 1)/2 and np/2
          int const nhalf[] = {(np + 1) >> 1, np >> 1};
          // MPIrank ranges are 0:{off, ..., off+nhalf[0]-1} and 1:{off+nhalf[0], ..., off+np-1}
          int const i01 = (rank >= rank_offset + nhalf[0]);
          // i01 == 0: this rank is part of the first (np + 1)/2 processes
          // i01 == 1: this rank is part of the second np/2 processes
          if (echo > 19) std::printf("# rank#%i divides %d into %d and %d\n", rank, np, nhalf[i01], nhalf[1 - i01]);
          assert(nhalf[0] + nhalf[1] == np);

          // determine the direction of the largest extent

          // determine the center of weight
          double cow[4];
          auto const w8sum = center_of_weight(cow, nuna, indirect.data(), xyzw);

          // determine the largest distance^2 from the center of weight
          double maxdist2{-1}; int64_t imax{-1}; // block with the largest distance
          for (size_t iuna = 0; iuna < nuna; ++iuna) { // serial due to special reduction
              auto const iall = indirect[iuna];
              auto const *const xyz = xyzw[iall];
              auto const dist2 = pow2(xyz[X] - cow[X]) + pow2(xyz[Y] - cow[Y]) + pow2(xyz[Z] - cow[Z]);
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
                  auto const f = xyz[X]*vec[X] // inner product
                               + xyz[Y]*vec[Y]
                               + xyz[Z]*vec[Z];
                  v[iuna].first  = f;
                  v[iuna].second = iall;
              } // iall

              auto lambda = [](fui_t i1, fui_t i2) { return i1.first < i2.first; };
              std::stable_sort(v.begin(), v.end(), lambda);

              auto const target_load0 = nhalf[0]*w8sum; // relative target for load0 multiplied with np
              { // scope: distribute according to target loads
                  auto const state0 = i01 ? ASSIGNED : UNASSIGNED,
                             state1 = i01 ? UNASSIGNED : ASSIGNED;
                  double load0{0}, load1{0}; // relative load multiplied with np to avoid floating point errors
                  size_t isrt{0};
                  for (isrt = 0; load0 < target_load0; ++isrt) { // serial
                      auto const iall = v[isrt].second;
                      load0 += w8s[iall]*np;
                      state[iall] = state0;
                  } // while load1 < target_load1
                  auto const isrt_middle = isrt;
                  for(; isrt < nuna; ++isrt) { // parallel reduction(+:load1)
                      auto const iall = v[isrt].second;
                      load1 += w8s[iall]*np;
                      state[iall] = state1;
                  } // isrt
                  assert(std::abs((load0 + load1) - (w8sum*np)) < epsilon*(w8sum*np) && "Maybe failed due to accuracy issues");
                  load_now = (i01 ? load1 : load0)/np;

                  if (echo > 29) std::printf("# plane level=%d %g %g %g isrt=%lu %d|%d\n", tree_level, vec[X], vec[Y], vec[Z], isrt_middle, nhalf[0],nhalf[1]);
#ifdef    LOAD_BALANCER_DRAW_SVG
                  if (echo > 29) { // show bisecting plane
                      // bisecting plane normal is the sorting vector vec, plane distance from the origin is ?
                      double pd{0}; int den{0};
                      if (isrt_middle < nuna) { pd += v[isrt_middle].first; ++den; }; // distance of the point that is closest to the plane and belongs to load1
                      if (isrt_middle > 0)    { pd += v[isrt_middle - 1].first; ++den; } // distance of the ... belongs to load0
                      if (rank == rank_offset) {
                          std::printf("plane level=%d %g %g %g  dist= %g  isrt=%lu %d|%d\n", tree_level, vec[X], vec[Y], vec[Z], pd/den, isrt_middle, nhalf[0],nhalf[1]);
                          // store the 2D plane in a global variable to be drawn into an SVG later
                          auto const s = draw2D.size();
                          if (s > 0) {
                              draw2D.resize(s + 4);
                              draw2D[s + 0] = vec[X];
                              draw2D[s + 1] = vec[Y];
                              draw2D[s + 2] = pd/den;
                              draw2D[s + 3] = tree_level;
                          } // s > 0
                      } // rank == rank_offset
                  } // echo
#endif // LOAD_BALANCER_DRAW_SVG
              } // scope

              if (echo > 19) std::printf("# rank#%i assign %g of %g (%.2f %%, target %.2f %%) to %d processes\n",
                                            rank, load_now, w8sum, load_now*100./w8sum, nhalf[i01]*100./np, nhalf[0]);

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

          ++tree_level;
      } // while np > 1

      if (nullptr != rank_center) {
          set(rank_center, 4, 0.0);
          if (load_now > 0) {
              // compute the center of weight again, for display and export
              double cow[4];
              auto const w8sum = center_of_weight(cow, nuna, indirect.data(), xyzw);
              if (echo > 13) std::printf("# rank#%i assign %.3f %% center %g %g %g, %g items\n",
                                            rank, w8sum*100/w8sum_all, cow[X], cow[Y], cow[Z], cow[W]);
              set(rank_center, 4, cow); // export
          } // load_now > 0
      } // rank_center

      if (echo > 9) std::printf("# rank#%i load target %.3f %%, assign %.3f %%\n",
                                   rank, 100./nprocs, load_now*100/w8sum_all);

      if (nullptr != owner_rank) {
          for (size_t iall = 0; iall < nall; ++iall) {
              if (UNASSIGNED == state[iall]) {
                  assert(no_owner == owner_rank[iall]);
                  owner_rank[iall] = rank;
                  assert(owner_rank[iall] == rank && "uint16_t too short for owner_ranks");
              } // unassigned
          } // iall
          // Beware: only the owned entries of owner_rank have been modified, so
          //         an MPI_MIN-Allreduce needs to be performed later. This is skipped here
          //         to maintain serial executability of this routine.
          if (echo > 99) {
              std::printf("# rank#%i owner_rank before MPI_MIN ", rank);
              printf_vector(" %i", owner_rank, nall);
          } // echo
      } // owner_rank

      if (1) { // parallelized consistency check
          for (size_t iuna = 0; iuna < nuna; ++iuna) { // parallel
              auto const iall = indirect[iuna];
              auto const w8 = xyzw[iall][W];
              assert(w8 == w8s[iall] && "weights inconsistent");
          } // iuna
      } // consistency check

      return load_now;
  } // plane_balancer


  double get(
        uint32_t const comm_size // number of MPI processes in this communicator
      , int32_t  const comm_rank // rank of this MPI process
      , uint32_t const nb[3] // number of blocks in X/Y/Z direction
      , int const echo // =0, log level
      , double rank_center[4] // =nullptr, export the rank center [0/1/2] and number of items [3]
      , uint16_t *const owner_rank // =nullptr, export the owner rank of each task, [nb[Z]*nb[Y]*nb[X]]
  ) {
      // distribute a rectangular box of nb[X] x nb[Y] x nb[Z] with all weights 1

      assert(comm_rank >= 0);
      assert(comm_rank < comm_size);

      auto const nall = nb[X]*size_t(nb[Y])*size_t(nb[Z]);
      assert(nall > 0);
      assert(nall <= (size_t(1) << 32) && "uint32_t is not long enough!");

      double w8sum_all{0};
      auto const xyzw = new float[nall][4];
      std::vector<double> w8s(nall, 0);

      for (uint32_t iz = 0; iz < nb[Z]; ++iz) {
      for (uint32_t iy = 0; iy < nb[Y]; ++iy) {
      for (uint32_t ix = 0; ix < nb[X]; ++ix) {
          auto const iall = (size_t(iz)*nb[Y] + iy)*nb[X] + ix;
//        assert(uint32_t(iall) == iall && "uint32_t is not long enough!");
          float const w8 = 1.f; // weight(ix,iy,iz); // WEIGHTS CAN BE INSERTED HERE
          w8s[iall]     = w8;
          w8sum_all    += w8;
          xyzw[iall][W] = w8;
          xyzw[iall][X] = ix;
          xyzw[iall][Y] = iy;
          xyzw[iall][Z] = iz;
      }}} // ix iy iz

      auto const load_now = plane_balancer(comm_size, comm_rank, nall, xyzw, w8s.data(), w8sum_all, echo
                                          , rank_center, owner_rank);
      delete[] xyzw;
      return load_now;
  } // get

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <typename real_t>
  real_t analyze_load_imbalance(real_t const load[], int const nprocs, int const echo=1) {
      simple_stats::Stats<> st(0);
      for (int rank = 0; rank < nprocs; ++rank) {
          st.add(load[rank]);
          if (echo > 19) std::printf("# myrank=%i owns %g\n", rank, load[rank]);
      } // rank
      auto const mean = st.mean(), max = st.max();
      if (echo > 0) std::printf("# %ld processes own %g blocks, per process [%g, %.2f +/- %.2f, %g]\n",
                                    st.tim(), st.sum(), st.min(), mean, st.dev(), max);
      return (mean > 0) ? max/mean : -1;
  } // analyze_load_imbalance


  template <typename real_t>
  real_t distance_squared(real_t const a[], real_t const b[]) {
      return pow2(a[X] - b[X]) + pow2(a[Y] - b[Y]) + pow2(a[Z] - b[Z]);
  } // distance_squared

  double intersect(double xy[2], double const v1[3], double const v2[3], float const threshold=1e-7) {
      auto const d1 = v1[2], d2 = v2[2];
      // Given two lines with normal vectors v1 and v2 and distances to the origin d1 and d2, respectively.
      // Compute their intersection (if any)
      double const det2x2 = v1[0]*v2[1] - v1[1]*v2[0]; // 2d determinant
      xy[0] = 0; xy[1] = 0; // init result (should not be used if false is returned)
      // we look for a point (x,y) on both lines, i.e.
      //         v1 dot xy == d1
      // and simultaneousy
      //         v2 dot xy == d2
      // set result
      double const denom = 1./det2x2;
      if (denom == denom) {
          xy[0] = (v2[1]*d1 - v1[1]*d2)*denom;
          xy[1] = (v1[0]*d2 - v2[0]*d1)*denom;
//        std::printf("# line1 (%g,%g)*(x,y) == %g intersects with line2 (%g,%g)*(x,y) == %g at (%g,%g)\n", v1[0],v1[1], d1, v2[0],v2[1], d2, xy[0],xy[1]);
      } // denom is not NaN
      return std::abs(det2x2);
  } // intersect

  status_t test_plane_balancer(int const nprocs, int const n[3], int const echo=0) {
      status_t stat(0);

      if (nprocs < 1) return stat;

      auto const nall = (n[X])*size_t(n[Y])*size_t(n[Z]);
      assert(nall <= (size_t(1) << 32) && "uint32_t is not long enough!");

      double w8sum_all{0};
      int constexpr W = 3;
      auto const xyzw = new float[nall][4];
      std::vector<double> w8s(nall, 0);

      int const holes = control::get("load_balancer.test.holes", 0.);
      auto const hole_radius_squared = pow2(control::get("load_balancer.test.holes.radius", 8.));

      for (int iz = 0; iz < n[Z]; ++iz) {
      for (int iy = 0; iy < n[Y]; ++iy) {
      for (int ix = 0; ix < n[X]; ++ix) {
          auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
//        assert(uint32_t(iall) == iall && "uint32_t is not long enough!");
          float h{1};
          for (int ih = 1-holes; ih < holes; ih += 2) {
              auto const x_hole = n[X]*ih/(2.*holes);
              auto const r2 = pow2(ix - .5*n[X] - x_hole) + pow2(iy - .5*n[Y]) + pow2(iz - .5*n[Z]);
              h *= (r2 > hole_radius_squared); // radius_squared
          } //
          float const w8 = 1.f*h; // weight(ix,iy,iz); // WEIGHTS CAN BE INSERTED HERE
          w8s[iall]     = w8;
          w8sum_all    += w8;
          xyzw[iall][W] = w8;
          xyzw[iall][X] = ix;
          xyzw[iall][Y] = iy;
          xyzw[iall][Z] = iz;
      }}} // ix iy iz
      double const longest_possible_distance = std::sqrt(pow2(n[X]) + pow2(n[Y]) + pow2(n[Z]));

      std::vector<double> load(nprocs, 0.0);
      auto const rank_center = new double[nprocs][4];
      bool constexpr compute_rank_centers = true;
      std::vector<uint16_t> owner_rank(nall, no_owner);

      if (echo > 0) std::printf("# %s: distribute %g blocks to %d processes\n\n", __func__, w8sum_all, nprocs);

#ifdef    LOAD_BALANCER_DRAW_SVG
      draw2D.resize(2); draw2D[0] = n[X]; draw2D[1] = n[Y]; // init
#endif // LOAD_BALANCER_DRAW_SVG

      int const echo_rank0 = control::get("load_balancer.test.echo.rank0", 0.); // increase the verbosity for rank0

      ProgressReport timer(__FILE__, __LINE__, 2.5, echo); // update the line every 2.5 seconds
      for (int rank = 0; rank < nprocs; ++rank) {
          load[rank] = plane_balancer(nprocs, rank, nall, xyzw, w8s.data(), w8sum_all, echo + (0 == rank)*echo_rank0
                                            , rank_center[rank], owner_rank.data());
          timer.report(rank, nprocs);
      } // rank

#ifdef    LOAD_BALANCER_DRAW_SVG
      { // scope: SVG
          // assume that draw2D is an array of sets of 4 doubles which results from a depth-first traversal of the bisection tree
          int const nplanes = draw2D.size()/4;
          assert(4*nplanes + 2 == draw2D.size());
          auto const nx = int(draw2D[0]), ny = int(draw2D[1]);
          if (echo > 2) std::printf("# found %ld planes for https://editsvgcode.com/\n\n", nplanes);
          if (echo > 2) std::printf("<!-- SVG code generated by %s -->\n", __FILE__);
          if (echo > 2) std::printf("<svg viewBox=\"%d %d %d %d\" xmlns=\"http://www.w3.org/2000/svg\">\n", -10, -10, nx + 20, ny + 20);
          double const frame[4][4] = {{1,0,0,0}, {0,1,0,0}, {1,0,1.*nx,0}, {0,1,1.*ny,0}};
          assert(0 == draw2D[5] && "the 1st plane must be the origin");
          std::vector<int> ancestor(32, -1);
          std::vector<int8_t>  side(32, -1);
          for (int ip = 0; ip < nplanes; ++ip) {
              double const *const v1 = &draw2D[ip*4 + 2];

              int const tree_level = v1[3];
              assert(tree_level >= 0 && tree_level < 32);
              side[tree_level] = (-1 == ancestor[tree_level]) ? 0 : 1;
              ancestor[tree_level] = ip;
              if (0) {
                std::printf("  <!-- I am plane #%i, level=%d, my ancestors are", ip, tree_level);
                for (int jp = 0; jp < tree_level; ++jp) {
                    std::printf(" %i", ancestor[jp]);
                } // jp 
                std::printf(" -->\n");
              } // 0

              double points[99][2];
              int npoints{0}; int ipoint[99];
              // determine who are my ancestors i.e. which lines are my parents and grandparents and so on.
              // The tree has been traversed depth-first
              for (int jp = tree_level - 1; jp >= -4; --jp) { // loops over ancestor lines and frame 
         //   for (int jp = -4; jp < tree_level; ++jp) { // loops over frame and ancestor lines
                  double const *const v2 = (jp < 0) ? frame[jp + 4] : &draw2D[ancestor[jp]*4 + 2];
                  if (true) {
                      // compute intersection of the lines
                      auto const intersects = intersect(points[npoints], v1, v2);
                      if (intersects > 1e-12) {
                          ipoint[npoints] = (jp < 0) ? jp : ancestor[jp];
                          auto const x = points[npoints][0], y = points[npoints][1];
                          double constexpr eps = 1e-9;
                          // check if they are within the border rect [0...nx, 0...ny]
                          bool const right_side = (jp < 4) || (int(x*v2[0] + y*v2[1] <= v2[2]) == side[jp]);
                          if (right_side && (x > -eps) && (x < nx + eps) && (y > -eps) && (y < ny + eps)) ++npoints; // accept the point
                      }
                  } // level index is higher
              } // jp
              if (npoints > 1) {
                  if (echo > 2) std::printf("  <!-- line #%i has %d points, take #%i and #%i -->\n", ip, npoints, ipoint[0], ipoint[1]);
                  if (echo > 2) std::printf("  <line x1=\"%g\" y1=\"%g\" x2=\"%g\" y2=\"%g\" stroke=\"black\" />\n",
                                                points[0][0], points[0][1], points[1][0], points[1][1]);
              } // npoints > 1
          } // ip
          if (echo > 2) std::printf("  <rect width=\"%d\" height=\"%d\" x=\"%d\" y=\"%d\" fill=\"none\" stroke=\"grey\" />\n", nx, ny, 0, 0);
          if (echo > 2) std::printf("</svg>\n\n");
      } // scope: SVG
#endif // LOAD_BALANCER_DRAW_SVG



      analyze_load_imbalance(load.data(), nprocs, echo);

      if (compute_rank_centers && nprocs > 1) {
          // analyze positions of rank centers
          double mindist2{9e300}; int ijmin[] = {-1, -1};
          double maxdist2{-1.0};  int ijmax[] = {-1, -1};
          simple_stats::Stats<> st2, st1;
          double const wbin = control::get("load_balancer.bin.width", 0.25), invbin = 1./wbin;
          int const nbin = 1 + int(longest_possible_distance/wbin);
          std::vector<uint32_t> hist(nbin, 0);
          int np{0}; // counter for the number of processes with a non-zero load
          for (int irank = 0; irank < nprocs; ++irank) {
              if (load[irank] > 0) {
                  ++np;
                  for (int jrank = 0; jrank < nprocs; ++jrank) { // self-avoiding triangular loop
                      if (load[jrank] > 0) {
                          auto const dist2 = distance_squared(rank_center[irank], rank_center[jrank]);
                          if (dist2 > 0 && dist2 < mindist2) { mindist2 = dist2; ijmin[0] = irank; ijmin[1] = jrank; }
                          if (dist2 > maxdist2) { maxdist2 = dist2; ijmax[0] = irank; ijmax[1] = jrank; }
                          auto const dist = std::sqrt(dist2);
                          if (echo > 15) std::printf("# distance-ij is %g\n", dist);
                          int const ibin = dist*invbin; // floor
                          ++hist[std::min(ibin, nbin - 1)];
                          st2.add(dist2);
                          st1.add(dist);
                      } // load
                  } // jrank
              } // load
          } // irank
          auto const maxdist = std::sqrt(maxdist2), mindist = std::sqrt(mindist2);
          if (echo > 1) std::printf("# shortest distance between centers is %g between rank#%i and #%i, longest is %g\n",
                                        mindist, ijmin[0], ijmin[1], maxdist);
          if (echo > 9) std::printf("# longest distance between centers is %g between rank#%i and #%i, shortest is %g\n",
                                        maxdist, ijmax[0], ijmax[1], mindist);
          if (echo > 7) {
              double const denom = 1./pow2(std::max(1, np));
              std::printf("## center-distance histogram, bin width %g\n", wbin);
              for (int ibin = 0; ibin < nbin; ++ibin) {
                  std::printf("%g %g\n", ibin*wbin, hist[ibin]*denom);
              } // ibin
              std::printf("\n\n");
          } // echo
          if (echo > 2) std::printf("# stats: distance [%g, %g +/- %g, %g]\n"
                                    "#        distance^2 [%g, %g +/- %g, %g]\n",
                                    st1.min(), st1.mean(), st1.dev(), st1.max(),
                                    st2.min(), st2.mean(), st2.dev(), st2.max());
      } // compute_rank_centers
      delete[] rank_center;

      // check masks
      if (1) {
          int strange{0};
          for (int iz = 0; iz < n[Z]; ++iz) {
            for (int iy = 0; iy < n[Y]; ++iy) {
              for (int ix = 0; ix < n[X]; ++ix) {
                  auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
                  auto const owner = owner_rank[iall];
                  strange += (no_owner == owner); // under-assignement
                  if (no_owner == owner) warn("work item %d %d %d has not been assigned to any rank", ix,iy,iz);
              } // ix
            } // iy
          } // iz
          if (strange) warn("strange: %d under-assignments", strange);
      }

      if (1 == n[Z] && echo > 5 && n[X] <= 300 && n[Y] <= 300) {
          std::printf("\n# visualize plane balancer %d x %d on %d processes:%s", n[Y], n[X], nprocs,
                              nprocs > 64 ? " (symbols are not unique!)" : "");
          int constexpr iz = 0;
          char const chars[65] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ<>";
          char constexpr mask_char = ' ';
          for (int iy = 0; iy < n[Y]; ++iy) {
              std::printf("\n# ");
              for(int ix = 0; ix < n[X]; ++ix) {
                  auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
                  auto const owner = owner_rank[iall];
                  auto const c = (xyzw[iall][W] < 1) ? mask_char : chars[owner & 0x3f];
                  std::printf("%c", c);
              } // ix
          } // iy
          std::printf("\n#\n\n");
          //
          // ToDo: to export the Voronoi diagrams, we need to access the plane normals
          //       and plane parameters at which the plane separates the two processes
          //       and the tree structure in order to know in which half space the plane is valid
          //
#ifdef    HAS_BITMAP_EXPORT
          // create a color map
          int const n3 = std::ceil(std::cbrt(nprocs));      assert(n3*n3*n3 >= nprocs);
          std::vector<uint8_t> bmp_color(n3*n3*n3*4, 0);
          { // scope: generate >=nprocs distinct colors
              int iproc{0};
              for (int cx = 1; cx <= n3; ++cx) {
                  for (int cy = 1; cy <= n3; ++cy) {            // start from 1 to avoid black(0x000000)
                      for (int cz = 1; cz <= n3; ++cz) {
                          bmp_color[iproc*4    ] = (250*cz)/n3;
                          bmp_color[iproc*4 + 1] = (250*cy)/n3; // use 250 instead of 255 to avoid white(0xffffff)
                          bmp_color[iproc*4 + 2] = (250*cx)/n3;
                          ++iproc;
              }}} // cx cy cz
              assert(iproc >= nprocs && "We need enough distinct colors");
          } // scope

          // color the owned region
          std::vector<uint8_t> bmp_data(n[Y]*n[X]*4, 0);
          for (int iy = 0; iy < n[Y]; ++iy) {
              for(int ix = 0; ix < n[X]; ++ix) {
                  auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
                  auto const owner = owner_rank[iall];
                  for (int rgb = 0; rgb < 3; ++rgb) {
                      bmp_data[iall*4 + rgb] = bmp_color[owner*4 + rgb];
                      if (xyzw[iall][W] < 1) bmp_data[iall*4 + rgb] = 255; // white(0xffffff)
                  } // rgb
              } // ix
          } // iy

          // store as image
          char filename[64]; std::snprintf(filename, 64, "load_balancer-n%d-%dx%d", nprocs, n[X], n[Y]);
          if (echo > 5) std::printf("# try to create bitmap file %s with %d colors for %d domains\n", filename, n3*n3*n3, nprocs);
          stat += bitmap::write_bmp_file(filename, bmp_data.data(), n[Y], n[X], -1, 1.f);
#endif // HAS_BITMAP_EXPORT
      } // visualize
      delete[] xyzw;
      return stat;
  } // test_plane_balancer

  // example 5 processes -->
  //        rank#0      rank#1      rank#2      rank#3      rank#4
  //  take    3/5         3/5         3/5         2/5         2/5
  //  take    2/3         2/3         1/3         1/2         1/2
  //  take    1/2         1/2

  inline double random_between_0_and_1() {
      double constexpr rand_denom = 1./RAND_MAX;
      return rand()*rand_denom;
  } // random_between_0_and_1


  status_t test_reference_point_cloud(int const nxyz[3], int const echo=9) {
      auto const maxdist_diagonal = std::sqrt(pow2(nxyz[X]) + pow2(nxyz[Y]) + pow2(nxyz[Z]));
      double const wbin = control::get("load_balancer.bin.width", 0.25), invbin = 1./wbin;
      int const nbin = int(maxdist_diagonal/wbin) + 1;
      if (echo > 0) std::printf("# reference point-cloud histogram for %d x %d x %d points, bin width %g, %d bins\n",
                                                                  nxyz[X], nxyz[Y], nxyz[Z],        wbin, nbin);
      std::vector<uint32_t> hist(nbin, 0);
      simple_stats::Stats<> st2, st1;
      for (int iz = 0; iz < nxyz[Z]; ++iz) {
      for (int iy = 0; iy < nxyz[Y]; ++iy) {
      for (int ix = 0; ix < nxyz[X]; ++ix) {
              for (int jz = 0; jz < nxyz[Z]; ++jz) {
              for (int jy = 0; jy < nxyz[Y]; ++jy) {
              for (int jx = 0; jx < nxyz[X]; ++jx) {
//                    double const rnd[] = {0.5*random_between_0_and_1() - 0.25,
//                                          0.5*random_between_0_and_1() - 0.25,
//                                          0.5*random_between_0_and_1() - 0.25};
                      int constexpr rnd[] = {0, 0, 0};
                      double const dist2 = pow2(ix - jx + rnd[X])
                                         + pow2(iy - jy + rnd[Y])
                                         + pow2(iz - jz + rnd[Z]);
                      auto const dist = std::sqrt(dist2);
                      st2.add(dist2);
                      st1.add(dist);
                      if (echo > 15) std::printf("# distance-ij is %g\n", dist);
                      int const ibin = dist*invbin; // floor
                      ++hist[ibin];
              }}} // jx jy jz
      }}} // ix iy iz
      if (echo > 5) {
          double const by_n = 1./(1.*nxyz[X]*nxyz[Y]*nxyz[Z]);
          double const denom = pow2(by_n);
          std::printf("## point-distance histogram, bin width %g\n", wbin);
          for (int ibin = 0; ibin < nbin; ++ibin) {
              auto const radius = ibin*wbin;
              // number of points inside a sphere shell of from radius r - wbin to r is 4/3*pi*(r^3 - (r - wbin)^3)*rho
              // so in the limit of small wbin, this becomes 4*pi*r^2*wbin*rho
              auto const analytical = 4*constants::pi*pow2(radius)*wbin*by_n;
              std::printf("%g %g %g\n", radius, hist[ibin]*denom, analytical);
          } // ibin
          std::printf("\n\n");
      } // echo
      if (echo > 2) std::printf("# stats: distance [%g, %g +/- %g, %g]\n"
                                "#        distance^2 [%g, %g +/- %g, %g]\n",
                                st1.min(), st1.mean(), st1.dev(), st1.max(),
                                st2.min(), st2.mean(), st2.dev(), st2.max());
      return 0;
  } // test_reference_point_cloud


  status_t all_tests(int const echo) {
      status_t stat(0);

      auto const niso = control::get("load_balancer.test.n", 19.);
      int const nxyz[] = {int(control::get("load_balancer.test.nx", niso)), // number of blocks
                          int(control::get("load_balancer.test.ny", niso)),
                          int(control::get("load_balancer.test.nz", 1.))};
      auto const nprocs = int(control::get("load_balancer.test.nprocs", 53.));
      if (echo > 0) std::printf("\n\n# %s start %d x %d x %d = %d with %d MPI processes\n",
                      __func__, nxyz[X], nxyz[Y], nxyz[Z], nxyz[X]*nxyz[Y]*nxyz[Z], nprocs);

      stat += test_plane_balancer(nprocs, nxyz, echo);
//    stat += test_reference_point_cloud(nxyz, echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace load_balancer
