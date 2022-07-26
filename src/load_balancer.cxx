
/*
 * Where does A43 need load balancing?
 *
 *
 */

#include <cassert> // assert
#include <cstdint> // int64_t, size_t, int32_t, uint8_t
#include <cstdio> // std::printf
#include <algorithm> // std::max, ::min, ::swap
#include <numeric> // std::iota
#include <vector> // std::vector<T>
#include <cmath> // std::ceil, ::pow, ::cbrt
#include <cstdlib> // rand

#include "load_balancer.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
// #include "mpi_parallel.hxx" // ...
#include "simple_stats.hxx" // ::Stats<>
#include "constants.hxx" // ::pi
#include "data_view.hxx" // view4D<T>
#include "inline_math.hxx" // pow2, scale, align<nBits>
#include "recorded_warnings.hxx" // error

namespace load_balancer {

  template <typename number_t, int sgn=1> // use sgn=-1 for smallest_of_3
  int largest_of_3(number_t const n[3]) {
      if (n[0]*sgn >= n[1]*sgn) {
          // answer is 0 or 2
          return (n[0]*sgn >= n[2]*sgn) ? 0 : 2;
      } else {
          // answer is 1 or 2
          return (n[1]*sgn >= n[2]*sgn) ? 1 : 2;
      }
  } // largest_of_3

} // namespace load_balancer


#ifndef  NO_UNIT_TESTS

  #include "control.hxx" // ::get

#define weight(x,y,z) 1.f

namespace load_balancer {

  template <typename real_t>
  inline real_t analyze_load_imbalance(real_t const load[], int const nprocs, int const echo=1) {
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

  int constexpr X=0, Y=1, Z=2;

  template <typename real_t>
  inline real_t distance_squared(real_t const a[], real_t const b[]) {
      return pow2(a[X] - b[X]) + pow2(a[Y] - b[Y]) + pow2(a[Z] - b[Z]);
  } // distance_squared


  // idea for a stable load balancer or at least a stable initial guess:
  //    for np processes, uses celing(log_2(np)) iterations
  //    in iteration #0, place the center at 0,0,0 and the diagonal opposite corner finding its longest extent
  //                     divide the work by a plane sweep into chunks proportional to ceiling(np/2) and floor(np/2)
  //    pass the new processor numbers to the next iteration ...
  //    last iteration: number of processors is 1 --> done or 2 --> as before

  // find the largest extent:
  // 1) find the center of weight (not accounting the load but only if there is a non-zero load)
  // 2) find the maximum distance^2 and its position
  // 3) assume that the largest extent is along the line through the center of weight and that position
  // --> a cuboid will always use the space diagonal

  inline status_t test_plane_balancer(int const nprocs, int const n[3], int const echo=0) {
      status_t stat(0);

      if (nprocs < 1) return stat;

      auto const nall = (n[X])*size_t(n[Y])*size_t(n[Z]);
      assert(nall <= (size_t(1) << 32) && "uint32_t is not long enough!");

      double w8sum_all{0};
      int constexpr W = 3;
      view2D<float> xyzw(nall, 4, 0.f);
      std::vector<float> w8s(nall, 0.f);

      for (int iz = 0; iz < n[Z]; ++iz) {
      for (int iy = 0; iy < n[Y]; ++iy) {
      for (int ix = 0; ix < n[X]; ++ix) {
          auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
//        assert(uint32_t(iall) == iall && "uint32_t is not long enough!");
          xyzw(iall,X) = ix;
          xyzw(iall,Y) = iy;
          xyzw(iall,Z) = iz;
          auto const w8 = weight(ix,iy,iz);
          w8s[iall]  = w8;
          w8sum_all += w8;
          xyzw(iall,W) = w8;
      }}} // ix iy iz

      std::vector<double> load(nprocs, 0.0);
      view2D<double> rank_center(nprocs, 4, 0.0);
      bool constexpr compute_rank_centers = true;

      if (echo > 0) std::printf("# %s: distribute %g blocks to %d processes\n\n", __func__, w8sum_all, nprocs);

for (int rank = 0; rank < nprocs; ++rank) {


      bool constexpr UNASSIGNED = 0, ASSIGNED = 1;
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
          double cow[] = {0, 0, 0, 0};
          double w8sum{0};
          for (size_t iuna = 0; iuna < nuna; ++iuna) {
              auto const iall = indirect[iuna];
              auto const *const xyz = xyzw[iall];
              double const w8 = xyz[W];
//            assert(w8 == w8s[iall]);
              w8sum += w8;
              // contributes if the weight is positive
              double const w8pos = double(w8 > 0);
              add_product(cow, 3, xyz, w8pos);
              cow[3] += w8pos;
          } // iall
          if (cow[3] > 0) scale(cow, 3, 1./cow[3]);

          // determine the largest distance^2 from the center
          double maxdist2{-1};
          int64_t imax{-1}; // index of the block with the largest distance
          for (size_t iuna = 0; iuna < nuna; ++iuna) {
              auto const iall = indirect[iuna];
              auto const *const xyz = xyzw[iall];
              auto const dist2 = pow2(xyz[X] - cow[X])
                               + pow2(xyz[Y] - cow[Y])
                               + pow2(xyz[Z] - cow[Z]);
              if (dist2 > maxdist2) {
                  maxdist2 = dist2;
                  imax = iall;
              } // dist2 > maxdist2
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

              for (size_t iuna = 0; iuna < nuna; ++iuna) {
                  auto const iall = indirect[iuna];
                  auto const *const xyz = xyzw[iall];
                  auto const f = xyz[X]*vec[X] + xyz[Y]*vec[Y] + xyz[Z]*vec[Z]; // inner product
                  v[iuna] = std::make_pair(float(f), uint32_t(iall));
              } // iall

              auto lambda = [](fui_t i1, fui_t i2) { return i1.first < i2.first; };
              std::stable_sort(v.begin(), v.end(), lambda);

              auto const target_frac = nhalf[0]/double(np), // the "larger" half
                         target_load = target_frac*w8sum;
              {
                  auto const state0 = i01 ? ASSIGNED : UNASSIGNED,
                             state1 = i01 ? UNASSIGNED : ASSIGNED;
                  double load0{0}, load1{0};
                  size_t isrt;
                  for (isrt = 0; load0 < target_load; ++isrt) {
                      auto const iall = v[isrt].second;
                      load0 += w8s[iall];
                      state[iall] = state0;
                  }
                  for(; isrt < nuna; ++isrt) {
                      auto const iall = v[isrt].second;
                      load1 += w8s[iall];
                      state[iall] = state1;
                  } // i
                  assert(load0 + load1 == w8sum && "Maybe failed due to accuracy issues");
                  load_now = i01 ? load1 : load0;
              }

              if (echo > 19) std::printf("# rank#%i assign %g of %g (%.2f %%, target %.2f %%) to %d processes\n", 
                                          rank, load_now, w8sum, load_now*100./w8sum, target_frac*100, nhalf[0]);

              // prepare for the next iteration
              rank_offset += i01*nhalf[0];
              np = nhalf[i01];
              // update the indirection list
              // determine which blocks are still unassigned
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

      if (compute_rank_centers && load_now > 0) {
          // compute the center of weight again, for display
          double cow[] = {0, 0, 0, 0}, w8sum{0};
          for (size_t iuna = 0; iuna < nuna; ++iuna) {
              auto const iall = indirect[iuna];
              auto const *const xyz = xyzw[iall];
              double const w8 = xyz[W];
              assert(w8 == w8s[iall]);
              w8sum += w8;
              // contributes if the weight is positive
              double const w8pos = double(w8 > 0);
              add_product(cow, 3, xyz, w8pos);
              cow[3] += w8pos;
          } // iall
          if (cow[3] > 0) scale(cow, 3, 1./cow[3]);
          set(rank_center[rank], 4, cow);
          if (echo > 13) std::printf("# rank#%i assign %.3f %% center %g %g %g\n", 
                                        rank, w8sum*100/w8sum_all, cow[X], cow[Y], cow[Z]);
      } // compute_rank_centers

      if (echo > 9) std::printf("# rank#%i assign %.3f %%, target %.3f %%\n\n", 
                                   rank, load_now*100/w8sum_all, 100./nprocs);
      load[rank] = load_now;

} // rank


      analyze_load_imbalance(load.data(), nprocs, echo);

      if (compute_rank_centers && nprocs > 1) {
          // analyze positions of rank centers
          double mindist2{9e300}; int ijmin[] = {-1, -1};
          double maxdist2{-1.0};  int ijmax[] = {-1, -1};
          for (int irank = 0; irank < nprocs; ++irank) {
              if (load[irank] > 0) {
                  for (int jrank = 0; jrank < irank; ++jrank) { // self-avoiding triangular loop
                      if (load[jrank] > 0) {
                          auto const dist2 = distance_squared(rank_center[irank], rank_center[jrank]);
                          if (dist2 < mindist2) { mindist2 = dist2; ijmin[0] = irank; ijmin[1] = jrank; }
                          if (dist2 > maxdist2) { maxdist2 = dist2; ijmax[0] = irank; ijmax[1] = jrank; }
                      } // load
                  } // jrank
              } // load
          } // irank
          auto const maxdist = std::sqrt(maxdist2);
          if (echo > 1) std::printf("# shortest distance between centers is %g between rank#%i and #%i, longest is %g\n", 
                                        std::sqrt(mindist2), ijmin[0], ijmin[1], maxdist);
          double const wbin = control::get("load_balancer.bin.width", 0.25), invbin = 1./wbin;
          int const nbin = int(maxdist/wbin) + 1;
          std::vector<uint32_t> hist(nbin, 0);
          int np{0};
          for (int irank = 0; irank < nprocs; ++irank) {
              if (load[irank] > 0) {
                  ++np;
                  for (int jrank = 0; jrank < irank; ++jrank) { // self-avoiding triangular loop
                      if (load[jrank] > 0) {
                          auto const dist2 = distance_squared(rank_center[irank], rank_center[jrank]);
                          auto const dist = std::sqrt(dist2);
                          if (echo > 15) std::printf("# distance-ij is %g\n", dist);
                          int const ibin = dist*invbin; // floor
                          ++hist[ibin];
                      } // load
                  } // jrank
              } // load
          } // irank
          if (echo > 5) {
              double const denom = 1./pow2(std::max(1, np));
              std::printf("## center-distance histogram, bin width %g\n", wbin);
              for (int ibin = 0; ibin < nbin; ++ibin) {
                  std::printf("%g %g\n", ibin*wbin, hist[ibin]*denom);
              } // ibin
              std::printf("\n\n");
          } // echo

      } // compute_rank_centers

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

  inline status_t test_reference_point_cloud(int const nxyz[3], int const echo=9) {
          if (echo > 0) std::printf("# reference point-cloud histogram for %d x %d x %d points\n",
                                                                      nxyz[X], nxyz[Y], nxyz[Z]);
          auto const maxdist_diagonal = std::sqrt(pow2(nxyz[X]) + pow2(nxyz[Y]) + pow2(nxyz[Z]));
          double const wbin = control::get("load_balancer.bin.width", 0.25), invbin = 1./wbin;
          int const nbin = int(maxdist_diagonal/wbin) + 1;
          std::vector<uint32_t> hist(nbin, 0);
          double mindist2{9e300}, maxdist2{-1.0};
          for (int iz = 0; iz < nxyz[Z]; ++iz) {
          for (int iy = 0; iy < nxyz[Y]; ++iy) {
          for (int ix = 0; ix < nxyz[X]; ++ix) {
                  for (int jz = 0; jz < nxyz[Z]; ++jz) {
                  for (int jy = 0; jy < nxyz[Y]; ++jy) {
                  for (int jx = 0; jx < nxyz[X]; ++jx) {
                          double const rnd[] = {0.5*random_between_0_and_1() - 0.25,
                                                0.5*random_between_0_and_1() - 0.25,
                                                0.5*random_between_0_and_1() - 0.25};
                          double const dist2 = pow2(ix - jx + rnd[X])
                                             + pow2(iy - jy + rnd[Y])
                                             + pow2(iz - jz + rnd[Z]);
                          mindist2 = std::min(mindist2, dist2);
                          maxdist2 = std::max(maxdist2, dist2);
                          auto const dist = std::sqrt(dist2);
                          if (echo > 15) std::printf("# distance-ij is %g\n", dist);
                          int const ibin = dist*invbin; // floor
                          ++hist[ibin];
                  }}}
          }}}
          auto const maxdist = std::sqrt(maxdist2);
          if (echo > 1) std::printf("# shortest distance between centers is %g, longest is %g\n", 
                                        std::sqrt(mindist2), maxdist);
          if (echo > 5) {
              double const by_n = 1./(1.*nxyz[X]*nxyz[Y]*nxyz[Z]);
              double const denom = 0.5*pow2(by_n);
              std::printf("## point-distance histogram, bin width %g\n", wbin);
              for (int ibin = 0; ibin < nbin; ++ibin) {
                  auto const radius = ibin*wbin;
                  // number of points inside a sphere shell of from radius r - wbin to r is 4/3*pi*(r^3 - (r - wbin)^3)*rho
                  // so in the limit of small wbin, this becomes 4*pi*r^2*wbin*rho
                  auto const analytical = 2*constants::pi*pow2(radius)*wbin*by_n;
                  std::printf("%g %g %g\n", radius, hist[ibin]*denom, analytical);
              } // ibin
              std::printf("\n\n");
          } // echo
          return 0;
  } // test_reference_point_cloud

#undef weight

  status_t all_tests(int const echo) { 
      status_t stat(0);

      int const nxyz[] = {int(control::get("load_balancer.test.nx", 17.)), // number of blocks
                          int(control::get("load_balancer.test.ny", 19.)),
                          int(control::get("load_balancer.test.nz", 23.))};
      int const nprocs = control::get("load_balancer.test.nprocs", 53.);
      if (echo > 0) std::printf("\n\n# %s start %d x %d x %d = %d with %d MPI processes\n", 
                      __func__, nxyz[X], nxyz[Y], nxyz[Z], nxyz[X]*nxyz[Y]*nxyz[Z], nprocs);

      stat += test_plane_balancer(nprocs, nxyz, echo);
//    stat += test_reference_point_cloud(nxyz, echo);
      return stat;
  } // all_tests

} // namespace load_balancer

#endif // NO_UNIT_TESTS
