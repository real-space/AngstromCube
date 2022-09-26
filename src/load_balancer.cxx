
/*
 * load balancing for A43
 *    main purpose is the distribution of right-hand-side blocks
 *    to MPI processes and GPUs
 *
 */

#ifndef NO_UNIT_TESTS

#include <cstdio> // std::printf
#include <cassert> // assert
#include <cstdint> // size_t, uint32_t
#include <algorithm> // std::max, ::min, ::swap
#include <cmath> // std::ceil, ::pow, ::cbrt
#include <cstdlib> // rand, RAND_MAX

#include "load_balancer.hxx" // ::plane_balancer, ::center_of_weight

#include "status.hxx" // status_t
#include "simple_stats.hxx" // ::Stats<>
#include "constants.hxx" // ::pi
#include "data_view.hxx" // view2D<T>
#include "control.hxx" // ::get

namespace load_balancer {

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


  status_t test_plane_balancer(int const nprocs, int const n[3], int const echo=0) {
      status_t stat(0);

      if (nprocs < 1) return stat;

      auto const nall = (n[X])*size_t(n[Y])*size_t(n[Z]);
      assert(nall <= (size_t(1) << 32) && "uint32_t is not long enough!");

      double w8sum_all{0};
      int constexpr W = 3;
      view2D<float> xyzw(nall, 4, 0.f);
      std::vector<double> w8s(nall, 0);

      for (int iz = 0; iz < n[Z]; ++iz) {
      for (int iy = 0; iy < n[Y]; ++iy) {
      for (int ix = 0; ix < n[X]; ++ix) {
          auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
//        assert(uint32_t(iall) == iall && "uint32_t is not long enough!");
          auto const w8 = 1.f; // weight(ix,iy,iz); // WEIGHTS CAN BE INSERTED HERE
          w8s[iall]    = w8;
          w8sum_all   += w8;
          xyzw(iall,W) = w8;
          xyzw(iall,X) = ix;
          xyzw(iall,Y) = iy;
          xyzw(iall,Z) = iz;
      }}} // ix iy iz

      std::vector<double> load(nprocs, 0.0);
      view2D<double> rank_center(nprocs, 4, 0.0);
      bool constexpr compute_rank_centers = true;
      uint16_t constexpr no_owner = (1 << 16) - 1;
      std::vector<uint16_t> owner_rank(nall, no_owner);

      if (echo > 0) std::printf("# %s: distribute %g blocks to %d processes\n\n", __func__, w8sum_all, nprocs);

      for (int rank = 0; rank < nprocs; ++rank) {
          load[rank] = plane_balancer(nprocs, rank, nall, xyzw, w8s.data(), w8sum_all, echo
                                            , rank_center[rank], owner_rank.data());
      } // rank

      analyze_load_imbalance(load.data(), nprocs, echo);

      if (compute_rank_centers && nprocs > 1) {
          // analyze positions of rank centers
          double mindist2{9e300}; int ijmin[] = {-1, -1};
          double maxdist2{-1.0};  int ijmax[] = {-1, -1};
          simple_stats::Stats<> st2, st1;
          for (int irank = 0; irank < nprocs; ++irank) {
              if (load[irank] > 0) {
                  for (int jrank = 0; jrank < nprocs; ++jrank) { // self-avoiding triangular loop
                      if (load[jrank] > 0) {
                          auto const dist2 = distance_squared(rank_center[irank], rank_center[jrank]);
                          if (dist2 > 0 && dist2 < mindist2) { mindist2 = dist2; ijmin[0] = irank; ijmin[1] = jrank; }
                          if (dist2 > maxdist2) { maxdist2 = dist2; ijmax[0] = irank; ijmax[1] = jrank; }
                          auto const dist = std::sqrt(dist2);
                          st2.add(dist2);
                          st1.add(dist);
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
          int np{0}; // counter for the number of processes with a non-zero load
          for (int irank = 0; irank < nprocs; ++irank) {
              if (load[irank] > 0) {
                  ++np;
                  for (int jrank = 0; jrank < nprocs; ++jrank) { // self-avoiding triangular loop
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


      // check masks
      if (1) {
          int strange{0};
          for (int iz = 0; iz < n[Z]; ++iz) {
            for (int iy = 0; iy < n[Y]; ++iy) {
              for(int ix = 0; ix < n[X]; ++ix) {
                  auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
                  auto const owner = owner_rank[iall];
                  strange += (no_owner == owner); // under-assignement
                  if (no_owner == owner) warn("work item %d %d %d has not been assigned to any rank", ix,iy,iz);
              } // ix
            } // iy
          } // iz
          if (strange) warn("strange: %d under-assignments", strange);
      }

      if (1 == n[Z] && echo > 5) {
          std::printf("\n# visualize plane balancer %d x %d on %d processes:%s", n[Y], n[X], nprocs,
                              nprocs > 64 ? " (symbols are not unique!)" : "");
          int constexpr iz = 0;
          char const chars[] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ<>";
          for (int iy = 0; iy < n[Y]; ++iy) {
              std::printf("\n# ");
              for(int ix = 0; ix < n[X]; ++ix) {
                  auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
                  auto const owner = owner_rank[iall];
                  auto const c = chars[owner & 0x3f];
                  std::printf("%c", c);
              } // ix
          } // iy
          std::printf("\n#\n\n");
          //
          // ToDo: to export the Voronoi diagrams, we need to access the plane normals
          //       and plane parameters at which the plane separates the two processes
          //
      } // visualize

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

      int const nxyz[] = {int(control::get("load_balancer.test.nx", 17.)), // number of blocks
                          int(control::get("load_balancer.test.ny", 19.)),
                          int(control::get("load_balancer.test.nz", 23.))};
      auto const nprocs = int(control::get("load_balancer.test.nprocs", 53.));
      if (echo > 0) std::printf("\n\n# %s start %d x %d x %d = %d with %d MPI processes\n", 
                      __func__, nxyz[X], nxyz[Y], nxyz[Z], nxyz[X]*nxyz[Y]*nxyz[Z], nprocs);

      stat += test_plane_balancer(nprocs, nxyz, echo);
//    stat += test_reference_point_cloud(nxyz, echo);
      return stat;
  } // all_tests

} // namespace load_balancer

#endif // NO_UNIT_TESTS
