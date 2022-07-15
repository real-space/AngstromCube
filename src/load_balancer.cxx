
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

#include "load_balancer.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
// #include "mpi_parallel.hxx" // ...
#include "simple_stats.hxx" // ::Stats<>
#include "data_view.hxx" // view4D<T>
#include "inline_math.hxx" // pow2, scale, align<nBits>
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

  inline status_t test_bisection_balancer(int const nprocs, int const n[3], int const echo=0) {
      status_t stat(0);

      // iterative process for balancing the computational load onto MPI processes.
      // The load can be the number and costs of compute tasks, memory consumption, etc.

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


  inline status_t test_diffusion_balancer(int const nprocs, int const n[3], int const echo=0) {
      status_t stat(0);
      int constexpr X=0, Y=1, Z=2;

      // iterative process for balancing the computational load onto MPI processes.
      // The load can be the number and costs of compute tasks, memory consumption, etc.

      if (nprocs > 1) return stat; // no surface
      auto const blockvolume = n[X]*size_t(n[Y])*size_t(n[Z]);
      auto const blocksperproc = blockvolume/double(nprocs);
      auto const cube_root = std::cbrt(blocksperproc);
      if (echo > 0) std::printf("# division of %d x %d x %d = %ld blocks onto %d processes, expected %.2f^3 = %.2f blocks per process\n",
                                  n[X], n[Y], n[Z], blockvolume, nprocs, cube_root, blocksperproc);
      assert(blockvolume > 0);

      view2D<float> pos(nprocs, 4, 0.f); // centroid positions
      view3D<int32_t> owner(n[Z], n[Y], n[X], -1); // owning ranks for each block
      view2D<float> load_data(2, align<3>(nprocs), 0.f);
      float *load = load_data[0], *last_load = load_data[1];

      { // scope: create an initial distribution using a Voronoi construction

          // approximately factorize nprocs
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

          if (echo > 0) std::printf("# division into %d x %d x %d = %d processes is closest to %d processes\n",
                                            npd[X], npd[Y], npd[Z], nprocs_more,             nprocs);


          // skip positions that are too much

          for (int myrank = 0; myrank < nprocs; ++myrank) {

              int const iproc_more = (myrank*nprocs_more)/nprocs;     assert(iproc_more < nprocs_more);
              int const coords[] = {iproc_more % npd[X], (iproc_more / npd[X]) % npd[Y], iproc_more / (npd[X]*npd[Y])};

              if (echo > 29) std::printf("# myrank=%i pseudorank=%i of %d, coords %d %d %d\n", myrank, iproc_more, nprocs_more, coords[X], coords[Y], coords[Z]);
              // place centroid on the map
              float const position[] = {(coords[X] + .5f)*by_npd[X],
                                        (coords[Y] + .5f)*by_npd[Y],
                                        (coords[Z] + .5f)*by_npd[Z]};

              set(pos[myrank], 3, position);
              if (echo > 9) std::printf("# myrank=%i initial coords %.1f %.1f %.1f\n", myrank, position[X], position[Y], position[Z]);

          } // myrank

          set(load, nprocs, 0.f); // clear
          for (int iz = 0; iz < n[Z]; ++iz) {
          for (int iy = 0; iy < n[Y]; ++iy) {
          for (int ix = 0; ix < n[X]; ++ix) {
              float largest{-3e33f};
              int l_rank{-1};
              for (int myrank = 0; myrank < nprocs; ++myrank) {
    //            auto const p = pop(myrank,iz,iy,ix);
                  auto const dist2 = pow2(ix - pos(myrank,X))
                                   + pow2(iy - pos(myrank,Y))
                                   + pow2(iz - pos(myrank,Z));
                  auto const p = -dist2;
                  if (p > largest) { largest = p; l_rank = myrank; }
              } // myrank
              assert(l_rank >= 0);
              owner(iz,iy,ix) = l_rank;
              float const w8 = 1; // weight(iz,iy,ix);
              load[l_rank] += w8;
          }}} // ix iy iz

          simple_stats::Stats<> st(0);
          for (int myrank = 0; myrank < nprocs; ++myrank) {
              st.add(load[myrank]);
              if (echo > 9) std::printf("# myrank=%i owns %g\n", myrank, load[myrank]);
          } // myrank
          if (echo > 0) std::printf("# %ld processes own %g blocks, per process [%g, %.2f +/- %.2f, %g]\n",
                                        st.tim(), st.sum(), st.min(), st.mean(), st.dev(), st.max());

          auto const load_balance = st.mean()/st.max(), load_imbalance = st.max()/st.min();
          if (echo > 0) std::printf("# load balance %.2f %%, imbalance %.2f\n", load_balance*100, load_imbalance);

      } // scope: initial distribution

      // now we have an initial distribution, certainly not perfectly balanced


      // compute forces onto the centroids
      // 1.) force/energy contribution is imbalance
      // 2.) force/energy contribution is surface

      view3D<float> f_surface(n[Z], n[Y], n[X], 0.f);

      double total_load{0};
      for (int iz = 0; iz < n[Z]; ++iz) {
      for (int iy = 0; iy < n[Y]; ++iy) {
      for (int ix = 0; ix < n[X]; ++ix) {
          float const w8 = 1; // weight(iz,iy,ix);
          total_load += w8;
      }}} // ix iy iz

      float const load_average = total_load/nprocs;

      double last_E_total;
      float last_maxfrc;

      for (int iteration = 0; iteration < 9; ++iteration) {

          std::swap(last_load, load); // swap pointers

          // E_imbalance ~ sum_r (load[r] - <load>)^2 = sum_r load[r]^2 - constant // square deviation
          // <load> = ngrid/nprocs is constant
          // load[r] = sum_xyz (owner[xyz] == r)
          //
          // E_imbalance = sum_r ( sum_xyz (owner[xyz] == r) )^2
          // E_surface   = sum_xyz (owner[xyz] != owner[xyz +/- 1])^2

          set(load, nprocs, 0.f); // clear
          for (int iz = 0; iz < n[Z]; ++iz) {
          for (int iy = 0; iy < n[Y]; ++iy) {
          for (int ix = 0; ix < n[X]; ++ix) {
              auto const owner_rank = owner(iz,iy,ix);    assert(owner_rank >= 0); assert(owner_rank < nprocs);
              float const w8 = 1; // weight(iz,iy,ix);

              load[owner_rank] += w8;
          }}} // ix iy iz

          double E_surface{0};
          for (int iz = 0; iz < n[Z]; ++iz) {
          for (int iy = 0; iy < n[Y]; ++iy) {
          for (int ix = 0; ix < n[X]; ++ix) {
              auto const owner_rank = owner(iz,iy,ix);
              float const w8 = 1; // weight(iz,iy,ix);

              int number_of_different_neighbors{0};
              float stress{0};
#define STAR
#ifdef  STAR
              int const ixyz[] = {ix, iy, iz};
              int jxyz[] = {ix, iy, iz}; // non-const
              // loop over Cartesian neighborhood star (6 points)
              for (int d = 0; d < 3; ++d) {
                  for (int mp1 = -1; mp1 <= 1; mp1 += 2) {
                      jxyz[d] = ixyz[d] + mp1;
                      if (jxyz[d] >= 0 && jxyz[d] < n[d]) {
                          auto const neighbor_rank = owner(jxyz[Z],jxyz[Y],jxyz[X]);
                          int  const differ = (owner_rank != neighbor_rank);
                          number_of_different_neighbors += differ;
                          stress += differ*(load[owner_rank] - load_average)
                                          *(load[neighbor_rank] - load_average); // we can combine these triple loops with the previous ones if we use last_load instead of load
                      } // in bounds
                  } // -/+1
                  jxyz[d] = ixyz[d]; // reset
              } // d
#else  // STAR
              // loop over 3^3-1 surrounding points, weighted with 4, 2 and 1
              // for the face(100), edge(110) and vertex(111) directions, respectively.
              for (int jz = iz-1; jz <= iz+1; ++jz) {  if (jz >= 0 && jz < n[Z]) {
              for (int jy = iy-1; jy <= iy+1; ++jy) {  if (jy >= 0 && jy < n[Y]) {
              for (int jx = ix-1; jx <= ix+1; ++jx) {  if (jx >= 0 && jx < n[X]) {
                  int const d2 = pow2(jx - ix) + pow2(jy - iy) + pow2(jz - iz);      assert(d2 <= 3);
                  number_of_different_neighbors += (owner_rank != owner(jz,jy,jx)) * (1 << (3 - d2));
              }}} }}} // jx jy jz
#endif // STAR
              E_surface += number_of_different_neighbors*w8;
              f_surface(iz,iy,ix) = stress;
          }}} // ix iy iz

          double E_imbalance{0}; // energy contribution
          for (int32_t owner_rank = 0; owner_rank < nprocs; ++owner_rank) {
              E_imbalance += pow2(load[owner_rank] - load_average);
          } // owner_rank
          auto const E_total = E_imbalance + E_surface;
          if (echo > 0) std::printf("# E_total= %g = %g + %g (E_imbalance + E_surface)\n",
                                      E_total,               E_imbalance,  E_surface);

          // now compute the gradients to suggest a change
          // d E_imbalance
          // ------------- = 2 load[r]
          //  d load[r]

          //  d load[r]
          // ------------ = 1 * (owner[xyz] != r) - 1 * (owner[xyz] == r)
          // d owner[xyz]

          float maxfrc{-9e9}, minfrc{9e9};
          int imax[3] = {-1, -1, -1};
          for (int iz = 0; iz < n[Z]; ++iz) {
          for (int iy = 0; iy < n[Y]; ++iy) {
          for (int ix = 0; ix < n[X]; ++ix) {
              auto const owner_rank = owner(iz,iy,ix);
              float const frc_imbalance = 2*(load[owner_rank] - load_average);
              auto const frc_total = frc_imbalance + f_surface(iz,iy,ix);
              if (frc_total > maxfrc) {
                  maxfrc = frc_total;
                  imax[X] = ix;  imax[Y] = iy;  imax[Z] = iz;
              } // is larger
              if (frc_total < minfrc) minfrc = frc_total;
          }}} // ix iy iz
          if (echo > 0) std::printf("# smallest force is %g\n", minfrc);

          auto const maxfrc_surf = f_surface(imax[Z],imax[Y],imax[X]);
          if (echo > 0) std::printf("# largest force is %g = %g + %g (imbalance + surface) at %i %i %i\n",
                                    maxfrc, maxfrc - maxfrc_surf, maxfrc_surf, imax[X], imax[Y], imax[Z]);


          if (iteration > 0) {
              if (echo > 0) std::printf("# iteration=%i new E_total= %g, last E_total= %g, change= %g, new frc= %g, last frc= %g\n",
                                            iteration, E_total, last_E_total, E_total - last_E_total, maxfrc, last_maxfrc);
          }


          // suggest to change towards a rank out of the neighborhood of imax
          int32_t suggested_rank{-1};
          { // scope
              int32_t neigh[26];
              int i6{0};
              int const ixyz[] = {imax[X], imax[Y], imax[Z]};
              auto const owner_rank = owner(ixyz[Z],ixyz[Y],ixyz[X]);
#ifdef  STAR
              int jxyz[] = {ixyz[X], ixyz[Y], ixyz[Z]}; // non-const
              // loop over Cartesian neighborhood star (6 points)
              for (int d = 0; d < 3; ++d) {
                  for (int mp1 = -1; mp1 <= 1; mp1 += 2) {
                      jxyz[d] = ixyz[d] + mp1;
                      if (jxyz[d] >= 0 && jxyz[d] < n[d]) {
                          auto const neighbor_rank = owner(jxyz[Z],jxyz[Y],jxyz[X]);
                          if (owner_rank != neighbor_rank) {
                              neigh[i6] = neighbor_rank;
                              ++i6;
                          } // differ
                      } // in bounds
                  } // -/+1
                  jxyz[d] = ixyz[d]; // reset
              } // d
#else  // STAR
  #error "not implemented"
#endif // STAR
              int const n6 = i6;
              assert(n6 > 0 && "we expect the largest force to appear on a surface");
              // see how many neighbors have the same rank
              uint8_t count[26] = {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0};
              int32_t ranks[26];
              int ndiff{0};
              for (int i6 = 0; i6 < n6; ++i6) {
                  auto const neighbor_rank = neigh[i6];
                  int listed{-1};
                  for (int id = 0; id < ndiff; ++id) {
                      if (neighbor_rank == ranks[id]) {
                          ++count[id];
                          listed = id;
                      }
                  } // id
                  if (-1 == listed) { // not listed
                      ranks[ndiff] = neighbor_rank;
                      count[ndiff] = 1;
                      ++ndiff;
                  } // not listed
              } // i6

              // find the highest count
              int maxcount{-1};
              for (int id = 0; id < ndiff; ++id) {
                  if (count[id] > maxcount) {
                      maxcount = count[id];
                      suggested_rank = ranks[id];
                  } // is larger
              } // id

              if (echo > 0) std::printf("# suggest to change owning rank for %i %i %i from %i to %i\n",
                                                imax[X], imax[Y], imax[Z], owner_rank, suggested_rank);
              assert(suggested_rank > -1);
          } // scope

          // do the suggested change
          owner(imax[Z],imax[Y],imax[X]) = suggested_rank;

          last_E_total = E_total;
          last_maxfrc  = maxfrc;

      } // while

      // compute new centroids
      for (int iz = 0; iz < n[Z]; ++iz) {
      for (int iy = 0; iy < n[Y]; ++iy) {
      for (int ix = 0; ix < n[X]; ++ix) {
          auto const owner_rank = owner(iz,iy,ix);
          float const w8 = 1; // weight(iz,iy,ix);
          int const xyz1[] = {ix, iy, iz, 1};
          add_product(pos[owner_rank], 4, xyz1, w8);
      }}} // ix iy iz

      for (int myrank = 0; myrank < nprocs; ++myrank) {
          auto *const position = pos[myrank];
          auto const denom = position[3];
          if (denom > 0) {
              scale(position, 3, 1.f/denom);
          } else {
              set(position, 4, 0.f);
              assert(0 == denom);
          } // owns some
          if (echo > 9) std::printf("# myrank=%i   final coords %.1f %.1f %.1f   %g\n", myrank, position[X], position[Y], position[Z], position[3]);
      } // myrank

      return stat;
  } // test_diffusion_balancer


  inline status_t test_RobinHood_balancer(int const nprocs, int const n[3], int const echo=0) {
      status_t stat(0);
      int constexpr X=0, Y=1, Z=2;

      // iterative process for balancing the computational load onto MPI processes.
      // The load can be the number and costs of compute tasks, memory consumption, etc.

      if (nprocs < 2) return stat; // no surface
      auto const blockvolume = n[X]*size_t(n[Y])*size_t(n[Z]);
      auto const blocksperproc = blockvolume/double(nprocs);
      auto const cube_root = std::cbrt(blocksperproc);
      if (echo > 0) std::printf("# division of %d x %d x %d = %ld blocks onto %d processes, expected %.2f^3 = %.2f blocks per process\n",
                                  n[X], n[Y], n[Z], blockvolume, nprocs, cube_root, blocksperproc);
      assert(blockvolume > 0);

      view2D<float> pos(nprocs, 4, 0.f); // centroid positions
      view3D<int32_t> owner(n[Z], n[Y], n[X], -1); // owning ranks for each block
//    view3D<float> weight(n[Z], n[Y], n[X], 1.f); // can be modified
      view2D<float> load_data(2, align<3>(nprocs), 0.f);
      float *load = load_data[0], *last_load = load_data[1];


      { // scope: create an initial distribution using a Voronoi construction

          // approximately factorize nprocs
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

          if (echo > 0) std::printf("# division into %d x %d x %d = %d processes is closest to %d processes\n",
                                            npd[X], npd[Y], npd[Z], nprocs_more,             nprocs);


          // skip positions that are too much

          for (int myrank = 0; myrank < nprocs; ++myrank) {

              int const iproc_more = (myrank*nprocs_more)/nprocs;     assert(iproc_more < nprocs_more);
              int const coords[] = {iproc_more % npd[X], (iproc_more / npd[X]) % npd[Y], iproc_more / (npd[X]*npd[Y])};

              if (echo > 29) std::printf("# myrank=%i pseudorank=%i of %d, coords %d %d %d\n", myrank, iproc_more, nprocs_more, coords[X], coords[Y], coords[Z]);
              // place centroid on the map
              float const position[] = {(coords[X] + .5f)*by_npd[X],
                                        (coords[Y] + .5f)*by_npd[Y],
                                        (coords[Z] + .5f)*by_npd[Z]};

              set(pos[myrank], 3, position);
              if (echo > 9) std::printf("# myrank=%i initial coords %.1f %.1f %.1f\n", myrank, position[X], position[Y], position[Z]);

          } // myrank

          set(load, nprocs, 0.f); // clear
          for (int iz = 0; iz < n[Z]; ++iz) {
          for (int iy = 0; iy < n[Y]; ++iy) {
          for (int ix = 0; ix < n[X]; ++ix) {
              float largest{-3e33f};
              int rank_of_largest{-1};
              for (int myrank = 0; myrank < nprocs; ++myrank) {
    //            auto const p = pop(myrank,iz,iy,ix);
                  auto const dist2 = pow2(ix - pos(myrank,X))
                                   + pow2(iy - pos(myrank,Y))
                                   + pow2(iz - pos(myrank,Z));
                  auto const p = -dist2;
                  if (p > largest) { largest = p; rank_of_largest = myrank; }
              } // myrank
              assert(rank_of_largest >= 0);
              owner(iz,iy,ix) = rank_of_largest;
              float const w8 = 1; // weight(iz,iy,ix);
              load[rank_of_largest] += w8;
          }}} // ix iy iz

          simple_stats::Stats<> st(0);
          for (int myrank = 0; myrank < nprocs; ++myrank) {
              st.add(load[myrank]);
              if (echo > 9) std::printf("# myrank=%i owns %g\n", myrank, load[myrank]);
          } // myrank
          if (echo > 0) std::printf("# %ld processes own %g blocks, per process [%g, %.2f +/- %.2f, %g]\n",
                                        st.tim(), st.sum(), st.min(), st.mean(), st.dev(), st.max());

          auto const load_balance = st.mean()/st.max(), load_imbalance = st.max()/st.min();
          if (echo > 0) std::printf("# load balance %.2f %%, imbalance %.2f\n", load_balance*100, load_imbalance);

      } // scope: initial distribution

      // now we have an initial distribution, certainly not perfectly balanced


      // compute forces onto the centroids
      // 1.) force/energy contribution is imbalance
      // 2.) force/energy contribution is surface

      double total_load{0};
      for (int iz = 0; iz < n[Z]; ++iz) {
      for (int iy = 0; iy < n[Y]; ++iy) {
      for (int ix = 0; ix < n[X]; ++ix) {
          float const w8 = 1; // weight(iz,iy,ix);
          total_load += w8;
      }}} // ix iy iz

      float const load_average = total_load/nprocs;
      double last_E_total;

      for (int iteration = 0; iteration < 9; ++iteration) {

          std::swap(last_load, load); // swap pointers

          // E_imbalance ~ sum_r (load[r] - <load>)^2 = sum_r load[r]^2 - constant // square deviation
          // <load> = ngrid/nprocs is constant
          // load[r] = sum_xyz (owner[xyz] == r)
          //
          // E_imbalance = sum_r ( sum_xyz (owner[xyz] == r) )^2
          // E_surface   = sum_xyz (owner[xyz] != owner[xyz +/- 1])^2

          set(load, nprocs, 0.f); // clear
          for (int iz = 0; iz < n[Z]; ++iz) {
          for (int iy = 0; iy < n[Y]; ++iy) {
          for (int ix = 0; ix < n[X]; ++ix) {
              auto const owner_rank = owner(iz,iy,ix);    assert(owner_rank >= 0); assert(owner_rank < nprocs);
              float const w8 = 1; // weight(iz,iy,ix);

              load[owner_rank] += w8;
          }}} // ix iy iz

          double E_surface{0};
          float maxstress{-9e9}, minstress{9e9};
          int32_t im[4] = {-1, -1, -1, -1}; // indices of the extremum
          for (int iz = 0; iz < n[Z]; ++iz) {
          for (int iy = 0; iy < n[Y]; ++iy) {
          for (int ix = 0; ix < n[X]; ++ix) {
              auto const owner_rank = owner(iz,iy,ix);
              float const w8 = 1; // weight(iz,iy,ix);

              int number_of_different_neighbors{0};

              int const ixyz[] = {ix, iy, iz};
              int jxyz[] = {ix, iy, iz}; // non-const
              auto const owner_wants_to_give = load[owner_rank] - load_average;
                         // we can combine these triple loops with the previous ones if we use last_load instead of load
              // loop over Cartesian neighborhood star (6 points)
              for (int d = 0; d < 3; ++d) {
                  for (int mp1 = -1; mp1 <= 1; mp1 += 2) {
                      jxyz[d] = ixyz[d] + mp1;
                      if (jxyz[d] >= 0 && jxyz[d] < n[d]) { // boundary
                          auto const neighbor_rank = owner(jxyz[Z],jxyz[Y],jxyz[X]);
                          int  const differ = (owner_rank != neighbor_rank); // is 1 or 0
                          auto const neighbor_wants_to_take = load_average - load[neighbor_rank];
                          auto const stress = differ*owner_wants_to_give*neighbor_wants_to_take;
                          if (stress > maxstress) {
                              maxstress = stress;
                          }
                          if (stress < minstress) {
                              minstress = stress;
                              im[X] = ix; im[Y] = iy; im[Z] = iz;
                              im[3] = neighbor_rank;
//                            im[4] = owner_rank;
//                            im[5] = d*2 + (1 + mp1)/2;
                          }
                          number_of_different_neighbors += differ;
                      } // boundary
                  } // -/+1
                  jxyz[d] = ixyz[d]; // reset
              } // d

              E_surface += number_of_different_neighbors*w8;
          }}} // ix iy iz

          auto const owner_of_largest = owner(im[Z],im[Y],im[X]);
          auto const neighbor_of_largest = im[3];
          if (echo > 0) std::printf("\n# smallest stress is %g, largest stress is %g at [%i %i %i], rank %i --> rank %i\n", 
                                     minstress, maxstress,  im[X], im[Y], im[Z], owner_of_largest, neighbor_of_largest);
          auto const weight_of_largest = 1.f; // weight(im[Z],im[Y],im[X]);

          double E_imbalance{0}, suggested_E_imbalance{0}; // energy contribution
          set(last_load, nprocs, load);
          last_load[owner_of_largest]    -= weight_of_largest;
          last_load[neighbor_of_largest] += weight_of_largest;
          for (int32_t owner_rank = 0; owner_rank < nprocs; ++owner_rank) {
              E_imbalance += pow2(load[owner_rank] - load_average);
              suggested_E_imbalance += pow2(last_load[owner_rank] - load_average);
          } // owner_rank
          auto const E_total = E_imbalance + E_surface;
          if (echo > 0) std::printf("# E_total= %g = %g + %g (E_imbalance + E_surface), suggested E_imbalance %g\n",
                                       E_total,               E_imbalance,  E_surface,  suggested_E_imbalance);


          if (iteration > 0) {
              if (echo > 0) std::printf("# iteration=%i new E_total= %g, last E_total= %g, change= %g\n",
                                            iteration, E_total, last_E_total, E_total - last_E_total);
          }

          // suggest to change towards a rank out of the neighborhood of imax
          int32_t const suggested_rank = im[3];

          // do the suggested change
          owner(im[Z],im[Y],im[X]) = suggested_rank;

          last_E_total = E_total;

      } // while

      // compute new centroids
      for (int iz = 0; iz < n[Z]; ++iz) {
      for (int iy = 0; iy < n[Y]; ++iy) {
      for (int ix = 0; ix < n[X]; ++ix) {
          auto const owner_rank = owner(iz,iy,ix);
          float const w8 = 1; // weight(iz,iy,ix);
          int const xyz1[] = {ix, iy, iz, 1};
          add_product(pos[owner_rank], 4, xyz1, w8);
      }}} // ix iy iz

      for (int myrank = 0; myrank < nprocs; ++myrank) {
          auto *const position = pos[myrank];
          auto const denom = position[3];
          if (denom > 0) {
              scale(position, 3, 1.f/denom);
          } else {
              set(position, 4, 0.f);
              assert(0 == denom);
          } // owns some
          if (echo > 9) std::printf("# myrank=%i   final coords %.1f %.1f %.1f   %g\n", myrank, position[X], position[Y], position[Z], position[3]);
      } // myrank

      return stat;
  } // test_RobinHood_balancer



  inline status_t test_pressure_balancer(int const nprocs, int const n[3], int const echo=0) {
      status_t stat(0);

      // iterative process for balancing the computational load onto MPI processes.
      // The load can be the number and costs of compute tasks, memory consumption, etc.

      if (nprocs < 2) return stat;
      auto const blockvolume = n[X]*size_t(n[Y])*size_t(n[Z]);
      auto const blocksperproc = blockvolume/double(nprocs);
      auto const cube_root = std::cbrt(blocksperproc);
      if (echo > 0) std::printf("# division of %d x %d x %d = %ld blocks onto %d processes, expected %.2f^3 = %.2f blocks per process\n",
                                  n[X], n[Y], n[Z], blockvolume, nprocs, cube_root, blocksperproc);
      assert(blockvolume > 0);

      view2D<float> pos(nprocs, 4, 0.f); // centroid positions
      view3D<int32_t> owner(n[Z], n[Y], n[X], -1); // owning ranks for each block
      auto const nprocs_aligned = align<3>(nprocs);
      view2D<float> load_data(2, nprocs_aligned, 0.f);
      float *load = load_data[0], *new_load = load_data[1];

      { // scope: create an initial distribution using a Voronoi construction

          // approximately factorize nprocs
          auto const inv_cbrt = 1./cube_root;
          int npd[3];
          for (int d = 0; d < 3; ++d) {
              assert(n[d] > 0);
              auto const npdf = n[d]*inv_cbrt;
              npd[d] = std::ceil(npdf);
              if (echo > 2) std::printf("# in %c-direction [%d, %.2f, %d] blocks per process\n", 'x' + d, npd[d] - 1, npdf, npd[d]);
              assert(npd[d] > 0);
          } // d
          int const nprocs_more = npd[X] * npd[Y] * npd[Z];
          assert(nprocs_more > 0);
          float const by_npd[] = {n[X]/float(npd[X]), n[Y]/float(npd[Y]), n[Z]/float(npd[Z])};

          if (echo > 0) std::printf("# division into %d x %d x %d = %d processes is closest to %d processes\n",
                                            npd[X], npd[Y], npd[Z], nprocs_more,             nprocs);


          // skip positions that are too much

          for (int myrank = 0; myrank < nprocs; ++myrank) {

              int const iproc_more = (myrank*nprocs_more)/nprocs;     assert(iproc_more < nprocs_more);
              int const coords[] = {iproc_more % npd[X], (iproc_more / npd[X]) % npd[Y], iproc_more / (npd[X]*npd[Y])};

              if (echo > 29) std::printf("# myrank=%i pseudorank=%i of %d, coords %d %d %d\n", myrank, iproc_more, nprocs_more, coords[X], coords[Y], coords[Z]);
              // place centroid on the map
              float const position[] = {(coords[X] + .5f)*by_npd[X],
                                        (coords[Y] + .5f)*by_npd[Y],
                                        (coords[Z] + .5f)*by_npd[Z]};
              if (echo > 9) std::printf("# myrank=%i initial coords %9.1f %9.1f %9.1f\n", myrank, position[X], position[Y], position[Z]);
              set(pos[myrank], 3, position);

          } // myrank

          set(load, nprocs, 0.f); // clear
          for (int iz = 0; iz < n[Z]; ++iz) {
          for (int iy = 0; iy < n[Y]; ++iy) {
          for (int ix = 0; ix < n[X]; ++ix) {
              float shortest_distance2{3e33f};
              int32_t owner_rank{-1};
              for (int myrank = 0; myrank < nprocs; ++myrank) {
                  auto const dist2 = pow2(ix - pos(myrank,X))
                                   + pow2(iy - pos(myrank,Y))
                                   + pow2(iz - pos(myrank,Z));
                  if (dist2 < shortest_distance2) {
                      shortest_distance2 = dist2;
                      owner_rank = myrank;
                  } // is shorter
              } // myrank
              assert(owner_rank >= 0);
              owner(iz,iy,ix) = owner_rank;

              // compute the initial load
              auto const w8 = weight(iz,iy,ix);
              load[owner_rank] += w8;
          }}} // ix iy iz

          simple_stats::Stats<> st(0);
          for (int myrank = 0; myrank < nprocs; ++myrank) {
              st.add(load[myrank]);
              if (echo > 9) std::printf("# myrank=%i owns %g\n", myrank, load[myrank]);
          } // myrank
          if (echo > 0) std::printf("# %ld processes own %g blocks, per process [%g, %.2f +/- %.2f, %g]\n",
                                        st.tim(), st.sum(), st.min(), st.mean(), st.dev(), st.max());

          auto const load_balance = st.mean()/st.max(), load_imbalance = st.max()/st.min();
          if (echo > 0) std::printf("# load balance %.2f %%, imbalance %.2f\n", load_balance*100, load_imbalance);

      } // scope: initial distribution

      // now we have an initial distribution, certainly not perfectly balanced


      // compute forces onto the centroids
      // 1.) force/energy contribution is imbalance
      // 2.) force/energy contribution is surface


      double total_load{0};
      for (int iz = 0; iz < n[Z]; ++iz) {
      for (int iy = 0; iy < n[Y]; ++iy) {
      for (int ix = 0; ix < n[X]; ++ix) {
          auto const w8 = weight(iz,iy,ix);
          total_load += w8;
      }}} // ix iy iz

      float const load_average = total_load/nprocs;

      double E_imbalance_initial;
      double new_E_imbalance;
      double last_E_total;

      view2D<float> contact(nprocs, nprocs_aligned); // count on how many surface points two processes are in contact

      int stop{0};
      for (int iteration = 0; iteration < 999 && 0 == stop; ++iteration) {
          int const echo_it = 7 + 4*(iteration > 0);

          std::swap(new_load, load); // swap pointers

          for (int32_t rank = 0; rank < nprocs; ++rank) {
              set(contact[rank], contact.stride(), 0.f); // clear
          } // rank
          set(load, nprocs, 0.f); // clear
          std::vector<uint32_t> n_surf(nprocs, 0);
          for (int iz = 0; iz < n[Z]; ++iz) {
          for (int iy = 0; iy < n[Y]; ++iy) {
          for (int ix = 0; ix < n[X]; ++ix) {
              auto const owner_rank = owner(iz,iy,ix);    assert(owner_rank >= 0); assert(owner_rank < nprocs);
              auto const w8        = weight(iz,iy,ix);

              load[owner_rank] += w8;

              int number_of_different_neighbors{0};

              int const ixyz[] = {ix, iy, iz};
              int jxyz[] = {ix, iy, iz}; // non-const
              // loop over Cartesian neighborhood star (6 points)
              for (int d = 0; d < 3; ++d) {
                  for (int mp1 = -1; mp1 <= 1; mp1 += 2) {
                      jxyz[d] = ixyz[d] + mp1;
                      if (jxyz[d] >= 0 && jxyz[d] < n[d]) {
                          auto const neighbor_rank = owner(jxyz[Z],jxyz[Y],jxyz[X]);
                          int  const differ = (owner_rank != neighbor_rank);
                          number_of_different_neighbors += differ;
                          contact(owner_rank,neighbor_rank) += w8;
                      } // in bounds
                  } // -/+1
                  jxyz[d] = ixyz[d]; // reset
              } // d

              n_surf[owner_rank] += number_of_different_neighbors;
          }}} // ix iy iz

          double E_imbalance{0}; // energy contribution
          for (int32_t owner_rank = 0; owner_rank < nprocs; ++owner_rank) {
              E_imbalance += pow2(load[owner_rank] - load_average);
          } // owner_rank

          if (iteration > 0) {
              if (echo > echo_it) std::printf("# iteration=%i new E_total= %g, last E_total= %g, change= %g\n",
                                            iteration, E_imbalance, last_E_total, E_imbalance - last_E_total);
          } // iteration
          if (0 == iteration) E_imbalance_initial = E_imbalance;

          // suggest a giving and a receiving rank
          int32_t giver{-1}, taker{-1};
          int imax[] = {-1, -1, -1};
          { // scope
              float maxstress{-9e9}, minstress{9e9};
              for (int32_t irank = 0; irank < nprocs; ++irank) {
                  auto const wants2give = load[irank] - load_average;
                  if (wants2give >= 0) {
                      for (int32_t jrank = 0; jrank < nprocs; ++jrank) {
                          if (irank != jrank) {
                              auto const total_contact = contact(irank,jrank);
                              if (total_contact > 0) {
                                  auto const wants2take = load_average - load[jrank];
                                  if (wants2take >= 0) {
                                      float const stress = total_contact*wants2give*wants2take;
                                      if (stress > maxstress) {
                                          maxstress = stress;
                                          giver = irank;
                                          taker = jrank;
                                      }
                                      if (stress < minstress) {
                                          minstress = stress;
                                      }
                                  } // wants2take
                              } // total_contact
                          } // i != j
                      } // jrank
                  } // wants2give
              } // irank
          } // scope
          if (echo > echo_it) std::printf("# iteration=%i suggested transfer rank %i --> rank %i\n", iteration, giver, taker);
          assert(giver != taker);

          // now find where exactly this transfer should happen
          float maxtension{-9e9}, mintension{9e9};
          {
              for (int iz = 0; iz < n[Z]; ++iz) {
              for (int iy = 0; iy < n[Y]; ++iy) {
              for (int ix = 0; ix < n[X]; ++ix) {
                  auto const owner_rank = owner(iz,iy,ix);
                  if (owner_rank == giver) {
//                    auto const w8 = weight(iz,iy,ix); // not needed
                      int tension_giver{0}, tension_taker{0};

                      int const ixyz[] = {ix, iy, iz};
                      int jxyz[] = {ix, iy, iz}; // non-const
                      // loop over Cartesian neighborhood star (6 points)
                      for (int d = 0; d < 3; ++d) {
                          for (int mp1 = -1; mp1 <= 1; mp1 += 2) {
                              jxyz[d] = ixyz[d] + mp1;
                              if (jxyz[d] >= 0 && jxyz[d] < n[d]) {
                                  auto const neighbor_rank = owner(jxyz[Z],jxyz[Y],jxyz[X]);
                                  tension_taker += (neighbor_rank == taker);
                                  tension_giver += (neighbor_rank == giver);
                              } // in bounds
                          } // -/+1
                          jxyz[d] = ixyz[d]; // reset
                      } // d

                      float const tension = tension_giver - tension_taker;
                      if (tension > maxtension) {
                          maxtension = tension;
                          imax[X] = ix; imax[Y] = iy; imax[Z] = iz;
                      }
                      if (tension < mintension) {
                          mintension = tension;
                      }

                  } // giver
              }}} // ix iy iz
              if (echo > echo_it) std::printf("# iteration=%i smallest, largest tension %g %g at %i %i %i\n",
                                                 iteration, mintension, maxtension, imax[X], imax[Y], imax[Z]);
          } // scope

          // we can make a prediction for E_imbalance
          new_E_imbalance = 0; // energy contribution
          {
              set(new_load, nprocs, load); // copy
              auto const w8 = weight(imax[Z],imax[Y],imax[X]);
              new_load[giver] -= w8;
              new_load[taker] += w8;
              for (int32_t rank = 0; rank < nprocs; ++rank) {
                  new_E_imbalance += pow2(new_load[rank] - load_average);
              } // rank
          }
          if (echo > echo_it) std::printf("# iteration=%i old E= %g, new E= %g\n", iteration, E_imbalance, new_E_imbalance);

          // now we have to estimate the change in surface energy
          // new_n_surf = n_surf.copy();
          // determine the local coordination number at imax w.r.t. giver and taker
          // new_n_surf[giver] -= 1;
          // new_n_surf[taker] += 1;
          // maybe that is already done in maxtension....

          last_E_total = E_imbalance;

          if (new_E_imbalance < E_imbalance - 0*maxtension) {
              // accept the suggested change
              assert(owner(imax[Z],imax[Y],imax[X]) == giver);
              owner(imax[Z],imax[Y],imax[X]) = taker;
              if (echo > 11) std::printf("# owner[%i,%i,%i] changes from %i to %i\n", imax[X], imax[Y], imax[Z], giver, taker);
          } else {
              // stop the iterative procedure
              stop = 1 + iteration;
              if (echo > 3) std::printf("# stop in iteration %i\n", iteration);
          } // energy is minimized

      } // while

      if (echo > 3) std::printf("# after %d iterations initial E= %g, new E= %g\n", stop, E_imbalance_initial, new_E_imbalance);

      // compute new centroids
      for (int rank = 0; rank < nprocs; ++rank) {
          set(pos[rank], 4, 0.f);
      } // rank
      for (int iz = 0; iz < n[Z]; ++iz) {
      for (int iy = 0; iy < n[Y]; ++iy) {
      for (int ix = 0; ix < n[X]; ++ix) {
          auto const owner_rank = owner(iz,iy,ix);
          auto const w8        = weight(iz,iy,ix);
          int  const xyz1[] = {ix, iy, iz, 1};
          add_product(pos[owner_rank], 4, xyz1, w8);
      }}} // ix iy iz

      for (int rank = 0; rank < nprocs; ++rank) {
          auto *const position = pos[rank];
          auto const denom = position[3];
          if (denom > 0) {
              scale(position, 3, 1.f/denom);
          } else {
              set(position, 4, 0.f);
              assert(0 == denom);
          } // owns some
          if (echo > 9) std::printf("# myrank=%i   final coords %9.1f %9.1f %9.1f   load= %g\n", 
                                         rank, position[X], position[Y], position[Z], position[3]);
      } // rank

      analyze_load_imbalance(load, nprocs, echo);

      return stat;
  } // test_pressure_balancer

  // continuum theory:
  //    define pop(r,xyz) between 0 and 1
  //    normalized to fulfill sum_r pop(r,xyz) == 1 for each xyz
  //    total_load = sum_r sum_xyz pop(r,xyz) * weights(xyz)
  //    E_imbalance = sum_r (sum_xyz pop(r,xyz) * weights(xyz) - load_average)^2
  //    E_surface = sum_r sum_xyz (sum_neighborsxyz pop(r,xyz))*(6 - sum_neighborsxyz pop(r,xyz))
  //    apply diffusion to let domains leak into each other and explore the phase space,
  //    apply sharpening e.g. with 3x^2-2x^3, to minimize the surface energy and find a mapping


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

      if (nprocs < 2) return stat;

      auto const all = (n[X])*size_t(n[Y])*size_t(n[Z]);

      double w8sum_all{0};
      view2D<float> xyzw(all, 4, 0.f);
      std::vector<double> w8s(all, 0.0);

      for (int iz = 0; iz < n[Z]; ++iz) {
      for (int iy = 0; iy < n[Y]; ++iy) {
      for (int ix = 0; ix < n[X]; ++ix) {
          auto const iall = size_t(iz*n[Y] + iy)*n[X] + ix;
          assert(uint32_t(iall) == iall);
          xyzw(iall,X) = ix;
          xyzw(iall,Y) = iy;
          xyzw(iall,Z) = iz;
          w8s[iall]  = weight(ix,iy,iz);
          w8sum_all += w8s[iall];
      }}} // ix iy iz

      std::vector<double> load(nprocs, 0.0);

      if (echo > 0) std::printf("# %s: distribute %g blocks to %d processes\n\n", __func__, w8sum_all, nprocs);

for (int rank = 0; rank < nprocs; ++rank) {


      bool constexpr UNASSIGNED = 0, ASSIGNED = 1;
      std::vector<bool> state(all, UNASSIGNED);

      std::vector<uint32_t> indirect(all, 0);
      std::iota(indirect.begin(), indirect.end(), 0); // initialize with 0,1,....,all-1

      size_t nuna{all}; // number of unassigned blocks

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
          //    determine the center of weight
          size_t iuna{0};
          double cow[3] = {0, 0, 0};
          double w8sum{0};

//           for (size_t iall = 0; iall < all; ++iall) {
          auto const nunk = nuna; // copy old number of unassigned blocks
          for (size_t iunk = 0; iunk < nunk; ++iunk) {
              auto const iall = indirect[iunk]; // read from array indirect

              if (UNASSIGNED == state[iall]) {
                  auto const w8 = w8s[iall];
                  w8sum += w8;
                  // only use as weight if the weight is non-zero
                  add_product(cow, 3, xyzw[iall], double(w8 > 0));
                  indirect[iuna] = iall; // write to array indirect
                  ++iuna;
              } // is_unassigned
          } // iall
          nuna = iuna; // set new number of unassigned blocks
          assert(nuna <= nunk);

          // determine the largest distance^2 from the center
          double maxdist2{-1};
          int64_t imax{-1}; // index of the block with the largest distance
//           for (size_t iall = 0; iall < all; ++iall) {
//               if (UNASSIGNED == state[iall]) {
          for (size_t iuna = 0; iuna < nuna; ++iuna) {
              auto const iall = indirect[iuna];
              {
                  auto const *const xyz = xyzw[iall];
                  auto const dist2 = pow2(xyz[X] - cow[X])
                                   + pow2(xyz[Y] - cow[Y])
                                   + pow2(xyz[Z] - cow[Z]);
                  if (dist2 > maxdist2) {
                      maxdist2 = dist2;
                      imax = iall;
                  }
              } // is_unassigned
          } // i

          if (imax < 0) {
              // this happens if the process is idle, i.e. there are no blocks assigned to this one
              load_now = 0;
              np = 1; // stop the loop
          } else {

              add_product(cow, 3, xyzw[imax], -1.);

              auto const len2 = pow2(cow[X]) + pow2(cow[Y]) + pow2(cow[Z]);
              auto const norm = (len2 > 0) ? 1./std::sqrt(len2) : 0.0;
              double const vec[] = {cow[X]*norm, cow[Y]*norm, cow[Z]*norm};
              if (echo > 19) std::printf("# rank#%i sort along the [%g %g %g] direction\n", rank, vec[X], vec[Y], vec[Z]);

              using fui_t = std::pair<float,uint32_t>;
              std::vector<fui_t> v(nuna);

    //           iuna = 0;
    //           for (size_t iall = 0; iall < all; ++iall) {
    //               if (UNASSIGNED == state[i]) {
              for (size_t iuna = 0; iuna < nuna; ++iuna) {
                  auto const iall = indirect[iuna];
                  {
                      auto const *const xyz = xyzw[iall];
                      auto const f = xyz[X]*vec[X] + xyz[Y]*vec[Y] + xyz[Z]*vec[Z]; // inner product
                      v[iuna] = std::make_pair(float(f), uint32_t(iall));
    //                   ++iuna;
                  } // is_unassigned
              } // iall
              assert(iuna == nuna);

              auto lambda = [](fui_t i1, fui_t i2) { return i1.first < i2.first; };
              std::stable_sort(v.begin(), v.end(), lambda);

              auto const target_frac = nhalf[0]/double(np); // the "larger" half
              auto const target_load = target_frac*w8sum;
              {
                  auto const state0 = i01 ? ASSIGNED : UNASSIGNED,
                            state1 = i01 ? UNASSIGNED : ASSIGNED;
                  double load0{0}, load1{0};
                  int isrt;
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
                  assert(load0 + load1 == w8sum && "Maybe failed due to accuarcy issues");
                  load_now = i01 ? load1 : load0;
              }

              if (echo > 19) std::printf("# rank#%i assign %g of %g (%.2f %%, target %.2f %%) to %d processes\n", 
                                          rank, load_now, w8sum, load_now*100./w8sum, target_frac*100, nhalf[0]);

              // prepare for the next iteration
              rank_offset += i01*nhalf[0];
              np = nhalf[i01];

          } // imax < 0

      } // while

      if (echo > 9) std::printf("# rank#%i assign %.3f %%, target %.3f %%\n\n", rank, load_now*100/w8sum_all, 100./nprocs);
      load[rank] = load_now;

} // rank


      analyze_load_imbalance(load.data(), nprocs, echo);

      return stat;
  } // test_plane_balancer

  // example 5 processes -->
  //        rank#0      rank#1      rank#2      rank#3      rank#4
  //  take    3/5         3/5         3/5         2/5         2/5
  //  take    2/3         2/3         1/3         1/2         1/2
  //  take    1/2         1/2

#undef weight

  status_t all_tests(int const echo) { 
      status_t stat(0);

      int const nxyz[] = {int(control::get("load_balancer.test.nx", 18.)), // number of blocks
                          int(control::get("load_balancer.test.ny", 19.)),
                          int(control::get("load_balancer.test.nz", 20.))};
      int const nprocs = control::get("load_balancer.test.nprocs", 53.);
      if (echo > 0) std::printf("\n\n# %s start %d x %d x %d = %d with %d MPI processes\n", __func__, nxyz[X], nxyz[Y], nxyz[Z], nxyz[X]*nxyz[Y]*nxyz[Z], nprocs);

//    stat += test_bisection_balancer(nprocs, nxyz, echo);
//    stat += test_diffusion_balancer(nprocs, nxyz, echo);
//    stat += test_RobinHood_balancer(nprocs, nxyz, echo);
//    stat += test_pressure_balancer(nprocs, nxyz, echo);
      stat += test_plane_balancer(nprocs, nxyz, echo);
      return stat;
  } // all_tests

} // namespace load_balancer

#endif // NO_UNIT_TESTS
