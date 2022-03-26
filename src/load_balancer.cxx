
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
#include "inline_math.hxx" // pow2, scale
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
      
      view2D<float> pos(nprocs, 4, 0.f); // centroid positions
      view3D<int32_t> owner(n[Z], n[Y], n[X], -1); // owning ranks for each block

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

          std::vector<float> n_own(nprocs, 0.f);
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
              n_own[l_rank] += w8;
          }}} // ix iy iz

          simple_stats::Stats<> st(0);
          for (int myrank = 0; myrank < nprocs; ++myrank) {
              st.add(n_own[myrank]);
              if (echo > 9) std::printf("# myrank=%i owns %g\n", myrank, n_own[myrank]);
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

//       view2D<float> frc(nprocs, 4, 0.f); // forces
      view3D<uint8_t> n_surface(n[Z], n[Y], n[X], 0.f);

      float const n_own_average = blocksperproc;
      double last_E_total;
      float last_maxfrc;

      for (int iteration = 0; iteration < 99; ++iteration) {

          double E_surface{0}, E_imbalance{0}; // energy
          // E_imbalance ~ sum_r (n_own[r] - <n_own>)^2 = sum_r n_own[r]^2 - constant // square deviation
          // <n_own> = ngrid/nprocs is constant
          // n_own[r] = sum_xyz (owner[xyz] == r)
          //
          // E_imbalance = sum_r ( sum_xyz (owner[xyz] == r) )^2
          // E_surface   = sum_xyz (owner[xyz] != owner[xyz +/- 1])^2
          
          std::vector<float> n_own(nprocs, 0.f);
          for (int iz = 0; iz < n[Z]; ++iz) {
          for (int iy = 0; iy < n[Y]; ++iy) {
          for (int ix = 0; ix < n[X]; ++ix) {
              auto const owner_rank = owner(iz,iy,ix);
              float const w8 = 1; // weight(iz,iy,ix);
              n_own[owner_rank] += w8;
              uint8_t number_of_different_neighbors{0};
#define STAR
#ifdef  STAR
              int const ixyz[] = {ix, iy, iz};
              int jxyz[] = {ix, iy, iz}; // non-const
              // loop over Cartesian neighborhood star (6 points)
              for (int d = 0; d < 3; ++d) {
                  for (int mp1 = -1; mp1 <= 1; mp1 += 2) {
                      jxyz[d] = ixyz[d] + mp1;
                      if (jxyz[d] >= 0 && jxyz[d] < n[d]) {
                          number_of_different_neighbors += (owner_rank != owner(jxyz[Z],jxyz[Y],jxyz[X]));
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
              n_surface(iz,iy,ix) = number_of_different_neighbors;
          }}} // ix iy iz
          for (int rank = 0; rank < nprocs; ++rank) {
              E_imbalance += pow2(n_own[rank]);
          } // rank
          auto const E_total = E_imbalance + E_surface;
          if (echo > 0) std::printf("# E_total= %g = %g + %g (E_imbalance + E_surface)\n",
                                      E_total,               E_imbalance,  E_surface);

          // now compute the gradients to suggest a change
          // d E_imbalance
          // ------------- = 2 n_own[r]
          //  d n_own[r]
          
          //  d n_own[r]
          // ------------ = 1 * (owner[xyz] != r) - 1 * (owner[xyz] == r)
          // d owner[xyz]

          float maxfrc{-9e9};
          int imax[3] = {-1, -1, -1};
          for (int iz = 0; iz < n[Z]; ++iz) {
          for (int iy = 0; iy < n[Y]; ++iy) {
          for (int ix = 0; ix < n[X]; ++ix) {
              auto const owner_rank = owner(iz,iy,ix);
              float const frc_imbalance = 2*(n_own[owner_rank] - n_own_average);
              auto const frc_total = frc_imbalance + n_surface(iz,iy,ix);
              if (frc_total > maxfrc) {
                  maxfrc = frc_total;
                  imax[X] = ix;  imax[Y] = iy;  imax[Z] = iz;
              } // is larger
          }}} // ix iy iz

          float const maxfrc_surf = n_surface(imax[Z],imax[Y],imax[X]);
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
                          auto const owner_neigh = owner(jxyz[Z],jxyz[Y],jxyz[X]);
                          if (owner_rank != owner_neigh) {
                              neigh[i6] = owner_neigh;
                              ++i6;
                          }
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


  status_t all_tests(int const echo) { 
      status_t stat(0);
      stat += test_diffusion_balancer(echo);
//    stat += test_bisection_balancer(echo);
      return stat;
  } // all_tests

} // namespace load_balancer

#endif // NO_UNIT_TESTS
