#pragma once

#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cmath> // std::sqrt
#include <algorithm> // std::max
#include <utility> // std::swap //, ::move
#include <vector> // std::vector<T>
#include <cstdio> // std::printf, ::snprintf

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "simple_timer.hxx" // SimpleTimer
#include "unit_system.hxx" // eV, _eV, Ang, _Ang
#include "print_tools.hxx" // printf_vector
#include "global_coordinates.hxx" // ::get
#include "green_input.hxx" // ::load_Hamitonian
#include "green_memory.hxx" // get_memory, free_memory, real_t_name
#include "green_sparse.hxx" // ::sparse_t<,>
#include "green_kinetic.hxx" // ::finite_difference_plan_t, index3D
#include "green_action.hxx" // ::plan_t, ::action_t
#include "sho_tools.hxx" // ::nSHO
#include "control.hxx" // ::get

/*
 *  ToDo plan:
 *    Implement CPU version of SHOprj and SHOadd
 *    Implement periodic boundary conditions in CPU version [optional, but would speed up the following item]
 *    Ensure that tfQMRgpu can invert CPU-only
 *    Verify that the density of states of the Green function matches that found by eigenvalues (explicit solver)
 *    Implement GPU versions of the Hamiltonian
 * 
 */

 /*
  *  Future plan:
  *   Support also density matrix purification scheme (McWeeney filter: x^2(3-2x)
  *   or with a non-trivial overlap operator S (3xSx - 2xSxSx from Phys. Rev. B 50, 17622)
  *   maybe better with norm-conserving PAW formalism --> S == 1
  */

namespace green_function {

  double const GByte = 1e-9; char const *const _GByte = "GByte";

  template <typename real_t, int R1C2=2, int Noco=1>
  void try_action(green_action::plan_t const & p, int const n_iterations=1, int const echo=9) {
      if (n_iterations >= 0) { // scope: try to apply the operator
          if (echo > 1) { std::printf("# %s<%s,R1C2=%d,Noco=%d>\n", __func__, real_t_name<real_t>(), R1C2, Noco); std::fflush(stdout); }
          green_action::action_t<real_t,R1C2,Noco,64> action(&p); // constructor

          uint32_t const nnzbX = p.colindx.size();
          int constexpr LM = Noco*64;
          auto x = get_memory<real_t[R1C2][LM][LM]>(nnzbX, echo, "x");
          auto y = get_memory<real_t[R1C2][LM][LM]>(nnzbX, echo, "y");
          for (size_t i = 0; i < nnzbX*R1C2*LM*LM; ++i) {
              x[0][0][0][i] = 0; // init x
          } // i
          
          bool const toy = control::get("green_function.benchmark.toy", 1.);

          // benchmark the action
          for (int iteration = 0; iteration < n_iterations; ++iteration) {
              SimpleTimer timer(__FILE__, __LINE__);
              if (echo > 5) { std::printf("# iteration #%i\n", iteration); std::fflush(stdout); }
              if (toy) {
                  action.toy_multiply(y, x, p.colindx.data(), nnzbX, p.nCols);
              } else {
                  action.multiply(y, x, p.colindx.data(), nnzbX, p.nCols);
              }
              std::swap(x, y);
          } // iteration

          free_memory(y);
          free_memory(x);
      } // n_iterations >= 0
  } // try_action


  inline status_t construct_Green_function(
        green_action::plan_t & p // result, create a plan how to apply the SHO-PAW Hamiltonian to a block-sparse truncated Green function
      , int const ng[3] // number of grid points of the unit cell in with the potential is defined
      , double const hg[3] // grid spacings
      , std::vector<double> const & Veff // [ng[2]*ng[1]*ng[0]]
      , std::vector<double> const & xyzZinso // [natoms*8]
      , std::vector<std::vector<double>> const & atom_mat // atomic hamiltonian and overlap matrix
      , int const echo=0 // log-level
      , std::complex<double> const *energy_parameter=nullptr // E in G = (H - E*S)^{-1}
      , int const Noco=1
  ) {
      int constexpr X=0, Y=1, Z=2;
      if (echo > 0) std::printf("\n#\n# %s(%i, %i, %i)\n#\n\n", __func__, ng[X], ng[Y], ng[Z]);
      
      assert(1 == Noco || 2 == Noco);

      std::complex<double> E_param(energy_parameter ? *energy_parameter : 0);

      int32_t const MANIPULATE = control::get("MANIPULATE", 0.);
      int32_t n_original_Veff_blocks[3] = {0, 0, 0};
      for (int d = 0; d < 3; ++d) {
          n_original_Veff_blocks[d] = (ng[d] >> 2); // divided by 4
          assert(ng[d] == 4*n_original_Veff_blocks[d] && "All grid dimensions must be a multiple of 4!");
          if (MANIPULATE) n_original_Veff_blocks[d] = std::min(n_original_Veff_blocks[d], MANIPULATE);
      } // d
      if (MANIPULATE && echo > 0) std::printf("\n# MANIPULATE=%d for n_original_Veff_blocks:\n", MANIPULATE);
      if (echo > 3) { std::printf("# n_original_Veff_blocks "); printf_vector(" %d", n_original_Veff_blocks, 3); }

      assert(Veff.size() == ng[Z]*ng[Y]*ng[X]);

//    view3D<block4<double>> Vtot_blocks(n_original_Veff_blocks[Z], n_original_Veff_blocks[Y], n_original_Veff_blocks[X]);
      size_t const n_all_Veff_blocks = n_original_Veff_blocks[Z]*n_original_Veff_blocks[Y]*n_original_Veff_blocks[X];

      // regroup effective potential into blocks of 4x4x4
//    auto Vtot_gpu = get_memory<double>(n_all_Veff_blocks*64, echo, "Vtot_gpu"); // does not really work flawlessly
//    view4D<double> Vtot(Vtot_gpu, n_original_Veff_blocks[Y], n_original_Veff_blocks[X], 64); // wrap
//    view4D<double> Vtot(n_original_Veff_blocks[Z], n_original_Veff_blocks[Y], n_original_Veff_blocks[X], 64); // get memory
      p.Veff = get_memory<double[64]>(n_all_Veff_blocks*Noco*Noco, echo, "Veff"); // in managed memory
      { // scope: reorder Veff into block-structured p.Veff

          for (size_t k = 0; k < n_all_Veff_blocks*Noco*Noco*64; ++k) {
              p.Veff[0][k] = 0; // clear
          } // k

          for (int ibz = 0; ibz < n_original_Veff_blocks[Z]; ++ibz) {
          for (int iby = 0; iby < n_original_Veff_blocks[Y]; ++iby) { // block index loops, parallel
          for (int ibx = 0; ibx < n_original_Veff_blocks[X]; ++ibx) {
              int const ibxyz[] = {ibx, iby, ibz};
              auto const Veff_index = index3D(n_original_Veff_blocks, ibxyz);
              for (int i4z = 0; i4z < 4; ++i4z) {
              for (int i4y = 0; i4y < 4; ++i4y) { // grid point index in [0, 4) loops
              for (int i4x = 0; i4x < 4; ++i4x) {
                  int const i64 = (i4z*4 + i4y)*4 + i4x;
                  size_t const izyx = (size_t(ibz*4 + i4z)
                              *ng[Y] + size_t(iby*4 + i4y)) 
                              *ng[X] + size_t(ibx*4 + i4x); // global grid point index
                  assert(izyx < size_t(ng[Z])*size_t(ng[Y])*size_t(ng[X]));
                  p.Veff[Veff_index*Noco*Noco + 0][i64] = Veff[izyx]; // copy potential value
              }}} // i4
          }}} // xyz
      } // scope

      // Cartesian cell parameters for the unit cell in which the potential is defined
      double const cell[3] = {ng[X]*hg[X], ng[Y]*hg[Y], ng[Z]*hg[Z]};

      // assume periodic boundary conditions and an infinite host crystal,
      // so there is no need to consider k-points

      // determine the largest and smallest indices of target blocks
      // given a max distance r_trunc between source blocks and target blocks

      // assume that the source blocks lie compact in space
      int32_t const nRHSs = n_original_Veff_blocks[Z] * n_original_Veff_blocks[Y] * n_original_Veff_blocks[X];
      if (echo > 0) std::printf("# total number of source blocks is %d\n", nRHSs);
      assert(nRHSs >= 0);
      view2D<uint32_t> global_source_coords(nRHSs, 4, 0); // unsigned since they must be in [0, 2^21) anyway
      p.global_source_indices.resize(nRHSs); // [nRHSs]
      p.nCols = nRHSs;
      double center_of_mass_RHS[] = {0, 0, 0};
      double center_of_RHSs[]     = {0, 0, 0};
      int32_t min_source_coords[] = {0, 0, 0}; // global coordinates
      int32_t max_source_coords[] = {0, 0, 0}; // global coordinates
      int32_t internal_global_offset[] = {0, 0, 0};
      double max_distance_from_comass{0};
      double max_distance_from_center{0};
      { // scope: determine min, max, center
          {   int iRHS{0};
              for (int ibz = 0; ibz < n_original_Veff_blocks[Z]; ++ibz) {
              for (int iby = 0; iby < n_original_Veff_blocks[Y]; ++iby) { // block index loops, serial
              for (int ibx = 0; ibx < n_original_Veff_blocks[X]; ++ibx) {
                  int const ibxyz[] = {ibx, iby, ibz, 0};
                  set(global_source_coords[iRHS], 4, ibxyz);
                  p.global_source_indices[iRHS] = global_coordinates::get(ibx, iby, ibz);
                  for (int d = 0; d < 3; ++d) {
                      int32_t const rhs_coord = global_source_coords(iRHS,d);
                      center_of_mass_RHS[d] += (rhs_coord*4 + 1.5)*hg[d];
                      min_source_coords[d] = std::min(min_source_coords[d], rhs_coord);
                      max_source_coords[d] = std::max(max_source_coords[d], rhs_coord);
                  } // d
                  ++iRHS;
              }}} // xyz
              assert(nRHSs == iRHS);
          } // iRHS

          for (int d = 0; d < 3; ++d) {
              center_of_mass_RHS[d] /= std::max(1, nRHSs);
              auto const middle2 = min_source_coords[d] + max_source_coords[d];
              internal_global_offset[d] = middle2/2;
              center_of_RHSs[d] = ((middle2*0.5)*4 + 1.5)*hg[d];
          } // d

          if (echo > 0) std::printf("# internal and global coordinates differ by %d %d %d\n",
              internal_global_offset[X], internal_global_offset[Y], internal_global_offset[Z]);

          p.source_coords = get_memory<int16_t[4]>(nRHSs, echo, "source_coords"); // internal coordinates
          { // scope: compute also the largest distance from the center or center of mass
              double max_d2m{0}, max_d2c{0};
              for (int iRHS = 0; iRHS < nRHSs; ++iRHS) {
                  double d2m{0}, d2c{0};
                  for (int d = 0; d < 3; ++d) {
                      p.source_coords[iRHS][d] = global_source_coords(iRHS,d) - internal_global_offset[d];
                      d2m += pow2((global_source_coords(iRHS,d)*4 + 1.5)*hg[d] - center_of_mass_RHS[d]);
                      d2c += pow2((global_source_coords(iRHS,d)*4 + 1.5)*hg[d] - center_of_RHSs[d]);
                  } // d
                  p.source_coords[iRHS][3] = 0; // not used
                  max_d2c = std::max(max_d2c, d2c);
                  max_d2m = std::max(max_d2m, d2m);
              } // iRHS
              max_distance_from_center = std::sqrt(max_d2c);
              max_distance_from_comass = std::sqrt(max_d2m);
          } // scope

      } // scope
      if (echo > 0) std::printf("# center of mass of RHS blocks is %g %g %g %s\n",
          center_of_mass_RHS[X]*Ang, center_of_mass_RHS[Y]*Ang, center_of_mass_RHS[Z]*Ang, _Ang);
      if (echo > 0) std::printf("# center of coords  RHS blocks is %g %g %g %s\n",
          center_of_RHSs[X]*Ang, center_of_RHSs[Y]*Ang, center_of_RHSs[Z]*Ang, _Ang);
      if (echo > 0) std::printf("# largest distance of RHS blocks from center of mass is %g, from center is %g %s\n",
                              max_distance_from_comass*Ang, max_distance_from_center*Ang, _Ang);


      // truncation radius
      double const r_trunc = control::get("green_function.truncation.radius", 10.);
      double const r_proj  = control::get("green_function.projection.radius", 6.); // in units of sigma
      if (echo > 0) std::printf("# green_function.truncation.radius=%g %s\n", r_trunc*Ang, _Ang);
      p.r_truncation = std::max(0., r_trunc);
      p.r_Vconfinement = std::min(std::max(0., r_trunc - 2.0), p.r_truncation);

// example:
//    truncation radius in Cu (fcc)
//    lattice constant = 3.522 Ang
//    volume per atom = alat^3 / 4 == 10.922 Ang^3 / atom
//    1000 atoms inside the truncation cluster
//    truncation sphere volume 10922 Ang^3
//    truncation radius = cbrt(3*V/(4*pi)) = 13.764 Ang = 26 Bohr
//    8000 atoms --> 52 Bohr
//    64000 atoms --> 104 Bohr

      double const r_block_circumscribing_sphere = 0.5*(4 - 1)*std::sqrt(pow2(hg[X]) + pow2(hg[Y]) + pow2(hg[Z]));
      if (echo > 0) std::printf("# circumscribing radius= %g %s\n", r_block_circumscribing_sphere*Ang, _Ang);
      

      // count the number of green function elements for each target block

      uint16_t num_target_coords[3] = {0, 0, 0};
      int16_t  min_target_coords[3] = {0, 0, 0}; // internal coordinates
      int16_t  max_target_coords[3] = {0, 0, 0}; // internal coordinates
      { // scope: create the truncated Green function block-sparsity pattern
          auto const rtrunc       = std::max(0., r_trunc);
          auto const rtrunc_plus  =              rtrunc + 2*r_block_circumscribing_sphere;
          auto const rtrunc_minus = std::max(0., rtrunc - 2*r_block_circumscribing_sphere);
          if (echo > 0) std::printf("# truncation radius %g, search within %g %s\n", rtrunc*Ang, rtrunc_plus*Ang, _Ang);
          if (echo > 0 && rtrunc_minus > 0) std::printf("# blocks with center distance below %g %s are fully inside\n", rtrunc_minus*Ang, _Ang);

          int16_t itr[3]; // 16bit, range [-32768, 32767] should be enough
          for (int d = 0; d < 3; ++d) {
              // how many blocks around the source block do we need to check
              auto const itrunc = std::floor(rtrunc_plus/(4*hg[d]));
              assert(itrunc < 32768 && "target coordinate type is int16_t!");
              itr[d] = int16_t(itrunc);
              assert(itr[d] >= 0);
              min_target_coords[d] = min_source_coords[d] - internal_global_offset[d] - itr[d];
              max_target_coords[d] = max_source_coords[d] - internal_global_offset[d] + itr[d];
              num_target_coords[d] = max_target_coords[d] + 1 - min_target_coords[d];
          } // d
          auto const product_target_blocks = size_t(num_target_coords[Z])*
                                             size_t(num_target_coords[Y])*
                                             size_t(num_target_coords[X]);
          if (echo > 0) std::printf("# all targets within (%i, %i, %i) and (%i, %i, %i) --> %d x %d x %d = %.3f k\n",
              min_target_coords[X], min_target_coords[Y], min_target_coords[Z],
              max_target_coords[X], max_target_coords[Y], max_target_coords[Z], 
              num_target_coords[X], num_target_coords[Y], num_target_coords[Z], product_target_blocks*.001);
          assert(product_target_blocks > 0);
          std::vector<std::vector<uint16_t>> column_indices(product_target_blocks);

          double const r2trunc        = pow2(rtrunc),
                       r2trunc_plus   = pow2(rtrunc_plus),
                       r2trunc_minus  = pow2(rtrunc_minus),
                       r2block_circum = pow2(r_block_circumscribing_sphere*3); // Maybe this can be reduced to *2

          std::vector<int32_t> tag_diagonal(product_target_blocks, -1);
          assert(nRHSs < 65536 && "the integer type of ColIndex is uint16_t!");
          std::vector<std::vector<bool>> sparsity_pattern(nRHSs); // std::vector<bool> is a memory-saving bit-array

          for (uint16_t iRHS = 0; iRHS < nRHSs; ++iRHS) {
              auto & sparsity_RHS = sparsity_pattern[iRHS]; // abbreviate
              sparsity_RHS.resize(product_target_blocks, false);
              auto const *const source_coords = p.source_coords[iRHS]; // internal source block coordinates
              simple_stats::Stats<> stats[3];
              int constexpr max_nci = 27;
              std::vector<int> hist(1 + max_nci, 0); // distribution of nci
              std::vector<simple_stats::Stats<>> stats_d2(1 + max_nci);
              for (int16_t bz = -itr[Z]; bz <= itr[Z]; ++bz) {
              for (int16_t by = -itr[Y]; by <= itr[Y]; ++by) {
              for (int16_t bx = -itr[X]; bx <= itr[X]; ++bx) {
                  int16_t const bxyz[3] = {bx, by, bz}; // block difference vector
                  int16_t target_coords[3]; // internal target block coordinates
                  for (int d = 0; d < 3; ++d) {
                      target_coords[d] = source_coords[d] + bxyz[d];
                      assert(target_coords[d] >= min_target_coords[d]);
                      assert(target_coords[d] <= max_target_coords[d]);
                  } // d
                  // d2 is the distance^2 of the block centers
                  double const d2 = pow2(bx*4*hg[X]) + pow2(by*4*hg[Y]) + pow2(bz*4*hg[Z]);

                  int nci{0}; // init number of corners inside
                  if (d2 < r2trunc_plus) { // potentially inside, check all 8 or 27 corner cases
//                    if (d2 < r2trunc_minus) { nci = max_nci; } else // skip the 8- or 27-corners test for inner blocks -> some speedup
                      { // scope: 8 or 27 corner test
                          int const far = (d2 > r2block_circum);
                          int const mci = far ? 8 : 27;
                          // i = i4 - j4 --> i in [-3, 3],
                          //     if two blocks are far from each other, we test only the 8 combinations of |{-3, 3}|^3
                          //     for blocks close to each other, we test all 27 combinations of |{-3, 0, 3}|^3
                          for (int iz = -3; iz <= 3; iz += 3 + 3*far) { double const d2z = pow2((bz*4 + iz)*hg[Z]);
                          for (int iy = -3; iy <= 3; iy += 3 + 3*far) { double const d2y = pow2((by*4 + iy)*hg[Y]) + d2z;
                          for (int ix = -3; ix <= 3; ix += 3 + 3*far) { double const d2c = pow2((bx*4 + ix)*hg[X]) + d2y;
#if 0
                              if (0 == iRHS && (d2c < r2trunc) && echo > 17) {
                                  std::printf("# %s: b= %i %i %i, i-j %i %i %i, d^2= %g %s\n", 
                                      __func__, bx,by,bz, ix,iy,iz, d2c, (d2c < r2trunc)?"in":"out");
                              }
#endif // 0
                              nci += (d2c < r2trunc); // add 1 if inside
                          }}} // ix iy iz
                          if (d2 < r2trunc_minus) assert(mci == nci); // for these, we could skip the 8-corners test
                          nci = (nci*27)/mci; // limit nci to [0, 27]
                      } // scope

                  } // d2 < r2trunc_plus

                  if (nci > 0) {
                      // any grid point in block (bx,by,bz) is closer than rtrunc to any grid point in the source block
                      int16_t idx[3];
                      for (int d = 0; d < 3; ++d) {
                          idx[d] = target_coords[d] - min_target_coords[d];
                          assert(idx[d] >= 0);
                          assert(idx[d] < num_target_coords[d]);
                          stats[d].add(target_coords[d]);
                      } // d
                      auto const idx3 = index3D(num_target_coords, idx);
                      assert(idx3 < product_target_blocks);
                      column_indices[idx3].push_back(iRHS);
                      sparsity_RHS[idx3] = true;
                      if (0 == bx && 0 == by && 0 == bz) {
                          tag_diagonal[idx3] = iRHS;
                      } // diagonal entry
                  } // inside
                  ++hist[nci];
                  stats_d2[nci].add(d2);

              }}} // xyz
              if (echo > 7) {
                  std::printf("# RHS at %i %i %i reaches from (%g, %g, %g) to (%g, %g, %g)\n",
                      global_source_coords(iRHS,X), global_source_coords(iRHS,Y), global_source_coords(iRHS,Z),
                      stats[X].min(), stats[Y].min(), stats[Z].min(),
                      stats[X].max(), stats[Y].max(), stats[Z].max());
              } // echo
              if (echo > 2) {
                  if (0 == iRHS) {
                      int total_checked{0};
                      // list in detail
                      for (int nci = 0; nci <= max_nci; ++nci) {
                          if (hist[nci] > 0) {
                              std::printf("# RHS has%9.3f k cases with %d corners inside, d2 stats: %g +/- %g in [%g, %g] Bohr^2\n",
                                  hist[nci]*.001, nci, stats_d2[nci].mean(), stats_d2[nci].dev(), stats_d2[nci].min(), stats_d2[nci].max());
                              total_checked += hist[nci];
                          } // hist[nci] > 0
                      } // nci
                      auto const partial = total_checked - hist[0] - hist[max_nci];
                      std::printf("# RHS has %.3f k inside, %.3f k partial and %.3f k outside (of %.3f k checked blocks)\n",
                                    hist[max_nci]*.001, partial*.001, hist[0]*.001, total_checked*.001);
                  } // RHS #0
              } // echo
          } // iRHS

          // create a histogram about the distribution of the number of columns per row
          std::vector<uint32_t> hist(1 + nRHSs, 0);
          for (size_t idx3 = 0; idx3 < column_indices.size(); ++idx3) {
              auto const n = column_indices[idx3].size();
              ++hist[n];
          } // idx3

          // eval the histogram
          size_t nall{0};
          size_t nnz{0}; // number of non-zero BSR entries
          for (int n = 0; n <= nRHSs; ++n) {
              nall += hist[n];
              nnz  += hist[n]*n;
          } // n
          if (echo > 5) { std::printf("# histogram total= %.3f k: ", nall*.001); printf_vector(" %d", hist.data(), nRHSs + 1); }
          assert(nall == product_target_blocks); // sanity check

          p.nRows = nall - hist[0]; // the target block entries with no RHS do not create a row
          if (echo > 0) std::printf("# total number of Green function blocks is %.3f k, "
                               "average %.1f per source block\n", nnz*.001, nnz/double(nRHSs));
          if (echo > 0) std::printf("# %.3f k (%.1f %% of %.3f k) target blocks are active\n", 
              p.nRows*.001, p.nRows/(product_target_blocks*.01), product_target_blocks*.001);

          assert(nnz < (uint64_t(1) << 32) && "the integer type or RowStart is uint32_t!");

          // resize BSR tables: (Block-compressed Sparse Row format)
          if (echo > 3) { std::printf("# memory of Green function is %.6f %s (float, twice for double)\n",
                              nnz*2.*64.*64.*sizeof(float)*GByte, _GByte); std::fflush(stdout); }
          p.colindx.resize(nnz);
          p.rowindx = get_memory<int32_t>(nnz, echo, "rowindx");
          p.RowStart = get_memory<uint32_t>(p.nRows + 1, echo, "RowStart");
          p.RowStart[0] = 0;
          p.veff_index = get_memory<int32_t>(p.nRows, echo, "veff_index"); // indirection list for the local potential
          set(p.veff_index, p.nRows, -1); // init as non-existing
          p.target_coords = get_memory<int16_t[3+1]>(p.nRows, echo, "target_coords");
          p.target_minus_source = get_memory<int16_t[3+1]>(nnz, echo, "target_minus_source");
          p.CubePos = get_memory<float[3+1]>(p.nRows, echo, "CubePos");

          p.global_target_indices.resize(p.nRows);
          p.subset.resize(p.nCols); // we assume columns of the unit operator as RHS

          view3D<int32_t> iRow_of_coords(num_target_coords[Z],
                                         num_target_coords[Y],
                                         num_target_coords[X], -1); // init as non-existing

          { // scope: fill BSR tables
              simple_stats::Stats<> st;
              uint32_t iRow{0};
              for (uint16_t z = 0; z < num_target_coords[Z]; ++z) { // serial
              for (uint16_t y = 0; y < num_target_coords[Y]; ++y) { // serial
              for (uint16_t x = 0; x < num_target_coords[X]; ++x) { // serial
                  int const idx[3] = {x, y, z};
                  auto const idx3 = index3D(num_target_coords, idx);
                  assert(idx3 < product_target_blocks);

                  auto const n = column_indices[idx3].size();
                  if (n > 0) {
                      st.add(n);
                      iRow_of_coords(idx[Z], idx[Y], idx[X]) = iRow; // set existing

                      p.RowStart[iRow + 1] = p.RowStart[iRow] + n;
                      // copy the column indices
                      set(p.colindx.data() + p.RowStart[iRow], n, column_indices[idx3].data());
                      // copy the target block coordinates
                      int32_t global_target_coords[3];
                      for (int d = 0; d < 3; ++d) {
                          p.target_coords[iRow][d] = idx[d] + min_target_coords[d];
                          p.CubePos[iRow][d] = p.target_coords[iRow][d];
                          global_target_coords[d] = p.target_coords[iRow][d] + internal_global_offset[d];
                      } // d
                      p.target_coords[iRow][3] = 0; // not used
                      p.CubePos[iRow][3] = 0.f; // not used

                      p.global_target_indices[iRow] = global_coordinates::get(global_target_coords); 
                      // global_target_indices is needed to gather the local potential data from other MPI processes

                      { // scope: determine the diagonal entry (source == target)
                          auto const iCol = tag_diagonal[idx3];
                          if (iCol > -1) {
                              assert(iCol < (1ul << 16)); // number range of uint16_t
                              for (int d = 0; d < 3; ++d) {
                                  // sanity check on internal coordinates
                                  assert(p.source_coords[iCol][d] == p.target_coords[iRow][d]);
                              } // d
                              { // search inz such that p.colindx[inz] == iCol
                                  int32_t inz_found{-1};
                                  for (int32_t inz = p.RowStart[iRow]; inz < p.RowStart[iRow + 1] && inz_found == -1; ++inz) {
                                      if (iCol == p.colindx[inz]) {
                                          inz_found = inz;
                                      }
                                  } // inz
                                  assert(-1 != inz_found);
                                  p.subset[iCol] = inz_found;
                              } // search
                          } // iCol valid
                      } // scope

                      for (int32_t inz = p.RowStart[iRow]; inz < p.RowStart[iRow + 1]; ++inz) {
                          auto const iCol = p.colindx[inz];
                          for (int d = 0; d < 3; ++d) {
                              p.target_minus_source[inz][d] = p.target_coords[iRow][d] - p.source_coords[iCol][d]; // ToDo: with periodic boundary conditions, we need to find the shortest distance vector accounting for periodic images of the cell
                          } // d
                          p.target_minus_source[inz][3] = 0; // not used
                          p.rowindx[inz] = iRow;
                      } // inz

                      if (1) { // scope: fill indirection table for having the local potential only defined in 1 unit cell and repeated periodically
                          int32_t mod[3];
                          for (int d = 0; d < 3; ++d) {
                              mod[d] = global_target_coords[d] % n_original_Veff_blocks[d];
                              mod[d] += (mod[d] < 0)*n_original_Veff_blocks[d];
                          } // d
                          p.veff_index[iRow] = index3D(n_original_Veff_blocks, mod);
                      } // scope
 
                      // count up the number of active rows
                      ++iRow;
                  } // n > 0
              }}} // idx
              assert(p.nRows == iRow);
              assert(nnz == p.RowStart[p.nRows]);
              std::printf("# source blocks per target block: average %.1f +/- %.1f in [%g, %g]\n", st.mean(), st.dev(), st.min(), st.max());
          } // scope
          column_indices.clear(); // not needed beyond this point

          if (echo > 1) { // measure the difference in the number of target blocks of each RHS
              std::vector<uint32_t> nt(nRHSs, 0);
              // traverse the BSR structure
              for (uint32_t iRow = 0; iRow < p.nRows; ++iRow) {
                  for (auto inz = p.RowStart[iRow]; inz < p.RowStart[iRow + 1]; ++inz) {
                      auto const iCol = p.colindx[inz];
                      ++nt[iCol];
                  } // inz
              } // iRow
              // analyze nt
              simple_stats::Stats<> st;
              for (uint16_t iRHS = 0; iRHS < nRHSs; ++iRHS) {
                  st.add(nt[iRHS]);
              } // iRHS
              std::printf("# target blocks per source block: average %.1f +/- %.1f in [%g, %g]\n", st.mean(), st.dev(), st.min(), st.max());
          } // echo

          // Green function is stored sparse 
          // as std::complex<real_t> green[nnz][64][64] 
          // or real_t green[nnz][2][64][64] for the GPU;

          for (int dd = 0; dd < 3; ++dd) { // derivate direction
              // create lists for the finite-difference derivatives
              p.fd_plan[dd] = green_kinetic::finite_difference_plan_t(dd
                , num_target_coords
                , p.RowStart, p.colindx.data()
                , iRow_of_coords
                , sparsity_pattern.data()
                , nRHSs, echo);
          } // dd

          // transfer grid spacing into managed GPU memory
          p.grid_spacing = get_memory<double>(4, echo, "grid_spacing");
          set(p.grid_spacing, 3, hg);
          p.grid_spacing[3] = r_proj; // radius in units of sigma at which the projection stops

      } // scope







      int const natoms = atom_mat.size();
      if (echo > 2) std::printf("\n#\n# %s: Start atom part, %d atoms\n#\n", __func__, natoms);

      // compute which atoms will contribute, the list of natoms atoms may contain a subset of all atoms
      double max_projection_radius{0};
      for (int ia = 0; ia < natoms; ++ia) {
          auto const sigma = xyzZinso[ia*8 + 6];
          auto const projection_radius = std::max(0.0, r_proj*sigma);
          max_projection_radius = std::max(max_projection_radius, projection_radius);
      } // ia
      if (echo > 3) std::printf("# largest projection radius is %g %s\n", max_projection_radius*Ang, _Ang);


      p.ApcStart = nullptr;
      p.natom_images = 0;
      { // scope:
          SimpleTimer atom_list_timer(__FILE__, __LINE__, "Atom part");
          
          auto const radius = r_trunc + max_distance_from_center + 2*max_projection_radius + 2*r_block_circumscribing_sphere;
          int iimage[3];
          size_t nimages{1};
          for (int d = 0; d < 3; ++d) { // parallel
              iimage[d] = std::ceil(radius/cell[d]);  // for periodic boundary conditions
              iimage[d] = 0;                          // for isolated boundary conditions
              nimages *= (iimage[d] + 1 + iimage[d]);
          } // d
          warn("used isolated boundary conditions for the generation of %ld atom images", nimages);
          auto const natom_images = natoms*nimages;
          if (echo > 3) std::printf("# replicate %d %d %d atom images, %.3f k images total\n",
                                        iimage[X], iimage[Y], iimage[Z], nimages*.001);

          std::vector<uint32_t> ApcStart(natom_images + 1, 0); // probably larger than needed, should call resize(nai + 1) later
          std::vector<green_action::atom_t> atom_data(natom_images);
          std::vector<uint16_t> atom_ncoeff(natoms, 0); // 0: atom does not contribute

          simple_stats::Stats<> nc_stats;
          double sparse{0}, dense{0}; // stats to assess how much memory can be saved using sparse storage

          std::vector<std::vector<uint32_t>> cubes; // stores the row indices of Green function rows
          cubes.reserve(natom_images); // maximum (needs 24 Byte per atom image)
          
          std::vector<std::vector<uint32_t>> imags; // stores the indices of atom images
          imags.resize(p.nRows);

          size_t iai{0};
          for (int z = -iimage[Z]; z <= iimage[Z]; ++z) { // serial
          for (int y = -iimage[Y]; y <= iimage[Y]; ++y) { // serial
          for (int x = -iimage[X]; x <= iimage[X]; ++x) { // serial
//            if (echo > 3) std::printf("# periodic shifts  %d %d %d\n", x, y, z);
              int const xyz_shift[] = {x, y, z};
              for (int ia = 0; ia < natoms; ++ia) { // loop over atoms in the unit cell, serial
                  // suggest a shifted atomic image position
                  double atom_pos[3];
                  for (int d = 0; d < 3; ++d) { // parallel
                      atom_pos[d] = xyzZinso[ia*8 + d] + xyz_shift[d]*cell[d];
                  } // d
                  auto const atom_id = int32_t(xyzZinso[ia*8 + 4]); 
                  auto const numax =       int(xyzZinso[ia*8 + 5]);
                  auto const sigma =           xyzZinso[ia*8 + 6] ;
//                   if (echo > 5) std::printf("# image of atom #%i at %g %g %g %s\n", atom_id, atom_pos[X]*Ang, atom_pos[Y]*Ang, atom_pos[Z]*Ang, _Ang);

                  double const r_projection = r_proj*sigma; // atom-dependent, precision dependent, assume float here
                  double const r2projection = pow2(r_projection);
//                double const r2projection_plus = pow2(r_projection + r_block_circumscribing_sphere);

                  // check all target blocks if they are inside the projection radius
                  uint32_t ntb{0}; // number of target blocks
                  for (uint32_t icube = 0; icube < p.nRows; ++icube) { // loop over blocks
                      auto const *const target_block = p.target_coords[icube];
//                    double d2{0};
//                    for (int d = 0; d < 3; ++d) { // serial
//                        double const center_of_block = (target_block[d]*4 + 1.5)*hg[d];
//                        d2 += pow2(center_of_block - atom_pos[d]);
//                    } // d
//                    if (d2 < r2projection_plus) {
                      if (1) {
                          // do more precise checking
//                           if (echo > 9) std::printf("# target block #%i at %i %i %i gets corner check\n",
//                                           icube, target_block[X], target_block[Y], target_block[Z]);
                          int nci{0}; // number of corners inside
                          // check 8 corners
                          for (int iz = 0; iz < 4; iz += 3) { // parallel, reduction on nci
                          for (int iy = 0; iy < 4; iy += 3) { // parallel, reduction on nci
                          for (int ix = 0; ix < 4; ix += 3) { // parallel, reduction on nci
                              int const ixyz[] = {ix, iy, iz};
                              double d2i{0}; // init distance^2 of the grid point from the center
                              for (int d = 0; d < 3; ++d) {
                                  double const grid_point = (target_block[d]*4 + ixyz[d] + 0.5)*hg[d];
                                  d2i += pow2(grid_point - atom_pos[d]);
                              } // d
                              if (d2i < r2projection) {
                                  ++nci; // at least one corner of the block is inside the projection radius of this atom
                              } // inside the projection radius
                          }}} // ix // iy // iz
                          // three different cases: 0, 1...7, 8
                          if (nci > 0) {
                              // atom image contributes
                              if (0 == ntb) {
                                  // this is the 1st cube to contribute
                                  std::vector<uint32_t> cube_list(0); // create an empty vector
                                  cubes.push_back(cube_list); // enlist the vector
                              } // 0 == ntb

                              cubes[iai].push_back(icube); // enlist
                              assert(cubes[iai][ntb] == icube); // check
                              imags[icube].push_back(iai); // also set up the transpose

                              ++ntb;
                              assert(cubes[iai].size() == ntb);
                              int const nCols = p.RowStart[icube + 1] - p.RowStart[icube];
                              sparse += nCols;
                              dense  += p.nCols; // all columns

//                               if (echo > 7) std::printf("# target block #%i at %i %i %i is inside\n",
//                                       icube, target_block[X], target_block[Y], target_block[Z]);
                          } else { // nci
//                               if (echo > 9) std::printf("# target block #%i at %i %i %i is outside\n",
//                                       icube, target_block[X], target_block[Y], target_block[Z]);
                          } // nci

                      } else { // d2 < r2projection_plus
//                           if (echo > 21) std::printf("# target block #%i at %i %i %i is far outside\n",
//                                       icube, target_block[X], target_block[Y], target_block[Z]);
                      } // d2 < r2projection_plus
                  } // icube

                  if (ntb > 0) {
                      // atom image contributes, mark in the list to have more than 0 coefficients
                      auto const nc = sho_tools::nSHO(numax); // number of coefficients for this atom
                      atom_ncoeff[ia] = nc; // atom image does contribute, so this atom does
                      assert(nc == atom_ncoeff[ia]); // conversion successful
                      nc_stats.add(nc);

                      // at least one target block has an intersection with the projection sphere of this atom image
                      auto & atom = atom_data[iai];
                      set(atom.pos, 3, atom_pos);
                      atom.sigma = sigma;
                      atom.gid = atom_id;
                      atom.ia = ia; // local atom index
                      set(atom.shifts, 3, xyz_shift);
                      atom.nc = nc; 
                      atom.numax = numax;

//                       if (echo > 5) std::printf("# image of atom #%i at %g %g %g %s contributes to %d target blocks\n",
//                                                    atom_id, atom_pos[X]*Ang, atom_pos[Y]*Ang, atom_pos[Z]*Ang, _Ang, ntb);
                      ApcStart[iai + 1] = ApcStart[iai] + nc;
                      ++iai;

                  } else {
//                       if (echo > 15) std::printf("# image of atom #%i at %g %g %g %s does not contribute\n",
//                                                     atom_id, atom_pos[X]*Ang, atom_pos[Y]*Ang, atom_pos[Z]*Ang, _Ang);
                  } // ntb > 0
              } // ia
          }}} // x // y // z

          auto const nai = iai; // corrected number of atomic images
          if (echo > 3) std::printf("# %ld of %lu (%.2f %%) atom images have an overlap with projection spheres\n",
                                    nai, natom_images, nai/(natom_images*.01));
          auto const napc = ApcStart[nai];

          if (echo > 3) std::printf("# sparse %g (%.2f %%) of dense %g\n", sparse, sparse/(dense*.01), dense);

          if (echo > 3) std::printf("# number of coefficients per image %.1f +/- %.1f in [%g, %g]\n",
                                    nc_stats.mean(), nc_stats.dev(), nc_stats.min(), nc_stats.max());

          if (echo > 3) std::printf("# %.3f k atomic projection coefficients, %.2f per atomic image\n", napc*.001, napc/double(nai));
          // projection coefficients for the non-local PAW operations are stored
          // as std::complex<real_t> apc[napc][nRHSs][Noco][Noco*64] or real_t apc[napc][nRHSs][2][Noco][Noco*64] on the GPU
          if (echo > 3) std::printf("# memory of atomic projection coefficients is %.6f %s (float, twice for double)\n",
                                                  napc*nRHSs*2.*64.*sizeof(float)*GByte, _GByte);

          p.natom_images = nai;
          assert(nai == p.natom_images); // verify

          { // scope: set up sparse tables
              using ::green_sparse::sparse_t;
              
              p.sparse_SHOadd = sparse_t<>(imags, false, "sparse_SHOadd", echo);
              // sparse_SHOadd: rows == Green function rows, cols == atom images
              if (echo > 29) {
                  std::printf("# sparse_SHOadd.rowStart(%p)= ", (void*)p.sparse_SHOadd.rowStart());
                  printf_vector(" %d", p.sparse_SHOadd.rowStart(), p.sparse_SHOadd.nRows());
                  std::printf("# sparse_SHOadd.colIndex(%p)= ", (void*)p.sparse_SHOadd.colIndex());
                  printf_vector(" %d", p.sparse_SHOadd.colIndex(), p.sparse_SHOadd.nNonzeros());
              } // echo
              imags.resize(0); // release host memory
             

              std::vector<std::vector<std::vector<uint32_t>>> vvv(p.nCols);
              for (uint16_t irhs = 0; irhs < p.nCols; ++irhs) {
                  vvv[irhs].resize(nai);
              } // irhs
              for (uint32_t iai = 0; iai < p.natom_images; ++iai) {
                  for (uint32_t itb = 0; itb < cubes[iai].size(); ++itb) {
                      auto const iRow = cubes[iai][itb];
                      for (auto inzb = p.RowStart[iRow]; inzb < p.RowStart[iRow + 1]; ++inzb) {
                          auto const irhs = p.colindx[inzb];
                          assert(irhs < p.nCols);
                          vvv[irhs][iai].push_back(inzb);
                      } // inzb
                  } // itb
              } // iai
              p.sparse_SHOprj = get_memory<sparse_t<>>(p.nCols, echo, "sparse_SHOprj");
              for (uint16_t irhs = 0; irhs < p.nCols; ++irhs) {
                  char name[32]; std::snprintf(name, 31, "sparse_SHOprj[irhs=%i]", irhs);
                  p.sparse_SHOprj[irhs] = sparse_t<>(vvv[irhs], false, name, echo);
              } // irhs

          } // scope: set up sparse tables

          assert(cubes.size() == p.natom_images);
          cubes.resize(0); // release host memory


          p.ApcStart = get_memory<uint32_t>(nai + 1, echo, "ApcStart");
          set(p.ApcStart, nai + 1, ApcStart.data()); // copy into GPU memory

          // get all info for the atomic matrices:
          std::vector<int32_t> global_atom_index(natoms); // translation table
          std::vector<int32_t>  local_atom_index(natoms, -1); // translation table
          int iac{0};
          for (int ia = 0; ia < natoms; ++ia) { // serial
              if (atom_ncoeff[ia] > 0) {
                  global_atom_index[iac] = ia;
                  local_atom_index[ia] = iac;
                  ++iac;
              } // atom contributes
          } // ia
          int const nac = iac; // number of contributing atoms
          global_atom_index.resize(nac);

          // now store the atomic positions in GPU memory
          p.atom_data = get_memory<green_action::atom_t>(nai, echo, "atom_data");
          p.AtomPos   = get_memory<double[3+1]>(nai, echo, "AtomPos");
          for (int iai = 0; iai < nai; ++iai) {
              p.atom_data[iai] = atom_data[iai]; // copy
              // translate index
              int const ia = atom_data[iai].ia;
              int const iac = local_atom_index[ia];
              assert(iac > -1);
              p.atom_data[iai].ia = iac;
              set(p.AtomPos[iai], 3, atom_data[iai].pos);
              p.AtomPos[iai][3] = 1./std::sqrt(atom_data[iai].sigma);
              // ToDo: transfer a list of numax for each atom image to GPU memory
          } // copy into GPU memory

          // get memory for the matrices and fill
          p.atom_mat = get_memory<double*>(nac, echo, "atom_mat");
          for (int iac = 0; iac < nac; ++iac) { // parallel
              int const ia = global_atom_index[iac];
              int const nc = atom_ncoeff[ia];
              assert(nc > 0); // the number of coefficients of contributing atoms must be non-zero
              char name[32]; std::snprintf(name, 31, "atom_mat[iac=%d/ia=%d]", iac, ia);
              p.atom_mat[iac] = get_memory<double>(Noco*Noco*2*nc*nc, echo, name);
              set(p.atom_mat[iac], Noco*Noco*2*nc*nc, 0.0); // clear
              // fill this with matrix values
              assert(2*nc*nc <= atom_mat[iac].size());
              // use MPI communication to find values in atom owner processes
              auto const hmt = atom_mat[ia].data();
              auto const ovl = hmt + nc*nc;
              for (int i = 0; i < nc; ++i) {
                  for (int j = 0; j < nc; ++j) {
                      int const ij = i*nc + j;
                      p.atom_mat[iac][ij]         = hmt[ij] - E_param.real() * ovl[ij]; // real part
                      p.atom_mat[iac][ij + nc*nc] =         - E_param.imag() * ovl[ij]; // imag part
                  } // j
              } // i
          } // iac
          p.number_of_contributing_atoms = nac;
          if (echo > 1) std::printf("# found %d contributing atoms and %ld atom images\n", nac, nai);
          assert(nac == nai && "So far, we do not distinguish between atoms and atom images");

      } // scope "Atom part"

      int const n_iterations = control::get("green_function.benchmark.iterations", -1.); 
                      // -1: no iterations, 0:run memory initialization only, >0: iterate
      if (n_iterations < 0) {
          if (echo > 2) std::printf("# green_function.benchmark.iterations=%d --> no benchmarks\n", n_iterations);
      } else { // n_iterations < 0
          // try one of the 6 combinations (strangely, we cannot run any two of these calls after each other, ToDo: find out what's wrong here)
          int const n3bit = control::get("green_function.benchmark.action", 0.);
          switch (n3bit) {
              case 0: try_action<float ,1,1>(p, n_iterations, echo); break; // real
              case 1: try_action<float ,2,1>(p, n_iterations, echo); break; // complex
              case 3: try_action<float ,2,2>(p, n_iterations, echo); break; // non-collinear
              case 4: try_action<double,1,1>(p, n_iterations, echo); break; // real
              case 5: try_action<double,2,1>(p, n_iterations, echo); break; // complex
              case 7: try_action<double,2,2>(p, n_iterations, echo); break; // non-collinear
              default: warn("green_function.benchmark.action must be in {0, 1, 3, 4, 5, 7} but found %d", n3bit);
          } // switch n3bit
      } // n_iterations < 0

      return 0;
  } // construct_Green_function



#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_Green_function(int const echo=0) {
      int    ng[3] = {0, 0, 0}; // grid sizes
      double hg[3] = {1, 1, 1}; // grid spacings
      std::vector<double> Veff(0); // local potential
      int natoms{0}; // number of atoms
      std::vector<double> xyzZinso(0); // atom info
      std::vector<std::vector<double>> atom_mat(0); // non-local potential

      auto const stat = green_input::load_Hamiltonian(ng, hg, Veff, natoms, xyzZinso, atom_mat, echo);
      if (stat) {
          warn("failed to load_Hamiltonian with status=%d", int(stat));
          return stat;
      } // stat

      green_action::plan_t p;
      return construct_Green_function(p, ng, hg, Veff, xyzZinso, atom_mat, echo);
  } // test_Green_function

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_Green_function(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_function
