#pragma once

#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cmath> // std::sqrt
#include <algorithm> // std::max
#include <utility> // std::swap //, std::move
#include <vector> // std::vector<T>
#include <cstdio> // std::printf


#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "simple_timer.hxx" // SimpleTimer

#ifdef  HAS_RAPIDXML
  // git clone https://github.com/dwd/rapidxml
  #include "rapidxml/rapidxml.hpp" // ::xml_document<>
  #include "rapidxml/rapidxml_utils.hpp" // ::file<>

  #include "xml_reading.hxx" // ::find_attribute, ::find_child
#endif // HAS_RAPIDXML

#include "xml_reading.hxx" // ::read_sequence

#include "unit_system.hxx" // eV, _eV, Ang, _Ang

#include "print_tools.hxx" // printf_vector
#include "global_coordinates.hxx" // ::get

#include "green_memory.hxx" // get_memory, free_memory
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

namespace green_function {

  double const GByte = 1e-9; char const *const _GByte = "GByte";

  inline status_t construct_Green_function(
        int const ng[3] // number of grid points of the unit cell in with the potential is defined
      , double const hg[3] // grid spacings
      , std::vector<double> const & Veff // [ng[2]*ng[1]*ng[0]]
      , std::vector<double> const & xyzZinso // [natoms*8]
      , std::vector<std::vector<double>> const & atom_mat // atomic hamiltonian and overlap matrix
      , int const echo=0 // log-level
      , std::complex<double> const *energy_parameter=nullptr // E in G = (H - E*S)^{-1}
  ) {
      int constexpr X=0, Y=1, Z=2;
      if (echo > 0) std::printf("\n#\n# %s(%i, %i, %i)\n#\n\n", __func__, ng[X], ng[Y], ng[Z]);

      std::complex<double> E_param(energy_parameter ? *energy_parameter : 0);

      auto const p = new green_action::plan_t(); // create a plan how to apply the SHO-PAW Hamiltonian to a block-sparse truncated Green function

      int32_t n_original_Veff_blocks[3] = {0, 0, 0};
      for (int d = 0; d < 3; ++d) {
          n_original_Veff_blocks[d] = (ng[d] >> 2); // divided by 4
          assert(ng[d] == 4*n_original_Veff_blocks[d] && "All grid dimensions must be a multiple of 4!");
      } // d

      if (echo > 3) { std::printf("# n_original_Veff_blocks "); printf_vector(" %d", n_original_Veff_blocks, 3); }

      assert(Veff.size() == ng[Z]*ng[Y]*ng[X]);

//    view3D<block4<double>> Vtot_blocks(n_original_Veff_blocks[Z], n_original_Veff_blocks[Y], n_original_Veff_blocks[X]);
      size_t const n_all_Veff_blocks = n_original_Veff_blocks[Z]*n_original_Veff_blocks[Y]*n_original_Veff_blocks[X];

//    auto Vtot_gpu = get_memory<double>(n_all_Veff_blocks*64); // does not really work flawlessly
//    view4D<double> Vtot(Vtot_gpu, n_original_Veff_blocks[Y], n_original_Veff_blocks[X], 64); // wrap
//    view4D<double> Vtot(n_original_Veff_blocks[Z], n_original_Veff_blocks[Y], n_original_Veff_blocks[X], 64); // get memory
      p->Veff = get_memory<double[64]>(n_all_Veff_blocks); // in managed memory
      { // scope: reorder Veff into block-structured p->Veff
          for (int ibz = 0; ibz < n_original_Veff_blocks[Z]; ++ibz) {
          for (int iby = 0; iby < n_original_Veff_blocks[Y]; ++iby) {
          for (int ibx = 0; ibx < n_original_Veff_blocks[X]; ++ibx) {
              int const ibxyz[3] = {ibx, iby, ibz};
              auto const Veff_index = index3D(n_original_Veff_blocks, ibxyz);
//            auto const Vtot_xyz = Vtot(ibz,iby,ibx); // get a view1D, i.e. double*
              auto const Vtot_xyz = p->Veff[Veff_index];
              for (int i4z = 0; i4z < 4; ++i4z) {
              for (int i4y = 0; i4y < 4; ++i4y) {
              for (int i4x = 0; i4x < 4; ++i4x) {
                  int const i64 = (i4z*4 + i4y)*4 + i4x;
                  size_t const izyx = ((ibz*4 + i4z)
                              *ng[Y] + (iby*4 + i4y)) 
                              *ng[X] + (ibx*4 + i4x);
                  assert(izyx < ng[Z]*ng[Y]*ng[X]);
                  Vtot_xyz[i64] = Veff[izyx];
              }}} // i4
          }}} // xyz
      } // scope

      // Cartesian cell parameters for the unit cell in which the potential is defined
      double const cell[3] = {ng[X]*hg[X], ng[Y]*hg[Y], ng[Z]*hg[Z]};

      // assume periodic boundary conditions and an infinite host crystal,
      // so there is no need to consider k-points
      
      // determine the largest and smallest indices of target blocks
      // given a max distance r_trunc between source blocks and target blocks

      double const r_block_circumscribing_sphere = 1.5*std::sqrt(pow2(hg[X]) + pow2(hg[Y]) + pow2(hg[Z]));
      if (echo > 0) std::printf("# circumscribing radius= %g %s\n", r_block_circumscribing_sphere*Ang, _Ang);

      // assume that the source blocks lie compact in space
      int const nRHSs = n_original_Veff_blocks[Z] * n_original_Veff_blocks[Y] * n_original_Veff_blocks[X];
      if (echo > 0) std::printf("# total number of source blocks is %d\n", nRHSs);
      view2D<uint32_t> global_source_coords(nRHSs, 4, 0); // unsigned since they must be in [0, 2^21) anyway
      p->global_source_indices.resize(nRHSs); // [nRHSs]
      p->nCols = nRHSs;
      double center_of_mass_RHS[3] = {0, 0, 0};
      double center_of_RHSs[3]     = {0, 0, 0};
      int32_t min_source_coords[3] = {0, 0, 0}; // global coordinates
      int32_t max_source_coords[3] = {0, 0, 0}; // global coordinates
      int32_t internal_global_offset[3] = {0, 0, 0};
      double max_distance_from_comass{0};
      double max_distance_from_center{0};
      { // scope: determine min, max, center
          {   int iRHS{0};
              for (int ibz = 0; ibz < n_original_Veff_blocks[Z]; ++ibz) {
              for (int iby = 0; iby < n_original_Veff_blocks[Y]; ++iby) {
              for (int ibx = 0; ibx < n_original_Veff_blocks[X]; ++ibx) {
                  global_source_coords(iRHS,X) = ibx;
                  global_source_coords(iRHS,Y) = iby;
                  global_source_coords(iRHS,Z) = ibz;
                  global_source_coords(iRHS,3) = 0; // unused
                  p->global_source_indices[iRHS] = global_coordinates::get(ibx, iby, ibz);
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

          p->source_coords = get_memory<int16_t[4]>(nRHSs); // internal coordindates
          { // scope: compute also the largest distance from the center or center of mass
              double max_d2m{0}, max_d2c{0};
              for (int iRHS = 0; iRHS < nRHSs; ++iRHS) {
                  double d2m{0}, d2c{0};
                  p->source_coords[iRHS][3] = 0; // not used
                  for (int d = 0; d < 3; ++d) {
                      p->source_coords[iRHS][d] = global_source_coords(iRHS,d) - internal_global_offset[d];
                      d2m += pow2((global_source_coords(iRHS,d)*4 + 1.5)*hg[d] - center_of_mass_RHS[d]);
                      d2c += pow2((global_source_coords(iRHS,d)*4 + 1.5)*hg[d] - center_of_RHSs[d]);
                  } // d
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
      double const r_trunc = control::get("green.function.truncation.radius", 10.);
      if (echo > 0) std::printf("# green.function.truncation.radius=%g %s\n", r_trunc*Ang, _Ang);
      p->r_truncation = std::max(0., r_trunc);
      p->r_Vconfinement = std::min(std::max(0., r_trunc - 2.0), p->r_truncation);

// example:
//    truncation radius in Cu (fcc)
//    lattice constant = 3.522 Ang
//    volume per atom = alat^3 / 4 == 10.922 Ang^3 / atom
//    1000 atoms inside the truncation cluster
//    truncation sphere volume 10922 Ang^3
//    truncation radius = cbrt(3*V/(4*pi)) = 13.764 Ang = 26 Bohr
//    8000 atoms --> 52 Bohr
//    64000 atoms --> 104 Bohr


      // count the number of green function elements for each target block
      size_t nnz{0}; // number of non-zero BSR entries

      uint16_t num_target_coords[3] = {0, 0, 0};
      int16_t  min_target_coords[3] = {0, 0, 0}; // internal coordinates
      int16_t  max_target_coords[3] = {0, 0, 0}; // internal coordinates
      { // scope: create the truncated Green function block-sparsity pattern
          auto const rtrunc       = std::max(0., r_trunc);
          auto const rtrunc_plus  =              rtrunc + 2*r_block_circumscribing_sphere;
          auto const rtrunc_minus = std::max(0., rtrunc - 2*r_block_circumscribing_sphere);
          if (echo > 0) std::printf("# truncation radius %g, search within %g %s\n", rtrunc*Ang, rtrunc_plus*Ang, _Ang);
          if (echo > 0 && rtrunc_minus > 0) std::printf("# blocks with center distance below %g %s are fully inside\n", rtrunc_minus*Ang, _Ang);

          int16_t itr[3]; // range [-32768, 32767] should be enough
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
          std::vector<std::vector<uint16_t>> column_indices(product_target_blocks);

          double const r2trunc        = pow2(rtrunc),
                       r2trunc_plus   = pow2(rtrunc_plus),
                       r2trunc_minus  = pow2(rtrunc_minus);
          double const r2block_circum = pow2(r_block_circumscribing_sphere*3);

          std::vector<int32_t> tag_diagonal(product_target_blocks, -1);
          assert(nRHSs < 65536 && "the integer type of ColIndex is uint16_t!");
          std::vector<std::vector<bool>> sparsity_pattern(nRHSs); // vector<bool> is a memory-saving bit-array

          for (uint16_t iRHS = 0; iRHS < nRHSs; ++iRHS) {
              sparsity_pattern[iRHS] = std::vector<bool>(product_target_blocks, false);
              auto & sparsity_RHS = sparsity_pattern[iRHS]; // abbreviate
              auto const *const source_coords = p->source_coords[iRHS]; // internal source block coordinates
              simple_stats::Stats<> stats[3];
              int constexpr max_nci = 27;
              std::vector<int> hist(1 + max_nci, 0); // distribution of nci
              simple_stats::Stats<> stats_d2[1 + max_nci];
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
                  // distance^2 of the block centers
                  double const d2 = pow2(bx*4*hg[X]) + pow2(by*4*hg[Y]) + pow2(bz*4*hg[Z]);

                  uint8_t nci{0}; // init number of corners inside
                  if (d2 < r2trunc_plus) { // potentially inside, check all 8 corner cases
                      int const far = (d2 > r2block_circum);
//                    if (d2 < r2trunc_minus) { nci = max_nci; } else // skip the 8-corners test for inner blocks -> some speedup
                      { // scope: 8 corner test
                         // i = i4 - j4 --> i in [-3, 3], 
                         //     if two blocks are far from each other, we test only the 8 combinations of {-3, 3}
                         //     for blocks close to each other, we test all 27 combinations of {-3, 0, 3}
                          for (int iz = -3; iz <= 3; iz += 3 + 3*far) {
                          for (int iy = -3; iy <= 3; iy += 3 + 3*far) {
                          for (int ix = -3; ix <= 3; ix += 3 + 3*far) {
                              double const d2c = pow2((bx*4 + ix)*hg[X])
                                               + pow2((by*4 + iy)*hg[Y])
                                               + pow2((bz*4 + iz)*hg[Z]);
//                               if (0 == iRHS && (d2c < r2trunc)) {
//                                   std::printf("# %s: b= %i %i %i, i-j %i %i %i, d^2= %g %s\n", 
//                                       __func__, bx,by,bz, ix,iy,iz, d2c, (d2c < r2trunc)?"in":"out");
//                               }
                              nci += (d2c < r2trunc);
                          }}}
                      } // scope

                      if (d2 < r2trunc_minus) assert(27 + far*(8 - 27) == nci); // for these, we could skip the 8-corners test
                  } // d2

                  if (nci > 0) { 
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
                                  hist[nci]*.001, nci, stats_d2[nci].avg(), stats_d2[nci].var(), stats_d2[nci].min(), stats_d2[nci].max());
                              total_checked += hist[nci];
                          } // hist[nci] > 0
                      } // nci
                      auto const partial = total_checked - hist[0] - hist[max_nci];
                      std::printf("# RHS has %.3f k inside, %.3f k partial and %.3f k outside (of %.3f k checked blocks)\n",
                                    hist[max_nci]*.001, partial*.001, hist[0]*.001, total_checked*.001);
                  } // RHS #0
              } // echo
          } // iRHS

          // create a histogram about the distribution of number of columns per row
          std::vector<uint32_t> hist(nRHSs + 1, 0);
          for (size_t idx3 = 0; idx3 < column_indices.size(); ++idx3) {
              auto const n = column_indices[idx3].size();
              ++hist[n];
          } // idx3

          // eval the histogram
          size_t nall{0};
          for (int n = 0; n <= nRHSs; ++n) {
              nall += hist[n];
              nnz  += hist[n]*n;
          } // n
          if (echo > 5) { std::printf("# histogram total=%.3f k: ", nall*.001); printf_vector(" %d", hist.data(), nRHSs + 1); }
          assert(nall == product_target_blocks); // sanity check

          p->nRows = nall - hist[0]; // the target block entries with no RHS do not create a row
          if (echo > 0) std::printf("# total number of Green function blocks is %.3f k, "
                               "average %.1f per source block\n", nnz*.001, nnz/double(nRHSs));
          if (echo > 0) std::printf("# %.3f k (%.1f %% of %.3f k) target blocks are active\n", 
              p->nRows*.001, p->nRows/(product_target_blocks*.01), product_target_blocks*.001);

          assert(nnz < (uint64_t(1) << 32) && "the integer type or RowStart is uint32_t!");

          // resize BSR tables: (Block-compressed Sparse Row format)
          if (echo > 3) { std::printf("# memory of Green function is %.6f %s (float, twice for double)\n",
                              nnz*2.*64.*64.*sizeof(float)*GByte, _GByte); std::fflush(stdout); }
          auto & ColIndex = p->colindx;
          ColIndex.resize(nnz);
          p->rowindx = get_memory<uint32_t>(nnz);
          auto & RowStart = p->RowStart;
          RowStart = get_memory<uint32_t>(p->nRows + 1);
          RowStart[0] = 0;
          p->veff_index = get_memory<uint32_t>(p->nRows); // indirection list for the local potential
          p->target_coords = get_memory<int16_t[4]>(p->nRows); // view2D<int16_t>(p->nRows, 4, 0);
          p->global_target_indices.resize(p->nRows);
          p->subset.resize(p->nCols); // we assume columns of the unit operator as RHS

          view3D<int32_t> iRow_of_coords(num_target_coords[Z],
                                         num_target_coords[Y],
                                         num_target_coords[X], -1);

          { // scope: fill BSR tables
              simple_stats::Stats<> st;
              uint32_t iRow{0};
              for (uint16_t z = 0; z < num_target_coords[Z]; ++z) { // serial
              for (uint16_t y = 0; y < num_target_coords[Y]; ++y) { // serial
              for (uint16_t x = 0; x < num_target_coords[X]; ++x) { // serial
                  uint16_t const idx[3] = {x, y, z};
                  auto const idx3 = index3D(num_target_coords, idx);
                  assert(idx3 < product_target_blocks);

                  auto const n = column_indices[idx3].size();
                  if (n > 0) {
                      st.add(n);
                      iRow_of_coords(idx[Z], idx[Y], idx[X]) = iRow;

                      RowStart[iRow + 1] = RowStart[iRow] + n;
                      // copy the column indices
                      set(ColIndex.data() + RowStart[iRow], n, column_indices[idx3].data());
                      set(p->rowindx      + RowStart[iRow], n, iRow);
                      // copy the target block coordinates
                      int32_t global_target_coords[3];
                      for (int d = 0; d < 3; ++d) {
                          p->target_coords[iRow][d] = idx[d] + min_target_coords[d];
                          global_target_coords[d] = p->target_coords[iRow][d] + internal_global_offset[d];
                      } // d
                      p->target_coords[iRow][3] = 0; // not used

                      p->global_target_indices[iRow] = global_coordinates::get(global_target_coords); 
                      // global_target_indices is needed to gather the local potential data from other MPI processes

                      { // scope: determine the diagonal entry (source == target)
                          auto const iCol = tag_diagonal[idx3];
                          if (iCol > -1) {
                              assert(iCol < (1ul << 16)); // number range if uint16_t
                              for (int d = 0; d < 3; ++d) {
                                  // sanity check onto internal coordinates
                                  assert(p->source_coords[iCol][d] == p->target_coords[iRow][d]);
                              } // d
                              auto inz = RowStart[iRow];
                              while(iCol != ColIndex[inz]) { ++inz; }
                              p->subset[iCol] = inz;
                          } // iCol valid
                      } // scope

                      if (1) { // scope: fill indirection table for having the local potential only defined in 1 unit cell and repeated periodically
                          int32_t mod[3];
                          for (int d = 0; d < 3; ++d) {
                              mod[d] = global_target_coords[d] % n_original_Veff_blocks[d];
                              mod[d] += (mod[d] < 0)*n_original_Veff_blocks[d];
                          } // d
                          p->veff_index[iRow] = index3D(n_original_Veff_blocks, mod);
                      } // scope
 
                      // count up the number of active rows
                      ++iRow;
                  } // n > 0
              }}} // idx
              assert(p->nRows == iRow);
              assert(nnz == RowStart[p->nRows]);
              std::printf("# source blocks per target block: average %.1f +/- %.1f in [%g, %g]\n", st.avg(), st.var(), st.min(), st.max());
          } // scope
          column_indices.clear(); // not needed beyond this point

          if (echo > 1) { // measure the difference in the number of target blocks of each RHS
              std::vector<uint32_t> nt(nRHSs, 0);
              // traverse the BSR structure
              for (uint32_t iRow = 0; iRow < p->nRows; ++iRow) {
                  for (auto inz = RowStart[iRow]; inz < RowStart[iRow + 1]; ++inz) {
                      auto const iCol = ColIndex[inz];
                      ++nt[iCol];
                  } // inz
              } // iRow
              // analyze nt
              simple_stats::Stats<> st;
              for (uint16_t iRHS = 0; iRHS < nRHSs; ++iRHS) {
                  st.add(nt[iRHS]);
              } // iRHS
              std::printf("# target blocks per source block: average %.1f +/- %.1f in [%g, %g]\n", st.avg(), st.var(), st.min(), st.max());
          } // echo

          // Green function is stored sparse 
          // as std::complex<real_t> green[nnz][64][64] 
          // or real_t green[nnz][2][64][64] for the GPU;

          for (int dd = 0; dd < 3; ++dd) { // derivate direction
              // create lists for the finite-difference derivatives
              p->fd_plan[dd] = green_kinetic::finite_difference_plan_t(dd
                , num_target_coords
                , RowStart, ColIndex.data()
                , iRow_of_coords
                , sparsity_pattern.data()
                , nRHSs, echo);
          } // dd

          
          // transfer grid spacing into GPU memory
          p->grid_spacing = get_memory<double>(4);
          set(p->grid_spacing, 3, hg);

      } // scope







      int const natoms = atom_mat.size();
      if (echo > 2) std::printf("\n#\n# %s: Start atom part, %d atoms\n#\n", __func__, natoms);

      // compute which atoms will contribute, the list of natoms atoms may contain a subset of all atoms
      double max_projection_radius{0};
      for (int ia = 0; ia < natoms; ++ia) {
          auto const sigma = xyzZinso[ia*8 + 6];
          auto const projection_radius = std::max(0.0, 6*sigma);
          max_projection_radius = std::max(max_projection_radius, projection_radius);
      } // ia
      if (echo > 3) std::printf("# largest projection radius is %g %s\n", max_projection_radius*Ang, _Ang);


      p->ApcStart = nullptr;
      p->natom_images = 0;
      { // scope:
          SimpleTimer atom_list_timer(__FILE__, __LINE__, "Atom part");
          
          auto const radius = r_trunc + max_distance_from_center 
                            + 2*max_projection_radius + 2*r_block_circumscribing_sphere;
          int iimage[3];
          size_t nimages{1};
          for (int d = 0; d < 3; ++d) { // parallel
              iimage[d] = int(std::ceil(radius/cell[d]));
              nimages *= (iimage[d]*2 + 1);
          } // d
          auto const natom_images = natoms*nimages;
          if (echo > 3) std::printf("# replicate %d %d %d atom images, %.3f k total\n", 
                                  iimage[X], iimage[Y], iimage[Z], nimages*.001);

          std::vector<uint32_t> ApcStart(natom_images + 1, 0); // probably larger than needed, should call resize(nai + 1) later
          std::vector<green_action::atom_t> atom_data(natom_images);
          std::vector<uint8_t> atom_ncoeff(natoms, 0); // 0: atom does not contribute

          simple_stats::Stats<> nc_stats;
          uint32_t constexpr COUNT = 0; // 0:count how many blocks are really involved
          double sparse{0}, dense{0}; // stats to assess how much memory can be saved using sparse storage

          size_t iai{0};
          for (int z = -iimage[Z]; z <= iimage[Z]; ++z) { // serial
          for (int y = -iimage[Y]; y <= iimage[Y]; ++y) { // serial
          for (int x = -iimage[X]; x <= iimage[X]; ++x) { // serial
//            if (echo > 3) std::printf("# periodic shifts  %d %d %d\n", x, y, z);
              int const xyz[3] = {x, y, z};
              for (int ia = 0; ia < natoms; ++ia) { // loop over atoms in the unit cell, serial
                  // suggest a new atomic image position
                  double pos[3];
                  for (int d = 0; d < 3; ++d) { // parallel
                      pos[d] = xyzZinso[ia*8 + d] + xyz[d]*cell[d];
                  } // d
                  auto const atom_id = int32_t(xyzZinso[ia*8 + 4]); 
                  auto const numax =       int(xyzZinso[ia*8 + 5]);
                  auto const sigma =           xyzZinso[ia*8 + 6];
//                   if (echo > 5) std::printf("# image of atom #%i at %g %g %g %s\n", atom_id, pos[X]*Ang, pos[Y]*Ang, pos[Z]*Ang, _Ang);

                  double const r_projection = pow2(6*sigma); // atom-dependent, precision dependent, assume float here
                  double const r2projection = pow2(r_projection);
//                double const r2projection_plus = pow2(r_projection + r_block_circumscribing_sphere);

                  // check all target blocks if they are inside the projection radius
                  uint32_t ntb{0}; // number of target blocks
                  for (uint32_t iRow = 0; (iRow < p->nRows) && (0 == COUNT*ntb); ++iRow) {
                      auto const *const target_block = p->target_coords[iRow];
                      double d2{0};
                      for (int d = 0; d < 3; ++d) { // serial
                          double const center_of_block = (target_block[d]*4 + 1.5)*hg[d];
                          d2 += pow2(center_of_block - pos[d]);
                      } // d
//                    if (d2 < r2projection_plus) {
                      if (1) {
                          // do more precise checking
//                           if (echo > 9) std::printf("# target block #%i at %i %i %i gets corner check\n",
//                                           iRow, target_block[X], target_block[Y], target_block[Z]);
                          int nci{0}; // number of corners inside
                          // check 8 corners
                          for (int iz = 0; iz < 4; iz += 3) { // parallel, reduction
                          for (int iy = 0; iy < 4; iy += 3) { // parallel, reduction
                          for (int ix = 0; ix < 4; ix += 3) { // parallel, reduction
                              int const ixyz[3] = {ix, iy, iz};
                              double d2i{0};
                              for (int d = 0; d < 3; ++d) {
                                  double const grid_point = (target_block[d]*4 + ixyz[d])*hg[d];
                                  d2i += pow2(grid_point - pos[d]);
                              } // d
                              if (d2i < r2projection) {
                                  ++nci; // at least one corner of the block 
                                  // ... is inside the projection radius of this atom
                              } // inside the projection radius
                          }}} // ixyz
                          // three different cases: 0, 1...7, 8
                          if (nci > 0) {
                              // atom image contributes
                              ++ntb; // stop outer loop if COUNT==1
                              int const nCols = p->RowStart[iRow + 1] - p->RowStart[iRow];
                              sparse += nCols;
                              dense  += p->nCols; // all columns
                              
//                               if (echo > 7) std::printf("# target block #%i at %i %i %i is inside\n",
//                                       iRow, target_block[X], target_block[Y], target_block[Z]);
                          } else { // nci
//                               if (echo > 9) std::printf("# target block #%i at %i %i %i is outside\n",
//                                       iRow, target_block[X], target_block[Y], target_block[Z]);
                          } // nci

                      } else { // d2
//                           if (echo > 21) std::printf("# target block #%i at %i %i %i is far outside\n",
//                                       iRow, target_block[X], target_block[Y], target_block[Z]);
                      } // d2
                  } // iRow

                  if (ntb > 0) {
                      // atom image contributes, mark in the list to have more than 0 coefficients
                      auto const nc = sho_tools::nSHO(numax);
                      atom_ncoeff[ia] = nc;
                      assert(nc == atom_ncoeff[ia]); // conversion successful
                      nc_stats.add(nc);

                      // at least one target block has an intersection with the projection sphere of this atom image
                      auto & atom = atom_data[iai];
                      set(atom.pos, 3, pos);
                      atom.sigma = sigma;
                      atom.gid = atom_id;
                      atom.ia = ia;
                      set(atom.shifts, 3, xyz);
                      atom.nc = nc; 
                      atom.numax = numax;

                      
//                       if (echo > 5) std::printf("# image of atom #%i at %g %g %g %s contributes to %d target blocks\n",
//                                                 atom_id, pos[X]*Ang, pos[Y]*Ang, pos[Z]*Ang, _Ang, ntb);
                      ApcStart[iai + 1] = ApcStart[iai] + atom.nc;
                      ++iai;
                  } else {
//                       if (echo > 15) std::printf("# image of atom #%i at %g %g %g %s does not contribute\n",
//                                                 atom_id, pos[X]*Ang, pos[Y]*Ang, pos[Z]*Ang, _Ang);
                  } // ntb > 0
              } // ia
          }}} // xyz

          auto const nai = iai; // corrected number of atomic images
          if (echo > 3) std::printf("# %ld of %lu (%.2f %%) atom images have an overlap with projection spheres\n",
                                  nai, natom_images, nai/(natom_images*.01));
          auto const napc = ApcStart[nai];

          if (echo > 3 && 0 == COUNT) std::printf("# sparse %g and dense %g\n", sparse, dense);

          if (echo > 3) std::printf("# number of coefficients per image %.1f +/- %.1f in [%g, %g]\n",
                                nc_stats.avg(), nc_stats.var(), nc_stats.min(), nc_stats.max());
          
          if (echo > 3) std::printf("# %.3f k atomic projection coefficients, %.2f per atomic image\n", napc*.001, napc/double(nai));
          // projection coefficients for the non-local PAW operations are stored
          // as std::complex<real_t> apc[napc][nRHSs][64]
          // or real_t apc[napc][nRHSs][2][64] for the GPU
          if (echo > 3) std::printf("# memory of atomic projection coefficients is %.6f %s (float, twice for double)\n",
                                  napc*nRHSs*2.*64.*sizeof(float)*GByte, _GByte);
          p->natom_images = nai;
          p->ApcStart = get_memory<uint32_t>(nai + 1);
          set(p->ApcStart, nai + 1, ApcStart.data()); // copy into GPU memory

          // get all info for the atomic matrices:
          std::vector<int32_t> global_atom_index(natoms, -1); // translation table
          std::vector<int32_t>  local_atom_index(natoms, -1); // translation table
          int iac{0};
          for (int ia = 0; ia < natoms; ++ia) { // serial
              if (atom_ncoeff[ia] > 0) {
                  global_atom_index[iac] = ia;
                  local_atom_index[ia] = iac;
                  ++iac;
              } // atom_contributes
          } // ia
          int const nac = iac;
          global_atom_index.resize(nac);

          // now store the atomic positions in GPU memory
          p->atom_data = get_memory<green_action::atom_t>(nai);
          for (int iai = 0; iai < nai; ++iai) {
              p->atom_data[iai] = atom_data[iai]; // copy
              // translate index
              int const ia = atom_data[iai].ia;
              int const iac = local_atom_index[ia];
              assert(iac > -1);
              p->atom_data[iai].ia = iac;
          } // copy into GPU memory

          // get memory for the matrices and fill
          p->atom_mat = get_memory<double(*)[2]>(nac);
          for (int iac = 0; iac < nac; ++iac) { // parallel
              int const ia = global_atom_index[iac];
              int const nc = atom_ncoeff[ia];
              assert(nc > 0); // the number of coefficients of contributing atoms must be non-zero
              p->atom_mat[iac] = get_memory<double[2]>(nc*nc);
              // fill this with matrix values
              // use MPI communication to find values in atom owner processes
              auto const hmt = atom_mat[ia].data();
              auto const ovl = hmt + nc*nc;
              for (int i = 0; i < nc; ++i) {
                  for (int j = 0; j < nc; ++j) {
                      int const ij = i*nc + j;
                      p->atom_mat[iac][ij][0] = hmt[ij] - E_param.real() * ovl[ij];
                      p->atom_mat[iac][ij][1] =         - E_param.imag() * ovl[ij];
                  } // j
              } // i
          } // iac

      } // scope

      int const n_iterations = control::get("green.function.benchmark.iterations", -1.); 
                      // -1: no iterations, 0:run memory initilziation only, >0: iterate
      if (n_iterations >= 0) { // scope: try out the operator
          typedef float real_t;
          int constexpr LM = 64;
          green_action::action_t<real_t,LM> action(p);
          auto const nnzbX = p->colindx.size();
          auto x = get_memory<real_t[2][LM][LM]>(nnzbX);
          auto y = get_memory<real_t[2][LM][LM]>(nnzbX);
          for (size_t i = 0; i < nnzbX*2*LM*LM; ++i) {
              x[0][0][0][i] = 0; // init x
          } // i

          // benchmark the action
          for (int iteration = 0; iteration < n_iterations; ++iteration) {
              if (echo > 5) { std::printf("# iteration #%i\n", iteration); std::fflush(stdout); }
              action.multiply(y, x, p->colindx.data(), nnzbX, p->nCols);
              std::swap(x, y);
          } // iteration

      } else {
          if (echo > 2) std::printf("# green.function.benchmark.iterations=%d --> no benchmarks\n", n_iterations);
      } // n_iterations >= 0

      return 0;
  } // construct_Green_function
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_Green_function(int const echo=0) {
#ifndef HAS_RAPIDXML
      warn("Unable to test Green function construction when compiled without -D HAS_RAPIDXML", 0);
      return STATUS_TEST_NOT_INCLUDED;
#else  // HAS_RAPIDXML
      char const *filename = "Hmt.xml";
      rapidxml::file<> infile(filename);

      rapidxml::xml_document<> doc;
      doc.parse<0>(infile.data());

      double hg[3] = {1, 1, 1};
      int ng[3] = {0, 0, 0};
      std::vector<double> Veff;
      std::vector<double> xyzZinso;
      std::vector<std::vector<double>> atom_mat;
      int natoms{0};

      auto const grid_Hamiltonian = doc.first_node("grid_Hamiltonian");
      if (grid_Hamiltonian) {
          auto const sho_atoms = xml_reading::find_child(grid_Hamiltonian, "sho_atoms", echo);
          if (sho_atoms) {
              auto const number = xml_reading::find_attribute(sho_atoms, "number", "0", echo);
              if (echo > 5) std::printf("# found number=%s\n", number);
              natoms = std::atoi(number);
              xyzZinso.resize(natoms*8);
              atom_mat.resize(natoms);
              int ia{0};
              for (auto atom = sho_atoms->first_node(); atom; atom = atom->next_sibling()) {
                  auto const gid = xml_reading::find_attribute(atom, "gid", "-1");
                  if (echo > 5) std::printf("# <%s gid=%s>\n", atom->name(), gid);
                  xyzZinso[ia*8 + 4] = std::atoi(gid);

                  double pos[3] = {0, 0, 0};
                  auto const position = xml_reading::find_child(atom, "position", echo);
                  for (int d = 0; d < 3; ++d) {
                      char axyz[] = {0, 0}; axyz[0] = 'x' + d; // "x", "y", "z"
                      auto const value = xml_reading::find_attribute(position, axyz);
                      if (*value != '\0') {
                          pos[d] = std::atof(value);
                          if (echo > 5) std::printf("# %s = %.15g\n", axyz, pos[d]);
                      } // value != ""
                      xyzZinso[ia*8 + d] = pos[d];
                  } // d

                  auto const projectors = xml_reading::find_child(atom, "projectors", echo);
                  int numax{-1};
                  {
                      auto const value = xml_reading::find_attribute(projectors, "numax", "-1");
                      if (*value != '\0') {
                          numax = std::atoi(value);
                          if (echo > 5) std::printf("# numax= %d\n", numax);
                      } // value != ""
                  }
                  double sigma{-1};
                  {
                      auto const value = xml_reading::find_attribute(projectors, "sigma", "-1");
                      if (*value != '\0') {
                          sigma = std::atof(value);
                          if (echo > 5) std::printf("# sigma= %g\n", sigma);
                      } // value != ""
                  }
                  xyzZinso[ia*8 + 5] = numax;
                  xyzZinso[ia*8 + 6] = sigma;
                  int const nSHO = sho_tools::nSHO(numax);
                  atom_mat[ia].resize(2*nSHO*nSHO);
                  for (int h0s1 = 0; h0s1 < 2; ++h0s1) {
                      auto const matrix_name = h0s1 ? "overlap" : "hamiltonian";
                      auto const matrix = xml_reading::find_child(atom, matrix_name, echo);
                      if (matrix) {
                          if (echo > 22) std::printf("# %s.values= %s\n", matrix_name, matrix->value());
                          auto const v = xml_reading::read_sequence<double>(matrix->value(), echo, nSHO*nSHO);
                          if (echo > 5) std::printf("# %s matrix has %ld values, expect %d x %d = %d\n",
                              matrix_name, v.size(), nSHO, nSHO, nSHO*nSHO);
                          assert(v.size() == nSHO*nSHO);
                          set(atom_mat[ia].data() + h0s1*nSHO*nSHO, nSHO*nSHO, v.data()); // copy
                      } else warn("atom with global_id=%s has no %s matrix!", gid, matrix_name);
                  } // h0s1
                  ++ia; // count up the number of atoms
              } // atom
              assert(natoms == ia); // sanity
          } else warn("no <sho_atoms> found in grid_Hamiltonian in file %s", filename);

          auto const spacing = xml_reading::find_child(grid_Hamiltonian, "spacing", echo);
          for (int d = 0; d < 3; ++d) {
              char axyz[] = {0, 0}; axyz[0] = 'x' + d; // "x", "y", "z"
              auto const value = xml_reading::find_attribute(spacing, axyz);
              if (*value != '\0') {
                  hg[d] = std::atof(value);
                  if (echo > 5) std::printf("# h%s = %.15g\n", axyz, hg[d]);
              } // value != ""
          } // d

          auto const potential = xml_reading::find_child(grid_Hamiltonian, "potential", echo);
          if (potential) {
              for (int d = 0; d < 3; ++d) {
                  char axyz[] = {'n', 0, 0}; axyz[1] = 'x' + d; // "nx", "ny", "nz"
                  auto const value = xml_reading::find_attribute(potential, axyz);
                  if (*value != '\0') {
                      ng[d] = std::atoi(value);
                      if (echo > 5) std::printf("# %s = %d\n", axyz, ng[d]);
                  } // value != ""
              } // d
              if (echo > 33) std::printf("# potential.values= %s\n", potential->value());
              Veff = xml_reading::read_sequence<double>(potential->value(), echo, ng[2]*ng[1]*ng[0]);
              if (echo > 5) std::printf("# potential has %ld values, expect %d x %d x %d = %d\n",
                  Veff.size(), ng[0], ng[1], ng[2], ng[2]*ng[1]*ng[0]);
              assert(Veff.size() == ng[2]*ng[1]*ng[0]);
          } else warn("grid_Hamiltonian has no potential in file %s", filename);

      } else warn("no grid_Hamiltonian found in file %s", filename);

      return construct_Green_function(ng, hg, Veff, xyzZinso, atom_mat, echo);
#endif // HAS_RAPIDXML
  } // test_Green_function

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_Green_function(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_function
