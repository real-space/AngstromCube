#pragma once

#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cmath> // std::sqrt
#include <algorithm> // std::max
#include <vector> // std::vector<T>

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "recorded_warnings.hxx" // warn
#include "simple_math.hxx" // ::random
#include "print_tools.hxx" // printf_vector
#include "data_view.hxx" // view4D<T>, view2D<T>
#include "display_units.h" // Ang, _Ang
#include "inline_math.hxx" // pow2
#include "control.hxx" // ::get
#include "simple_stats.hxx" // ::Stats
#include "sho_tools.hxx" // ::nSHO

namespace green_function {
  
  double const GByte = 1e-9; char const *const _GByte = "GByte";
  
  inline int64_t global_coordinates(int32_t x, int32_t y, int32_t z) {
      // interleaved bit pattern of the lowest 21 bits to 63 bits
      int64_t i63{0};
      for(int b21 = 0; b21 < 21; ++b21) {
          int32_t const i3 = (x & 0x1) + 2*(y & 0x1) + 4*(z & 0x1);
          x >>= 1; y >>= 1; z >>= 1; // delete the least significant bit of the coordinates
          int64_t const i3_shifted = int64_t(i3) << (b21*3);
          i63 |= i3_shifted;
      } // b21
      return i63; // when printing i63, use format %22.22lo
  } // global_coordinates

  template <typename int_t>
  inline int64_t global_coordinates(int_t const xyz[3]) {
      return global_coordinates(xyz[0], xyz[1], xyz[2]);
  } // global_coordinates

  inline status_t global_coordinates(uint32_t xyz[3], int64_t i63) {
      // retrieve global_coordinates from i63 identifyer
      // mind that the coordinates will be cast into the half-open range [0, 2^21)
      uint32_t x{0}, y{0}, z{0};
      for(int b21 = 0; b21 < 21; ++b21) {
          x |= (i63 & 0x1) << b21;   i63 >>= 1;
          y |= (i63 & 0x1) << b21;   i63 >>= 1;
          z |= (i63 & 0x1) << b21;   i63 >>= 1;
      } // b21
      xyz[0] = x; xyz[1] = y; xyz[2] = z;
      return status_t(i63 & 0x1); // this is 0 if i63 >= 0
  } // global_coordinates

  template <typename T>
  class block4 {
  public:
      block4(T const initial_value=T(0)) {
          for(int i = 0; i < 64; ++i) { _data[0][0][i] = initial_value; }
      } // constructor

      // access operators (i2,i1,i0) with i0,i1,i2 in [0, 3]
      T const & operator () (int i2, int i1, int i0) const { return _data[i2 & 0x3][i1 & 0x3][i0 & 0x3]; }
      T       & operator () (int i2, int i1, int i0)       { return _data[i2 & 0x3][i1 & 0x3][i0 & 0x3]; }
      // access operators [i] with i in [0, 63]
      T const & operator [] (int i64) const { return _data[0][0][i64 & 0x3f]; }
      T       & operator [] (int i64)       { return _data[0][0][i64 & 0x3f]; }

  private:
      T _data[4][4][4];
  }; // class

//   typedef struct {
//       double pos[3];
//       int32_t gid;
//       int8_t shifts[3];
//       int8_t iZ;
//   } atom_t;

  template <typename uint_t, typename int_t>
  inline size_t index3D(uint_t const n[3], int_t const i[3]) { return size_t(i[2]*n[1] + i[1])*n[0] + i[0]; }

  inline status_t construct_Green_function(
        int const ng[3]
      , double const hg[3]
      , std::vector<double> const & Veff // [ng[2]*ng[1]*ng[0]]
      , std::vector<double> const & xyzZinso // [natoms*8]
      , std::vector<std::vector<double>> const & atom_mat // atomic hamiltonian and overlap matrix
      , int const natoms=0
      , int const echo=0
  ) {
      int constexpr X=0, Y=1, Z=2;
      if (echo > 0) printf("\n#\n# %s(%i, %i, %i)\n#\n\n", __func__, ng[X], ng[Y], ng[Z]);

      int32_t n_original_Veff_blocks[3] = {0, 0, 0};
      for(int d = 0; d < 3; ++d) {
          assert((0 == (ng[d] & 0x3)) && "All grid dimensions must be a multiple of 4!");
          n_original_Veff_blocks[d] = (ng[d] >> 2); // divided by 4
      } // d

      if (echo > 3) { printf("# n_original_Veff_blocks "); printf_vector(" %d", n_original_Veff_blocks, 3); }

      assert(Veff.size() == ng[Z]*ng[Y]*ng[X]);

//    view3D<block4<double>> Vtot_blocks(n_original_Veff_blocks[2], n_original_Veff_blocks[1], n_original_Veff_blocks[0]);
      
      view4D<double> Vtot(n_original_Veff_blocks[2], n_original_Veff_blocks[1], n_original_Veff_blocks[0], 64, 0.0);
      { // scope: reorder Veff into block-structured Vtot
          for(int ibz = 0; ibz < n_original_Veff_blocks[Z]; ++ibz) {
          for(int iby = 0; iby < n_original_Veff_blocks[Y]; ++iby) {
          for(int ibx = 0; ibx < n_original_Veff_blocks[X]; ++ibx) {
              auto const Vtot_xyz = Vtot(ibz,iby,ibx); // get a view1D, i.e. double*
              size_t const of[3] = {ibx*4ul, iby*4ul, ibz*4ul}; // offsets
              for(int i64 = 0; i64 < 64; ++i64) {
                  size_t const coords[3] = {(i64 & 0x3) + of[X],
                                     ((i64 >> 2) & 0x3) + of[Y],
                                             (i64 >> 4) + of[Z]};
                  size_t const izyx = (coords[Z]*ng[Y] + coords[Y])*ng[X] + coords[X];
                  assert(izyx < ng[Z]*ng[Y]*ng[X]);
                  Vtot_xyz[i64] = Veff[izyx];
              } // i64
          }}} // xyz
      } // scope

      // Cartesian cell parameters for the unit cell in which the potential is defined
      double const cell[3] = {ng[0]*hg[0], ng[1]*hg[1], ng[2]*hg[2]};

      // assume periodic boundary conditions and an infinite host crystal,
      // so there is no need to consider k-points
      
      // determine the largest and smallest indices of target blocks
      // given a max distance r_trunc between source blocks and target blocks
      
      double const r_block_circumscribing_sphere = 2*std::sqrt(pow2(hg[X]) + pow2(hg[Y]) + pow2(hg[Z]));

      // assume that the source blocks lie compact in space
      int const nRHSs = n_original_Veff_blocks[2] * n_original_Veff_blocks[1] * n_original_Veff_blocks[0];
      if (echo > 0) printf("# total number of source blocks is %d\n", nRHSs);
      view2D<int16_t> RHS_coords(nRHSs, 4, 0);
      double center_of_mass_RHS[3] = {0, 0, 0};
      double center_of_RHSs[3]     = {0, 0, 0};
      int16_t min_source_coords[3] = {0, 0, 0};
      int16_t max_source_coords[3] = {0, 0, 0};
      double max_distance_from_comass{0};
      double max_distance_from_center{0};
      double farest_from_center{0};
      { // scope:
          int iRHS{0};
          for(int ibz = 0; ibz < n_original_Veff_blocks[Z]; ++ibz) {
          for(int iby = 0; iby < n_original_Veff_blocks[Y]; ++iby) {
          for(int ibx = 0; ibx < n_original_Veff_blocks[X]; ++ibx) {
              RHS_coords(iRHS,X) = ibx;
              RHS_coords(iRHS,Y) = iby;
              RHS_coords(iRHS,Z) = ibz;
              for(int d = 0; d < 3; ++d) {
                  int16_t const rhs_coord = RHS_coords(iRHS,d);
                  center_of_mass_RHS[d] += (rhs_coord*4 + 1.5)*hg[d];
                  min_source_coords[d] = std::min(min_source_coords[d], rhs_coord);
                  max_source_coords[d] = std::max(max_source_coords[d], rhs_coord);
              } // d
              ++iRHS;
          }}} // xyz
          assert(nRHSs == iRHS);
          for(int d = 0; d < 3; ++d) {
              center_of_mass_RHS[d] /= std::max(1, nRHSs);
              center_of_RHSs[d] = ((0.5*(min_source_coords[d] + max_source_coords[d]))*4 + 1.5)*hg[d];
          } // d

          // compute also the largest distance from the center or center of mass
          double max_d2_from_comass{0}, max_d2_from_center{0};
          for(int iRHS = 0; iRHS < nRHSs; ++iRHS) {
              double d2m{0}, d2c{0};
              for(int d = 0; d < 3; ++d) {
                  d2m += pow2((RHS_coords(iRHS,d)*4 + 1.5)*hg[d] - center_of_mass_RHS[d]);
                  auto const dv = (RHS_coords(iRHS,d)*4 + 1.5)*hg[d] - center_of_RHSs[d];
                  d2c += pow2(dv);
                  farest_from_center = std::max(farest_from_center, std::abs(dv));
              } // d
              max_d2_from_center = std::max(max_d2_from_center, d2c);
              max_d2_from_comass = std::max(max_d2_from_comass, d2m);
          } // iRHS
          max_distance_from_center = std::sqrt(max_d2_from_center);
          max_distance_from_comass = std::sqrt(max_d2_from_comass);
      } // scope
      if (echo > 0) printf("# center of mass of RHS blocks is %g %g %g %s\n",
          center_of_mass_RHS[0]*Ang, center_of_mass_RHS[1]*Ang, center_of_mass_RHS[2]*Ang, _Ang);
      if (echo > 0) printf("# center of of RHS blocks is %g %g %g %s\n",
          center_of_RHSs[0]*Ang, center_of_RHSs[1]*Ang, center_of_RHSs[2]*Ang, _Ang);
      if (echo > 0) printf("# largest distance of RHS blocks from center of mass is %g, from center is %g %s\n",
                              max_distance_from_comass*Ang, max_distance_from_center*Ang, _Ang);

      // truncation radius
      double const r_trunc = control::get("green.function.truncation.radius", 10.);
      if (echo > 0) printf("# green.function.truncation.radius=%g %s\n", r_trunc*Ang, _Ang);

      
      
      // count the number of green function elements for each target block
      std::vector<std::vector<uint16_t>> column_indices; // [product_target_blocks];
      size_t nnz{0}; // number of non-zero BRS entries
      uint32_t nRows{0}; // number of rows
      std::vector<uint16_t> ColIndex; // [nnz]
      std::vector<uint32_t> RowStart; // [nRows + 1]

      view2D<int16_t> target_coords; // [nRows][4]
      std::vector<int64_t> target_indices; // [nRows]
      
      
      uint16_t num_target_coords[3] = {0, 0, 0};
      int16_t  min_target_coords[3] = {0, 0, 0};
      int16_t  max_target_coords[3] = {0, 0, 0};
      size_t product_target_blocks{0};
      { // scope:
          auto const radius = std::max(0.0, r_trunc) + r_block_circumscribing_sphere + farest_from_center;
          auto const rtrunc = std::max(0.0, r_trunc) + r_block_circumscribing_sphere;
          int16_t const itr[3] = {int16_t(std::ceil(rtrunc/(4*hg[X]))),
                                  int16_t(std::ceil(rtrunc/(4*hg[Y]))),
                                  int16_t(std::ceil(rtrunc/(4*hg[Z])))};
          auto const r2trunc = pow2(rtrunc);

          for(int d = 0; d < 3; ++d) {
//            double const inv_hg = 1./(4*hg[d]);
              min_target_coords[d] = min_source_coords[d] - itr[d];
              max_target_coords[d] = max_source_coords[d] + itr[d];
              num_target_coords[d] = max_target_coords[d] + 1 - min_target_coords[d];
          } // d
          product_target_blocks = num_target_coords[Z]*num_target_coords[Y]*num_target_coords[X];
          if (echo > 0) printf("# all targets within (%i, %i, %i) and (%i, %i, %i) --> %d x %d x %d = %ld\n",
              min_target_coords[X], min_target_coords[Y], min_target_coords[Z],
              max_target_coords[X], max_target_coords[Y], max_target_coords[Z], 
              num_target_coords[X], num_target_coords[Y], num_target_coords[Z], product_target_blocks);
          column_indices.resize(product_target_blocks);

          assert(nRHSs < 65536 && "the integer type of ColIndex is uint16_t!");
          std::vector<std::vector<bool>> sparsity_pattern(nRHSs);
          for(uint16_t iRHS = 0; iRHS < nRHSs; ++iRHS) {
              sparsity_pattern[iRHS] = std::vector<bool>(product_target_blocks, false);
              int16_t const *const source_coords = RHS_coords[iRHS];
              simple_stats::Stats<> stats[3];
              int inside{0}, outside{0};
              for(int16_t bz = -itr[Z]; bz <= itr[Z]; ++bz) {
              for(int16_t by = -itr[Y]; by <= itr[Y]; ++by) {
              for(int16_t bx = -itr[X]; bx <= itr[X]; ++bx) {
                  int16_t const bxyz[3] = {bx, by, bz}; // difference vector
                  int16_t target_coords[3];
                  for(int d = 0; d < 3; ++d) {
                      target_coords[d] = source_coords[d] + bxyz[d];
                      assert(target_coords[d] >= min_target_coords[d]);
                      assert(target_coords[d] <= max_target_coords[d]);
                  } // d
                  double const d2 = pow2(bx*4*hg[X]) + pow2(by*4*hg[Y]) + pow2(bz*4*hg[Z]);
                  if (d2 < r2trunc) { // inside
                      int16_t idx[3];
                      for(int d = 0; d < 3; ++d) {
                          idx[d] = target_coords[d] - min_target_coords[d];
                          assert(idx[d] >= 0);
                          assert(idx[d] < num_target_coords[d]);
                          stats[d].add(target_coords[d]);
                      } // d
                      auto const idx3 = index3D(num_target_coords, idx);
                      assert(idx3 < product_target_blocks);
                      column_indices[idx3].push_back(iRHS);
                      sparsity_pattern[iRHS][idx3] = true;
                      ++inside;
                  } else { // inside
                      ++outside;
                  } // inside
              }}} // xyz
              if (echo > 7) printf("# RHS at %i %i %i reaches from (%g, %g, %g) to (%g, %g, %g)\n",
                      RHS_coords(iRHS,X), RHS_coords(iRHS,Y), RHS_coords(iRHS,Z),
                      stats[X].min(), stats[Y].min(), stats[Z].min(),
                      stats[X].max(), stats[Y].max(), stats[Z].max());
//               if (echo > 7) printf("# RHS has %d and %d of %d blocks inside and outside, respectively\n",
//                                               inside, outside, inside + outside);
          } // iRHS

          // create a histogram
          std::vector<uint32_t> hist(nRHSs + 1, 0);
          view2D<int16_t> box_target_coords(product_target_blocks, 4);
          for(uint16_t z = 0; z < num_target_coords[Z]; ++z) {
          for(uint16_t y = 0; y < num_target_coords[Y]; ++y) {
          for(uint16_t x = 0; x < num_target_coords[X]; ++x) {
              uint16_t const idx[3] = {x, y, z};
              auto const idx3 = index3D(num_target_coords, idx);
              assert(idx3 < product_target_blocks);
              
              int32_t jdx[3];
              set(jdx, 3, idx);
              add_product(jdx, 3, min_target_coords, 1);
              set(box_target_coords[idx3], 3, jdx);

              int const n = column_indices[idx3].size();
              ++hist[n];
          }}} // xyz
          size_t nall{0};
          for(int n = 0; n <= nRHSs; ++n) {
              nall += hist[n];
              nnz  += hist[n]*n;
          } // n
          if (echo > 0) { printf("# histogram total=%ld  ", nall); printf_vector(" %d", hist.data(), nRHSs + 1); }
          assert(nall == product_target_blocks);
          nRows = nall - hist[0]; // the target block ntries with no RHS do not create a row
          assert(nRows >= 0);
          if (echo > 0) printf("# total number of Green function blocks is %.3f k, average %.1f per source block\n", nnz*.001, nnz/double(nRHSs));
          if (echo > 0) printf("# %d of %d (%.1f %%) target blocks are active\n", nRows, product_target_blocks, nRows/(product_target_blocks*.01));

          assert(nnz < (uint64_t(1) << 32) && "the integer type or RowStart is uint32_t!");

          // resize BSR tables: (Block-compressed Sparse Row format)
          if (echo > 3) { printf("# memory of Green function is %.6f %s (float, twice for double)\n",
                              nnz*2.*64.*64.*sizeof(float)*GByte, _GByte); std::fflush(stdout); }
          ColIndex.resize(nnz);
          RowStart.resize(nRows + 1);
          RowStart[0] = 0;
          target_coords = view2D<int16_t>(nRows, 4, 0);
          target_indices.resize(nRows);
          
          view3D<int32_t> iRow_of_coords(num_target_coords[Z], num_target_coords[Y], num_target_coords[X], -1);

          { // scope: fill BSR tables
              uint32_t iRow{0};
              for(size_t idx3 = 0; idx3 < column_indices.size(); ++idx3) {
                  auto const n = column_indices[idx3].size();
                  if (n > 0) {
                      RowStart[iRow + 1] = RowStart[iRow] + n;
                      // copy the column indices
                      set(&ColIndex[RowStart[iRow]], n, column_indices[idx3].data());
                      // copy the target block coordinates
                      set(target_coords[iRow], 3, box_target_coords[idx3]);
                      target_indices[iRow] = global_coordinates(target_coords[iRow]);

                      int32_t idx[3];
                      set(idx, 3, box_target_coords[idx3]);
                      add_product(idx, 3, min_target_coords, -1);
                      for(int d = 0; d < 3; ++d) { assert(idx[d] >= 0); assert(idx[d] < num_target_coords[d]); }
                      iRow_of_coords(idx[Z], idx[Y], idx[X]) = iRow;

                      // count up the number of active rows
                      ++iRow;
                  } // n > 0
              } // idx3
              assert(nRows == iRow);
              assert(nnz == RowStart[nRows]);
          } // scope
          column_indices.clear(); // not needed beyond this point

          if (echo > 1) { // measure the difference in the number of target blocks of each RHS
              std::vector<uint32_t> nt(nRHSs, 0);
              // traverse BSR structure
              for(uint32_t iRow = 0; iRow < nRows; ++iRow) {
                  for(auto inz = RowStart[iRow]; inz < RowStart[iRow + 1]; ++inz) {
                      auto const iCol = ColIndex[inz];
                      ++nt[iCol];
                  } // inz
              } // iRow
              // analyze nt
              simple_stats::Stats<> st;
              for(uint16_t iRHS = 0; iRHS < nRHSs; ++iRHS) {
                  st.add(nt[iRHS]);
              } // iRHS
              printf("# target blocks per source block: average %.1f +/- %.1f\n", st.avg(), st.var());
          } // echo
          
          // Green function is stored sparse 
          // as std::complex<real_t> green[nnz][64][64] 
          // or real_t green[nnz][2]64][64] for the GPU;

          
          // prepare the finite-difference sequence lists
          for(int dd = 0; dd < 3; ++dd) { // direction of derivative
              int num[3];
              set(num, 3, num_target_coords);
              int const num_dd = num[dd];
              num[dd] = 1;
              if (echo > 0) printf("# FD lists for the %c-direction %d %d %d\n", 'x' + dd, num[X], num[Y], num[Z]);
              simple_stats::Stats<> length_stats;
              std::vector<std::vector<int32_t>> list;
              int const max_lists = nRHSs*num[Z]*num[Y]*num[X];
              list.resize(max_lists);
              int ilist{0};
              for(int iRHS = 0; iRHS < nRHSs; ++iRHS) {
//                if (echo > 0) printf("# FD list for RHS #%i\n", iRHS);
                  auto const & sparsity_RHS = sparsity_pattern[iRHS];
                  for(int iz = 0; iz < num[Z]; ++iz) { //  
                  for(int iy = 0; iy < num[Y]; ++iy) { //   only 2 of these 3 loops have a range > 1
                  for(int ix = 0; ix < num[X]; ++ix) { // 
                      int idx[3] = {ix, iy, iz};
                      for(int id = 0; id < num_dd; ++id) {
                          idx[dd] = id;
//                           if (echo > 0) printf("# FD list for RHS #%i test coordinates %i %i %i\n",
//                                                   iRHS, idx[X], idx[Y], idx[Z]);
                          auto const idx3 = index3D(num_target_coords, idx);
                          if (sparsity_RHS[idx3]) {
                              auto const iRow = iRow_of_coords(idx[Z], idx[Y], idx[X]);
                              assert(iRow >= 0);

                              int inz_found{-1};
                              for(int inz = RowStart[iRow]; inz < RowStart[iRow + 1]; ++inz) {
                                  if (ColIndex[inz] == iRHS) {
                                      inz_found = inz; // store where it was found
                                      inz = RowStart[iRow + 1]; // stop search loop
                                  } // found
                              } // search
                              assert(inz_found >= 0);

                              assert(ilist < max_lists);
                              list[ilist].push_back(inz_found);
                          } // sparsity pattern
                      } // id
                      int const list_length = list[ilist].size();
                      if (list_length > 0) {
                          length_stats.add(list_length);
//                           if (echo > 0) printf("# FD list of length %d for the %c-direction %i %i %i\n",
//                                                   list_length, 'x' + dd, idx[X], idx[Y], idx[Z]);
                          // add end markers
                          list[ilist].push_back(-1);
                          list[ilist].push_back(-1);
                          
                          ++ilist; // create a new list index
                      } // list_length > 0
                  }}} // ixyz
              } // iRHS
              int const number_of_lists = ilist;
              if (echo > 0) printf("# %d FD lists for the %c-direction (%.2f %%), length %.3f +/- %.3f, min %g max %g\n",
                                    number_of_lists, 'x' + dd, number_of_lists/(max_lists*.01),
                                    length_stats.avg(), length_stats.var(), length_stats.min(), length_stats.max());
          } // dd
          
      } // scope
      
      
      
      
      
      // compute which atoms will contribute
      double max_projection_radius{0};
      for(int ia = 0; ia < natoms; ++ia) {
          auto const sigma = xyzZinso[ia*8 + 6];
          auto const projection_radius = 9*std::max(0.0, sigma);
          max_projection_radius = std::max(max_projection_radius, projection_radius);
      } // ia
      if (echo > 3) printf("# largest projection radius is %g %s\n", max_projection_radius*Ang, _Ang);
      
      view2D<double> atom_coords;
      uint32_t nai{0}; // corrected number of all relevant atoms
      uint32_t napc{0}; // number of all projection coefficients
      std::vector<uint32_t> apc_start; // [nai + 1]
      { // scope:
          auto const radius = r_trunc + max_distance_from_center + 2*max_projection_radius;
          int const iimage[3] = {int(std::ceil(radius/cell[X])),
                                 int(std::ceil(radius/cell[Y])),
                                 int(std::ceil(radius/cell[Z]))};
          auto const nimages = size_t(iimage[Z]*2 + 1)*size_t(iimage[Y]*2 + 1)*size_t(iimage[X]*2 + 1);
          auto const natom_images = natoms*nimages;
          apc_start.resize(natom_images + 1); // probably larger than needed
          apc_start[0] = 0;
          atom_coords = view2D<double>(natom_images, 4);
          size_t iai{0};
          for(int z = -iimage[Z]; z <= iimage[Z]; ++z) {
          for(int y = -iimage[Y]; y <= iimage[Y]; ++y) {
          for(int x = -iimage[X]; x <= iimage[X]; ++x) {
//            if (echo > 3) printf("# periodic shifts  %d %d %d\n", x, y, z);
              int const xyz[3] = {x, y, z};
              for(int ia = 0; ia < natoms; ++ia) { // loop over atoms in the unit cell
                  // suggest a new atomic images position
                  double pos[3];
                  for(int d = 0; d < 3; ++d) {
                      pos[d] = xyzZinso[ia*8 + d] + xyz[d]*cell[d];
                  } // d
                  auto const atom_id = int32_t(xyzZinso[ia*8 + 4]); 
                  auto const numax =       int(xyzZinso[ia*8 + 5]);
                  auto const sigma =           xyzZinso[ia*8 + 6];
//                   if (echo > 5) printf("# image of atom #%i at %g %g %g %s\n", atom_id, pos[X]*Ang, pos[Y]*Ang, pos[Z]*Ang, _Ang);

                  double const r_projection = pow2(9*sigma); // atom-dependent
                  double const r2projection = pow2(r_projection);
                  double const r2projection_plus_block = pow2(r_projection + r_block_circumscribing_sphere);
                  
                  // check all target blocks if they are inside the projection radius
                  uint32_t ntb{0}; // number of target blocks
                  for(uint32_t iRow = 0; (iRow < nRows) && (0 == ntb); ++iRow) {
                      int16_t const *const target_block = target_coords[iRow];
                      double d2{0};
                      for(int d = 0; d < 3; ++d) {
                          double const center_of_block = (target_block[d]*4 + 1.5)*hg[d];
                          d2 += pow2(center_of_block - pos[d]);
                      } // d
                      if (d2 < r2projection_plus_block) {
                          // do more precise checking
//                           if (echo > 9) printf("# target block #%i at %i %i %i gets corner check\n",
//                                           iRow, target_block[X], target_block[Y], target_block[Z]);
                          bool any_corner{false};
                          for(int iz = 0; iz < 4; iz += 3) {
                          for(int iy = 0; iy < 4; iy += 3) {
                          for(int ix = 0; ix < 4; ix += 3) {
                              int const ixyz[3] = {ix, iy, iz};
                              double d2i{0};
                              for(int d = 0; d < 3; ++d) {
                                  double const grid_point = (target_block[d]*4 + ixyz[d])*hg[d];
                                  d2i += pow2(grid_point - pos[d]);
                              } // d
                              if (d2i < r2projection) {
                                  any_corner = true; // at least one corner of the block 
                                  // ... is inside the projection radius of this atom
                              } // inside the projection radius
                          }}} // ixyz
                          if (any_corner) {
                              // atom image contributes
                              ++ntb; // stop outer loop
//                               if (echo > 7) printf("# target block #%i at %i %i %i is inside\n",
//                                       iRow, target_block[X], target_block[Y], target_block[Z]);
                          } else { // any_corner
//                               if (echo > 9) printf("# target block #%i at %i %i %i is outside\n",
//                                       iRow, target_block[X], target_block[Y], target_block[Z]);
                          } // any_corner
                      } else { // d2
//                           if (echo > 21) printf("# target block #%i at %i %i %i is far outside\n",
//                                       iRow, target_block[X], target_block[Y], target_block[Z]);
                      } // d2
                  } // iRow

                  if (ntb > 0) {
                      // atom image contributes
                      // at least one target block has an intersection with the projection sphere of this atom image
                      for(int d = 0; d < 3; ++d) {
                          atom_coords(iai,d) = pos[d];
                      } // d
                      atom_coords(iai,3) = atom_id; // global atom index
//                       if (echo > 5) printf("# image of atom #%i at %g %g %g %s contributes to %d target blocks\n",
//                                                 atom_id, pos[X]*Ang, pos[Y]*Ang, pos[Z]*Ang, _Ang, ntb);
                      int const nc = sho_tools::nSHO(numax);
                      apc_start[iai + 1] = apc_start[iai] + nc;
                      ++iai;
                  } else {
//                       if (echo > 15) printf("# image of atom #%i at %g %g %g %s does not contribute\n",
//                                                 atom_id, pos[X]*Ang, pos[Y]*Ang, pos[Z]*Ang, _Ang);
                  } // ntb > 0
              } // ia
          }}} // xyz
          nai = iai; // corrected number of atomic images
          if (echo > 3) printf("# %d of %d (%.2f %%) atom images have an overlap with projection spheres\n",
                                  nai, natom_images, nai/(natom_images*.01));
          napc = apc_start[nai];
          if (echo > 3) printf("# %.3f k atomic projection coefficients, %.2f per atomic image\n", napc*.001, napc/double(nai));
          // projection coefficients for the non-local PAW operations are stored
          // as std::complex<real_t> apc[napc][nRHSs][64]
          // or real_t apc[napc][nRHSs][2][64] for the GPU
          if (echo > 3) printf("# memory of atomic projection coefficients is %.6f %s (float, twice for double)\n",
                                  napc*nRHSs*2.*64.*sizeof(float)*GByte, _GByte);
      } // scope
      
      // Action of the SHO-PAW Hamiltonian H onto a trial Green function G:
      // indicies named after H_{ij} G_{jk} = HG_{ik}
      //
      //      projection:
      //      apc.set(0); // clear
      //      for(jRow < nRows) // reduction
      //        for(RowStart[jRow] <= jnz < Rowstart[jRow + 1]) // reduction
      //          for(jai < nai) // parallel (with data-reuse on Green function element)
      //             for(jprj < nSHO(jai))
      //                for(k64 < 64) // vector parallel
      //                  for(ri < 2) // vector parallel
      //                    for(j64 < 64) // reduction
      //                    {
      //                        apc[apc_start[jai] + jprj][ColIndex[jnz]][ri][k64] +=
      //                        sho_projector[jprj][j64] *
      //                        G[jnz][ri][j64][k64];
      //                    }
      //
      //      meanwhile: kinetic energy including the local potential
      //      HG = finite_difference_kinetic_energy(G);
      //
      //      when projection is done: multiply atomic matrices to projection coefficients
      //      aac.set(0); // clear
      //      for(iai < nai) // parallel
      //        for(iRHS < nRHSs) // parallel
      //          for(iprj < nSHO(iai)) // parallel
      //            for(jprj < nSHO(iai)) // reduction
      //              for(kprj < nSHO(iai)) // parallel
      //                for(k64 < 64) // vector parallel
      //                {
      //                    aac[apc_start[iai] + iprj][iRHS][complex][k64] +=
      //                    atom_matrix[ia[iai]][iprj][jprj][complex] *
      //                    apc[apc_start[iai] + jprj][iRHS][complex][k64];
      //                 }
      //                 // with indirection list ia[iai]
      //
      //      after both, addition:
      //      for(iRow < nRows) // parallel
      //        for(RowStart[iRow] <= inz < Rowstart[iRow + 1]) // parallel
      //          for(iai < nai) // reduction
      //            for(iprj < nSHO(iai)) // reduction
      //              for(i64 < 64)
      //                for(ri < 2) // vector parallel
      //                  for(k64 < 64) // vector parallel
      //                  {
      //                      HG[inz][ri][i64][k64] += 
      //                      aac[apc_start[iai] + iprj][ColIndex[inz]][ri][k64] *
      //                      sho_projector[iprj][i64];
      //                  }
      //
      //      and the finite_difference_kinetic_energy reads:
      //      for(d < 3) // serial
      //        for(istick < nstick[d]) // parallel
      //          while(jnz >= 0, jnz from index_list[d][istick])
      //            for(j4 < 4) // serial (data-reuse)
      //              for(-8 <= imj <= 8) // reduction
      //                for(j16 < 16) // parallel
      //                for(ri < 2) // vector parallel
      //                  for(k64 < 64) // vector parallel
      //                  {
      //                      HG[inz][ri][i4,j16][k64] +=
      //                       G[jnz][ri][j4,j16][k64] * T_FD[|imj|]
      //                  }
      //

      return 0;
  } // construct_Green_function
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_global_coordinates(int const echo=0) {
      status_t stat(0);
      int const n_tested = (1 << 12);
      for(int i = 0; i < n_tested; ++i) {
          int32_t ixyz[3]; // input coordinates
          for(int d = 0; d < 3; ++d) {
              ixyz[d] = simple_math::random(0, 0x1fffff); // [0, 2^21 - 1]
          } // d
          int64_t const i63 = global_coordinates(ixyz);
          uint32_t oxyz[3]; // output coordinates
          status_t is = global_coordinates(oxyz, i63);
          // the output format for i63 indices is octal with 22 digits.
          // exactly 21 octal digits represent the bit-interleaved 3D coordinates
          // and one leading 0 indices the octal notation for non-negative i63
          if (echo > 11) printf("# global_coordinates(%i, %i, %i)\t--> %22.22lo --> (%i, %i, %i)\n",
                                  ixyz[0], ixyz[1], ixyz[2], i63, oxyz[0], oxyz[1], oxyz[2]);
          for(int d = 0; d < 3; ++d) {
              is += (ixyz[d] != oxyz[d]);
          } // d
          if (is != 0) {
              if (echo > 1) printf("# global_coordinates(%i, %i, %i)\t--> %22.22lo --> (%i, %i, %i) failed!\n",
                                  ixyz[0], ixyz[1], ixyz[2], i63, oxyz[0], oxyz[1], oxyz[2]);
              ++stat;
          } // is
      } // i
      if (echo > 3) printf("# %s tested %.3f k random coordinate tuples, found %d errors\n", 
                              __func__, n_tested*.001, int(stat));
      if (echo > 9) {
          int64_t const i63 = -1; // show how invalid i63 indices are displayed
          printf("# global_coordinates(impossible coordinates)\t--> %22.22lo == %d\n", i63, i63);
      } // echo
      return stat;
  } // test_global_coordinates

  inline status_t test_Green_function(int const echo=0) {
      status_t stat(0);
      stat = STATUS_TEST_NOT_INCLUDED; // ToDo: not implemented
      return stat;
  } // test_Green_function

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_global_coordinates(echo);
//    stat += test_Green_function(echo); // deactivated for now
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_function
