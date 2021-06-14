#pragma once

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cmath> // std::sqrt
#include <algorithm> // std::max
#include <utility> // std::swap //, std::move
#include <vector> // std::vector<T>

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "green_memory.hxx" // get_memory, free_memory
#include "data_view.hxx" // view3D<T>
#include "inline_math.hxx" // set
#include "simple_stats.hxx" // ::Stats<>

  // ToDo: move to green_utils.hxx or similar
  template <typename uint_t, typename int_t> inline
  size_t index3D(uint_t const n[3], int_t const i[3]) {
      // usual 3D indexing
      return size_t(i[2]*n[1] + i[1])*n[0] + i[0];
  } // index3D


namespace green_kinetic {

  class finite_difference_plan_t {
  private:
      uint32_t *prefix; // in managed memory
      int32_t *fd_list; // in managed memory
      uint32_t n_lists;

  public:
      finite_difference_plan_t() : prefix(nullptr), fd_list(nullptr), n_lists(0) {}

      finite_difference_plan_t(
            int const dd // direction of derivative
          , uint16_t const num_target_coords[3]
          , uint32_t const RowStart[]
          , uint16_t const ColIndex[]
          , view3D<int32_t> const & iRow_of_coords // (Z,Y,X)
          , std::vector<bool> const sparsity_pattern[]
          , unsigned const nRHSs=1
          , int const echo=0
      ) {
          int constexpr X=0, Y=1, Z=2;
          // prepare the finite-difference sequence lists
          char const direction = 'x' + dd;
          assert(X == dd || Y == dd || Z == dd); 
          int num[3];
          set(num, 3, num_target_coords);
          int const num_dd = num[dd];
          num[dd] = 1; // replace number of target blocks in derivative direction
          if (echo > 0) std::printf("# FD lists for the %c-direction %d %d %d\n", direction, num[X], num[Y], num[Z]);
          simple_stats::Stats<> length_stats;
          std::vector<std::vector<int32_t>> list;
          size_t const max_lists = nRHSs*size_t(num[Z])*size_t(num[Y])*size_t(num[X]);
          list.resize(max_lists);
          size_t ilist{0};
          for (int iRHS = 0; iRHS < nRHSs; ++iRHS) {
//                if (echo > 0) std::printf("# FD list for RHS #%i\n", iRHS);
              auto const & sparsity_RHS = sparsity_pattern[iRHS];
              for (int iz = 0; iz < num[Z]; ++iz) { //  
              for (int iy = 0; iy < num[Y]; ++iy) { //   only 2 of these 3 loops have a range > 1
              for (int ix = 0; ix < num[X]; ++ix) { // 
                  int idx[3] = {ix, iy, iz};
                  for (int id = 0; id < num_dd; ++id) { // loop over direction to derive
                      idx[dd] = id; // replace index in the derivate direction
//                           if (echo > 0) std::printf("# FD list for RHS #%i test coordinates %i %i %i\n",
//                                                   iRHS, idx[X], idx[Y], idx[Z]);
                      auto const idx3 = index3D(num_target_coords, idx);
                      if (sparsity_RHS[idx3]) {
                          auto const iRow = iRow_of_coords(idx[Z], idx[Y], idx[X]);
                          assert(iRow >= 0);

                          int32_t inz_found{-1};
                          for (auto inz = RowStart[iRow]; inz < RowStart[iRow + 1]; ++inz) {
                              if (ColIndex[inz] == iRHS) {
                                  inz_found = inz; // store where it was found
                                  inz = RowStart[iRow + 1]; // stop search loop
                              } // found
                          } // search
                          assert(inz_found >= 0); // fails at inconsistency between sparsity_pattern and the BSR tables

                          assert(ilist < max_lists);
                          list[ilist].push_back(inz_found);
                      } // sparsity pattern
                  } // id
                  int const list_length = list[ilist].size();
                  if (list_length > 0) {
                      length_stats.add(list_length);
//                           if (echo > 0) std::printf("# FD list of length %d for the %c-direction %i %i %i\n",
//                                                   list_length, direction, idx[X], idx[Y], idx[Z]);
                      // add end-of-sequence markers
                      list[ilist].push_back(-1);
                      
                      ++ilist; // create a new list index
                  } // list_length > 0
              }}} // ixyz
          } // iRHS
          n_lists = ilist; assert(n_lists == ilist && "too many lists, max. 2^32-1");
          if (echo > 0) std::printf("# %d FD lists for the %c-direction (%.2f %%), length %.3f +/- %.3f, min %g max %g\n",
                                n_lists, direction, n_lists/(max_lists*.01),
                                length_stats.avg(), length_stats.var(), length_stats.min(), length_stats.max());

          // store index starts in managed memory
          prefix = get_memory<uint32_t>(n_lists + 1); // create in GPU memory
          prefix[0] = 0;
          for (int ilist = 0; ilist < n_lists; ++ilist) {
              int const n = list[ilist].size();
              prefix[ilist + 1] = prefix[ilist] + n;
          } // ilist
          size_t const ntotal = prefix[n_lists];
          if (echo > 0) std::printf("# FD lists for the %c-direction require %ld uint32_t, i.e. %.3f kByte\n",
                                  direction, ntotal, ntotal*sizeof(uint32_t)*1e-3);
          // store indices in managed memory
          fd_list = get_memory<int32_t>(ntotal); // create in GPU memory
          { // scope: copy indices into managed memory
              size_t ntotal_check{0};
              for (int ilist = 0; ilist < n_lists; ++ilist) {
                  int const n = list[ilist].size();
                  assert(ntotal_check == prefix[ilist]); // sanity
                  ntotal_check += n;
                  set(&fd_list[prefix[ilist]], n, list[ilist].data()); // copy into GPU memory
              } // ilist
              assert(ntotal == ntotal_check); // sanity
          } // scope
 
      } // constructor

      finite_difference_plan_t& operator= (finite_difference_plan_t && rhs) {
          // std::printf("# finite_difference_plan_t& operator= (finite_difference_plan_t && rhs);\n");
          std::swap(fd_list, rhs.fd_list);
          std::swap(prefix , rhs.prefix);
          std::swap(n_lists, rhs.n_lists);
          return *this;
      } // move assignment

      finite_difference_plan_t(finite_difference_plan_t && rhs) = delete;
//       {   std::printf("# finite_difference_plan_t(finite_difference_plan_t && rhs);\n");
//           *this = std::move(rhs);
//       } // move constructor

      finite_difference_plan_t(finite_difference_plan_t const & rhs) = delete; // copy constructor

      finite_difference_plan_t& operator= (finite_difference_plan_t const & rhs) = delete; // move assignment

      ~finite_difference_plan_t() {
          std::printf("# destruct %s, pointers= %p and %p\n", __func__, (void*)fd_list, (void*)prefix); std::fflush(stdout);
          free_memory(fd_list);
          free_memory(prefix);
      } // destructor

      uint32_t size() const { return n_lists; } // number of lists
      int32_t const * list(uint32_t const i) const { return fd_list + prefix[i]; }

  }; // class finite_difference_plan_t

  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_finite_difference(int const echo=0) {
      return STATUS_TEST_NOT_INCLUDED;
  } // test_finite_difference

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_finite_difference(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_kinetic
