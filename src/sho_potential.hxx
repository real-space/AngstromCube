#pragma once

#include <vector> // std::vector<T>
#include <fstream> // std::ifstream

#include "status.hxx" // status_t
#include "inline_math.hxx" // set
#include "recorded_warnings.hxx" // warn

namespace sho_potential {

  template <typename real_t>
  status_t load_local_potential(std::vector<real_t> & vtot, int dims[3], char const *filename, int const echo=0) {
      status_t stat(0);
      set(dims, 3, 0); // clear
      vtot.clear();
      { // scope: read in the potential from a file
          std::ifstream infile(filename);
          size_t npt = 0; 
          if (!infile.is_open()) {
              if (echo > 1) printf("# %s failed to open file '%s'\n",  __func__, filename);
              return 1; // failure
          }
          for(int d = 2; d >= 0; --d) {
              char sep;
              infile >> sep >> dims[d];
              if (echo > 3) printf("# found dim %c with %i grid points\n", 120+d, dims[d]);
          } // d
          size_t const all = dims[0]*dims[1]*dims[2];
          vtot.reserve(all);
          size_t idx;
          double val;
          while (infile >> idx >> val) {
              assert(idx < all);
              ++npt;
              vtot.push_back(val);
          } // while
          if (echo > 3) printf("# %s use %i values from file '%s'\n", __func__, npt, filename);
          if (npt < all) {
              warn("when loading local potential from file '%s' found only %ld entries while expected %ld", filename, vtot.size(), all);
              ++stat;
          } // not enough entries
      } // scope
      return stat;
  } // load_local_potential
  
  status_t all_tests(int const echo=0);
  
} // namespace sho_potential
