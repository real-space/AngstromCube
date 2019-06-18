#pragma once

namespace sho_tools {
  

  inline constexpr int nSHO(int const numax) { return ((1 + numax)*(2 + numax)*(3 + numax))/6; }

  inline int lnm_index(int const ell, int const nrn, int const emm) {
      return 0; // ToDo
  } // lnm_index

  inline int ln_index(int const ell, int const nrn) {
      return 0; // ToDo
  } // ln_index
  
  inline int all_tests() {
//       int constexpr numax = 9;
//       int lnm = 0, ln = 0;
//       for(int ell = 0; ell <= numax; ++ell) {
//           for(int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
//               assert(ln_index(ell, nrn) == ln);
//               ++ln;
//               for(int emm = -ell; emm <= ell; ++emm) {
//                   assert(lnm_index(ell, nrn, emm) == lnm);
//                   ++lnm;
//               } // emm
//           } // nrn
//       } // ell
      return 0;
  }
  
} // namespace sho_tools
