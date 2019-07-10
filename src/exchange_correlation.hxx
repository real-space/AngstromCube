#pragma once

typedef int status_t;

namespace exchange_correlation {

  template<typename real_t>
  real_t lda_PZ81_kernel(     // returns the XC-energy
      real_t const rho,       // [in] density
      real_t &Vup,            // [out] potential (if Vdn is given)
      real_t const mag=0,     // [optional in] magnetization
      real_t *Vdn=nullptr);   // [optional out] down-spin potential

  status_t all_tests();

} // namespace exchange_correlation
