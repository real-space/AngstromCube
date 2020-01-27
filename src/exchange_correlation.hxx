#pragma once

typedef int status_t;

namespace exchange_correlation {

  template<typename real_t>
  real_t lda_PZ81_kernel(     // returns the XC-energy
      real_t const rho,       // [in] density
      real_t &Vdn,            // [out] potential (down-spin if Vup is given)
      real_t const mag=0,     // [optional in] magnetization
      real_t *Vup=nullptr);   // [optional out] up-spin potential

  status_t all_tests(int const echo=0);

} // namespace exchange_correlation
