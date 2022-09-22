#pragma once

#include "status.hxx" // status_t

namespace exchange_correlation {

  template <typename real_t>
  real_t lda_PZ81_kernel(     // returns the XC-energy
      real_t const rho,       // [in] density
      real_t &Vdn,            // [out] potential (down-spin if Vup is given)
      real_t const mag=0,     // [optional in] magnetization
      real_t *Vup=nullptr);   // [optional out] up-spin potential

  template <typename real_t>
  real_t lda_PW91_kernel(     // returns the XC-energy
      real_t const rho,       // [in] density
      real_t &Vdn,            // [out] potential (down-spin if Vup is given)
      real_t const mag=0,     // [optional in] magnetization
      real_t *Vup=nullptr);   // [optional out] up-spin potential


  // define here which one is the default LDA flavor
  char const default_LDA[] = "PW";

  template <typename real_t> inline
  real_t LDA_kernel(          // returns the XC-energy
      real_t const rho,       // [in] density
      real_t &Vdn,            // [out] potential (down-spin if Vup is given)
      real_t const mag=0,     // [optional in] magnetization
      real_t *Vup=nullptr
  ) {
      return ('W' == default_LDA[1]) ? 
          lda_PW91_kernel(rho, Vdn, mag, Vup): // "Perdew-Wang 1991"
          lda_PZ81_kernel(rho, Vdn, mag, Vup); // "Perdew-Zunger 1981"
  } // LDA_kernel

  status_t all_tests(int const echo=0); // declaration only

} // namespace exchange_correlation
