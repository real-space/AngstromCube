#pragma once

#include "quantum_numbers.h" // enn_QN_t, ell_QN_t

typedef int status_t;

namespace sho_radial {
  
  inline int constexpr nSHO(int const numax) { return ((3 + numax)*(2 + numax)*(1 + numax))/6; } // number of SHO eigenstates
  inline int constexpr nSHO_radial(int const numax) { return (numax*(numax + 4) + 4)/4; } // number of different radial eigenstates
  
  status_t all_tests();
  
} // namespace sho_radial
