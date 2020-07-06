#pragma once

#include "status.hxx" // status_t
#include "grid_operators.hxx" // ::grid_operator_t

namespace conjugate_gradients {

  template<typename real_t>
  status_t eigensolve(real_t eigenstates[] // on entry start wave functions, on exit improved eigenfunctions
    , int const nbands // number of bands
    , grid_operators::grid_operator_t<real_t,real_t> const & op // grid operator descriptor
    , int const echo=9 // log output level
    , float const threshold=1e-8f
    , double *const eigenvalues=nullptr); // export results
  
  template<typename real_t> inline double tiny();
  template<> inline double tiny<double>() { return 2.25e-308; }
  template<> inline double tiny<float> () { return 1.18e-38; }

  status_t all_tests(int const echo=0);

} // namespace conjugate_gradients
