#pragma once

#include "grid_operators.hxx" // ::grid_operator_t

#include "status.hxx" // status_t

namespace davidson_solver {

  template<typename real_t, typename real_fd_t=real_t>
  status_t eigensolve(real_t waves[] // on entry start wave functions, on exit improved eigenfunctions
    , double *const energies // export eigenenergies
    , int const nbands // number of bands
    , grid_operators::grid_operator_t<real_t,real_fd_t> const &op
    , int const echo=0 // log output level
    , float const mbasis=2 // factor enlarging the space of trial functions
    , unsigned const niterations=9 // number of Davidson iterations
  );

  template<typename real_t, typename real_fd_t=real_t> inline
  status_t rotate(real_t waves[] // on entry start wave functions, on exit improved eigenfunctions
    , double *const energies // export eigenenergies
    , int const nbands // number of bands
    , grid_operators::grid_operator_t<real_t,real_fd_t> const &op
    , int const echo=0)
      { return eigensolve(waves, energies, nbands, op, echo, 1, 1); }

  status_t all_tests(int const echo=0);

} // namespace davidson_solver
