#pragma once
// This file is part of AngstromCube under MIT License

// #doc
// This is the C++ interface to the main routine of the libliveatom.so shared library.
//

#include <cstdint> // int32_t

#include "status.hxx" // status_t

namespace single_atom {

  status_t atom_update(
        char const *what    // decide, what to import or export, only 1st char relevant (w), unless what="#...", then also the 2nd char
      , int32_t natoms      // number of atoms, most relevant when (i) ("initialize"), this is the C++-interface so here call-by-value
      , double  *dp=nullptr // (i)in atomic number, (s)out sigma_cmp, (l)in decide if _cmp or _pot
      , int32_t *ip=nullptr // (i)in numax, (n)out numax, (l)out lmax_cmp or lmax_pot
      , float   *fp=nullptr // (i)in ionization, (c/v/z)in ar2, (u)opt-in {mix_pot, mix_rho}
      , double *const *dpp=nullptr // (c/v/z)out quantities on r2-grid, (u)in vlm, (q)out qlm, (h)out: aHm, aSm, (a)in aDm
  ); // declaration only

  // set the internal variable "version.atom" to a git key passed as preprocessor macro
  status_t set_version(int const echo=0); // declaration only

  // self-tests
  status_t all_tests(int const echo=0); // declaration only

} // namespace single_atom
