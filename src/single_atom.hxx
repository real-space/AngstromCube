#pragma once

#include <cstdint> // int32_t

#include "radial_grid.h" // radial_grid_t

#include "status.hxx" // status_t

namespace single_atom {

  status_t atom_update(
      char const *what,     // decide, what to import or export, only 1st char relevant (w)
      int natoms,           // number of atoms, most relevant when (i) ("initialize")
      double  *dp=nullptr,  // (i)in atomic number, (s)out sigma_cmp, (l)in decide if _cmp or _pot
      int32_t *ip=nullptr,  // (i)in numax, (n)out numax, (l)out lmax_cmp or lmax_pot
      float   *fp=nullptr,  // (i)in ionization, (c/v/z)in ar2, (u)opt-in {mix_pot, mix_rho}
      double**dpp=nullptr); // (c/v/z)out quantities on r2-grid, (u)in vlm, (q)out qlm, (h)out: aHm, aSm, (a)in aDm

  status_t all_tests(int const echo=0);

} // namespace single_atom

