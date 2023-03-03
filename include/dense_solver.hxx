#pragma once

#include <cstdio> // std::printf

#include "status.hxx" // status_t
#include "data_view.hxx" // view3D<T>
#include "inline_math.hxx" // pow2
#include "recorded_warnings.hxx" // warn

namespace dense_solver {

  template <typename real_t>
  inline real_t Lorentzian(real_t const re, real_t const im) { return -im/(pow2(im) + pow2(re)); }

  template <typename real_t>
  inline void display_spectrum(real_t const eigvals[], int const nB, char const *x_axis
      , double const u=1, char const *_u="", char const *matrix_name="", int const mB=32) {
      if (nB < 2) return;
      std::printf("%s%s", x_axis, matrix_name);
      // show at most the (mB - 2) lowest + 2 highest eigenvalues
      for (int iB = 0; iB < std::min(nB - 2, mB - 2); ++iB) {
          std::printf(" %g", eigvals[iB]*u);
      } // iB
      if (nB > mB) std::printf(" ..."); // there are more eigenvalues than we display
      std::printf(" %g %g %s\n", eigvals[nB - 2]*u, eigvals[nB - 1]*u, _u); // last two
  } // display_spectrum

  template <typename complex_t>
  status_t solve(
        view3D<complex_t> & HSm // Hamiltonian and Overlap, both[nB][stride]
      , char const *x_axis
      , int const echo=0 // log-level
      , int const nbands=0 // number of bands
      , double *eigenenergies=nullptr // export nbands eigenvalues
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace dense_solver
