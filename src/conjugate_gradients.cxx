#include <cstdio> // printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <algorithm> // std::swap<T>
#include <cmath> // std::sin

#include "conjugate_gradients.hxx" // ::eigensolve

#include "real_space.hxx" // ::grid_t
#include "control.hxx" // ::get
#include "complex_tools.hxx" // complex_name
#include "simple_math.hxx" // random
#include "grid_operators.hxx" // ::grid_operator_t

namespace conjugate_gradients {
  // solve iteratively for the lowest eigenstates of an implicitly given Hamiltonian using the conjugate gradients method

#ifndef NO_UNIT_TESTS
  
  template <typename complex_t>
  status_t test_solver(int const echo=9) {
      int const nbands = std::min(8, int(control::get("conjugate_gradients.test.num.bands", 4)));
      if (echo > 3) printf("\n# %s %s<%s> with %d bands\n", __FILE__, __func__, complex_name<complex_t>(), nbands);
      status_t stat{0};
      int const dims[] = {8, 8, 8};
      // particle in a box: lowest mode: sin(xyz*pi/L)^3 --> k_x=k_y=k_z=pi/L
      // --> ground state energy = 3*(pi/9)**2 Rydberg = 0.182 Hartree
      //                           3*(pi/8.78)**2      = 0.192 (found)
      //                           3*(pi/8)**2         = 0.231
      // first excitation energies should be 2*(pi/9)**2 + (2*pi/9)**2 = 0.384 Hartree (3-fold degenerate)
      real_space::grid_t g(dims);
      std::vector<complex_t> psi(nbands*g.all(), 0.0);

      char const swm = *control::get("conjugate_gradients.test.start.waves", "a"); // 'a':good, 'r':random
      if ('a' == swm) { // scope: create good start wave functions
          double const k = constants::pi/8.78; // ground state wave vector
          double wxyz[8] = {1, 0,0,0, 0,0,0, 0};
          for(int iz = 0; iz < dims[2]; ++iz) { wxyz[3] = iz - 3.5; double const cos_z = std::cos(k*wxyz[3]);
          for(int iy = 0; iy < dims[1]; ++iy) { wxyz[2] = iy - 3.5; double const cos_y = std::cos(k*wxyz[2]);
          for(int ix = 0; ix < dims[0]; ++ix) { wxyz[1] = ix - 3.5; double const cos_x = std::cos(k*wxyz[1]);
              if (nbands > 4) {
                  wxyz[4] = wxyz[1]*wxyz[2]; // x*y (ell=2)
                  wxyz[5] = wxyz[2]*wxyz[3]; // y*z (ell=2)
                  wxyz[6] = wxyz[3]*wxyz[1]; // z*x (ell=2)
                  wxyz[7] = wxyz[1]*wxyz[2]*wxyz[3]; // x*y*z (ell=3)
              } // nbands > 4
              int const ixyz = (iz*dims[1] + iy)*dims[0] + ix;
              for(int iband = 0; iband < nbands; ++iband) {
                  psi[iband*g.all() + ixyz] = wxyz[iband]*cos_x*cos_y*cos_z; // good start wave functions
              } // iband
          }}} // ix iy iz
          if (echo > 2) printf("\n# %s: use cosine solutions as start vectors\n", __func__);
      } else if ('r' == swm) {
          for(size_t i = 0; i < nbands*g.all(); ++i) {
              psi[i] = simple_math::random(-1., 1.); // random wave functions (most probably not very good)
          } // i
          if (echo > 2) printf("\n# %s: use random values as start vectors\n", __func__);
      } else {
          for(int iband = 0; iband < nbands; ++iband) {
              psi[iband*g.all() + iband] = 1; // bad start wave functions
          } // iband
          if (echo > 2) printf("\n# %s: use as start vectors some delta functions at the boundary\n", __func__);
      } // swm (start wave method)

      // construct Hamiltonian and Overlap operator
      using real_t = decltype(std::real(complex_t(1)));
      grid_operators::grid_operator_t<complex_t,real_t> const op(g); 

      std::vector<double> eigenvalues(nbands, 0.0);
      int const nit = control::get("conjugate_gradients.test.max.iterations", 1.);
      for(int it = 0; it < nit; ++it) {
          stat += eigensolve(psi.data(), eigenvalues.data(), nbands, op, echo);
      } // it
      return stat;
  } // test_solver

  status_t all_tests(int const echo) {
      status_t stat{0};
      stat += test_solver<std::complex<double>>(echo);
      stat += test_solver<std::complex<float>> (echo);
      stat += test_solver<double>(echo);
      stat += test_solver<float> (echo); // test complation and convergence
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace conjugate_gradients
