#pragma once

#include "status.hxx" // status_t
#include "complex_tools.hxx" // complex_name, is_complex, conjuagte, to_complex_t
#include "complex_tools.hxx" // Lorentzian
#include "data_view.hxx" // view2D<T>

namespace dense_solver {
  
  template<typename complex_t, typename real_t>
  inline status_t solve(view3D<complex_t> & SHm
          , char const *x_axis
          , int const echo=0) { // log-level

      int constexpr S=0, H=1; // static indices for S:overlap matrix, H:Hamiltonian matrix

      status_t stat(0);
      
      int const nB  = SHm.dim1();
      int const nBa = SHm.stride();
      assert(nBa >= nB); // stride may not be smaller than the dimension of the two square matrices

      std::vector<real_t> eigvals(nB, 0.0);
      auto const ovl_eig = int(control::get("dense_solver.test.overlap.eigvals", 0.));
      char const hermitian = *control::get("dense_solver.test.hermitian", "both") | 32; // 'n':none, 's':overlap, 'h':Hamiltonian, 'b':both
      real_t const E_imag = control::get("electronic.temperature", 9.765625e-4);
      double const f_dos = -1./constants::pi;
      double const energy_range[2] = {-1, 1}; // ToDo: external input
      
      auto const nE_Green = int(control::get("dense_solver.test.green.function", 0.));
      if (nE_Green > 0) {
          auto const nE = nE_Green;
          complex_t const minus1(-1);
          if (!is_complex(minus1)) {
              warn("# Green functions can only be computed in complex versions"); return stat;
          } // is not complex

          double const dE = (energy_range[1] - energy_range[0])/nE;
          if (echo > 0) printf("\n## E_real (%s), DoS:\n", _eV);
          view2D<complex_t> ESmH(nB, nBa, complex_t(0)); // get memory
          for(int iE = 0; iE <= nE; ++iE) {
              real_t const E_real = iE*dE + energy_range[0];
              auto const E = to_complex_t<complex_t,real_t>(std::complex<real_t>(E_real, E_imag));
              if (echo > 99) printf("# Green function for energy point (%g %s, %g %s)\n",
                                        std::real(E)*eV,_eV, std::imag(E)*Kelvin,_Kelvin);
              int constexpr check = 1;
              view2D<complex_t> ESmH_copy(check*nB, nBa, complex_t(0)); // get memory
              view2D<complex_t> Sinv(nB, nBa, complex_t(0)); // get memory
              // construct matrix to be inverted: E*S - H
              for(int iB = 0; iB < nB; ++iB) {
                  set(Sinv[iB], nB, SHm(S,iB)); // copy S
                  set(ESmH[iB], nB, SHm(H,iB), minus1); // -H
                  add_product(ESmH[iB], nB, SHm(S,iB), E); // +E*S
                  if (check) set(ESmH_copy[iB], nB, ESmH[iB]); // copy
              } // iB

              // Green function G = (E*S - H)^{-1}, mind that with a non-orthogonal basis set, we have to multiply S before interpreting the trace
              auto const stat_inv = linear_algebra::inverse(nB, ESmH.data(), ESmH.stride());
              view2D<complex_t> Green(ESmH.data(), ESmH.stride()); // wrap
              if (int(stat_inv)) warn("inversion failed for %s, E= %g %s", x_axis, E_real*eV,_eV);
              stat += stat_inv;

#ifdef DEVEL
              if (check) {
                  real_t devN(0), devT(0);
                  for(int iB = 0; iB < nB; ++iB) {
                      for(int jB = 0; jB < nB; ++jB) {
                          complex_t cN(0), cT(0);
                          for(int kB = 0; kB < nB; ++kB) {
                              cN += Green(iB,kB) * ESmH_copy(kB,jB); 
                              cT += ESmH_copy(iB,kB) * Green(kB,jB); 
                          } // kB
                          complex_t const diag = (iB == jB);
                          devN = std::max(devN, std::abs(cN - diag));
                          devT = std::max(devT, std::abs(cT - diag));
                          if (echo > 99) printf("# iB=%i jB=%i cN= %g %g, dev= %.1e,\tcT= %g %g, dev= %.1e\n", iB, jB,
                              std::real(cN), std::imag(cN), std::abs(cN - diag), std::real(cT), std::imag(cT), std::abs(cT - diag));
                      } // jB
                  } // iB
                  if ((echo > 19) || ((echo > 0) && (devN + devT > 1e-7))) {
                      printf("# deviation of G * (ES-H) from unity is %.2e and %.2e transposed\n", devN, devT);
                  }
              } // check
#endif

              view2D<complex_t> GreenS(nB, nBa, complex_t(0)); // get memory
              for(int iB = 0; iB < nB; ++iB) {
                  for(int jB = 0; jB < nB; ++jB) {
                      complex_t c(0);
                      for(int kB = 0; kB < nB; ++kB) {
                          c += Green(iB,kB) * SHm(S,kB,jB); 
                      } // kB
                      GreenS(iB,jB) = c;
                  } // jB
              } // iB

              // density of states
              double density{0};
              for(int iB = 0; iB < nB; ++iB) {
                  density += f_dos*std::imag(GreenS(iB,iB));
              } // iB
              if (echo > 0) printf("%g %g\n", std::real(E)*eV, density);
              if (echo > 99) printf("# DoS for energy point (%g %s, %g %s) %g\n",
                  std::real(E)*eV,_eV, std::imag(E)*Kelvin,_Kelvin, density);
            
          } // iE
          if (echo > 0) printf("\n");
      } // green_function
      
      
      for(int s0h1 = 0; s0h1 < 2; ++s0h1) { // loop must run forward and serial
          if (echo > 0) printf("\n");
          auto const matrix_name = s0h1 ? "Hamiltonian" : "overlap";
          real_t const  u = s0h1 ?  eV :  1; // output unit conversion factor 
          auto   const _u = s0h1 ? _eV : ""; // unit symbol

#ifdef DEVEL
          // display S and H
          if (echo > 9 - s0h1) {
              printf("\n# %s matrix (%s) for Bloch phase", matrix_name, _u);
              for(int d = 0; d < 3; ++d) {
//                   printf(" %g+i*%g ", std::real(Bloch_phase[d]), std::imag(Bloch_phase[d]));
              } // d
              printf(":\n");
              for(int iB = 0; iB < nB; ++iB) {
                  printf("# row%3i ", iB);
                  for(int jB = 0; jB < nB; ++jB) {
                      printf("%8.3f", SHm(s0h1,iB,jB)*u);
                  } // jB
                  printf("\n");
              } // iB
              printf("\n");
          } // echo

          // check if the matrix is symmetric/Hermitian
          if (('b' == hermitian) || ((s0h1 ? 'h' : 's') == hermitian)) {
              real_t diag{0}; // imaginary part of diagonal elements
              real_t offr{0}; // off-diagonal elements real part
              real_t offi{0}; // off-diagonal elements imaginary part
              for(int iB = 0; iB < nB; ++iB) {
                  diag = std::max(diag, std::abs(std::imag(SHm(s0h1,iB,iB))));
                  for(int jB = iB + 1; jB < nB; ++jB) { // triangular loop
                      auto const dev = SHm(s0h1,iB,jB) - conjugate(SHm(s0h1,jB,iB));
                      offr = std::max(offr, std::abs(std::real(dev)));
                      offi = std::max(offi, std::abs(std::imag(dev)));
                  } // jB
              } // iB
              if (echo > 1) printf("# %s deviates from hermitian: imag(diag)= %.1e  imag(off)= %.1e  real(off)= %.1e\n",
                                      matrix_name, diag, offi, offr);
          } // check hermiticity
#endif
          
          // diagonalize matrix
          status_t stat_eig(0);
          if (1 == s0h1) {
              stat_eig = linear_algebra::eigenvalues(eigvals.data(), nB, SHm(H,0), nBa, SHm(S,0), nBa);
          } else if (ovl_eig) {
              view2D<complex_t> S_copy(nB, nBa); // get memory
              set(S_copy.data(), nB*nBa, SHm(S,0)); // copy overlap matrix S into work array W
              stat_eig = linear_algebra::eigenvalues(eigvals.data(), nB, S_copy.data(), nBa);
          } // ovl_eig

          // show result
          if ((H == s0h1) || ovl_eig) {
              if (stat_eig) {
                  warn("diagonalizing the %s matrix failed, status= %i", matrix_name, int(stat_eig));
                  stat += stat_eig;
              } else if (nB > 0) {
                  double const lowest_eigenvalue = eigvals[0], highest_eigenvalue = eigvals[nB - 1];
                  if (echo > 2) {
                      printf("%s%s", x_axis, s0h1 ? "" : matrix_name);
                      int constexpr mB = 8; // show at most the 6 lowest + 2 highest eigenvalues
                      for(int iB = 0; iB < std::min(nB - 2, mB - 2); ++iB) {
                          printf(" %g", eigvals[iB]*u);
                      } // iB
                      if (nB > mB) printf(" ..."); // there are more eigenvalues than we display
                      printf(" %g %g %s\n", eigvals[nB - 2]*u, eigvals[nB - 1]*u, _u); // last two
                  } // echo
                  if (echo > 4) printf("# lowest and highest eigenvalue of the %s matrix are %g and %g %s, respectively\n", 
                                            matrix_name, lowest_eigenvalue*u, highest_eigenvalue*u, _u);
                  if (S == s0h1) {
                      if (lowest_eigenvalue <= 0) {
                          warn("overlap matrix has instable eigenvalues, lowest= %g", lowest_eigenvalue);
                      } else if (lowest_eigenvalue < .1) {
                          warn("overlap matrix has critical eigenvalues, lowest= %g", lowest_eigenvalue);
                      } // lowest_eigenvalue
#ifdef DEVEL
                  } else {
                      // experiment: show Lorentzians where the eigenenergies are
                      auto const nE = int(control::get("dense_solver.test.green.lehmann", 0.));
                      if (nE > 0) {
                          double const dE = (energy_range[1] - energy_range[0])/nE;
                          // from the eigenvalues, we can plot the density of states via the Lehmann representation
                          if (echo > 0) printf("\n## E_real (%s), DoS (from Lehmann representation):\n", _eV);
                          for(int iE = 0; iE <= nE; ++iE) {
                              real_t const E_real = iE*dE + energy_range[0];
                              double density{0};
                              for(int iB = 0; iB < nB; ++iB) {
                                  density += f_dos*Lorentzian(E_real - eigvals[iB], E_imag);
                              } // density
                              if (echo > 0) printf("%g %g\n", E_real*eV, density);
                          } // iE
                          if (echo > 0) printf("\n");
                      } // nE > 0
#endif
                  } // S == s0h1
              } // stat_eig
          } // H or ovl_eig

          if (1 == s0h1) {
              if (0 == stat_eig) {
                  // ToDo: export the entire spectrum and the lowest n_occupied_bands eigenvectors
                  //       for this k-point to be later used for the density generation.
                  //       n_occupied_bands should be chosen large so that also the tail 
                  //       of the Fermi-Dirac distribution is captured. We can check that
                  //       by ensuring that the occupation numbers of the next higher state
                  //          f_FD(eigvals[n_occupied_bands] - E_Fermi) < 10^{-15}
                  //       are small.
              } // success
          } // H

      } // s0h1
            
      return stat;
  } // solve
  
  inline status_t all_tests(int const echo=0) { return -1; }

} // namespace dense_solver
