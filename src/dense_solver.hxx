#pragma once

#include <cstdio> // printf


#include "status.hxx" // status_t
#include "complex_tools.hxx" // complex_name, is_complex, conjugate, to_complex_t
#include "linear_algebra.hxx" // ::eigenvalues, ::generalized_eigenvalues, ::inverse
#include "constants.hxx" // ::pi
#include "data_view.hxx" // view2D<T>
#include "display_units.h" // eV, _eV, Kelvin, _Kelvin
#include "control.hxx" // ::get
#include "recorded_warnings.hxx" // warn

#ifndef NO_UNIT_TESTS
  #include "simple_math.hxx" // ::random<real_t>
#endif

namespace dense_solver {

  template <typename real_t>
  inline real_t Lorentzian(real_t const re, real_t const im) { return -im/(pow2(im) + pow2(re)); }

  template <typename real_t>
  inline void display_spectrum(real_t const eigvals[], int const nB, char const *x_axis
      , double const u=1, char const *_u="", char const *matrix_name="", int const mB=32) {
      if (nB < 2) return;
      printf("%s%s", x_axis, matrix_name);
      // show at most the (mB - 2) lowest + 2 highest eigenvalues
      for(int iB = 0; iB < std::min(nB - 2, mB - 2); ++iB) {
          printf(" %g", eigvals[iB]*u);
      } // iB
      if (nB > mB) printf(" ..."); // there are more eigenvalues than we display
      printf(" %g %g %s\n", eigvals[nB - 2]*u, eigvals[nB - 1]*u, _u); // last two
  } // display_spectrum
  
  
  template <typename complex_t>
  inline status_t solve(
        view3D<complex_t> & HSm // Hamiltonian and Overlap both[nB,stride]
      , char const *x_axis
      , int const echo=0 // log-level
      , int const nbands=0
      , double *eigenenergies=nullptr // export nbands eigenvalues
  ) {
            
      using real_t = decltype(std::real(complex_t(1))); // base type

      int constexpr H=0, S=1; // static indices for H:Hamiltonian matrix, S:overlap matrix

      status_t stat(0);
      
      int const nB  = HSm.dim1(); // number of basis functions
      int const nBa = HSm.stride();
      assert(nBa >= nB && "stride may not be smaller than the dimension of the two square matrices");

      std::vector<real_t> eigvals(nB, 0.0);
      auto const ovl_eig = int(control::get("dense_solver.test.overlap.eigvals", 0.));
#ifdef DEVEL
      char const hermitian = *control::get("dense_solver.test.hermitian", "none") | 32; // 'n':none, 's':overlap, 'h':Hamiltonian, 'b':both
#endif
      real_t const E_imag = control::get("electronic.temperature", 9.765625e-4);
      double const f_dos = -1./constants::pi;
      double const energy_range[2] = {-1, 1}; // ToDo: external input

      auto const nE = int(control::get("dense_solver.test.green.function", 0.));
      if (nE > 0) {
          if (!is_complex<complex_t>()) {
              warn("# Green functions can only be computed in complex versions", 0); return stat;
          } // is not complex

          double const dE = (energy_range[1] - energy_range[0])/nE;
          if (echo > 0) printf("\n## E_real (%s), DoS:\n", _eV);
          view2D<complex_t> ESmH(nB, nBa, complex_t(0)); // get memory
          for(int iE = 0; iE <= nE; ++iE) {
              real_t const E_real = iE*dE + energy_range[0];
              auto const E = to_complex_t<complex_t,real_t>(std::complex<real_t>(E_real, E_imag));
              if (echo > 9) printf("# Green function for energy point (%g %s, %g %s)\n",
                                        std::real(E)*eV,_eV, std::imag(E)*Kelvin,_Kelvin);
              int constexpr check = 1; // 1:do checks, 0:no checks
              view2D<complex_t> ESmH_copy(check*nB, nBa, complex_t(0)); // get memory
              view2D<complex_t> Sinv(nB, nBa, complex_t(0)); // get memory
              // construct matrix to be inverted: E*S - H
              for(int iB = 0; iB < nB; ++iB) {
                  set(Sinv[iB], nB, HSm(S,iB)); // copy S
                  set(ESmH[iB], nB, HSm(H,iB), real_t(-1)); // -H
                  add_product(ESmH[iB], nB, HSm(S,iB), E); // +E*S
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
                          c += Green(iB,kB) * HSm(S,kB,jB); 
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
      
      
      status_t stat_eig(0);
      for(int h0s1 = 1; h0s1 >= 0; --h0s1) { // loop must run down and serial
          if (echo > 0) printf("\n");
          auto const matrix_name = h0s1 ? "overlap" : "Hamiltonian";
          real_t const  u = h0s1 ?  1 :  eV; // output unit conversion factor 
          auto   const _u = h0s1 ? "" : _eV; // unit symbol
#ifdef DEVEL
          // display H and S
          if (echo > 28 + h0s1) {
              printf("\n# %s matrix (%s)", matrix_name, _u);
//               printf(" for Bloch phase");
//               for(int d = 0; d < 3; ++d) {
//                   printf(" %g+i*%g ", std::real(Bloch_phase[d]), std::imag(Bloch_phase[d]));
//               } // d
              printf(":\n");
              for(int iB = 0; iB < nB; ++iB) {
                  printf("# row%3i ", iB);
                  for(int jB = 0; jB < nB; ++jB) {
                      printf(" %7.3f %g", std::real(HSm(h0s1,iB,jB))*u, std::imag(HSm(h0s1,iB,jB))*u);
                  } // jB
                  printf("\n");
              } // iB
              printf("\n");
          } // echo

          // check if the matrix is symmetric/Hermitian
          if (('b' == hermitian) || ((h0s1 ? 's' : 'h') == hermitian)) {
              real_t diag{0}; // imaginary part of diagonal elements
              real_t offr{0}; // off-diagonal elements real part
              real_t offi{0}; // off-diagonal elements imaginary part
              for(int iB = 0; iB < nB; ++iB) {
                  diag = std::max(diag, std::abs(std::imag(HSm(h0s1,iB,iB))));
                  for(int jB = iB + 1; jB < nB; ++jB) { // triangular loop
                      auto const dev = HSm(h0s1,iB,jB) - conjugate(HSm(h0s1,jB,iB));
                      offr = std::max(offr, std::abs(std::real(dev)));
                      offi = std::max(offi, std::abs(std::imag(dev)));
                  } // jB
              } // iB
              if (echo > 1) printf("# %s deviates from hermitian: imag(diag)= %.1e  imag(off)= %.1e  real(off)= %.1e\n",
                                      matrix_name, diag, offi, offr);
          } // check hermiticity
#endif // DEVEL
          
          // diagonalize matrix
          if (H == h0s1) {
              stat_eig = linear_algebra::eigenvalues(eigvals.data(), nB, HSm(H,0), nBa, HSm(S,0), nBa);
          } else if (ovl_eig) {
              view2D<complex_t> S_copy(nB, nBa); // get memory
              set(S_copy.data(), nB*nBa, HSm(S,0)); // copy overlap matrix S into work array S_copy
              stat_eig = linear_algebra::eigenvalues(eigvals.data(), nB, S_copy.data(), nBa);
          } // ovl_eig

          // show result
          if ((H == h0s1) || ovl_eig) {
              if (stat_eig) {
                  warn("diagonalizing the %s matrix failed, status= %i", matrix_name, int(stat_eig));
                  stat += stat_eig;
              } else if (nB > 0) {
                  double const lowest_eigenvalue = eigvals[0], highest_eigenvalue = eigvals[nB - 1];
                  if (echo > 2) {
                      display_spectrum(eigvals.data(), nB, x_axis, u, _u, h0s1?matrix_name:"");
                  } // echo
                  if (echo > 4) printf("# lowest and highest eigenvalue of the %s matrix is %g and %g %s, respectively\n", 
                                            matrix_name, lowest_eigenvalue*u, highest_eigenvalue*u, _u);
                  if (S == h0s1) {
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
                  } // S == h0s1
              } // stat_eig
          } // H or ovl_eig

      } // h0s1
      
      if (0 == stat_eig) {
                  // ToDo: export the entire spectrum and the lowest n_occupied_bands eigenvectors
                  //       for this k-point to be later used for the density generation.
                  //       n_occupied_bands should be chosen large so that also the tail 
                  //       of the Fermi-Dirac distribution is captured. We can check that
                  //       by ensuring that the occupation numbers of the next higher state
                  //          f_FD(eigvals[n_occupied_bands] - E_Fermi) < 10^{-15}
                  //       are small.
                
                  // Idea: exporting all n_occupied_bands eigenvectors may require
                  //       a lot of memory capacity and memory traffic during copying.
                  //       Alternatively, we could think about this scheme with higher data locality:
                  //       In iteration 0
                  //            we don't know the Fermi level yet, so we use bisection to find a 
                  //            E_Fermi_suggested for each kpoint separately
                  //            After iteration 0 and outside the k-loop, we can average
                  //            E_Fermi_suggested with the kpoint integration weights
                  //       In all further iterations
                  //            the Fermi level is an input.
                  //       With the Fermi level,
                  //       we compute the occupation number f_FD(eigvals[:] - E_Fermi)
                  //       to add the band densities to the total density
                  //       Furthermore, we also aggregate those
                  //       band densities for which d/d E_Fermi f_FD is nonzero
                  //       to a response density.
                  //       The response density is necessary because probably
                  //       the input Fermi level was not the one that produces the correct 
                  //       number of electrons requested. However, in particular close to
                  //       SCF convergence, we can linearize that.
                  //       The response density norm estimates the number of states at the (wrong)
                  //       Fermi level. So we can extract the necessary correction
                  //       to the Fermi level for the next iteration and
                  //       add as much of the response density to the total density
                  //       that the norm of the corrected total density matches the charge requested.
                
                  // Implementation:
#if 0               
                  // ToDo: test this
                  if (scf_iteration < 1) {
                      // in the first iteration, we assume that E_Fermi is unknown but q_electrons is given
                      assert( q_electrons >= 0 );
                      int const n_electrons = std::floor(q_electrons); // how many fully occupied states were there for kT=0
                      assert( n_electrons < nB );
                      int const n1 = std::min(n_electrons + 1, nB - 1);
                      double const w1 = q_electrons - n_electrons, w0 = 1 - w1; // linear interpolation weights
                      E_Fermi = w0*eigvals[n_electrons] + w1*eigvals[n1]; // start guess for bisection iterations
                      // find Fermi level by bisection (inside of this k-point spectrum only)
                      stat += fermi_distribution::Fermi_level(nB, eigvals.data(), kT, q_electrons, 
                                                              E_Fermi, nullptr, echo);
                  } // first SCF iteration

                  // assume that (now) E_Fermi is given
                  double q_density{0}, q_response{0};
                  for(int iB = 0; iB < nB; ++iB) {
                      double dfdE;
                      double const x = (eigval[iB] - E_Fermi)*beta; // beta is the inverse temperature
                      double const f_FD = fermi_distribution::FermiDirac(x, dfdE);
                      if (f_FD > 2e-16) {
                          q_density += f_FD;
                          // ToDO: add k_weight*f_FD*|\psi[iB]|^2 to the total density
                      }
                      if (std::abs(dfdE) > 2e-16) {
                          q_response += dfdE;
                          // ToDO: add k_weight*dfdE*|\psi[iB]|^2 to the response density
                      }
                  } // iB
                  weight_sum += k_weight; // to check the normalization of weights
                  density_sum += k_weight*q_density;
                  response_sum += k_weight*q_response;
                  E_Fermi_sum += k_weight*E_Fermi;
                  // export weight_sum, density_sum, response_sum, E_Fermi_sum
                  // 
                  // Why do we compute E_Fermi_sum?
                  //    After k-point summation E_Fermi_sum/weightsum will be a good guess
                  //    for the Fermi level at the end of the first SCF iteration
                  // But we could also just compute the correct Fermi level from all eigenvalues! Probably better
#endif
      } // success

      if (eigenenergies) {
          set(eigenenergies, std::min(size_t(nbands), eigvals.size()), eigvals.data()); // export
      } // eigenenergies
      
      return stat;
  } // solve
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template<typename real_t>
  status_t test_inverse(int const echo=0) {
      status_t status(0);
      int constexpr N = 5;
      real_t dev{0};
      view2D<std::complex<real_t>> mat(N, N, 0), inv(N, N);
      for(int n = 1; n <= N; ++n) { // dimension
          // fill with random values
          for(int i = 0; i < n; ++i) {
              for(int j = 0; j < n; ++j) {
                  auto const Re = simple_math::random<real_t>(-1, 1);
                  auto const Im = simple_math::random<real_t>(-1, 1);
                  mat(i,j) = std::complex<real_t>(Re, Im);
                  inv(i,j) = mat(i,j); // copy
              } // j
          } // i
          
          auto const stat = linear_algebra::inverse(n, inv.data(), inv.stride());
          if (stat) warn("inversion failed with status= %i", stat);
          status += stat;

          real_t devN{0}, devT{0};
          for(int i = 0; i < n; ++i) {
              for(int j = 0; j < n; ++j) {
                  std::complex<real_t> cN(0), cT(0);
                  for(int k = 0; k < n; ++k) {
                      cN += mat(i,k) * inv(k,j);
                      cT += inv(i,k) * mat(k,j);
                  } // k
                  std::complex<real_t> const diag = (i == j);
                  devN = std::max(devN, std::abs(cN - diag));
                  devT = std::max(devT, std::abs(cT - diag));
                  if (echo > 9) printf("# i=%i j=%i a=%g %g \tinv=%g %g \tcN=%g %g \tcT=%g %g\n", i, j,        
                      std::real(mat(i,j)), std::imag(mat(i,j)), std::real(inv(i,j)), std::imag(inv(i,j)),
                      std::real(cN), std::imag(cN), std::real(cT), std::imag(cT) );
              } // j
          } // i
          if (echo > 3) printf("# %s n= %d deviations from unity are %.2e and %.2e transposed\n",
                                __func__, n, devN, devT);
          dev = std::max(dev, std::max(devN, devT));
      } // n
      if (echo > 0) printf("# %s up to N= %d deviations from unity are %.2e\n\n", __func__, N, dev);
      return status;
  } // test_inverse
  
  inline status_t all_tests(int const echo=0) {
      status_t status(0);
      status += test_inverse<double>(echo);
      status += test_inverse<float>(echo);
      return status;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace dense_solver
