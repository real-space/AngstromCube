#pragma once

#include <cstdio> // printf
#include <cassert> // assert
#include <cstdint> // uint8_t

// #include "radial_integrator.hxx" // 
#include "radial_grid.hxx" // create_equidistant_radial_grid
#include "bessel_transform.hxx" // transform_s_function
#include "finite_difference.hxx" // set_Laplacian_coefficients
#include "sho_radial.hxx" // radial_eigenstates, radial_normalization, expand_poly
#include "display_units.h" // eV, Ang
#include "inline_tools.hxx" // align<nBits>
#include "inline_math.hxx" // set
#include "sho_tools.hxx" // num_ln_indices

typedef int status_t;

  extern "C" {
      // double symmetric generalized eigenvalue problem
      void dsygv_(int const* ITYPE, char const* JOBZ, char const* UPLO, int const* N, 
                  double* A, int const* LDA, double* B, int const* LDB, 
                  double* W, double* WORK, int const* LWORK, int* INFO);
      // double symmetric standard eigenvalue problem
      void dsyev_(char const* JOBZ, char const* UPLO, int const* N, double* A, int const* LDA,
                  double* W, double* WORK, int const* LWORK, int* INFO);
  } // LAPACK

namespace scattering_test {

  inline status_t logarithmic_derivative() {
      
  } // logarithmic_derivative
  
  inline status_t eigenstate_analysis(radial_grid_t const& gV // grid descriptor for Vsmt
              , double const Vsmt[] // smooth potential given on radial grid
              , double const sigma // sigma spread of SHO projectors
              , int const lmax // ellmax or numax of SHO projectors
              , uint8_t const nn[] // number of projectors per ell
              , double const aHm[] // non-local Hamiltonian elements in ln_basis
              , double const aSm[] // non-local overlap matrix elements
              , int const nr=384 // number of radial grid points in equidistance mesh
              , int const echo=2
    ) {
      status_t stat = 0;
      auto const g = *radial_grid::create_equidistant_radial_grid(nr + 1, gV.rmax);
      auto const dr = g.dr[0]; // in an equidistant grid, the grid spacing is constant and, hence, indepent of ir
      if (echo > 1) printf("\n# %s %s dr=%g nr=%i rmax=%g %s\n", __FILE__, __func__, dr*Ang, nr, dr*nr*Ang, _Ang); 
      
      auto const Vloc = new double[g.n];
      { // scope: interpolate to the equidistant grid by bessel transform
          int const nq = nr/2; double const dq = .125;
          auto const Vq = new double[nq];
          stat += bessel_transform::transform_s_function(Vq, Vsmt, gV, nq, dq, false, echo);
          stat += bessel_transform::transform_s_function(Vloc, Vq, g, nq, dq, true, echo); // back=true
          if (echo > 8) {
              printf("\n## Vsmt:\n"); for(int ir = 1; ir < gV.n; ++ir) printf("%g %g\n", gV.r[ir], Vsmt[ir]); printf("\n\n");
              printf("\n## Vq:\n");   for(int iq = 0; iq < nq; ++iq) printf("%g %g\n", iq*dq, Vq[iq]); printf("\n\n");
              printf("\n## Vloc:\n"); for(int ir = 1; ir < nr; ++ir) printf("%g %g\n", g.r[ir], Vloc[ir]); printf("\n\n");
          } // echo
          delete[] Vq;
      } // scope
      
      int const stride = align<2>(nr);
      // allocate Hamiltonian and overlap matrix
      auto const Ham = new double[(3*nr + 1)*stride], Ovl = &Ham[1*nr*stride], 
                           work = &Ham[2*nr*stride], eigs = &Ham[3*nr*stride];

      int nln = 0; int mprj = 0; 
      for(int ell = 0; ell <= lmax; ++ell) {
          nln += nn[ell];
          mprj = std::max(mprj, (int)nn[ell]);
      } // ell
      auto const rprj = new double[mprj*stride]; // projector functions*r
      auto const Gaussian = new double[nr]; // envelope function
      double const siginv = 1./sigma;
      for(int ir = 0; ir < nr; ++ir) {
          double const r = g.r[ir + 1];
          Gaussian[ir] = std::exp(-0.5*pow2(siginv*r));
      } // ir
      auto const poly = new double[mprj];

      int const nFD = 1; double cFD[1 + nFD];
      stat += finite_difference::set_Laplacian_coefficients(cFD, 1, dr);
      
      int ln_off = 0;
      for(int ell = 0; ell <= lmax; ++ell) {
          set(Ham, 2*nr*stride, 0.0); // clear Hamiltonian and overlap matrix

          // setup the local Hamiltonian
          for(int ir = 0; ir < nr; ++ir) {
              double const r = g.r[ir + 1];
              // local potential and repulsive angular part of the kinetic energy
              Ham[ir*stride + ir] = Vloc[ir + 1] + 0.5*(ell*(ell + 1.)/pow2(r)); 
              Ovl[ir*stride + ir] = 1.;
              for(int jr = std::max(0, ir - 1); jr <= std::min(ir + 1, nr - 1); ++jr) {
                  Ham[ir*stride + jr] += -0.5*cFD[std::abs(ir - jr)]; // finite_difference kinetic energy operator
              } // jr
              // in the radial representation, the usual Laplacian d^2/dr^2 can be applied 
              // for the radial component if the wave functions are in r*phi(r) representation
          } // ir
          
          // generate normalized SHO projectors
          int const nprj = nn[ell];
          for(int nrn = 0; nrn < nprj; ++nrn) {
              stat += sho_radial::radial_eigenstates(poly, nrn, ell);
              auto const f = sho_radial::radial_normalization(poly, nrn, ell) * std::sqrt(pow3(siginv));
              double norm2 = 0;
              if (echo > 7) printf("\n## %s ell=%i projector nrn=%i:\n", __func__, ell, nrn);
              for(int ir = 0; ir < nr; ++ir) {
                  double const r = g.r[ir + 1];
                  double &p = rprj[nrn*stride + ir];
                  rprj[nrn*stride + ir] = f*sho_radial::expand_poly(poly, 1 + nrn, pow2(siginv*r));
                  // now multiply what is missing: Gaussian envelope and r^{ell + 1} factors
                  rprj[nrn*stride + ir] *= Gaussian[ir] * intpow(siginv*r, ell) * r;
                  norm2 += pow2(rprj[nrn*stride + ir]);
                  if (echo > 7) printf("%g %g\n", r, rprj[nrn*stride + ir]);
              } // ir
              if (echo > 7) printf("\n\n");
              if (echo > 6) printf("# SHO projector for ell=%i nrn=%i is normalized to %g\n", ell, nrn, norm2*dr);
          } // nrn

          // add the non-local dyadic operators to the Hamiltonian and overlap
          for(int ir = 0; ir < nr; ++ir) {
              for(int jr = 0; jr < nr; ++jr) {
                  for(int nrn = 0; nrn < nprj; ++nrn) {
                      for(int mrn = 0; mrn < nprj; ++mrn) {
                          int const ijln = (ln_off + nrn)*nln + (ln_off + mrn);
                          Ham[ir*stride + jr] += rprj[nrn*stride + ir]*aHm[ijln]*rprj[mrn*stride + jr]*dr;
                          Ovl[ir*stride + jr] += rprj[nrn*stride + ir]*aSm[ijln]*rprj[mrn*stride + jr]*dr;
                      } // mrn
                  } // nrn
              } // jr
          } // ir
          
          {   // scope: diagonalize
              int const lwork = nr*stride, itype = 1; // for the LAPACK solver
              char const uplo = 'u', jobv = 'v'; int info = 0;
              
              // solve the generalized eigenvalue problem
              dsygv_(&itype, &jobv, &uplo, &nr, Ham, &stride, Ovl, &stride, eigs, work, &lwork, &info);
//            dsyev_(&jobv, &uplo, &nr, Ham, &stride, eigs, work, &lwork, &info);  // standard EV problem

              if (0 == info) {
                
                  int const nev = 5 - ell/2; // show less eigenvalues for higher ell-states
                  if (echo > 1) {
                      printf("# lowest eigenvalues for ell=%i  ", ell);
                      for(int iev = 0; iev < nev; ++iev) {
                          printf("  %.6f", eigs[iev]*eV);
                      }   printf("  %s\n", _eV);
                  } // echo
                  
                  if (echo > 2) {
                      // projection analysis for the lowest nev eigenvectors
                      auto const evec = Ham; // eigenvectors are store in the memory space that was used to input the Hamiltonian
                      for(int iev = 0; iev < nev; ++iev) {
                          // plot eigenvectors
                          if (echo > 8) {
                              printf("\n## %s ell=%i eigenvalue %.6f %s %i-th eigenvector:\n", __func__, ell, eigs[iev]*eV, _eV, iev);
                              for(int ir = 0; ir < nr; ++ir) {
                                  printf("%g %g\n", g.r[ir + 1], evec[iev*stride + ir]);
                              }   printf("\n\n");
                          } // echo
                          
                          printf("# projection analysis for ell=%i eigenvalue (#%i) %.6f %s  coefficients ", ell, iev, eigs[iev]*eV, _eV);
                          for(int nrn = 0; nrn < nprj; ++nrn) {
                              printf("%12.6f", dot_product(nr, &evec[iev*stride], &rprj[nrn*stride])*std::sqrt(dr));
                          }   printf("\n");
                      } // iev
                  } // echo
                  
              } else { // info
                  if (echo > 2) printf("# diagonalization for ell=%i returned info=%i\n", ell, info);
                  ++stat;
              } // info
          } // scope: diagonalize
          
          ln_off += nn[ell]; // forward
      } // ell

      delete[] rprj;
      delete[] poly;
      delete[] Vloc;
      delete[] Gaussian;
      delete[] Ham;
      
      return stat;
  } // eigenstate_analysis
  
  
  template<typename real_t>
  status_t emm_average(real_t Mln[], real_t const Mlmn[], int const lmax, uint8_t const nn[], int const stride=-1) 
  {
      int const echo = 9;
      if (echo > 4) printf("# %s: ellmax = %i\n", __func__, lmax);
      int iln = 0, ilmn = 0;
      for(int ell = 0; ell <= lmax; ++ell) {
          iln += nn[ell];
          ilmn += (2*ell + 1)*nn[ell];
      } // ell
      int const nln = iln, nlmn = ilmn;
      int const M_stride = (stride < 0) ? nlmn : stride;
      auto ln_list = std::vector<int>(nlmn, 0);
      auto lf_list = std::vector<real_t>(nln, 0);
      { // scope: fill the ln_list
          int ilmn = 0, ln_off = 0;
          if (echo > 6) printf("# %s: ln_list ", __func__);
          for(int ell = 0; ell <= lmax; ++ell) {
              for(int emm = -ell; emm <= ell; ++emm) {
                  for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                      ln_list[ilmn] = ln_off + nrn;
                      if (echo > 6) printf(" %i", ln_list[ilmn]);
                      ++ilmn;
                  } // nrn
              } // emm
              for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                  lf_list[ln_off] = real_t(1)/std::sqrt(2*ell + real_t(1));
                  ++ln_off;
              } // nrn
          } // ell
          if (echo > 6) printf("\n");
          assert(nlmn == ilmn);
          if (echo > 5) {
              printf("# %s: l2p1_list ", __func__);
              for(int iln = 0; iln < nln; ++iln) {
                  printf(" %.1f", 1./pow2(lf_list[iln]));
              }   printf("\n");
          } // echo
      } // scope
      
      { // scope: now average
          set(Mln, nln*nln, real_t(0)); // clear
          for(int ilmn = 0; ilmn < nlmn; ++ilmn) {      int const iln = ln_list[ilmn];
              for(int jlmn = 0; jlmn < nlmn; ++jlmn) {  int const jln = ln_list[jlmn];
                  Mln[iln*nln + jln] += Mlmn[ilmn*M_stride + jlmn];
              } // jlmn
          } // ilmn
          // apply (2*ell + 1) denominators
          for(int iln = 0; iln < nln; ++iln) {
              for(int jln = 0; jln < nln; ++jln) {
                  Mln[iln*nln + jln] *= lf_list[iln]*lf_list[jln];
              } // jln
          } // iln
      } // scope
      
  } // emm_average
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS
  
  inline status_t test_eigenstate_analysis(int const echo=3, int const lmax=7) {
      if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
      auto const rg = *radial_grid::create_default_radial_grid(0);
      int const nln = sho_tools::num_ln_indices(lmax);
      auto V = std::vector<double>(rg.n, 0.0);
      double const sigma = 1.0; // if the rmax ~= 10, lmax = 7, sigma <= 1.5, otherwise projectors leak out
      for(int ir = 0; ir < rg.n; ++ir) V[ir] = 0.5*(pow2(rg.r[ir]) - pow2(rg.rmax))/pow4(sigma); // harmonic potential 
      auto const aHm = std::vector<double>(nln*nln, 0.0);
      uint8_t nn[1 + lmax]; for(int ell = 0; ell <= lmax; ++ell) nn[ell] = 1 + (lmax - ell)/2;
      return eigenstate_analysis(rg, V.data(), sigma, lmax, nn, aHm.data(), aHm.data(), 128, echo);
      // expected result: eigenstates at E0 + (2*nrn + ell)*sigma^-2 Hartree
  } // test_eigenstate_analysis

  inline status_t all_tests(int const echo=3) {
      if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
      auto status = 0;
      status += test_eigenstate_analysis(echo);
      return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace scattering_test
