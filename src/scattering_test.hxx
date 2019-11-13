#pragma once

#include <cstdio> // printf
#include <cassert> // assert
#include <cstdint> // uint8_t
// #include <math.h> // fmod

#include "radial_grid.hxx" // create_equidistant_radial_grid, find_grid_index
#include "bessel_transform.hxx" // transform_s_function
#include "finite_difference.hxx" // set_Laplacian_coefficients
#include "sho_radial.hxx" // radial_eigenstates, radial_normalization, expand_poly
#include "display_units.h" // eV, Ang
#include "inline_tools.hxx" // align<nBits>
#include "inline_math.hxx" // set
#include "sho_tools.hxx" // num_ln_indices
#include "radial_integrator.hxx" // integrate_outwards<SRA>
#include "constants.hxx" // pi
#include "linear_algebra.hxx" // linear_solve, generalized_eigenvalues

typedef int status_t;

namespace scattering_test {

  int constexpr TRU=0, SMT=1, TRU_AND_SMT=2;

  template<typename real_t>
  inline real_t arcus_tangent(real_t const nom, real_t const den) {
      return (den < 0) ? std::atan2(-nom, -den) : std::atan2(nom, den); }
  
  template<typename real_t>
  inline size_t count_nodes(size_t const n, real_t const f[]) {
      size_t nnodes = 0;
      for(auto i = 1u; i < n; ++i) {
          nnodes += (f[i - 1]*f[i] < 0);
      };
      return nnodes;
  } // count_nodes

  inline status_t expand_sho_projectors(double prj[], int const stride, radial_grid_t const &rg, double const sigma,
       int const n, int const nrns[], int const ells[], int const rpow=0, int const echo=9) {
      
      status_t stat = 0;
      int maxpoly = 0;
      for(int i = 0; i < n; ++i) {
          maxpoly = std::max(maxpoly, 1 + nrns[i]);
      } // i
      double const siginv = 1./sigma;
      double const sigma_m23 = std::sqrt(pow3(siginv));
      auto const poly = new double[maxpoly];
      for(int i = 0; i < n; ++i) {
          int const ell = ells[i];
          int const nrn = nrns[i];
          stat += sho_radial::radial_eigenstates(poly, nrn, ell);
          double const f = sho_radial::radial_normalization(poly, nrn, ell) * sigma_m23;
          double norm = 0;
          if (echo > 8) printf("\n## ell=%i nrn=%i projectors on radial grid: r, p_ln(r):\n", ell, nrn);
          for(int ir = 0; ir < rg.n; ++ir) {
              double const r = rg.r[ir], dr = rg.dr[ir];
//            double const dr = 0.03125, r = dr*ir; // equidistant grid
              double const x = siginv*r, x2 = pow2(x);
              double const Gaussian = (x2 < 160) ? std::exp(-0.5*x2) : 0;
              prj[i*stride + ir] = f * sho_radial::expand_poly(poly, 1 + nrn, x2) * Gaussian * intpow(x, ell);
              if (echo > 8) printf("%g %g\n", r, prj[i*stride + ir]);
              norm += pow2(r*prj[i*stride + ir]) * dr;
              prj[i*stride + ir] *= intpow(r, rpow);
          } // ir
          if (echo > 8) printf("\n\n");
          if (echo > 5) printf("# projector normalization of ell=%i nrn=%i is %g\n", ell, nrn, norm);
      } // i
      delete[] poly;
      return stat;
  } // expand_sho_projectors
  

// #define _SELECTED_ENERGIES_LOGDER
  
  inline status_t logarithmic_derivative(
                radial_grid_t const *const rg[TRU_AND_SMT] // radial grid descriptors for Vtru, Vsmt
              , double        const *const rV[TRU_AND_SMT] // true and smooth potential given on the radial grid *r
              , double const sigma // sigma spread of SHO projectors
              , int const lmax // ellmax or numax of SHO projectors
              , uint8_t const nn[] // number of projectors per ell
              , double const aHm[] // non-local Hamiltonian elements in ln_basis
              , double const aSm[] // non-local overlap matrix elements
              , double const energy_range[3] // {lower, step, upper}
#ifdef  _SELECTED_ENERGIES_LOGDER
              , int const echo_level=2) {
                int const echo = echo_level + 10; // turn a lot of output on
#else
              , int const echo=2) {
#endif
    
      if (echo > 1) printf("\n# %s %s lmax=%i\n", __FILE__, __func__, lmax); 
      double const one_over_pi = 1./constants::pi;
      status_t stat = 0;
      
      double const dE = std::max(1e-9, energy_range[1]);

#ifdef  _SELECTED_ENERGIES_LOGDER
      int const nen = 6;
//    double const energy_list[nen] = {-0.221950, 0.047733, -0.045238,  0.120905, -0.359751,  0.181009};
//    int const ell_start = 0;
      double const energy_list[nen] = {-0.359751, -0.221950, 0.047733, -0.045238,  0.120905, 0.181009}; // check the Cu d-eigenstate
      int const ell_start = 2;
#else
      int const nen = (int)std::ceil((energy_range[2] - energy_range[0])/dE);
      if (echo > 0) printf("\n## logarithmic_derivative from %.3f to %.3f in %i steps of %g %s\n", 
                        energy_range[0]*eV, (energy_range[0] + nen*dE)*eV, nen + 1, dE*eV, _eV);
      int const ell_start = 0; // ToDo: delete this variable after debugging
#endif
      
      int const nr_diff = rg[TRU]->n - rg[SMT]->n; assert(nr_diff >= 0);
      int const mr = align<2>(rg[TRU]->n);
      auto const gg = new double[2*mr], ff = &gg[mr]; // greater and smaller component, TRU grid
      int ir_stop[TRU_AND_SMT];

      auto linsolfail = std::vector<size_t>(1 + lmax, 0);
      
      int nln = 0; for(int ell = 0; ell <= lmax; ++ell) nln += nn[ell];
      int const stride = mr;

      int const node_count = 1; // 1 or 0 switch
      auto const rphi = new double[node_count*9*rg[SMT]->n + 1];
      auto const rtru = new double[node_count*rg[TRU]->n];

      auto const rprj = new double[nln*stride]; // mr might be much larger than needed since mr is taken from the TRU grid
      { // scope: preparation for the projector functions
          int ells[nln], nrns[nln], iln = 0;
          for(int ell = 0; ell <= lmax; ++ell) {
              for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                  ells[iln] = ell;
                  nrns[iln] = nrn;
                  ++iln;
              } // nrn
          } // ell
          assert(nln == iln);
          stat += expand_sho_projectors(rprj, stride, *rg[SMT], sigma, nln, nrns, ells, 1, 0);
      } // scope

      ir_stop[SMT] = std::min(radial_grid::find_grid_index(*rg[SMT], 9*sigma), rg[SMT]->n - 2);
      double const Rlog = rg[SMT]->r[ir_stop[SMT]];
      if (echo > 0) printf("# %s check at radius %g %s\n", __func__, Rlog*Ang, _Ang);
      ir_stop[TRU] = ir_stop[SMT] + nr_diff;

      for(int ien = 0; ien <= nen; ++ien) {
#ifdef  _SELECTED_ENERGIES_LOGDER
          auto const energy = energy_list[ien];
          if (echo > 0) printf("# %s energy = ", __func__);
#else
          auto const energy = energy_range[0] + ien*dE;
#endif
//        if (echo > 0) printf("# node-count at %.6f %s", energy*eV, _eV);
          
          if (echo > 0) printf("%.6f", energy*eV);
          int iln_off = 0;
          for(int ell = ell_start; ell <= lmax; ++ell) 
          { // ell-loop
              double dg[TRU_AND_SMT], vg[TRU_AND_SMT];
              double deriv[9], value[9];
              double mat[99], gfp[99], x[9];
              double mat2[99];
              int nnodes[9];
              int const n = nn[ell];
              assert(n < 9);
              for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                  int success = 1;

                  for(int jrn = 0; jrn <= n*ts; ++jrn) {
                      bool const inhomgeneous = ((SMT == ts) && (jrn > 0));
//                       printf("# find %s %shomgeneous solution for ell=%i E=%g %s\n", (TRU == ts)?"true":"smooth",
//                              (inhomgeneous)?"in":"", ell, energy*eV,_eV); // DEBUG
                      int const nrn = jrn - 1, iln = iln_off + nrn;
                      double const *const rp = inhomgeneous ? &rprj[iln*stride] : nullptr;
                      int constexpr SRA = 1; // use the scalar relativistic approximation
                      nnodes[ts + jrn] = radial_integrator::integrate_outwards<SRA>(*rg[ts], rV[ts], ell, energy, 
                                                                     gg, ff, ir_stop[ts], &deriv[jrn], rp);
//                    if (TRU == ts) nnodes[TRU] = count_nodes(ir_stop[TRU] + 1, gg);
                      if ((TRU == ts) && node_count) set(rtru, ir_stop[TRU] + 1, gg);
                      value[jrn] = gg[ir_stop[ts]]; // value of the greater component at Rlog
                      if (SMT == ts) {
                          if (node_count) set(&rphi[jrn*rg[SMT]->n], ir_stop[ts] + 1, gg); // store radial solution
                          for(int krn = 0; krn < n; ++krn) {
                              // compute the inner products of the  projectors rprj with the solution gg
                              gfp[jrn*n + krn] = dot_product(ir_stop[SMT] + 1, &rprj[(iln_off + krn)*stride], gg, rg[SMT]->dr);
                              if (echo > 8) printf("# scattering solution for ell=%i E=%g %s <%i|%i> %g\n", 
                                                          ell, energy*eV,_eV, jrn, 1+krn, gfp[jrn*n + krn]);
                          } // krn

                          if (echo > 99) {
                              printf("\n## %shomogeneous scattering solution for ell=%i E=%g %s (r,phi):\n", 
                                    (jrn)?"in":"", ell, energy*eV,_eV);
                              for(int ir = 1; ir <= ir_stop[SMT]; ++ir) {
                                  auto const r = rg[SMT]->r[ir];
                                  printf("%g %g %g %g\n", r, gg[ir], (energy*r - rV[SMT][ir])*gg[ir], (jrn > 0)?rprj[iln*stride + ir]:0);
                              }   printf("\n\n");
                          } // echo
                          
                      } // SMT == ts
                  } // jrn

                  set(x, 1 + n, 0.0); // clear
                  x[0] = 1.0;
                  int nx = 1;
                  if ((SMT == ts) && n > 0) {
                      int const n0 = 1 + n;
                      set(mat, n0*n0, 0.0); // clear
                      for(int i = 0; i < n0; ++i) {
                          mat[i*n0 + i] = 1.0; // unity matrix
                      } // i
                      for(int irn = 0; irn < n; ++irn) {
                          for(int jrn = 0; jrn < n0; ++jrn) {
                              for(int krn = 0; krn < n; ++krn) {
                                  mat[jrn*n0 + (1 + irn)] += gfp[jrn*n + krn] * ( aHm[(irn + iln_off)*nln + (krn + iln_off)]
                                                                       - energy * aSm[(irn + iln_off)*nln + (krn + iln_off)] );
                              } // krn
                          } // jrn
                      } // irn

                      set(mat2, n0*n0, mat); // copy
                      auto const solving_status = linear_algebra::linear_solve(n0, mat, n0, x, n0);
                      stat += solving_status;
                      nnodes[SMT] = 0;
                      if (0 == solving_status) {
                          nx = n0; // successful --> linear combination with all n0 elements

                          if (echo > 8) { // show the problem
                              for(int i = 0; i < n0; ++i) {
                                  printf("# scattering_test mat[%i] ", i);
                                  for(int j = 0; j < n0; ++j) {
                                      printf("\t%g", mat2[i*n0 + j]);
                                  }   
                                  printf(" \t x[%i] = %g", i, x[i]);
                                  if (i > 0) {
                                      printf(" \t HmES[%i] ", i);
                                      for(int j = 0; j < n; ++j) printf("\t%g", aHm[(i-1 + iln_off)*nln + (j + iln_off)] 
                                                                     - energy * aSm[(i-1 + iln_off)*nln + (j + iln_off)]);
                                  } // i > 0
                                  printf("\n");
                              } // i
                          } // echo
                          
                          if (echo > 7) { // scope: verify
                              printf("# scattering_test linear solve test: ");
                              for(int i = 0; i < n0; ++i) {
                                  double bi = 0; for(int j = 0; j < n0; ++j) bi += x[j]*mat2[j*n0 + i];
                                  printf(" %g", bi);
                              }   printf("\n");
                          } // scope

                          if (node_count) {
                              scale(rphi, ir_stop[SMT] + 1, x[0]); // scale the homogeneous solution with x[0] (which is usually == 1)
                              if (echo > 9) printf("# scattering solution for ell=%i E=%g %s coefficients %g", ell, energy*eV,_eV, x[0]);
                              for(int i = 1; i < n0; ++i) {
                                  if (echo > 9) printf(" %g", x[i]);
                                  add_product(rphi, ir_stop[SMT] + 1, &rphi[i*rg[SMT]->n], x[i]); // add inhom. solutions
                              } // i
                              if (echo > 9) printf("\n");
                              nnodes[SMT] = count_nodes(ir_stop[SMT] + 1, rphi);
                              if (echo > 9) {
                                  auto const scal = rtru[ir_stop[TRU]]/rphi[ir_stop[SMT]]; // match in value at end point
                                  auto const s = 1./std::sqrt(dot_product(ir_stop[TRU] + 1, rtru, rtru, rg[TRU]->dr)); // normalize
                                  printf("\n## scattering solution for ell=%i E=%g %s (r,tru,smt,diff):\n", ell, energy*eV,_eV);
                                  for(int ir = 1; ir <= ir_stop[SMT]; ++ir) {
                                      printf("%g\t%g %g %g\n", rg[SMT]->r[ir], s*rtru[ir + nr_diff] , rphi[ir]*scal*s, 
                                                                                (rtru[ir + nr_diff] - rphi[ir]*scal)*s);
//                                    printf("\t%g %g\n", rV[TRU][ir + nr_diff], rV[SMT][ir]); // compare potentials
                                  }   printf("\n\n");
                              } // echo
                          } // node_count

                      } else {
                          ++linsolfail[ell];
                          success = 0;
                      }
                  } // SMT == ts

                  vg[ts] = dot_product(nx, x, value); // value of the greater component at Rlog
                  dg[ts] = dot_product(nx, x, deriv); // derivative

                  double const generalized_node_count = success*(0*nnodes[ts] + 0.5 - one_over_pi*arcus_tangent(dg[ts], vg[ts]));
#ifdef  _SELECTED_ENERGIES_LOGDER
                  if (echo > 0) printf("# %cL(ell=%i) =", ts?'~':' ', ell);
#endif
                  if (echo > 0) printf("%c%.6f", (ts)?' ':'\t', generalized_node_count);
#ifdef  _SELECTED_ENERGIES_LOGDER
                  if (echo > 0) printf("\n");
#endif
              } // ts
//            if (echo > 0) printf("# %i %g %g %g %g\n", ell, dg[TRU], vg[TRU], dg[SMT], vg[SMT]);
#ifdef  _SELECTED_ENERGIES_LOGDER
              printf("\n\n\n# EXIT at %s line %i \n\n", __FILE__, __LINE__); exit(__LINE__); // DEBUG
#endif              
              iln_off += n;
          } // ell
          if (echo > 0) printf("\n");
      } // ien

      for(int ell = 0; ell <= lmax; ++ell) {
          if (linsolfail[ell]) {
              printf("# %s linear solve failed %ld times for ell=%i\n", __func__, linsolfail[ell], ell);
          } // fail
      } // ell

      delete[] rphi;
      return stat;
  } // logarithmic_derivative
  
  
  
  
  
  inline status_t eigenstate_analysis(radial_grid_t const& gV // grid descriptor for Vsmt
              , double const Vsmt[] // smooth potential given on radial grid
              , double const sigma // sigma spread of SHO projectors
              , int const lmax // ellmax or numax of SHO projectors
              , uint8_t const nn[] // number of projectors per ell
              , double const aHm[] // non-local Hamiltonian elements in ln_basis
              , double const aSm[] // non-local overlap matrix elements
              , int const nr=384 // number of radial grid points in equidistance mesh
              , int const Vshift=0 // potential shift
              , int const echo=2
    ) {
      status_t stat = 0;
      auto const g = *radial_grid::create_equidistant_radial_grid(nr + 1, gV.rmax);
      auto const dr = g.dr[0]; // in an equidistant grid, the grid spacing is constant and, hence, indepent of ir
      if (echo > 1) printf("\n# %s %s dr=%g nr=%i rmax=%g %s\n", __FILE__, __func__, dr*Ang, nr, dr*nr*Ang, _Ang); 
      
      auto const Vloc = new double[g.n];
      { // scope: interpolate to the equidistant grid by Bessel-transform
          int const nq = nr/2; double const dq = .125;
          auto const Vq = new double[nq];
          stat += bessel_transform::transform_s_function(Vq, Vsmt, gV, nq, dq, false, 0);
          stat += bessel_transform::transform_s_function(Vloc, Vq, g, nq, dq, true, 0); // back=true
          for(int ir = 0; ir < g.n; ++ir) {
              Vloc[ir] += Vshift;
          } // ir
          if (echo > 8) {
              printf("\n## Vsmt:\n"); for(int ir = 1; ir < gV.n; ++ir) printf("%g %g\n", gV.r[ir], Vsmt[ir]); printf("\n\n");
              printf("\n## Vq:\n");   for(int iq = 0; iq < nq; ++iq) printf("%g %g\n", iq*dq, Vq[iq]); printf("\n\n");
              printf("\n## Vloc:\n"); for(int ir = 1; ir < nr; ++ir) printf("%g %g\n", g.r[ir], Vloc[ir]); printf("\n\n");
          } // echo
          delete[] Vq;
      } // scope
      
      int const stride = align<2>(nr); assert(stride >= nr);
      // allocate Hamiltonian and overlap matrix
      auto const Ham = new double[(3*nr + 1)*stride], Ovl = &Ham[1*nr*stride], 
                   Ovl_copy = &Ham[2*nr*stride],     eigs = &Ham[3*nr*stride];

      int nln = 0; int mprj = 0; 
      for(int ell = 0; ell <= lmax; ++ell) {
          nln += nn[ell];
          mprj = std::max(mprj, (int)nn[ell]);
      } // ell
      auto rprj = new double[mprj*stride]; // projector functions*r

      int const nFD = 4; double cFD[1 + nFD]; set(cFD, 1 + nFD, 0.0);
      stat += finite_difference::set_Laplacian_coefficients(cFD, nFD, dr);
      if (echo > 3) printf("# %s %s finite difference with %i neighbors\n", __FILE__, __func__, nFD); 

      int ln_off = 0;
      for(int ell = 0; ell <= lmax; ++ell) {
          set(Ham, 2*nr*stride, 0.0); // clear Hamiltonian and overlap matrix

          // setup the local Hamiltonian
          for(int ir = 0; ir < nr; ++ir) {
              double const r = g.r[ir + 1];
              // local potential and repulsive angular part of the kinetic energy
              Ham[ir*stride + ir] = Vloc[ir + 1] + 0.5*(ell*(ell + 1.)/pow2(r)); 
              Ovl[ir*stride + ir] = 1.;
              for(int jr = 0; jr < nr; ++jr) {
                  int const dij = std::abs(ir - jr);
                  if (dij <= nFD) Ham[ir*stride + jr] -= 0.5*cFD[dij]; // finite_difference kinetic energy operator
              } // jr
              // in the radial representation, the usual Laplacian d^2/dr^2 can be applied 
              // for the radial component if the wave functions are in r*phi(r) representation
          } // ir
          
          // generate normalized SHO projectors
          int const nprj = nn[ell];
          {
              int ells[nprj], nrns[nprj];
              for(int nrn = 0; nrn < nprj; ++nrn) {
                  nrns[nrn] = nrn;
                  ells[nrn] = ell;
              } // nrn
              stat += expand_sho_projectors(rprj, stride, g, sigma, nprj, nrns, ells, 1, echo);
          } // scope
          auto const rprj1 = rprj + 1; // forward the rprj-pointer by one so that ir=0 will access the first non-zero radius

          // add the non-local dyadic operators to the Hamiltonian and overlap
          for(int ir = 0; ir < nr; ++ir) {
              for(int jr = 0; jr < nr; ++jr) {
                  for(int nrn = 0; nrn < nprj; ++nrn) {
                      for(int mrn = 0; mrn < nprj; ++mrn) {
                          int const ijln = (ln_off + nrn)*nln + (ln_off + mrn);
                          Ham[ir*stride + jr] += rprj1[nrn*stride + ir]*aHm[ijln]*rprj1[mrn*stride + jr]*dr;
                          Ovl[ir*stride + jr] += rprj1[nrn*stride + ir]*aSm[ijln]*rprj1[mrn*stride + jr]*dr;
                      } // mrn
                  } // nrn
              } // jr
          } // ir
          
          set(Ovl_copy, nr*stride, Ovl); // copy
          
          { // scope: diagonalize
              // solve the generalized eigenvalue problem
              auto const info = linear_algebra::generalized_eigenvalues(nr, Ham, stride, Ovl, stride, eigs);
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
                      auto const evec = Ham; // eigenvectors are stored in the memory space that was used to input the Hamiltonian
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
                              printf("%12.6f", dot_product(nr, &evec[iev*stride], &rprj1[nrn*stride])*std::sqrt(dr));
                          }   printf("\n");
                      } // iev
                  } // echo
                  
              } else { // info
                  if (echo > 2) printf("# diagonalization for ell=%i returned info=%i\n", ell, info);
                  
                  // diagonalize Ovl_copy, standard eigenvalue problem
                  linear_algebra::eigenvalues(nr, Ovl_copy, stride, eigs);
                  if (eigs[0] <= 0) { // warn
                      if (echo > 0) printf("# %s ell=%i lowest eigenvalue of the overlap matrix is non-positive! %g\n", __func__, ell, eigs[0]);
                  } // overlap matrix is not positive definite
                  if (echo > 1) {
                      printf("# %s ell=%i eigenvalues of the overlap matrix", __func__, ell);
                      for(int iev = 0; iev < 8; ++iev) {
                          printf(" %g", eigs[iev]);
                      }   printf("\n");
                  } // echo
                  ++stat;
              } // info
          } // scope: diagonalize

          ln_off += nn[ell]; // forward
      } // ell

      delete[] rprj;
      delete[] Vloc;
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
      return 0;
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
      auto nn = std::vector<uint8_t>(1 + lmax, 0); for(int ell = 0; ell <= lmax; ++ell) nn[ell] = 1 + (lmax - ell)/2;
      return eigenstate_analysis(rg, V.data(), sigma, lmax, nn.data(), aHm.data(), aHm.data(), 128, echo);
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
