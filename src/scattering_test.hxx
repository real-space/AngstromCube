#pragma once

#include <cstdio> // printf
#include <cassert> // assert
#include <cstdint> // uint8_t
#include <cmath> // std::copysign

#include "radial_grid.hxx" // ::create_equidistant_radial_grid, ::find_grid_index
#include "bessel_transform.hxx" // ::transform_s_function
#include "finite_difference.hxx" // ::set_Laplacian_coefficients
#include "sho_radial.hxx" // ::radial_eigenstates, ::radial_normalization, ::expand_poly
#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_tools.hxx" // align<nBits>
#include "inline_math.hxx" // set
#include "sho_tools.hxx" // ::nSHO, ::nSHO_radial
#include "radial_integrator.hxx" // ::integrate_outwards<SRA>
#include "constants.hxx" // ::pi
#include "linear_algebra.hxx" // ::linear_solve, ::eigenvalues
#include "data_view.hxx" // view2D
#ifdef DEVEL
    #include "control.hxx" // ::get // ToDo: remove this again
#endif

#define DEBUG
#ifdef  DEBUG
  #include "debug_output.hxx" // here
#endif

#include "status.hxx" // status_t

namespace scattering_test {

  int constexpr TRU=0, SMT=1, TRU_AND_SMT=2;

  auto const ellchar = "spdfghijklmno"; // ToDo: by convention 'i' is not the proper char for ell=6
  
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

  inline status_t expand_sho_projectors(double prj[], int const stride, radial_grid_t const &rg, 
                    double const sigma, int const numax, int const rpow=0, int const echo=9) {
      
      status_t stat = 0;
      double const siginv = 1./sigma;
      double const sigma_m23 = std::sqrt(pow3(siginv));
      int const maxpoly = align<2>(1 + numax/2);
      int const nln = sho_tools::nSHO_radial(numax);
      view2D<double> poly(nln, maxpoly);
      std::vector<double> norm(nln, 0.0);
      std::vector<int> nrns(nln), ells(nln);
      for(int ell = 0; ell <= numax; ++ell) {
          for(int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
              int const iln = sho_tools::ln_index(numax, ell, nrn);
              stat += sho_radial::radial_eigenstates(poly[iln], nrn, ell);
              auto const norm_factor = sho_radial::radial_normalization(poly[iln], nrn, ell) * sigma_m23;
              scale(poly[iln], nrn + 1, norm_factor);
          } // nrn
      } // ell

      assert(rg.n <= stride);
      if (echo > 18) printf("\n## SHO projectors on radial grid: r, p_00(r), p_01, ... :\n");
      for(int ir = 0; ir < rg.n; ++ir) {
          double const r = rg.r[ir], r2dr = rg.r2dr[ir];
//        double const dr = 0.03125, r = dr*ir, r2dr = pow2(r)*dr; // equidistant grid
          double const x = siginv*r, x2 = pow2(x);
          double const Gaussian = (x2 < 160) ? std::exp(-0.5*x2) : 0;
          if (echo > 18) printf("%g ", r);
          auto const r_pow_rpow = intpow(r, rpow);
          for(int ell = 0; ell <= numax; ++ell) {
              auto const x_pow_ell = intpow(x, ell);
              for(int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
                  int const iln = sho_tools::ln_index(numax, ell, nrn);
                  double const projector_value = sho_radial::expand_poly(poly[iln], 1 + nrn, x2) * Gaussian * x_pow_ell;
                  if (echo > 18) printf(" %g", projector_value);
                  norm[iln] += pow2(projector_value) * r2dr;
                  prj[iln*stride + ir] = projector_value * r_pow_rpow;
              } // nrn
          } // ell
          if (echo > 18) printf("\n");
      } // ir
      if (echo > 18) printf("\n\n");
      if (echo > 9) {
          printf("# projector normalizations are ");
          for(int iln = 0; iln < nln; ++iln) {
              printf(" %g", norm[iln]);
          }   printf("\n");
      } // echo
      return stat;
  } // expand_sho_projectors
  
// #define _SELECTED_ENERGIES_LOGDER
  
  inline double find_outwards_solution(radial_grid_t const & rg, double const rV[] // effective potential r*V(r)
      , int const ell, double const energy // // angular moment quantum number and energy
      , double gg[], double ff[] // work arrays, greater and smaller component
      , int const ir_stop=-1 // radius where to stop
      , double const *inh=nullptr) // r*inhomogeneiety
  {
      int constexpr SRA = 1; // use the scalar relativistic approximation
      double deriv{0};
      radial_integrator::integrate_outwards<SRA>(rg, rV, ell, energy, gg, ff, ir_stop, &deriv, inh);
      return deriv;
  } // find_outwards_solution

  template<bool NodeCount=false>
  inline double generalized_node_count_TRU(radial_grid_t const & rg, double const rV[] // effective potential r*V(r)
      , int const ell, double const energy // // angular moment quantum number and energy
      , double gg[], double ff[] // work arrays, greater and smaller component
      , int const ir_stop // radius where to stop
      , int const echo=0) {

      if (echo > 9) printf("# find true homgeneous solution for ell=%i E=%g %s\n", ell, energy*eV,_eV); // DEBUG
      double const deriv = find_outwards_solution(rg, rV, ell, energy, gg, ff, ir_stop, nullptr);
      double const value = gg[ir_stop]; // value of the greater component at Rlog
      int nnodes{0}; if (NodeCount) nnodes = count_nodes(ir_stop + 1, gg);
      double constexpr one_over_pi = 1./constants::pi;
      return (nnodes + 0.5 - one_over_pi*arcus_tangent(deriv, value));
  } // generalized_node_count_TRU

  template <bool NodeCount=false, int const nLIM=8>
  inline double generalized_node_count_SMT(radial_grid_t const & rg, double const rV[] // effective potential r*V(r)
      , int const ell, double const energy // // angular moment quantum number and energy
      , double gg[], double ff[] // work arrays, greater and smaller component
      , int const ir_stop // radius where to stop
      , view2D<double> const & rprj // pointer to inhomogeneieties
      , int const n=0 // number of projector >= 0
      , double const *aHm=nullptr // pointer to Hamiltonian, can be zero if n=0
      , double const *aSm=nullptr // pointer to overlap, can be zero if n=0
      , int const stride=0 // stride for Hamiltonian and overlap
      , int const echo=0) {

      double deriv[nLIM], value[nLIM], gfp[nLIM*nLIM]; assert(n < nLIM);

      view2D<double> waves(NodeCount ? (1 + n) : 0, align<2>(ir_stop + 1)); // get memory

      for(int jrn = n; jrn >= 0; --jrn) { // >0:inhomogeneous, 0: homogeneous
          if (echo > 9) printf("# find smooth %shomgeneous solution for ell=%i E=%g %s (%i)\n", (jrn > 0)?"in":"", ell, energy*eV,_eV, jrn-1);
          auto const ggj = NodeCount ? waves[jrn] : gg;
          deriv[jrn] = find_outwards_solution(rg, rV, ell, energy, ggj, ff, ir_stop, (jrn > 0) ? rprj[jrn - 1] : nullptr);
          value[jrn] = ggj[ir_stop]; // value of the greater component at Rlog
          for(int krn = 0; krn < n; ++krn) {
              // compute the inner products of the projectors rprj with the solution ggj
              gfp[jrn*n + krn] = dot_product(ir_stop + 1, rprj[krn], ggj, rg.dr);
          } // krn
      } // jrn

      int const n0 = 1 + n;
      double x[nLIM];
      set(x, n0, 0.0); // clear
      x[0] = 1.0;
      if (n > 0) {
          double mat[nLIM*nLIM];
          set(mat, n0*n0, 0.0); // clear
          for(int i = 0; i < n0; ++i) {
              mat[i*n0 + i] = 1.0; // unit matrix
          } // i
          for(int irn = 0; irn < n; ++irn) {
              for(int jrn = 0; jrn < n0; ++jrn) {
                  for(int krn = 0; krn < n; ++krn) {
                      mat[jrn*n0 + (1 + irn)] += gfp[jrn*n + krn] * ( aHm[irn*stride + krn]
                                                           - energy * aSm[irn*stride + krn] );
                  } // krn
              } // jrn
          } // irn
          auto const solving_status = linear_algebra::linear_solve(n0, mat, n0, x, n0);
          if (solving_status) {
              warn("linear solve failed for E= %g %s", energy*eV, _eV);
              return 0; // failed
          } // solving_status is non-zero
      } // n > 0
      double const val = dot_product(n0, x, value); // value of the greater component at Rlog
      double const der = dot_product(n0, x, deriv); // derivative
      int nnodes{0}; 
      if (NodeCount) {
          set(gg, rg.n, 0.0); // clear
          for(int jrn = 0; jrn <= n; ++jrn) {
              add_product(gg, ir_stop + 1, waves[jrn], x[jrn]);
          } // jrn
          nnodes = count_nodes(ir_stop + 1, gg); 
      } // NodeCount
      double constexpr one_over_pi = 1./constants::pi;
      return (nnodes + 0.5 - one_over_pi*arcus_tangent(der, val));
  } // generalized_node_count_SMT

  inline status_t logarithmic_derivative(
                radial_grid_t const *const rg[TRU_AND_SMT] // radial grid descriptors for Vtru, Vsmt
              , double        const *const rV[TRU_AND_SMT] // true and smooth potential given on the radial grid *r
              , double const sigma // sigma spread of SHO projectors
              , int const lmax // ellmax up to which the analysis should go
              , uint8_t const nn_input[] // number of projectors per ell
              , int const numax
              , double const aHm[] // non-local Hamiltonian elements in ln_basis
              , double const aSm[] // non-local overlap matrix elements
              , double const energy_range[3] // {lower, step, upper}
              , char const *label=""
              , int const echo=2 // log-level
              , float const Rlog_over_sigma=6) {

      status_t stat(0);
      
      double const dE = std::copysign(std::max(1e-9, std::abs(energy_range[1])), energy_range[1]);
      int const nen = int(std::ceil((energy_range[2] - energy_range[0])/dE));

      if (echo > 1) printf("\n# %s %s energy range from %g to %g %s in %d steps of %g %s, lmax=%d\n", 
          label, __func__, energy_range[0]*eV, energy_range[2]*eV, _eV, 1 + nen, dE*eV, _eV, lmax);

      if (nen < 0) return stat; // empty range

      int const nr_diff = rg[TRU]->n - rg[SMT]->n; assert(nr_diff >= 0);
      int const mr = align<2>(rg[TRU]->n);
      std::vector<double> gg(mr), ff(mr); // greater and smaller component, TRU grid
      int ir_stop[TRU_AND_SMT];

      int const nln = sho_tools::nSHO_radial(numax);
      int const stride = mr;

      int nn[12]; for(int ell = 0; ell < 12; ++ell) { nn[ell] = (numax + 2 - ell)/2; } // == nnmax(numax, ell)
      if (echo > 9) printf("# %s nn = %d %d %d %d %d %d %d %d\n", __func__, nn[0], nn[1], nn[2], nn[3], nn[4], nn[5], nn[6], nn[7]);
      
      int constexpr node_count = 0; // 1 or 0 switch
      view2D<double> rphi(node_count*9, align<2>(rg[SMT]->n));
      std::vector<double> rtru(node_count*rg[TRU]->n);

      view2D<double> rprj(nln, stride); // mr might be much larger than needed since mr is taken from the TRU grid
      // preparation for the projector functions
      stat += expand_sho_projectors(rprj.data(), rprj.stride(), *rg[SMT], sigma, numax, 1, 0);
      
      ir_stop[SMT] = std::min(radial_grid::find_grid_index(*rg[SMT], Rlog_over_sigma*sigma), rg[SMT]->n - 2);
      double const Rlog = rg[SMT]->r[ir_stop[SMT]];
      if (echo > 0) printf("# %s %s check at radius %g %s\n", label, __func__, Rlog*Ang, _Ang);
      ir_stop[TRU] = ir_stop[SMT] + nr_diff;

      view3D<float> gncs(1 + nen, 1 + lmax, TRU_AND_SMT);
      for(int ien = 0; ien <= nen; ++ien) {
          auto const energy = energy_range[0] + ien*dE;
//        if (echo > 0) printf("# node-count at %.6f %s", energy*eV, _eV);

          if (echo > 22) printf("%.6f", energy*eV);
          for(int ell = 0; ell <= lmax; ++ell) 
          { // ell-loop
              int const iln_off = sho_tools::ln_index(numax, ell, 0);
              for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                  double const gnc = (TRU == ts) ?
                     generalized_node_count_TRU(*rg[ts], rV[ts], ell, energy, gg.data(), ff.data(), ir_stop[ts], echo) :
                     generalized_node_count_SMT(*rg[ts], rV[ts], ell, energy, gg.data(), ff.data(), ir_stop[ts],
                                                view2D<double>((iln_off < nln)?rprj[iln_off]:nullptr, rprj.stride()), nn[ell],
                                                &aHm[iln_off*nln + iln_off],
                                                &aSm[iln_off*nln + iln_off], nln, echo);
                  gncs(ien,ell,ts) = gnc;
                  if (echo > 22) printf("%c%.6f", (ts)?' ':'\t', gnc);
              } // ts
//            if (echo > 0) printf("# %i %g %g %g %g\n", ell, dg[TRU], vg[TRU], dg[SMT], vg[SMT]);

          } // ell
          if (echo > 22) printf("\n");
      } // ien

      // analyze and show a compressed summary
      if (echo > 1) { // show summary
          if (dE > 0) {
              if (echo > 2) printf("\n# %s logder at %.3f %s, summary for %.3f to %.3f %s, dE= %.1f m%s\n",
                  label, Rlog*Ang, _Ang, energy_range[0]*eV, energy_range[2]*eV, _eV, dE*1000*eV, _eV);
              int constexpr mres = 8;
              double at_energy[TRU_AND_SMT][mres];
              double max_abs_diff{0}, max_diff{0}; int ell_max_diff{-1}; // check the lowest resonances only
              for(int ell = 0; ell <= lmax; ++ell) {
                  char more{0}; // will be set to '+' if there are more than mres
                  int nres[TRU_AND_SMT] = {0, 0}; // number of resonances stored
                  for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                      set(at_energy[ts], mres, 0.0); // clear
                      float prev_gnc{-9};
                      for(int ien = 1; ien <= nen; ++ien) {
                          float const gnc = gncs(ien,ell,ts);
                          if (gnc < prev_gnc) { // possible resonance
                              auto const energy = energy_range[0] + (ien - .5)*dE;
                              if (nres[ts] < mres) { at_energy[ts][nres[ts]++] = energy; } else { more = '+'; }
                          } // not increasing
                          prev_gnc = gnc;
                      } // ien
                      int const nshow = nres[ts];
                      if (nshow > 0 && echo > 3) {
                          printf("# %s %c-%s resonances at", label, ellchar[ell], ts?"smt":"tru");
                          for(int ires = 0; ires < nshow; ++ires) {
                              printf("\t%.3f", at_energy[ts][ires]*eV);
                          } // ires
                          printf(" %s\n", _eV);
                      } // echo
                  } // ts
                  double const res_diff = at_energy[SMT][0] - at_energy[TRU][0]; // check the lowest resonance only
                  if (std::abs(res_diff) > max_abs_diff) {
                      ell_max_diff = ell;
                      max_diff = res_diff;
                      max_abs_diff = std::abs(res_diff);
                  } // new maximum found
                  int const ndiff = std::min(nres[TRU], nres[SMT]);
                  if (ndiff > 0 && echo > 2) {
                      printf("# %s %c-resonance differences", label, ellchar[ell]);
                      for(int ires = 0; ires < ndiff; ++ires) {
                          printf("\t%.1f", (at_energy[SMT][ires] - at_energy[TRU][ires])*1000*eV);
                      } // ires
                      printf("%s m%s\n", more?" ...":"", _eV); // should come out as " meV" or " mHa" or " mRy"
                      fflush(stdout);
                  } // echo
              } // ell

              if (ell_max_diff >= 0) { // now the most useful summary line:
                  printf("# %s absolute largest logder difference is %g %s = %d * %g m%s found in %c-channel\n\n",
                      label, max_diff*eV, _eV, int(std::round(max_diff/dE)), dE*1000*eV, _eV, ellchar[ell_max_diff]);
              } // ell valid

          } else { // dE > 0
              printf("# %s logarithmic_derivative summary needs a positive logder.step, found %g %s\n\n", label, dE*eV, _eV);
          } // dE > 0
      } // echo

      if (echo > 5) { // display logarithmic derivatives or generalized node counts
          printf("\n\n## %s logarithmic_derivative from %.3f to %.3f in %i steps of %g %s\n", 
                label, energy_range[0]*eV, (energy_range[0] + nen*dE)*eV, nen + 1, dE*eV, _eV);
          for(int ien = 0; ien <= nen; ++ien) {
              auto const energy = energy_range[0] + ien*dE;
              printf("%.6f", energy*eV);
              for(int ell = 0; ell <= lmax; ++ell) { // ell-loop
                  printf("\t%.6f %.6f", gncs(ien,ell,TRU), gncs(ien,ell,SMT));
              } // ell
              printf("\n");
          } // ien
          printf("\n\n");
      } // echo
      
      return stat;
  } // logarithmic_derivative
  
  
  
  
  
  inline status_t eigenstate_analysis(radial_grid_t const& gV // grid descriptor for Vsmt
              , double const Vsmt[] // smooth potential given on radial grid
              , double const sigma // sigma spread of SHO projectors
              , int const lmax // ellmax up to which the analysis should go
              , uint8_t const nn_input[] // number of projectors per ell
              , int const numax
              , double const aHm[] // non-local Hamiltonian elements in ln_basis, assume stride nln
              , double const aSm[] // non-local overlap matrix elements, assume stride nln
              , int const nr=384 // number of radial grid points in equidistance mesh
              , double const Vshift=0 // potential shift
              , char const *label=""
              , int const echo=2
    ) {
      status_t stat(0);
      auto const g = *radial_grid::create_equidistant_radial_grid(nr + 1, gV.rmax);
      auto const dr = g.dr[0]; // in an equidistant grid, the grid spacing is constant and, hence, indepent of ir
      if (echo > 1) printf("\n# %s %s %s dr=%g nr=%i rmax=%g %s\n", label, __FILE__, __func__, dr*Ang, nr, dr*nr*Ang, _Ang); 
      
      int nn[12]; for(int ell = 0; ell < 12; ++ell) { nn[ell] = (numax + 2 - ell)/2; } // == nnmax(numax, ell)
      
      auto Vloc = std::vector<double>(g.n);
      { // scope: interpolate to the equidistant grid by Bessel-transform
          int const nq = nr/2; double const dq = .125;
          auto Vq = std::vector<double>(nq);
          stat += bessel_transform::transform_s_function(Vq.data(), Vsmt, gV, nq, dq, false, 0);
          stat += bessel_transform::transform_s_function(Vloc.data(), Vq.data(), g, nq, dq, true, 0); // back=true
          for(int ir = 0; ir < g.n; ++ir) {
              Vloc[ir] += Vshift;
          } // ir
          if (echo > 8) {
              printf("\n## Vsmt:\n"); for(int ir = 1; ir < gV.n; ++ir) printf("%g %g\n", gV.r[ir], Vsmt[ir]); printf("\n\n");
              printf("\n## Vq:\n");   for(int iq = 0; iq < nq; ++iq) printf("%g %g\n", iq*dq, Vq[iq]); printf("\n\n");
              printf("\n## Vloc:\n"); for(int ir = 1; ir < nr; ++ir) printf("%g %g\n", g.r[ir], Vloc[ir]); printf("\n\n");
          } // echo
      } // scope

      int const stride = align<2>(nr + 1); assert(stride >= nr);
      // allocate Hamiltonian and overlap matrix
      view2D<double> Ham(nr, stride); // get memory
      view2D<double> Ovl(nr, stride); // get memory
      view2D<double> Ovl_copy(nr, stride); // get memory
      std::vector<double> eigs(nr); // get memory

      int const nln = sho_tools::nSHO_radial(numax);

      // generate normalized SHO projectors
      view2D<double> rprj(nln, stride); // projector functions*r
      // preparation for the projector functions
      stat += expand_sho_projectors(rprj.data(), rprj.stride(), g, sigma, numax, 1, echo);

      int constexpr nFD = 4; double cFD[1 + nFD]; set(cFD, 1 + nFD, 0.0);
      stat += (nFD != finite_difference::set_Laplacian_coefficients(cFD, nFD, dr, 'r'));
      if (echo > 3) printf("# %s %s finite difference with %i neighbors\n", __FILE__, __func__, nFD); 

      for(int ell = 0; ell <= lmax; ++ell) {
          int const ln_off = sho_tools::ln_index(numax, ell, 0);

          set(Ham, nr, 0.0); // clear Hamiltonian
          set(Ovl, nr, 0.0); // clear overlap matrix
          view2D<double const> aHm_ell(&aHm[ln_off*nln + ln_off], nln);
          view2D<double const> aSm_ell(&aSm[ln_off*nln + ln_off], nln);

          // setup the local Hamiltonian
          for(int ir = 0; ir < nr; ++ir) {
              double const r = g.r[ir + 1];
              // diagonal: local potential and repulsive angular part of the kinetic energy
              Ham(ir,ir) = Vloc[ir + 1] + 0.5*(ell*(ell + 1.)/pow2(r)); 
              Ovl(ir,ir) = 1.;
              for(int jr = 0; jr < nr; ++jr) {
                  int const dij = std::abs(ir - jr);
                  if (dij <= nFD) Ham(ir,jr) -= 0.5*cFD[dij]; // finite_difference kinetic energy operator
              } // jr
              // in the radial representation, the usual Laplacian d^2/dr^2 can be applied 
              // for the radial component if the wave functions are in r*phi(r) representation
          } // ir
          
          view2D<double const> rprj1((ln_off < nln)?(rprj[ln_off] + 1):nullptr, rprj.stride()); 
          // forward the rprj-pointer by one so that ir=0 will access the first non-zero radius

          // add the non-local dyadic operators to the Hamiltonian and overlap
          double projector_norm[8] = {0,0,0,0,0,0,0,0};
          for(int ir = 0; ir < nr; ++ir) {
              for(int jr = 0; jr < nr; ++jr) {
                  for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                      for(int mrn = 0; mrn < nn[ell]; ++mrn) {
                          Ham(ir,jr) += rprj1(nrn,ir) * aHm_ell(nrn,mrn) * rprj1(mrn,jr) * dr;
                          Ovl(ir,jr) += rprj1(nrn,ir) * aSm_ell(nrn,mrn) * rprj1(mrn,jr) * dr;
                      } // mrn
                  } // nrn
              } // jr
              for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                  projector_norm[nrn] += pow2(rprj1(nrn,ir)) * dr; 
              } // nrn
          } // ir

#ifdef DEVEL
          if (echo > 19) { // debug
              printf("# %s: projector norm of ell=%i is ", __func__, ell);
              for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                  printf(" %g", projector_norm[nrn]);
              }   printf("\n");
          } // echo

          if (nn[ell] > 0 && echo > 9) { // debug
              printf("# %s: charge deficits for ell=%i are ", __func__, ell);
              for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                  for(int mrn = nrn; mrn < nn[ell]; ++mrn) { // triangular loop
                      printf(" %g", aSm_ell(nrn,mrn));
                  } // mrn
              }   printf("\n");
          } // echo
#endif
          set(Ovl_copy.data(), nr*stride, Ovl.data()); // copy

          { // scope: diagonalize
              bool diagonalize_overlap_matrix = true; // true:always, false:if GEP fails (GEP=generalized eigenvalue problem)
              // solve the generalized eigenvalue problem
              auto const info = linear_algebra::eigenvalues(eigs.data(), nr, Ham.data(), Ham.stride(), Ovl.data(), Ovl.stride());
              if (0 == int(info)) {

                  int const nev = 5 - ell/2; // show less eigenvalues for higher ell-states
                  if (echo > 1) {
                      printf("# %s lowest %c-eigenvalues ", label, ellchar[ell]);
                      for(int iev = 0; iev < nev; ++iev) {
                          printf("  %.6f", eigs[iev]*eV);
                      }   printf("  %s\n", _eV);
                  } // echo
                  
                  if (echo > 2) {
                      // projection analysis for the lowest nev eigenvectors
                      view2D<double> evec(Ham.data(), Ham.stride()); // eigenvectors are stored in the memory space that was used to input the Hamiltonian
                      for(int iev = 0; iev < nev; ++iev) {
                          // plot eigenvectors
                          if (echo > 9) {
                              printf("\n## %s %s ell=%i eigenvalue%10.6f %s %i-th eigenvector:\n", label, __func__, ell, eigs[iev]*eV, _eV, iev);
                              for(int ir = 0; ir < nr; ++ir) {
                                  printf("%g %g\n", g.r[ir + 1], evec(iev,ir));
                              }   printf("\n\n");
                          } // echo

                          if (echo > 7) {
                              printf("# %s projection analysis for ell=%i eigenvalue (#%i)%10.6f %s  coefficients", label, ell, iev, eigs[iev]*eV, _eV);
                              for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                                  printf("%12.6f", dot_product(nr, evec[iev], rprj1[nrn])*std::sqrt(dr));
                              }   printf("\n");
                          } // echo
                      } // iev
                  } // echo
                  
              } else { // info
                  if (echo > 2) printf("# %s diagonalization for ell=%i returned info=%i\n", label, ell, int(info));
                  diagonalize_overlap_matrix = true;
                  ++stat;
              }

              if (diagonalize_overlap_matrix) {
                  // diagonalize Ovl_copy, standard eigenvalue problem
                  linear_algebra::eigenvalues(eigs.data(), nr, Ovl_copy.data(), Ovl_copy.stride());
                  if (eigs[0] <= 0) { // warn
                      if (echo > 0) printf("# %s %s ell=%i lowest eigenvalue of the overlap matrix is non-positive! %g\n", label, __func__, ell, eigs[0]);
                  } // overlap matrix is not positive definite
                  if (echo > 8) {
                      printf("# %s %s ell=%i eigenvalues of the overlap matrix", label, __func__, ell);
                      for(int iev = 0; iev < 8; ++iev) {
                          printf(" %g", eigs[iev]);
                      }   printf("\n");
                  } // echo
              } // info
          } // scope: diagonalize

      } // ell

      return stat;
  } // eigenstate_analysis
  
  
  template<typename real_t>
  status_t emm_average(real_t Mln[], real_t const Mlmn[], int const numax, uint8_t const nn[], int const stride=-1) {
      int const echo = 1;
      if (echo > 4) printf("# %s: numax = %i\n", __func__, numax);
      int const nln = sho_tools::nSHO_radial(numax);
      int const nlmn = sho_tools::nSHO(numax);
      int const M_stride = (stride < 0) ? nlmn : stride;
      auto ln_list = std::vector<int>(nlmn, 0);
      auto lf_list = std::vector<real_t>(nln, 0);
      { // scope: fill the ln_list
          
          if (echo > 6) printf("# %s: ln_list ", __func__);
          for(int ell = 0; ell <= numax; ++ell) {
              for(int emm = -ell; emm <= ell; ++emm) {
                  for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                      int const ilmn = sho_tools::lmn_index(numax, ell, emm, nrn);
                      int const iln = sho_tools::ln_index(numax, ell, nrn);
                      ln_list[ilmn] = iln;
                      if (echo > 6) printf(" %i", ln_list[ilmn]);
                  } // nrn
              } // emm
              for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                  int const iln = sho_tools::ln_index(numax, ell, nrn);
                  lf_list[iln] = real_t(1)/std::sqrt(2*ell + real_t(1));
              } // nrn
          } // ell
          if (echo > 6) printf("\n");
          if (echo > 5) {
              printf("# %s: l2p1_list ", __func__);
              for(int iln = 0; iln < nln; ++iln) {
                  printf(" %.1f", (lf_list[iln] > 0) ? 1./pow2(lf_list[iln]) : 0);
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
  inline status_t all_tests(int const echo=0) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS
  
  inline status_t test_eigenstate_analysis(int const echo=3, int const lmax=7) {
      if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
      // test the eigenstate analysis with a harmonic potential with projectors but zero non-local matrices
      auto const rg = *radial_grid::create_default_radial_grid(0);
      int const nln = sho_tools::nSHO_radial(lmax);
      std::vector<double> V(rg.n, 0.0);
      double const sigma = 1.0; // if the rmax ~= 10, lmax = 7, sigma <= 1.5, otherwise projectors leak out
      for(int ir = 0; ir < rg.n; ++ir) V[ir] = 0.5*(pow2(rg.r[ir]) - pow2(rg.rmax))/pow4(sigma); // harmonic potential 
      std::vector<double> const aHm(nln*nln, 0.0);
      std::vector<uint8_t> nn(1 + lmax, 0); for(int ell = 0; ell <= lmax; ++ell) nn[ell] = 1 + (lmax - ell)/2;
      return eigenstate_analysis(rg, V.data(), sigma, lmax, nn.data(), lmax, aHm.data(), aHm.data(), 128, 0.0, "", echo);
      // expected result: eigenstates at E0 + (2*nrn + ell)*sigma^-2 Hartree
  } // test_eigenstate_analysis

  inline status_t all_tests(int const echo=0) {
      if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
      status_t status(0);
      status += test_eigenstate_analysis(echo);
      return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace scattering_test
