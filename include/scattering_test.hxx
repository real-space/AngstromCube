#pragma once

#include <cstdio> // std::printf, std::fflush, stdout
#include <cassert> // assert
#include <cstdint> // uint8_t
#include <vector> // std::vector<T>
#include <cmath> // std::copysign, ::atan2, ::sqrt, ::ceil, ::abs, ::exp, ::round
#include <algorithm> // std::min, std::max

#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // ::create_radial_grid, ::equation_equidistant,
                           // ::create_default_radial_grid, ::find_grid_index
                           // ::destroy_radial_grid
#include "bessel_transform.hxx" // ::transform_s_function
#include "finite_difference.hxx" // ::set_Laplacian_coefficients
#include "sho_radial.hxx" // ::radial_eigenstates, ::radial_normalization, ::expand_poly
#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set, add_product, dot_product, align<nBits>, product, pow4
#include "sho_tools.hxx" // ::nSHO, ::nSHO_radial
#include "radial_integrator.hxx" // ::integrate_outwards<SRA>
#include "constants.hxx" // ::pi
#include "linear_algebra.hxx" // ::linear_solve, ::eigenvalues
#include "data_view.hxx" // view2D<T>
#include "print_tools.hxx" // printf_vector
#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

namespace scattering_test {

  int constexpr TRU=0, SMT=1, TRU_AND_SMT=2;

  char const ellchar[] = "spdfghijkl0123456789"; // ToDo: by convention 'i' is not the proper char for ell=6

  template <typename real_t> inline
  real_t arcus_tangent(real_t const nom, real_t const den) {
      return (den < 0) ? std::atan2(-nom, -den) : std::atan2(nom, den);
  } // arcus_tangent

  template <typename real_t>
  size_t count_nodes(real_t const f[], int const n) {
      size_t nnodes{0};
      for (int i = 1; i < n; ++i) {
          nnodes += ((f[i - 1]*f[i]) < 0);
      } // i
      return nnodes;
  } // count_nodes

 
  template <typename real_t>
  status_t expand_sho_projectors(
        real_t prj[]  // output projectors [nln*stride]
      , int const stride // stride >= rg.n
      , radial_grid_t const & rg // radial grid descriptor
      , double const sigma // SHO basis spread
      , int const numax // SHO basis size
      , int const rpow=0 // return r^rpow * p(r)
      , int const echo=0 // log-level
      , real_t dprj[]=nullptr // derivative of projectors w.r.t. sigma
  ) {
      assert(sigma > 0);
      double const sigma_inv = 1./sigma;
      double const sigma_m23 = std::sqrt(pow3(sigma_inv)); // == sigma^{-3/2}
      double const dds_sigma_m23 = -1.5*sigma_inv; // d/dsigma sigma^{-3/2} = -3/2 sigma^{-5/2} =  -3/2 * 1./sigma * sigma^{-3/2} 
      int const maxpoly = align<2>(1 + numax/2);
      int const nln = sho_tools::nSHO_radial(numax);
      view2D<double> poly(nln, maxpoly, 0.0);
      view2D<double> ddx2_poly(nln, maxpoly, 0.0);
      for (int ell = 0; ell <= numax; ++ell) {
          for (int nrn = 0; nrn < sho_tools::nn_max(numax, ell); ++nrn) { // number of radial nodes
              int const iln = sho_tools::ln_index(numax, ell, nrn);
              sho_radial::radial_eigenstates(poly[iln], nrn, ell); //  a polynomial in x^2
              double const norm_factor = sho_radial::radial_normalization(poly[iln], nrn, ell) * sigma_m23;
              scale(poly[iln], nrn + 1, norm_factor);
              for (int p = 1; p <= nrn; ++p) {
                  ddx2_poly(iln,p - 1) = p*poly(iln,p); // derive polynomial w.r.t. its argument x^2
              } // p
              // Beware for the derivative: the norm_factor depends on sigma
          } // nrn
      } // ell

      assert(rg.n <= stride);
#ifdef DEVEL
      std::vector<double> norm(nln, 0.0);
#endif // DEVEL
      if (echo > 18) std::printf("\n## SHO projectors on radial grid: r, p_00(r), p_01, ... :\n");
      for (int ir = 0; ir < rg.n; ++ir) { // parallel over ir
          double const r = rg.r[ir];
#ifdef DEVEL
          double const r2dr = rg.r2dr[ir]; // load grid from radial grid descriptor
#endif // DEVEL
//        double const dr = 0.03125, r = dr*ir, r2dr = r*r*dr; // use an equidistant grid
          double const r_pow_rpow = intpow(r, rpow);

          double const x = r*sigma_inv; // x == r/sigma
          double const x2 = x*x; // x^2
          double const Gaussian = (x2 < 160) ? std::exp(-0.5*x2) : 0;

          // prepare for the derivative w.r.t. sigma
          double const ddx_Gaussian = -x*Gaussian; // d/dx Gaussian(x)
          double const ddx_x2 = 2*x; // d/dx x^2 = 2x
          double const dds_x = -x*sigma_inv; // d/dsigma x = d/dsigma r/sigma = -r/sigma^2 = -x/sigma

          if (echo > 18) std::printf("%g ", r);
          double x_pow_ell{1}; // x^ell
          double ddx_x_pow_ell{0}; // ell*x^(ell - 1)
          for (int ell = 0; ell <= numax; ++ell) { // serial loop, must run forward 
              for (int nrn = 0; nrn < sho_tools::nn_max(numax, ell); ++nrn) {
                  int const iln = sho_tools::ln_index(numax, ell, nrn);
                  double const poly_value = sho_radial::expand_poly(poly[iln], 1 + nrn, x2);
                  double const projector_value = poly_value * Gaussian * x_pow_ell;
#ifdef DEVEL
                  if (echo > 18) std::printf(" %g", projector_value);
                  norm[iln] += pow2(projector_value) * r2dr;
#endif // DEVEL
                  prj[iln*stride + ir] = projector_value * r_pow_rpow; // store

                  if (nullptr != dprj) {
                      // construct the derivative w.r.t. sigma
                      double const ddx2_poly_value = sho_radial::expand_poly(ddx2_poly[iln], nrn, x2);
                      double const ddx_poly_value = ddx_x2 * ddx2_poly_value;
                      double const dds_projector_value =
                                   dds_sigma_m23 * projector_value + // derivate of the prefactor is an additional factor here
                                   dds_x * (ddx_poly_value * Gaussian * x_pow_ell +
                                            poly_value * ddx_Gaussian * x_pow_ell +
                                            poly_value * Gaussian * ddx_x_pow_ell);

                      dprj[iln*stride + ir] = dds_projector_value * r_pow_rpow; // store derivative
                  } // compute also the derivative w.r.t. sigma
              } // nrn

              // update x^ell and its derivative for the next ell-iteration
              ddx_x_pow_ell = (ell + 1)*x_pow_ell;
              x_pow_ell *= x;
          } // ell
          if (echo > 18) std::printf("\n");
      } // ir
      if (echo > 18) std::printf("\n\n");
#ifdef DEVEL
      if (echo > 9) {
          std::printf("# projector normalizations are ");
          printf_vector(" %g", norm.data(), nln);
      } // echo
#endif // DEVEL
      return 0;
  } // expand_sho_projectors


  inline double find_outwards_solution( // returns derivative at maximum
        radial_grid_t const & rg // radial grid descriptor
      , double const rV[] // effective potential r*V(r)
      , int const ell // angular moment quantum number
      , double const energy // energy parameter
      , double gg[] // work arrays, greater component
      , double ff[] // work arrays, smaller component
      , int const ir_stop=-1 // radius where to stop
      , double const *inh=nullptr // r*inhomogeneiety
  ) {
      int constexpr SRA = 1; // use the scalar relativistic approximation
      double deriv{0};
      radial_integrator::integrate_outwards<SRA>(gg, ff, rg, rV, ell, energy, ir_stop, &deriv, inh);
      return deriv;
  } // find_outwards_solution

  template <bool NodeCount=false> inline
  double generalized_node_count_TRU(
        radial_grid_t const & rg // radial grid descriptor
      , double const rV[] // effective potential r*V(r)
      , int const ell // angular moment quantum number
      , double const energy // energy parameter
      , double gg[] // work arrays, greater component
      , double ff[] // work arrays, smaller component
      , int const ir_stop // radius where to stop
      , int const echo=0 // log-level
  ) {
      if (echo > 19) std::printf("# find true homogeneous solution for ell=%i E=%g %s\n", ell, energy*eV,_eV); // DEBUG
      double const deriv = find_outwards_solution(rg, rV, ell, energy, gg, ff, ir_stop, nullptr);
      double const value = gg[ir_stop]; // value of the greater component at Rlog
      int nnodes{0}; if (NodeCount) nnodes = count_nodes(gg, ir_stop + 1);
      double constexpr one_over_pi = 1./constants::pi;
      return (nnodes + 0.5 - one_over_pi*arcus_tangent(deriv, value));
  } // generalized_node_count_TRU

  template <bool NodeCount=false, int const nLIM=9> inline
  double generalized_node_count_SMT(
        radial_grid_t const & rg // radial grid descriptor
      , double const rV[] // effective potential r*V(r)
      , int const ell // angular moment quantum number
      , double const energy // energy parameter
      , double gg[] // work arrays, greater component
      , double ff[] // work arrays, smaller component
      , int const ir_stop // radius where to stop
      , view2D<double> const & rprj // pointer to inhomogeneieties
      , int const n=0 // number of projector >= 0
      , double const *aHm=nullptr // pointer to Hamiltonian, can be zero if n=0
      , double const *aSm=nullptr // pointer to overlap, can be zero if n=0
      , int const stride=0 // stride for Hamiltonian and overlap
      , int const echo=0 // log-level
  ) {
      double deriv[nLIM], value[nLIM], gfp[nLIM*nLIM]; assert(n < nLIM);

      view2D<double> waves(NodeCount ? (1 + n) : 0, align<2>(ir_stop + 1)); // get memory

      for (int jrn = n; jrn >= 0; --jrn) { // 0:homogeneous, >0:inhomogeneous
          if (echo > 19) std::printf("# find smooth %shomogeneous solution for ell=%i E=%g %s (%i)\n", (jrn > 0)?"in":"", ell, energy*eV,_eV, jrn-1);
          auto const ggj = NodeCount ? waves[jrn] : gg;
          deriv[jrn] = find_outwards_solution(rg, rV, ell, energy, ggj, ff, ir_stop, (jrn > 0) ? rprj[jrn - 1] : nullptr);
          value[jrn] = ggj[ir_stop]; // value of the greater component at Rlog
          for (int krn = 0; krn < n; ++krn) {
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
          for (int i = 0; i < n0; ++i) {
              mat[i*n0 + i] = 1.0; // unit matrix
          } // i
          for (int irn = 0; irn < n; ++irn) {
              for (int jrn = 0; jrn < n0; ++jrn) {
                  for (int krn = 0; krn < n; ++krn) {
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
          for (int jrn = 0; jrn <= n; ++jrn) {
              add_product(gg, ir_stop + 1, waves[jrn], x[jrn]);
          } // jrn
          nnodes = count_nodes(gg, ir_stop + 1); 
      } // NodeCount
      double constexpr one_over_pi = 1./constants::pi;
      return (nnodes + 0.5 - one_over_pi*arcus_tangent(der, val));
  } // generalized_node_count_SMT

  inline status_t logarithmic_derivative(
        radial_grid_t const rg[TRU_AND_SMT] // radial grid descriptors for Vtru, Vsmt
      , double        const *const rV[TRU_AND_SMT] // true and smooth potential given on the radial grid *r
      , double const sigma // sigma spread of SHO projectors
      , int const lmax // ellmax up to which the analysis should go
      , int const numax
      , double const aHm[] // non-local Hamiltonian elements in ln_basis
      , double const aSm[] // non-local overlap matrix elements
      , double const energy_range[3] // {lower, step, upper}
      , char const *label="" // log-prefix
      , int const echo=0 // log-level
      , float const Rlog_over_sigma=6.f
  ) {
      status_t stat(0);
      
      double const dE = std::copysign(std::max(1e-9, std::abs(energy_range[1])), energy_range[1]);
      int const nen = int(std::ceil((energy_range[2] - energy_range[0])/dE));

      if (echo > 1) std::printf("# %s %s energy range from %g to %g %s in %d steps of %g %s, lmax=%d%c\n", 
          label, __func__, energy_range[0]*eV, energy_range[2]*eV, _eV, 1 + nen, dE*eV, _eV, lmax, (nen < 0)?'\n':' ');

      if (nen < 0) return stat; // empty range

      int const nr_diff = rg[TRU].n - rg[SMT].n; assert(nr_diff >= 0);
      int const mr = align<2>(rg[TRU].n);
      std::vector<double> gg(mr), ff(mr); // greater and smaller component, TRU grid
      int ir_stop[TRU_AND_SMT];

      int const nln = sho_tools::nSHO_radial(numax);
      int const stride = mr;

      int constexpr node_count = 0; // 1 or 0 switch
      view2D<double> rphi(node_count*9, align<2>(rg[SMT].n));
      std::vector<double> rtru(node_count*rg[TRU].n);

      view2D<double> rprj(nln, stride); // mr might be much larger than needed since mr is taken from the TRU grid
      // preparation for the projector functions
      stat += expand_sho_projectors(rprj.data(), rprj.stride(), rg[SMT], sigma, numax, 1, 0);

      ir_stop[SMT] = std::min(radial_grid::find_grid_index(rg[SMT], Rlog_over_sigma*sigma), rg[SMT].n - 2);
      double const Rlog = rg[SMT].r[ir_stop[SMT]];
      if (echo > 1) std::printf("# %s %s check at radius %g %s\n", label, __func__, Rlog*Ang, _Ang);
      ir_stop[TRU] = ir_stop[SMT] + nr_diff;

      view3D<float> gncs(1 + nen, 1 + lmax, TRU_AND_SMT);
      for (int ien = 0; ien <= nen; ++ien) {
          auto const energy = energy_range[0] + ien*dE;
//        if (echo > 0) std::printf("# node-count at %.6f %s", energy*eV, _eV);

          if (echo > 22) std::printf("%.6f", energy*eV);
          for (int ell = 0; ell <= lmax; ++ell) 
          { // ell-loop
              int const nn = (numax + 2 - ell)/2;
              int const iln_off = sho_tools::ln_index(numax, ell, 0);
              for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                  double const gnc = (TRU == ts) ?
                     generalized_node_count_TRU(rg[ts], rV[ts], ell, energy, gg.data(), ff.data(), ir_stop[ts], echo/2) :
                     generalized_node_count_SMT(rg[ts], rV[ts], ell, energy, gg.data(), ff.data(), ir_stop[ts],
                                                view2D<double>((iln_off < nln)?rprj[iln_off]:nullptr, rprj.stride()), nn,
                                                &aHm[iln_off*nln + iln_off],
                                                &aSm[iln_off*nln + iln_off], nln, echo/2);
                  gncs(ien,ell,ts) = gnc;
                  if (echo > 22) std::printf("%c%.6f", (ts)?' ':'\t', gnc);
              } // ts
//            if (echo > 0) std::printf("# %i %g %g %g %g\n", ell, dg[TRU], vg[TRU], dg[SMT], vg[SMT]);

          } // ell
          if (echo > 22) std::printf("\n");
      } // ien

      // analyze and show a compressed summary
      if (echo > 1) { // show summary
          if (dE > 0) {
              if (echo > 2) std::printf("\n# %s logder at %.3f %s, summary for %.3f to %.3f %s, dE= %.1f m%s\n",
                  label, Rlog*Ang, _Ang, energy_range[0]*eV, energy_range[2]*eV, _eV, dE*1000*eV, _eV);
              int constexpr mres = 8;
              double at_energy[TRU_AND_SMT][mres];
              double max_abs_diff{0}, max_diff{0}; int ell_max_diff{-1}; // check the lowest resonances only
              for (int ell = 0; ell <= lmax; ++ell) {
                  char more{0}; // will be set to '+' if there are more than mres
                  int nres[TRU_AND_SMT] = {0, 0}; // number of resonances stored
                  for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                      set(at_energy[ts], mres, 0.0); // clear
                      float prev_gnc{-9};
                      for (int ien = 1; ien <= nen; ++ien) {
                          float const gnc = gncs(ien,ell,ts);
                          if (gnc < prev_gnc) { // possible resonance
                              auto const energy = energy_range[0] + (ien - .5)*dE;
                              if (nres[ts] < mres) { at_energy[ts][nres[ts]++] = energy; } else { more = '+'; }
                          } // not increasing
                          prev_gnc = gnc;
                      } // ien
                      int const nshow = nres[ts];
                      if (nshow > 0 && echo > 3) {
                          std::printf("# %s %c-%s resonances at", label, ellchar[ell], ts?"smt":"tru");
//                           for (int ires = 0; ires < nshow; ++ires) {
//                               std::printf("\t%.3f", at_energy[ts][ires]*eV);
//                           } // ires
                          printf_vector("\t%.3f", at_energy[ts], nshow, "", eV);
                          std::printf(" %s\n", _eV);
                      } // echo
                  } // ts
                  double const res_diff = at_energy[SMT][0] - at_energy[TRU][0]; // check the lowest resonance only
                  if (std::abs(res_diff) > max_abs_diff) {
                      ell_max_diff = ell;
                      max_diff = res_diff;
                      max_abs_diff = std::abs(res_diff);
                  } // new maximum found
                  int const ndiff = std::min(nres[TRU], nres[SMT]);
                  if (echo > 2) {
                      if (ndiff > 0) {
                          std::printf("# %s %c-resonance differences", label, ellchar[ell]);
                          for (int ires = 0; ires < ndiff; ++ires) {
                              std::printf("\t%.1f", (at_energy[SMT][ires] - at_energy[TRU][ires])*1000*eV);
                          } // ires
                          std::printf("%s m%s\n", more?" ...":"", _eV); // should come out as " meV" or " mHa" or " mRy"
                      } else {
                          std::printf("# %s no %c-resonances found below %g %s\n", label, ellchar[ell], energy_range[2]*eV, _eV);
                      } // ndiff > 0
                      std::fflush(stdout);
                  } // echo
              } // ell

              if (ell_max_diff >= 0) { // now the most useful summary line:
                  std::printf("# %s absolute largest logder difference is %g %s = %d * %g m%s found in %c-channel\n",
                      label, max_diff*eV, _eV, int(std::round(max_diff/dE)), dE*1000*eV, _eV, ellchar[ell_max_diff]);
              } // ell valid
              std::printf("\n");

          } else { // dE > 0
              std::printf("# %s logarithmic_derivative summary needs a positive logder.step, found %g %s\n\n", label, dE*eV, _eV);
          } // dE > 0
      } // echo

      if (echo > 7) { // display logarithmic derivatives or generalized node counts
          std::printf("\n\n## %s logarithmic_derivative from %.3f to %.3f in %i steps of %g %s\n", 
                label, energy_range[0]*eV, (energy_range[0] + nen*dE)*eV, nen + 1, dE*eV, _eV);
          for (int ien = 0; ien <= nen; ++ien) {
              auto const energy = energy_range[0] + ien*dE;
              std::printf("%.6f", energy*eV);
              for (int ell = 0; ell <= lmax; ++ell) { // ell-loop
                  std::printf("\t%.6f %.6f", gncs(ien,ell,TRU), gncs(ien,ell,SMT));
              } // ell
              std::printf("\n");
          } // ien
          std::printf("\n\n");
      } // echo
      
      return stat;
  } // logarithmic_derivative
  
  
  
  
  
  inline status_t eigenstate_analysis(
        radial_grid_t const & gV // radial grid descriptor for Vsmt
      , double const Vsmt[] // smooth potential given on radial grid
      , double const sigma // sigma spread of SHO projectors
      , int const lmax // ellmax up to which the analysis should go
      , int const numax // SHO basis size
      , double const aHm[] // non-local Hamiltonian elements in ln_basis, assume stride nln
      , double const aSm[] // non-local overlap matrix elements, assume stride nln
      , int const nr=384 // number of radial grid points in equidistance mesh
      , double const Vshift=0 // potential shift
      , char const *label="" // log-prefix
      , int const echo=0 // log-level
      , float const reference[3][4]=nullptr // up to 3 reference eigenvalues for s,p,d,f
      , double const warning_threshold=3e-3
  ) {
      status_t stat(0);
      auto g = *radial_grid::create_radial_grid(nr + 1, gV.rmax, radial_grid::equation_equidistant);
      auto const dr = g.dr[0]; // in an equidistant grid, the grid spacing is constant and, hence, indepent of ir
      if (echo > 1) std::printf("\n# %s %s nr=%i dr=%g rmax=%g %s\n", label, __func__, nr, dr*Ang, dr*nr*Ang, _Ang); 
//    if (echo > 1) std::printf("# %s %s rmax=%g --> rmax=%g %s\n", label, __func__, gV.rmax*Ang, g.rmax*Ang, _Ang); 

      std::vector<double> Vloc(g.n);
      { // scope: interpolate to the equidistant grid by Bessel-transform
          int const nq = nr/2; double const dq = .125;
          std::vector<double> Vq(nq);
          stat += bessel_transform::transform_s_function(Vq.data(), Vsmt, gV, nq, dq, false, 0);
          stat += bessel_transform::transform_s_function(Vloc.data(), Vq.data(), g, nq, dq, true, 0); // back=true
          for (int ir = 0; ir < g.n; ++ir) {
              Vloc[ir] += Vshift;
          } // ir
          if (echo > 8) {
              std::printf("\n## Vsmt:\n"); for (int ir = 1; ir < gV.n; ++ir) std::printf("%g %g\n", gV.r[ir], Vsmt[ir]); std::printf("\n\n");
              std::printf("\n## Vq:\n");   for (int iq = 0; iq < nq; ++iq) std::printf("%g %g\n", iq*dq, Vq[iq]); std::printf("\n\n");
              std::printf("\n## Vloc:\n"); for (int ir = 1; ir < nr; ++ir) std::printf("%g %g\n", g.r[ir], Vloc[ir]); std::printf("\n\n");
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
      if (echo > 3) std::printf("# %s %s finite difference with %i neighbors\n", label, __func__, nFD); 

      double max_dev_from_reference{-1}; double energy_of_reference{0}; char ellchar_of_reference{'?'};

      for (int ell = 0; ell <= lmax; ++ell) {
          int const nn = (numax + 2 - ell)/2;
          int const ln_off = sho_tools::ln_index(numax, ell, 0);

          set(Ham, nr, 0.0); // clear Hamiltonian
          set(Ovl, nr, 0.0); // clear overlap matrix
          view2D<double const> aHm_ell(&aHm[ln_off*nln + ln_off], nln);
          view2D<double const> aSm_ell(&aSm[ln_off*nln + ln_off], nln);

          // setup the local Hamiltonian
          for (int ir = 0; ir < nr; ++ir) {
              double const r = g.r[ir + 1];
              // diagonal: local potential and repulsive angular part of the kinetic energy
              Ham(ir,ir) = Vloc[ir + 1] + 0.5*(ell*(ell + 1.)/pow2(r)); 
              Ovl(ir,ir) = 1.;
              for (int jr = 0; jr < nr; ++jr) {
                  int const dij = std::abs(ir - jr);
                  if (dij <= nFD) Ham(ir,jr) -= 0.5*cFD[dij]; // finite_difference kinetic energy operator
              } // jr
              // in the radial representation, the usual Laplacian d^2/dr^2 can be applied 
              // for the radial component if the wave functions are in r*phi(r) representation
          } // ir

          view2D<double const> rprj1((ln_off < nln)?(rprj[ln_off] + 1):nullptr, rprj.stride()); 
          // forward the rprj-pointer by one so that ir=0 will access the first non-zero radius

          // add the non-local dyadic operators to the Hamiltonian and overlap
          double projector_norm[8] = {0,0,0,0, 0,0,0,0}; assert(nn <= 8);
          for (int ir = 0; ir < nr; ++ir) {
              for (int jr = 0; jr < nr; ++jr) {
                  for (int nrn = 0; nrn < nn; ++nrn) {
                      for (int mrn = 0; mrn < nn; ++mrn) {
                          Ham(ir,jr) += rprj1(nrn,ir) * aHm_ell(nrn,mrn) * rprj1(mrn,jr) * dr;
                          Ovl(ir,jr) += rprj1(nrn,ir) * aSm_ell(nrn,mrn) * rprj1(mrn,jr) * dr;
                      } // mrn
                  } // nrn
              } // jr
              for (int nrn = 0; nrn < nn; ++nrn) {
                  projector_norm[nrn] += pow2(rprj1(nrn,ir)) * dr; 
              } // nrn
          } // ir

#ifdef DEVEL
          if (echo > 19) { // debug
              std::printf("# %s: projector norm of ell=%i is ", __func__, ell);
              printf_vector(" %g", projector_norm, nn);
          } // echo

          if (nn > 0 && echo > 9) { // debug
              std::printf("# %s: charge deficits for ell=%i are ", __func__, ell);
              for (int nrn = 0; nrn < nn; ++nrn) {
                  for (int mrn = nrn; mrn < nn; ++mrn) { // triangular loop
                      std::printf(" %g", aSm_ell(nrn,mrn));
                  } // mrn
              } // nrn
              std::printf("\n");
          } // echo
#endif // DEVEL
          set(Ovl_copy.data(), nr*stride, Ovl.data()); // copy

          { // scope: diagonalize
              bool diagonalize_overlap_matrix{true}; // true:always, false:if GEP fails (GEP=generalized eigenvalue problem)
              // solve the generalized eigenvalue problem
              auto const info = linear_algebra::eigenvalues(eigs.data(), nr, Ham.data(), Ham.stride(), Ovl.data(), Ovl.stride());
              if (0 == int(info)) {

                  int const nev = 5 - ell/2; // show less eigenvalues for higher ell-states
                  if (echo > 1 && nev > 0) {
                      std::printf("# %s lowest %c-eigenvalues ", label, ellchar[ell]);
                      printf_vector("  %.6f", eigs.data(), nev, "", eV);
                      std::printf("  %s\n", _eV);
                  } // echo

                  if (echo > 2) {
                      // projection analysis for the lowest nev eigenvectors
                      double const sqrt_dr = std::sqrt(dr);
                      view2D<double> evec(Ham.data(), Ham.stride()); // eigenvectors are stored in the memory space that was used to input the Hamiltonian
                      for (int iev = 0; iev < nev; ++iev) {
                          // plot eigenvectors
                          if (echo > 9) {
                              std::printf("\n## %s %s ell=%i eigenvalue%10.6f %s %i-th eigenvector:\n", label, __func__, ell, eigs[iev]*eV, _eV, iev);
                              for (int ir = 0; ir < nr; ++ir) {
                                  std::printf("%g %g\n", g.r[ir + 1], evec(iev,ir));
                              } // ir 
                              std::printf("\n\n");
                          } // echo

                          if (echo > 7) {
                              std::printf("# %s projection analysis for ell=%i eigenvalue (#%i) %10.6f %s  coefficients", label, ell, iev, eigs[iev]*eV, _eV);
                              for (int nrn = 0; nrn < nn; ++nrn) {
                                  std::printf("%12.6f", dot_product(nr, evec[iev], rprj1[nrn])*sqrt_dr);
                              } // nrn
                              std::printf("\n");
                          } // echo
                      } // iev
                  } // echo
                  
                  if (nullptr != reference && ell < 4) {
                      for (int iev = 0; iev < 3; ++iev) {
                          double const eig_ref = reference[iev][ell];
                          if (echo > 21) std::printf("# %s compare %c-eigenvalues %g to %g %s\n",
                                            label, ellchar[ell], eigs[iev]*eV, eig_ref*eV, _eV);
                          if (0 != eig_ref) {
                              auto const dev = eigs[iev] - eig_ref;
                              auto const absdev = std::abs(dev);
                              if (absdev > max_dev_from_reference) {
                                  ellchar_of_reference = ellchar[ell];
                                  energy_of_reference = eig_ref;
                                  max_dev_from_reference = absdev;
                              }
                              if (absdev > warning_threshold) {
                                  if (echo > 2) std::printf("# %s %c-eigenvalue #%i %g deviates %g from %g %s\n",
                                                    label, ellchar[ell], iev, eigs[iev]*eV, dev*eV, eig_ref*eV, _eV);
                                  warn("%s %c-eigenvalue #%i deviates %g m%s", label, ellchar[ell], iev, dev*eV*1e3, _eV);
                              } // deviates
                          } // eig_ref is not exactly zero
                      } // iev
                  } // compare to reference eigenvalues
                  
              } else { // info
                  if (echo > 2) std::printf("# %s diagonalization for ell=%i returned info=%i\n", label, ell, int(info));
                  diagonalize_overlap_matrix = true;
                  ++stat;
              } // info

              if (diagonalize_overlap_matrix) {
                  // diagonalize Ovl_copy, standard eigenvalue problem
                  linear_algebra::eigenvalues(eigs.data(), nr, Ovl_copy.data(), Ovl_copy.stride());
                  if (eigs[0] <= 0) { // warn
                      if (echo > 0) std::printf("# %s %s ell=%i lowest eigenvalue of the overlap matrix is non-positive! %g\n", label, __func__, ell, eigs[0]);
                  } // overlap matrix is not positive definite
                  if (echo > 8) {
                      std::printf("# %s %s ell=%i eigenvalues of the overlap matrix", label, __func__, ell);
                      printf_vector(" %g", eigs.data(), 8);
                  } // echo
              } // diagonalize_overlap_matrix

          } // scope: diagonalize

      } // ell

      if (max_dev_from_reference > 0 && echo > 3) {
          std::printf("# %s %s an eigenvalue deviates %g m%s from its %c-reference %g %s\n",
                          label, __func__, max_dev_from_reference*eV*1e3, _eV, 
                          ellchar_of_reference, energy_of_reference*eV, _eV);
      } // deviates from reference

      // destroy the equidistant radial grid descriptor
      radial_grid::destroy_radial_grid(&g);

      return stat;
  } // eigenstate_analysis


  template <typename real_t>
  status_t emm_average(
        real_t Mln[] // output emm-averaged or emm-summed array[nln*nln]
      , real_t const Mlmn[] // input array[nlmn*nlmn]
      , int const numax // SHO basis size
      , int const avg2sum0=2 // 2:average over emm-values withing each ell-channel, 0:sum only
      , int const echo=0 // log-level
      , int const stride=-1 // optional stride for Mlmn, defaults to nlmn
  ) {
      if (echo > 4) std::printf("# %s: numax = %i\n", __func__, numax);
      int const nln = sho_tools::nSHO_radial(numax);
      int const nlmn = sho_tools::nSHO(numax);
      int const M_stride = (stride < 0) ? nlmn : stride;
      std::vector<int> ln_list(nlmn, 0);
      std::vector<real_t> lf_list(nln, 0);
      // fill the ln_list

      if (echo > 6) std::printf("# %s: ln_list ", __func__);
      for (int ell = 0; ell <= numax; ++ell) {
          for (int emm = -ell; emm <= ell; ++emm) {
              for (int nrn = 0; nrn < sho_tools::nn_max(numax, ell); ++nrn) {
                  int const ilmn = sho_tools::lmn_index(numax, ell, emm, nrn);
                  int const iln = sho_tools::ln_index(numax, ell, nrn);
                  ln_list[ilmn] = iln;
                  if (echo > 6) std::printf(" %i", ln_list[ilmn]);
              } // nrn
          } // emm
          for (int nrn = 0; nrn < sho_tools::nn_max(numax, ell); ++nrn) {
              int const iln = sho_tools::ln_index(numax, ell, nrn);
              lf_list[iln] = real_t(1)/std::sqrt(avg2sum0*ell + real_t(1));
          } // nrn
      } // ell
      if (echo > 6) std::printf("\n");
      if (echo > 5) {
          std::printf("# %s: l2p1_list ", __func__);
          for (int iln = 0; iln < nln; ++iln) {
              std::printf(" %.1f", (lf_list[iln] > 0) ? 1./pow2(lf_list[iln]) : 0);
          } // iln
          std::printf("\n");
      } // echo

      set(Mln, nln*nln, real_t(0)); // clear
      // accumulate
      for (int ilmn = 0; ilmn < nlmn; ++ilmn) {
          int const iln = ln_list[ilmn];
          for (int jlmn = 0; jlmn < nlmn; ++jlmn) {
              int const jln = ln_list[jlmn];
              Mln[iln*nln + jln] += Mlmn[ilmn*M_stride + jlmn];
          } // jlmn
      } // ilmn

      // apply (2*ell + 1) denominators, can be omitted for (0 == avg2sum0)
      for (int iln = 0; iln < nln; ++iln) {
          for (int jln = 0; jln < nln; ++jln) {
              Mln[iln*nln + jln] *= lf_list[iln]*lf_list[jln];
          } // jln
      } // iln
      return 0;
  } // emm_average

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS
  
  inline status_t test_eigenstate_analysis(int const echo=3, int const lmax=7) {
      if (echo > 0) std::printf("\n# %s %s\n", __FILE__, __func__);
      // test the eigenstate analysis with a harmonic potential with projectors but zero non-local matrices
      auto const rg = *radial_grid::create_default_radial_grid();
      int const nln = sho_tools::nSHO_radial(lmax);
      double const sigma = 1.0; // if the rmax ~= 10, lmax = 7, sigma <= 1.5, otherwise projectors leak out
      std::vector<double> const aHm(nln*nln, 0.0); // dummy non-local matrices (constant at zero)
      std::vector<double> V(rg.n, 0.0);
      product(V.data(), rg.n, rg.r, rg.r, 0.5/pow4(sigma)); // harmonic potential
      return eigenstate_analysis(rg, V.data(), sigma, lmax, lmax, aHm.data(), aHm.data(), 128, 0.0, "", echo);
      // expected result: eigenstates at (1.5 + 2*nrn + ell)*sigma^-2 Hartree
      // needs to be checked by human, ToDo: how to export the result?
  } // test_eigenstate_analysis

  inline status_t test_expand_sho_projectors_derivative(int const echo=0
      , int const numax=9
      , double const sigma=1.0
  ) {
      status_t stat(0);
      auto const rg = *radial_grid::create_default_radial_grid();
      int const nr = align<2>(rg.n);
      int const nln = sho_tools::nSHO_radial(numax);
      view3D<double> prj(5, nln, nr, 0.0); // get memory for {sigma, sigma+delta, sigma-delta, ...
      // the numerical derivative (sigma+delta - sigma-delta)/(2 delta), and the analytical d/dsigma}

      stat += expand_sho_projectors(prj(0,0), prj.stride(), rg, sigma, numax, 0, echo, prj(4,0)); // derivative into prj[4]
      double constexpr delta = 1e-9;
      stat += expand_sho_projectors(prj(1,0), prj.stride(), rg, sigma*(1 - delta), numax, 0, echo);
      stat += expand_sho_projectors(prj(2,0), prj.stride(), rg, sigma*(1 + delta), numax, 0, echo);

#ifdef DEVEL
      // check how much <sho_ell_irn|sho_ell_jrn> deviates from a unit matrix
      for (int k = 0; k < 3; ++k) { // loop over 3 different sigma-values: sigma, sigma*(1-delta), sigma*(1+delta)
          double max_dev{0};
          for (int ell = 0; ell <= numax; ++ell) { // ell-block diagonal
              for (int irn = 0; irn < sho_tools::nn_max(numax, ell); ++irn) {
                  if (echo > 7) std::printf("# %s %c%d ", __func__, ellchar[ell], irn);
                  int const iln = sho_tools::ln_index(numax, ell, irn);
                  for (int jrn = 0; jrn < sho_tools::nn_max(numax, ell); ++jrn) {
                      int const jln = sho_tools::ln_index(numax, ell, jrn);
                      auto const aij = dot_product(rg.n, prj(k,iln), prj(k,jln), rg.r2dr);
                      auto const dev = aij - (irn == jrn);
                      max_dev = std::max(max_dev, std::abs(dev));
                      if (echo > 7) std::printf(" %8.1e", dev);
                  } // jrn
                  if (echo > 7) std::printf("\n");
              } // irn
          } // ell
          if (echo > 4) std::printf("# %s: largest deviation from unit matrix for numax= %d is %.1e\n", __func__, numax, max_dev);
          stat += (max_dev > 5e-13);
      } // k
#endif // DEVEL

      // construct a finite-difference derivative w.r.t. sigma in set#3
      add_product(prj(3,0), nln*prj.stride(), prj(1,0), -.5/delta);
      add_product(prj(3,0), nln*prj.stride(), prj(2,0),  .5/delta);

      { // scope: check the difference between the analytically derived projectors (#4) and a finite-difference derived set (#3)
          double max_dev{0};
          for (int ell = 0; ell <= numax; ++ell) {
              for (int irn = 0; irn < sho_tools::nn_max(numax, ell); ++irn) {
                  if (echo > 5) std::printf("# %s: %c%d ", __func__, ellchar[ell], irn);
                  int const iln = sho_tools::ln_index(numax, ell, irn);
                  for (int jrn = 0; jrn < sho_tools::nn_max(numax, ell); ++jrn) {
                      int const jln = sho_tools::ln_index(numax, ell, jrn);
                      auto const ana_ij = dot_product(rg.n, prj(0,iln), prj(4,jln), rg.r2dr);
                      auto const num_ij = dot_product(rg.n, prj(0,iln), prj(3,jln), rg.r2dr);
                      if (echo > 5) std::printf(" %11.6f", ana_ij);
//                       if (echo > 5) std::printf(" %11.6f", num_ij);
                      auto const dev = ana_ij - num_ij;
                      max_dev = std::max(max_dev, std::abs(dev));
                      if (echo > 5) std::printf(" %8.1e", dev);
                  } // jrn
                  if (echo > 5) std::printf("\n");
              } // irn
          } // ell
          if (echo > 4) std::printf("# %s: largest deviation |analytical - numerical| is %.1e\n", __func__, max_dev);
          stat += (max_dev > 2e-6);
      } // scope

#ifdef DEVEL
      if (echo > 21) {
          std::printf("\n## %s: plot numerical and analytical derivative:\n", __func__);
          for (int ir = 0; ir < rg.n; ++ir) {
              std::printf("%g", rg.r[ir]);
              for (int irn = 0; irn < nln; ++irn) {
                  std::printf("  %g %g", prj(3,irn,ir), prj(4,irn,ir));
              } // irn
              std::printf("\n");
          } // ir
          std::printf("\n\n");
      } // echo
#endif // DEVEL

      return stat;
  } // test_expand_sho_projectors_derivative

  inline status_t all_tests(int const echo=0) {
      if (echo > 0) std::printf("\n# %s %s\n", __FILE__, __func__);
      status_t stat(0);
      stat += test_eigenstate_analysis(echo);
      stat += test_expand_sho_projectors_derivative(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace scattering_test
