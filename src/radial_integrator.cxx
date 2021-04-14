#include <cstdio> // std::printf, std::fflush, stdout
#include <vector> // std::vector
#include <cstdlib> // std::abs
#include <cmath> // std::sqrt, std::pow
#include <algorithm> // std::min, std::max

#include "radial_integrator.hxx"

#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // ::create_exponential_radial_grid
#include "inline_math.hxx" // sgn, pow2, align<nBits>
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
#include "display_units.h" // eV, _eV, Ang, _Ang

// #define FULL_DEBUG
// #define DEBUG

#ifdef  DEBUG
  #include "debug_output.hxx" // dump_to_file
#endif

#ifdef FULL_DEBUG
  #define full_debug(print) print 
#else
  #define full_debug(print)
#endif

#ifdef DEBUG
  #define debug(print) print 
#else
  #define debug(print)
#endif

namespace radial_integrator {
  // integrates the radial SRA equation inwards and outwards,
  // accepts Z_effective(r) instead of a spherical potential V(r) or r*V(r)

  double relativistic_power_series( // returns exponent
        double const Z // number of protons
      , ell_QN_t const ell // angular momentum quantum number
      , double const V0mE // V(r=0) - E
      , double ps[][2] // large and small component
  ) {
      int constexpr G=0, F=1;
      double constexpr c0 = 137.0359895; // speed of light
      auto const llp1 = ell*(ell + 1.);
      auto const aa = -std::max(Z, 1.)/c0;
      // gf[G:F] = (sum_{i=0}^{10} ps[i][G:F] * r^(i+s))
      double const s = std::sqrt(llp1 + 1 - aa*aa);
      auto const bb = llp1 - aa*aa;
      auto const p21 = V0mE/c0;
      auto const p12 = 2*c0 - p21;
      ps[0][G] = 1.;          // for exponent 0
      ps[0][F] = (1. - s)/aa; // for exponent 0
      full_debug(printf("# %s: Z=%g l=%d V=%g %s coeff G,F  %g,%g", __func__, Z, ell, V0mE*eV, _eV, ps[0][G], ps[0][F]));
      for (auto k = 1; k < 9; ++k) {
          auto const alfa = p12*ps[k - 1][F];
          auto beta =    aa*p21*ps[k - 1][G];
          if (ell > 0) {
            beta -= p12*(aa*ps[k - 1][G] - (k + s)*ps[k - 1][F]);
            if (k > 1) beta -= p12*p21*ps[k - 2][G];
          } // ell > 0
          auto const detinv = 1./(k*(k + 2*s)*aa);
          ps[k][G] = (alfa*(k + s + 1) - beta)*aa*detinv;
          ps[k][F] = (beta*(k + s - 1) - bb*alfa)*detinv;
          full_debug(printf("  %g,%g", ps[k][G], ps[k][F]));
      } // k
      full_debug(printf("\n"));
      for (auto k = 9; k < 11; ++k) {
          ps[k][G] = 0;
          ps[k][F] = 0;
      } // k
      return s;
  } // relativistic_power_series
  

  void power_series_by_Horner_scheme(
        double & gg // result: large component*r
      , double & ff // result: small component*r
      , double const ps[][2] // polynomial coefficients
      , double const rpow // power for radius
      , double const r // radius
  ) {
        int constexpr G=0, F=1; // constants
        gg = ps[10][G];
        ff = ps[10][F];
        for (auto ip = 9; ip >= 0; --ip) {
            gg = gg*r + ps[ip][G];
            ff = ff*r + ps[ip][F];
        } // ip
        auto const r_power = std::pow(r, rpow);
        gg *= r_power;
        ff *= r_power;
  } // power_series_by_Horner_scheme


  //   !! derivation operator for the Schr\""odinger equation
  //   !! The original Schr\""odinger equation
  //   !!      1                l(l+1)
  //   !!   - -- d^2g(r)/dr^2 + ------ + [ V_Hxc(r) - Ze^2/r - E ]*g(r) = 0
  //   !!     2m                2m r^2
  //   !! becomes (by reduction of order):
  //   !! | dg |      | 1/r   m  | | g |
  //   !! |    | = dr |          | |   |
  //   !! | df |      | 2W  -1/r | | f |
  //   !!
  //   !!            l(l+1)
  //   !! where  W = ------ + (V_Hxc(r) - Ze^2/r - E)
  //   !!            2m r^2
  //   !!                   !!! Hartree units, m=1, e^2=1
  //   !!
  //   !! so       dg/dr = g/r + mf   <==>  mf = dg/dr - g/r
  //   !! and      df/dr = 2Wg - f/r
  //   !!
  //   !! then d^2g/dr^2 = dg/dr/r - g/r^2 + mdf/dr
  //   !!                = dg/dr/r - g/r^2 + m[2Wg - f/r]
  //   !!                = dg/dr/r - g/r^2 + 2mWg - [dg/dr - g/r]/r
  //   !!                = 2mWg
  //   !! <==> -1/(2m) d^2g/dr^2 + Wg = 0
  //   !!
  //   !! derivation operator for the
  //   !! scalar relativistic approximation of the Dirac equation
  //   !! Similar to the Schr\""odinger equation (above)
  //   !! | dg |      | 1/r  m(r) | | g |
  //   !! |    | = dr |           | |   |
  //   !! | df |      | 2W   -1/r | | f |
  //   !!
  //   !! but with the relativistic mass
  //   !!                   E-V_Hxc(r)+Ze^2/r
  //   !!  m(r) = m_0 (1 + -----------------)
  //   !!                      2 m_0 c_0^2
  //   !! also affecting the potential term W
  //   !!              l(l+1)
  //   !! where  W = --------- + (V_Hxc(r) - Ze^2/r - E)
  //   !!            2m(r) r^2
  //   !! Hartree atomic units, m_0=1, e^2=1, c_0=1/alpha=137.036

  template <int SRA> // 0:non_relativistic, 1:scalar_relativistic, 2:sqrt_approximated
  void sra(
        double s[2][2] // result see above
      , float const llp1 // ell*(ell + 1)
      , double const rV
      , double const E // energy
      , double const r // radius
      , double const dr // radial grid spacing
  ) {
      double constexpr m0 = 1; // the rest mass of the electron
      double const Ekin = E - rV/r; // = (E - (V_Hxc(r) - e^2*Z/r)) kinetic energy
      
      // double const c0 = 137.0359895; // speed of light
      double constexpr c0m2 = 5.325136192159324e-05; // 1/c0^2
      double const mrel = m0 * (
              (0 == SRA) ?  1 :                 // non-relativistic mass
            ((2 == SRA) ? (1 + Ekin*c0m2*0.5)  // approximation for sqrt(1 + Ekin/(m0 c0^2))
                : std::sqrt(1 + Ekin*c0m2))     // scalar relativistic mass
                                );
      s[0][0] = dr/r;      // 1/r
      s[0][1] = dr*mrel;   // m(r)
      s[1][0] = dr*(llp1/(mrel*r*r) - 2*Ekin); // 2W with m(r)
      s[1][1] = -s[0][0];  // = -dr/r  !! -1/r
//    full_debug(printf("# %s: \t%g  \t%g\n#\t%g  \t%g\n", __func__, s[0][0], s[0][1], s[1][0], s[1][1]));
  } // sra

  template <int SRA>
  int integrate_inwards( // return the number of nodes
        double gg[] // large component*r
      , double ff[] // small component*r
      , radial_grid_t const & g // radial grid descriptor
      , double const rV[] // radial potential r*V_Hxc(r) - e^2*Z
      , ell_QN_t const ell // angular momentum quantum number
      , double const E // energy (eigen-)value in Hartree
      , double const valder[2] // defaults to nullptr
      , double *dg // =nullptr // derivative at end point
      , int *ir_stopped // =nullptr // index at which the inwards integration stopped
      , int const ir_start // start index, -1:(g.n - 1), default is -1
      , int const ir_stop // latest stop index, default is 4
  ) {
      int constexpr Step = -1; // negative
      float const llp1 = ell*(ell + 1.f);
      double dG[3], dF[3]; // Adam''s Moulton multistep derivatives
      double bG, bF;
//    int sign_gg_prev{0};
      double s[2][2];

      int nnodes = 0; // init result

      // evaluate the first 3 derivatives of g, f and the value of g, f
      // upt to the 3rd grid point by power series expansion

      int const ir0 = (g.n + ir_start)%g.n; // init the inward integration
      int ir = ir0;

      {   // init scope
          double vd[2];
          if (nullptr != valder) {
              vd[0] = valder[0];
              vd[1] = valder[1];
          } else {
              vd[0] = 1e-42;
              vd[1] = -std::sqrt(std::max(1., 2*(rV[ir]/g.r[ir] - E))) * vd[0];
          } // valder
          
          sra<SRA>(s, llp1, rV[ir], E, g.r[ir], g.dr[ir]);
          gg[ir] =  vd[0];
          ff[ir] = (vd[1]*g.dr[ir] + s[1][1]*vd[0]) / s[0][1];

          double aG[4], aF[4];
          aG[0] = vd[1] * g.dr[ir];
          aF[0] = s[1][0] * gg[ir] + s[1][1] * ff[ir];

        for (auto i3 = 0; i3 < 3; ++i3) { // three steps
          ir += Step;

#ifdef FULL_DEBUG
          std::printf("# %s: loop ir=%d r= %g %s i3=%d\n", __func__, ir, g.r[ir]*Ang, _Ang, i3);
#endif
          double const dri12 = (3*g.dr[ir + 1] + 6*g.dr[ir] - g.dr[ir - 1])/8.;
          double const  ri12 = (3* g.r[ir + 1] + 6* g.r[ir] -  g.r[ir - 1])/8.;
          double const rVi12 = (3*  rV[ir + 1] + 6*  rV[ir] -   rV[ir - 1])/8.;

          sra<SRA>(s, llp1, rVi12, E, ri12, dri12); // half step

          bG = gg[ir - Step] - 0.5*aG[0];
          bF = ff[ir - Step] - 0.5*aF[0];

          aG[1] = s[0][0] * bG + s[0][1] * bF;
          aF[1] = s[1][0] * bG + s[1][1] * bF;

          bG = gg[ir - Step] - 0.5*aG[1];
          bF = ff[ir - Step] - 0.5*aF[1];

          aG[2] = s[0][0] * bG + s[0][1] * bF;
          aF[2] = s[1][0] * bG + s[1][1] * bF;

          sra<SRA>(s, llp1, rV[ir], E, g.r[ir], g.dr[ir]);

          bG = gg[ir - Step] - 1.0*aG[2];
          bF = ff[ir - Step] - 1.0*aF[2];

          aG[3] = s[0][0] * bG + s[0][1] * bF;
          aF[3] = s[1][0] * bG + s[1][1] * bF;

          bG = gg[ir - Step] - (aG[0] + 2*aG[1] + 2*aG[2] + aG[3])/6.;
          bF = ff[ir - Step] - (aF[0] + 2*aF[1] + 2*aF[2] + aF[3])/6.;

          gg[ir] = bG; ff[ir] = bF;

  //         auto const sign_gg = sgn(gg[ir]);
  //         nnodes += (sign_gg * sign_gg_prev < 0); // count nodes of the greater component
  //         sign_gg_prev = sign_gg; // for the next iteration

          dG[2 - i3] = s[0][0] * bG + s[0][1] * bF;
          dF[2 - i3] = s[1][0] * bG + s[1][1] * bF;

        } // i3

      } // init scope

      if (nullptr != ir_stopped) *ir_stopped = ir_stop;

      // main loop
      for (ir = ir0 - 4; ir >= ir_stop; ir += Step) {

#ifdef  FULL_DEBUG
          if ((ir < ir_stop + 4) || (ir > ir0 - 7))
          std::printf("# %s: loop ir=%d r= %g %s\n", __func__, ir, g.r[ir]*Ang, _Ang);
#endif

          sra<SRA>(s, llp1, rV[ir], E, g.r[ir], g.dr[ir]);

          // b = (24/h*gf(:,ir-h) +19*d(:,1) -5*d(:,2) + d(:,3) ) / 9.
          // Integrator weights for Adams Moulton multistep (3-step)
          bG = ((24*Step)*gg[ir - Step] + 19*dG[0] - 5*dG[1] + dG[2])/9.;
          bF = ((24*Step)*ff[ir - Step] + 19*dF[0] - 5*dF[1] + dF[2])/9.;

          // determinant(8/3h - s)
          double const eight3rds = 8./(3.*Step);
          double const det = (eight3rds - s[0][0])*(eight3rds - s[1][1]) - s[0][1]*s[1][0];
          double const detinv = 1./det;

          // a special property of s: inv(8/3h-s) = (8/3h+s)/det(8/3h-s)
          // gf(:,ir) = (matmul(s, b) + 8./(3.*h)*b)/det
          gg[ir] = ((eight3rds + s[0][0]) * bG + s[0][1] * bF)*detinv;
          ff[ir] = ((eight3rds + s[1][1]) * bF + s[1][0] * bG)*detinv;

  //         auto const sign_gg = sgn(gg[ir]);
  //         nnodes += (sign_gg * sign_gg_prev < 0); // count nodes of the greater component
  //         sign_gg_prev = sign_gg; // for the next iteration

          // shift the derivatives of gg and ff to next grid point
          dG[2] = dG[1]; dG[1] = dG[0];
          dF[2] = dF[1]; dF[1] = dF[0];

          bG = gg[ir];
          bF = ff[ir];

          // compute the new derivative
          dG[0] = s[0][0] * bG + s[0][1] * bF;
          dF[0] = s[1][0] * bG + s[1][1] * bF;

          if (sgn(gg[ir]) * sgn(dG[0]) >= 0) {
              if (nullptr != ir_stopped) *ir_stopped = ir;
#ifdef FULL_DEBUG
              std::printf("# %s: extremum at ir=%d r= %g %s\n", __func__, ir, g.r[ir]*Ang, _Ang);
#endif
              ir = -8; // stop the loop
          } // extremum reached

      } // ir

#ifdef DEBUG
      if (ir >= 0) std::printf("# %s: extremum not found %d \n", __func__, ir);
#endif

      if (nullptr != dg) *dg = dG[0]; // store derivative at the cut radius to evaluate the eigenvalue

      return nnodes;
  } // integrate_inwards

  template <int SRA>
  int integrate_outwards( // return the number of nodes
        double gg[] // result: large component*r
      , double ff[] // result: small component*r
      , radial_grid_t const & g // radial grid descriptor
      , double const rV[] // radial potential r*V_Hxc(r) - e^2*Z
      , ell_QN_t const ell // angular momentum quantum number
      , double const E // energy (eigen-)value in Hartree
      , int const ir_stop // latest stop index, -1:(g.n - 1), default is -1
      , double *dg // =nullptr, derivative at end point
      , double const *rp // =nullptr // optional input: inhomogeneity*r
  ) {
      int constexpr Step = 1; // positive
      float const llp1 = ell*(ell + 1.f);
      double dG[3], dF[3]; // Adam''s Moulton multistep derivatives
      double bG, bF;
      int sign_gg_prev{0};
      double s[2][2];

      int nnodes{0}; // init result
      int ir{0}; // radial grid index

      { // init scope
        double ps[11][2], rpow;
        if (nullptr == rp) {
          rpow = relativistic_power_series(-rV[0], ell, (rV[1] - rV[0])/g.r[1] - E, ps);
        } // rp

        // evaluate the first 3 derivatives of g, f and the value of g, f
        // upt to the 3rd grid point by power series expansion

        gg[ir] = 0; ff[ir] = 0; // r*wave function is zero at r=0

        for (auto i3 = 0; i3 < 3; ++i3) { // three steps
          ir += Step;

          if (nullptr != rp) {
              double const denom = rV[ir] - g.r[ir]*E;
              if (fabs(denom) < 1e-14) exit(__LINE__); // Energy must differ from Potential value!
              bG = rp[ir]*g.r[ir]/denom; // the behaviour of p(ir) is assumed to be p0*r^(ell+1), too
              bF = 0.; // b(G) * m0 * ell / r(ir) ! ToDo: derive properly
          } else {
              power_series_by_Horner_scheme(bG, bF, ps, rpow, g.r[ir]);
              bG = g.r[ir]; // cos
              bF = 0;       // sin
          } // rp

          gg[ir] = bG; ff[ir] = bF;

          auto const sign_gg = sgn(gg[ir]);
          nnodes += (sign_gg * sign_gg_prev < 0); // count nodes of the greater component
          sign_gg_prev = sign_gg; // for the next iteration

          sra<SRA>(s, llp1, rV[ir], E, g.r[ir], g.dr[ir]);

          dG[2 - i3] = s[0][0] * bG + s[0][1] * bF;
          dF[2 - i3] = s[1][0] * bG + s[1][1] * bF;

          if (nullptr != rp) {
              dF[2 - i3] -= 2*rp[ir]*g.dr[ir];  // inhomogeneity comes in here
          } // rp

#ifdef FULL_DEBUG
          std::printf("# %s: loop ir=%d r= %g %s i3=%d\n", __func__, ir, g.r[ir]*Ang, _Ang, i3);
#endif
        } // i3
      } // init scope

      int const ir_max = (-1 == ir_stop)? (g.n - 1) : std::min(g.n - 1, std::abs(ir_stop));

      // main loop
      for (ir = 4; ir <= ir_max; ++ir) {
#ifdef  FULL_DEBUG
          if ((ir < 7) || (ir > ir_max - 3))
          std::printf("# %s: loop ir=%d r= %g %s\n", __func__, ir, g.r[ir]*Ang, _Ang);
#endif

          // b = (24/h*gf(:,ir-h) +19*d(:,1) -5*d(:,2) + d(:,3) ) / 9.
          // Integrator weights for Adams Moulton multistep (3-step)
          bG = ((24*Step)*gg[ir - Step] + 19*dG[0] - 5*dG[1] + dG[2])/9.;
          bF = ((24*Step)*ff[ir - Step] + 19*dF[0] - 5*dF[1] + dF[2])/9.;

          if (nullptr != rp) {
              bF -= (2*8./3.)*g.dr[ir]*rp[ir]; // inhomogeneity comes in here
          } // rp

          sra<SRA>(s, llp1, rV[ir], E, g.r[ir], g.dr[ir]);

          // determinant(8/3h - s)
          double const eight3rds = 8./(3.*Step);
          double const det = (eight3rds - s[0][0])*(eight3rds - s[1][1]) - s[0][1]*s[1][0];
          double const detinv = 1./det;

          // a special property of s: inv(8/3h-s) = (8/3h+s)/det(8/3h-s)
          // gf(:,ir) = (matmul(s, b) + 8./(3.*h)*b)/det
          gg[ir] = ((eight3rds + s[0][0]) * bG + s[0][1] * bF)*detinv;
          ff[ir] = (s[1][0] * bG + (eight3rds + s[1][1]) * bF)*detinv;
          
          bG = gg[ir];
          bF = ff[ir];

          auto const sign_gg = sgn(gg[ir]);
          nnodes += (sign_gg * sign_gg_prev < 0); // count nodes of the greater component
#ifdef  FULL_DEBUG
          if (sign_gg * sign_gg_prev < 0) {
              std::printf("# %s: found node #%d at r[%d]= %g %s for E= %.9f %s\n", __func__, nnodes, ir, g.r[ir]*Ang, _Ang, E*eV, _eV);
          } // node found
#endif
          sign_gg_prev = sign_gg; // for the next iteration

          // shift the derivatives of gg and ff to next grid point
          dG[2] = dG[1]; dG[1] = dG[0];
          dF[2] = dF[1]; dF[1] = dF[0];
          
          // compute the new derivative
          dG[0] = s[0][0] * bG + s[0][1] * bF; // optimization: subexpression of line __LINE__-12
          dF[0] = s[1][0] * bG + s[1][1] * bF;

      } // ir

#ifdef DEBUG
      std::printf("# %s: up to r[%d]= %g %s for E= %.9f %s, %d nodes\n", __func__, ir, g.r[ir]*Ang, _Ang, E*eV, _eV, nnodes);
#endif

      if (dg) *dg = dG[0]; // store derivative at the cut radius to evaluate the eigenvalue

      return nnodes;
  } // integrate_outwards

  template // explicit template instantiation needed for logarithmic derivatives
  int integrate_outwards<1>(double gg[], double ff[],
        radial_grid_t const & g, double const rV[], ell_QN_t const ell, double const E,
        int const ir_stop=-1, double *dg=nullptr, double const *rp=nullptr);

  template <int SRA>
  double shoot_sra( // returns the kink
        radial_grid_t const & g // radial grid descriptor
      , double const rV[] // radial potential r*V_Hxc(r) - e^2*Z
      , ell_QN_t const ell // angular momentum quantum number
      , double const E // energy (eigen-)value in Hartree
      , int & nnodes // number of nodes
      , double* rf=nullptr // radial wave function*r
      , double* r2rho=nullptr // density of that wave function*r^2
  ) {
#ifdef  DEBUG
      std::printf("# %s: l=%d E= %.9f %s\n", __func__, ell, E*eV, _eV);
#endif
      int ir_x{0}; // radial grid index of the extremum
      double dg_inw{0};
      int const N_aligned = align<2>(g.n);
//    std::printf("# %s: N_aligned = %i\n", __func__, N_aligned); std::fflush(stdout); // DEBUG

      std::vector<double>  gf(4*N_aligned);
      auto const gg_inw = &gf[0*N_aligned];
      auto const ff_inw = &gf[1*N_aligned];

      // integrate the Schrodinger equation or SRA equation inwards until
      // the first extremum (maximum/minimum) is reached
      nnodes = integrate_inwards<SRA>(gg_inw, ff_inw, g, rV, ell, E, nullptr, &dg_inw, &ir_x);
#ifdef  DEBUG
      std::printf("# %s: extremum at r[%d]= %g %s\n", __func__, ir_x, g.r[ir_x]*Ang, _Ang);
#endif

      double dg_out{0};
      auto const gg_out = &gf[2*N_aligned];
      auto const ff_out = &gf[3*N_aligned];

      // integrate the Schrodinger equation or SRA equation outwards
      // stopping at the extremum at r(ir)
      nnodes += integrate_outwards<SRA>(gg_out, ff_out, g, rV, ell, E, ir_x, &dg_out);

      // match and scale the tail
      auto const s = gg_out[ir_x] / gg_inw[ir_x];

      if (nullptr != rf) {
          // export the radial wave function*r
          for (auto ir = 0; ir < ir_x; ++ir) {
              rf[ir] = gg_out[ir];
          } // ir
          for (auto ir = ir_x; ir < g.n; ++ir) {
              rf[ir] = gg_inw[ir] * s;
          } // ir
      } // rf

      if (nullptr != r2rho) {
          // create partial density*r^2 of this wave function
          for (auto ir = 0; ir < ir_x; ++ir) {
              r2rho[ir] = pow2(gg_out[ir]);
              // r2rho[ir] += pow2(ff_out[ir]); // relativistic norm correction
          } // ir
          for (auto ir = ir_x; ir < g.n; ++ir) {
              r2rho[ir] = pow2(gg_inw[ir] * s);
              // r2rho[ir] += pow2(ff_inw[ir] * s); // relativistic norm correction
          } // ir
      } // r2rho

      // a step in first derivative of the wave function corresponds to
      // an additional delta-function in the potential with the area:
      double const kink = (dg_out - dg_inw*s) / (gg_out[ir_x] + (0 == gg_out[ir_x]));

#ifdef  DEBUG
      std::printf("# %s: kink of %g at r= %g %s for E= %.9f %s, %d nodes\n", __func__, kink, g.r[ir_x]*Ang, _Ang, E*eV, _eV, nnodes);
#endif

      return kink;
  } // shoot_sra

  double shoot( // returns the kink
        int const sra // 1:scalar relativistic approximation, 0:Schroedinger equation
      , radial_grid_t const & g // radial grid descriptor
      , double const rV[] // radial potential r*V_Hxc(r) - e^2*Z
      , ell_QN_t const ell // angular momentum quantum number
      , double const E // energy (eigen-)value in Hartree
      , int & nn // number of nodes
      , double* rf // [optional, out] radial wave function*r
      , double* r2rho // [optional, out] density of that wave function*r^2
  ) {
    if (0 == sra) return shoot_sra<0>(g, rV, ell, E, nn, rf, r2rho); // non relativistic mass
    if (1 == sra) return shoot_sra<1>(g, rV, ell, E, nn, rf, r2rho); // scalar relativistic mass
    if (2 == sra) return shoot_sra<2>(g, rV, ell, E, nn, rf, r2rho); // approx relativistic mass 
    std::printf("Error in %s, sra=%d not implemented!\n", __func__, sra); exit(__LINE__);
  } // shoot

    
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_hydrogen_atom(
        radial_grid_t const & g // radial grid descriptor
      , double const Z // number of protons in the nucleus
      , int const echo=0 // log-level
  ) {
      // this plots the kink and number of nodes as a function of energy, see doc/fig/20190313_kink_of_energy.*
      std::vector<double> rV(g.n, -Z); // fill all potential values with r*V(r) == -Z
      int nnn_prev{-1};
      for (int iE = 1; iE < 1e6; ++iE) {
          int nnn, ell{0};
          double const E = -.6e-12*pow2(iE*Z);
          auto const kink = shoot(0, g, rV.data(), ell, E, nnn);
          if (nnn_prev != nnn) std::printf("\n");
          if (echo > 3) std::printf("%g %g %d\n", -E, kink, nnn); // shows that the kink is a falling function of E 
          // with poles at the energies where the number of nodes changes, we plot -E for log axis
          nnn_prev = nnn; // for the next iteration
      } // iE
      return 0;
  } // test_hydrogen_atom

  status_t test_hydrogen_wave_functions(
        radial_grid_t const & g // radial grid descriptor
      , double const Z // number of protons in the nucleus
      , int const echo=0 // log-level
  ) {
      std::vector<double> rf(g.n), rV(g.n, -Z); // fill all potential values with r*V(r) == -Z
      int nnn;
      auto const kink = shoot(0, g, rV.data(), 0, -0.5, nnn, rf.data());
      if (echo > 5) debug(dump_to_file("H1s_radial_wave_function.dat", g.n, rf.data(), g.r));
      return (std::abs(kink) < 1e-3);
  } // test_hydrogen_wave_functions

  
  status_t test_Bessel_functions(int const echo=0) {
      // unit test for the outwards integration
      auto g = *radial_grid::create_exponential_radial_grid(512); // radial grid descriptor
      std::vector<double> rV(g.n, 0); // fill all potential values with r*V(r) == 0 everywhere
      std::vector<double> rf(g.n), rg(g.n); // rf = large component of the radial function
      double const k = 1; // wave number of the Bessel function
      integrate_outwards<0>(rf.data(), rg.data(), g, rV.data(), 0, 0.5*k*k);
      if (echo > 3) std::printf("\n## %s: x, sin(x), f(x), sin(x) - f(x):\n", __func__);
      double dev{0}, norm{0};
      for (auto ir = 1; ir < g.n; ++ir) {
          auto const x = k*g.r[ir];
          double const ref = std::sin(x), diff = ref - rf[ir];
          if (echo > 4) std::printf("%g %g %g %g\n", x, ref, rf[ir], diff); // compare result to x*j0(x)==sin(x)
          norm += pow2(ref)*g.dr[ir]; // integrate
          dev += pow2(diff)*g.dr[ir]; // integrate
      } // ir
      if (echo > 2) std::printf("\n# %s: integral (sin(%g*r) - r*f(r))^2 dr= %g\n", __func__, k, dev);
      if (echo > 5) std::printf("# %s: integral sin^2(%g*r) dr= %g\n", __func__, k, norm);
      radial_grid::destroy_radial_grid(&g);
      return (dev > 7e-11*norm); // return error if deviations become large
  } // test_Bessel_functions

  // unit test for the inhomogeneous outwards integration
  status_t test_inhomogeneous(
        radial_grid_t const & g // radial grid descriptor
      , double const Z=0 // number of protons in the nucleus
      , int const echo=0 // log-level
  ) {
      std::vector<double> mem(4*g.n);
      auto const gg = &mem[0], ff = &mem[g.n], rp = &mem[2*g.n], rV = &mem[3*g.n];
      for (int ir = 0; ir < g.n; ++ir) {
          rp[ir] = g.r[ir]*exp(-0.5*pow2(g.r[ir])); // init inhomogeneity*r
          rV[ir] = -Z; // bare hydrogen potential
      } // ir
      double E{0.5}, dg{0};
      for (ell_QN_t ell = 0; ell < 4; ++ell) { // loop must run serial and forward
          integrate_outwards<0>(gg, ff, g, rV, ell, E/(ell + 1.), -1, &dg, rp);
          if (echo > 3) std::printf("\n## %s for ell=%d: r, f(r), inhomogeneity(r):\n", __func__, ell);
          for (int ir = 0; ir < g.n; ++ir) {
              if (echo > 3 && ir > 0) std::printf("%g %g %g\n", g.r[ir], gg[ir], rp[ir]); //, ff[ir], rV[ir]);
              rp[ir] *= g.r[ir]; // udate rp for the next ell-iteration
          } // ir
      } // ell
      return 0;
  } // test_inhomogeneous

  status_t all_tests(int const echo) {
      status_t stat(0);
//       auto const g = *radial_grid::create_exponential_radial_grid(512);
//    double const Z = 1;
//    stat += test_hydrogen_atom(*radial_grid::create_exponential_radial_grid(256), Z, echo);
//    stat += test_hydrogen_wave_functions(*radial_grid::create_exponential_radial_grid(2610), Z, echo);
      stat += test_Bessel_functions(echo);
//    stat += test_inhomogeneous(g, Z, echo); // solutions need to be inspected manually
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace radial_integrator
