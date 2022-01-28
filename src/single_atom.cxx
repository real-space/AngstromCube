#ifndef DEVEL

/*
 *  Make sure to modify this file only in the development branch
 */
#endif // not DEVEL

#include <cstdio> // std::printf, ::snprintf, ::fflush, stdout
#include <cstring> // std::strncpy
#include <cmath> // std::sqrt, ::abs
#include <cassert> // assert
#include <algorithm> // std::max, ::min
#include <fstream> // std::ofstream

#include "single_atom.hxx"

#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // ::create_default_radial_grid, ::destroy_radial_grid
#include "radial_eigensolver.hxx" // ::shooting_method
#include "radial_potential.hxx" // ::Hartree_potential
#include "angular_grid.hxx" // ::transform, ::get_grid_size
#include "radial_integrator.hxx" // ::integrate_outwards
#include "exchange_correlation.hxx" // ::lda_PZ81_kernel

#include "inline_math.hxx" // pow2, pow3, set, scale, product, add_product, intpow, dot_product, align<nBits>

#include "sho_unitary.hxx" // ::Unitary_SHO_Transform<real_t>
#include "sho_tools.hxx" // ...

#include "solid_harmonics.hxx" // ::lm_index, ::Y00, ::Y00inv
#include "atom_core.hxx" // ::initial_density, ::rad_pot, ::nl_index
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t, emm_Degenerate, spin_QN_t, spin_Degenerate
#include "energy_level.hxx" // TRU, SMT, TRU_AND_SMT, TRU_ONLY, partial_wave_t
#include "spherical_state.hxx" // spherical_state_t, core, semicore, valence, csv_name, show_state_analysis, show_state
#include "display_units.h" // eV, _eV, Ang, _Ang
#include "unit_system.hxx" // ::energy_unit
#include "simple_math.hxx" // ::invert
#include "simple_timer.hxx" // SimpleTimer
#include "bessel_transform.hxx" // ::transform_to_r2grid
#include "scattering_test.hxx" // ::eigenstate_analysis, ::logarithmic_derivative, ::emm_average
#include "linear_algebra.hxx" // ::eigenvalues
#include "data_view.hxx" // view4D<T>, view3D<T>, view2D<T>, transpose, gemm
#include "lossful_compression.hxx" // print_compressed
#include "control.hxx" // ::get
#include "chemical_symbol.hxx" // ::get
#include "sigma_config.hxx" // ::get, element_t
#include "bisection_tools.hxx" // bisector_t
#include "print_tools.hxx" // printf_vector<T>(fmt, vec, n, final="\n", scale=1, add=0)
#include "energy_contribution.hxx" // ::TOTAL, ::KINETIC, ::ELECTROSTATIC, ...

#include "recorded_warnings.hxx" // warn, error, abort
#include "pseudo_tools.hxx" // ::pseudize_local_potential, ::pseudize_function, ...
                            // ::pseudize_spherical_density, ::perform_Gram_Schmidt
#define   HAS_ELEMENT_CONFIG
#ifdef    HAS_ELEMENT_CONFIG
    #include "element_config.hxx" // ::get
#endif // HAS_ELEMENT_CONFIG


#ifdef HAS_RAPIDXML
    #include "pawxml_import.hxx" // ::parse_pawxml, pawxml_t, pawxmlstate_t
    // #include "radial_grid.hxx" // ::create_radial_grid, ::equation_reciprocal
#endif // HAS_RAPIDXML

#include "pawxml_export.hxx" // ::write_to_file

// #define DEBUG
#ifdef  DEBUG
    #include "debug_output.hxx" // here
#else  // DEBUG
    #define here
#endif // DEBUG


namespace single_atom {

  int constexpr ELLMAX=15;
  char const ellchar[] = "spdfghijklmnopq";

  char const ts_name[][8] = {"true", "smooth"};  // TRU:true, SMT:smooth

  double constexpr Y00    = solid_harmonics::Y00; // == 1./sqrt(4*pi)
  double constexpr Y004pi = solid_harmonics::Y00inv; // == sqrt(4*pi)

  template <typename real_t>
  inline void symmetrize(real_t &left, real_t &right) {
      // given two elements, set both of them to their common arithmetic average
      left = (left + right)/2; right = left;
  } // symmetrize

#ifdef DEVEL
  status_t minimize_curvature(
        int const n // dimension
      , view2D<double> & A // kinetic energy operator
      , view2D<double> & B // overlap operator
      , double* lowest=nullptr // optional export coefficient
  ) {
      // solve a generalized eigenvalue problem to find the lowest eigenvector coefficient
      if (n < 1) return 0;
      std::vector<double> eigs(n);
      auto const info = linear_algebra::eigenvalues(eigs.data(), n, A.data(), A.stride(), B.data(), B.stride());
      if (lowest) *lowest = eigs[0];
      return info;
  } // minimize_curvature
#endif // DEVEL

  template <typename int_t>
  int display_delimiter( // returns mln (or mlmn if resolve=='m')
        int const numax // size of the SHO basis
      , int_t const nn[] // numbers of partial waves per ell-channel
      , char const resolve='\0' // 'm':emm_Resolved, otherwise emm_Degenerate
  ) {
      // we want to display only the non-zero SHO contributions and skip higher ell entries

      int const m2 = 2*('m' == resolve);

      int lmax{-1}; // find the largest ell for which there are partial waves
      for (int ell = 0; ell <= numax; ++ell) {
          if (nn[ell] > 0) lmax = ell;
      } // ell

      int mln{0}; // count the number of [emm-degenerate] radial SHO states
      for (int ell = 0; ell <= lmax; ++ell) {
          mln += sho_tools::nn_max(numax, ell)*(m2*ell + 1);
      } // ell
      return mln;
  } // display_delimiter


  template <int ADD0_or_PROJECT1> // 0:add, 1:project, 2:expand, 3:test_normalization
  void add_or_project_compensators(
        view2D<double> & Alm // ADD0_or_PROJECT1 == 0 or 2 result
      , double qlm[]         // ADD0_or_PROJECT1 == 1 or 3 result
      , radial_grid_t const & rg // radial grid descriptor
      , int const lmax       // cutoff for angular momentum expansion
      , double const sigma   // spread of generalized Gaussian
      , int const echo=0     // log-level
  ) {
      // compensation charge densities on the radial grid

      int const nr = rg.n; // number of radial grid points
      if (echo > 0) std::printf("# sigma = %g\n", sigma);
      auto const sig2inv = .5/(sigma*sigma);
      std::vector<double> rl(nr), rlGauss(nr);
      for (int ell = 0; ell <= lmax; ++ell) { // loop must run forward and serial!
          double norm{0};
          for (int ir = 0; ir < nr; ++ir) {
              auto const r = rg.r[ir];
              if (0 == ell) {
                  rl[ir] = 1; // init r^ell as r^0
                  rlGauss[ir] = (ADD0_or_PROJECT1 < 2) ? std::exp(-sig2inv*r*r) : 1; // r^0*Gaussian
              } else {
                  rl[ir]      *= r; // construct r^ell
                  rlGauss[ir] *= r; // construct r^ell*Gaussian
              }
              norm += rlGauss[ir] * rl[ir] * rg.r2dr[ir];
              if (echo > 8) std::printf("# ell=%i norm=%g ir=%i rlGauss=%g rl=%g r2dr=%g\n",
                                      ell, norm, ir, rlGauss[ir], rl[ir], rg.r2dr[ir]);
          } // ir
          if (echo > 1) std::printf("# ell=%i norm=%g nr=%i\n", ell, norm, nr);
          assert(norm > 0);
          auto const scal = 1./norm;
          for (int emm = -ell; emm <= ell; ++emm) {
              int const lm = solid_harmonics::lm_index(ell, emm);
              if (0 == ADD0_or_PROJECT1) {
                  add_product(Alm[lm], nr, rlGauss.data(), qlm[lm]*scal); // add normalized compensator to augmented density
              } else if (2 == ADD0_or_PROJECT1) {
                  add_product(Alm[lm], nr, rl.data(), qlm[lm]); // add q_{\ell m} * r^\ell to electrostatic potential
                  // Note: setup of normalized rlGauss is not necessary in this version
              } else if (3 == ADD0_or_PROJECT1) {
                  qlm[lm] = dot_product(nr, Alm[lm], rl.data(), rg.r2dr); // dot_product with metric
                  // Note: setup of normalized rlGauss is not necessary in this version
              } else {
                  qlm[lm] = dot_product(nr, Alm[lm], rlGauss.data(), rg.r2dr) * scal; // dot_product with metric
              } // add or project
          } // emm
      } // ell

  } // add_or_project_compensators


  

  template <typename int_t>
  void get_valence_mapping(
        int_t ln_index_list[] // ln-index of ilmn (retrieves the emm_Degenerate index)
      , int_t lm_index_list[] // lm-index of ilmn (retrieves the (ell,emm)-combined spherical harmonic quantum number)
      , int_t lmn_begin[] // begin of ilmn indices wich have a given lm-index
      , int_t lmn_end[]   //   end of ilmn indices wich have a given lm-index
      , int const numax // SHO basis size
      , char const *label="" // log-prefix
      , int const echo=0 // log-level
  ) {
      // create various index tables used in the construction of matrix elements

      int const mlm = pow2(1 + numax);
      set(lmn_begin, mlm, int_t(-1));
      for (int ell = 0; ell <= numax; ++ell) {
          for (int emm = -ell; emm <= ell; ++emm) {
              for (int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
                  int const ilmn      = sho_tools::lmn_index(numax, ell, emm, nrn);
                  ln_index_list[ilmn] = sho_tools::ln_index(numax, ell, nrn); // valence state index
                  int const ilm       = sho_tools::lm_index(ell, emm);
                  lm_index_list[ilmn] = ilm;
                  if (lmn_begin[ilm] < 0) lmn_begin[ilm] = ilmn; // store the first index of this lm
                  lmn_end[ilm] = ilmn + 1; // store the last index of this lm
              } // nrn
          } // emm
      } // ell
#ifdef DEVEL
      if (echo > 3) { // display
          std::printf("# %s ln_index_list ", label);
          printf_vector(" %i", ln_index_list, sho_tools::nSHO(numax));
          std::printf("# %s lmn_begin-lmn_end ", label);
          for (int ilm = 0; ilm < mlm; ++ilm) {
              std::printf(" %i", lmn_begin[ilm]);
              if (lmn_begin[ilm] < lmn_end[ilm] - 1) std::printf("-%i", lmn_end[ilm] - 1); // do not display e.g. " 7-7" but only " 7"
          } // ilm
          std::printf("\n");
      } // echo
#endif // DEVEL
  } // get_valence_mapping

   
  double expand_numerical_functions_in_SHO_basis( // return the quality of "occupied" projectors
        double & gradient // second result, gradient of the quality w.r.t. sigma
      , double const sigma // spread of the SHO basis
      , int const numax_basis // size of the SHO basis
      , int const numax // size of the projector set
      , radial_grid_t const & rg // radial grid descriptor, typically the smooth grid
      , view2D<double> const & rfunc // r*functions(numerical), rfunc(nln,>= rg.n)
      , double const weight_ln[] // weights, usually partial wave occupations [nln], negative weights for skip
      , char const *label="" // log-prefix
      , int const echo=0 // log-level
      , view2D<double> *func_coeff=nullptr // optional result
  ) {

      int const nln = sho_tools::nSHO_radial(numax_basis);
      view3D<double> sho_basis(2, nln, align<2>(rg.n), 0.0); // get memory for projector(r) and the derivative w.r.t. sigma

      // expand the normalized radial SHO basis functions for this value of sigma
      scattering_test::expand_sho_projectors(sho_basis(0,0), sho_basis.stride(), rg, sigma, numax_basis, 0, echo >> 1, sho_basis(1,0));

      double weighted_quality{0};
      gradient = 0;
      for (int ell = 0; ell <= numax; ++ell) {
          int const nn_max = sho_tools::nn_max(numax_basis, ell);
          if (func_coeff) assert(nn_max <= func_coeff->stride());
          std::vector<double> denom_sho(nn_max, 0.0);
          { // scope: fill norm2 of the SHO basis functions
              for (int jrn = 0; jrn < nn_max; ++jrn) { // smooth number or radial nodes
                  int const jln = sho_tools::ln_index(numax_basis, ell, jrn);
                  denom_sho[jrn] = dot_product(rg.n, sho_basis(0,jln), sho_basis(0,jln), rg.r2dr); // should be close to 1.0
    //            std::printf("# for sigma= %g %s radial SHO function #%i normalized %g\n", sigma*Ang,_Ang, jln, denom_sho);
              } // jrn
          } // scope

          int const mrn = sho_tools::nn_max(numax, ell);
          for (int nrn = 0; nrn < mrn; ++nrn) { // number of the partial wave
              int const iln = sho_tools::ln_index(numax, ell, nrn); // index of partial wave
              if (weight_ln[iln] >= 0) {
                  auto const denom_num = dot_product(rg.n, rfunc[iln], rfunc[iln], rg.dr); // norm^2 of the numerically given projectors
                  if (denom_num <= 0) {
                      if (weight_ln[iln] > 0) {
                          if (echo > 1) std::printf("# %s for sigma= %g %s cannot normalize %c%i proj: num %g\n", 
                                                      label, sigma*Ang,_Ang, ellchar[ell], nrn, denom_num);
                      } // weight
                  } // cannot be normalized

                  double quality_ln{0}, gradient_ln{0};
                  for (int jrn = 0; jrn < nn_max; ++jrn) { // smooth number or radial nodes
                      int const jln = sho_tools::ln_index(numax_basis, ell, jrn); // index of the radial SHO state

    //                if (echo > 1) std::printf("# %s in iteration %i norm of %c%i projectors: sho %g classical %g\n",
    //                                             label, iter, ellchar[ell], nrn, denom_sho, denom_num);
                      if (denom_sho[jrn]*denom_num > 0) {
                          auto const inner   = dot_product(rg.n, rfunc[iln], sho_basis(0,jln), rg.rdr); // metric is rdr since numerical projectors enter as r*prj
                          auto const d_inner = dot_product(rg.n, rfunc[iln], sho_basis(1,jln), rg.rdr); // derivative w.r.t. sigma
                          auto const denoms  = 1./(denom_sho[jrn]*denom_num);
                          auto const quality = pow2(inner)*denoms;
                          gradient_ln += d_inner * 2 * inner * denoms;
                          quality_ln  += quality;
                          if (echo > 13) std::printf("# %s quality for %c%i with sigma= %g %s is %g\n",
                                                        label, ellchar[ell], nrn, sigma*Ang, _Ang, quality);

                          if (func_coeff) (*func_coeff)(iln,jrn) = inner / std::sqrt(denom_sho[jrn]); // store
                      } else {
                          if (func_coeff) (*func_coeff)(iln,jrn) = 0;
                          if (weight_ln[iln] > 0 && echo > 3) {
                              std::printf("# %s for sigma= %g %s cannot normalize %c%i proj: sho %g num %g\n", 
                                            label, sigma*Ang,_Ang, ellchar[ell], nrn, denom_sho[jrn], denom_num);
                          } // weight echo
                      } // denom > 0
                  } // jrn
                  if (echo > 11) std::printf("# %s quality for %c nrn=%d with sigma= %g %s is %g (weight %.2f)\n", 
                                                label, ellchar[ell], nrn, sigma*Ang, _Ang, quality_ln, weight_ln[iln]);
                  weighted_quality += weight_ln[iln]*quality_ln;
                  gradient         += weight_ln[iln]*gradient_ln;
              } else { // weight
                  if (echo > 11) std::printf("# %s skip %c nrn=%d\n", label, ellchar[ell], nrn);
              } // weight skip
          } // nrn
      } // ell
      if (echo > 9) std::printf("# %s weighted quality with sigma= %g %s is %g\n", 
                                   label, sigma*Ang, _Ang, weighted_quality);

      return weighted_quality; // and gradient
  } // expand_numerical_functions_in_SHO_basis

  
  
  
    double fit_function_set( // returns optimized sigma
        view2D<double> & prj_coeff // (nln,nn_max(numax_basis, 0))
      , int const numax
      , double const weight_ln[] // weights
      , radial_grid_t const & rg // radial grid descriptor, typically the smooth grid
      , view2D<double> const & rfunc // r*functions(numerical), rfunc(nln,>= rg.n)
      , double const sigma_in // input
      , int const numax_basis
      , double const range=1 // 1:no optimization, >1: result in [sigma_in/range, sigma_in*range]
      , char const *label="" // log-prefix    
      , int const echo=0 // log-level
    ) {
        // create smooth partial waves according to Bloechl/GPAW and corresponding
        // numerically defined projector functions. On demand, fit the best sigma
        // to represent those projectors in a SHO basis of given size numax.

        if (echo > 1) std::printf("\n# %s %s\n", label, __func__);

        // We can define tphi(nn[ell], nr_tru) inside the ell-loop
        // Mind that without the experimental section activated by single_atom.suggest.local.potential
        // we could also define sphi(nn[ell], nr_smt) and skin inside the ell-loop, however, 
        // this does not apply to the projectors since they will all be needed later in the fitting.
//         int const nln = sho_tools::nSHO_radial(numax); // == (numax*(numax + 4) + 4)/4

//         std::vector<double> occ_ln(nln, 0.0);
        
        // get weights only
        double total_weight{0};
        for (int ell = 0; ell <= numax; ++ell) {
            for (int nrn = 0; nrn < sho_tools::nn_max(numax, ell); ++nrn) { // partial waves
                int const iln = sho_tools::ln_index(numax, ell, nrn);
                total_weight += std::max(0.0, weight_ln[iln]);
            } // nrn
        } // ell

//         char const *optimized{""};
//         int const optimize_sigma = control::get("single_atom.optimize.sigma", 0.); // 0:no optimization,
          // 1:optimize + use, -1:optimize + display, -10: optimize + display + abort

        double const sigma_old = sigma_in;
        double sigma_opt{sigma_old};
        if (range > 1) {
            if (echo > 2) std::printf("# %s start optimization of sigma at %g %s\n", label, sigma_old*Ang, _Ang);
            // search for the optimal sigma value to best represent the projectors in a SHO basis of size numax
            double const sigma_range[] = {sigma_old/range, sigma_old*range};
            bisection_tools::bisector_t<double> bisection(sigma_range[0], sigma_range[1], 1e-12, '*');
            double sigma_now, weighted_quality{0}, gradient{0}; // sigma_now does not need initialization
            double const original_quality = expand_numerical_functions_in_SHO_basis(gradient,
                    sigma_old, numax_basis, numax, rg, rfunc, weight_ln, label, 0);
            while (bisection.root(sigma_now, gradient, echo/2))
            {
                weighted_quality = expand_numerical_functions_in_SHO_basis(gradient,
                    sigma_now, numax_basis, numax, rg, rfunc, weight_ln, label, echo);
            } // while
            double const best_weighted_quality = weighted_quality;
            sigma_opt = sigma_now;
            if (echo > 5) std::printf("# %s after %d bisection steps found optimized sigma= %g %s (gradient= %.1e)\n", 
                                  label, bisection.get_iterations_needed(), sigma_opt*Ang, _Ang, gradient);

            // the upper limit for the weighted quality is the number of valence electrons
            if (echo > 2) std::printf("\n# %s  original sigma= %g %s has quality %g of max. %g, %.3f %%\n", label, sigma_old*Ang, _Ang,
                                      original_quality, total_weight, original_quality*100/std::max(1., total_weight));
            if (echo > 2) std::printf("# %s optimized sigma= %g %s for numax= %d with quality %g of max. %g, %.3f %%\n\n", label, sigma_opt*Ang, _Ang, numax,
                                      best_weighted_quality, total_weight, best_weighted_quality*100/std::max(1., total_weight));
            if (sigma_range[0]*1.001 > sigma_opt) warn("%s optimal sigma is at the lower end of the analyzed range!", label);
            if (sigma_range[1]*0.999 < sigma_opt) warn("%s optimal sigma is at the upper end of the analyzed range!", label);

//             optimized = "optimized ";
        } // range > 1
        double const sigma_out = sigma_opt; // return value

        // show expansion coefficients of the projectors
        { // scope: fill prj_coeff
            double gradient{0};
            auto const weighted_quality = expand_numerical_functions_in_SHO_basis(gradient,
                sigma_out, numax_basis, numax, rg, rfunc, weight_ln, label, echo, &prj_coeff);
            if (echo > 3) std::printf("\n# %s sigma= %g %s with quality %g of max. %g, %.3f %%\n\n", label, sigma_out*Ang, _Ang, 
                                          weighted_quality, total_weight, weighted_quality*100/std::max(1., total_weight));
        } // scope: fill prj_coeff

        for (int ell = 0; ell <= numax; ++ell) {
            int const nn_max = sho_tools::nn_max(numax, ell);
            int const mn_max = sho_tools::nn_max(numax_basis, ell);
            assert(mn_max <= 8 && "numax_basis > hard limit 15");
            for (int nrn = 0; nrn < nn_max; ++nrn) {
                int const iln = sho_tools::ln_index(numax, ell, nrn);
                if (echo > 0) {
                    std::printf("# %s fitted coefficient ell=%d nrn=%d are", label, ell, nrn);
                    printf_vector(" %g", prj_coeff[iln], mn_max);
                } // echo
            } // nrn
        } // ell

        return sigma_out;
    } // fit_function_set
  



  class LiveAtom {
  public:
      // ToDo: separate everything which is energy-parameter-set dependent 
      //        and group it into a class valence_window_t (or some better name)
      //        so that single_atom can treat eigenvalue equations and
      //        Green functions, i.e. for a fixed energy or interpolated between energies
      //              valence_window_t members:
//                         std::vector<char> partial_wave_active;
//                         view2D<double> hamiltonian, overlap; // matrices [nSHO][>=nSHO]
//                         view2D<double> density_matrix; // atomic density matrix [nSHO][>=nSHO]
//                         view3D<double> kinetic_energy; // tensor [TRU_AND_SMT][nln][nln]
//                         view4D<double> charge_deficit; // tensor [1 + ellmax_cmp][TRU_AND_SMT][nln][nln]
//                         view2D<double> projectors; // [nln][rg[SMT].n] r*projectors, depend only on sigma and numax (MAYBE not needed, only temporary)
//                         view2D<double> projector_coeff[1 + ELLMAX]; // [ell][nn[ell]][nn_max(numax,ell)] coeff of projectors in the SHO basis
//                         view3D<double> partial_wave_radial_part[TRU_AND_SMT]; // matrix [wave0_or_wKin1][nln or less][nr], valence states point into this
//                         std::vector<double> zero_potential; // PAW potential shape correction, POTENTIALLY energy-parameter-set dependent
      //
      //        maybe also separate the spherical states into a class spherical_atom_t
      //        from which the TRU core and semicore density is exported (+valence on demand)
      //        This sperical_atom_t could at some point even be unified with or
      //        replace the atom_core module

      // general config
      int32_t atom_id; // global atom identifier
      double Z_core; // number of protons in the core
      char label[16]; // label string
      radial_grid_t rg[TRU_AND_SMT]; // radial grid descriptor for the true and smooth grid:
              // SMT may point to TRU, but at least both radial grids must have the same tail
      int nr_diff; // how many more radial grid points are in *rg[TRU] compared to *rg[SMT]
      ell_QN_t ellmax_pot; // limit ell for full_potential
      ell_QN_t ellmax_rho; // limit ell for full_density
      ell_QN_t ellmax_cmp; // limit ell for the charge deficit compensators
      double sigma_compensator; // Gaussian spread for the charge deficit compensators
      std::vector<double> qlm_compensator; // coefficients for the charge deficit compensators, (1 + ellmax_cmp)^2

      // the following quantities are energy-parameter-set dependent
      double r_cut; // classical augmentation radius for potential and core density
      int   ir_cut[TRU_AND_SMT]; // classical augmentation radius index for potential and core density
//    float r_match; // radius for matching of true and smooth partial wave, usually 6--9*sigma
      ell_QN_t numax; // limit of the SHO projector quantum numbers
      double sigma; // spread of the SHO projectors and its inverse
      uint8_t nn[1 + ELLMAX]; // number of projectors and partial waves used in each ell-channel
      std::vector<partial_wave_t> partial_wave;
      // the following quantities are energy-parameter-set dependent and spin-resolved (nspins=1 or =2)
      std::vector<char> partial_wave_active; // used as std::vector<bool> but bool vectors have no .data() member
      view2D<double> hamiltonian, overlap; // matrices [nSHO][>=nSHO]
      view3D<double> kinetic_energy; // tensor [TRU_AND_SMT][nln][nln]
      view4D<double> charge_deficit; // tensor [1 + ellmax_cmp][TRU_AND_SMT][nln][nln]
      view2D<double> projectors; // [nln][rg[SMT].n] r*projectors, depend only on sigma and numax
      view2D<double> projector_coeff[1 + ELLMAX]; // [ell][nn[ell]][nn_max(numax,ell)] coeff of projectors in the SHO basis
      view3D<double> partial_wave_radial_part[TRU_AND_SMT]; // matrix [wave0_or_wKin1][nln or less][nr], valence states point into this
      // end of energy-parameter-set dependent members
      view3D<double> true_core_waves; // matrix [wave0_or_wKin1][nln][nr], core states point into this
      std::vector<double> zero_potential; // PAW potential shape correction, potentially energy-parameter-set dependent
      view2D<double> density_matrix; // atomic density matrix [nSHO][>=nSHO]

      view2D<double> aug_density; // augmented density, core + valence + compensation, (1+ellmax_rho)^2 radial functions
      // int nspins; // 1 or 2 (order 0z) or 4 (order 0zxy), so far not used

      // spin-resolved members
      double csv_charge[3]; // [csv], in units of electrons
      double take_spherical_density[3]; // [csv], 1: use the spherical density only, 0: use the density from partial waves, mixtures possible.
      std::vector<spherical_state_t> spherical_state; // 20 core states are the usual max., 32 core states are enough if spin-orbit-interaction is on
      view2D<double> spherical_density[TRU_AND_SMT]; // spherical densities*4pi, no Y00 factor, for {core, semicore, valence}

      view2D<double> full_density[TRU_AND_SMT]; // total density, core + valence, (1 + ellmax_rho)^2 radial functions
      view2D<double> full_potential[TRU_AND_SMT]; // (1 + ellmax_pot)^2 radial functions
      std::vector<double> potential[TRU_AND_SMT]; // spherical potential r*V(r), no Y00 factor, used for the generation of partial waves

      double spherical_charge_deficit[3]; // [csv], in units of electrons

      double logder_energy_range[3]; // [start, increment, stop]
      std::vector<char> partial_wave_char; // [iln]
      double partial_wave_energy_split[1 + ELLMAX]; // [ell]
      int    local_potential_method;

      bool freeze_partial_waves;
      bool regenerate_partial_waves;

      view2D<double> unitary_zyx_lmn; // unitary sho transformation matrix [order_Ezyx][order_lmn], stride=nSHO(numax)

      int8_t gaunt_init;
      std::vector<gaunt_entry_t> gaunt;

      std::vector<int16_t> ln_index_list; // iln = ln_index_list[ilmn]
      std::vector<int16_t> lm_index_list; // ilm = lm_index_list[ilmn]
      std::vector<int16_t> lmn_begin; // lm_index_list[lmn_begin[ilm]] == ilm
      std::vector<int16_t> lmn_end; // lm_index_list[lmn_end[ilm] - 1] == ilm


      // energy contributions
      double energy_xc[TRU_AND_SMT]; // exchange-correlation[ts]
      double energy_dc[TRU_AND_SMT]; // xc double counting[ts]
      double energy_dm; // atomic density matrix double counting, SMT only
      double energy_kin[TRU_AND_SMT];
      double energy_kin_csvn[4][TRU_AND_SMT]; // kinetic[csv][ts], csv=3 --> non-spherical
      double energy_es[TRU_AND_SMT]; // electrostatic[ts]
      double energy_tot[TRU_AND_SMT]; // total energy[ts]
      double energy_ref[TRU_ONLY]; // reference energy

      static int constexpr SRA = 1; // 1: Scalar Relativistic Approximation

      float reference_spdf[3][4]; // max. 3 valence eigenvalues for s,p,d,f to compare in eigenstate_analysis

  public:

    LiveAtom( // constructor method
          double const Z_protons // number of protons in the nucleus
        , bool const atomic_valence_density=false // this option allows to make an isolated atom scf calculation
        , int32_t const global_atom_id=-1 // global atom identifier
        , int const echo=0 // logg level for this constructor method only
        , double const ionization=0 // charged valence configurations
    )
        : // initializations
          atom_id{global_atom_id} 
        , Z_core{Z_protons}
        , gaunt_init{-1}
    {
        // constructor routine

        char chemical_symbol[4];
        chemical_symbol::get(chemical_symbol, Z_core);
        set_label(chemical_symbol);

        if (echo > 0) std::printf("\n\n#\n# %s LiveAtom with %g protons\n", label, Z_core);
#ifndef   HAS_ELEMENT_CONFIG
        if (0 != ionization) error("ionization inactive!");
#endif // HAS_ELEMENT_CONFIG

        take_spherical_density[core]     = 1; // must always be 1 since we can represent the true core density only on the radial grid
        take_spherical_density[semicore] = 1;
        take_spherical_density[valence]  = atomic_valence_density ? 1 : 0;


        // here use the preliminary Z_core, may be adjusted
        rg[TRU] = *radial_grid::create_default_radial_grid(Z_core);
        // create a radial grid descriptor which has less points at the origin
        rg[SMT] = *radial_grid::create_pseudo_radial_grid(rg[TRU], control::get("single_atom.smooth.radial.grid.from", 1e-3));
        // Warning: *rg[TRU] and *rg[SMT] need an explicit destructor call

        int const nr[] = {int(align<2>(rg[TRU].n)), int(align<2>(rg[SMT].n))}; // optional memory access alignment
        if (echo > 0) std::printf("# %s radial grid up to %g %s\n", label, rg[TRU].rmax*Ang, _Ang);
        if (echo > 0) std::printf("# %s radial grid numbers are %d and %d\n", label, rg[TRU].n, rg[SMT].n);
        if (echo > 0) std::printf("# %s radial grid numbers are %d and %d (padded to align)\n", label, nr[TRU], nr[SMT]);

        // allocate spherically symmetric quantities
        for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            spherical_density[ts] = view2D<double>(3, nr[ts], 0.0); // get memory for core, semicore, valence
            potential[ts]       = std::vector<double>(nr[ts], 0.0); // get memory for r*V(r)
        } // true and smooth

        { // scope to load potential[TRU] from file pot/Zeff.00Z
            auto const load_stat = atom_core::read_Zeff_from_file(potential[TRU].data(), rg[TRU], Z_core, "pot/Zeff", -1, echo, label);
            if (0 != load_stat) {
                if ('g' == (*control::get("single_atom.start.potentials", "generate") | 32)) {
                    if (echo > 0) std::printf("\n# %s generate self-consistent atomic potential for Z= %g\n", label, Z_core);
                    auto const gen_stat = atom_core::solve(Z_core, echo/2, 'c', &rg[TRU]);
                    if (0 != gen_stat) error("failed to generate a self-consistent atomic potential for Z= %g", Z_core);
                    // should be able to read from file now
                    auto const s = atom_core::read_Zeff_from_file(potential[TRU].data(), rg[TRU], Z_core, "pot/Zeff", -1, echo, label);
                    if (0 != s) error("loading of potential file failed for Z= %g failed although generated", Z_core);
                } else { // generate
                    error("loading of potential file failed for Z= %g\n#   run -t atom_core +atom_core.test.Z=%g", Z_core, Z_core);
                } // start.potentials=generate
            } // load_stat
        } // scope

#ifdef DEVEL
        // show the loaded Zeff(r) == -r*V(r)
        if (echo > 33) {
           std::printf("\n## loaded Z_eff(r) function:\n");
           for (int ir = 0; ir < rg[TRU].n; ++ir) {
               std::printf("%.15g %.15g\n", rg[TRU].r[ir], -potential[TRU][ir]);
           } // ir
           std::printf("\n\n");
        } // echo
#endif // DEVEL

        std::vector<double> occ_custom(36, 0.); // customized occupation numbers for the radial states
        std::vector<int8_t> csv_custom(36, csv_undefined);
        set(nn, 1 + ELLMAX, uint8_t(0)); // clear

        double const core_state_localization = control::get("single_atom.core.state.localization", -1.); // in units of electrons, -1:inactive
        char const *custom_configuration{"auto"};
        char local_potential_method_name[16];

        { // scope: initialize from sigma_config [or optionally from element_config]
            if (echo > 8) std::printf("# %s get PAW configuration data for Z=%g\n", label, Z_core);
            auto const ec =
#ifdef    HAS_ELEMENT_CONFIG
                            ('e' == (*control::get("single_atom.config", "sigma") | 32)) ? // s:sigma_config, e:element_config
                            element_config::get(Z_core, ionization, rg[TRU], potential[TRU].data(), 
                                             chemical_symbol, core_state_localization, echo) :
#endif // HAS_ELEMENT_CONFIG
                            sigma_config::get(Z_core, echo - 4, &custom_configuration);
            if (echo > 0) std::printf("# %s got PAW configuration data for Z=%g: rcut=%g sigma=%g %s\n",
                                         label, ec.Z, ec.rcut*Ang, ec.sigma*Ang, _Ang);

            if (ec.Z != Z_core) warn("%s number of protons adjusted from %g to %g", label, Z_core, ec.Z);
            Z_core = ec.Z;
            sigma = std::abs(ec.sigma);
            r_cut = std::abs(ec.rcut);
            // get numax from nn[]
            int nu_max{ec.numax}; // -1: choose the minimum automatically
            int const nn_limiter = control::get("single_atom.nn.limit", 2);
            for (int ell = 0; ell < 8; ++ell) {
                if (ec.nn[ell] > 0) {
                    nn[ell] = std::min(int(ec.nn[ell]), nn_limiter); // take a smaller number of partial waves
                    int const numax_min = 2*nn[ell] - 2 + ell; // this is the minimum numax that we need to get ec.nn[ell]
                    nu_max = std::max(nu_max, numax_min);
                } // nn > 0
            } // ell
            numax = nu_max;
            assert(numax >= 0);
            for (int inl = 0; inl < 32; ++inl) {
                occ_custom[inl] = ec.occ[inl][0] + ec.occ[inl][1]; // spin polarization is neglected so far
                csv_custom[inl] = ec.csv[inl];
            } // inl enn-ell all-electron orbital index
            std::strncpy(local_potential_method_name, ec.method, 15);

        } // initialize from sigma_config

        // after this, Z_core may not change any more

        // neutral atom reference energy
        energy_ref[TRU] = atom_core::neutral_atom_total_energy(Z_core);


        if (echo > 1) std::printf("# %s projectors are expanded up to numax= %d\n", label, numax);
        ellmax_rho = 2*numax; // could be smaller than 2*numax
        ellmax_pot = ellmax_rho;
        if (echo > 1) std::printf("# %s radial density and potentials are expanded up to lmax= %d and %d, respectively\n", label, ellmax_rho, ellmax_pot);
        ellmax_cmp = std::min(4, int(ellmax_rho));
        if (echo > 1) std::printf("# %s compensation charges are expanded up to lmax= %d\n", label, ellmax_cmp);

        sigma_compensator = r_cut/std::sqrt(20.); // Gaussian decays to exp(-10)=4.54e-5 at r=r_cut

        nr_diff = rg[TRU].n - rg[SMT].n; // how many more radial grid points are in the radial for TRU quantities compared to the SMT grid
        ir_cut[SMT] = radial_grid::find_grid_index(rg[SMT], r_cut);
        ir_cut[TRU] = ir_cut[SMT] + nr_diff;
        if (echo > 1) std::printf("# %s pseudize the core density at r[%i or %i]= %.6f, requested %.3f %s\n",
                          label, ir_cut[TRU], ir_cut[SMT], rg[SMT].r[ir_cut[SMT]]*Ang, r_cut*Ang, _Ang);
        assert(rg[SMT].r[ir_cut[SMT]] == rg[TRU].r[ir_cut[TRU]]); // should be exactly equal, otherwise r_cut may be to small

        if (echo > 1) {
            std::printf("# %s number of projectors per ell ", label);
            printf_vector(" %d", nn, 1 + numax);
        } // echo
        assert( numax <= ELLMAX );

        // allocate quantities with lm-resolution:
        for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            full_potential[ts] = view2D<double>(pow2(1 + ellmax_pot), nr[ts], 0.0); // get memory
            full_density[ts]   = view2D<double>(pow2(1 + ellmax_rho), nr[ts], 0.0); // get memory
        } // true and smooth



        //
        // Spherical States
        // 
        // A set of atomic eigenstates of the spherical part of the potential is computed
        // to serve several purposes:
        //    - core level energies and the true core density
        //    - energy parameters for valence states
        //    - a spherical valence density to initialize a full calculation
        //    - enable automatic analysis
        //

        std::vector<int8_t> as_valence(36, -1);
        enn_QN_t enn_core_ell[16]; // energy quantum number of the highest occupied core level
        set(enn_core_ell, 16, enn_QN_t(0));

        set(csv_charge, 3, 0.); // clear numbers of electrons for {core, semicore, valence}

        { // scope: initialize the spherical states and true spherical densities
            // typically, there are at most 20 spherical states without spin-orbit coupling
            spherical_state = std::vector<spherical_state_t>(32);
            int highest_occupied_state_index{-1};
            int ics{0}; // init counter of spherical states
            for (enn_QN_t enn = 1; enn < 9; ++enn) { // principal quantum number n
                for (ell_QN_t ell = 0; ell < enn; ++ell) { // angular momentum character l
                    int const inl = atom_core::nl_index(enn, ell);
                    if (echo > 15) std::printf("# %s inl= %d, ics= %d\n", label, inl, ics);
                    assert(inl < 36);
                    if (0 != occ_custom[inl]) {
                        auto & cs = spherical_state[ics]; // abbreviate "core state"

                        cs.energy = atom_core::guess_energy(Z_core, enn);
                        std::snprintf(cs.tag, 7, "%d%c", enn, ellchar[ell]); // create a state label
                        cs.nrn[TRU] = enn - ell - 1; // true number of radial nodes
                        cs.enn = enn;
                        cs.ell = ell;
                        cs.emm = emm_Degenerate;
                        cs.spin = spin_Degenerate;
                        cs.csv = csv_custom[inl];
                        if (valence == cs.csv) as_valence[inl] = ics; // mark as good for the valence band, store the core state index

                        double const occ = std::abs(occ_custom[inl]);
                        cs.occupation = occ;
                        if (occ > 0) {
                            csv_charge[cs.csv] += occ;
                            highest_occupied_state_index = ics; // store the index of the highest occupied core state
                            if (echo > 7) std::printf("# %s %-9s%-4s%6.1f\n", label, csv_name(cs.csv), cs.tag, occ);
                            if (as_valence[inl] < 0) {
                                enn_core_ell[ell] = std::max(enn, enn_core_ell[ell]); // find the largest enn-quantum number of the occupied core states
                            } // not as valence
                        } // occupied

                        ++ics; // add a spherical state
                    } // occupied
                } // ell
            } // enn
            int const n_spherical_states = highest_occupied_state_index + 1; // correct the number of spherical states to those occupied
            spherical_state.resize(n_spherical_states); // destruct unused core states

            // get memory for the true radial wave functions and kinetic waves
            true_core_waves = view3D<double>(2, n_spherical_states, nr[TRU], 0.0);
            for (int ics = 0; ics < n_spherical_states; ++ics) {
                auto & cs = spherical_state[ics]; // abbreviate "core state"
                cs.wave[TRU] = true_core_waves(0,ics); // the true radial function
                cs.wKin[TRU] = true_core_waves(1,ics); // the kinetic energy wave
            } // ics

        } // scope

        if (echo > 5) {
            std::printf("# %s enn_core_ell ", label);
            printf_vector(" %d", enn_core_ell, 4);
        } // echo

        double const total_n_electrons = csv_charge[core] + csv_charge[semicore] + csv_charge[valence];
        if (echo > 2) std::printf("# %s initial occupation with %g electrons: %g core, %g semicore and %g valence electrons\n", 
                                    label, total_n_electrons, csv_charge[core], csv_charge[semicore], csv_charge[valence]);

        { // scope: initialize the spherical states and spherical densities
            float const mixing[3] = {1, 1, 1}; // take 100% of the new density (old densities are zero)
            update_spherical_states(mixing, echo, true, core_state_localization);
        } // scope








        // 
        // Partial Waves
        // 
        // The set of partial waves describes the valence energy window
        //

        int nlnn{0};
        for (int ell = 0; ell <= numax; ++ell) { 
            nlnn += nn[ell]; // count active partial waves
        } // ell
        for (int ts = TRU; ts <= SMT; ++ts) {
            partial_wave_radial_part[ts] = view3D<double>(2, nlnn, nr[ts], 0.0); // get memory for the true/smooth radial wave function and kinetic wave
        } // ts
        
        double constexpr energy_derivative = -8.0;
        double const excited_energy = control::get("single_atom.partial.wave.energy", 1.0); // 1.0 Hartree higher, special function for -8.0:energy derivative
#ifndef DEVEL
        if (energy_derivative == excited_energy) warn("energy derivative for partial wave feature only with -D DEVEL", 0);
#endif // DEVEL
        set(partial_wave_energy_split, 1 + ELLMAX, excited_energy);

        auto const energy_parameter = control::get("single_atom.partial.wave.energy.parameter", -9e9);
        int const prev_energy_parameter = control::get("single_atom.previous.energy.parameter", 0.); // bit array, -1:all, 0:none, 2:p, 6:p+d, 14:p+d+f, ...
        bool const use_energy_parameter = (energy_parameter > -9e9);
        if (use_energy_parameter) {
            if (echo > 5) std::printf("# %s use energy parameter %g %s for all ell-channels\n", label, energy_parameter*eV, _eV);
        } // use_energy_parameter

        int const nln = sho_tools::nSHO_radial(numax); // == (numax*(numax + 4) + 4)/4
        partial_wave_char = std::vector<char>(nln + 1, '\0'); // init partial wave characteristics, +1 for a zero-terminated string
        partial_wave = std::vector<partial_wave_t>(nln);
        partial_wave_active = std::vector<char>(nln, 0);
        { // scope: initialize partial waves
            int ilnn{0}; // index of the existing partial waves
            for (int ell = 0; ell <= numax; ++ell) {
                for (int nrn = 0; nrn <= (numax - ell)/2; ++nrn) { // smooth number or radial nodes

                    int const iln = sho_tools::ln_index(numax, ell, nrn);
                    assert(iln < nln);
                    auto & vs = partial_wave[iln]; // abbreviate "valence state"
                    int const enn = std::max(ell + 1, enn_core_ell[ell] + 1) + nrn;
//                  if (echo > 0) std::printf(" %d%c", enn, ellchar[ell]);
                    vs.nrn[TRU] = enn - ell - 1; // number of radial nodes in the true wave
                    vs.nrn[SMT] = nrn;           // number of radial nodes in the smooth wave
                    vs.occupation = 0.0;
                    vs.enn = enn; // principal quantum number
                    vs.ell = ell; // angular momentum quantum number
                    vs.emm = emm_Degenerate; // angular z-component quantum number
                    vs.spin = spin_Degenerate; // spin quantum number
                    vs.csv = valence; // classifier for {core, semicore, valence}

                    if (nrn < nn[ell]) {
                        partial_wave_active[iln] = 1;
                        assert(ilnn < nlnn);
                        for (int ts = TRU; ts <= SMT; ++ts) {
                            vs.wave[ts] = partial_wave_radial_part[ts](0,ilnn); // pointer for the {true, smooth} radial function
                            vs.wKin[ts] = partial_wave_radial_part[ts](1,ilnn); // pointer for the {true, smooth} kinetic energy
                        } // ts
                        ++ilnn; // add one partial wave
                    } else {
                        for (int ts = TRU; ts <= SMT; ++ts) {
                            vs.wave[ts] = nullptr;
                            vs.wKin[ts] = nullptr;
                        } // ts
                    }

                    int const inl = atom_core::nl_index(enn, ell); // global emm-degenerate orbital index
                    if (partial_wave_active[iln]) {
                        double E{atom_core::guess_energy(Z_core, enn)}; // init with a guess
                        if (nrn > 0) E = std::max(E, partial_wave[iln - 1].energy); // higher than previous energy
                        radial_eigensolver::shooting_method(SRA, rg[TRU], potential[TRU].data(), enn, ell, E, vs.wave[TRU]);
                        vs.energy = E;
                        std::snprintf(vs.tag, 7, "%d%c", enn, ellchar[ell]); // create a state label

                        double occ{occ_custom[inl]};
                        { // scope: transfer valence occupation
                            int const ics = as_valence[inl]; // index of the corresponding spherical state
//                          std::printf("# as_valence[nl_index(enn=%d, ell=%d) = %d] = %d\n", enn, ell, inl, ics);
                            if (ics >= 0) { // atomic eigenstate was marked as valence
                                occ = spherical_state[ics].occupation; //
                                vs.occupation = occ;
                                if (occ > 0) {
//                                  if (echo > 3) std::printf("# %s transfer %.1f electrons from %s-spherical state #%d"
//                                      " to the %s-partial wave #%d\n", label, occ, spherical_state[ics].tag, ics, vs.tag, iln);
                                    if (echo > 3) std::printf("# %s transfer %.1f electrons from %s-spherical state"
                                        " to the %s-partial wave\n", label, occ, spherical_state[ics].tag, vs.tag);
                                } // occupation transfer
                            } // ics
                        } // scope

                        if (echo > 0) show_state(label, csv_name(valence), vs.tag, vs.occupation, E);

                        partial_wave_char[iln] = '0' + enn; // eigenstate: '1', '2', '3', ...
                        if (use_energy_parameter) {
                            vs.energy = energy_parameter;
                            partial_wave_char[iln] = 'e'; // use energy parameter
                            std::snprintf(vs.tag, 7, "%c%.5g", ellchar[ell], vs.energy*eV); // create a state label
                        } else // use_energy_parameter
                        if (occ > 0) {
                            // eigenstate unchanged, see above
                        } else { // occ > 0
                            if (nrn > 0) {
                                char asterisk[8] = "*******"; asterisk[nrn] = '\0';
                                std::snprintf(vs.tag, 7, "%c%s", ellchar[ell], asterisk); // create a state label
                                partial_wave_char[iln] = '*'; // set as excited
#ifdef DEVEL
                                if (energy_derivative == partial_wave_energy_split[ell]) partial_wave_char[iln] = 'D';
#endif // DEVEL
                                if ('*' == partial_wave_char[iln] && std::abs(partial_wave_energy_split[ell]) < 1e-3) {
                                    warn("%s partial wave energy split for ell=%c is small (%g %s)",
                                          label, ellchar[ell], partial_wave_energy_split[ell]*eV, _eV);
                                } // |dE| small
                            } else
                            if (prev_energy_parameter & (1 << ell)) { // nrn > 0
                                if (0 == ell) warn("%s cannot make use of the energy parameter of the previous ell-channel when ell=0", label);
                                if (ell > 0) partial_wave_char[iln] = 'p'; // use base energy of previous ell-channel (polarization)
                            } else {
                                // eigenstate unchanged, see above
                            } // nrn > 0
                        } // occ > 0

                    } else { // active
                        partial_wave_char[iln] = '_'; // partial wave inactive
                        std::snprintf(vs.tag, 7, "%c?", ellchar[ell]); // create a label for inactive states
                    } // partial_wave_active

                    if (echo > 19) std::printf("# %s created a state descriptor with tag=\'%s\'\n", label, vs.tag);
                } // nrn
            } // ell
            assert(nlnn == ilnn && "count of partial waves incorrect");
        } // scope: initialize partial waves

        set(reference_spdf[0], 3*4, 0.f); // clear reference valence eigenvalues

        int const nlm_aug = pow2(1 + std::max(ellmax_rho, ellmax_cmp));
        aug_density = view2D<double>(nlm_aug, full_density[SMT].stride(), 0.0); // get memory
        int const nlm_cmp = pow2(1 + ellmax_cmp);
        qlm_compensator = std::vector<double>(nlm_cmp, 0.0); // get memory
        charge_deficit = view4D<double>(1 + ellmax_cmp, TRU_AND_SMT, nln, nln, 0.0); // get memory
        kinetic_energy = view3D<double>(TRU_AND_SMT, nln, nln, 0.0); // get memory
        zero_potential = std::vector<double>(nr[SMT], 0.0); // get memory

        projectors = view2D<double>(nln, nr[SMT], 0.0); // get memory, radial representation
        for (int ell = 0; ell <= numax; ++ell) {
            projector_coeff[ell] = view2D<double>(nn[ell], sho_tools::nn_max(numax, ell), 0.0); // get memory, block-diagonal in ell
            for (int nrn = 0; nrn < nn[ell]; ++nrn) {
                projector_coeff[ell](nrn,nrn) = 1.0; // Kronecker
            } // nrn
        } // ell

        int const nSHO = sho_tools::nSHO(numax);
        int const matrix_stride = align<2>(nSHO); // 2^<2> doubles = 32 Byte alignment
        if (echo > 0) std::printf("# %s matrix size for hamiltonian and overlap: dim= %d, stride= %d\n", label, nSHO, matrix_stride);
        density_matrix = view2D<double>(nSHO, matrix_stride, 0.0); // get memory
        hamiltonian    = view2D<double>(nSHO, matrix_stride, 0.0); // get memory
        overlap        = view2D<double>(nSHO, matrix_stride, 0.0); // get memory

        unitary_zyx_lmn = view2D<double>(nSHO, nSHO, 0.0);
        { // scope: fill the Unitary_SHO_Transform with values from a file
            sho_unitary::Unitary_SHO_Transform<double> const u(numax);
            auto const stat = u.construct_dense_matrix(unitary_zyx_lmn.data(), numax, nSHO, sho_tools::order_zyx, sho_tools::order_lmn);
            assert(0 == int(stat));
        } // scope

        int const mlm = pow2(1 + numax);
        ln_index_list.resize(nSHO);
        lm_index_list.resize(nSHO);
        lmn_begin.resize(mlm);
        lmn_end.resize(mlm);
        get_valence_mapping(ln_index_list.data(), lm_index_list.data(), lmn_begin.data(), lmn_end.data(), numax, label, echo - 5);

        { // scope: specify the range of the logarithmic derivative analysis
            char const *_eu;
            auto const eu = unit_system::energy_unit(control::get("logder.unit", "Ha"), &_eu);
            auto const in_eu = 1./eu;
            logder_energy_range[0] = control::get("logder.start", -2.0*eu)*in_eu;
            logder_energy_range[1] = control::get("logder.step",  1e-2*eu)*in_eu; // ToDo: these getter calls should be moved to the main function
            logder_energy_range[2] = control::get("logder.stop",   1.0*eu)*in_eu;
            if (echo > 3) std::printf("# %s logder.start=%g logder.step=%g logder.stop=%g %s\n",
                label, logder_energy_range[0]*eV, logder_energy_range[1]*eV, logder_energy_range[2]*eV, _eV);
            if (eu != eV) {
                if (echo > 4) std::printf("# %s logder.start=%g logder.step=%g logder.stop=%g %s\n",
                    label, logder_energy_range[0]*eu, logder_energy_range[1]*eu, logder_energy_range[2]*eu, _eu);
            } // logder.unit != output.energy.unit
        } // scope

        { // scope: pseudize the local potential r*V(r)
            local_potential_method = ('p' == (*local_potential_method_name | 32)) ? 0 : //  0: parabola fit
                              int(control::get("single_atom.lagrange.derivative", 7.)); // >0: sinc fit
            pseudo_tools::pseudize_local_potential<1>(potential[SMT].data(),
                    potential[TRU].data(), rg, ir_cut, local_potential_method, label, echo);
            // ToDo: use a sinc to pseudize an attractive potential. Result should be very similar to parabola
        } // scope

        regenerate_partial_waves = true; // must be true at start to generate the partial waves at least once
        freeze_partial_waves = (control::get("single_atom.relax.partial.waves", 1.) < 1);

        float const potential_mixing = 0.25; // this controls the percentage of the full_potential that is ...
        // ... taken to relax the spherical potential from which spherical states and true partial wave are computed
        float const density_mixing[3] = {.25, .25, .25}; // [csv]

        update_density(density_mixing, echo); // run the density update at least once to generate partial waves

        int const maxit_scf = control::get("single_atom.init.max.scf", 0.);
        if (echo > 1) std::printf("# %s run single_atom.init.max.scf=%d initial SCF-iterations\n", label, maxit_scf);
        int const echo_minimal = std::min(echo, int(control::get("single_atom.init.scf.echo", 1.)));
        for (int scf = 0; scf < maxit_scf; ++scf) {
            int const echo_scf = (scf > maxit_scf - 2) ? echo : echo_minimal; // turn on output only for the last iteration
            if (echo_scf > echo_minimal) std::printf("\n\n"); // spacer
            if (echo > 1) std::printf("# %s SCF-iteration %d%c", label, scf, (echo_scf > echo_minimal) ? '\n': '\t');
            update_density(density_mixing, echo_scf);
            double const *const ves_multipoles = nullptr; // no potential shifts
            update_potential(potential_mixing, ves_multipoles, echo_scf);
        } // self-consistency iterations

#ifdef DEVEL

        int const export_xml = control::get("single_atom.export.xml", 0.);
        if (export_xml) {
            update_potential(potential_mixing, nullptr, echo); // compute the zero_potential and total energy contributions
            if (echo > 0) std::printf("\n\n# %s export configuration to file\n", label);
            auto const stat = pawxml_export::write_to_file(Z_core, rg, 
                partial_wave, partial_wave_active.data(),
                kinetic_energy, csv_charge, spherical_density, projectors,
                r_cut, sigma_compensator, zero_potential.data(), echo,
                energy_kin_csvn[core][TRU],
                energy_kin[TRU], energy_xc[TRU], energy_es[TRU], energy_tot[TRU],
                custom_configuration,
                exchange_correlation::default_LDA,
                control::get("single_atom.export.path", "."));
            if (stat) warn("pawxml_export::write_to_file returned status= %i", int(stat));
            if (control::get("single_atom.optimize.sigma", 0.) > 0) {
                warn("optimized sigma may differ from config string for Z=%g", Z_core);
            }
            if (echo > 0) std::printf("# %s exported configuration to file\n", label);
            if (export_xml < 0) abort("single_atom.export.xml=%d (negative leads to a stop)", export_xml);
        } // export_xml

        // show the smooth and true potential
        if (false && echo > 0) {
            std::printf("\n## %s spherical parts: r in Bohr, "
            "Zeff_tru(r), Zeff_smt(r)"
            ", r^2*rho_tot_tru(r), r^2*rho_tot_smt(r)"
            ", r^2*rho_cor_tru(r), r^2*rho_cor_smt(r)"
            ", r^2*rho_val_tru(r), r^2*rho_val_smt(r)"
            ", zero_potential(r) in Ha"
            ":\n", label);
            for (int irs = 0; irs < rg[SMT].n; irs += 1) {
                int const irt = irs + nr_diff;
                auto const r = rg[SMT].r[irs], r2 = pow2(r);
                auto const ct = spherical_density[TRU](core,irt)*r2*Y00*Y00;
                auto const cs = spherical_density[SMT](core,irs)*r2*Y00*Y00;
                auto const vt = spherical_density[TRU](valence,irt)*r2*Y00*Y00;
                auto const vs = spherical_density[SMT](valence,irs)*r2*Y00*Y00;
                std::printf("%g %g %g %g %g %g %g %g %g %g\n", r
//                         , -full_potential[TRU](00,irt)*Y00*r // for comparison, should be the same as Z_eff(r)
                        , -potential[TRU][irt] // Z_eff(r)
                        , -potential[SMT][irs] // \tilde Z_eff(r)
                        , ct + vt, cs + vs
                        , ct, cs
                        , vt, vs
                        , zero_potential[irs]*Y00
                      );
            } // ir
            std::printf("\n\n");

            if (dot_product(rg[TRU].n, spherical_density[TRU][semicore], rg[TRU].r2dr) > 0) {
                warn("%s semicore density was not plotted", label);
            } // non-vanishing semicore density

        } // echo

        // show the smooth and true partial waves
        if (false && echo > 0) {
            std::printf("\n\n\n# %s valence partial waves:\n", label);
            for (int ir = 0; ir < rg[SMT].n; ir += 1) {
                auto const r = rg[SMT].r[ir];
                auto const f = r; // scale with r? scale with Y00?
                std::printf("%g ", r);
                for (int ics = 0; ics < spherical_state.size(); ++ics) {
                    std::printf("   %g", spherical_state[ics].wave[TRU][ir + nr_diff]*f);
                } // ics
                for (int iln = 0; iln < nln; ++iln) {
                    if (partial_wave_active[iln]) {
                        std::printf("   %g %g", partial_wave[iln].wave[TRU][ir + nr_diff]*f
                                              , partial_wave[iln].wave[SMT][ir]*f);
                    } // active
                } // iln
                std::printf("\n");
            } // ir
            std::printf("\n\n\n\n");
        } // echo

#endif // DEVEL

    } // constructor



    
    
#ifdef HAS_RAPIDXML

    LiveAtom( // constructor method from pawxml
          double const Z_protons // number of protons in the nucleus
        , int const echo=0 // logg level for this constructor method only
    )
        : // initializations
          atom_id{-1}
        , gaunt_init{-1}
    {
        // constructor routine

        char chemical_symbol[4];
        chemical_symbol::get(chemical_symbol, Z_protons);
        set_label(chemical_symbol);

        auto const pawpath = control::get("single_atom.pawxml.path", "gpaws");
        char xmlfilename[512]; std::snprintf(xmlfilename, 511, "%s/%s.LDA", pawpath, chemical_symbol);
        if (echo > 0) std::printf("\n\n#\n# %s LiveAtom loads \'%s\'\n", label, xmlfilename);
        auto const p = pawxml_import::parse_pawxml(xmlfilename, echo);
        if (0 != p.parse_status) error("%s parsing \'%s\' returned status=%d", label, xmlfilename, int(p.parse_status));

        Z_core = p.Z;
        if (echo > 3) std::printf("\n\n#\n# %s loading of \'%s\' successful, %g protons\n", label, xmlfilename, Z_core);
        if (Z_protons != Z_core) warn("%s number of protons adjusted from %g to %g", label, Z_protons, Z_core);

        rg[TRU] = *radial_grid::create_radial_grid(p.n, p.n*p.radial_grid_a, radial_grid::equation_reciprocal);
        rg[SMT] = rg[TRU]; rg[SMT].memory_owner = false; // same grid for true and smooth quantities, shallow copy


        take_spherical_density[core]     = 1; // must always be 1 since we can represent the true core density only on the radial grid
        take_spherical_density[semicore] = 1;
        take_spherical_density[valence]  = 0; // use a synthetic density matrix instead

        int const nr[] = {int(align<2>(rg[TRU].n)), int(align<2>(rg[SMT].n))}; // optional memory access alignment
        if (echo > 0) std::printf("# %s radial grid up to %g %s\n", label, rg[TRU].rmax*Ang, _Ang);
        if (echo > 0) std::printf("# %s radial grid numbers are %d and %d\n", label, rg[TRU].n, rg[SMT].n);
        if (echo > 0) std::printf("# %s radial grid numbers are %d and %d (padded to align)\n", label, nr[TRU], nr[SMT]);

        // allocate spherically symmetric quantities
        for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            spherical_density[ts] = view2D<double>(3, nr[ts], 0.0); // get memory for core, semicore, valence
            potential[ts]       = std::vector<double>(nr[ts], 0.0); // get memory for r*V(r)
        } // true and smooth

        std::vector<double> occ_custom(36, 0.); // customized occupation numbers for the radial states
        std::vector<int8_t> csv_custom(36, csv_undefined);
        std::vector<int8_t> ist_custom(36, -1); // state indices
        set(nn, 1 + ELLMAX, uint8_t(0)); // clear

        int ncmx[4]; // largest enn of the core electrons
        sigma_config::set_default_core_shells(ncmx, Z_core);
        if (echo > 6) std::printf("# %s preliminary core states up to %ds %dp %dd %df\n", label, ncmx[0], ncmx[1], ncmx[2], ncmx[3]);

        r_cut = rg[SMT].rmax; // init at maximum
        std::vector<int8_t> enn_ell(1 + ELLMAX, 0);
        int ist{0};
        for (auto & s : p.states) {
            int const enn = s.n, ell = s.l;
            if (echo > 4) std::printf("# %s valence state %d%c E= %.6f %s\n", label, enn, ellchar[ell], s.e*eV, _eV);
            assert(0 <= ell); assert(ell <= ELLMAX);
            ++nn[ell]; // increase the number of partial waves for this ell
            if (enn > 0) {
                if (ell < 4) ncmx[ell] = std::min(ncmx[ell], enn - 1);
                enn_ell[ell] = enn;
                int const inl = atom_core::nl_index(enn, ell);
                assert(inl < 36);
                ist_custom[inl] = ist;
                assert(0 <= s.f); assert(s.f <= 2*(ell + 1 + ell)); // sanity check for occupation numbers
                occ_custom[inl] = s.f; // copy occupation numbers
                csv_custom[inl] = valence;
            } else {
                int const enn_prime = std::max(ell + 1, enn_ell[ell] + 1);
                int const inl = atom_core::nl_index(enn_prime, ell);
                assert(inl < 36);
                ist_custom[inl] = ist;
                ++enn_ell[ell];
            } // enn valid
            r_cut = std::min(r_cut, double(s.rc));
            ++ist;
        } // valence states
        assert(p.states.size() == ist && "fatal counting error");
        if (echo > 3) std::printf("# %s core states up to %ds %dp %dd %df\n", label, ncmx[0], ncmx[1], ncmx[2], ncmx[3]);
        if (echo > 3) std::printf("# %s smallest cutoff radius is %g %s\n", label, r_cut*Ang, _Ang);

        { // scope: determine numax
            int nu_max{-1};
            for (int ell = 0; ell <= ELLMAX; ++ell) {
                if (nn[ell] > 0) {
                    int const numax_min = 2*nn[ell] - 2 + ell; // this is the minimum numax that we need to get for nn[ell]
                    nu_max = std::max(nu_max, numax_min);
                } // nn > 0
            } // ell
            numax = nu_max;
            assert(numax >= 0);
            if (echo > 3) std::printf("# %s numax is %d\n", label, numax);
        } // scope

        // fill the core
        for (int ell = 0; ell < 4; ++ell) {
            for (int enn = ell + 1; enn <= ncmx[ell]; ++enn) {
                int const inl = atom_core::nl_index(enn, ell);
                if (csv_undefined == csv_custom[inl]) {
                    csv_custom[inl] = core;
                    occ_custom[inl] = 2*(ell + 1 + ell);
                    if (echo > 8) std::printf("# %s %d%c is a fully occupied core state\n", label, enn, ellchar[ell]);
                } // new core state
            } // enn
        } // ell

        sigma_compensator = p.shape_function_rc*std::sqrt(.5); // exp(-(r/rc)^2) == exp(-1/2*(r/sigma_compensator)^2)

        energy_ref[TRU] = p.ae_energy[3]; // reference total energy of neutral atom calculation


        if (echo > 1) std::printf("# %s projectors are expanded up to numax= %d\n", label, numax);
        ellmax_rho = 2*numax; // could be smaller than 2*numax
        ellmax_pot = ellmax_rho;
        if (echo > 1) std::printf("# %s radial density and potentials are expanded up to lmax= %d and %d, respectively\n", label, ellmax_rho, ellmax_pot);
        ellmax_cmp = std::min(2, int(ellmax_rho)); // see https://wiki.fysik.dtu.dk/gpaw/algorithms.html at "Compensation charges"
        if (echo > 1) std::printf("# %s compensation charges are expanded up to lmax= %d\n", label, ellmax_cmp);

        nr_diff = rg[TRU].n - rg[SMT].n; // how many more radial grid points are in the radial for TRU quantities compared to the SMT grid
        ir_cut[SMT] = radial_grid::find_grid_index(rg[SMT], r_cut);
        assert(0 == nr_diff && "same radial grid for true and smooth quantities");
        ir_cut[TRU] = ir_cut[SMT] + nr_diff;
        if (echo > 1) std::printf("# %s pseudize the core density at r[%i or %i]= %.6f, requested %.3f %s\n",
                          label, ir_cut[TRU], ir_cut[SMT], rg[SMT].r[ir_cut[SMT]]*Ang, r_cut*Ang, _Ang);
        assert(rg[SMT].r[ir_cut[SMT]] == rg[TRU].r[ir_cut[TRU]]); // should be exactly equal, otherwise r_cut may be to small

        if (echo > 1) {
            std::printf("# %s number of projectors per ell ", label);
            printf_vector(" %d", nn, 1 + numax);
        } // echo
        assert( numax <= ELLMAX );

        // allocate quantities with lm-resolution:
        for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            full_potential[ts] = view2D<double>(pow2(1 + ellmax_pot), nr[ts], 0.0); // get memory
            full_density[ts]   = view2D<double>(pow2(1 + ellmax_rho), nr[ts], 0.0); // get memory
        } // true and smooth

        //
        // Spherical States are inactive
        //

        set(csv_charge, 3, 0.); // clear numbers of electrons for {core, semicore, valence}
        for (int inl = 0; inl < 36; ++inl) {
            if (csv_undefined != csv_custom[inl]) {
                csv_charge[csv_custom[inl]] += occ_custom[inl];
            } // csv defined
        } // inl

        double const total_n_electrons = csv_charge[core] + csv_charge[semicore] + csv_charge[valence];
        if (echo > 2) std::printf("# %s initial occupation with %g electrons: %g core, %g semicore and %g valence electrons\n", 
                                    label, total_n_electrons, csv_charge[core], csv_charge[semicore], csv_charge[valence]);

        // 
        // Partial Waves
        // 
        // The set of partial waves describes the valence energy window
        //

        int nlnn{0};
        for (int ell = 0; ell <= numax; ++ell) { 
            nlnn += nn[ell]; // count active partial waves
        } // ell
        for (int ts = TRU; ts <= SMT; ++ts) {
            partial_wave_radial_part[ts] = view3D<double>(1, nlnn, nr[ts], 0.0); // get memory for the true/smooth radial wave function, no kinetic wave
        } // ts

        int const nln = sho_tools::nSHO_radial(numax); // == (numax*(numax + 4) + 4)/4

        projectors = view2D<double>(nln, nr[SMT], 0.0); // get memory, radial representation
        for (int ell = 0; ell <= numax; ++ell) {
            projector_coeff[ell] = view2D<double>(nn[ell], sho_tools::nn_max(numax, ell), 0.0); // get memory, block-diagonal in ell
            for (int nrn = 0; nrn < nn[ell]; ++nrn) {
                projector_coeff[ell](nrn,nrn) = 1.0; // Kronecker
            } // nrn
        } // ell

        std::vector<int> ist_index(nln, -1); // translation table from iln indices to ist indices in p.states
        partial_wave_char = std::vector<char>(nln, '\0'); // init partial wave characteristics
        partial_wave = std::vector<partial_wave_t>(nln);
        partial_wave_active = std::vector<char>(nln, 0);
        { // scope: initialize partial waves
            int ilnn{0}; // index of the existing partial waves
            for (int ell = 0; ell <= numax; ++ell) {
                for (int nrn = 0; nrn <= (numax - ell)/2; ++nrn) { // smooth number or radial nodes

                    int const iln = sho_tools::ln_index(numax, ell, nrn);
                    assert(iln < nln);
                    auto & vs = partial_wave[iln]; // abbreviate "valence state"
                    int const enn = std::max(ell + 1, ((ell < 4)?ncmx[ell]:0) + 1) + nrn;
                    if (echo > 15) std::printf("# %s enn= %d, ell= %d, nrn= %d, iln= %d\n", label, enn, ell, nrn, iln);
                    vs.nrn[TRU] = enn - ell - 1; // number of radial nodes in the true wave
                    vs.nrn[SMT] = nrn;           // number of radial nodes in the smooth wave
                    vs.occupation = 0.0;
                    vs.energy = 0.0;
                    vs.enn = enn; // principal quantum number
                    vs.ell = ell; // angular momentum quantum number
                    vs.emm = emm_Degenerate; // angular z-component quantum number
                    vs.spin = spin_Degenerate; // spin quantum number
                    vs.csv = valence; // classifier for {core, semicore, valence}

                    if (nrn < nn[ell]) {
                        partial_wave_active[iln] = 1;
                        assert(ilnn < nlnn);
                        for (int ts = TRU; ts <= SMT; ++ts) {
                            vs.wave[ts] = partial_wave_radial_part[ts](0,ilnn); // pointer for the {true, smooth} radial function
                            vs.wKin[ts] = nullptr; // _radial_part[ts](1,ilnn); // pointer for the {true, smooth} kinetic energy
                        } // ts
                        ++ilnn; // add one partial wave
                    } else {
                        for (int ts = TRU; ts <= SMT; ++ts) {
                            vs.wave[ts] = nullptr;
                            vs.wKin[ts] = nullptr;
                        } // ts
                    } // active

                    if (partial_wave_active[iln]) {
                        std::snprintf(vs.tag, 7, "%d%c", enn, ellchar[ell]); // create a state label

                        int const inl = atom_core::nl_index(enn, ell); // global emm-degenerate orbital index
                        int const ist = ist_custom[inl]; // state index
                        if (ist >= 0) {
                            assert(ist < p.states.size());
                            ist_index[iln] = ist;
                            vs.energy = p.states[ist].e;
                            vs.occupation = occ_custom[inl];
                        } // ist valid

                        partial_wave_char[iln] = '0' + enn; // eigenstate: '1', '2', '3', ...
                        if (vs.occupation > 0) {
                            // keep eigenstate notation
                        } else { // occ > 0
                            if (nrn > 0) {
                                char asterisk[8] = "*******"; asterisk[nrn] = '\0';
                                std::snprintf(vs.tag, 7, "%c%s", ellchar[ell], asterisk); // create a state label
                                partial_wave_char[iln] = '*'; // set as excited
                            } else { // nrn > 0
                                // keep eigenstate notation
                            } // nrn > 0
                        } // occ > 0

                        if (ist >= 0) {
                            // copy the true and smmoth partial wave function
                            auto const & tsp = p.states[ist].tsp;
                            for (int ts = TRU; ts <= SMT; ++ts) {
                                if (tsp[ts].size() != rg[ts].n) error("%s %s partial wave %s has %ld grid points but expected %d",
                                                                      label, ts_name[ts], vs.tag, tsp[ts].size(), rg[ts].n);
                                assert(tsp[ts].size() == rg[ts].n);
                                set(vs.wave[ts], rg[ts].n, tsp[ts].data());
                            } // ts
                            // copy the numerical projector function
                            {   assert(tsp[2].size() == rg[SMT].n);
                                set(projectors[iln], rg[SMT].n, tsp[2].data());
                            }
                        } // ist >= 0
                        
                        // add to inital valence density
                        if (vs.occupation > 0) {
                            assert(vs.wave[TRU]); // make sure the pointer is not null
                            double const norm2 = dot_product(nr[TRU], vs.wave[TRU], vs.wave[TRU], rg[TRU].r2dr);
                            if (echo > 9) std::printf("# %s state \'%s\' is normalized as %g\n", label, vs.tag, norm2);
                            assert(norm2 > 0);
                            assert(vs.wave[SMT]); // make sure the pointer is not null
                            for (int ts = TRU; ts <= SMT; ++ts) {
                                add_product(spherical_density[ts][valence], nr[ts], vs.wave[ts], vs.wave[ts], vs.occupation/norm2);
                            } // ts
                        } // occupied

                        if (echo > 0) show_state(label, csv_name(valence), vs.tag, vs.occupation, vs.energy);

                    } else { // active
                        partial_wave_char[iln] = '_'; // partial wave inactive
                        std::snprintf(vs.tag, 7, "%c?", ellchar[ell]); // create a label for inactive states
                    } // partial_wave_active

                    if (echo > 19) std::printf("# %s created a state descriptor with tag=\'%s\'\n", label, vs.tag);
                } // nrn
            } // ell
            assert(nlnn == ilnn && "count of partial waves incorrect");
        } // scope: initialize partial waves

        set(reference_spdf[0], 3*4, 0.f); // clear reference valence eigenvalues

        int const nlm_aug = pow2(1 + std::max(ellmax_rho, ellmax_cmp));
        aug_density = view2D<double>(nlm_aug, full_density[SMT].stride(), 0.0); // get memory
        int const nlm_cmp = pow2(1 + ellmax_cmp);
        qlm_compensator = std::vector<double>(nlm_cmp, 0.0); // get memory
        charge_deficit = view4D<double>(1 + ellmax_cmp, TRU_AND_SMT, nln, nln, 0.0); // get memory
        kinetic_energy = view3D<double>(TRU_AND_SMT, nln, nln, 0.0); // get memory
        zero_potential = std::vector<double>(nr[SMT], 0.0); // get memory

        int const nSHO = sho_tools::nSHO(numax);
        int const matrix_stride = align<2>(nSHO); // 2^<2> doubles = 32 Byte alignment
        if (echo > 0) std::printf("# %s matrix size for hamiltonian and overlap: dim= %d, stride= %d\n", label, nSHO, matrix_stride);
        density_matrix = view2D<double>(nSHO, matrix_stride, 0.0); // get memory
        hamiltonian    = view2D<double>(nSHO, matrix_stride, 0.0); // get memory
        overlap        = view2D<double>(nSHO, matrix_stride, 0.0); // get memory

        unitary_zyx_lmn = view2D<double>(nSHO, nSHO, 0.0);
        { // scope: fill the Unitary_SHO_Transform with values from a file
            sho_unitary::Unitary_SHO_Transform<double> const u(numax);
            auto const stat = u.construct_dense_matrix(unitary_zyx_lmn.data(), numax, nSHO, sho_tools::order_zyx, sho_tools::order_lmn);
            assert(0 == int(stat));
        } // scope

        int const mlm = pow2(1 + numax);
        ln_index_list.resize(nSHO);
        lm_index_list.resize(nSHO);
        lmn_begin.resize(mlm);
        lmn_end.resize(mlm);
        get_valence_mapping(ln_index_list.data(), lm_index_list.data(), lmn_begin.data(), lmn_end.data(), numax, label, echo - 5);

        { // scope: specify the range of the logarithmic derivative analysis
            char const *_eu;
            auto const eu = unit_system::energy_unit(control::get("logder.unit", "Ha"), &_eu);
            auto const in_eu = 1./eu;
            logder_energy_range[0] = control::get("logder.start", -2.0*eu)*in_eu;
            logder_energy_range[1] = control::get("logder.step",  1e-2*eu)*in_eu; // ToDo: these getter calls should be moved to the main function
            logder_energy_range[2] = control::get("logder.stop",   1.0*eu)*in_eu;
            if (echo > 3) std::printf("# %s logder.start=%g logder.step=%g logder.stop=%g %s\n",
                label, logder_energy_range[0]*eV, logder_energy_range[1]*eV, logder_energy_range[2]*eV, _eV);
            if (eu != eV) {
                if (echo > 4) std::printf("# %s logder.start=%g logder.step=%g logder.stop=%g %s\n",
                    label, logder_energy_range[0]*eu, logder_energy_range[1]*eu, logder_energy_range[2]*eu, _eu);
            } // logder.unit != output.energy.unit
        } // scope

        local_potential_method = -1; // 0: parabola fit, -1: do not recompute zero_potential

        regenerate_partial_waves = false; // must be true at start to generate the partial waves at least once
        freeze_partial_waves = true; // do not try to regenerate_partial_waves on a radial_exponential_grid
        
        { // scope: copy remaining radial functions
            auto const sq4pi = std::sqrt(4*constants::pi);
            // true core density
            assert(p.func[0].size() == rg[TRU].n);
            set(spherical_density[TRU][core], rg[TRU].n, p.func[0].data(), sq4pi);
            // smmoth core density
            assert(p.func[1].size() == rg[SMT].n);
            set(spherical_density[SMT][core], rg[SMT].n, p.func[1].data(), sq4pi);
            // zero_potential
            assert(p.func[4].size() == rg[SMT].n);
            set(zero_potential.data(), rg[SMT].n, p.func[4].data());
        } // scope
        if (echo > 0) std::printf("# %s zero_potential at origin %g %s\n", label, zero_potential[0]*Y00*eV, _eV);

        
        { // scope: copy kinetic energy difference matrix of partial waves
            int const nstates = p.states.size();
            assert(nstates*nstates == p.dkin.size());
            for (int iln = 0; iln < nln; ++iln) {
                int const ist = ist_index[iln];
                if (ist >= 0) {
                    assert(ist < nstates);
                    for (int jln = 0; jln < nln; ++jln) {
                        int const jst = ist_index[jln];
                        if (jst >= 0) {
                            auto const dkin = p.dkin[ist*nstates + jst]; // import kinetic energy deficit
                            // dkin[][] is the difference kinetic_energy[TRU][][] - kinetic_energy[SMT][][] but we do not know the exact numbers
                            kinetic_energy(TRU,iln,jln) = dkin;
                        } // jst
                    } // jln
                } // ist
            } // iln
        } // scope

        sigma = 0.5*r_cut; // estimate ToDo: sigma from optimizing the projector representation in SHO basis
        
        // ToDo: Gram-Schmidt orthogonalize projectors against lower partial waves and higher partial waves agains projectors
        // ToDo: rotate the kinetic_energy deficit matrix acccordingly

        // kinetic energy of core electrons
        set(energy_kin_csvn[0], 4*2, 0.0);
        energy_kin_csvn[core][TRU] = p.core_energy_kinetic;

        update_charge_deficit(echo);
        // scope: compute charge deficit of spherical densities
        for (int csv = core; csv <= valence; ++csv) { 
            double const cd[] = {dot_product(rg[TRU].n, spherical_density[TRU][csv], rg[TRU].r2dr),
                                 dot_product(rg[SMT].n, spherical_density[SMT][csv], rg[SMT].r2dr)};
            spherical_charge_deficit[csv] = cd[TRU] - cd[SMT];
            if (echo > 0 && csv_charge[csv] > 0) {
                std::printf("# %s true and smooth %s densities have %g and %g electrons, respectively\n",
                               label, csv_name(csv), cd[TRU], cd[SMT]);
            } // echo
        } // scope

        float const density_mixing[] = {0, 0, 0};
        bool const synthetic_density_matrix = true;
        update_density(density_mixing, echo, synthetic_density_matrix);
        update_potential(0.f, nullptr, echo);
        

    } // constructor from pawxml

#endif // HAS_RAPIDXML




    

    ~LiveAtom() { // destructor
        radial_grid::destroy_radial_grid(&rg[SMT], ts_name[SMT]);
        radial_grid::destroy_radial_grid(&rg[TRU], ts_name[TRU]);
    } // destructor

    

    void initialize_Gaunt(int const echo=0) {
        if (gaunt_init < 0) {
            gaunt = angular_grid::create_numerical_Gaunt(6, echo);
            gaunt_init = 6;
        } // not initialized
    } // initialize_Gaunt

    
    
    void update_spherical_states(
          float const mixing[3]
        , int const echo=0
        , bool const norm_warning=false
        , double const core_state_localization=-1
    ) {
        if (echo > 2) std::printf("\n# %s %s Z=%g\n", label, __func__, Z_core);
        int const n_spherical_states = spherical_state.size();
        if (n_spherical_states < 1) return;
        // core states are feeling the spherical part of the hamiltonian only
        auto const & g = rg[TRU]; // the radial grid for true quantities
        std::vector<double> r2rho(g.n);
        double const *const rV_tru = potential[TRU].data();
        view2D<double> new_r2density(3, align<2>(g.n), 0.0); // new TRU densities for {core, semicore, valence}
        double nelectrons[]     = {0, 0, 0}; // core, semicore, valence
        double band_energy[]    = {0, 0, 0}; // core, semicore, valence
        double kinetic_energy[] = {0, 0, 0}; // core, semicore, valence
        double Coulomb_energy[] = {0, 0, 0}; // core, semicore, valence
#ifdef DEVEL
        if (echo > 17) {
            std::printf("# %s %s: solve for eigenstates of the full radially symmetric potential\n", label, __func__);
            std::printf("\n## %s %s: r, -Zeff(r)\n", label, __func__);
            print_compressed(g.r, rV_tru, g.n);
        } // echo
#endif // DEVEL

        for (int ics = 0; ics < n_spherical_states; ++ics) { // private r2rho, reduction(+:new_r2density)
            auto & cs = spherical_state[ics]; // abbreviate "core state"
            double const occ = std::abs(cs.occupation);
            radial_eigensolver::shooting_method(SRA, g, rV_tru, cs.enn, cs.ell, cs.energy, cs.wave[TRU], r2rho.data());
            int const csv = cs.csv; assert(core <= csv && csv <= valence); // {core, semicore, valence}
            auto const norm = dot_product(g.n, r2rho.data(), g.dr);
            auto const norm_factor = (norm > 0)? 1./std::sqrt(norm) : 0;
            if (norm_warning && norm <= 0) warn("%s spherical %s-state cannot be normalized, Z= %g", cs.tag, label, Z_core);
            // warn if occupied and not normalizable
            auto const scal = pow2(norm_factor)*occ; // scaling factor for the density contribution of this state
            // transform r*wave(r) as produced by the radial_eigensolver to wave(r)
            // and normalize the core level wave function to one
            scale(cs.wave[TRU], g.n, g.rinv, norm_factor);
            // create wKin for the computation of the kinetic energy density
            product(cs.wKin[TRU], g.n, rV_tru, cs.wave[TRU], -1.); // start as wKin = -r*V(r)*wave(r)
            add_product(cs.wKin[TRU], g.n, g.r, cs.wave[TRU], cs.energy); // now wKin = r*(E - V(r))*wave
            cs.kinetic_energy = dot_product(g.n, cs.wKin[TRU], cs.wave[TRU], g.rdr);
            // for more precision, we could eval
            //    auto const E_pot = dot_product(g.n, rho, rV_tru, g.rdr);
            //    cs.kinetic_energy = cs.energy - Epot;
            // with rho=pow2(norm_factor)*r2rho/r^2

            if (scal > 0) {
                add_product(new_r2density[csv], g.n, r2rho.data(), scal);
                nelectrons[csv] += occ;
                band_energy[csv] += occ*cs.energy;
                kinetic_energy[csv] += occ*cs.kinetic_energy;
            } // scal > 0
        } // ics

        if (true) { // scope: warnings

            double constexpr x = 9e9; // extreme value
            double max_energy[] = {-x, -x, -x}, min_energy[] = {x, x, x};
            double E_hfos{-9e9}, E_hpos{-9e9}; // energy of the highest fully/partially occupied state
            int  ics_hfos{-1}, ics_hpos{-1}; // indices
            int iref[4] = {0, 0, 0, 0};

            for (int ics = 0; ics < n_spherical_states; ++ics) { // serial loop
                auto const & cs = spherical_state[ics]; // abbreviate "core state"
                int const csv = cs.csv;
                double const E = cs.energy;
                double const occ = std::abs(cs.occupation);
                double const max_occ = 2.*(2*cs.ell + 1);
                if (std::abs(occ - max_occ) < 1e-15) { // spherical state is fully occupied
                    if (E > E_hfos) { E_hfos = E; ics_hfos = ics; }
                } else if (occ > 1e-15) { // this state is partically occupied (considerably)
                    if (E > E_hpos) { E_hpos = E; ics_hpos = ics; }
                }
                min_energy[csv] = std::min(min_energy[csv], E);
                max_energy[csv] = std::max(max_energy[csv], E);

                // display state tag, occupation and energy
                auto const charge_outside = show_state_analysis(echo - 5, label, g, cs.wave[TRU],
                                    cs.tag, cs.occupation, cs.energy, csv_name(csv), ir_cut[TRU]);
                if (core_state_localization > 0) {
                    auto const csv_auto = (charge_outside > core_state_localization) ? valence : core;
                    if (csv_auto !=    csv) warn("%s the spherical %s state is %s, but has %g %% of its charge outside",
                                                  label, cs.tag, csv_name(csv), charge_outside*100);
                    if (echo > 11) std::printf("# %s the spherical %s %s state has %g %% charge outside the sphere, suggest %s\n",
                                                  label, cs.tag, csv_name(csv), charge_outside*100, csv_name(csv_auto));
                } // check core state criterion

                // create reference eigenvalues to compare them automatically
                if (valence == cs.csv && cs.ell < 4) {
                    reference_spdf[iref[cs.ell]][cs.ell] = cs.energy;
                    ++iref[cs.ell];
                } // valence
            } // ics

            if (E_hpos < E_hfos && ics_hfos >= 0 && ics_hpos >= 0) {
                warn("%s Energy of the partially occupied %s-state is %g < %g %s, the energy of the fully occupied %s-state", label,
                    spherical_state[ics_hpos].tag, E_hpos*eV, E_hfos*eV, _eV, spherical_state[ics_hfos].tag);
            } // the state occupation is not the ground state, partially occupied states are not at the Fermi level

            for (int csl = core; csl <= semicore; ++csl) {
                for (int svh = csl + 1; svh <= valence; ++svh) {
                    if (max_energy[csl] > min_energy[svh]) {
                        warn("%s some %s states are higher than %s states", label, csv_name(csl), csv_name(svh));
                    } // band overlap
                } // svh
            } // csl

        } // scope: warnings

        // report integrals
        for (int csv = 0; csv < 3; ++csv) { // {core, semicore, valence}
            double const old_charge = dot_product(g.n, g.r2dr, spherical_density[TRU][csv]);
            double const new_charge = dot_product(g.n, g.dr, new_r2density[csv]);

            double mix_new = mixing[csv], mix_old = 1.0 - mix_new;
            // rescale the mixing coefficients such that the desired number of electrons comes out
            auto const mixed_charge = mix_old*old_charge + mix_new*new_charge;
            if (0.0 != mixed_charge) {
                auto const rescale = nelectrons[csv]/mixed_charge;
                mix_old *= rescale;
                mix_new *= rescale;
            } // rescale

            // generate a mixed density
            double density_change{0}, Coulomb_change{0}, check_charge{0}; // some stats
            for (int ir = 0; ir < g.n; ++ir) {
                auto const new_rho = new_r2density(csv,ir)*pow2(g.rinv[ir]); // *r^{-2}
                auto const old_rho = spherical_density[TRU](csv,ir);
                density_change  += std::abs(new_rho - old_rho)*g.r2dr[ir];
                Coulomb_change  +=         (new_rho - old_rho)*g.rdr[ir]; // Coulomb integral change
                Coulomb_energy[csv] +=      new_rho           *g.rdr[ir]; // Coulomb integral
                auto const mixed_rho = mix_new*new_rho + mix_old*old_rho;
                spherical_density[TRU](csv,ir) = mixed_rho;
                check_charge    += mixed_rho*g.r2dr[ir];
            } // ir
            Coulomb_energy[csv] *= -Z_core;
            Coulomb_change      *= -Z_core; // to get an energy estimate

            if (nelectrons[csv] > 0) {
#ifdef DEVEL
                if (echo > 7) std::printf("# %s previous %s density has %g electrons, expected %g\n"
                                          "# %s new %s density has %g electrons\n",
                                          label, csv_name(csv), old_charge, nelectrons[csv], 
                                          label, csv_name(csv), new_charge);
                auto const dev = check_charge - nelectrons[csv];
                if (std::abs(dev) > 2e-15*Z_core) warn("%s New spherical %s density has %g electrons, expected %g, diff %.1e e",
                                                        label, csv_name(csv), check_charge, nelectrons[csv], dev);
#endif // DEVEL
                if (echo > 3) std::printf("# %s %-8s density change %g e, Coulomb energy change %g %s\n", 
                    label, csv_name(csv), density_change, Coulomb_change*eV, _eV);
                if (echo > 5) std::printf("# %s %-8s E_Coulomb= %.9f  E_band= %.9f  E_kinetic= %.9f %s\n", 
                    label, csv_name(csv), Coulomb_energy[csv]*eV, band_energy[csv]*eV, kinetic_energy[csv]*eV, _eV);
            } // output only for contributing densities

            energy_kin_csvn[csv][SMT] = 0;
            energy_kin_csvn[csv][TRU] = kinetic_energy[csv]; // store

            spherical_charge_deficit[csv] = pseudo_tools::pseudize_spherical_density(
                spherical_density[SMT][csv],
                spherical_density[TRU][csv], rg, ir_cut, csv_name(csv), label, echo*(nelectrons[csv] > 0));
        } // csv

        if (echo > 2) {
            double const w111[] = {1, 1, 1};
            std::printf("# %s total    E_Coulomb= %.9f  E_band= %.9f  E_kinetic= %.9f %s\n",
                    label, dot_product(3, w111, Coulomb_energy)*eV,
                           dot_product(3, w111, band_energy)*eV, 
                           dot_product(3, w111, kinetic_energy)*eV, _eV);
        } // echo

    } // update_spherical_states



    void show_ell_block_diagonal(
          view2D<double> const & matrix_ln
        , char const *title=""
        , double const unit=1
        , bool const all_i=false // otherwise only active partial waves
        , bool const all_j=false // otherwise only active partial waves
    ) const {
        int const nlnr = sho_tools::nSHO_radial(numax);
        view2D<char> ln_label(nlnr, 4);
        sho_tools::construct_label_table<4>(ln_label.data(), numax, sho_tools::order_ln);
        int const mlnr = display_delimiter(numax, nn);
        for (int iln = 0; iln < mlnr; ++iln) {
            if (partial_wave_active[iln] || all_i) {
                std::printf("# %s %s %-4s ", label, title, all_i ? ln_label[iln] : partial_wave[iln].tag);
                for (int jln = 0; jln < mlnr; ++jln) {
                    if (partial_wave_active[jln] || all_j) {
                        if (partial_wave[iln].ell == partial_wave[jln].ell) {
                            std::printf(" %11.6f", matrix_ln(iln,jln)*unit);
                        } else {
                            std::printf("            ");
                        } // ells match
                    } // only active or all
                } // jln
                std::printf("\n");
            } // only active or all
        } // iln
        std::printf("\n");
    } // show_ell_block_diagonal

    
    
    
    
    
    void show_projector_coefficients(char const *attribute="", double const sigma=1, int const echo=0) {
        // Total energy of SHO state is (nu + 3/2) sigma^-2 in Hartree, nu=ell+2*nrn
        // Kinetic energy is half of the total energy (due to x <--> p symmetry)
        double const kinetic_pref = (sigma > 0) ? 0.5/pow2(sigma) : 0; // in Hartree units
        double kinetic_coeff[8];
        for (int ell = 0; ell <= numax; ++ell) {
            int const nmx = sho_tools::nn_max(numax, ell);
            assert(nmx <= 8);
            for (int irn = 0; irn < nn[ell]; ++irn) {
                int const iln = sho_tools::ln_index(numax, ell, irn);
                auto const norm2 = dot_product(nmx, projector_coeff[ell][irn], projector_coeff[ell][irn]);
                for (int i = 0; i < nmx; ++i) kinetic_coeff[i] = projector_coeff[ell][irn][i] * (ell + 2*i + 1.5); // nu + 3/2
                auto const E_kin = kinetic_pref * dot_product(nmx, projector_coeff[ell][irn], kinetic_coeff);
                if (norm2 > 0) {
                    if (echo > 0) {
                        auto const length = std::sqrt(norm2); assert(length > 0);
                        std::printf("# %s %scoefficients of the %s-projector: %.6e * [", label, attribute, partial_wave[iln].tag, length);
                        printf_vector(" %9.6f", projector_coeff[ell][irn], nmx, " ]", 1./length);
                        if (sigma > 0) std::printf(", %.3f %s", 2*E_kin/norm2, unit_system::_Rydberg); // display the minimum required plane wave cutoff
                        std::printf("\n");
                    } // echo
                } else {
                    warn("%s failed to normalize %s-projector coefficients", label, partial_wave[iln].tag);
                } // norm2 > 0
            } // irn
        } // ell
    } // show_projector_coefficients 








    void update_energy_parameters(int const echo=0) {
        // determine the energy parameters for the partial waves

        if (echo > 2) std::printf("\n# %s %s Z=%g partial wave characteristics=\"%s\"\n",
                                       label, __func__, Z_core, partial_wave_char.data());

        double previous_ell_energy{0};
        for (int ell = 0; ell <= numax; ++ell) { // loop runs serial, loop-carried dependency on previous_ell_energy
            for (int nrn = 0; nrn < nn[ell]; ++nrn) { // smooth number or radial nodes, serial due to partial_wave[iln - 1].energy
                int const iln = sho_tools::ln_index(numax, ell, nrn);
                auto & vs = partial_wave[iln]; // abbreviate ('vs' stands for valence state)

                char const c = partial_wave_char[iln];
                assert ('_' != c && "energy parameters only need to be set for active partial waves");
                if ('e' == c) {
                    // energy parameter, vs.energy has already been set earlier
                    if (echo > 4) std::printf("# %s the %s partial wave uses energy parameter E= %g %s\n", 
                                                 label, vs.tag, vs.energy*eV, _eV);
                } else if ('D' == c) {
#ifndef DEVEL
                    error("%s partial_wave_char may only be 'D' with -D DEVEL", label);
#else  // DEVEL
                    // energy derivative at the energy of the lower partial wave
                    assert(nrn > 0);
                    vs.energy = partial_wave[iln - 1].energy;
                    if (echo > 4) std::printf("# %s the %s partial wave is an energy derivative at E= %g %s\n", 
                                                 label, vs.tag, vs.energy*eV, _eV);
#endif // DEVEL
                } else if ('*' == c) {
                    assert(nrn > 0);
                    vs.energy = partial_wave[iln - 1].energy + partial_wave_energy_split[ell];
                    if (echo > 4) std::printf("# %s the %s partial wave is at E= %g %s\n",
                                                 label, vs.tag, vs.energy*eV, _eV);
                } else if ('p' == c) {
                    if (0 == ell) warn("%s energy parameter for ell=%i nrn=%i undetermined", label, ell, nrn);
                    vs.energy = previous_ell_energy;
                    if (echo > 4) std::printf("# %s the %s partial wave is at E= %g %s, copy %c-energy\n", 
                                                 label, vs.tag, vs.energy*eV, _eV, (ell > 0)?ellchar[ell - 1]:'?');
                } else {
                    assert(c == '0' + vs.enn);
                    // find the eigenenergy of the TRU spherical potential
                    assert(vs.wave[TRU] != nullptr);
                    radial_eigensolver::shooting_method(SRA, rg[TRU], potential[TRU].data(), vs.enn, ell, vs.energy, vs.wave[TRU]);
                    if (echo > 4) std::printf("# %s the %s partial wave is at E= %g %s, the %i%c-eigenvalue\n",
                                                 label, vs.tag, vs.energy*eV,_eV, vs.enn, ellchar[ell]);
                } // c
                if (0 == nrn) previous_ell_energy = vs.energy;
            } // nrn
        } // ell

    } // update_energy_parameters










    double update_projector_coefficients( // returns optimized sigma if optimization is active
        int const echo=0 // log-level
    ) {
        // create smooth partial waves according to Bloechl/GPAW and corresponding
        // numerically defined projector functions. On demand, fit the best sigma
        // to represent those projectors in a SHO basis of given size numax.

        if (echo > 1) std::printf("\n# %s %s Z=%g\n", label, __func__, Z_core);

        // We can define tphi(nn[ell], nr_tru) inside the ell-loop
        // Mind that without the experimental section activated by single_atom.suggest.local.potential
        // we could also define sphi(nn[ell], nr_smt) and skin inside the ell-loop, however, 
        // this does not apply to the projectors since they will all be needed later in the fitting.

        // create partial waves with polynomials like GPAW
        // and preliminary projectors as suggested by Bloechl
        update_partial_waves(echo - 6, 'C', 2); // suppress output, 'C':classical_scheme, 2xGramSchmidt

        int const nln = sho_tools::nSHO_radial(numax); // == (numax*(numax + 4) + 4)/4

        std::vector<double> occ_ln(nln, -1.); // negative occupations for inactive projectors
        double total_occ{0};

        // get weights only
        for (int ell = 0; ell <= numax; ++ell) {
            for (int nrn = 0; nrn < nn[ell]; ++nrn) { // active partial waves
                int const iln = sho_tools::ln_index(numax, ell, nrn);
                assert(partial_wave_active[iln]);
                auto const occ = partial_wave[iln].occupation;
                // the total deviation is weighted with the ell-channel-summed occupation numbers
                occ_ln[iln] = occ;
                total_occ  += occ;
            } // nrn
        } // ell        
        
        int const optimize_sigma = control::get("single_atom.optimize.sigma", 0.); // 0:no optimization,
                            // 1:optimize + use, -1:optimize + display, -10: optimize + display + abort
#if 0
        char const *optimized{""};
#ifdef DEVEL
        double const sigma_old = sigma; // copy member variable
        double sigma_opt{sigma_old};
        if (optimize_sigma) {
            if (echo > 2) std::printf("# %s start optimization of sigma at %g %s\n", label, sigma_old*Ang, _Ang);
            // search for the optimal sigma value to best represent the projectors in a SHO basis of size numax
            double const sigma_range[] = {0.5*sigma_old, 2.0*sigma_old};
            bisection_tools::bisector_t<double> bisection(sigma_range[0], sigma_range[1], 1e-12, '*');
            double sigma_now, weighted_quality{0}, gradient{0}; // sigma_now does not need initialization
            double const original_quality = expand_numerical_functions_in_SHO_basis(gradient,
                    sigma_old, numax, nn, rg[SMT], projectors, occ_ln.data(), label, 0);
            while (bisection.root(sigma_now, gradient, echo >> 1))
            {
                weighted_quality = expand_numerical_functions_in_SHO_basis(gradient,
                    sigma_now, numax, nn, rg[SMT], projectors, occ_ln.data(), label, echo);
            } // while
            double const best_weighted_quality = weighted_quality;
            sigma_opt = sigma_now;
            if (echo > 5) std::printf("# %s after %d bisection steps found optimized sigma= %g %s (gradient= %.1e)\n", 
                                  label, bisection.get_iterations_needed(), sigma_opt*Ang, _Ang, gradient);

            // the upper limit for the weighted quality is the number of valence electrons
            if (echo > 2) std::printf("\n# %s  original sigma= %g %s has quality %g of max. %g, %.3f %%\n", label, sigma_old*Ang, _Ang,
                                      original_quality, total_occ, original_quality*100/std::max(1., total_occ));
            if (echo > 2) std::printf("# %s optimized sigma= %g %s for numax= %d with quality %g of max. %g, %.3f %%\n\n", label, sigma_opt*Ang, _Ang, numax,
                                      best_weighted_quality, total_occ, best_weighted_quality*100/std::max(1., total_occ));
            if (sigma_range[0]*1.01 > sigma_opt) warn("%s optimal sigma is at the lower end of the analyzed range!", label);
            if (sigma_range[1]*0.99 < sigma_opt) warn("%s optimal sigma is at the upper end of the analyzed range!", label);

            optimized = "optimized ";
        } // optimize_sigma
        double const sigma_out = (optimize_sigma > 0 || optimize_sigma < -9) ? sigma_opt : sigma_old; // return value
#else  // DEVEL
        double const sigma_out = sigma; // return value
        if (optimize_sigma) warn("single_atom.optimize.sigma=%i active only with -D DEVEL", optimize_sigma);
#endif // DEVEL

        // show expansion coefficients of the projectors
        view2D<double> prj_coeff(nln, 8, 0.0); // get memory, plain memory layout [nln][8]
        { // scope: fill prj_coeff
            double gradient{0};
            auto const weighted_quality = expand_numerical_functions_in_SHO_basis(gradient, 
                sigma_out, numax, nn, rg[SMT], projectors, occ_ln.data(), label, echo, prj_coeff.data());
            if (echo > 3) std::printf("\n# %s sigma= %g %s with quality %g of max. %g, %.3f %%\n\n", label, sigma_out*Ang, _Ang, 
                                          weighted_quality, total_occ, weighted_quality*100/std::max(1., total_occ));
        } // scope: fill prj_coeff

#else  // 0


//     double fit_function_set( // returns optimized sigma
//         view2D<double> & prj_coeff // (nln, 8)
//       , int const numax
//       , double const weight_ln[] // weights
//       , radial_grid_t const & rg // radial grid descriptor, typically the smooth grid
//       , view2D<double> const & rfunc // r*functions(numerical), rfunc(nln,>= rg.n)
//       , double const sigma_in // input
//       , int const numax_basis
//       , char const *label="" // log-prefix    
//       , int const echo=0 // log-level
//     ) {

        view2D<double> prj_coeff(nln, 8, 0.0); // get memory, plain memory layout [nln][8]
        double const sigma_out = fit_function_set( // returns optimized sigma
            prj_coeff // (nln, 8)
          , numax
          , occ_ln.data()
          , rg[SMT] // radial grid descriptor, typically the smooth grid
          , projectors // r*functions(numerical), rfunc(nln,>= rg.n)
          , sigma // input
          , numax // same numax
          , optimize_sigma ? 2.0 : 1.0 // range
          , label // log-prefix    
          , echo); // log-level
#endif // 0

        for (int ell = 0; ell <= numax; ++ell) {
            int const nn_max = sho_tools::nn_max(numax, ell);
            assert(nn_max <= 8 && "numax > hard limit 15");

            for (int irn = 0; irn < nn[ell]; ++irn) {
                int const iln = sho_tools::ln_index(numax, ell, irn);

                // display
                assert(nn_max == projector_coeff[ell].stride());
                set(projector_coeff[ell][irn], nn_max, prj_coeff[iln]); // copy into member variable

            } // irn

        } // ell

        show_projector_coefficients(optimize_sigma?"optimized ":"", sigma_out, echo - 6); // and warn if not normalizable

#ifdef DEVEL
        if (optimize_sigma < -9) abort("%s after sigma optimization, sigma= %g %s", label, sigma_out*Ang, _Ang);
#endif // DEVEL

        return sigma_out;
    } // update_projector_coefficients










    void update_partial_waves(
          int const echo=0
        , char const generation_method='?'
        , int const GS_iterations=-1
    ) {
        // generate new smooth partial waves given fixed projector functions

        // classical schemes:
        char constexpr classical_scheme = 'C';          // find partial waves by fitting a polynomial and evaluate preliminary projector functions according to Bloechl
        char constexpr recreate_second = 'r';           // construct the higher projector orthogonal to the lowest smooth partial wave
#ifdef DEVEL
        char constexpr classical_partial_waves = 'c';   // find partial waves by fitting a polynomial, SHO-filtered projector functions unchanged
                                                        // works ok but does not give perfect agreement in logder, instable s-charge deficit eigenvalues in Cu if numax=4
                                                        // but stable for numax=3, maybe limit the number of radial SHO basis functions used to nn[ell]
        // revPAW schemes:
        char constexpr minimize_radial_curvature = 'm'; // as suggested by Morian Sonnet: minimize the radial curvature of the smooth partial wave
        char constexpr energy_ordering = 'e';           // as suggested by Baumeister+Tsukamoto in PASC19 proceedings
        char constexpr orthogonalize_second = '2';      // take the same lowest partial wave as for nn==1 and use the freedom of the second
        char constexpr orthogonalize_first  = '1';      // ToDo
#endif // DEVEL

        char const method = ('?' != generation_method) ? generation_method :
                      *control::get("single_atom.partial.wave.method", "r");

        int const Gram_Schmidt_iterations = (GS_iterations >= 0) ? GS_iterations :
                              control::get("single_atom.gram.schmidt.repeat", 2.);

        double r_match{r_cut}; // init classical

        view2D<double> radial_sho_basis; // r*SHO_{ln}(r)
        if (method != classical_scheme) {
            sigma = update_projector_coefficients(echo); // and potentially run optimization on sigma

            int const nln = sho_tools::nSHO_radial(numax);
            radial_sho_basis = view2D<double>(nln, align<2>(rg[SMT].n), 0.0); // get memory
            // create a set of radial SHO basis function for given sigma
            scattering_test::expand_sho_projectors(radial_sho_basis.data(),
                radial_sho_basis.stride(), rg[SMT], sigma, numax, 1, echo >> 1);
#ifdef DEVEL
            if (echo > 5) { // show normalization and orthogonality of the radial SHO basis
                double max_dev{0};
                for (int ell = 0; ell <= numax; ++ell) {
                    for     (int irn = 0; irn < sho_tools::nn_max(numax, ell); ++irn) {   int const iln = sho_tools::ln_index(numax, ell, irn);
                        for (int jrn = 0; jrn < sho_tools::nn_max(numax, ell); ++jrn) {   int const jln = sho_tools::ln_index(numax, ell, jrn);
                            auto const ovl_ij = dot_product(rg[SMT].n, radial_sho_basis[iln], radial_sho_basis[jln], rg[SMT].dr);
                            auto const dev = ovl_ij - (irn == jrn);
                            if (echo > 7) std::printf("# %s radial SHO basis <%c%d|%c%d> = %i %c %.1e sigma=%g %s\n", label,
                                ellchar[ell],irn, ellchar[ell],jrn, irn == jrn, (dev < 0)?'-':'+', std::abs(dev), sigma*Ang, _Ang);
                            max_dev = std::max(max_dev, std::abs(dev));
                        } // jrn
                    } // irn
                } // ell
                std::printf("# %s radial SHO basis deviates %.1e from diagonal, sigma= %g %s\n", label, max_dev, sigma*Ang, _Ang);
            } // echo

            if (echo > 9) {
                std::printf("\n## %s show the local potentials (r, r*Vtru, r*Vsmt):\n", label);
                for (int ir = 1; ir < rg[SMT].n; ++ir) {
                    std::printf("%g %g %g\n", rg[SMT].r[ir], potential[TRU][ir + nr_diff], potential[SMT][ir]);
                } // ir
                std::printf("\n\n");
            } // echo
#endif // DEVEL

            r_match = 9*sigma; // exp(-9^2/2) = 2.6e-18, all projectors are zero beyond this r_match
        } // not classical_scheme

        if (echo > 2) std::printf("\n# %s %s Z=%g method=\'%c\'\n", label, __func__, Z_core, method);
        // the basis for valence partial waves is generated from the spherical part of the hamiltonian
        if ('?' == generation_method && classical_scheme == method) 
            warn("%s Classical scheme leads to invalid potentials, option C only for internal use!", label);


        int const ir_match[] = {radial_grid::find_grid_index(rg[TRU], r_match),
                                radial_grid::find_grid_index(rg[SMT], r_match)};
        if (echo > 3) std::printf("# %s matching radius %g %s at radial indices %i and %i\n",
                                     label, r_match*Ang, _Ang, ir_match[TRU], ir_match[SMT]);

        int const nr = rg[TRU].n;
        std::vector<double> r2rho(nr); // must be thread-private if OMP parallel over ell


        for (int ell_ = 0; ell_ <= numax; ++ell_) { int const ell = ell_;
            int const ln_off = sho_tools::ln_index(numax, ell, 0); // offset where to start indexing emm-degenerate partial waves
            int const n = nn[ell]; // abbreviation for the number of active partial waves in this ell
            int const nmx = sho_tools::nn_max(numax, ell);

            if (echo > 3 && n > 0) std::printf("\n# %s %s for ell=%i\n\n", label, __func__, ell);

            view2D<double> projectors_ell(projectors[ln_off], projectors.stride()); // sub-view of the member array
            if (method != classical_scheme) {
                // construct projector functions as linear combination of the radial SHO basis functions
                for (int nrn = 0; nrn < n; ++nrn) {
                    int const iln = ln_off + nrn;
                    // construct the projectors from the linear combinations of the radial SHO basis
                    set(projectors_ell[nrn], projectors_ell.stride(), 0.0); // clear
                    for (int mrn = 0; mrn < nmx; ++mrn) { // radial SHO basis size
                        auto const c = projector_coeff[ell](nrn,mrn);
                        if (0.0 != c) {
                            add_product(projectors_ell[nrn], rg[SMT].n, radial_sho_basis[ln_off + mrn], c);
                            if (echo > 5 && (0 == nrn || recreate_second != method)) {
                                std::printf("# %s construct %s projector by taking %9.6f of the %c%i radial SHO basis function\n",
                                               label, partial_wave[iln].tag, c, ellchar[ell],mrn);
                            } // echo
                        } // coefficient non-zero
                    } // mrn
                } // nrn
                if (echo > 4 && n > 0) std::printf("# %s %c-projectors prepared\n\n", label, ellchar[ell]);
            } // not classical

            for (int nrn_ = 0; nrn_ < n; ++nrn_) { int const nrn = nrn_; // smooth number or radial nodes
                int const iln = ln_off + nrn;
                auto & vs = partial_wave[iln]; // abbreviate ('vs' stands for valence state)

                // create the true partial wave
                set(vs.wave[TRU], nr, 0.0); // clear
                double normalize{1};
#ifdef DEVEL
                bool orthogonalize{false};
                bool const use_energy_derivative = ('D' == partial_wave_char[iln]);
                if (use_energy_derivative) {
                    // choose Psi_1 as energy derivative:
                    // solve inhomogeneous equation (H - E) Psi_1 = Psi_0
                    assert(iln > 0); assert(partial_wave[iln - 1].ell == ell);
                    vs.energy = partial_wave[iln - 1].energy; // same energy parameter as for Psi_0
                    std::vector<double> ff(rg[TRU].n), inh(rg[TRU].n);
                    product(inh.data(), rg[TRU].n, partial_wave[iln - 1].wave[TRU], rg[TRU].r);
                    double dg;
                    if (echo > 1) std::printf("# %s for ell=%i use energy derivative at E= %g %s\n", label, ell, vs.energy*eV, _eV);
                    radial_integrator::integrate_outwards<SRA>(vs.wave[TRU], ff.data(),
                        rg[TRU], potential[TRU].data(), ell, vs.energy, -1, &dg, inh.data());
                    // and T Psi_1 = (E - V) Psi_1 + Psi_0 for the kinetic part later
                    normalize = 0; // do not normalize
                    orthogonalize = true; // orthogonalize the higher true partial waves
                } else
#endif // DEVEL
                if (partial_wave_char[iln] >= '1' && partial_wave_char[iln] <= '9') { // regular quantum number
                    // if vs.energy matches an eigenenergy of the spherical potential, we can use this integrator and normalize
                    int nnodes{0};
                    // integrate SRA equation (^T + V(r) - E) phi(r) = 0 outwards homogeneously
                    radial_integrator::shoot(SRA, rg[TRU], potential[TRU].data(), ell, vs.energy, nnodes, vs.wave[TRU], r2rho.data());
                    normalize = 1; // normalize to unit charge
                } else {
                    view2D<double> ff(1, rg[TRU].n);
                    if (echo > 1) std::printf("# %s for the %s-partial wave integrate outwards at E= %g %s\n", label, vs.tag, vs.energy*eV, _eV);
                    radial_integrator::integrate_outwards<SRA>(vs.wave[TRU], ff[0], rg[TRU], potential[TRU].data(), ell, vs.energy);
                    product(r2rho.data(), rg[TRU].n, vs.wave[TRU], vs.wave[TRU]); // ToDo: maybe ff needs to be addded for a correct norm
                    normalize = 1; // dot_product(ir_cut[TRU], r2rho.data(), rg[TRU].dr); // normalize to have unit charge inside rcut
                } // valence eigenstate?

#ifdef DEVEL
                if (nrn > 0 && orthogonalize) {
                    auto const psi0 = partial_wave[iln - 1].wave[TRU];
                    auto const d  = dot_product(nr, vs.wave[TRU], psi0, rg[TRU].rdr);
                    auto const d2 = dot_product(nr, psi0, psi0, rg[TRU].r2dr);
                    auto const p  = -d/d2; // projection part
                    if (echo > 1) std::printf("# %s for ell=%i orthogonalize energy derivative with coefficient %g\n", label, ell, p);
                    add_product(vs.wave[TRU], nr, psi0, rg[TRU].r, p);
                } // orthogonalize
#endif // DEVEL
                { // scope: divide by r and normalize the partial waves
                    double norm_factor{1};
                    if (normalize > 0) {
                        auto const norm_wf2 = dot_product(nr, r2rho.data(), rg[TRU].dr);
                        if (norm_wf2 > 0) {
                            norm_factor = std::sqrt(normalize/norm_wf2);
                        } else {
                            warn("%s cannot normalize true partial %s-wave", label, vs.tag);
                            norm_factor = 0;
                        }
                    } // normalize
                    scale(vs.wave[TRU], nr, rg[TRU].rinv, norm_factor); // transform r*wave(r) as produced by the radial_eigensolver to wave(r)
                } // scope

                // create wKin for the computation of the kinetic energy density
                product(vs.wKin[TRU], nr, potential[TRU].data(), vs.wave[TRU], -1.); // start as wKin = -r*V(r)*wave(r)
                add_product(vs.wKin[TRU], nr, rg[TRU].r, vs.wave[TRU], vs.energy); // now wKin = r*(E - V(r))*wave(r)
#ifdef DEVEL
                if (use_energy_derivative && nrn > 0) {
                    add_product(vs.wKin[TRU], nr, rg[TRU].r, partial_wave[iln - 1].wave[TRU]); // now wKin = r*(E - V(r))*wave(r) + r*psi0
                } // kinetic energy wave in the case of energy derivative has an extra term
#endif // DEVEL
//                 auto const tru_norm = dot_product(ir_cut[TRU], r2rho.data(), rg[TRU].dr)/norm_wf2; // integrate only up to rcut
//                 auto const work = r2rho.data();
//                 scale(work, nr, potential[TRU].data());
//                 auto const tru_Epot = dot_product(ir_cut[TRU], work, rg[TRU].dr)/norm_wf2;
//                 auto const tru_kinetic_E = vs.energy*tru_norm - tru_Epot; // kinetic energy contribution up to r_cut

//                 if (echo > 1) std::printf("# valence %2d%c%6.1f E=%16.6f %s\n", vs.enn, ellchar[ell], vs.occupation, vs.energy*eV,_eV);
                show_state_analysis(echo - 5, label, rg[TRU], vs.wave[TRU], vs.tag, vs.occupation, vs.energy, csv_name(valence), ir_cut[TRU]);
                
                if (recreate_second == method && '*' == partial_wave_char[iln] && 1 == nrn) {
                    if (echo > 7) std::printf("\n# %s recreate_second for %s\n", label, vs.tag);
                    int const iln0 = ln_off + (nrn - 1); // index of the previous partial wave
                    // project the previous smooth partial wave iln0 onto the radial SHO basis
                    double co[8] = {0,0,0,0, 0,0,0,0};
                    double cp[8] = {0,0,0,0, 0,0,0,0};
                    assert(nmx < 8);
                    for (int irn = 0; irn < nmx; ++irn) {
                        co[irn] = dot_product(rg[SMT].n, partial_wave[iln0].wave[SMT], radial_sho_basis[ln_off + irn], rg[SMT].rdr);
                        cp[irn] = projector_coeff[ell](nrn - 1,irn); // lower projector
                    } // irn
                    // the new projector for the excited partial wave needs to have a projector
                    // orthogonal to the lower partial wave and linear independent of the lower projector
                    // lucky shot: orthogonalize the projector coefficients of the lower partial wave, cp against co
                    double const cocp = dot_product(nmx, co, cp);
                    double const coco = dot_product(nmx, co, co);
                    if (echo > 3) std::printf("# %s recreate_second for %s: inner products %g and %g\n", label, vs.tag, cocp, coco);
                    add_product(cp, 8, co, -cocp/coco); // orthogonalize
                    double const cocp_check = dot_product(8, co, cp);
                    if (echo > 7) std::printf("# %s recreate_second for %s: inner product after orthogonalization %.1e\n", label, vs.tag, cocp_check);
                    for (int irn = 0; irn < nmx; ++irn) {
                        projector_coeff[ell](nrn,irn) = cp[irn]; // store in member variable
                        if (echo > 7) std::printf("# %s recreate_second for %s coefficient #%i is %g\n", label, vs.tag, irn, cp[irn]);
                    } // irn

                    // reconstruct the second projector
                    set(projectors_ell[nrn], projectors_ell.stride(), 0.0);
                    for (int mrn = 0; mrn < nmx; ++mrn) { // radial SHO basis size
                        auto const c = projector_coeff[ell](nrn,mrn);
                        if (0.0 != c) {
                            add_product(projectors_ell[nrn], rg[SMT].n, radial_sho_basis[ln_off + mrn], c);
                            if (echo > 5) std::printf("# %s re-construct %s projector by taking %9.6f of the %c%i radial SHO basis function\n",
                                                         label, partial_wave[iln].tag, c, ellchar[ell],mrn);
                        } // coefficient non-zero
                    } // mrn
                    
                } // recreate_second

                if (classical_scheme == method
#ifdef DEVEL
                    || classical_partial_waves == method
#endif // DEVEL
                   ) {
                    bool const modify_projectors = (classical_scheme == method);

                    // Bloechl scheme
                    
                    // Caveat: for method == classical_partial_waves
                    // We can use the usual pseudized partial waves together with SHO-filtered projectors,
                    // however, the PAW equation is not fullfilled and we can observe small deviations in
                    // the logarithmic_derivative, e.g. for C -.1 eV in the s- and +.3 eV in the p-channel
                    // furthermore, eigenstate_analysis is not prepared for numerically given projectors

                    // copy the tail of the true wave function into the smooth wave function
                    set(vs.wave[SMT], rg[SMT].n, vs.wave[TRU] + nr_diff);
                    set(vs.wKin[SMT], rg[SMT].n, vs.wKin[TRU] + nr_diff);
                    if (modify_projectors) {
                        set(projectors_ell[nrn], projectors_ell.stride(), 0.0); // clear
                    } // modify_projectors

                    // pseudize the true partial wave by matching a polynomial r^ell*(c_0 + c_1 r^2 + c_2 r^4 + c_3 r^6)
                    double coeff[4]; // matching polynomial coefficients
                    auto const stat = pseudo_tools::pseudize_function(vs.wave[SMT], rg[SMT].r, ir_cut[SMT], 4, ell, coeff);
                    // ToDo: ir_cut could be ell-dependent here
                    assert(0 == stat);

                    // construct a preliminary projector according to the Bloechl scheme:
                    //    ~p(r) = (^T + ~V - E) ~phi(r)
                    double ckin[3]; // polynomial coefficients of the kinetic operator times the smooth wave
                    // The kinetic energy operator is
                    //    ^T = .5*( ell(ell+1)/r^2 - d^2/dr^2 ) when acting onto r*wave(r)
                    for (int k = 1; k < 4; ++k) {
                        ckin[k - 1] = 0.5*( ell*(ell + 1) - (ell + 2*k)*(ell + 1 + 2*k) )*coeff[k];
                    } // k

                    if (echo > 19) std::printf("\n## %s classical method for %c%i: r, smooth r*wave, smooth r*Twave, "
                                          "true r*wave, true r*Twave, projector:\n", label, ellchar[ell], nrn);
                    // expand smooth kinetic wave up to r_cut
                    for (int ir = 0; ir <= ir_cut[SMT]; ++ir) {
                        double const r = rg[SMT].r[ir], r2 = pow2(r), rl1 = intpow(r, ell + 1);
 
                        vs.wKin[SMT][ir] = rl1*(ckin[0] + r2*(ckin[1] + r2*ckin[2])); // expand polynomial ...
                        // ... describing the kinetic energy operator times the smooth wave function times r

                        if (modify_projectors) {
                            // Here, we construct the preliminary projector functions accoding to the Bloechl scheme:
                            //    ~p(r) := (^T + V_loc - E) ~phi(r)
                            //    i.e. the projector functions are only given numerically and are non-zero only inside the sphere
                            projectors_ell(nrn,ir) = vs.wKin[SMT][ir] + (potential[SMT][ir] - r*vs.energy)*vs.wave[SMT][ir];
                        } // modify_projectors

                        if (echo > 19) std::printf("%g %g %g %g %g %g\n", r, vs.wave[SMT][ir], vs.wKin[SMT][ir],
                            vs.wave[TRU][ir + nr_diff], vs.wKin[TRU][ir + nr_diff], projectors_ell(nrn,ir)*rg[SMT].rinv[ir]);
                    } // ir

                    if (echo > 19) {
                        // beyond the cutoff radius show only SMT since smooth and true are identical by construction
                        std::printf("\n## %s classical method for %c%i: r, smooth/true r*wave, smooth/true r*Twave, projector:\n", label, ellchar[ell], nrn);
                        for (int ir = ir_cut[SMT]; ir < rg[SMT].n; ++ir) {
                            std::printf("%g %g %g %g\n", rg[SMT].r[ir], vs.wave[SMT][ir], vs.wKin[SMT][ir], projectors_ell(nrn,ir)*rg[SMT].rinv[ir]);
                        } // ir
                        std::printf("\n\n");
                    } // echo

                    // Beware: the numerical projectors generated are only preliminary and not yet SHO filtered!
                    // please run the sigma optimization after using this method.

                } // classical_scheme or classical_partial_waves        
                else
                { // scope: generate smooth partial waves from projectors, revPAW scheme
                    assert(ir_match[TRU] > 0);
                    int const stride = align<2>(rg[SMT].n);
                    view2D<double> rphi(1 + n, stride);
                    view2D<double> Tphi(1 + n, stride);
                    std::vector<double> ff(rg[SMT].n);
                    auto const vgtru = vs.wave[TRU][ir_match[TRU]]    *rg[TRU].r[ir_match[TRU]];
                    auto const dgtru = vs.wave[TRU][ir_match[TRU] - 1]*rg[TRU].r[ir_match[TRU] - 1] - vgtru;
                    double vghom{1}, dghom{1};

                    int constexpr HOM = 0; // homogeneous solution index

                    for (int krn = HOM; krn <= n; ++krn) { // loop must run serial and forward
                        // krn == 0 generates the homogeneous solution in the first iteration
                        auto projector = (krn > HOM) ? projectors_ell[krn-1] : nullptr;
#ifdef DEVEL
                        std::vector<double> rhs;
                        if (use_energy_derivative && nrn > 0 && krn > HOM) {
                            rhs = std::vector<double>(stride);
                            product(rhs.data(), rg[SMT].n, rg[SMT].r, partial_wave[iln - 1].wave[SMT]); // r*phi_0
                            projector = rhs.data(); // use as inhomogeneity
                            // comment: this will be 2x the same solution, for krn==1 and krn==2, code results from a quick fix
                        }
#endif // DEVEL

                        double dg; // derivative at end point
                        // solve inhomgeneous equation and match true wave in value and derivative
                        radial_integrator::integrate_outwards<SRA>(rphi[krn], ff.data(),
                            rg[SMT], potential[SMT].data(), ell, vs.energy, -1, &dg, projector);

                        auto const vginh = rphi(krn,ir_match[SMT]);
                        auto const dginh = rphi(krn,ir_match[SMT] - 1) - vginh;
                        if (HOM == krn) {
                            vghom = vginh; dghom = dginh;
                        } else {
                            // matching coefficient - how much of the homogeneous solution do we need to add to match logder
                            auto const denom = vgtru*dghom - dgtru*vghom; // ToDo: check if denom is not too small
#ifdef DEVEL
                            if (echo > 17) { std::printf("# %s LINE= %d    %g * %g - %g * %g = %g\n", label,
                                             __LINE__, vgtru, dghom, dgtru, vghom, denom); std::fflush(stdout); }
#endif // DEVEL
                            auto const inv_denom = (std::abs(denom) > 1e-300) ? 1./denom : 0;
                            auto const c_hom = - (vgtru*dginh - dgtru*vginh) * inv_denom;

                            // combine inhomogeneous and homogeneous solution to match true solution outside r_match
                            add_product(rphi[krn], rg[SMT].n, rphi[HOM], c_hom);

                            auto const denom2 = vginh + c_hom*vghom;
                            auto const scal = (std::abs(denom2) > 1e-300) ? vgtru/denom2 : 0; // == vgtru/(vginh + c_hom*vghom) match not only in logder but also in value
                            assert(scal == scal);
                            scale(rphi[krn], rg[SMT].n, rg[SMT].rinv, scal); // scale and divide by r
                            // ToDo: extrapolate lowest point?

#ifdef DEVEL
                            if (use_energy_derivative && nrn > 0) {
                                // (T + V - E) phi_0 = linear combination of projectors
                                // (T + V - E) phi_1 = phi_0 
                                // -. T phi_1 = (E - V) phi_1 + phi_0
                                for (int ir = 0; ir < rg[SMT].n; ++ir) {
                                    Tphi(krn,ir) = (vs.energy*rg[SMT].r[ir] - potential[SMT][ir])*rphi(krn,ir) + scal*rhs[ir];
                                } // ir
                                // ToDo: check these equations and normalization factors
                                if (echo > 1) std::printf("# %s generate Tphi with inhomogeneity\n", label);
                                // seems like the tails of TRU and SMT wave and wKin are deviating slightly beyond r_match

                            } else 
#endif // DEVEL
                            { // scope: kinetic energy operator times wave
                                // Tphi = (E - V)*phi + projector
                                for (int ir = 0; ir < rg[SMT].n; ++ir) {
                                    Tphi(krn,ir) = (vs.energy*rg[SMT].r[ir] - potential[SMT][ir])*rphi(krn,ir) + scal*projector[ir];
                                } // ir
                            } // scope

#ifdef DEVEL
                            // now visually check that the matching of value and derivative of rphi is ok.
                            if (echo > 29) {
                                std::printf("\n## %s check matching of rphi for ell=%i nrn=%i krn=%i (r, phi_tru,phi_smt, prj, rTphi_tru,rTphi_smt):\n",
                                                  label, ell, nrn, krn-1);
                                for (int ir = 1; ir < rg[SMT].n; ++ir) {
                                    std::printf("%g  %g %g  %g  %g %g\n", rg[SMT].r[ir],     vs.wave[TRU][ir + nr_diff], rphi(krn,ir),
                                                 projectors_ell(krn-1,ir)*rg[SMT].rinv[ir],  vs.wKin[TRU][ir + nr_diff], Tphi(krn,ir));
                                } // ir
                                std::printf("\n\n");
                            } // echo

                            // check that the matching of value and derivative of rphi is ok by comparing value and derivative
                            if (echo > 9) {
                                std::printf("# %s check matching of vg and dg for ell=%i nrn=%i krn=%i: %g == %g ? and %g == %g ?\n",
                                    label, ell, nrn, krn-1, vgtru, scal*(vginh + c_hom*vghom), dgtru, scal*(dginh + c_hom*dghom));
                            } // echo
#endif // DEVEL
                        } // krn > 0
                    } // krn


                    std::vector<double> evec(n, 0.0);
                    if (recreate_second == method) {
                        evec[0] = 1; // for method==recreate_second
                    } else {
                        evec[nrn] = 1.0; // e.g. for method==energy_ordering
                    }

                    if (n > 1) {
#ifdef DEVEL
                        if (minimize_radial_curvature == method) {
                            view2D<double> Ekin(     3*n , n);
                            view2D<double> Olap(Ekin[1*n], n);
                            view2D<double> Vkin(Ekin[2*n], n);
                            int const nr_max = ir_match[SMT];
                            for (int krn = 0; krn < n; ++krn) {
                                for (int jrn = 0; jrn < n; ++jrn) {
                                    // definition of Tphi contains factor *r
                                    Ekin(krn,jrn) = dot_product(nr_max, Tphi[1 + krn], rphi[1 + jrn], rg[SMT].rdr);
                                    Olap(krn,jrn) = dot_product(nr_max, rphi[1 + krn], rphi[1 + jrn], rg[SMT].r2dr);
                                    Vkin(krn,jrn) = dot_product(nr_max, rphi[1 + krn], rphi[1 + jrn], rg[SMT].dr)
                                                    * 0.5*(ell*(ell + 1)); // kinetic energy in angular direction, Hartree units
                                    // minimize only the radial kinetic energy
                                    Ekin(krn,jrn) -= Vkin(krn,jrn); // subtract the angular part, suggested by Morian Sonnet
                                } // jrn
                            } // krn

                            // symmetrize
                            for (int krn = 0; krn < n; ++krn) {
                                for (int jrn = 0; jrn < krn; ++jrn) { // triangular loop
                                    symmetrize(Ekin(krn,jrn), Ekin(jrn,krn));
                                    symmetrize(Olap(krn,jrn), Olap(jrn,krn));
                                } // jrn
                            } // krn

                            if (echo > 7) {
    //                          std::printf("\n");
                                for (int krn = 0; krn < n; ++krn) {
                                    std::printf("# %s curvature (%s) and overlap for i=%i ", label, _eV, krn);
                                    printf_vector(" %g", Ekin[krn], n, "\t\t", eV);
                                    printf_vector(" %g", Olap[krn], n);
                                } // krn
    //                          std::printf("\n");
                            } // echo

                            { // scope: minimize the radial curvature of the smooth partial wave
                                double lowest_eigenvalue{0};
                                auto const info = minimize_curvature(n, Ekin, Olap, &lowest_eigenvalue);
                                if (info) {
                                    warn("%s generalized eigenvalue problem in minimize_curvature failed, info= %i", label, int(info));
                                } else {
                                    set(evec.data(), n, Ekin[0]);
                                    if (echo > 6) {
                                        std::printf("# %s lowest eigenvalue of the radial curvature is %g %s", label, lowest_eigenvalue*eV, _eV);
                                        if (echo > 8) { std::printf(", coefficients"); printf_vector(" %g", Ekin[0], n, nullptr); }
                                        std::printf("\n");
                                    } // echo
                                } // info
                            } // scope

                        } else // method
                        if (orthogonalize_second == method) {
                         
                            if (nrn > 0) { // only the second
                                assert(1 == nrn); // we do not know how to treat a 3rd wave yet
                                
                                // minimize the matrix element <Psi_1|p_0>
                                double c[2];
                                for (int krn = 0; krn < 2; ++krn) {
                                    c[krn] = dot_product(rg[SMT].n, rphi[1 + krn], projectors_ell[0], rg[SMT].rdr);
                                } // krn
                                // now find the angle phi such that cos(phi)*c[1] + sin(phi)*c[0] is zero;
#if 0
                                for (int iang = -18; iang <= 18; ++iang) {
                                    double const angle = iang * constants::pi / 18;
                                    double const ovl10 = std::cos(angle)*c[1] + std::sin(angle)*c[0];
                                    std::printf("# method=orthogonalize_second angle=%g\t<Psi_1|p_0>= %g\n", angle, ovl10);
                                } // iang
#endif // 0
                                double const angle = std::atan2(-c[1], c[0]);
                                evec[0] = std::sin(angle);
                                evec[1] = std::cos(angle);
                                if (echo > 8) {
                                    auto const ovl10 = evec[0]*c[0] + evec[1]*c[1];
                                    std::printf("# %s method=orthogonalize_second angle=%g\t<Psi_1|p_0>= %g coeffs= %g %g\n",
                                                   label, angle, ovl10, evec[0], evec[1]);
                                } // echo

                            } // if nrn > 0

                        } else // method
                        if (orthogonalize_first == method) {

                            if (0 == nrn) { // only the first

                                // minimize the matrix element <Psi_0|p_1>
                                double c[2];
                                for (int krn = 0; krn < 2; ++krn) {
                                    c[krn] = dot_product(rg[SMT].n, rphi[1 + krn], projectors_ell[1], rg[SMT].rdr);
                                } // krn
                                // now find the angle phi such that cos(phi)*c[1] + sin(phi)*c[0] is zero;
                                double const angle = std::atan2(-c[1], c[0]);
                                evec[0] = std::sin(angle);
                                evec[1] = std::cos(angle);
                                if (echo > 1) {
                                    auto const ovl01 = evec[0]*c[0] + evec[1]*c[1];
                                    std::printf("# %s method=orthogonalize_first angle=%g\t<Psi_0|p_1>= %g coeffs= %g %g\n", 
                                                   label, angle, ovl01, evec[0], evec[1]);
                                } // echo

                            } // if nrn > 0

                        } else // method

                        if (energy_ordering == method) {
                            assert(1 == evec[nrn]);
                        } // energy_ordering
#endif // DEVEL
                        if (recreate_second == method) {
                            assert(1 == evec[0]);
                        } // recreate_second

                    } // n > 1 more than one partial wave

                    set(vs.wave[SMT], rg[SMT].n, 0.0);
                    set(vs.wKin[SMT], rg[SMT].n, 0.0);
                    double sum{0}; for (int krn = 0; krn < n; ++krn) { sum += evec[krn]; }
                    assert(std::abs(sum) > 1e-300);
                    double const scal = 1.0/sum; // scaling such that the sum is 1
                    for (int krn = 0; krn < n; ++krn) {
                        auto const coeff = scal*evec[krn];
                        add_product(vs.wave[SMT], rg[SMT].n, rphi[1 + krn], coeff);
                        add_product(vs.wKin[SMT], rg[SMT].n, Tphi[1 + krn], coeff);
                    } // krn
#ifdef DEVEL
                    if (echo > 19) {
                        std::printf("\n## %s check matching of partial waves ell=%i nrn=%i (r, phi_tru,phi_smt, rTphi_tru,rTphi_smt):\n",
                                          label, ell, nrn);
                        for (int ir = 0; ir < rg[SMT].n; ++ir) {
                            std::printf("%g  %g %g  %g %g\n",  rg[SMT].r[ir],
                                vs.wave[TRU][ir + nr_diff], vs.wave[SMT][ir],
                                vs.wKin[TRU][ir + nr_diff], vs.wKin[SMT][ir]);
                        } // ir
                        std::printf("\n\n");
                    } // echo
#endif // DEVEL
                } // not classical, revPAW

            } // nrn


            { // scope: establish dual orthgonality with projectors
#ifdef DEVEL
                if (echo > 24) { // show normalization and orthogonality of projectors
                    for (int irn = 0; irn < n; ++irn) {
                        for (int jrn = 0; jrn < n; ++jrn) {
                            int const diag = (irn == jrn);
                            std::printf("# %s %c-projector <#%d|#%d> = %i + %.1e sigma=%g %s\n", label, ellchar[ell], irn, jrn, diag,
                                dot_product(rg[SMT].n, projectors_ell[irn], projectors_ell[jrn], rg[SMT].dr) - diag, sigma*Ang, _Ang);
                        } // jrn
                    } // irn
                    std::printf("# %s Mind: unlike in previous versions, projectors are not necessarily orthonormalized!\n", label);
                } // echo
#endif // DEVEL
                view2D<double> ovl(n, n); // get memory
                view3D<double> LU_inv(2, n, n, 0.0); // get memory

                for (int iGS = 0; iGS < Gram_Schmidt_iterations; ++iGS) {
                    int const echo_GS = echo - 4*iGS; // most output in first iteration

                    for     (int irn = 0; irn < n; ++irn) { // number of partial waves
                        for (int jrn = 0; jrn < n; ++jrn) { // number of partial waves
                            ovl(irn,jrn) = dot_product(rg[SMT].n, projectors_ell[irn], partial_wave[ln_off + jrn].wave[SMT], rg[SMT].rdr);
                            if (echo > 4) std::printf("# %s %c-projector #%i with partial wave #%i has overlap %g\n",
                                                         label, ellchar[ell], irn, jrn, ovl(irn,jrn));
                        } // jrn
                    } // irn

                    pseudo_tools::perform_Gram_Schmidt(n, LU_inv, ovl, label, ellchar[ell], echo_GS);

                    { // scope: orthogonalize the projectors
                        int const mr = projectors_ell.stride();
                        view2D<double> proj(n, mr, 0.0); // temporary storage for projectors on radial grid
                        int const nmx = sho_tools::nn_max(numax, ell); // number of radial SHO basis functions
                        view2D<double> proj_coeff(n, nmx, 0.0);
                        for (int irn = 0; irn < n; ++irn) { // number of partial waves
                            set(proj[irn], mr, projectors_ell[irn]); // copy into temporary
                            set(proj_coeff[irn], nmx, projector_coeff[ell][irn]); // copy into temporary
                        } // irn
                        // reconstruct projectors
                        for (int irn = 0; irn < n; ++irn) {
                            set(projectors_ell[irn], mr, 0.0); // clear
                            set(projector_coeff[ell][irn], nmx, 0.0); // clear
                            for (int jrn = 0; jrn < n; ++jrn) {
                                auto const p_coeff = LU_inv(0,irn,jrn);
                                if (p_coeff != 0) {
                                    if (echo_GS > 3) {
                                        std::printf("# %s create orthogonalized %c-projector #%i with %g * projector #%i\n",
                                                       label, ellchar[ell], irn, p_coeff, jrn);
                                    } // echo
                                    add_product(projectors_ell[irn], mr, proj[jrn], p_coeff);
                                    add_product(projector_coeff[ell][irn], nmx, proj_coeff[jrn], p_coeff);
                                } // coeff nonzero
                            } // jrn
                        } // irn
                    } // scope 
                    
                    // make a new linear combination of the partial waves
                    for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                        int const nr = rg[ts].n;
                        view3D<double> waves(2, n, align<2>(nr), 0.0); // temporary storage for pairs {wave, wKin}
                        for (int irn = 0; irn < n; ++irn) {
                            int const iln = ln_off + irn;
                            set(waves(0,irn), nr, partial_wave[iln].wave[ts]); // copy into temporary
                            set(waves(1,irn), nr, partial_wave[iln].wKin[ts]); // copy into temporary
                        } // irn
                        // reconstruct partial waves
                        for (int irn = 0; irn < n; ++irn) {
                            int const iln = ln_off + irn;
                            set(partial_wave[iln].wave[ts], nr, 0.0); // clear
                            set(partial_wave[iln].wKin[ts], nr, 0.0); // clear
                            for (int jrn = 0; jrn < n; ++jrn) {
                                auto const w_coeff = LU_inv(1,jrn,irn);
                                if (w_coeff != 0) {
                                    if (ts == TRU && echo_GS > 3) {
                                        std::printf("# %s create orthogonalized partial %c-wave #%i with %g * wave #%i\n",
                                                       label, ellchar[ell], irn, w_coeff, jrn);
                                    } // echo
                                    add_product(partial_wave[iln].wave[ts], nr, waves(0,jrn), w_coeff);
                                    add_product(partial_wave[iln].wKin[ts], nr, waves(1,jrn), w_coeff);
                                } // coeff nonzero
                            } // jrn
                        } // irn
                    } // ts in {tru, smt}

                } // iGS

#ifdef DEVEL
                if (n > 0) { // scope: check the overlap again
                    view2D<double> ovl_new(n, n); // get memory
                    double dev{0}; // largest deviations from unity
                    for     (int irn = 0; irn < n; ++irn) { // smooth number or radial nodes
                        for (int jrn = 0; jrn < n; ++jrn) { // smooth number or radial nodes
                            int const jln = ln_off + jrn;
                            ovl_new(irn,jrn) = dot_product(rg[SMT].n, projectors_ell[irn], partial_wave[jln].wave[SMT], rg[SMT].rdr);
                            int const Kronecker = (irn == jrn);
                            auto const deviate = ovl_new(irn,jrn) - Kronecker;
                            if (echo > 7) std::printf("# %s %c-projector #%d with partial %c-wave #%d with new overlap= %i + %g\n",
                                                         label, ellchar[ell], irn, ellchar[ell], jrn, Kronecker, deviate);
                            dev = std::max(dev, std::abs(deviate));
                        } // jrn
                    } // irn
                    if (echo > 2) std::printf("# %s after %dx orthogonalization %c-<projectors|partial waves> deviates max. %.1e from unity matrix\n",
                                                 label, Gram_Schmidt_iterations, ellchar[ell], dev);
                    if (dev > 1e-12) {
                        warn("%s %c-duality violated, deviates %g from unity", label, ellchar[ell], dev);
                        if (echo > 0) {
                            std::printf("\n# %s %c-<projectors|partial waves> matrix deviates %.1e:\n", label, ellchar[ell], dev);
                            for (int i = 0; i < n; ++i) {
                                std::printf("# %s irn=%2i ", label, i);
                                printf_vector(" %11.6f", ovl_new[i], n);
                            } // i
                            std::printf("\n");
                        } // echo
                    } // dev
                } // scope

                if (echo > 15) {
                    for (int irn = 0; irn < n; ++irn) { // smooth number or radial nodes
                        auto const & vs = partial_wave[ln_off + irn];
                        std::printf("\n## %s show orthogonalized partial %c-waves for irn=%i (r, phi_tru, phi_smt, rTphi_tru, rTphi_smt):\n",
                                          label, ellchar[ell], irn);
                        for (int ir = 1; ir < rg[SMT].n; ++ir) {
                            std::printf("%g  %g %g  %g %g\n",  rg[SMT].r[ir],
                                vs.wave[TRU][ir + nr_diff], vs.wave[SMT][ir],
                                vs.wKin[TRU][ir + nr_diff], vs.wKin[SMT][ir]);
                        } // ir
                        std::printf("\n\n");
                    } // irn
                } // echo
#endif // DEVEL

            } // scope: establish dual orthogonality with [SHO] projectors

        } // ell

        if (method != classical_scheme) {
            show_projector_coefficients("orthogonalized ", sigma, echo - 6); // and warn if not normalizable
        } // not classical

#ifdef DEVEL
        char const *export_as_json = control::get("single_atom.export.json", "");
        if ((nullptr != export_as_json) && ('\0' != *export_as_json)) {

            std::ofstream json(export_as_json, std::ios::out);
            json << "{"; // begin json
            json <<  "\n\t\"atomic number\": " << Z_core;
            json << ",\n\t\"length unit\": \"Bohr\"";
            json << ",\n\t\"energy unit\": \"Hartree\"";

            json << ",\n\t\"partial waves\": ";
            for (int i = 0, ell = 0; ell <= numax; ++ell) {
                int const ln_off = sho_tools::ln_index(numax, ell, 0); // offset where to start indexing emm-degenerate partial waves
                for (int nrn = 0; nrn < nn[ell]; ++nrn) { // smooth number or radial nodes
                    auto const & vs = partial_wave[ln_off + nrn]; // abbreviate ('vs' stands for valence state)
                    if (vs.ell != ell) warn("%s inconsistency in partial waves!", label);
                    
                    json << (i?',':'{') << "\n\t\t\"" << int(vs.enn) << ellchar[vs.ell] << "\": {"; // begin state
                    json <<  "\n\t\t\t\"energy\": " << vs.energy;
                    json << ",\n\t\t\t\"occupation\": " << vs.occupation;
                    for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                        int const ir_off = (TRU == ts)*nr_diff;
                        json << ",\n\t\t\t\"" << ts_name[ts] << " wave\": ";
                        for (int ir = 0; ir < rg[SMT].n; ++ir) {
                            json << (ir?((ir&7)?", ":",\n\t\t\t\t"):"[") << vs.wave[ts][ir + ir_off];
                        } // ir
                        json << "\n\t\t\t]"; // end wave
                    } // ts
                    json << "\n\t\t}"; // end state
                    ++i;
                } // nrn
            } // ell
            json << "\n\t}"; // end partial waves

            for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                int const ir_off = (TRU == ts)*nr_diff;
                json << ",\n\t\"r*" << ts_name[ts] << " potential\": ";
                for (int ir = 0; ir < rg[SMT].n; ++ir) {
                    json << (ir?((ir&7)?", ":",\n\t\t"):"[") << potential[ts][ir + ir_off];
                } // ir
                json << "\n\t]"; // end r*potential
            } // ts

            json << ",\n\t\"radial grid\": {";
            json <<  "\n\t\t\"number\": " << rg[SMT].n - 1;
            json << ",\n\t\t\"values\": ";
            for (int ir = 0; ir < rg[SMT].n; ++ir) {
                json << (ir?((ir&7)?", ":",\n\t\t\t"):"[") << rg[SMT].r[ir];
            } // ir
            json << "\n\t\t]"; // end radial_grid.values
            json << "\n\t}"; // end radial_grid

            json << "\n}"; // end json

        } // export_as_json
#endif // DEVEL        

    } // update_partial_waves

    
    void update_kinetic_energy_deficit(int const echo=0) {
        
        for (int ell = 0; ell <= numax; ++ell) {
            int const ln_off = sho_tools::ln_index(numax, ell, 0); // offset where to start indexing emm-degenerate partial waves
            int const n = nn[ell]; // abbreviation

            // compute kinetic energy difference matrix from wKin
            for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                int const nr_cut = rg[ts].n; // integration over the entire grid, diagonal elements then appear positive
                for     (int iln = 0 + ln_off; iln < n + ln_off; ++iln) {
                    for (int jln = 0 + ln_off; jln < n + ln_off; ++jln) {
                        kinetic_energy(ts,iln,jln) = dot_product(nr_cut,
                            partial_wave[iln].wKin[ts],
                            partial_wave[jln].wave[ts], rg[ts].rdr); // we only need rdr here since wKin is defined as r*(E - V(r))*wave(r)
                    } // j
                } // i
            } // ts

#ifdef DEVEL       
            // display
            for     (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    auto const E_kin_tru = kinetic_energy(TRU,i+ln_off,j+ln_off);
                    auto const E_kin_smt = kinetic_energy(SMT,i+ln_off,j+ln_off);
                    if (echo > 19) std::printf("# %s %c-channel <%d|T|%d> kinetic energy [unsymmetrized] (true) %g and (smooth) %g (diff) %g %s\n",
                        label, ellchar[ell], i, j, E_kin_tru*eV, E_kin_smt*eV, (E_kin_tru - E_kin_smt)*eV, _eV);
                } // j
            } // i
#endif // DEVEL

            if (1) { // symmetrize the kinetic energy tensor
                for (int iln = 0 + ln_off; iln < n + ln_off; ++iln) {
                    for (int jln = 0 + ln_off; jln < iln; ++jln) { // triangular loop excluding the diagonal elements
                        if (1) { // usual
                            for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                                symmetrize(kinetic_energy(ts,iln,jln), kinetic_energy(ts,jln,iln));
                            } // ts
                        } else {
                            // actually we do not need to symmetrize each contribution kinetic_energy[ts]
                            // but it would be enough to symmetrize the difference matrix kinetic_energy[TRU] - kinetic_energy[SMT]
                            // -. symmetrize only the difference
                            for (int twice = 0; twice < 2; ++twice) { // do this twice for numerical accuracy
                                auto const dij = kinetic_energy(TRU,iln,jln) - kinetic_energy(SMT,iln,jln);
                                auto const dji = kinetic_energy(TRU,jln,iln) - kinetic_energy(SMT,jln,iln);
                                auto const asy = 0.5*(dij - dji); // asymmetry
                                for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                                    int const sign = (TRU == ts) ? -1 : 1;
                                    kinetic_energy(ts,iln,jln) += 0.5*sign*asy;
                                    kinetic_energy(ts,jln,iln) -= 0.5*sign*asy;
                                } // ts
                            } // twice
                            // however, it should at the end not make a difference since only the difference enters the Hamiltonian
                            // and the energy contributions are computed with the density matrix which is symmetric.
                        } // symmetrize only the difference?
                    } // j
                } // i
            } // symmetrization active?

#ifdef DEVEL                
            if (echo > 19) { // display
                for     (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        auto const E_kin_tru = kinetic_energy(TRU,i+ln_off,j+ln_off);
                        auto const E_kin_smt = kinetic_energy(SMT,i+ln_off,j+ln_off);
                        std::printf("# %s %c-channel <%d|T|%d> kinetic energy [symmetrized] (true) %g and (smooth) %g (diff) %g %s\n",
                            label, ellchar[ell], i, j, E_kin_tru*eV, E_kin_smt*eV, (E_kin_tru - E_kin_smt)*eV, _eV);
                    } // j
                } // i
            } // echo
#endif // DEVEL

        } // ell
        
    } // update_kinetic_energy_deficit


    void update_charge_deficit(int const echo=0) {
        int const nln = sho_tools::nSHO_radial(numax);

        for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            int const nr = rg[ts].n; // integrate over the full radial grid
            std::vector<double> rl(nr, 1.0); // init as r^0
            std::vector<double> wave_r2rl_dr(nr);
            if (echo > 4) std::printf("\n# %s charges for %s partial waves\n", label, ts_name[ts]);
            for (int ell = 0; ell <= ellmax_cmp; ++ell) { // loop-carried dependency on rl, run forward, run serial!
                bool const echo_l = (echo > 4 + 4*(ell > 0));
                if (echo_l) std::printf("# %s charges for ell=%i\n", label, ell);
                if (ell > 0) scale(rl.data(), nr, rg[ts].r); // create r^{\ell}
                for (int iln = 0; iln < nln; ++iln) {
                    if (partial_wave_active[iln]) {
                        auto const *const wave_i = partial_wave[iln].wave[ts]; assert(wave_i);
                        if (echo_l) std::printf("# %s %s %-4s", label, ts?"smt":"tru", partial_wave[iln].tag);
                        product(wave_r2rl_dr.data(), nr, wave_i, rl.data(), rg[ts].r2dr); // product of three arrays
                        for (int jln = 0; jln < nln; ++jln) {
                            if (partial_wave_active[jln]) {
                                auto const *const wave_j = partial_wave[jln].wave[ts]; assert(wave_j);
                                auto const cd = dot_product(nr, wave_r2rl_dr.data(), wave_j);
                                charge_deficit(ell,ts,iln,jln) = cd;
                                if (echo_l) std::printf("\t%10.6f", cd);
//                              if (SMT == ts && echo > 1) std::printf("\t%10.6f", charge_deficit(ell,TRU,iln,jln) - cd);
                            } // active j
                        } // jln
                        if (echo_l) std::printf("\n");
                    } // active i
                } // iln
                if (echo_l) std::printf("\n");
            } // ell
        } // ts
    } // update_charge_deficit


    void transform_SHO(
          double out[] // assumed shape [N][out_stride]
        , int const out_stride
        , double const in[] // assumed shape [N][in_stride]
        , int const in_stride
        , bool const in_Cartesian // true: Cartesian to Radial, false: Radial to Cartesian
    ) const {

        int const N = unitary_zyx_lmn.stride(); // assumed shape [N][N]
        view2D<double const> const inp(in, in_stride); // wrap input
        view2D<double> const uni(unitary_zyx_lmn.data(), N); // wrap transformation matrix
        view2D<double> res(out, out_stride); // wrap
        view2D<double> tmp(N, N); // get memory

        auto const uni_transposed = transpose(uni, N);

        if (in_Cartesian) {
            // transform Cartesian input to Radial output

            gemm(tmp, N, inp, N, uni); // *u
            gemm(res, N, uni_transposed, N, tmp); // u^T*

        } else { // in_Cartesian
            // transform Radial input to Cartesian output

            gemm(tmp, N, uni, N, inp); // u*
            gemm(res, N, tmp, N, uni_transposed); // *u^T

        } // in_Cartesian

    } // transform_SHO




    view2D<double> unfold_projector_coefficients(
          emm_QN_t const m=-1
    ) const {

        int const nln = sho_tools::nSHO_radial(numax);
        std::vector<int16_t> ln_offset(nln, -1);
        std::vector<int8_t>   ell_list(nln, -1);
        for (int ell = 0; ell <= numax; ++ell) {
            for (int nrn = 0; nrn < sho_tools::nn_max(numax, ell); ++nrn) {
                int const iln  = sho_tools::ln_index(numax, ell, nrn);
                ln_offset[iln] = sho_tools::ln_index(numax, ell, 0);
                ell_list[iln] = ell;
            } // nrn
        } // ell
        for (int iln = 0; iln < nln; ++iln) { 
            assert(ln_offset[iln] > -1); // assert coverage
            assert(ell_list[iln]  > -1); // assert coverage
        } // iln

        if (emm_Degenerate == m) { 
            view2D<double> u_proj(nln, nln, 0.0);
            for (int iln = 0; iln < nln; ++iln) {
                int const ell = ell_list[iln];
                int const nrn = iln - ln_offset[iln];
                if (nrn < nn[ell]) {
                    assert(nrn >= 0);
                    for (int jln = 0; jln < nln; ++jln) {
                        if (ell == ell_list[jln]) {
                            int const mrn = jln - ln_offset[jln];
                            assert(mrn >= 0);
                            u_proj(iln,jln) = projector_coeff[ell](nrn,mrn);
                        } // ell-diagonal
                    } // jln
                } // active
            } // iln
            return u_proj;
        } // emm_Degenerate

        int const nlmn = sho_tools::nSHO(numax);
        view2D<double> u_proj(nlmn, nlmn, 0.0);
        for (int ilmn = 0; ilmn < nlmn; ++ilmn) {
            int const iln = ln_index_list[ilmn];
            int const ilm = lm_index_list[ilmn];
            int const ell = ell_list[iln];
            int const nrn = iln - ln_offset[iln];
            if (nrn < nn[ell]) {
                assert(nrn >= 0);
                for (int jlmn = 0; jlmn < nlmn; ++jlmn) {
                    int const jln = ln_index_list[jlmn];
                    if (ilm == lm_index_list[jlmn]) {
                        int const mrn = jln - ln_offset[jln];
                        assert(mrn >= 0);
                        u_proj(ilmn,jlmn) = projector_coeff[ell](nrn,mrn);
                    } // ell-diagonal
                } // jlmn
            } // active
        } // ilmn

        return u_proj;
    } // unfold_projector_coefficients

    
    void create_synthetic_density_matrix(view2D<double> & radial_density_matrix, int const echo=0) {
        if (echo > 1) std::printf("\n# %s create a synthetic radial density matrix\n", label);
        int const nSHO = sho_tools::nSHO(numax);
        for (int ilmn = 0; ilmn < nSHO; ++ilmn) {
            int const iln = ln_index_list[ilmn];
            if (partial_wave_active[iln]) {
                auto   const ell = partial_wave[iln].ell;
                double const occ = partial_wave[iln].occupation/(2*ell + 1.);
                for (int jlmn = 0; jlmn < nSHO; ++jlmn) {
                    int const jln = ln_index_list[jlmn];
                    if (partial_wave_active[jln]) {
                        radial_density_matrix(ilmn,jlmn) = (ilmn == jlmn)*occ;
                    } // active_j
                } // jlmn
            } // active_i
        } // ilmn
    } // create_synthetic_density_matrix


    void get_rho_tensor(
          view3D<double> & rho_tensor // result tensor with indices (lm,iln,jln) where iln,jln are in radial SHO basis
        , view2D<double> const & rho_matrix
        , sho_tools::SHO_order_t const order=sho_tools::order_zyx
        , int const echo=0
        , bool const synthetic_density_matrix=false
    ) {
        // from a density matrix (various representations possible) create a density tensor

        int const nSHO = sho_tools::nSHO(numax);

        int const lmax = std::max(ellmax_rho, ellmax_cmp);
        int const nlm = pow2(1 + lmax);  // limit for L=(ell,emm) for the density
        int const mlm = pow2(1 + numax); // limit for L=(ell,emm) of the partial waves

        view2D<double> radial_density_matrix;

        if (sho_tools::order_zyx == order) {
          
            energy_dm = 0; // double counting correction (with old matrix elements, is this a problem?)
            for (int izyx = 0; izyx < nSHO; ++izyx) {
                energy_dm += dot_product(nSHO, rho_matrix[izyx], hamiltonian[izyx]);
            } // izyx

            // we need to transform from Cartesian to Radial first

            // now transform the rho_matrix[izyx][jzyx]
            //     into a radial_density_matrix[ilmn][jlmn]
            //     using the unitary transform from left and right
            radial_density_matrix = view2D<double>(nSHO, nSHO); // get memory
            transform_SHO(radial_density_matrix.data(), radial_density_matrix.stride(),
                  rho_matrix.data(), rho_matrix.stride(), true);
#ifdef DEVEL
            if (1) { // debugging
                view2D<double> check_matrix(nSHO, nSHO);
                transform_SHO(check_matrix.data(), check_matrix.stride(),
                              radial_density_matrix.data(), radial_density_matrix.stride(), false);
                double maxdev{0};
                for (int i = 0; i < nSHO; ++i) {
                    for (int j = 0; j < nSHO; ++j) {
                        maxdev = std::max(maxdev, std::abs(check_matrix(i,j) - rho_matrix(i,j)));
                    } // j
                } // i
                if (maxdev > 1e-12) warn("%s found max deviation %.1e when backtransforming the density matrix", label, maxdev);
                if (maxdev > 1e-9) error("%s found max deviation %.1e when backtransforming the density matrix", label, maxdev);
            } // debugging
#endif // DEVEL
        } else if (sho_tools::order_lmn == order) {
            radial_density_matrix = view2D<double>(rho_matrix.data(), rho_matrix.stride()); // wrap
        } else {
            error("%s should be in either order_zyx or order_lmn", label);
        } // order

#ifdef DEVEL
        if (echo > 7 && take_spherical_density[valence] < 1) {
            std::printf("# %s Radial SHO density matrix in %s-order:\n", label, SHO_order2string(sho_tools::order_lmn).c_str());
            view2D<char> labels(sho_tools::nSHO(numax), 8, '\0');
            sho_tools::construct_label_table(labels.data(), numax, sho_tools::order_lmn);
            for (int ilmn = 0; ilmn < nSHO; ++ilmn) {
                std::printf("# %s %-8s ", label, labels[ilmn]);
                printf_vector(" %11.6f", radial_density_matrix[ilmn], nSHO);
            } // ilmn
            std::printf("\n");
        } // echo

        if (echo > 12) {
            int const nlnr = sho_tools::nSHO_radial(numax);
            view2D<double> density_matrix_ln(nlnr,nlnr);
            scattering_test::emm_average(density_matrix_ln.data(), radial_density_matrix.data(), int(numax), 0);
            show_ell_block_diagonal(density_matrix_ln, "emm-summed SHO density matrix", 1, true, true);
        } // echo
#endif // DEVEL

        if (synthetic_density_matrix) {
            create_synthetic_density_matrix(radial_density_matrix, echo);
        } else {
            // bring the iln,jln indices of rho_tensor into the partial-wave space
            auto const u_proj = unfold_projector_coefficients();
            int const N = nSHO;
            view2D<double> tmp_mat(N, N);
            auto const uT_proj = transpose(u_proj, N);
            gemm(tmp_mat, N, radial_density_matrix, N, uT_proj); // *u^T
            gemm(radial_density_matrix, N, u_proj, N, tmp_mat); // u*
        } // scope

#ifdef DEVEL
        auto const n_modify = 1 + int(control::get("single_atom.synthetic.density.matrix", 0.));
        for (int i_modify = 0; i_modify < n_modify; ++i_modify) {

            if (i_modify > 0) create_synthetic_density_matrix(radial_density_matrix, echo);

            // display radial density matrix in partial waves
            if (echo > 7) std::printf("# %s Radial density matrix in partial waves:\n", label);
            for (int ilmn = 0; ilmn < nSHO; ++ilmn) {
                int const iln = ln_index_list[ilmn];
                if (partial_wave_active[iln]) {
                    if (echo > 7) std::printf("# %s %-8s ", label, partial_wave[iln].tag);
                    for (int jlmn = 0; jlmn < nSHO; ++jlmn) {
                        int const jln = ln_index_list[jlmn];
                        if (partial_wave_active[jln]) {
                            if (echo > 7) std::printf(" %11.6f", radial_density_matrix(ilmn,jlmn));
                        } // active_j
                    } // jlmn
                    if (echo > 7) std::printf("\n");
                } // active_i
            } // ilmn
            if (echo > 7) std::printf("\n");
            
        } // i_modify

        if (echo > 2) {
            int const nlnr = sho_tools::nSHO_radial(numax);
            view2D<double> density_matrix_ln(nlnr,nlnr);
            scattering_test::emm_average(density_matrix_ln.data(), radial_density_matrix.data(), int(numax), 0);
            show_ell_block_diagonal(density_matrix_ln, "emm-summed density matrix"); // in partial waves
        } // echo
#endif // DEVEL

        // compute the kinetic energy contribution
        double E_kin[TRU_AND_SMT] = {0, 0};
        for (int ilmn = 0; ilmn < nSHO; ++ilmn) {
            int const iln = ln_index_list[ilmn];
            if (partial_wave_active[iln]) {
                for (int jlmn = 0; jlmn < nSHO; ++jlmn) {
                    int const jln = ln_index_list[jlmn];
                    if (partial_wave_active[jln]) {
                        E_kin[TRU] += radial_density_matrix(ilmn,jlmn) * kinetic_energy(TRU,iln,jln);
                        E_kin[SMT] += radial_density_matrix(ilmn,jlmn) * kinetic_energy(SMT,iln,jln);
                    } // active_j
                } // jlmn
            } // active_i
        } // ilmn
        set(energy_kin_csvn[3], TRU_AND_SMT, E_kin); // non-spherical kinetic energy contribution
        // energy_kin_csvn[3] and energy_dm have the problem that they use the mixed density matrix
        // while they should use the new density matrix to compute the term integral V_ij*D_ij
        

        //   Now, contract with the Gaunt tensor over m_1 and m_2
        initialize_Gaunt(); // make sure the Gaunt tensor has been precomputed
        //   rho_tensor[lm][iln][jln] =
        //     G_{lm l_1m_1 l_2m_2} * radial_density_matrix{il_1m_1n_1 jl_2m_2n_2}

        set(rho_tensor, nlm, 0.0); // clear
        for (auto gnt : gaunt) {
            int const lm = gnt.lm, ilm = gnt.lm1, jlm = gnt.lm2; auto G = gnt.G;
            // if (0 == lm) assert(std::abs(G - Y00*(ilm == jlm)) < 1e-12); // make sure that G_00ij = delta_ij*Y00
            if (0 == lm) G = Y00*(ilm == jlm); // make sure that G_00ij = delta_ij*Y00
            if (lm < nlm && ilm < mlm && jlm < mlm) {
                for (int ilmn = lmn_begin[ilm]; ilmn < lmn_end[ilm]; ++ilmn) {
                    int const iln = ln_index_list[ilmn];
                    for (int jlmn = lmn_begin[jlm]; jlmn < lmn_end[jlm]; ++jlmn) {
                        int const jln = ln_index_list[jlmn];
                        rho_tensor(lm,iln,jln) += G * radial_density_matrix(ilmn,jlmn);
#ifdef FULL_DEBUG
//                         auto const rho_ij = rho_tensor[lm][iln][jln];
//                         if (std::abs(rho_ij) > 1e-9)
//                             std::printf("# LINE=%d rho_ij = %g for lm=%d iln=%d jln=%d\n", __LINE__, rho_ij*Y004pi, lm, iln, jln);
#endif // FULL_DEBUG
                    } // jlmn
                } // ilmn
            } // limits
        } // gnt

#ifdef DEVEL
// #ifdef FULL_DEBUG
        if (echo > 0) {
            int const nln = sho_tools::nSHO_radial(numax);
            for (int lm = 0; lm < 1; ++lm) {
                for (int iln = 0; iln < nln; ++iln) {
                    for (int jln = 0; jln < nln; ++jln) {
                        auto const rho_ij = rho_tensor(lm,iln,jln);
                        if (std::abs(rho_ij) > 2e-16 && echo > 11) {
                            std::printf("# %s LINE=%d rho_ij = %g for lm=%d iln=%d jln=%d\n", 
                                label, __LINE__, rho_ij*Y004pi, lm, iln, jln);
                        } // rho_ij > 0
                    } // jln
                } // iln
            } // lm
        } // echo
// #endif // FULL_DEBUG
#endif // DEVEL
    } // get_rho_tensor


    void update_full_density(
          view3D<double> const & density_tensor // density tensor rho_{lm iln jln}
        , int const echo=0 // log-level
    ) {
        // from a density tensor and partial waves expand the non-spherical densities 

        int const nlm = pow2(1 + ellmax_rho);
        int const nln = sho_tools::nSHO_radial(numax);

        double const mix_valence_density = 1.0 - take_spherical_density[valence]; // fraction of valence density from partial waves
        for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            size_t const nr = rg[ts].n;
            assert(full_density[ts].stride() >= nr);
            for (int lm = 0; lm < nlm; ++lm) {
                set(full_density[ts][lm], full_density[ts].stride(), 0.0); // clear
                if (00 == lm) {
                    // add spherical densities, in particular the core density
                    for (int csv = core; csv <= valence; ++csv) {
                        if (take_spherical_density[csv] > 0 && csv_charge[csv] > 0) {
                            add_product(full_density[ts][00], nr, spherical_density[ts][csv], Y00*take_spherical_density[csv]); 
                            // needs scaling with Y00 since core_density has a factor 4*pi
                            if (echo > 2 + ts) std::printf("# %s %s density has %g electrons after adding the spherical %s density\n",
                                label, ts_name[ts], dot_product(nr, full_density[ts][00], rg[ts].r2dr)*Y004pi, csv_name(csv));
                        } // take
                    } // csv, spherical {core, semicore, valence} densities
                } // 00 == lm
                if (mix_valence_density > 0) {
                    for (int iln = 0; iln < nln; ++iln) {
                        if (partial_wave_active[iln]) {
                            auto const wave_i = partial_wave[iln].wave[ts];
                            for (int jln = 0; jln < nln; ++jln) {
                                if (partial_wave_active[jln]) {
                                    auto const wave_j = partial_wave[jln].wave[ts];
                                    double const rho_ij = density_tensor(lm,iln,jln) * mix_valence_density;
#ifdef FULL_DEBUG
                                    if (echo > 0 && 00 == lm && ts == TRU && std::abs(rho_ij) > 2e-16) 
                                        std::printf("# %s rho_ij = %g for lm=%d iln=%d jln=%d\n",
                                            label, rho_ij*Y004pi, lm, iln, jln);
#endif // FULL_DEBUG
                                    add_product(full_density[ts][lm], nr, wave_i, wave_j, rho_ij);
                                } // active
                            } // jln
                        } // active
                    } // iln
                } // mix_valence_density
            } // lm
            if (echo > 2) std::printf("# %s %s density has %g electrons\n",
                      label, ts_name[ts], dot_product(nr, full_density[ts][00], rg[ts].r2dr)*Y004pi);

        } // ts: true and smooth

        // determine the compensator charges
        int const nlm_cmp = pow2(1 + ellmax_cmp);
        for (int ell = 0; ell <= ellmax_cmp; ++ell) {
            for (int emm = -ell; emm <= ell; ++emm) {
                int const lm = solid_harmonics::lm_index(ell, emm);
                double rho_lm{0}, tru_lm{0}, smt_lm{0};
                if (mix_valence_density > 0 || echo > 0) {
                    for (int iln = 0; iln < nln; ++iln) {
                        for (int jln = 0; jln < nln; ++jln) {
                            double const rho_ij = density_tensor(lm,iln,jln);
#ifdef FULL_DEBUG
                            if (echo > 0 && std::abs(rho_ij) > 1e-9)
                                std::printf("# %s rho_ij = %g for ell=%d emm=%d iln=%d jln=%d\n",
                                    label, rho_ij*Y004pi, ell, emm, iln, jln);
#endif // FULL_DEBUG
                            rho_lm += rho_ij * ( charge_deficit(ell,TRU,iln,jln)
                                               - charge_deficit(ell,SMT,iln,jln) );
                            tru_lm += rho_ij *   charge_deficit(ell,TRU,iln,jln)  ; // only for display
                            smt_lm += rho_ij *   charge_deficit(ell,SMT,iln,jln)  ; // only for display
                        } // jln
                    } // iln
                } // mix_valence_density
                assert(lm >= 0);
                assert(lm < nlm_cmp);
                qlm_compensator[lm] = rho_lm * mix_valence_density;
                if (0 == ell && echo > 2) std::printf("# %s valence density matrix proposes %g true, %g smooth electrons\n", label, tru_lm*Y004pi, smt_lm*Y004pi);

            } // emm
        } // ell

        // account for Z protons in the nucleus and the missing charge in the spherical densities
        assert(1 == take_spherical_density[core]); // must always be 1 since we can compute the core density only on the radial grid
        double const spherical_charge_deficits = dot_product(3, spherical_charge_deficit, take_spherical_density);
        qlm_compensator[00] += (spherical_charge_deficits - Z_core)*Y00;
        if (echo > 5) std::printf("# %s compensator monopole charge is %g electrons\n", label, qlm_compensator[00]*Y004pi);
        
        
        { // scope: measure the ionization inside the sphere, should be small
            double sph_charge{0};
            for (int csv = core; csv <= valence; ++csv) {
                sph_charge += dot_product(ir_cut[TRU] + 1, rg[TRU].r2dr, spherical_density[TRU][csv]);
            } // csv
            if (echo > 2) std::printf("# %s true spherical density has %g electrons inside the sphere\n", label, sph_charge);
            double const s00_charge = dot_product(ir_cut[TRU] + 1, rg[TRU].r2dr, full_density[TRU][00])*Y004pi;
            if (echo > 2) std::printf("# %s true full density has %g electrons inside the sphere\n", label, s00_charge);
            if (echo > 2) std::printf("# %s density shows an ionization of %g electrons inside the sphere\n", label, s00_charge - sph_charge);
        } // scope
            

        { // scope: construct the augmented density
            int const nlm_aug = pow2(1 + std::max(ellmax_rho, ellmax_cmp));
            auto const mr = full_density[SMT].stride(); // on the smooth grid
            assert(aug_density.stride() == mr);
            set(aug_density, nlm_aug, 0.0); // clear all entries
            set(aug_density.data(), nlm*mr, full_density[SMT].data()); // copy smooth full_density, need spin summation?
            add_or_project_compensators<0>(aug_density, qlm_compensator.data(), rg[SMT], ellmax_cmp, sigma_compensator);
             
#ifdef DEVEL
            if (echo > 19) {
                std::printf("\n## r, aug_density, true_density, smooth_density, spherical_density:\n");
                for (int ir = 0; ir < rg[SMT].n; ++ir) {
                    int const ir_tru = ir + nr_diff;
                    std::printf("%g %g %g %g %g\n", rg[SMT].r[ir], aug_density(00,ir),
                          full_density[TRU](00,ir_tru), full_density[SMT](00,ir),
                          (spherical_density[TRU](core,ir_tru) + spherical_density[TRU](valence,ir_tru))*Y00);
                } // ir
                std::printf("\n\n");
            } // echo
#endif // DEVEL

            double const tru_charge = dot_product(rg[TRU].n, rg[TRU].r2dr, full_density[TRU][00]); // only full_density[0==lm]
            if (echo > 3) std::printf("# %s true density has %g electrons\n", label, tru_charge*Y004pi); // this value can differ ...
            // ... from Z_core since we integrate over the entire grid
        } // scope

    } // update_full_density










    status_t update_full_potential(
          float const mixing // how much of the true potential is transferred to the spherical potential
        , double const ves_multipole[]
        , int const echo=0
    ) {
        // from non-spherical densities create a non-spherical potential V_{ell emm}(r)
        status_t stat(0);

        // due to its non-linearity, evaluate the exchange-correlation contributions on an angular grid
        int const nlm = pow2(1 + ellmax_rho);
        assert(ellmax_rho == ellmax_pot);
        int const npt = angular_grid::get_grid_size(ellmax_rho);

        std::vector<double> vlm(nlm, 0.0); // potential shifts
        for (int ts = SMT; ts >= TRU; --ts) { // smooth quantities first, so we can determine vlm
            int const nr = rg[ts].n;
            auto const mr = full_density[ts].stride();
            assert(mr >= nr); // stride

            { // scope: quantities on the angular grid
                if (echo > 6) std::printf("# %s quantities on the angular grid are %i * %li = %li\n", label, npt, mr, npt*mr);
                view2D<double> on_grid(2, npt*mr);
                auto const rho_on_grid = on_grid[0];
                auto const vxc_on_grid = on_grid[0];
                auto const exc_on_grid = on_grid[1];

                // transform the lm-index into real-space
                // using an angular grid quadrature, e.g. Lebedev-Laikov grids
                if (echo > 6 && SMT == ts) std::printf("# %s local smooth density at origin %g a.u.\n",
                                                          label, full_density[ts](00,0)*Y00);
                stat += angular_grid::transform(rho_on_grid, full_density[ts].data(), mr, ellmax_rho, false);
                // envoke the exchange-correlation potential (acts in place)
//              if (echo > 7) std::printf("# envoke the exchange-correlation on angular grid\n");
                for (size_t ip = 0; ip < npt*mr; ++ip) {
                    double const rho = rho_on_grid[ip];
                    double vxc{0};
                    exc_on_grid[ip] = exchange_correlation::lda_PZ81_kernel(rho, vxc); // only needed if we want to compute the total energy
                    vxc_on_grid[ip] = vxc; // overwrites rho_on_grid[ip] due to pointer aliasing
                } // ip
                // transform back to lm-index
                assert(full_potential[ts].stride() == mr);
                stat += angular_grid::transform(full_potential[ts].data(), vxc_on_grid, mr, ellmax_pot, true);
                { // scope: transform also the exchange-correlation energy density
                    view2D<double> exc_lm(nlm, mr);
                    stat += angular_grid::transform(exc_lm.data(), exc_on_grid, mr, ellmax_rho, true);
                    if (SMT == ts) {
                        if (echo > 7) std::printf("# %s local smooth exchange-correlation potential at origin is %g %s\n",
                                                            label, full_potential[SMT](00,0)*Y00*eV,_eV);
                    } // SMT only
                    
                    double E_dc{0}, E_xc{0};

                    auto const Edc00 = dot_product(nr, full_potential[ts][00], full_density[ts][00], rg[ts].r2dr); // dot_product with diagonal metric
                    E_dc += Edc00;
                    if (echo > 5) std::printf("# %s double counting correction  in %s 00 channel %.9f %s\n",
                                            label, ts_name[ts], Edc00*eV,_eV);
                    auto const Exc00 = dot_product(nr, exc_lm[00], full_density[ts][00], rg[ts].r2dr); // dot_product with diagonal metric
                    E_xc += Exc00;
                    if (echo > 5) std::printf("# %s exchange-correlation energy in %s 00 channel %.9f %s\n",
                                            label, ts_name[ts], Exc00*eV,_eV);
                    for (int ell = 1; ell <= std::min(ellmax_rho, ellmax_pot); ++ell) {
                        double Edc_L{0}, Exc_L{0};
                        for (int emm = -ell; emm <= ell; ++emm) {
                            int const ilm = sho_tools::lm_index(ell, emm);
                            Edc_L += dot_product(nr, full_potential[ts][ilm], full_density[ts][ilm], rg[ts].r2dr); // dot_product with diagonal metric
                            Exc_L += dot_product(nr, exc_lm[ilm], full_density[ts][ilm], rg[ts].r2dr); // dot_product with diagonal metric
                        } // emm
                        if (echo > 5 + ell) std::printf("# %s double counting correction  in %s ell=%i channel %g %s\n",
                                        label, ts_name[ts], ell, Edc_L*eV,_eV);
                        if (echo > 5 + ell) std::printf("# %s exchange-correlation energy in %s ell=%i channel %g %s\n",
                                        label, ts_name[ts], ell, Exc_L*eV,_eV);
                        E_dc += Edc_L;
                        E_xc += Exc_L;
                    } // ell

                    energy_dc[ts] = E_dc;
                    energy_xc[ts] = E_xc;

                } // scope
            } // scope: quantities on the angular grid

            // solve electrostatics inside the spheres
            view2D<double> Ves(nlm, mr);
            double const q_nucleus = (TRU == ts) ? -Z_core*Y00 : 0; // Z_core = number of protons in the nucleus
            auto   const & rho_aug = (TRU == ts) ? full_density[TRU] : aug_density;
            // solve electrostatics with singularity (q_nucleus) // but no outer boundary conditions (v_lm)
            radial_potential::Hartree_potential(Ves.data(), rg[ts], rho_aug.data(), rho_aug.stride(), ellmax_pot, q_nucleus);

            if (SMT == ts) {
                add_or_project_compensators<1>(Ves, vlm.data(), rg[SMT], ellmax_cmp, sigma_compensator); // project Ves to compensators
                if (echo > 7) std::printf("# %s inner integral between normalized compensator and smooth Ves(r) = %g %s\n", label, vlm[0]*Y00*eV,_eV);

                if (echo > 21) {
                    std::printf("\n## %s smooth spherical electrostatic potential (a.u.):\n", label);
                    print_compressed(rg[ts].r, Ves[00], rg[ts].n);                    
                } // echo

                if (echo > 21) {
                    std::printf("\n## %s smooth spherical exchange-correlation potential (a.u.):\n", label);
                    print_compressed(rg[ts].r, full_potential[ts][00], rg[ts].n);                    
                } // echo

                // but the solution of the 3D problem found that these integrals should be ves_multipole, therefore
                if (nullptr == ves_multipole) {
                    set(vlm.data(), nlm, 0.); // no correction of the electrostatic potential heights for isolated atoms
                } else {
                    if (echo > 6) std::printf("# %s v_00 found %g but expected %g, shift by %g %s\n", // report monopole shift
                        label, vlm[00]*Y00*eV, ves_multipole[00]*Y00*eV, ves_multipole[00]*Y00*eV - vlm[00]*Y00*eV, _eV);
                    scale(vlm.data(), nlm, -1.); add_product(vlm.data(), nlm, ves_multipole, 1.); // vlm := ves_multipole - vlm
                } // no ves_multipole given
            } // smooth only

            if (SMT == ts) {
                if (echo > 7) std::printf("# %s local smooth electrostatic potential at origin is %g %s\n", label, Ves(00,0)*Y00*eV,_eV);
            } // SMT only

            add_or_project_compensators<2>(Ves, vlm.data(), rg[ts], ellmax_cmp, sigma_compensator);

            
            double E_Hartree{0};
            for (int ilm = 0; ilm < pow2(1 + ellmax_pot); ++ilm) {
                E_Hartree += 0.5*dot_product(nr, Ves[ilm], rho_aug[ilm], rg[ts].r2dr);
            } // ilm
            energy_es[ts] = E_Hartree;
            
            if (SMT == ts) { // debug: project again to see if the correction worked out for the ell=0 channel
                std::vector<double> v_test(pow2(1 + ellmax_cmp), 0.0);
                add_or_project_compensators<1>(Ves, v_test.data(), rg[SMT], ellmax_cmp, sigma_compensator); // project to compensators with ellmax_cmp=0
                if (echo > 7) {
                    std::printf("# %s after correction v_00 is %g %s\n", label, v_test[00]*Y00*eV,_eV);
                    if (1) {
                        std::printf("# %s after correction v_lm (%s Bohr^-ell) is", label, _eV);
                        printf_vector(" %.6f", v_test.data(), v_test.size(), "\n", Y00*eV); // Warning, not consistent with _Ang != "Bohr"
                    } // 1
                    std::printf("# %s local smooth electrostatic potential at origin is %g %s\n", label, Ves(00,0)*Y00*eV,_eV);
                    std::printf("# %s local smooth augmented density at origin is %g a.u.\n", label, aug_density(00,0)*Y00);
                    if (echo > 8) {
                        std::printf("\n## %s local smooth electrostatic potential and augmented density in a.u.:\n", label);
                        for (int ir = 0; ir < rg[SMT].n; ++ir) {
                            std::printf("%g %g %g\n", rg[SMT].r[ir], Ves(00,ir)*Y00, aug_density(00,ir)*Y00);
                        } // ir
                        std::printf("\n\n");
                    } // show radial function of Ves[00]*Y00 to be compared to projections of the 3D electrostatic potential
                } // echo
                if (echo > 5) std::printf("# %s smooth Hartree energy %.9f %s\n", label, energy_es[SMT]*eV, _eV);
            } else {
                // TRU
                if (echo > 8) std::printf("# %s local true electrostatic potential*r at origin is %g (should match -Z=%.1f)\n",
                                        label, Ves(00,1)*(rg[TRU].r[1])*Y00, -Z_core);
                auto const E_Coulomb = -Z_core*dot_product(nr, full_density[TRU][00], rg[TRU].rdr)*Y004pi;
                if (echo > 5) std::printf("# %s true Hartree energy %.9f Coulomb energy %.9f %s\n",
                                        label, (energy_es[TRU] - 0.5*E_Coulomb)*eV, E_Coulomb*eV, _eV);
                energy_es[TRU] += 0.5*E_Coulomb; // the other half is included in E_Hartree
            } // ts

            add_product(full_potential[ts].data(), nlm*mr, Ves.data(), 1.0); // add the electrostatic potential, scale_factor=1.0

        } // ts smooth and true
        if (0 != stat) warn("%s angular grid transformations failed with status sum = %i", label, int(stat));

        
        auto const df = Y00*eV; assert(df > 0); // display factor
        if (local_potential_method >= 0) {
            // construct the zero_potential V_bar

            if (echo > 5) std::printf("# %s old zero potential: V_bar(0) = %g, V_bar(R_cut) = %g, V_bar(R_max) = %g %s\n",
                label, zero_potential[0]*df, zero_potential[ir_cut[SMT]]*df, zero_potential[rg[SMT].n - 1]*df, _eV);

            set(zero_potential.data(), zero_potential.size(), 0.0); // init zero
            
            std::vector<double> V_smt(rg[SMT].n);
            auto const stat_pseudo = pseudo_tools::pseudize_local_potential<0>(V_smt.data(),
                                    full_potential[TRU][00], rg, ir_cut, local_potential_method, label, echo, Y00);
            if (stat_pseudo) {
                warn("%s matching procedure for the local potential failed, status= %i", label, int(stat_pseudo));
                stat += stat_pseudo;
            } else {
                if (echo > 5) std::printf("# %s smooth potential: V_smt(0) = %g, V_smt(R_cut) = %g %s\n",
                                        label, V_smt[0]*df, V_smt[ir_cut[SMT]]*df, _eV);
                for (int ir = 0; ir < rg[SMT].n; ++ir) {
                    zero_potential[ir] = V_smt[ir] - full_potential[SMT](00,ir);
                } // ir
#ifdef DEVEL
                if (echo > 21) {
                    std::printf("\n## %s pseudized total potential (a.u.):\n", label);
                    print_compressed(rg[SMT].r, V_smt.data(), rg[SMT].n);                    
                } // echo
#endif // DEVEL
            } // pseudization successful
        } else {
            if (echo > 6) std::printf("# %s zero_potential stays unchanged!\n", label);
        } // modify zero_potential

#ifdef DEVEL
        if (echo > 31) {
            std::printf("# %s local smooth zero_potential:\n", label);
            for (int ir = 0; ir < rg[SMT].n; ++ir) {
                std::printf("%g %g\n", rg[SMT].r[ir], zero_potential[ir]*Y00);
            } // ir
            std::printf("\n\n");
        } // echo

        { // scope: analyze zero potential
            double vol = 0, Vint = 0, r1Vint = 0, r2Vint = 0;
            for (int ir = ir_cut[SMT]; ir < rg[SMT].n; ++ir) {
                auto const r  = rg[SMT].r[ir];
                auto const dV = rg[SMT].r2dr[ir];
                vol    +=                    dV;
                Vint   += zero_potential[ir]*dV;
                r1Vint += zero_potential[ir]*dV*r;
                r2Vint += zero_potential[ir]*dV*r*r;
            } // ir
            if (echo > 5) std::printf("# %s zero_potential statistics = %g %g %g %s\n",
                        label, Vint/vol*eV, r1Vint/(vol*r_cut)*eV, r2Vint/(vol*pow2(r_cut))*eV, _eV);
                // these numbers should be small since they indicate that V_bar is localized inside the sphere
                // and how much V_smt deviates from V_tru outside the sphere
        } // scope: analyze zero potential

        if (echo > 21) {
            std::printf("\n## %s zero potential (a.u.):\n", label);
            print_compressed(rg[SMT].r, zero_potential.data(), rg[SMT].n);                    
        } // echo
#endif // DEVEL
        if (echo > 5) std::printf("# %s zero_potential: V_bar(0) = %g, V_bar(R_cut) = %g, V_bar(R_max) = %g %s\n",
            label, zero_potential[0]*df, zero_potential[ir_cut[SMT]]*df, zero_potential[rg[SMT].n - 1]*df, _eV);

        // add spherical zero potential for SMT==ts and 00==lm
        add_product(full_potential[SMT][00], rg[SMT].n, zero_potential.data(), 1.0);
#ifdef DEVEL
        if (echo > 21) {
            std::printf("\n## %s smooth total potential (a.u.):\n", label);
            print_compressed(rg[SMT].r, full_potential[SMT][00], rg[SMT].n);                    
        } // echo
#endif // DEVEL

        // feed back spherical part of the full potentials into the spherical potentials r*V
        // which determines core states and true partial waves
        if (mixing > 0) {
            for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                scale(potential[ts].data(), rg[ts].n, 1. - mixing);
                add_product(potential[ts].data(), rg[ts].n, full_potential[ts][00], rg[ts].r, Y00*mixing);
            } // ts true and smooth
        } // mixing > 0

        // fix the true spherical potential at the origin r=0
        potential[TRU][0] = -Z_core; // r*V(r) is finite due to V(r) diverging to -infinity as -Z/r

#ifdef DEVEL
        if (0) { // scope: test: use the spherical routines from atom_core::rad_pot(output=r*V(r), input=rho(r)*4pi)
            std::vector<double> rho4pi(rg[TRU].n);
            set(rho4pi.data(), rg[TRU].n, full_density[TRU][00], Y004pi);
            std::printf("\n# WARNING: use rad_pot to construct the r*V_tru(r) [for DEBUGGING]\n\n");
            atom_core::rad_pot(potential[TRU].data(), rg[TRU], rho4pi.data(), Z_core);
        } // scope
#endif // DEVEL
        return stat;
    } // update_full_potential









    void update_matrix_elements(
          int const echo=0 // log-level
    ) {
        // create matrix elements for the PAW correction of the Hamiltonian and overlap 
        // from the current partial waves and the non-sperical potential

        int const nlm = pow2(1 + ellmax_pot);
        int const nln = sho_tools::nSHO_radial(numax);
        int const nSHO = sho_tools::nSHO(numax);
        int const nlmn = nSHO;

        // first generate the matrix elements in the valence basis
        //    overlap[iln][jln] and potential_lm[iln][jln]
        // then, generate the matrix elements in the radial representation
        //    overlap[ilmn][jlmn] and hamiltonian[ilmn][jlmn]
        //    where hamiltonian[ilmn][jlmn] = kinetic_energy_deficit_{iln jln}
        //            + sum_lm Gaunt_{lm ilm jlm} * potential_{lm iln jln}
        // then, transform the matrix elements into the Cartesian representation using sho_unitary
        //    overlap[iSHO][jSHO] and hamiltonian[iSHO][jSHO]
#ifdef DEVEL
        if (echo > 19) {
            std::printf("\n\n## %s compare potentials (Bohr, Ha, Ha, Ha, Ha) r, r*V_tru[00](r), r*V_smt[00](r), r*V_tru(r), r*V_smt(r):\n", label);
            for (int ir = 1; ir < rg[SMT].n; ++ir) {
                double const r = rg[SMT].r[ir];
                std::printf("%g %g %g %g %g\n", r, r*full_potential[TRU](00,ir + nr_diff)*Y00,
                    r*full_potential[SMT](00,ir)*Y00, potential[TRU][ir + nr_diff], potential[SMT][ir]);
            } // ir
            std::printf("\n\n");
        } // echo
#endif // DEVEL

        // construct the potential tensor in terms of the emm_Degenerate partial waves
        view4D<double> potential_ln(nlm, TRU_AND_SMT, nln, nln, 0.0); // get memory, // emm1-emm2-degenerate
        for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            int const nr = rg[ts].n;
            std::vector<double> wave_pot_r2dr(nr);
            for (int ell = 0; ell <= ellmax_pot; ++ell) {
                for (int emm = -ell; emm <= ell; ++emm) {
                    int const lm = solid_harmonics::lm_index(ell, emm);
                    assert(lm < nlm);
                    // similar but not the same with check_spherical_matrix_elements
                    for (int iln = 0; iln < nln; ++iln) {
                        if (partial_wave_active[iln]) {
                            auto const wave_i = partial_wave[iln].wave[ts];
                            product(wave_pot_r2dr.data(), nr, wave_i, full_potential[ts][lm], rg[ts].r2dr);
                            for (int jln = 0; jln < nln; ++jln) {
                                if (partial_wave_active[jln]) {
                                    auto const wave_j = partial_wave[jln].wave[ts];
                                    potential_ln(lm,ts,iln,jln) = dot_product(nr, wave_pot_r2dr.data(), wave_j);
                                } // active
                            } // jln
                        } // active
                    } // iln
                } // emm
            } // ell
        } // ts true and smooth

        view3D<double> matrices_lmn(2, nlmn, nlmn, 0.0); // get memory
        auto hamiltonian_lmn = matrices_lmn[0],
                 overlap_lmn = matrices_lmn[1];

        initialize_Gaunt(); // make sure the Gaunt tensor is precomputed

        // start by adding the contribution of the non-spherical potential
        int const mlm = pow2(1 + numax);
        for (auto gnt : gaunt) {
            int const lm = gnt.lm, lm1 = gnt.lm1, lm2 = gnt.lm2; auto G = gnt.G;
            if (00 == lm) G = Y00*(lm1 == lm2); // make sure that G_00ij = delta_ij*Y00
            if (lm1 < mlm && lm2 < mlm) {
                if (lm < nlm) {
                    for (int ilmn = lmn_begin[lm1]; ilmn < lmn_end[lm1]; ++ilmn) {
                        int const iln = ln_index_list[ilmn];
                        for (int jlmn = lmn_begin[lm2]; jlmn < lmn_end[lm2]; ++jlmn) {
                            int const jln = ln_index_list[jlmn];
                            hamiltonian_lmn(ilmn,jlmn) +=
                              G * ( potential_ln(lm,TRU,iln,jln)
                                  - potential_ln(lm,SMT,iln,jln) );
                        } // jlmn
                    } // ilmn
                } // lm
            } // limits
        } // gnt

        // add the kinetic_energy deficit to the hamiltonian
        if (echo > 7) std::printf("\n# %s Hamiltonian elements %s-ordered in %s:\n",
                        label, sho_tools::SHO_order2string(sho_tools::order_lmn).c_str(), _eV);
        for (int ilmn = 0; ilmn < nlmn; ++ilmn) {
            int const iln = ln_index_list[ilmn];
            int const ilm = lm_index_list[ilmn];
            for (int jlmn = 0; jlmn < nlmn; ++jlmn) {
                int const jln = ln_index_list[jlmn];
                int const jlm = lm_index_list[jlmn];
                if (ilm == jlm) {
                    hamiltonian_lmn(ilmn,jlmn) += ( kinetic_energy(TRU,iln,jln) 
                                                  - kinetic_energy(SMT,iln,jln) );
                    overlap_lmn(ilmn,jlmn) = ( charge_deficit(0,TRU,iln,jln)
                                             - charge_deficit(0,SMT,iln,jln) ); // ell=0
                } // diagonal in lm, offdiagonal in nrn
            } // jlmn
            if (echo > 7) {
                std::printf("# %s hamiltonian elements for ilmn=%3i  ", label, ilmn);
                printf_vector(" %7.3f", hamiltonian_lmn[ilmn], nlmn, "\n", eV);
            } // echo
        } // ilmn

        auto const u_proj = unfold_projector_coefficients();
#ifdef DEVEL
        if (echo > 8) { // display
            std::printf("\n# %s lmn-based projector matrix:\n", label);
            int const mlmn = display_delimiter(numax, nn, 'm');
            for (int ilmn = 0; ilmn < nlmn; ++ilmn) {
                if (partial_wave_active[ln_index_list[ilmn]]) {
                    std::printf("# %s u_proj for ilmn=%3i  ", label, ilmn);
                    printf_vector(" %g", u_proj[ilmn], mlmn);
                } // active
            } // ilmn
        } // echo
#endif // DEVEL

        { // scope: wrap the matrices with the projector_coeff from left and right
            auto const uT_proj = transpose(u_proj, nlmn);
            view2D<double> tmp_mat(nlmn, nlmn, 0.0);
            for (int iHS = 0; iHS < 2; ++iHS) {
                auto matrix_lmn = matrices_lmn[iHS];
                gemm(tmp_mat, nlmn, matrix_lmn, nlmn, u_proj); // *u
                gemm(matrix_lmn, nlmn, uT_proj, nlmn, tmp_mat); // u^T*
            } // iHS
        } // scope

#ifdef DEVEL
        if (echo > 8) { // display
            std::printf("\n# %s lmn-based Overlap elements:\n", label);
            for (int ilmn = 0; ilmn < nlmn; ++ilmn) {
                std::printf("# %s overlap elements for ilmn=%3i  ", label, ilmn);
                printf_vector(" %g", overlap_lmn[ilmn], nlmn);
            } // ilmn
            std::printf("\n");
        } // echo

        if (1) { // scope: check if averaging over emm gives back the same operators in the case of a spherical potential

            view3D<double> matrices_ln(2, nln, nln); // get memory
            for (int iHS = 0; iHS < 2; ++iHS) {
                scattering_test::emm_average(matrices_ln(iHS,0), matrices_lmn(iHS,0), int(numax));
            } // iHS
            if (echo > 4) {
                std::printf("\n");
                show_ell_block_diagonal(matrices_ln[0], "emm-averaged hamiltonian", eV, true, true);
                show_ell_block_diagonal(matrices_ln[1], "emm-averaged overlap    ",  1, true, true);
            } // echo

#ifdef    NEVER
            if (1) {
                // Warning: can only produce the same eigenenergies if potentials are converged:
                //          i.e.      Y00*r*full_potential[ts][00](r) == potential[ts](r)
                if (echo > 1) std::printf("\n\n# %s perform a diagonalization of the pseudo Hamiltonian\n\n", label);
                // prepare a smooth local potential which goes to zero at Rmax
                std::vector<double> Vsmt(rg[SMT].n, 0);
                double const V_rmax = full_potential[SMT](00,rg[SMT].n - 1)*Y00;
                for (int ir = 0; ir < rg[SMT].n; ++ir) {
                    Vsmt[ir] = full_potential[SMT](00,ir)*Y00 - V_rmax;
                } // ir
                if (echo > 1) std::printf("\n# %s %s eigenstate_analysis\n\n", label, __func__);
                if (echo > 5) std::printf("# local potential for eigenstate_analysis is shifted by %g %s\n", V_rmax*eV,_eV);
                scattering_test::eigenstate_analysis // find the eigenstates of the spherical Hamiltonian
                  (rg[SMT], Vsmt.data(), sigma, int(numax + 1), numax, matrices_ln(0,0), matrices_ln(1,0), 384, V_rmax, label, echo);
                std::fflush(stdout);
            } else if (echo > 0) std::printf("\n# eigenstate_analysis deactivated for now! %s %s:%i\n\n", __func__, __FILE__, __LINE__);

            if (1) {
                // Warning: can only produce the same eigenenergies if potentials are converged:
                //          i.e.      Y00*r*full_potential[ts][00](r) == potential[ts](r)
                std::vector<double> rV_vec[TRU_AND_SMT];
                for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                    rV_vec[ts] = std::vector<double>(rg[ts].n);
                    product(rV_vec[ts].data(), rg[ts].n, full_potential[ts][00], rg[ts].r, Y00); // rV(r) = V_00(r)*r*Y00
                } // ts
                double const *const rV[TRU_AND_SMT] = {rV_vec[TRU].data(), rV_vec[SMT].data()};
                if (echo > 1) std::printf("\n# %s %s logarithmic_derivative\n\n", label, __func__);
                scattering_test::logarithmic_derivative // scan the logarithmic derivatives
                  (rg, rV, sigma, int(numax + 1), numax, matrices_ln(0,0), matrices_ln(1,0), logder_energy_range, label, echo);
                std::fflush(stdout);
            } else if (echo > 0) std::printf("\n# logarithmic_derivative deactivated for now! %s %s:%i\n\n", __func__, __FILE__, __LINE__);
#endif // NEVER

        } // scope
#endif // DEVEL

        set(hamiltonian, nSHO, 0.0); // clear
        set(overlap,     nSHO, 0.0); // clear
        // Now transform _lmn quantities to Cartesian representations using sho_unitary
        transform_SHO(hamiltonian.data(), hamiltonian.stride(), hamiltonian_lmn.data(), hamiltonian_lmn.stride(), false);
        transform_SHO(    overlap.data(),     overlap.stride(),     overlap_lmn.data(),     overlap_lmn.stride(), false);
        // Mind that this transform is unitary and assumes square-normalized SHO-projectors
        // ... which might require proper normalization factors f(i)*f(j) to be multiplied in, see sho_projection::sho_prefactor
#ifdef DEVEL
        if (echo > 7) { // display
            std::printf("\n# %s SHO-transformed Hamiltonian elements (%s-order) in %s:\n",
                        label, sho_tools::SHO_order2string(sho_tools::order_zyx).c_str(), _eV);
            view2D<char> zyx_label(nSHO, 8);
            sho_tools::construct_label_table<8>(zyx_label.data(), numax, sho_tools::order_zyx);
            for (int iSHO = 0; iSHO < nSHO; ++iSHO) {
                std::printf("# %s hamiltonian elements for %-6s", label, zyx_label[iSHO]);
                printf_vector(" %11.6f", hamiltonian[iSHO], nSHO, "\n", eV);
            } // iSHO
        } // echo
#endif // DEVEL

    } // update_matrix_elements










    void check_spherical_matrix_elements(
          int const echo // log-level
    ) const {
        // construct the matrix elements for the spherical potentials
        // and compute logarithmic derivatives
        // and eigenstates on a radial equidistant grid
        // no member variable is affected

        // check if the emm-averaged Hamiltonian elements produce the same scattering properties
        // as the spherical part of the full potential
        int const nln = sho_tools::nSHO_radial(numax);
        std::vector<int8_t> ells(nln, -1);
        for (int ell = 0; ell <= numax; ++ell) {
            for (int nrn = 0; nrn < nn[ell]; ++nrn) {
                int const iln = sho_tools::ln_index(numax, ell, nrn);
                ells[iln] = ell;
            } // nrn
        } // ell

        view3D<double> potential_ln(TRU_AND_SMT, nln, nln); // fully emm-degenerate potential matrix elements
        for (int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            int const nr = rg[ts].n;
            std::vector<double> wave_pot_r2dr(nr); // temporary product of partial wave, potential and metric
            for (int ell = 0; ell <= numax; ++ell) {
                for (int nrn = 0; nrn < nn[ell]; ++nrn) {
                    int const iln = sho_tools::ln_index(numax, ell, nrn);
                    if (partial_wave_active[iln]) {
                        auto const wave_i = partial_wave[iln].wave[ts];
                        // potential is defined as r*V(r), so we only need r*dr to get r^2*dr as integration weights
                        product(wave_pot_r2dr.data(), nr, wave_i, potential[ts].data(), rg[ts].rdr);
                        for (int mrn = 0; mrn < nn[ell]; ++mrn) {
                            int const jln = sho_tools::ln_index(numax, ell, mrn);
                            if (partial_wave_active[jln]) {
                                auto const wave_j = partial_wave[jln].wave[ts];
                                potential_ln(ts,iln,jln) = dot_product(nr, wave_pot_r2dr.data(), wave_j);
                            } // active
                        } // mrn
                    } // active
                } // nrn
            } // ell
        } // ts: true and smooth

        view3D<double> aHSm(2, nln, nln, 0.0);
        auto hamiltonian_ln = aHSm[0], overlap_ln = aHSm[1];
        { // scope
            for (int iln = 0; iln < nln; ++iln) {
                for (int jln = 0; jln < nln; ++jln) {
                    if (ells[iln] == ells[jln]) {
                        hamiltonian_ln(iln,jln) = ( kinetic_energy(TRU,iln,jln)
                                                  - kinetic_energy(SMT,iln,jln) )
                                                + ( potential_ln(TRU,iln,jln)
                                                  - potential_ln(SMT,iln,jln) );
                        int const ell = 0;
                        overlap_ln(iln,jln) = ( charge_deficit(ell,TRU,iln,jln)
                                              - charge_deficit(ell,SMT,iln,jln) );
                    } // ells match
                } // jln
            } // iln
        } // scope


#ifdef DEVEL
        if (echo > 4) { // display
            show_ell_block_diagonal(aHSm[0], "spherical hamiltonian", eV);
            show_ell_block_diagonal(aHSm[1], "spherical overlap    ");
        } // echo
#endif // DEVEL

        auto const u_proj = unfold_projector_coefficients(emm_Degenerate);

#ifdef DEVEL
        if (echo > 2) { // display
            int const mln = display_delimiter(numax, nn);
            std::printf("\n# %s ln-based projector matrix:\n", label);
            if (1) { // show a legend
                view2D<char> ln_labels(nln, 8);
                sho_tools::construct_label_table<8>(ln_labels.data(), numax, sho_tools::order_ln);
                std::printf("# %s              radial SHO:      ", label);
                for (int jln = 0; jln < mln; ++jln) {
                    std::printf("%-12s", ln_labels[jln]);
                } // jln
                std::printf("\n");
            } // show a legend
            show_ell_block_diagonal(u_proj, "ln-based projector", 1, false, true);
        } // echo
#endif // DEVEL

        { // scope: wrap the matrices with the projector_coeff from left and right, in-place
            auto const uT_proj = transpose(u_proj, nln);
            view2D<double> tmp_mat(nln, nln, 0.0); // get temporary memory
            for (int iHS = 0; iHS < 2; ++iHS) {
                auto matrix_ln = aHSm[iHS]; // get a non-const sub-view
                gemm(tmp_mat, nln, matrix_ln, nln, u_proj); // *u
                gemm(matrix_ln, nln, uT_proj, nln, tmp_mat); // u^T*
#ifdef DEVEL
                if (echo > 4) { // display
                    std::printf("\n# %s spherical %s in radial SHO basis:\n", label, iHS?"overlap":"hamiltonian");
                    show_ell_block_diagonal(matrix_ln, iHS?"overlap    ":"hamiltonian" , iHS?1:eV, true, true);
                } // echo
#endif // DEVEL
            } // iHS
        } // scope

#ifdef DEVEL
        int const check_overlap_eigenvalues = control::get("single_atom.check.overlap.eigenvalues", 1.); // 0:false, 1:true, -10:true and abort
#else  // DEVEL
        int const check_overlap_eigenvalues = 1; // 1:true
#endif // DEVEL
        if (check_overlap_eigenvalues) {
            auto const overlap_ln = aHSm[1];
            for (int ell = 0; ell <= numax; ++ell) {
                if (nn[ell] > 0) {
                 
                    int const nmx = sho_tools::nn_max(numax, ell);
                    assert(nmx >= 1);
                    std::vector<double> eigvals(nmx + 2, 0.0);
                    // create a copy since the eigenvalue routine modifies the memory regions
                    view2D<double> ovl_ell(nmx, nmx); 
                    for (int irn = 0; irn < nmx; ++irn) {
                        int const iln = sho_tools::ln_index(numax, ell, irn);
                        int const jln = sho_tools::ln_index(numax, ell, 0);
                        set(ovl_ell[irn], nmx, &overlap_ln(iln,jln)); // copy ell-diagonal blocks
                    } // irn

//                      if (echo > 0) std::printf("%s charge deficit operator for ell=%c is [%g %g, %g %g]\n", // display 2x2
//                          label, ellchar[ell], ovl_ell(0,0), ovl_ell(0,nmx-1), ovl_ell(nmx-1,0), ovl_ell(nmx-1,nmx-1));

                    // the eigenvalues of the non-local part of the overlap operator may not be <= -1
                    auto const info = linear_algebra::eigenvalues(eigvals.data(), nmx, ovl_ell.data(), ovl_ell.stride());
                    if (0 != info) warn("%s when trying to diagonalize %dx%d charge deficit operator info= %i", label, nmx, nmx, info);

                    if (eigvals[0] <= -1.0) {
                        warn("%s instable eigenvalues of %c-overlap operator %g , %g and %g", label, ellchar[ell], 1 + eigvals[0], 1 + eigvals[1], 1 + eigvals[2]);
                    } else if (eigvals[0] < -0.9) {
                        warn("%s critical eigenvalues of %c-overlap operator %g , %g and %g", label, ellchar[ell], 1 + eigvals[0], 1 + eigvals[1], 1 + eigvals[2]);
                    } // warnings

                    if (echo > 1 + check_overlap_eigenvalues) {
                        std::printf("# %s eigenvalues of the %c-overlap operator are %g , %g and %g\n",
                                      label, ellchar[ell], 1 + eigvals[0], 1 + eigvals[1], 1 + eigvals[2]);
                    } // echo
#ifdef DEVEL
                    if (echo > 11 + check_overlap_eigenvalues) {
                        std::printf("# %s eigenvalues of the %c-charge deficit operator are %g , %g and %g\n",
                                      label, ellchar[ell], eigvals[0], eigvals[1], eigvals[2]);
                    } // echo
#endif // DEVEL

                } // nn[ell] > 0
            } // ell
            if (check_overlap_eigenvalues < -9) abort("# single_atom.check.overlap.eigenvalues=%d < -9", check_overlap_eigenvalues);
        } // check_overlap_eigenvalues

        if (echo > 0) { // eigenstate_analysis only writes to the log, so we can save the efforts
            double const V_rmax = potential[SMT][rg[SMT].n - 1]*rg[SMT].rinv[rg[SMT].n - 1];
            std::vector<double> Vsmt(rg[SMT].n);
            // prepare a smooth local potential which goes to zero at Rmax
            for (int ir = 0; ir < rg[SMT].n; ++ir) {
                Vsmt[ir] = potential[SMT][ir]*rg[SMT].rinv[ir] - V_rmax;
            } // ir
            if (echo > 1) std::printf("\n\n# %s %s: eigenstate_analysis\n", label, __func__);
            if (echo > 5) std::printf("# local potential for eigenstate_analysis is shifted by %g %s\n", V_rmax*eV,_eV);
            scattering_test::eigenstate_analysis // find the eigenstates of the spherical Hamiltonian
              (rg[SMT], Vsmt.data(), sigma, int(numax + 1), numax, hamiltonian_ln.data(), overlap_ln.data(), 384,
               V_rmax, label, echo, reference_spdf, control::get("eigenstate.analysis.warn", 3.675e-3)); // threshold= 0.1 eV
            std::fflush(stdout);
        } // echo

        if (echo > 0) { // logarithmic_derivative only writes to the log, so we can save the efforts
            if (echo > 1) std::printf("\n\n# %s %s: logarithmic_derivative\n", label, __func__);
            double const *const rV[TRU_AND_SMT] = {potential[TRU].data(), potential[SMT].data()};
            scattering_test::logarithmic_derivative // scan the logarithmic derivatives
              (rg, rV, sigma, int(numax + 1), numax, hamiltonian_ln.data(), overlap_ln.data(), logder_energy_range, label, echo);
            std::fflush(stdout);
        } // echo

    } // check_spherical_matrix_elements





    void total_energy_contributions(
          double E_tot[TRU_AND_SMT] // result
        , double E_kin[TRU_AND_SMT] // side-result
        , int const echo=0 // log-level
    ) const {

        set(E_kin, TRU_AND_SMT, 0.0); // init
        for (int csv = core; csv <= valence; ++csv) {
            add_product(E_kin, TRU_AND_SMT, energy_kin_csvn[csv], take_spherical_density[csv]);
        } // csv
        double const non_spherical = 1. - take_spherical_density[valence];
        // energy_kin_csvn[3] contains the kinetic energy of the density matrix
        add_product(E_kin, TRU_AND_SMT, energy_kin_csvn[3], non_spherical);
        // energy_dm contains the contraction of density matrix and hamiltonian correction
        E_kin[SMT] += energy_dm*non_spherical;

        set(E_tot, TRU_AND_SMT, E_kin);
        add_product(E_tot, TRU_AND_SMT, energy_es, 1.0); // electrostatic
        add_product(E_tot, TRU_AND_SMT, energy_xc, 1.0); // exchange_correlation

        if (echo > 1) {
            std::printf("\n# %s\n", label);
            std::printf("# %s Total energy (%s) true contribution and smooth contribution:\n", label, _eV);
            for (int csv = core; csv <= valence; ++csv) {
                if (csv_charge[csv] > 0) { // do not display core or semicore electrons if there are none
                    std::printf("# %s kinetic   %20.9f  spherical %-8s %6.1f %%\n",
                        label, energy_kin_csvn[csv][TRU]*eV, csv_name(csv), take_spherical_density[csv]*100);
                } // csv_charge
            } // csv
            // kinetic valence energy correction from the atomic density matrix:
            std::printf("# %s kinetic   %20.9f %20.9f%6.1f %%\n", label, energy_kin_csvn[3][TRU]*eV,
                                                                         energy_kin_csvn[3][SMT]*eV, non_spherical*100);
            std::printf("# %s dm                             %20.9f%6.1f %%\n", label, energy_dm*eV, non_spherical*100);
            std::printf("# %s kineticsum%20.9f %20.9f\n", label, E_kin[TRU]*eV, E_kin[SMT]*eV);
            std::printf("# %s es        %20.9f %20.9f\n", label, energy_es[TRU]*eV, energy_es[SMT]*eV);
            std::printf("# %s xc        %20.9f %20.9f\n", label, energy_xc[TRU]*eV, energy_xc[SMT]*eV);
            if (echo > 7) // so far, the xc-dc term does not contribute
            std::printf("# %s dc        %20.9f %20.9f\n", label, energy_dc[TRU]*eV, energy_dc[SMT]*eV); 
            std::printf("# %s\n", label);
            std::printf("# %s reference %20.9f\n", label, energy_ref[TRU]*eV);
        } // echo
        if (echo > 0) {
            std::printf("# %s Total     %20.9f %20.9f %s\n", label, E_tot[TRU]*eV, E_tot[SMT]*eV, _eV);
        } // echo
        if (echo > 1) {
            std::printf("# %s diff      %20.9f\n", label, (E_tot[TRU] - energy_ref[TRU])*eV);
            std::printf("# %s\n\n", label);
        } // echo
        
    } // total_energy_contributions




    void update_density(
          float const density_mixing[3] // mixing of the spherical {core,semicore,valence} density
        , int const echo=0 // log-level
        , bool const synthetic_density_matrix=false
    ) {
        if (echo > 2) std::printf("\n# %s %s Z=%g\n", label, __func__, Z_core);
        update_spherical_states(density_mixing, echo); // compute core level shifts

        if (regenerate_partial_waves) {
            // create new partial waves for the valence description
            update_energy_parameters(echo);
            update_partial_waves(echo);
            
            // update quantities derived from the partial waves
            update_kinetic_energy_deficit(echo);
            update_charge_deficit(echo);
            
            // check scattering properties for the reference potential
            check_spherical_matrix_elements(echo); 
            if (freeze_partial_waves) {
                regenerate_partial_waves = false;
                if (echo > 0) std::printf("# %s single_atom.relax.partial.waves=0\n", label);
            } // freeze_partial_waves
        } // regenerate_partial_waves

        int const nln = sho_tools::nSHO_radial(numax);
        int const lmax = std::max(ellmax_rho, ellmax_cmp);
        int const mlm = pow2(1 + lmax);
        view3D<double> rho_tensor(mlm, nln, nln, 0.0); // get memory
        get_rho_tensor(rho_tensor, density_matrix, sho_tools::order_zyx, echo, synthetic_density_matrix);
        update_full_density(rho_tensor, echo);
    } // update_density

    // ==============
    // after update_density we need to export qlm_compensator, 
    // then solve the 3D electrostatic problem,
    // and finally return here with ves_multipoles
    // ==============

    void update_potential(
          float const potential_mixing // mixing of the potential
        , double const ves_multipoles[] // multipoles of the electrostatic potential found in the 3D Poisson solver
        , int const echo=0 // log-level
    ) {
        if (echo > 2) std::printf("\n# %s %s\n", label, __func__);
        update_full_potential(potential_mixing, ves_multipoles, echo);
        total_energy_contributions(energy_tot, energy_kin, echo);
        update_matrix_elements(echo); // this line does not compile with icpc (ICC) 19.0.2.187 20190117
        perturbation_theory(echo);
    } // update_potential

    double get_total_energy(
          double E[]=nullptr
    ) const {
        if (nullptr != E) {
            E[energy_contribution::TOTAL]              = energy_tot[TRU] - energy_tot[SMT];
            E[energy_contribution::KINETIC]            = energy_kin[TRU] - energy_kin[SMT];
            E[energy_contribution::REFERENCE]          = energy_ref[TRU];
            E[energy_contribution::EXCHANGE_CORRELATION] = energy_xc[TRU] - energy_xc[SMT];
            E[energy_contribution::DOUBLE_COUNTING]      = energy_dc[TRU] - energy_dc[SMT];
            E[energy_contribution::ELECTROSTATIC]        = energy_es[TRU] - energy_es[SMT];
        } // nullptr != E

        // the reference energy is substracted to make the numeric values ...
        //                ... smaller and the total energy human readable
        return energy_tot[TRU] - energy_ref[TRU] - energy_tot[SMT];
    } // get_total_energy

    status_t get_smooth_spherical_quantity(
          double qnt[] // result array: function on an r2-grid
        , float const ar2 // r2-grid parameter factor
        , int const nr2   // r2-grid parameter size
        , char const what // which quantity to export
        , int const echo=1 // log-level
    ) const {
        // export smooth spherical functions on r2-grids

        std::vector<double> mixed_spherical;
        char   const *qnt_name{nullptr};   
        double const *qnt_vector{nullptr};
        if ('z' == what) { qnt_name = "zero_potential";  qnt_vector = zero_potential.data(); } else
        if ('c' == what) { qnt_name = "core_density";    qnt_vector = spherical_density[SMT][core]; } else
        if ('v' == what) { qnt_name = "valence_density"; qnt_vector = spherical_density[SMT][valence]; } else
        {   // create a csv-mixture of spherical_density[SMT] using take_spherical_density weights
            mixed_spherical = std::vector<double>(rg[SMT].n, 0.0);
            for (int csv = 0; csv < 3; ++csv) {
                add_product(mixed_spherical.data(), rg[SMT].n, spherical_density[SMT][csv], take_spherical_density[csv]);
            } // csv
            qnt_name = "mixed_density"; qnt_vector = mixed_spherical.data();
        }

        if (echo > 8) std::printf("# %s call transform_to_r2grid(%p, %.1f, %d, %s=%p, rg=%p)\n",
                        label, (void*)qnt, ar2, nr2, qnt_name, (void*)qnt_vector, (void*)&rg[SMT]);
        double const minval = ('z' == what) ? -9e307 : 0.0; // zero potential may be negative, densities should not
#ifdef DEVEL
        double const Y00s = Y00*(('c' == what) ? Y00 : 1); // Y00 for zero_pot and Y00^2 for densities
        if (echo > 8) {
            std::printf("\n## %s %s before filtering:\n", label, qnt_name);
            for (int ir = 0; ir < rg[SMT].n; ++ir) {
                std::printf("%g %g\n", rg[SMT].r[ir], qnt_vector[ir]*Y00s);
            } // ir
            std::printf("\n\n");
        } // echo
#endif // DEVEL

        // divide input by mask function
        std::vector<double> inp(rg[SMT].n);
        auto const r2cut = pow2(rg[SMT].rmax)*1.1, r2inv = 1./r2cut;
        for (int ir = 0; ir < rg[SMT].n; ++ir) {
            auto const r2 = pow2(rg[SMT].r[ir]);
            auto const mask = pow8(1. - pow8(r2*r2inv));
            inp[ir] = qnt_vector[ir] / mask;
        } // ir

        auto const stat = bessel_transform::transform_to_r2grid(qnt, ar2, nr2, inp.data(), rg[SMT], echo);

        // multiply output with mask function
        auto const inv_ar2 = 1./ar2;
        for (int ir2 = 0; ir2 < nr2; ++ir2) {
            auto const r2 = ir2*inv_ar2;
            auto const mask = (r2 < r2cut) ? pow8(1. - pow8(r2*r2inv)) : 0;
            qnt[ir2] = std::max(minval, qnt[ir2]) * mask;
        } // ir2

#ifdef DEVEL
        if (echo > 8) {
            std::printf("\n## %s %s  after filtering:\n", label, qnt_name);
            for (int ir2 = 0; ir2 < nr2; ++ir2) {
                std::printf("%g %g\n", std::sqrt(ir2*inv_ar2), qnt[ir2]*Y00s);
            } // ir2
            std::printf("\n\n");
        } // echo
#endif // DEVEL
        return stat;
    } // get_smooth_spherical_quantity

    radial_grid_t const* get_smooth_radial_grid(int const echo=0) const { return &rg[SMT]; }

    double get_number_of_electrons(
          char const csv='v'
    ) const {
        if ('c' == (csv | 32)) return csv_charge[core];
        if ('s' == (csv | 32)) return csv_charge[semicore];
        if ('v' == (csv | 32)) return csv_charge[valence];
        // otherwise total number of electrons
        return csv_charge[core] + csv_charge[semicore] + csv_charge[valence];
    } // get_number_of_electrons

    double const get_sigma() const { return sigma; }
    int    const get_numax() const { return numax; }

  private:
    
    void set_label(char const *chemical_symbol) {
        if (atom_id < 0) std::snprintf(label, 15, "%s", chemical_symbol); 
        else std::snprintf(label, 15, "%s#%i", chemical_symbol, atom_id); 
    } // set_label
    

    status_t perturbation_theory(int const echo=0) {
        SimpleTimer timer(__FILE__, __LINE__, __func__, echo);

        auto const lambda = control::get("single_atom.perturbation.strength", 1.);
        int        ellmax = control::get("single_atom.perturbation.ellmax", -1.);
        ellmax = std::min(ellmax, int(ellmax_pot));

        status_t status(0);
        if (ellmax < 0) return status; // no perturbation theory
        auto const & g = rg[TRU];

        std::vector<double> rV(g.n); // r*V(r)
        product(rV.data(), g.n, g.r, full_potential[TRU][00], Y00);
        rV[0] = -Z_core; // correct for singularity

        if (echo > 5) {
            std::printf("\n## r rV (a.u.):\n");
            for (int ir = 0; ir < 9; ++ir) {
                std::printf("%g %g\n", g.r[ir], rV[ir]);
            } // ir
            std::printf("\n\n");
        } // echo

        int const nlm = pow2(1 + ellmax);
        int const nln = sho_tools::nSHO_radial(numax);

        // compute eigenstates of the spherical (unperturbed) Hamiltonian

        int const nr_aligned = align<2>(g.n);
        view2D<double> rwave0(nln, nr_aligned, 0.0);
        std::vector<double> energy0(nln);
        for (int iln = 0; iln < nln; ++iln) {
            int const enn = partial_wave[iln].enn;
            int const ell = partial_wave[iln].ell;
            double E      = partial_wave[iln].energy; // alternatively take the corresponding energy from spherical_states
            auto const stat = radial_eigensolver::shooting_method(SRA, g, rV.data(), enn, ell, E, rwave0[iln]);
            if (stat) {
                if (echo > 0) std::printf("# %s failed to solve spherical problem for ell=%d\n", label, ell);
                status += std::abs(stat);
            } else {
                if (echo > 0) std::printf("# %s found %d%c-energy at%12.6f %s\n", label, enn, ellchar[ell], E*eV, _eV);
            } // stat
            // ToDo: normalize states rwave0
            auto const norm2 = dot_product(g.n, rwave0[iln], rwave0[iln], g.dr);
            if (norm2 > 0) {
                scale(rwave0[iln], g.n, 1./std::sqrt(norm2));
            } else {
                status += 1;
                warn("%s failed to normalize %d%c-state\n", label, enn, ellchar[ell]);
            } // norm2
            energy0[iln] = E; // store
        } // ell

        // start perturbation theory
        // < rwave0 | V_nonspherical | rwave0' >
        // in analogy to the construction of the true Hamiltonian elements in update_matrix_elements of single_atom.cxx

        // construct the potential tensor in terms of the emm_Degenerate partial waves
        view3D<double> potential_ln(nlm, nln, nln, 0.0); // get memory, // emm1-emm2-degenerate
        { // scope: fill potential_ln
            std::vector<double> wave_pot_rdr(g.n);
            for (int lm = 00; lm < nlm; ++lm) {
                // similar but not the same with check_spherical_matrix_elements
                for (int iln = 0; iln < nln; ++iln) {
                    product(wave_pot_rdr.data(), g.n, rwave0[iln], full_potential[TRU][lm], g.dr);
                    for (int jln = 0; jln < nln; ++jln) {
                        potential_ln(lm,iln,jln) = dot_product(g.n, wave_pot_rdr.data(), rwave0[jln]);
                    } // jln
                } // iln
            } // lm
        } // scope

//      initialize_Gaunt(); // make sure the Gaunt tensor is precomputed

        int const nlmn = sho_tools::nSHO(numax);
        view2D<double> hamiltonian_lmn(nlmn, nlmn, 0.0);
        for (int ilmn = 0; ilmn < nlmn; ++ilmn) {
            int const iln = ln_index_list[ilmn];
            hamiltonian_lmn(ilmn,ilmn) = energy0[iln]; // diagonal elements
        } // ilmn


        // add contribution of the non-spherical potential
        int const mlm = pow2(1 + numax);
        for (auto gnt : gaunt) {
            int const lm = gnt.lm, lm1 = gnt.lm1, lm2 = gnt.lm2; auto const G = gnt.G * lambda;
            if (lm > 00) { // the spherical part is already included in energy0[]
                if (lm1 < mlm && lm2 < mlm) {
                    if (lm < nlm) {
                        for (int ilmn = lmn_begin[lm1]; ilmn < lmn_end[lm1]; ++ilmn) {
                            int const iln = ln_index_list[ilmn];
                            for (int jlmn = lmn_begin[lm2]; jlmn < lmn_end[lm2]; ++jlmn) {
                                int const jln = ln_index_list[jlmn];
                                hamiltonian_lmn(ilmn,jlmn) += G * potential_ln(lm,iln,jln);
                            } // jlmn
                        } // ilmn
                    } // lm
                } // limits
            } // lm > 00
        } // gnt

        std::vector<double> energy1(nlmn, 0.0);
        // diagonalize hamiltonian_lmn
        auto const stat = linear_algebra::eigenvalues(energy1.data(), nlmn, hamiltonian_lmn.data(), hamiltonian_lmn.stride());
        if (0 == stat) {
            if (echo > 0) {
                std::printf("# %s local spectrum(%s) ", label, _eV);
                printf_vector(" %g", energy1.data(), nlmn, "\n", eV);
            } // echo

            // ToDo:
            //  - export eigenenergies
            //  - store eigenvector coefficients for density generation

        } else {
            status += std::abs(stat);
            warn("%s diagonalization failed with status=%d", label, int(stat));
        } // stat

        return status;
    } // perturbation_theory
    
    
    
    
  }; // class LiveAtom

  // instead of having a switch only onto  the first char of a string (as in atom_update),
  // we could use a switch onto int or long with these functions
  
// #define __IS_A_CXX14_COMPILER__

  inline uint64_t constexpr string2long(char const s[8]) {
#ifdef  __IS_A_CXX14_COMPILER__
      // this version allows also null-delimited strings shorter than 8 chars
      uint64_t ui64{0};
      if (nullptr != s) {
          #define instr(n,__more__) if (s[n]) { ui64 |= (uint64_t(s[n]) << (8*n)); __more__ }
          instr(0,instr(1,instr(2,instr(3,instr(4,instr(5,instr(6,instr(7,;))))))))
          #undef  instr
      } // s is a valid pointer
      return ui64;
#else  // __IS_A_CXX14_COMPILER__
      return s[0] | (uint64_t(s[1]) <<  8)
                  | (uint64_t(s[2]) << 16) 
                  | (uint64_t(s[3]) << 24)
                  | (uint64_t(s[4]) << 32)
                  | (uint64_t(s[5]) << 40)
                  | (uint64_t(s[6]) << 48)
                  | (uint64_t(s[7]) << 56);
#endif // __IS_A_CXX14_COMPILER__   
  } // string2long

  inline uint32_t constexpr string2int(char const s[4]) {
#ifdef  __IS_A_CXX14_COMPILER__
      // this version allows also null-delimited strings shorter than 4 chars
      uint32_t ui32{0};
      if (nullptr != s) {
          #define instr(n,__more__) if (s[n]) { ui32 |= (uint32_t(s[n]) << (8*n)); __more__ }
          instr(0,instr(1,instr(2,instr(3,;))))
          #undef  instr
      } // s is a valid pointer
      return ui32;
#else  // __IS_A_CXX14_COMPILER__
      return s[0] | (uint32_t(s[1]) <<  8) 
                  | (uint32_t(s[2]) << 16) 
                  | (uint32_t(s[3]) << 24);
#endif // __IS_A_CXX14_COMPILER__
  } // string2int

  // from https://hbfs.wordpress.com/2017/01/10/strings-in-c-switchcase-statements/
  inline uint64_t constexpr __mix_hash_(char const m, uint64_t const s) { return ((s << 7) + ~(s >> 3)) + ~m; }
  inline uint64_t constexpr string2hash(char const * m) { return (*m) ? __mix_hash_(*m, string2hash(m + 1)) : 0; }

  status_t test_string_switch(char const *const what, int const echo=0) {
      #define str2int string2hash
      switch ( str2int(what) ) {
          case str2int("initialize"):
          case str2int("memory cleanup"):
          case str2int("lmax qlm"):
          case str2int("lmax vlm"): // does not work with uint32_t
          case str2int("sigma cmp"):
          case str2int("x densities"):
          case str2int("core densities"):
          case str2int("valence densities"):
          case str2int("#valence electrons"):
          case str2int("#semicore electrons"):
          case str2int("#core electrons"):
          case str2int("projectors"):
          case str2int("energies"): // reserve for the export of atomic energy contributions
          case str2int("qlm charges"):
          case str2int("update"): // needs a trailing ' ' if str2int==string2long
          case str2int("hamiltonian"):
          case str2int("zero potentials"):
          case str2int("atomic density matrices"):
          case str2int("radial grids"):
            if (echo > 0) std::printf("# %s found selector what=\"%s\".\n", __func__, what);
          break; default:
            if (echo > 0) std::printf("# %s unknown selector what=\"%s\"!\n", __func__, what);
            return 1; // error
      } // switch (what)
      #undef  str2int
      return 0; // no error
  } // test_string_switch

  status_t atom_update(
        char const *const what    // selector string
      , int const natoms          // number of atoms
      , double  *const dp         // quantities (input/output) double  dp[natoms]
      , int32_t *const ip         // quantities (input/output) integer ip[natoms]
      , float   *const fp         // quantities (input)        float   fp[natoms or less]
      , double  *const *const dpp // quantities (input/output) double* dpp[natoms]
  ) {
      // single interface: the LiveAtom class is not exposed outside this compilation unit
      // instead the interface of this function alone and all functionality can be
      // handeled via the various use options of this interface.
      // A vector of LiveAtom instances is stored in a static variable.

      static std::vector<LiveAtom*> a; // internal state, so this function may not be templated!!
      static std::vector<bool> echo_mask;
      static int echo = -9;
      if (-9 == echo) echo = int(control::get("single_atom.echo", 0.)); // initialize only on the 1st call to atom_update()

      if (nullptr == what) return -1;
      
      char const how = what[0] | 32; // first char, convert to lowercase
      if (echo > 4) std::printf("\n# %s %s what=\"%s\" --> \'%c\'\n\n", __FILE__, __func__, what, how);

      float constexpr ar2_default = 16.f;
      int   constexpr nr2_default = 1 << 12;
      float const mix_defaults[] = {.5f, .5f, .5f, .5f}; // {mix_pot, mix_rho_core, mix_rho_semicore, mix_rho_valence}

      int na{natoms};

      status_t stat(0);
      stat += test_string_switch(what); // muted

      switch (how) {

          case 'i': // interface usage: atom_update("initialize", natoms, dp=Za[], ip=numax[], ion[]=0, dpp=null);
          {
              double const *Za = dp; assert(nullptr != Za); // may not be nullptr as it holds the atomic core charge Z[ia]
              a.resize(na);
              echo_mask.resize(na);
              bool const atomic_valence_density = (nullptr != dpp); // global control for all atoms
              auto const echo_init = int(control::get("single_atom.init.echo", double(echo))); // log-level for the LiveAtom constructor
              auto const bmask = int64_t(control::get("single_atom.echo.mask", -1.)); // log-level mask, -1:all
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  float const ion = (fp) ? fp[ia] : 0;
                  echo_mask[ia] = (-1 == bmask) ? 1 : ((bmask >> ia) & 0x1);
                  int type{0}; if (ip) type = (-9 == ip[ia]);
                  if (0 == type) {
                      a[ia] = new LiveAtom(Za[ia], atomic_valence_density, ia, echo_mask[ia]*echo_init, ion);
                  } else {
                      a[ia] = new LiveAtom(Za[ia], echo_mask[ia]*echo_init); // load from pawxml files
                  }
                  if (ip) ip[ia] = a[ia]->get_numax(); // export numax, optional
              } // ia
          }
          break;

          case 'm': // interface usage: atom_update("memory cleanup", natoms);
          {
              if (a.size() != na) warn("what='%s' for %d atoms, but only %ld atoms active!", what, na, a.size());
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  a[ia]->~LiveAtom(); // envoke destructor
              } // ia
              a.clear();
              angular_grid::cleanup(echo);
              na = a.size(); // set na to fulfill consistency check at the end of this routine
              assert(!dp); assert(!ip); assert(!fp); assert(!dpp); // all other arguments must be nullptr (by default)
          }
          break;

          case '#': // interface usage: atom_update("#core electrons", natoms, dp=ne[]);
                    // interface usage: atom_update("#semicore electrons", natoms, dp=ne[]);
                    // interface usage: atom_update("#valence electrons", natoms, dp=ne[]);
                    // interface usage: atom_update("#all electrons", natoms, dp=ne[]);
          {
              double *ne = dp; assert(nullptr != ne);
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  ne[ia] = a[ia]->get_number_of_electrons(what[1]);
              } // ia
              assert(!ip); assert(!fp); assert(!dpp); // all other arguments must be nullptr (by default)
          }
          break;

          case 'c': // interface usage: atom_update("core densities",    natoms, null, nr2=2^12, ar2=16.f, qnt=rho_c);
          case 'v': // interface usage: atom_update("valence densities", natoms, null, nr2=2^12, ar2=16.f, qnt=rho_v);
          case 'z': // interface usage: atom_update("zero potentials",   natoms, null, nr2=2^12, ar2=16.f, qnt=v_bar);
          {
              double *const *const qnt = dpp; assert(nullptr != qnt);
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  assert(nullptr != qnt[ia]);
                  int   const nr2 = ip ? ip[ia] : nr2_default;
                  float const ar2 = fp ? fp[ia] : ar2_default;
                  stat += a[ia]->get_smooth_spherical_quantity(qnt[ia], ar2, nr2, how);
              } // ia
              assert(!dp);
          }
          break;

          case 'r': // interface usage: atom_update("radial grid", natoms, null, null, null, (double**)rg_ptr);
          {
#ifdef  DEVEL
              assert(nullptr != dpp);
              double const **const dnc = const_cast<double const**>(dpp);
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  dnc[ia] = reinterpret_cast<double const*>(a[ia]->get_smooth_radial_grid()); // pointers to smooth radial grids
              } // ia
#else  // DEVEL
              stat = -1;
              warn("%s only available with -D DEVEL", what);
#endif // DEVEL
              assert(!dp); assert(!ip); assert(!fp); // all other arguments must be nullptr (by default)
          }
          break;

          case 's': // interface usage: atom_update("sigma compensator", natoms, sigma);
          {
              double *const sigma = dp; assert(nullptr != sigma);
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  sigma[ia] = a[ia]->sigma_compensator; // spreads of the compensators // ToDo: use a getter function
              } // ia
              assert(!ip); assert(!fp); assert(!dpp); // all other arguments must be nullptr (by default)
          }
          break;

          case 'p': // interface usage: atom_update("projectors", natoms, sigma, numax);
          {
              double  *const sigma = dp; assert(nullptr != sigma);
              int32_t *const numax = ip; assert(nullptr != numax);
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  sigma[ia] = a[ia]->get_sigma(); // spreads of the projectors
                  numax[ia] = a[ia]->get_numax(); //  number of SHO-projectors
              } // ia
              assert(!fp); assert(!dpp); // all other arguments must be nullptr (by default)
          }
          break;

          case 'u': // interface usage: atom_update("update", natoms, null, null, mix[1]={mix_pot}, vlm);
          {
              double const *const *const vlm = dpp; assert(nullptr != vlm);
              float const mix_pot = fp ? fp[0] : mix_defaults[0];
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  a[ia]->update_potential(mix_pot, vlm[ia], echo_mask[ia]*echo); // set electrostatic multipole shifts
              } // ia
              assert(!dp); assert(!ip); // all other arguments must be nullptr (by default)
          }
          break;

          case 'a': // interface usage: atom_update("atomic density matrix", natoms, null, null, mix_rho[3], atom_rho);
          {
              double const *const *const atom_rho = dpp; assert(nullptr != atom_rho);
              float const *const mix_rho = fp ? fp : &mix_defaults[1];
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  assert(nullptr != atom_rho[ia]);
                  int const numax = a[ia]->get_numax();
                  int const ncoeff = sho_tools::nSHO(numax);
                  for (int i = 0; i < ncoeff; ++i) {
                      set(a[ia]->density_matrix[i], ncoeff, &atom_rho[ia][i*ncoeff + 0]);
                  } // i
                  a[ia]->update_density(mix_rho, echo_mask[ia]*echo);                  
              } // ia
              assert(!dp); assert(!ip); // all other arguments must be nullptr (by default)
          }
          break;
              
          case 'q': // interface usage: atom_update("qlm charges", natoms, null, null, null, qlm);
          {
              double *const *const qlm = dpp; assert(nullptr != qlm);
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  int const nlm = pow2(1 + a[ia]->ellmax_cmp);
                  set(qlm[ia], nlm, a[ia]->qlm_compensator.data()); // copy compensator multipoles
              } // ia
              assert(!dp); assert(!ip); assert(!fp); // all other arguments must be nullptr (by default)
          }
          break;

          case 'l': // interface usage: atom_update("lmax qlm", natoms, dp=null, lmax, fp=mix_spherical=null);
                    // interface usage: atom_update("lmax vlm", natoms, dp= ~0 , lmax, fp=mix_spherical=null);
          {
              int32_t *const lmax = ip; assert(nullptr != lmax);
              float const mix_spherical = fp ? std::min(std::max(0.f, fp[0]), 1.f) : 0;
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  lmax[ia] = dp ? a[ia]->ellmax_pot : a[ia]->ellmax_cmp;
                  // fine-control take_spherical_density[valence] any float in [0, 1], NOT atom-resolved! consumes only fp[0]
                  if (fp) a[ia]->take_spherical_density[valence] = mix_spherical;
              } // ia
              assert(!dpp); // last argument must be nullptr (by default)
          }
          break;

          case 'n': // interface usage: atom_update("numax", natoms, null, numax);
          {
              int32_t *const numax = ip; assert(nullptr != numax);
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  numax[ia] = a[ia]->get_numax();
              } // ia
              assert(!dp); assert(!fp); assert(!dpp); // all other arguments must be nullptr (by default)
          }
          break;
          
          case 'h': // interface usage: atom_update("hamiltonian and overlap", natoms, null, nelements, null, atom_mat);
          {
              double *const *const atom_mat = dpp; assert(nullptr != atom_mat);
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  assert(nullptr != atom_mat[ia]);
                  int const numax = a[ia]->get_numax();
                  int const ncoeff = sho_tools::nSHO(numax);
                  int const nelements = ip ? ip[ia] : 2*pow2(ncoeff);
                  int const h1hs2 = nelements/pow2(ncoeff); // ToDo: make this more elegant
                  for (int i = 0; i < ncoeff; ++i) {
                      if (h1hs2 > 1)
                      set(&atom_mat[ia][(1*ncoeff + i)*ncoeff + 0], ncoeff, a[ia]->overlap[i]);
                      set(&atom_mat[ia][(0*ncoeff + i)*ncoeff + 0], ncoeff, a[ia]->hamiltonian[i]);
                  } // i
              } // ia
              assert(!dp); assert(!fp); // all other arguments must be nullptr (by default)
          }
          break;

          case 'e': // interface usage: atom_update("energies", natoms, delta_Etot, null, null, dpp=atom_ene=null);
          {
              double *const *const atom_ene = dpp;
              for (size_t ia = 0; ia < a.size(); ++ia) {
                  dp[ia] = a[ia]->get_total_energy(atom_ene ? atom_ene[ia] : nullptr);
              } // ia
              assert(!ip); assert(!fp); // all other arguments must be nullptr (by default)
          }
          break;
          
          default:
          {
              if (echo > 0) std::printf("# %s: first argument \'%s\' undefined, no action!\n", __func__, what);
              stat = how;
          }
          break;

      } // switch(how)

      if (a.size() != na) warn("inconsistency: internally %ld atoms active but natoms=%d", a.size(), na);
      if (stat) warn("what='%s' returns status = %i", what, int(stat));
      return stat;
  } // atom_update


#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_compensator_normalization(int const echo=5) {
      double maxdev{0};
#ifdef DEVEL
      if (echo > 1) std::printf("\n# %s: %s\n", __FILE__, __func__);
      auto & rg = *radial_grid::create_radial_grid(512, 2.f);
      int const nr = rg.n, lmax = 0, nlm = pow2(1 + lmax);
      std::vector<double> qlm(nlm, 0.0);
      view2D<double> cmp(1, nr);
      for (double sigma = 0.5; sigma < 2.25; sigma *= 1.1) {
          set(qlm.data(), nlm, 0.0); qlm[0] = 1.0;
          set(cmp, 1, 0.0); // clear
          add_or_project_compensators<0>(cmp, qlm.data(), rg, lmax, sigma, 0); // add normalized compensator
//        add_or_project_compensators<1>(cmp, qlm.data(), rg, lmax, sigma, 0); // project
//        if (echo > 0) std::printf("# %s: square-norm of normalized compensator with sigma = %g is %g\n", __func__, sigma, qlm[0]);
          add_or_project_compensators<3>(cmp, qlm.data(), rg, lmax, sigma, 0); // test normalization
          maxdev = std::max(maxdev, std::abs(qlm[0] - 1.0));
          if (echo > 7) std::printf("# %s: for sigma = %g is 1 + %.1e\n", __func__, sigma, qlm[0] - 1);
      } // sigma
      if (echo > 2) std::printf("# %s: largest deviation is %.1e\n", __func__, maxdev);
      radial_grid::destroy_radial_grid(&rg);
#endif // DEVEL
      return (maxdev > 2e-15);
  } // test_compensator_normalization

  status_t test_LiveAtom(int const echo=9) {
      bool const avd     = control::get("single_atom.test.atomic.valence.density", 1.) > 0; // atomic valence density
      auto const ionic   = control::get("single_atom.test.ion", 0.); // only with HAS_ELEMENT_CONFIG
      auto const Z_begin = control::get("single_atom.test.Z", 29.); // default copper
      auto const Z_inc   = control::get("single_atom.test.Z.inc", 1.); // default: sample only integer values
      auto const Z_end   = control::get("single_atom.test.Z.end", Z_begin + Z_inc); // default: only one atom
      for (double Z = Z_begin; Z < Z_end; Z += Z_inc) {
          if (echo > 3) std::printf("\n\n\n\n\n\n\n\n\n\n\n\n");
          if (echo > 1) std::printf("\n# %s: Z = %g\n", __func__, Z);
          LiveAtom a(Z, avd, -1, echo, ionic); // envoke constructor
          if (echo > 3) std::printf("\n\n\n\n\n\n\n\n\n\n\n\n");
      } // Z
      if (Z_begin >= Z_end) warn("Empty range for Z in [%g, %g]", Z_begin, Z_end);
      return 0;
  } // test_LiveAtom

  status_t test_pawxml_constructor(int const echo=9) {
      auto const Z = control::get("single_atom.test.Z", 29.); // default copper
      LiveAtom a(Z, echo); // envoke constructor from pawxml
      return 0;
  } // test_pawxml_constructor

  status_t all_tests(int const echo) {
      status_t stat(0);
      int n{0}; int const t = control::get("single_atom.select.test", -1.); // -1:all
      if (t & (1 << n++)) stat += test_compensator_normalization(echo);
      if (t & (1 << n++)) stat += test_pawxml_constructor(echo);
      if (t & (1 << n++)) stat += test_LiveAtom(echo);
      assert(3 == n);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace single_atom


#define LIVE_ATOM_AS_LIBRARY
#ifdef  LIVE_ATOM_AS_LIBRARY
  // C and Fortran interface
  extern "C" {
      #include <cstdint> // and <stdint.h> in C
      #define SINGLE_ATOM_SOURCE
        #include "single_atom.h"
      #undef  SINGLE_ATOM_SOURCE
  } // extern "C"
#endif // LIVE_ATOM_AS_LIBRARY
