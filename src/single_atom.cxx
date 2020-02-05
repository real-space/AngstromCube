#include <cstdio> // printf
#include <cstdlib> // abs
#include <cmath> // sqrt
#include <cassert> // assert
#include <algorithm> // max

#include "single_atom.hxx"

#include "radial_grid.hxx" // create_default_radial_grid, destroy_radial_grid, dot_product
#include "radial_eigensolver.hxx" // shooting_method
#include "radial_potential.hxx" // Hartree_potential
#include "angular_grid.hxx" // ::transform, ::Lebedev_grid_size
#include "radial_integrator.hxx" // integrate_outwards
#include "exchange_correlation.hxx" // lda_PZ81_kernel
#include "inline_tools.hxx" // align<nbits>
#include "sho_unitary.hxx" // Unitary_SHO_Transform<real_t>
#include "solid_harmonics.hxx" // lm_index, Y00, Y00inv
#include "atom_core.hxx" // initial_density, rad_pot, nl_index
#include "sho_tools.hxx" // lnm_index, nSHO
#include "sho_radial.hxx" // nSHO_radial
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t, emm_Degenerate, spin_QN_t, spin_Degenerate
#include "atom_core.hxx" // ::ellchar
#include "energy_level.hxx" // TRU, SMT, TRU_AND_SMT, core_level_t, valence_level_t
#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // pow2, pow3, set, scale, product, add_product, intpow
#include "simple_math.hxx" // invert
#include "simple_timer.hxx" // SimpleTimer
#include "bessel_transform.hxx" // transform_to_r2_grid
#include "scattering_test.hxx" // eigenstate_analysis, emm_average
#include "linear_algebra.hxx" // linear_solve, generalized_eigenvalues
#include "data_view.hxx" // view2D<T>, view3D<T>
#include "control.hxx" // get

// #define FULL_DEBUG
#define DEBUG

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

extern "C" {
   // BLAS interface to matrix matrix multiplication
  void dgemm_(const char*, const char*, const int*, const int*, const int*, const double*,
              const double*, const int*, const double*, const int*, const double*, double*, const int*);
} // extern "C"

  status_t minimize_curvature(int const n, view2D<double> & A, view2D<double> & B, double* lowest=nullptr) {
      // solve a generalized eigenvalue problem
      auto eigs = std::vector<double>(n);
      auto const info = linear_algebra::generalized_eigenvalues(n, A.data(), A.stride(),
                                                                   B.data(), B.stride(), eigs.data());
      if (lowest) *lowest = eigs[0];
      return info;
  } // minimize_curvature


  int constexpr ELLMAX=7;
  char const ellchar[] = "spdfghijklmno";
  char const c0s1v2_char[] = "csv?";  // 0:core, 1:semicore, 2:valence
  char const c0s1v2_name[][8] = {"core", "semi", "valence"};  // 0:core, 1:semicore, 2:valence

  double constexpr Y00 = solid_harmonics::Y00; // == 1./sqrt(4*pi)
  double constexpr Y00inv = solid_harmonics::Y00inv; // == sqrt(4*pi)

  template<typename real_t>
  inline void symmetrize(real_t &left, real_t &right) { left = (left + right)/2; right = left; }

  status_t pseudize_function(double fun[], radial_grid_t const *rg, int const irc,
              int const nmax=4, int const ell=0, double *coeff=nullptr) {
      int constexpr echo = 0;
      // match a radial function with an even-order polynomial inside r[irc]
      double Amat[4][4], bvec[4] = {0,0,0,0};
      set(Amat[0], 4*4, 0.0);
      int const nm = std::min(std::max(1, nmax), 4);
      for(int i4 = 0; i4 < nm; ++i4) {
          // use up to 4 radial indices [irc-2,irc-1,irc+0,irc+1]
          int const ir = irc + i4 - nm/2; // not fully centered around irc
          double const r = rg->r[ir];
          double rl = intpow(r, ell);
          // set up a basis of 4 functions: r^0, r^0, r^4, r^6 at up to 4 neighboring grid points around r(irc)
          for(int j4 = 0; j4 < nm; ++j4) {
              Amat[j4][i4] = rl;
              rl *= pow2(r);
          } // j4
          // b is the inhomogeneus right side of the set of linear equations
          bvec[i4] = fun[ir];
      } // i4

      if (echo > 7) {
          for(int i4 = 0; i4 < nm; ++i4) {
              printf("# %s Amat i=%i ", __func__, i4);
              for(int j4 = 0; j4 < nm; ++j4) {
                  printf(" %16.9f", Amat[j4][i4]);
              } // j4
              printf(" bvec %16.9f\n", bvec[i4]);
          } // i4
      } // echo

      double* x = bvec;
      auto const info = linear_algebra::linear_solve(nm, Amat[0], 4, bvec, 4, 1);

      if (echo > 7) {
          printf("# %s xvec     ", __func__);
          for(int j4 = 0; j4 < nm; ++j4) {
              printf(" %16.9f", x[j4]);
          }   printf("\n");
      } // echo

      set(x + nm, 4 - nm, 0.0); // clear the unused coefficients
      // replace the inner part of the function by the even-order polynomial
      for(int ir = 0; ir < irc; ++ir) {
          double const r = rg->r[ir];
          double const rr = pow2(r);
          double const rl = intpow(r, ell);
          fun[ir] = rl*(x[0] + rr*(x[1] + rr*(x[2] + rr*x[3])));
      } // ir
      if (nullptr != coeff) set(coeff, nm, x); // export coefficients
      return info;
  } // pseudize_function



    template<int ADD0_or_PROJECT1>
    void add_or_project_compensators(
    	  view2D<double> & Alm // ADD0_or_PROJECT1 == 0 or 2 result
    	, double qlm[]         // ADD0_or_PROJECT1 == 1 or 3 result
    	, int const lmax	   // cutoff for anuular momentum expansion
    	, radial_grid_t const *rg // radial grid despriptor
    	, double const sigma   // spread_compensator
    	, int const echo=0     // log output level
    ) {
        int const nr = rg->n;
        auto const sig2inv = -.5/(sigma*sigma);
        if (echo > 0) printf("# sigma = %g\n", sigma);
        std::vector<double> rlgauss(nr);
        std::vector<double> rl(nr);
        for(int ell = 0; ell <= lmax; ++ell) { // serial!
            double norm = 0;
            for(int ir = 0; ir < nr; ++ir) {
                auto const r = rg->r[ir];
                if (0 == ell) {
                    rl[ir] = 1; // start with r^0
                    rlgauss[ir] = std::exp(sig2inv*r*r);
                } else {
                    rl[ir]      *= r; // construct r^ell
                    rlgauss[ir] *= r; // construct r^ell*gaussian
                }
                norm += rlgauss[ir] * rl[ir] * rg->r2dr[ir];
                if (echo > 8) printf("# ell=%i norm=%g ir=%i rlgauss=%g rl=%g r2dr=%g\n",
                						ell, norm, ir, rlgauss[ir], rl[ir], rg->r2dr[ir]);
            } // ir
            if (echo > 1) printf("# ell=%i norm=%g nr=%i\n", ell, norm, nr);
            assert(norm > 0);
            auto const scal = 1./norm;
            for(int emm = -ell; emm <= ell; ++emm) {
                int const lm = solid_harmonics::lm_index(ell, emm);
                if (0 == ADD0_or_PROJECT1) {
                    add_product(Alm[lm], nr, rlgauss.data(), qlm[lm]*scal); // add normalized compensator to augmented density
                } else if (2 == ADD0_or_PROJECT1) {
                    add_product(Alm[lm], nr, rl.data(), qlm[lm]); // add q_{\ell m} * r^\ell to electrostatic potential
                    // ToDo: setup of normalized rlgauss is not necessary in this version
                } else if (3 == ADD0_or_PROJECT1) {
                    qlm[lm] = dot_product(nr, Alm[lm], rl.data(), rg->r2dr); // dot_product with metric
                    // ToDo: setup of normalized rlgauss is not necessary in this version
                } else {
                    qlm[lm] = dot_product(nr, Alm[lm], rlgauss.data(), rg->r2dr) * scal; // dot_product with metric
                } // add or project
            } // emm
        } // ell
    } // add_or_project_compensators


  class LiveAtom {
  public:
      // ToDo: separate everything which is energy-parameter-set dependent and group it into a class valence_t

      // general config
      int32_t atom_id; // global atom identifyer
      double Z_core; // number of nucleons in the core
      char label[16]; // label string
      radial_grid_t* rg[TRU_AND_SMT]; // radial grid descriptor for the true and smooth grid:
              // SMT may point to TRU, but at least both radial grids must have the same tail
      int nr_diff; // how many more radial grid points are in *rg[TRU] compared to *rg[SMT]
      ell_QN_t ellmax; // limit ell for full_potential and full_density
      ell_QN_t ellmax_compensator; // limit ell for the charge deficit compensators
      double sigma_compensator; // Gaussian spread for the charge deficit compensators
      std::vector<double> qlm_compensator; // coefficients for the charge deficit compensators, (1+ellmax_compensator)^2

      // the following quantities are energy-parameter-set dependent
      double r_cut; // classical augmentation radius for potential and core density
      int ir_cut[TRU_AND_SMT]; // classical augmentation radius index for potential and core density
      float r_match; // radius for matching of true and smooth partial wave, usually 6--9*sigma
      ell_QN_t numax; // limit of the SHO projector quantum numbers
      double  sigma, sigma_inv; // spread of the SHO projectors and its inverse
      uint8_t nn[1+ELLMAX+2]; // number of projectors and partial waves used in each ell-channel
      std::vector<valence_level_t> partial_wave;
      // the following quantities are energy-parameter-set dependent and spin-resolved (nspins=1 or =2)
      std::vector<double> zero_potential; // PAW potential shape correction
      view2D<double> hamiltonian, overlap; // matrices [nSHO][>=nSHO]
      view3D<double> kinetic_energy; // tensor [TRU_AND_SMT][nln][nln]
      view4D<double> charge_deficit; // tensor [1 + ellmax_compensator][TRU_AND_SMT][nln][nln]
      view2D<double> projectors; // [nln][nr_smt] r*projectors
      view3D<double> partial_waves[TRU_AND_SMT]; // matrix [wave0_or_wKin1][nln][nr], valence states point into this
      view3D<double> true_core_waves; // matrix [wave0_or_wKin1][nln][nr], core states point into this
      std::vector<double> true_norm; // vector[nln] for display of partial wave results
      // end of energy-parameter-set dependent members

      view2D<double> aug_density; // augmented density, core + valence + compensation, (1+ellmax)^2 radial functions
      int ncorestates; // for emm-Degenerate representations, 20 (or 32 with spin-orbit) core states are maximum
      int nspins; // 1 or 2 or 4, order: 0zxy

      view2D<double> unitary_Ezyx_lmn; // unitary sho transformation matrix [order_Ezyx][order_lmn], stride=nSHO(numax)

      double logder_energy_range[3]; // [start, increment, stop]

      // spin-resolved members
      std::vector<core_level_t> core_state; // 20 core states are the usual max., 32 core states are enough if spin-orbit-interaction is on
      std::vector<double> core_density[TRU_AND_SMT]; // spherical core density*4pi, no Y00 factor
      std::vector<double> spherical_valence_density[TRU_AND_SMT]; // spherical valence density*4pi, no Y00 factor
      view2D<double> full_density[TRU_AND_SMT]; // total density, core + valence, (1+ellmax)^2 radial functions
      view2D<double> full_potential[TRU_AND_SMT]; // (1+ellmax)^2 radial functions
      std::vector<double> potential[TRU_AND_SMT]; // spherical potential r*V(r), no Y00 factor, used for the generation of partial waves


      double core_charge_deficit; // in units of electrons

      bool gaunt_init;
      std::vector<gaunt_entry_t> gaunt;
      // ToDo: check if all of these 4 lists are used
      std::vector<int16_t> ln_index_list;
      std::vector<int16_t> lm_index_list;
      std::vector<int16_t> lmn_begin;
      std::vector<int16_t> lmn_end;

  public:

    // constructor method:
    LiveAtom(double const Z_nucleons
            , int const nu_max=3
            , bool const transfer2valence=true // depending on transfer2valence results look
    // slightly different in the shape of smooth potentials but matrix elements are the same
            , double const ionization=0
            , int const global_atom_id=-1
            , int const echo=0 // logg level for this contructor method only
            ) : atom_id{global_atom_id} 
              , Z_core{Z_nucleons}
              , gaunt_init{false}
        { // constructor

        if (atom_id >= 0) { sprintf(label, "a#%d", atom_id); } else { label[0]=0; }
        if (echo > 0) printf("\n\n#\n# %s LiveAtom with %g nucleons, ionization=%g\n", label, Z_core, ionization);

        rg[TRU] = radial_grid::create_default_radial_grid(Z_core);

        if (0) { // flat copy, true and smooth quantities live on the same radial grid
            rg[SMT] = rg[TRU]; rg[SMT]->memory_owner = false; // avoid double free
        } else { // create a radial grid descriptor which has less points at the origin
            auto const rc = control::get("smooth.radial.grid.from", 1e-4);
            rg[SMT] = radial_grid::create_pseudo_radial_grid(*rg[TRU], rc);
        } // use the same number of radial grid points for true and smooth quantities
        // Warning: *rg[TRU] and *rg[SMT] need an explicit destructor call

        int const nrt = align<2>(rg[TRU]->n), // optional memory access alignment
                  nrs = align<2>(rg[SMT]->n);
        if (echo > 0) printf("# %s radial grid numbers are %d and %d\n", label, rg[TRU]->n, rg[SMT]->n);
        if (echo > 0) printf("# %s radial grid numbers are %d and %d (padded to align)\n", label, nrt, nrs);

        numax = nu_max; // 3; // 3:up to f-projectors
        if (echo > 0) printf("# %s projectors and partial waves are expanded up to numax = %d\n", label,  numax);
        ellmax = 2*numax; // can be smaller than 2*numax
        if (echo > 0) printf("# %s radial density and potentials are expanded up to lmax = %d\n", label, ellmax);
        ellmax_compensator = std::min(4, (int)ellmax);
        if (echo > 0) printf("# %s compensation charges are expanded up to lmax = %d\n", label, ellmax_compensator);
        r_cut = 2.0; // Bohr
        sigma_compensator = r_cut/std::sqrt(20.); // Bohr
        sigma_compensator = 1.5*r_cut/std::sqrt(20.); // Bohr
//      sigma_compensator *= 3; // 3x as large as usually taken in GPAW, much smoother
        sigma = 0.61; // Bohr, spread for projectors (Cu)
        sigma_inv = 1./sigma;
        r_match = 9*sigma;
        set(nn, 1+ELLMAX+2, uint8_t(0));
        if (echo > 0) printf("# %s numbers of projectors ", label);
        int const nn_limiter = control::get("single_atom.nn.limit", 9);
        for(int ell = 0; ell <= ELLMAX; ++ell) {
            nn[ell] = std::max(0, (numax + 2 - ell)/2);
            nn[ell] = std::min((int)nn[ell], nn_limiter); // take a smaller numer of partial waves
            if (echo > 0) printf(" %d", nn[ell]);
        } // ell
        if (echo > 0) printf("\n");


        int const nlm = pow2(1 + ellmax);
        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            int const nr = (TRU == ts)? nrt : nrs;
            // spherically symmetric quantities:
            core_density[ts] = std::vector<double>(nr, 0.0); // get memory
            potential[ts]    = std::vector<double>(nr, 0.0); // get memory
            // quantities with lm-resolution:
            int const mr = align<2>(nr);
            full_density[ts]   = view2D<double>(nlm, mr, 0.0); // get memory
            full_potential[ts] = view2D<double>(nlm, mr, 0.0); // get memory
        } // true and smooth

        atom_core::read_Zeff_from_file(potential[TRU].data(), *rg[TRU], Z_core, "pot/Zeff", -1, echo);

//         // show the loaded Zeff(r)
//         if (echo > 0) {
//            printf("\n## loaded Z_eff(r) function:\n");
//            for(int ir = 0; ir < rg[TRU]->n; ++ir) {
//                printf("%.15g %.15g\n", rg[TRU]->r[ir], -potential[TRU][ir]);
//            } // ir
//         } // echo

        std::vector<int8_t> c0s1v2(96, -1); // 0:core, 1:semicore, 2:valence
        std::vector<int8_t> as_valence(96, -1);

        int const inl_core_hole = control::get("core.hole.index", -1.);;
        double const core_hole_charge = std::min(std::max(0.0, control::get("core.hole.charge", 1.)), 1.0);
        double core_hole_charge_used{0};

        double csv_charge[] = {0, 0, 0};
        
        double const core_valence_separation  = control::get("core.valence.separation", -2.0);
        double       core_semicore_separation = control::get("core.semicore.separation",    core_valence_separation);
        double const semi_valence_separation  = control::get("semicore.valence.separation", core_valence_separation);
        if (core_semicore_separation > semi_valence_separation) {
            warn("%s core.semicore.separation=%g may not be higher than semicore.valence.separation=%g, correct", 
                 label, core_semicore_separation, semi_valence_separation);
            core_semicore_separation = semi_valence_separation; // correct --> no semicore states possible
        } // warning
        int enn_core_ell[12] = {0,0,0,0, 0,0,0,0, 0,0,0,0};
        ncorestates = 20;
        true_core_waves = view3D<double>(2, ncorestates, nrt, 0.0); // get memory for the true radial wave functions and kinetic waves
        int constexpr SRA = 1;
        core_state = std::vector<core_level_t>(ncorestates);
        {
            int ics = 0, jcs = -1; 
            double n_electrons = Z_core - ionization; // init number of electrons to be distributed
            for(int nq_aux = 0; nq_aux < 8; ++nq_aux) { // auxiliary quantum number
                int enn = (nq_aux + 1)/2; // init principal quantum number
                for(int ell = nq_aux/2; ell >= 0; --ell) { // angular momentum character
                    ++enn; // update principal quantum number
                    for(int jj = 2*ell; jj >= 2*ell; jj -= 2) {
                        auto &cs = core_state[ics]; // abbreviate
                        cs.wave[TRU] = true_core_waves(0,ics); // the true radial function
                        cs.wKin[TRU] = true_core_waves(1,ics); // the kinetic energy wave
                        double E = atom_core::guess_energy(Z_core, enn);
                        std::vector<double> r2rho(nrt, 0.0);
                        radial_eigensolver::shooting_method(SRA, *rg[TRU], potential[TRU].data(),
                                 enn, ell, E, cs.wave[TRU], r2rho.data());
                        cs.energy = E;

                        int const inl = atom_core::nl_index(enn, ell);
                        if (E > semi_valence_separation) {
                            as_valence[inl] = ics; // mark as good for the valence band
                            c0s1v2[inl] = 2; // mark as valence state
//                          printf("# as_valence[nl_index(enn=%d, ell=%d) = %d] = %d\n", enn, ell, inl, ics);
                        } else // move to the valence band
                        if (E > core_semicore_separation) {
                            c0s1v2[inl] = 1; // mark as semicore state
                        } else { // move to the semicore band
                            c0s1v2[inl] = 0; // mark as core state
                        } // stay in the core

                        cs.nrn[TRU] = enn - ell - 1; // true number of radial nodes
                        cs.enn = enn;
                        cs.ell = ell;
                        cs.emm = emm_Degenerate;
                        cs.spin = spin_Degenerate;
                        double const max_occ = 2*(jj + 1);
                        double const core_hole = core_hole_charge*(inl_core_hole == inl);
                        double const occ_no_core_hole = std::min(std::max(0., n_electrons), max_occ);
                        double const occ = std::min(std::max(0., occ_no_core_hole - core_hole), max_occ) ;
                        double const real_core_hole_charge = occ_no_core_hole - occ;
                        if (real_core_hole_charge > 0) {
                            core_hole_charge_used = real_core_hole_charge;
                            if (echo > 1) printf("# %s use a core.hole.index=%i (%d%c) missing core.hole.charge=%g electrons\n", 
                                                    label, inl, enn, ellchar[ell], real_core_hole_charge);
                        } // core hole active
                        cs.occupation = occ;
                        if (occ > 0) {
                            jcs = ics; // store the index of the highest occupied core state
                            if (echo > 0) printf("# %s %s %2d%c%6.1f E= %g %s\n",
                                                    label, c0s1v2_name[c0s1v2[inl]], enn, ellchar[ell], occ, E*eV,_eV);
                            if (as_valence[inl] < 0) {
                                enn_core_ell[ell] = std::max(enn, enn_core_ell[ell]);
                            } else {
                                if (transfer2valence) cs.occupation = -occ; // mark as valence state by negative occupations
                            } // not as valence
                            csv_charge[c0s1v2[inl]] += occ;
                        } // occupied
                        if (cs.occupation > 0) {
                            double const norm = occ/dot_product(rg[TRU]->n, r2rho.data(), rg[TRU]->dr);
                            add_product(core_density[TRU].data(), rg[TRU]->n, r2rho.data(), norm);
                        } // occupied

                        n_electrons -= max_occ;
                        ++ics;
                    } // jj
                } // ell
            } // nq_aux
            ncorestates = jcs + 1; // correct the number of core states to those occupied
        } // core states
        double const total_n_electrons = csv_charge[0] + csv_charge[1] + csv_charge[2];
        if (echo > 2) printf("# %s initial occupation with %g electrons: %g core, %g semi-core and %g valence electrons\n", 
                                    label, total_n_electrons, csv_charge[0], csv_charge[1], csv_charge[2]);
        
        if ((inl_core_hole >= 0) && (std::abs(core_hole_charge_used - core_hole_charge) > 5e-16)) {
            warn("%s core.hole.charge=%g requested in core.hole.index=%i but used %g electrons (diff %.1e)",
                  label, core_hole_charge, inl_core_hole, core_hole_charge_used, core_hole_charge - core_hole_charge_used);
        } // warning when deviates

        if (echo > 3) {
            if (core_semicore_separation >= semi_valence_separation) {
                printf("# %s core.valence.separation at %g %s\n", label, core_valence_separation*eV,_eV);
            } else {
                printf("# %s core.semicore.separation at %g %s\n", label, core_semicore_separation*eV, _eV);
                printf("# %s semicore.valence.separation at %g %s\n", label, semi_valence_separation*eV, _eV);
            }
        } // echo
        
        return;
        
        scale(core_density[TRU].data(), rg[TRU]->n, rg[TRU]->rinv); // initial_density produces r^2*rho --> reduce to r*rho
        scale(core_density[TRU].data(), rg[TRU]->n, rg[TRU]->rinv); // initial_density produces r^2*rho --> reduce to   rho
        if (echo > 2) printf("# %s initial core density has %g electrons\n", label, dot_product(rg[TRU]->n, core_density[TRU].data(), rg[TRU]->r2dr));
        
        if (echo > 5) printf("# %s enn_core_ell  %i %i %i %i\n", label, enn_core_ell[0], enn_core_ell[1], enn_core_ell[2], enn_core_ell[3]);

        int const nln = sho_radial::nSHO_radial(numax); // == (numax*(numax + 4) + 4)/4
        partial_waves[TRU] = view3D<double>(2, nln, nrt, 0.0); // get memory for the true   radial wave function and kinetic wave
        partial_waves[SMT] = view3D<double>(2, nln, nrs, 0.0); // get memory for the smooth radial wave function and kinetic wave
        partial_wave = std::vector<valence_level_t>(nln);
        {
//          if (echo > 0) printf("# valence "); // no new line, compact list follows
            for(int ell = 0; ell <= numax; ++ell) {
//              for(int nrn = 0; nrn < nn[ell]; ++nrn) { // smooth number or radial nodes
                for(int nrn = 0; nrn <= (numax - ell)/2; ++nrn) { // smooth number or radial nodes

                    int const iln = sho_tools::ln_index(numax, ell, nrn);
                    assert(iln < nln);
                    auto &vs = partial_wave[iln]; // abbreviate
                    int const enn = std::max(ell + 1, enn_core_ell[ell] + 1) + nrn;
//                  if (echo > 0) printf(" %d%c", enn, ellchar[ell]);
                    vs.nrn[TRU] = enn - ell - 1; // true number of radial nodes
                    vs.nrn[SMT] = nrn;
                    vs.occupation = 0;
                    vs.enn = enn;
                    vs.ell = ell;
                    vs.emm = emm_Degenerate;
                    vs.spin = spin_Degenerate;

                    vs.wave[SMT] = partial_waves[SMT](0,iln); // the smooth radial function
                    vs.wave[TRU] = partial_waves[TRU](0,iln); // the true radial function
                    vs.wKin[SMT] = partial_waves[SMT](1,iln); // the smooth kinetic energy
                    vs.wKin[TRU] = partial_waves[TRU](1,iln); // the true kinetic energy

                    if (nrn < nn[ell]) {
                        double E = std::max(atom_core::guess_energy(Z_core, enn), core_valence_separation);
                        if (nrn > 0) E = std::max(E, partial_wave[iln - 1].energy); // higher than previous energy
                        radial_eigensolver::shooting_method(SRA, *rg[TRU], potential[TRU].data(), enn, ell, E, vs.wave[TRU]);
                        vs.energy = E;

                        {
                            int const inl = atom_core::nl_index(enn, ell);
                            int const ics = as_valence[inl];
    //                      printf("# as_valence[nl_index(enn=%d, ell=%d) = %d] = %d\n", enn, ell, inl, ics);
                            if (ics >= 0) { // atomic eigenstate was marked as valence
                                if (transfer2valence) {
                                    auto const occ = -core_state[ics].occupation; //
                                    vs.occupation = occ;
//                                  core_state[ics].occupation = 0;
                                    if (occ > 0) printf("# %s transfer %.1f electrons from %d%c-core state #%d"
                                            " to valence state #%d\n", label, occ, enn, ellchar[ell], ics, iln);
                                } // transfer2valence
                            } // ics
                        }
                        if (echo > 0) printf("# %s valence %2d%c%6.1f E = %g %s\n", label, enn, ellchar[ell], vs.occupation, E*eV,_eV);
                    } // nrn < nn[ell]
                } // nrn
            } // ell
        } // valence states

        nr_diff = rg[TRU]->n - rg[SMT]->n;
        ir_cut[SMT] = radial_grid::find_grid_index(*rg[SMT], r_cut);
        ir_cut[TRU] = ir_cut[SMT] + nr_diff;
        if (echo > 0) printf("# %s pseudize the core density at r[%i or %i] = %.6f, requested %.3f %s\n",
                          label, ir_cut[TRU], ir_cut[SMT], rg[SMT]->r[ir_cut[SMT]]*Ang, r_cut*Ang, _Ang);
        assert(rg[SMT]->r[ir_cut[SMT]] == rg[TRU]->r[ir_cut[TRU]]); // should be exactly equal

        int const nlm_aug = pow2(1 + std::max(ellmax, ellmax_compensator));
        aug_density = view2D<double>(nlm_aug, full_density[SMT].stride(), 0.0); // get memory
        int const nlm_cmp = pow2(1 + ellmax_compensator);
        qlm_compensator = std::vector<double>(nlm_cmp, 0.0); // get memory
        charge_deficit = view4D<double>(1 + ellmax_compensator, TRU_AND_SMT, nln, nln, 0.0); // get memory
        kinetic_energy = view3D<double>(TRU_AND_SMT, nln, nln, 0.0); // get memory
        zero_potential = std::vector<double>(nrs, 0.0); // get memory
        true_norm      = std::vector<double>(nln, 1.0); // get memory

        projectors = view2D<double>(nln, align<2>(rg[SMT]->n), 0.0); // get memory

        int const nSHO = sho_tools::nSHO(numax);
        int const matrix_stride = align<2>(nSHO); // 2^<2> doubles = 32 Byte alignment
        if (echo > 0) printf("# %s matrix size for hamiltonian and overlap: dim = %d, stride = %d\n", label, nSHO, matrix_stride);
        hamiltonian = view2D<double>(nSHO, matrix_stride, 0.0); // get memory
        overlap     = view2D<double>(nSHO, matrix_stride, 0.0); // get memory

        unitary_Ezyx_lmn = view2D<double>(nSHO, nSHO, 0.0);
        {   sho_unitary::Unitary_SHO_Transform<double> const u(numax);
            auto const stat = u.construct_dense_matrix(unitary_Ezyx_lmn.data(), numax, nSHO, sho_tools::order_Ezyx, sho_tools::order_lmn);
            assert(0 == stat);
        } // scope to fill unitary

        int const mlm = pow2(1 + numax);
        ln_index_list.resize(nSHO);
        lm_index_list.resize(nSHO);
        lmn_begin.resize(mlm);
        lmn_end.resize(mlm);
        get_valence_mapping(ln_index_list.data(), lm_index_list.data(), lmn_begin.data(), lmn_end.data(), echo);

        logder_energy_range[0] = control::get("logder.start", -2.0);
        logder_energy_range[1] = control::get("logder.step",  1e-2); // ToDo: these should be moved to the main function
        logder_energy_range[2] = control::get("logder.stop",   1.0);
        if (echo > 3) printf("# %s logder.start=%g logder.step=%g logder.stop=%g %s\n",
            label, logder_energy_range[0]*eV, logder_energy_range[1]*eV, logder_energy_range[2]*eV,_eV);


        {   // construct initial smooth spherical potential
            set(potential[SMT].data(), rg[SMT]->n, potential[TRU].data() + nr_diff); // copy the tail of the spherical part of r*V_tru(r)
            if (echo > 2) printf("\n# %s construct initial smooth spherical potential as parabola\n", label);
            pseudize_function(potential[SMT].data(), rg[SMT], ir_cut[SMT], 2, 1); // replace by a parabola
        }

        int const maxit_scf = control::get("single_atom.init.scf.maxit", 1.);
        float const mixing = 0.5f;

        for(int scf = 0; scf < maxit_scf; ++scf) {
            if (echo > 1) printf("\n\n# %s SCF-iteration %d\n\n", label, scf);
//             update((scf >= maxit_scf - 333)*9); // switch full echo on in the last 3 iterations
            update_density(mixing, echo);
            double* const ves_multipoles = nullptr;
            update_potential(mixing, ves_multipoles, echo);
        } // self-consistency iterations

        // show the smooth and true potential
        if (false && (echo > 0)) {
            printf("\n## %s spherical parts: r, "
            "Zeff_tru(r), Zeff_smt(r)"
            ", r^2*rho_tru(r), r^2*rho_smt(r)"
            ", zero_potential(r) in Ha"
            ":\n", label);
            for(int ir = 0; ir < rg[SMT]->n; ir += 1) {
                auto const r = rg[SMT]->r[ir];
                printf("%g %g %g %g %g %g\n", r
//                         , -full_potential[TRU][00][ir + nr_diff]*Y00*r // for comparison, should be the same as Z_eff(r)
                        , -potential[TRU][ir + nr_diff] // Z_eff(r)
                        , -potential[SMT][ir] //    \tilde Z_eff(r)
                        , core_density[TRU][ir + nr_diff]*r*r*Y00*Y00
                        , core_density[SMT][ir]*r*r*Y00*Y00
                        , zero_potential[ir]*Y00 // not converted to eV
                      );
            } // ir
            printf("\n\n");
        } // echo

        // show the smooth and true partial waves
        if (false && (echo > 0)) {
            printf("\n\n\n# %s valence partial waves:\n", label);
            for(int ir = 0; ir < rg[SMT]->n; ir += 1) {
                auto const r = rg[SMT]->r[ir];
                auto const f = r; // scale with r? scale with Y00?
                printf("%g ", r);
                for(int ics = 0; ics < ncorestates; ++ics) {
                    printf("   %g", core_state[ics].wave[TRU][ir + nr_diff]*f);
                } // ics
                for(int iln = 0; iln < nln; ++iln) {
                    printf("   %g %g", partial_wave[iln].wave[TRU][ir + nr_diff]*f
                                     , partial_wave[iln].wave[SMT][ir]*f);
                } // iln
                printf("\n");
            } // ir
            printf("\n\n\n\n");
        } // echo

    } // constructor

    ~LiveAtom() { // destructor
        radial_grid::destroy_radial_grid(rg[TRU]);
        radial_grid::destroy_radial_grid(rg[SMT]);
    } // destructor

    status_t initialize_Gaunt() {
      if (gaunt_init) return 0; // success
      auto const stat = angular_grid::create_numerical_Gaunt<6>(gaunt);
      gaunt_init = (0 == stat);
      return stat;
    } // initialize_Gaunt

    void show_state_analysis(int const echo, radial_grid_t const *rg, double const wave[],
            int const enn, int const ell, float const occ, double const energy, int8_t const c0s1v2, int const ir_cut=-1) const {
        if (echo < 1) return; // this function only prints to the logg
        
        double stats[] = {0,0,0,0,0};
        for(int ir = 0; ir < rg->n; ++ir) {
            double const rho_wf = pow2(wave[ir]);
            double const dV = rg->r2dr[ir], r = rg->r[ir], r_inv = rg->rinv[ir];
            stats[0] += dV;
            stats[1] += rho_wf*dV;
            stats[2] += rho_wf*r*dV;
            stats[3] += rho_wf*r*r*dV;
            stats[4] += rho_wf*r_inv*dV; // Coulomb integral without -Z
        } // ir

        double charge_outside{0};
        if (ir_cut >= 0) {
            for(int ir = ir_cut; ir < rg->n; ++ir) {
                double const rho_wf = pow2(wave[ir]);
                double const dV = rg->r2dr[ir];
                charge_outside += rho_wf*dV;
            } // ir
        } // ir_cut

        printf("# %s %s %2d%c%6.1f E=%16.6f %s  <r>=%g rms=%g %s <r^-1>=%g %s q_out=%.3g e\n", label,
               c0s1v2_name[c0s1v2], enn, ellchar[ell], occ, energy*eV,_eV, 
               stats[2]/stats[1]*Ang, std::sqrt(std::max(0., stats[3]/stats[1]))*Ang,_Ang, 
               stats[4]/stats[1]*eV,_eV, charge_outside/stats[1]);
    } // show_state_analysis

    status_t pseudize_spherical_density(double smooth_density[], double const true_density[]
                                    , char const *quantity="core", int const echo=0) const {
        int const nrs = rg[SMT]->n;
        set(smooth_density, nrs, true_density + nr_diff); // copy the tail of the true density into the smooth density

        auto const stat = pseudize_function(smooth_density, rg[SMT], ir_cut[SMT], 3); // 3: use r^0, r^2 and r^4
        // alternatively, pseudize_function(smooth_density, rg[SMT], ir_cut[SMT], 3, 2); // 3, 2: use r^2, r^4 and r^6
        if (stat && (echo > 0)) printf("# %s Matching procedure for the smooth core density failed! info = %d\n", label, stat);

        if (echo > 8) { // plot the densities
            printf("\n## %s radius, {smooth, true} for %s density:\n", label, quantity);
            for(int ir = 0; ir < nrs; ir += 2) { // plot only every second entry
                printf("%g %g %g\n", rg[SMT]->r[ir], smooth_density[ir], true_density[ir + nr_diff]);
            }   printf("\n\n");
        } // plot
        return stat;
    } // pseudize_spherical_density
    
    
    void update_core_states(float const mixing, int const echo=0) {
        if (echo > 1) printf("\n# %s %s Z=%g\n", label, __func__, Z_core);
        // core states are feeling the spherical part of the hamiltonian only
        int const nr = rg[TRU]->n;
        std::vector<double> r2rho(nr);
        std::vector<double> new_r2core_density(nr, 0.0);
        double nelectrons{0};
        for(int ics = 0; ics < ncorestates; ++ics) {
            auto & cs = core_state[ics]; // abbreviate
            int constexpr SRA = 1;
            radial_eigensolver::shooting_method(SRA, *rg[TRU], potential[TRU].data(), cs.enn, cs.ell, cs.energy, cs.wave[TRU], r2rho.data());
            auto const norm = dot_product(nr, r2rho.data(), rg[TRU]->dr);
            auto const norm_factor = (norm > 0)? 1./std::sqrt(norm) : 0;
            auto const scal = pow2(norm_factor)*cs.occupation; // scaling factor for the density contribution of this state
            nelectrons += std::max(0.f, cs.occupation);
            // transform r*wave(r) as produced by the radial_eigensolver to wave(r)
            // and normalize the core level wave function to one
            scale(cs.wave[TRU], nr, rg[TRU]->rinv, norm_factor);
            // create wKin for the computation of the kinetic energy density
            product(cs.wKin[TRU], nr, potential[TRU].data(), cs.wave[TRU], -1.); // start as wKin = -r*V(r)*wave(r)
            add_product(cs.wKin[TRU], nr, rg[TRU]->r, cs.wave[TRU], cs.energy); // now wKin = r*(E - V(r))*wave

            if (scal > 0) add_product(new_r2core_density.data(), nr, r2rho.data(), scal);
            show_state_analysis(echo, rg[TRU], cs.wave[TRU], cs.enn, cs.ell, cs.occupation, cs.energy, 2*(scal < 0), ir_cut[TRU]);
        } // ics

        // report integrals
        auto const old_core_charge = dot_product(nr, rg[TRU]->r2dr, core_density[TRU].data());
        auto const new_core_charge = dot_product(nr, rg[TRU]->dr, new_r2core_density.data());
        if (echo > 0) printf("# %s expect a core density with %g electrons\n", label, nelectrons);
        if (echo > 0) printf("# %s previous core density has %g electrons\n", label, old_core_charge);
        if (echo > 0) printf("# %s new core density has %g electrons\n",      label, new_core_charge);
        double mix_new = mixing, mix_old = 1 - mixing;
        // can we rescale the mixing coefficients such that the desired number of core electrons comes out?
        auto const mixed_charge = mix_old*old_core_charge + mix_new*new_core_charge;
        if (mixed_charge != 0) {
            auto const rescale = nelectrons/mixed_charge;
            mix_old *= rescale;
            mix_new *= rescale;
        } // rescale

        double core_density_change{0}, core_density_change2{0}, core_nuclear_energy{0};
        for(int ir = 0; ir < nr; ++ir) {
            auto const new_rho = new_r2core_density[ir]*pow2(rg[TRU]->rinv[ir]); // *r^{-2}
            core_density_change  += std::abs(new_rho - core_density[TRU][ir])*rg[TRU]->r2dr[ir];
            core_density_change2 +=     pow2(new_rho - core_density[TRU][ir])*rg[TRU]->r2dr[ir];
            core_nuclear_energy  +=         (new_rho - core_density[TRU][ir])*rg[TRU]->rdr[ir]; // Coulomb integral change
            core_density[TRU][ir] = mix_new*new_rho + mix_old*core_density[TRU][ir];
        } // ir
        core_nuclear_energy *= -Z_core;
        if (echo > 0) printf("# %s core density change %g e (rms %g e) energy change %g %s\n", label,
            core_density_change, std::sqrt(std::max(0.0, core_density_change2)), core_nuclear_energy*eV,_eV);

        pseudize_spherical_density(core_density[SMT].data(), core_density[TRU].data(), "core", echo - 1); 

        // report integrals
        auto const tru_core_charge = dot_product(rg[TRU]->n, rg[TRU]->r2dr, core_density[TRU].data());
        auto const smt_core_charge = dot_product(rg[SMT]->n, rg[SMT]->r2dr, core_density[SMT].data());
        core_charge_deficit = tru_core_charge - smt_core_charge;
        if (echo > 1) printf("# %s true and smooth core density have %g and %g electrons\n", label, tru_core_charge, smt_core_charge);

    } // update_core_states

    void update_partial_waves(int const echo=0) {
        if (echo > 1) printf("\n# %s %s Z=%g\n", label, __func__, Z_core);
        // the basis for valence partial waves is generated from the spherical part of the hamiltonian
//      auto const small_component = new double[rg[TRU]->n];
        int const nr = rg[TRU]->n;
        std::vector<double> r2rho(nr);

        // ToDo: projectors only depend on sigma, numax and the radial grid --> move this to the contructor
        scattering_test::expand_sho_projectors(projectors.data(), projectors.stride(), *rg[SMT], sigma, numax, 1, echo/2);


        int const ir_match[] = {radial_grid::find_grid_index(*rg[TRU], r_match),
                                radial_grid::find_grid_index(*rg[SMT], r_match)};
        if (echo > 3) printf("# %s matching radius %g %s at radial indices %i and %i\n",
                               label, r_match*Ang, _Ang, ir_match[TRU], ir_match[SMT]);

        if (echo > 9) { printf("\n## %s show the local potentials (r, r*Vtru, r*Vsmt):\n", label);
            for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                printf("%g %g %g\n", rg[SMT]->r[ir], potential[TRU][ir + nr_diff], potential[SMT][ir]);
            }   printf("\n\n");
        } // echo

        for(int ell = 0; ell <= numax; ++ell) {
            int const ln_off = sho_tools::ln_index(numax, ell, 0); // offset where to start indexing valence states

            if (echo > 3) printf("\n# %s %s for ell=%i\n\n", label, __func__, ell);

            view2D<double> projectors_ell(projectors[ln_off], projectors.stride()); // sub-view

            int const n = nn[ell];
            for(int nrn = 0; nrn < n; ++nrn) { // smooth number or radial nodes
                int const iln = ln_off + nrn;
                auto & vs = partial_wave[iln]; // abbreviate
                int constexpr SRA = 1;

                // solve for a true partial wave
    //          radial_integrator::integrate_outwards<SRA>(*rg[TRU], potential[TRU], ell, vs.energy, vs.wave[TRU], small_component);
                set(vs.wave[TRU], nr, 0.0); // clear

                if (true) {
                    // solve for a true valence eigenstate
                    radial_eigensolver::shooting_method(SRA, *rg[TRU], potential[TRU].data(), vs.enn, ell, vs.energy, vs.wave[TRU], r2rho.data());
                    if (echo > 7) printf("# %s %s found a true %i%c-eigenstate of the spherical potential at E=%g %s\n",
                                          label, __func__, vs.enn, atom_core::ellchar(ell), vs.energy*eV,_eV);
                } else {
                    assert(nrn > 0);
                    vs.energy = partial_wave[iln - 1].energy + 1.0; // copy energy from lower state and add 1.0 Hartree
                    int nnodes = 0;
                    // integrate outwards
                    radial_integrator::shoot(SRA, *rg[TRU], potential[TRU].data(), ell, vs.energy, nnodes, vs.wave[TRU], r2rho.data());
                } // bound state?

                // normalize the partial waves
                auto const norm_wf2 = dot_product(nr, r2rho.data(), rg[TRU]->dr);
                auto const norm_factor = 1./std::sqrt(norm_wf2);
                scale(vs.wave[TRU], nr, rg[TRU]->rinv, norm_factor); // transform r*wave(r) as produced by the radial_eigensolver to wave(r)

                // create wKin for the computation of the kinetic energy density
                product(vs.wKin[TRU], nr, potential[TRU].data(), vs.wave[TRU], -1.); // start as wKin = -r*V(r)*wave(r)
                add_product(vs.wKin[TRU], nr, rg[TRU]->r, vs.wave[TRU], vs.energy); // now wKin = r*(E - V(r))*wave(r)

//                 auto const tru_norm = dot_product(ir_cut[TRU], r2rho.data(), rg[TRU]->dr)/norm_wf2; // integrate only up to rcut
//                 auto const work = r2rho.data();
//                 scale(work, nr, potential[TRU].data());
//                 auto const tru_Epot = dot_product(ir_cut[TRU], work, rg[TRU]->dr)/norm_wf2;
//                 auto const tru_kinetic_E = vs.energy*tru_norm - tru_Epot; // kinetic energy contribution up to r_cut

//                 if (echo > 1) printf("# valence %2d%c%6.1f E=%16.6f %s\n", vs.enn, ellchar[ell], vs.occupation, vs.energy*eV,_eV);
                show_state_analysis(echo, rg[TRU], vs.wave[TRU], vs.enn, ell, vs.occupation, vs.energy, 0, ir_cut[TRU]);


                // idea: make this module flexible enough so it can load a potential and
                //       generate PAW data in XML format (GPAW, ABINIT) using the SHO method

                { // scope: generate smooth partial waves from projectors, revPAW scheme
                    int const stride = align<2>(rg[SMT]->n);
                    view2D<double> rphi(1 + n, stride);
                    view2D<double> Tphi(1 + n, stride);
                    std::vector<double> ff(rg[SMT]->n);
                    auto const vgtru = vs.wave[TRU][ir_match[TRU]]    *rg[TRU]->r[ir_match[TRU]];
                    auto const dgtru = vs.wave[TRU][ir_match[TRU] - 1]*rg[TRU]->r[ir_match[TRU] - 1] - vgtru;
                    double vghom, dghom;

                    int constexpr HOM = 0;
                    for(int krn = HOM; krn <= n; ++krn) { // loop must run serial and forward
                        // krn == 0 generates the homogeneous solution in the first iteration
                        auto const projector = (krn > HOM) ? projectors_ell[krn - 1] : nullptr;
                        double dg; // derivative at end point, not used

                        // solve inhomgeneous equation and match true wave in value and derivative
                        radial_integrator::integrate_outwards<SRA>(*rg[SMT], potential[SMT].data(), ell, vs.energy,
                                                                  rphi[krn], ff.data(), -1, &dg, projector);

                        auto const vginh = rphi[krn][ir_match[SMT]];
                        auto const dginh = rphi[krn][ir_match[SMT] - 1] - vginh;
                        if (HOM == krn) {
                            vghom = vginh; dghom = dginh;
                        } else {
                            // matching coefficient - how much of the homogeneous solution do we need to add to match logder
                            auto const denom =    vgtru*dghom - dgtru*vghom; // ToDo: check if denom is not too small
                            auto const c_hom = - (vgtru*dginh - dgtru*vginh) / denom;

                            // combine inhomogeneous and homogeneous solution to match true solution outside r_match
                            add_product(rphi[krn], rg[SMT]->n, rphi[HOM], c_hom);

                            auto const scal = vgtru / (vginh + c_hom*vghom); // match not only in logder but also in value
                            scale(rphi[krn], rg[SMT]->n, rg[SMT]->rinv, scal); // and divide by r
                            // ToDo: extrapolate lowest point?

                            // Tphi = projector + (E - V)*phi
                            for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                                Tphi[krn][ir] = scal*projector[ir] + (vs.energy*rg[SMT]->r[ir] - potential[SMT][ir])*rphi[krn][ir];
                            } // ir

                            // now visually check that the matching of value and derivative of rphi is ok.
                            if (echo > 9) {
                                printf("\n## %s check matching of rphi for ell=%i nrn=%i krn=%i (r, phi_tru, phi_smt, rTphi_tru, rTphi_smt):\n",
                                        label, ell, nrn, krn-1);
                                for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                                    printf("%g  %g %g  %g %g\n", rg[SMT]->r[ir],
                                        vs.wave[TRU][ir + nr_diff], rphi[krn][ir],
                                        vs.wKin[TRU][ir + nr_diff], Tphi[krn][ir]);
                                } // ir
                                printf("\n\n");
                            } // echo

                            // check that the matching of value and derivative of rphi is ok by comparing value and derivative
                            if (echo > 9) {
                                printf("# %s check matching of vg and dg for ell=%i nrn=%i krn=%i: %g == %g ? and %g == %g ?\n",
                                       label, ell, nrn, krn-1,  vgtru, scal*(vginh + c_hom*vghom),  dgtru, scal*(dginh + c_hom*dghom));
                            } // echo

                        } // krn > 0
                    } // krn

                    std::vector<double> evec(n, 0.0);
                    evec[nrn] = 1.0;    // as suggested by Baumeister+Tsukamoto in PASC19 proceedings
                    if (true && (n > 1)) { // as suggested by Morian Sonnet: minimize the radial curvature of the smooth partial wave
                        view2D<double> Ekin(     3*n , n);
                        view2D<double> Olap(Ekin[1*n], n);
                        view2D<double> Vkin(Ekin[2*n], n);
                        int const nr_max = ir_match[SMT];
                        for(int krn = 0; krn < n; ++krn) {
                            for(int jrn = 0; jrn < n; ++jrn) {
                                // definition of Tphi contains factor *r
                                Ekin[krn][jrn] = dot_product(nr_max, Tphi[1 + krn], rphi[1 + jrn], rg[SMT]->rdr);
                                Olap[krn][jrn] = dot_product(nr_max, rphi[1 + krn], rphi[1 + jrn], rg[SMT]->r2dr);
                                Vkin[krn][jrn] = dot_product(nr_max, rphi[1 + krn], rphi[1 + jrn], rg[SMT]->dr)
                                                 * 0.5*(ell*(ell + 1)); // kinetic energy in angular direction, Hartree units
                                // minimize only the radial kinetic energy
                                Ekin[krn][jrn] -= Vkin[krn][jrn]; // subtract the angular part, suggested by Morian Sonnet
                            } // jrn
                        } // krn

                        // symmetrize
                        for(int krn = 0; krn < n; ++krn) {
                            for(int jrn = 0; jrn < krn; ++jrn) { // triangular loop
                                symmetrize(Ekin[krn][jrn], Ekin[jrn][krn]);
                                symmetrize(Olap[krn][jrn], Olap[jrn][krn]);
                            } // jrn
                        } // krn

                        if (echo > 7) {
//                          printf("\n");
                            for(int krn = 0; krn < n; ++krn) {
                                printf("# %s curvature|overlap for i=%i ", label, krn);
                                for(int jrn = 0; jrn < n; ++jrn) {
                                    printf(" %g|%g", Ekin[krn][jrn], Olap[krn][jrn]);
                                } // jrn
                                printf("\n");
                            } // krn
//                          printf("\n");
                        } // echo

                        { // scope: minimize the radial curvature of the smooth partial wave
                            double lowest_eigenvalue = 0;
                            auto const info = minimize_curvature(n, Ekin, Olap, &lowest_eigenvalue);
                            if (info) {
                                printf("# %s generalized eigenvalue problem failed info=%i\n", label, info);
                            } else {
                                set(evec.data(), n, Ekin[0]);
                                if (echo > 6) {
                                    printf("# %s lowest eigenvalue %g %s", label, lowest_eigenvalue*eV, _eV);
                                    if (echo > 8) {
                                        printf(", coefficients"); for(int krn = 0; krn < n; ++krn) printf(" %g", Ekin[0][krn]);
                                    } // high echo
                                    printf("\n");
                                } // echo
                            } // info
                        } // scope

                    } // minimize radial kinetic energy of partial wave

                    set(vs.wave[SMT], rg[SMT]->n, 0.0);
                    set(vs.wKin[SMT], rg[SMT]->n, 0.0);
                    double sum{0}; for(int krn = 0; krn < n; ++krn) sum += evec[krn];
                    double const scal = 1.0/sum; // scaling such that the sum is 1
                    for(int krn = 0; krn < n; ++krn) {
                        auto const coeff = scal*evec[krn];
                        add_product(vs.wave[SMT], rg[SMT]->n, rphi[1 + krn], coeff);
                        add_product(vs.wKin[SMT], rg[SMT]->n, Tphi[1 + krn], coeff);
                    } // krn

                    // ToDo: now check that the matching of value and derivative of wave and wKin are ok.
                    if (echo > 9) {
                        printf("\n## %s check matching of partial waves ell=%i nrn=%i (r, phi_tru, phi_smt, rTphi_tru, rTphi_smt):\n",
                                label, ell, nrn);
                        for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                            printf("%g  %g %g  %g %g\n", rg[SMT]->r[ir],
                                vs.wave[TRU][ir + nr_diff], vs.wave[SMT][ir],
                                vs.wKin[TRU][ir + nr_diff], vs.wKin[SMT][ir]);
                        } // ir
                        printf("\n\n");
                    } // echo

//                     if (echo > 8) {
//                         printf("# fflush(stdout) in line %i\n", __LINE__ + 1);
//                         fflush(stdout);
//                     } // echo

                } // scope


            } // nrn



            { // scope: establish dual orthgonality with [SHO] projectors
                int const nr = rg[SMT]->n; //, mr = align<2>(nr);
                // int const stride = align<2>(rg[SMT]->n);
                int const msub = (1 + numax/2); // max. size of the subspace

                if (echo > 4) { // show normalization and orthogonality of projectors
                    for(int nrn = 0; nrn < n; ++nrn) {
                        for(int krn = 0; krn < n; ++krn) {
                            printf("# %s %c-projector <#%d|#%d> = %i + %g  sigma=%g %s\n", label, ellchar[ell], nrn, krn,
                                (nrn == krn), dot_product(nr, projectors_ell[nrn], projectors_ell[krn], rg[SMT]->dr) - (nrn==krn), sigma*Ang,_Ang);
                        } // krn
                    } // nrn
                } // echo

                view2D<double> ovl(msub, msub); // get memory
                for(int nrn = 0; nrn < n; ++nrn) { // smooth number or radial nodes
                    int const iln = ln_off + nrn;
                    auto const wave = partial_wave[iln].wave[SMT];
                    for(int krn = 0; krn < n; ++krn) { // smooth number or radial nodes
                        ovl[nrn][krn] = dot_product(nr, wave, projectors_ell[krn], rg[SMT]->rdr);
                        if (echo > 2) printf("# %s smooth partial %c-wave #%d with %c-projector #%d has overlap %g\n",
                                               label, ellchar[ell], nrn, ellchar[ell], krn, ovl[nrn][krn]);
                    } // krn
                } // nrn

                view2D<double> inv(msub, msub); // get memory
                if (n > 4) exit(__LINE__); // error not implemented
                double const det = simple_math::invert(n, inv.data(), inv.stride(), ovl.data(), ovl.stride());
                if (echo > 2) printf("# %s determinant for %c-projectors %g\n", label, ellchar[ell], det);

                // make a new linear combination
                for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                    int const nrts = rg[ts]->n, mrts = align<2>(nrts);
                    view3D<double> waves(2, n, mrts, 0.0); // temporary storage for pairs {wave, wKin}
                    for(int nrn = 0; nrn < n; ++nrn) {
                        int const iln = ln_off + nrn;
                        set(waves[0][nrn], nrts, partial_wave[iln].wave[ts]); // copy
                        set(waves[1][nrn], nrts, partial_wave[iln].wKin[ts]); // copy
                    } // nrn
                    for(int nrn = 0; nrn < n; ++nrn) {
                        int const iln = ln_off + nrn;
                        set(partial_wave[iln].wave[ts], nrts, 0.0); // clear
                        set(partial_wave[iln].wKin[ts], nrts, 0.0); // clear
                        for(int krn = 0; krn < n; ++krn) {
                            add_product(partial_wave[iln].wave[ts], nrts, waves[0][krn], inv[nrn][krn]);
                            add_product(partial_wave[iln].wKin[ts], nrts, waves[1][krn], inv[nrn][krn]);
                        } // krn
                    } // nrn
                } // ts - tru and smt

                // check again the overlap, seems ok
                for(int nrn = 0; nrn < n; ++nrn) { // smooth number or radial nodes
                    int const iln = ln_off + nrn;
                    auto const wave = partial_wave[iln].wave[SMT];
                    for(int krn = 0; krn < n; ++krn) { // smooth number or radial nodes
                        ovl[nrn][krn] = dot_product(nr, wave, projectors_ell[krn], rg[SMT]->rdr);
                        if (echo > 2) printf("# %s smooth partial %c-wave #%d with %c-projector #%d new overlap %g\n",
                                               label, ellchar[ell], nrn, ellchar[ell], krn, ovl[nrn][krn]);
                    } // krn
                } // nrn

                // compute kinetic energy difference matrix from wKin
                for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                    int const nr_cut = rg[ts]->n; // integration over the entire grid --> diagonal elements then appear positive.
                    for(int iln = 0 + ln_off; iln < n + ln_off; ++iln) {
                        for(int jln = 0 + ln_off; jln < n + ln_off; ++jln) {
                            kinetic_energy(ts,iln,jln) = dot_product(nr_cut,
                            	partial_wave[iln].wKin[ts],
                              partial_wave[jln].wave[ts], rg[ts]->rdr); // we only need rdr here since wKin is defined as r*(E - V(r))*wave(r)
                        } // j
                    } // i
                } // ts

                // display
                for(int i = 0; i < n; ++i) {
                    for(int j = 0; j < n; ++j) {
                        auto const E_kin_tru = kinetic_energy(TRU,i+ln_off,j+ln_off);
                        auto const E_kin_smt = kinetic_energy(SMT,i+ln_off,j+ln_off);
                        if (echo > 0) printf("# %s %c-channel <%d|T|%d> kinetic energy [unsymmetrized] (true) %g and (smooth) %g (diff) %g %s\n",
                          label, ellchar[ell], i, j, E_kin_tru*eV, E_kin_smt*eV, (E_kin_tru - E_kin_smt)*eV, _eV);
                    } // j
                } // i

                if (1) { // symmetrize the kinetic energy tensor
                    for(int iln = 0 + ln_off; iln < n + ln_off; ++iln) {
                        for(int jln = 0 + ln_off; jln < iln; ++jln) { // triangular loop excluding the diagonal elements
                            if (1) { // usual
                                for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                                    auto &aij = kinetic_energy(ts,iln,jln);
                                    auto &aji = kinetic_energy(ts,jln,iln);
                                    auto const avg = 0.5*(aij + aji);
                                    aij = avg;
                                    aji = avg;
                                } // ts
                            } else {
                                // actually we do not need to symmetrize each contribution kinetic_energy[ts]
                                // but it would be enough to symmetrize the difference matrix kinetic_energy[TRU] - kinetic_energy[SMT]
                                // --> symmetrize only the difference
                                for(int twice = 0; twice < 2; ++twice) { // do this twice for numerical accuracy
                                    auto const dij = kinetic_energy(TRU,iln,jln) - kinetic_energy(SMT,iln,jln);
                                    auto const dji = kinetic_energy(TRU,jln,iln) - kinetic_energy(SMT,jln,iln);
                                    auto const asy = 0.5*(dij - dji); // asymmetry
                                    for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
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

                // display
                for(int i = 0; i < n; ++i) {
                    for(int j = 0; j < n; ++j) {
                        auto const E_kin_tru = kinetic_energy(TRU,i+ln_off,j+ln_off);
                        auto const E_kin_smt = kinetic_energy(SMT,i+ln_off,j+ln_off);
                        if (echo > 0) printf("# %s %c-channel <%d|T|%d> kinetic energy [symmetrized] (true) %g and (smooth) %g (diff) %g %s\n",
                          label, ellchar[ell], i, j, E_kin_tru*eV, E_kin_smt*eV, (E_kin_tru - E_kin_smt)*eV, _eV);
                    } // j
                } // i

            } // scope: establish dual orthgonality with [SHO] projectors

        } // ell

    } // update_partial_waves

    void update_charge_deficit(int const echo=0) {
        int const nln = sho_radial::nSHO_radial(numax);
        { // scope: generate a vector true_norm such that
          // true_norm[iln]*charge_deficit[0][TRU][iln][jln]*true_norm[jln] is normalized to 1
            int const ts = TRU;
            int const nr = rg[ts]->n; // entire radial grid
            for(int iln = 0; iln < nln; ++iln) {
                if (0) {
                    auto const wave_i = partial_wave[iln].wave[ts];
                    auto const norm2 = dot_product(nr, wave_i, wave_i, rg[ts]->r2dr);
                    true_norm[iln] = (norm2 > 1e-99) ? 1./std::sqrt(norm2) : 0;
                } else true_norm[iln] = 1.;
            } // iln
        } // scope

        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            int const nr = rg[ts]->n; // integrate over the full radial grid
            std::vector<double> rl(nr, 1.0); // init as r^0
            std::vector<double> wave_r2rl_dr(nr);
            if (echo > 4) printf("\n# %s charges for %s partial waves\n", label, (TRU==ts)?"true":"smooth");
            for(int ell = 0; ell <= ellmax_compensator; ++ell) { // loop-carried dependency on rl, run forward, run serial!
                bool const echo_l = (echo > 4 + 2*(ell > 0));
                if (echo_l) printf("# %s charges for ell=%d, jln = 0, 1, ...\n", label, ell);
                if (ell > 0) scale(rl.data(), nr, rg[ts]->r); // create r^{\ell}
                for(int iln = 0; iln < nln; ++iln) {
                    if (echo_l) printf("# %s iln = %d ", label, iln);
                    auto const wave_i = partial_wave[iln].wave[ts];
                    product(wave_r2rl_dr.data(), nr, wave_i, rl.data(), rg[ts]->r2dr); // product of three arrays
                    for(int jln = 0; jln < nln; ++jln) {
                        auto const wave_j = partial_wave[jln].wave[ts];
                        auto const cd = dot_product(nr, wave_r2rl_dr.data(), wave_j);
                        charge_deficit(ell,ts,iln,jln) = cd;
                        if (echo_l) printf("\t%10.6f", true_norm[iln]*cd*true_norm[jln]);
//                      if ((SMT==ts) && (echo > 1)) printf("\t%10.6f", charge_deficit(ell,TRU,iln,jln) - cd);
                    } // jln
                    if (echo_l) printf("\n");
                } // iln
                if (echo_l) printf("\n");
            } // ell
        } // ts
    } // update_charge_deficit

    template<typename int_t>
    void get_valence_mapping(int_t ln_index_list[], int_t lm_index_list[],
                             int_t lmn_begin[], int_t lmn_end[],
                             int const echo=0) const {
        int const mlm = pow2(1 + numax);
        set(lmn_begin, mlm, int_t(-1));
        for(int ell = 0; ell <= numax; ++ell) {
            for(int emm = -ell; emm <= ell; ++emm) {
                for(int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
                    int const ilmn      = sho_tools::lmn_index(numax, ell, emm, nrn);
                    ln_index_list[ilmn] = sho_tools::ln_index(numax, ell, nrn); // valence state index
                    int const lm        = sho_tools::lm_index(ell, emm);
                    lm_index_list[ilmn] = lm;
                    if (lmn_begin[lm] < 0) lmn_begin[lm] = ilmn; // store the first index of this lm
                    lmn_end[lm] = ilmn + 1; // store the last index of this lm
                } // nrn
            } // emm
        } // ell

        if (echo > 3) {
            printf("# %s ln_index_list ", label);
            for(int i = 0; i < sho_tools::nSHO(numax); ++i) {
                printf(" %i", ln_index_list[i]);
            }   printf("\n");
            printf("# %s lmn_begin-lmn_end ", label);
            for(int i = 0; i < mlm; ++i) {
                if (lmn_begin[i] == lmn_end[i] - 1) {
                    printf(" %i", lmn_begin[i]);
                } else {
                    printf(" %i-%i", lmn_begin[i], lmn_end[i] - 1);
                }
            }   printf("\n");
        } // echo

    } // get_valence_mapping


    void transform_SHO(double out[], int const out_stride,
                  double const in[], int const in_stride,
                  bool const in_Cartesian, double const alpha=1) const {

        int const N = unitary_Ezyx_lmn.stride();
        view2D<double const> uni(unitary_Ezyx_lmn.data(), N);
        view2D<double> tmp(N, N); // get memory
        view2D<double const> inp(in, in_stride);
        view2D<double> res(out, out_stride);

//         double const beta = 0;
//         char const nn = 'n', tn = in_Cartesian?'n':'t', nt = in_Cartesian?'t':'n';

        if (in_Cartesian) {
            // transform Cartesian input to Radial output

//             tmp(=C)[n_C][m_R] = in(=B)[n_C][k_C] * unitary(=A)[k_C][m_R]    // step 1
//             out(=C)[n_R][m_R] = unitary(=B^T)[n_R][k_C] * tmp(=A)[k_C][m_R] // step 2

            for(int nC = 0; nC < N; ++nC) {
                for(int mR = 0; mR < N; ++mR) {
                    double tij = 0;
                    for(int kC = 0; kC < N; ++kC) {
                        tij += inp[nC][kC] * uni[kC][mR]; // *u
                    } // kC
                    tmp[nC][mR] = alpha*tij;
                } // mR
            } // nC
            for(int nR = 0; nR < N; ++nR) {
                for(int mR = 0; mR < N; ++mR) {
                    double tij = 0;
                    for(int kC = 0; kC < N; ++kC) {
                        tij += uni[kC][nR] * tmp[kC][mR]; // u^T*
                    } // kC
                    res[nR][mR] = alpha*tij;
                } // mR
            } // nR

        } else { // in_Cartesian
            // transform Radial input to Cartesian output

//             tmp(=C)[n_R][m_C] = in(=B)[n_R][k_R] * unitary(=A^T)[k_R][m_C] // step 1
//             out(=C)[n_C][m_C] = unitary(=B)[n_C][k_R] * tmp(=A)[k_R][m_C]  // step 2

            for(int nC = 0; nC < N; ++nC) {
                for(int mR = 0; mR < N; ++mR) {
                    double tij = 0;
                    for(int kR = 0; kR < N; ++kR) {
                        tij += uni[nC][kR] * inp[kR][mR]; // u*
                    } // kR
                    tmp[nC][mR] = alpha*tij;
                } // mR
            } // nC
            for(int nC = 0; nC < N; ++nC) {
                for(int mC = 0; mC < N; ++mC) {
                    double tij = 0;
                    for(int kR = 0; kR < N; ++kR) {
                        tij += tmp[nC][kR] * uni[mC][kR]; // *u^T
                    } // kR
                    res[nC][mC] = alpha*tij;
                } // mC
            } // nC

        } // in_Cartesian

// now using BLAS            ToDo!
//         void dgemm_(const char* tA, const char* tB, const int* M, const int* N, const int* K, const double* alpha,
//               const double* A, const int* sA, const double* B, const int* sB, const double* beta, double* C, const int* sC);
//         performs C[n][m] := beta*C[n][m] + alpha*sum_k B[n][k] * A[k][m];

//         dgemm_(&tn, &nn, &N, &N, &N, &alpha, uni.data(), &N, in, &in_stride, &beta, tmp, &N);
//         dgemm_(&nn, &nt, &N, &N, &N, &alpha, tmp, &N, uni.data(), &N, &beta, out, &out_stride);

    } // transform_SHO

    void get_rho_tensor(view3D<double> & density_tensor
        , view2D<double> const & density_matrix, sho_tools::SHO_order_t const order=sho_tools::order_Ezyx
        , int const echo=0) {
        int const nSHO = sho_tools::nSHO(numax);
        int const stride = nSHO;
        assert(stride >= nSHO);

        initialize_Gaunt();

        int const lmax = std::max(ellmax, ellmax_compensator);
        int const nlm = pow2(1 + lmax);
        int const mlm = pow2(1 + numax);

        view2D<double> radial_density_matrix;
        
        if (sho_tools::order_Ezyx == order) {
            // now transform the density_matrix[izyx][jzyx]
            //     into a radial_density_matrix[ilmn][jlmn]
            //     using the unitary transform from left and right
            radial_density_matrix = view2D<double>(nSHO, stride);
            transform_SHO(radial_density_matrix.data(), stride,
                  density_matrix.data(), density_matrix.stride(), true);
            
            if (0) { // debugging
                view2D<double> check_matrix(nSHO, nSHO);
                transform_SHO(check_matrix.data(), check_matrix.stride(),
                              radial_density_matrix.data(), radial_density_matrix.stride(), false);
                double maxdev{0};
                for(int i = 0; i < nSHO; ++i) {
                    for(int j = 0; j < nSHO; ++j) {
                          maxdev = std::max(maxdev, std::abs(check_matrix[i][j] - density_matrix[i][j]));
                    } // j
                } // i
                printf("# %s found max deviation %.1e when backtransforming the density matrix\n\n", label, maxdev);
                assert(maxdev < 1e-9);
            } // debugging
            
        } else
        if (sho_tools::order_lmn == order) {
            radial_density_matrix = view2D<double>(density_matrix.data(), density_matrix.stride()); // wrap
        }

        if (0) {
            printf("# %s Radial density matrix\n", label);
            for(int i = 0; i < nSHO; ++i) {
                for(int j = 0; j < nSHO; ++j) {
                    printf("\t%.1f", radial_density_matrix[i][j]);
                }   printf("\n");
            }   printf("\n");
        } // plot

        //   Then, contract with the Gaunt tensor over m_1 and m_2
        //   rho_tensor[lm][iln][jln] =
        //     G_{lm l_1m_1 l_2m_2} * radial_density_matrix[il_1m_1n_1][jl_2m_2n_2]

        set(density_tensor, nlm, 0.0); // clear
        for(auto gnt : gaunt) {
            int const lm = gnt.lm, lm1 = gnt.lm1, lm2 = gnt.lm2; auto G = gnt.G;
            // if (0 == lm) assert(std::abs(G - Y00*(lm1 == lm2)) < 1e-12); // make sure that G_00ij = delta_ij*Y00
            if (0 == lm) G = Y00*(lm1 == lm2); // make sure that G_00ij = delta_ij*Y00
            if ((lm < nlm) && (lm1 < mlm) && (lm2 < mlm)) {
                for(int ilmn = lmn_begin[lm1]; ilmn < lmn_end[lm1]; ++ilmn) {
                    int const iln = ln_index_list[ilmn];
                    for(int jlmn = lmn_begin[lm2]; jlmn < lmn_end[lm2]; ++jlmn) {
                        int const jln = ln_index_list[jlmn];
                        density_tensor(lm,iln,jln) += G * radial_density_matrix[ilmn][jlmn];

//                         auto const rho_ij = rho_tensor[lm][iln][jln];
//                         if (std::abs(rho_ij) > 1e-9)
//                             printf("# LINE=%d rho_ij = %g for lm=%d iln=%d jln=%d\n", __LINE__, rho_ij*Y00inv, lm, iln, jln);
                    } // jlmn
                } // ilmn
            } // limits
        } // gnt

#ifdef FULL_DEBUG
        for(int lm = 0; lm < nlm; ++lm) {
            for(int iln = 0; iln < nln; ++iln) {
                for(int jln = 0; jln < nln; ++jln) {
                    auto const rho_ij = density_tensor(lm,iln,jln);
                    if (std::abs(rho_ij) > 1e-9)
                        printf("# %s LINE=%d rho_ij = %g for lm=%d iln=%d jln=%d\n", label, __LINE__, rho_ij*Y00inv, lm, iln, jln);
                } // jln
            } // iln
        } // lm
#endif
    } // get_rho_tensor


    void update_full_density(view3D<double> const & density_tensor, int const echo=0) { // density tensor rho_{lm iln jln}
        int const nlm = pow2(1 + ellmax);
        int const nln = sho_radial::nSHO_radial(numax);
        // view3D<double const> density_tensor(rho_tensor, nln, nln); // rho_tensor[mln][nln][nln]

        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            size_t const nr = rg[ts]->n;
            assert(full_density[ts].stride() >= nr);
            set(full_density[ts], nlm, 0.0); // clear
            for(int lm = 0; lm < nlm; ++lm) {
                if (00 == lm) {
                    set(full_density[ts][lm], nr, core_density[ts].data(), Y00); // needs scaling with Y00 since core_density has a factor 4*pi
                    if (echo > 0) printf("# %s %s density has %g electrons after adding the core density\n", label,
                     (TRU == ts)?"true":"smooth", dot_product(nr, full_density[ts][lm], rg[ts]->r2dr)*Y00inv);
                } // 00 == lm
                for(int iln = 0; iln < nln; ++iln) {
                    auto const wave_i = partial_wave[iln].wave[ts];
                    assert(nullptr != wave_i);
                    for(int jln = 0; jln < nln; ++jln) {
                        auto const wave_j = partial_wave[jln].wave[ts];
                        assert(nullptr != wave_j);
                        double const rho_ij = density_tensor(lm,iln,jln);
//                         if (std::abs(rho_ij) > 1e-9)
//                             printf("# %s rho_ij = %g for lm=%d iln=%d jln=%d\n", label, rho_ij*Y00inv, lm, iln, jln);
                        add_product(full_density[ts][lm], nr, wave_i, wave_j, rho_ij);
                    } // jln
                } // iln
            } // lm
            if (echo > 0) printf("# %s %s density has %g electrons after adding the valence density\n",  label,
                    (TRU == ts)?"true":"smooth", dot_product(nr, full_density[ts][00], rg[ts]->r2dr)*Y00inv);

        } // true and smooth

        int const nlm_cmp = pow2(1 + ellmax_compensator);
        for(int ell = 0; ell <= ellmax_compensator; ++ell) {
            for(int emm = -ell; emm <= ell; ++emm) {
                int const lm = solid_harmonics::lm_index(ell, emm);
                double rho_lm = 0;
                for(int iln = 0; iln < nln; ++iln) {
                    for(int jln = 0; jln < nln; ++jln) {
                        double const rho_ij = density_tensor(lm,iln,jln);
                        if (std::abs(rho_ij) > 1e-9)
                            printf("# %s rho_ij = %g for ell=%d emm=%d iln=%d jln=%d\n", label, rho_ij*Y00inv, ell, emm, iln, jln);
                        rho_lm += rho_ij * ( charge_deficit(ell,TRU,iln,jln)
                                           - charge_deficit(ell,SMT,iln,jln) );
                    } // jln
                } // iln
                assert(lm >= 0);
                assert(lm < nlm_cmp);
                qlm_compensator[lm] = rho_lm;
            } // emm
        } // ell

        // account for Z protons in the nucleus and the missing charge in the smooth core density
        qlm_compensator[0] += Y00*(core_charge_deficit - Z_core);
        if (echo > 5) printf("# %s compensator monopole charge is %g electrons\n", label, qlm_compensator[0]*Y00inv);

        { // scope: construct the augmented density
            int const nlm_aug = pow2(1 + std::max(ellmax, ellmax_compensator));
            auto const mr = full_density[SMT].stride(); // on the smooth grid
            assert(aug_density.stride() == mr);
            set(aug_density, nlm_aug, 0.0); // clear all entries
            set(aug_density.data(), nlm*mr, full_density[SMT].data()); // copy smooth full_density, need spin summation?
            add_or_project_compensators<0>(aug_density, qlm_compensator.data(), ellmax_compensator, rg[SMT], sigma_compensator);
            double const aug_charge = dot_product(rg[SMT]->n, rg[SMT]->r2dr, aug_density[00]); // only aug_density[00==lm]
            if (echo > 5) printf("# %s augmented density shows an ionization of %g electrons\n", label, aug_charge*Y00inv); // this value should be small

            double const tru_charge = dot_product(rg[TRU]->n, rg[TRU]->r2dr, full_density[TRU][00]); // only full_density[0==lm]
            if (echo > 5) printf("# %s true density has %g electrons\n", label, tru_charge*Y00inv); // this value should be of the order of Z
        } // scope

    } // update_full_density


    void update_full_potential(float const mixing, double const ves_multipole[], int const echo=0) {
        int const nlm = pow2(1 + ellmax);
        int const npt = angular_grid::Lebedev_grid_size(ellmax);
        auto vlm = new double[nlm];
        for(int ts = SMT; ts >= TRU; --ts) { // smooth quantities first, so we can determine vlm
            int const nr = rg[ts]->n;
            auto const mr = full_density[ts].stride();

            { // scope: quantities on the angular grid
                if (echo > 6) printf("# %s quantities on the angular grid are %i * %li = %li\n", label, npt, mr, npt*mr);
                view2D<double> on_grid(2, npt*mr);
                auto const rho_on_grid = on_grid[0];
                auto const vxc_on_grid = on_grid[0];
                auto const exc_on_grid = on_grid[1];

                // transform the lm-index into real-space
                // using an angular grid quadrature, e.g. Lebedev-Laikov grids
                if ((echo > 6) && (SMT == ts)) printf("# %s local smooth density at origin %g a.u.\n",
                                                          label, full_density[ts][00][0]*Y00);
                angular_grid::transform(rho_on_grid, full_density[ts].data(), mr, ellmax, false);
                // envoke the exchange-correlation potential (acts in place)
    //          printf("# envoke the exchange-correlation on angular grid\n");
                for(size_t ip = 0; ip < npt*mr; ++ip) {
                    double const rho = rho_on_grid[ip];
                    double vxc = 0, exc = 0;
                    exc = exchange_correlation::lda_PZ81_kernel(rho, vxc);
                    vxc_on_grid[ip] = vxc; // overwrites rho_on_grid[ip] due to pointer aliasing
                    exc_on_grid[ip] = exc; // only needed if we want to compute the total energy
                } // ip
                // transform back to lm-index
                assert(full_potential[ts].stride() == mr);
                angular_grid::transform(full_potential[ts].data(), vxc_on_grid, mr, ellmax, true);
                { // scope: transform also the exchange-correlation energy
                    view2D<double> exc_lm(nlm, mr);
                    angular_grid::transform(exc_lm.data(), exc_on_grid, mr, ellmax, true);
                    if ((echo > 7) && (SMT == ts)) printf("# %s local smooth exchange-correlation potential at origin is %g %s\n",
                                                            label, full_potential[SMT][00][0]*Y00*eV,_eV);
                    if (echo > 5) {
                        auto const Edc00 = dot_product(nr, full_potential[ts][00], full_density[ts][00], rg[ts]->r2dr); // dot_product with diagonal metric
                        printf("# %s double counting correction  in %s 00 channel %.12g %s\n",
                                  label, (TRU == ts)?"true":"smooth", Edc00*eV,_eV);
                        auto const Exc00 = dot_product(nr, exc_lm[00], full_density[ts][00], rg[ts]->r2dr); // dot_product with diagonal metric
                        printf("# %s exchange-correlation energy in %s 00 channel %.12g %s\n",
                                  label, (TRU == ts)?"true":"smooth", Exc00*eV,_eV);
                    } // echo
                } // scope
            } // scope: quantities on the angular grid

            // solve electrostatics inside the spheres
            view2D<double> Ves(nlm, mr);
            double const q_nucleus = (TRU == ts) ? -Z_core*Y00 : 0; // Z_core = number of protons in the nucleus
            auto   const & rho_aug = (TRU == ts) ? full_density[TRU] : aug_density;
            // solve electrostatics with singularity (q_nucleus) // but no outer boundary conditions (v_lm)
            radial_potential::Hartree_potential(Ves.data(), *rg[ts], rho_aug.data(), rho_aug.stride(), ellmax, q_nucleus);

            if (SMT == ts) {
                add_or_project_compensators<1>(Ves, vlm, ellmax_compensator, rg[SMT], sigma_compensator); // project Ves to compensators
                if (echo > 7) printf("# %s inner integral between normalized compensator and smooth Ves(r) = %g %s\n", label, vlm[0]*Y00*eV,_eV);
                // this seems to high by a factor sqrt(4*pi), ToDo: find out why
//                 scale(vlm, nlm, Y00);

                // but the solution of the 3D problem found that these integrals should be v_lm, therefore
                if (nullptr == ves_multipole) {
                    set(vlm, nlm, 0.); // no correction of the electrostatic potential heights for isolated atoms
                } else {
                    if (echo > 6) printf("# %s v_00 found %g but expected %g %s\n", label, vlm[0]*Y00*eV, ves_multipole[0]*Y00*eV,_eV);
                    scale(vlm, nlm, -1.); add_product(vlm, nlm, ves_multipole, 1.); // vlm := ves_multipole - vlm
                } // no ves_multipole given
            } // smooth only

            if (echo > 7) {
                if (SMT == ts) printf("# %s local smooth electrostatic potential at origin is %g %s\n", label, Ves[00][0]*Y00*eV,_eV);
            }

            add_or_project_compensators<2>(Ves, vlm, ellmax_compensator, rg[ts], sigma_compensator);

            if (SMT == ts) { // debug: project again to see if the correction worked out for the ell=0 channel
                double v_[1];
                add_or_project_compensators<1>(Ves, v_, 0, rg[SMT], sigma_compensator); // project to compensators
                if (echo > 7) {
                    printf("# %s after correction v_00 is %g %s\n", label, v_[00]*Y00*eV,_eV);
                    printf("# %s local smooth electrostatic potential at origin is %g %s\n", label, Ves[00][0]*Y00*eV,_eV);
                    printf("# %s local smooth augmented density at origin is %g a.u.\n", label, aug_density[00][0]*Y00);
                    if (echo > 8) {
                        printf("\n## %s local smooth electrostatic potential and augmented density in a.u.:\n", label);
                        for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                            printf("%g %g %g\n", rg[SMT]->r[ir], Ves[00][ir]*Y00, aug_density[00][ir]*Y00);
                        }   printf("\n\n");
                    } // show radial function of Ves[00]*Y00 to be compared to projections of the 3D electrostatic potential
                } // echo
            } else {
                if (echo > 8) printf("# %s local true electrostatic potential*r at origin is %g (should match -Z=%.1f)\n",
                                label, Ves[00][1]*(rg[TRU]->r[1])*Y00, -Z_core);
            } // ts

            add_product(full_potential[ts].data(), nlm*mr, Ves.data(), 1.0); // add the electrostatic potential, scale_factor=1.0

        } // true and smooth

        // construct the zero_potential V_bar
        std::vector<double> V_smt(rg[SMT]->n);
        set(V_smt.data(), rg[SMT]->n, full_potential[TRU][00] + nr_diff); // copy the tail of the spherical part of the true potential
        set(zero_potential.data(), rg[SMT]->n, 0.0); // init zero
        auto const df = Y00*eV; assert(df > 0); // display factor
        if (echo > 5) printf("# %s match local potential to parabola at R_cut = %g %s, V_tru(R_cut) = %g %s\n",
                    label, rg[SMT]->r[ir_cut[SMT]]*Ang, _Ang, full_potential[TRU][00][ir_cut[TRU]]*df, _eV);
        auto const stat = pseudize_function(V_smt.data(), rg[SMT], ir_cut[SMT], 2);
        if (stat) {
            if (echo > 0) printf("# %s matching procedure for the potential parabola failed! info = %d\n", label, stat);
        } else {
//             if (echo > -1) printf("# local smooth zero_potential:\n");
            for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                zero_potential[ir] = V_smt[ir] - full_potential[SMT][00][ir];
//                 if (echo > -1) printf("%g %g\n", rg[SMT]->r[ir], zero_potential[ir]*Y00);
            } // ir
//             if (echo > -1) printf("\n\n");
            if (echo > 5) printf("# %s potential parabola: V_smt(0) = %g, V_smt(R_cut) = %g %s\n",
                                  label, V_smt[0]*df, V_smt[ir_cut[SMT]]*df, _eV);
            // analyze the zero potential
            double vol = 0, Vint = 0, r1Vint = 0, r2Vint = 0;
            for(int ir = ir_cut[SMT]; ir < rg[SMT]->n; ++ir) {
                auto const r  = rg[SMT]->r[ir];
                auto const dV = rg[SMT]->r2dr[ir];
                vol    +=                    dV;
                Vint   += zero_potential[ir]*dV;
                r1Vint += zero_potential[ir]*dV*r;
                r2Vint += zero_potential[ir]*dV*pow2(r);
            } // ir
            if (echo > 5) printf("# %s zero potential statistics = %g %g %g %s\n",
                        label, Vint/vol*eV, r1Vint/(vol*r_cut)*eV, r2Vint/(vol*pow2(r_cut))*eV, _eV);
                // these numbers should be small since they indicate that V_bar is localized inside the sphere
                // and how much V_smt deviates from V_tru ouside the sphere
        } // pseudization successful
        if (echo > 5) printf("# %s zero potential: V_bar(0) = %g, V_bar(R_cut) = %g, V_bar(R_max) = %g %s\n",
            label, zero_potential[0]*df, zero_potential[ir_cut[SMT]]*df, zero_potential[rg[SMT]->n - 1]*df, _eV);

        // add spherical zero potential for SMT==ts and 0==lm
        add_product(full_potential[SMT][00], rg[SMT]->n, zero_potential.data(), 1.0);

        // feed back spherical part of the true potential into the spherical true potential r*V
        // which determines core states and true partial waves
        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            scale(potential[ts].data(), rg[ts]->n, 1. - mixing);
            add_product(potential[ts].data(), rg[ts]->n, full_potential[ts][00], rg[ts]->r, Y00*mixing);
        } // ts true and smooth

        if (0) { // scope: test: use the spherical routines from atom_core::rad_pot(output=r*V(r), input=rho(r)*4pi)
            std::vector<double> rho4pi(rg[TRU]->n);
            set(rho4pi.data(), rg[TRU]->n, full_density[TRU][00], Y00inv);
            printf("\n# WARNING: use rad_pot to construct the r*V_tru(r)\n\n");
            atom_core::rad_pot(potential[TRU].data(), *rg[TRU], rho4pi.data(), Z_core);
        } // scope

    } // update_full_potential

    void update_matrix_elements(int const echo=0) {
        int const nlm = pow2(1 + ellmax);
        int const mlm = pow2(1 + numax);
        int const nln = sho_radial::nSHO_radial(numax);
        int const nSHO = sho_tools::nSHO(numax);
        int const nlmn = nSHO;
        initialize_Gaunt();

        // first generate the matrix elemnts in the valence basis
        //    overlap[iln][jln] and potential_lm[iln][jln]
        // then, generate the matrix elements in the radial representation
        //    overlap[ilmn][jlmn] and hamiltonian[ilmn][jlmn]
        //    where hamiltonian[ilmn][jlmn] = kinetic_energy_deficit_{iln jln}
        //            + sum_lm Gaunt_{lm ilm jlm} * potential_{lm iln jln}
        // then, transform the matrix elements into the Cartesian representation using sho_unitary
        //    overlap[iSHO][jSHO] and hamiltonian[iSHO][jSHO]

        view4D<double> potential_ln(nlm, TRU_AND_SMT, nln, nln); // get memory, // emm1-emm2-degenerate
        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            int const nr = rg[ts]->n;
            std::vector<double> wave_pot_r2dr(nr);
            for(int ell = 0; ell <= ellmax; ++ell) {
                for(int emm = -ell; emm <= ell; ++emm) {
                    int const lm = solid_harmonics::lm_index(ell, emm);
                    assert(lm < nlm);
                    for(int iln = 0; iln < nln; ++iln) {
                        auto const wave_i = partial_wave[iln].wave[ts];
                        product(wave_pot_r2dr.data(), nr, wave_i, full_potential[ts][lm], rg[ts]->r2dr);
                        for(int jln = 0; jln < nln; ++jln) {
                            auto const wave_j = partial_wave[jln].wave[ts];
                            potential_ln(lm,ts,iln,jln) = dot_product(nr, wave_pot_r2dr.data(), wave_j);
                        } // jln
                    } // iln
                } // emm
            } // ell
        } // ts: true and smooth

        view2D<double> hamiltonian_lmn(nlmn, nlmn); // get memory
        view2D<double> overlap_lmn(nlmn, nlmn);
        set(hamiltonian_lmn, nlmn, 0.0); // clear
        set(overlap_lmn,     nlmn, 0.0); // clear
        for(auto gnt : gaunt) {
            int const lm = gnt.lm, lm1 = gnt.lm1, lm2 = gnt.lm2; auto G = gnt.G;
            if (00 == lm) G = Y00*(lm1 == lm2); // make sure that G_00ij = delta_ij*Y00
            if (lm1 < mlm && lm2 < mlm) {
                if (lm < nlm) {
                    for(int ilmn = lmn_begin[lm1]; ilmn < lmn_end[lm1]; ++ilmn) {
                        int const iln = ln_index_list[ilmn];
                        for(int jlmn = lmn_begin[lm2]; jlmn < lmn_end[lm2]; ++jlmn) {
                            int const jln = ln_index_list[jlmn];
                            hamiltonian_lmn[ilmn][jlmn] +=
                              G * ( potential_ln(lm,TRU,iln,jln)
                              	  - potential_ln(lm,SMT,iln,jln) );
                        } // jlmn
                    } // ilmn
                } // lm
            } // limits
        } // gnt

        // add the kinetic_energy deficit to the hamiltonian
        if (echo > 7) printf("\n# %s Hamiltonian elements %s-ordered in %s:\n",
                        label, sho_tools::SHO_order2string(sho_tools::order_lmn).c_str(), _eV);
        for(int ilmn = 0; ilmn < nlmn; ++ilmn) {
            int const iln = ln_index_list[ilmn];
            int const ilm = lm_index_list[ilmn];
            if (echo > 7) printf("# %s hamiltonian elements for ilmn=%3i  ", label, ilmn);
            for(int jlmn = 0; jlmn < nlmn; ++jlmn) {
                int const jlm = lm_index_list[jlmn];
                int const jln = ln_index_list[jlmn];
                if (ilm == jlm) {
                    hamiltonian_lmn[ilmn][jlmn] += ( kinetic_energy(TRU,iln,jln) - kinetic_energy(SMT,iln,jln) );
                    overlap_lmn[ilmn][jlmn] =
                        ( charge_deficit(0,TRU,iln,jln)
                        - charge_deficit(0,SMT,iln,jln) ); // ell=0
                } // diagonal in lm, offdiagonal in nrn
                if ((echo > 7)) printf(" %7.3f", true_norm[iln]*hamiltonian_lmn[ilmn][jlmn]*true_norm[jln]*eV);
            } // jlmn
            if ((echo > 7)) printf("\n");
        } // ilmn

        if (echo > 8) {
            printf("\n# %s lmn-based Overlap elements:\n", label);
            for(int ilmn = 0; ilmn < nlmn; ++ilmn) {      int const iln = ln_index_list[ilmn];
                printf("# %s overlap elements for ilmn=%3i  ", label, ilmn);
                for(int jlmn = 0; jlmn < nlmn; ++jlmn) {  int const jln = ln_index_list[jlmn];
                    printf(" %g", true_norm[iln]*true_norm[jln]*overlap_lmn[ilmn][jlmn]);
                }   printf("\n");
            } // ilmn
        } // echo

        if (1) { // scope: check if averaging over emm gives back the same operators

            view2D<double> hamiltonian_ln(nln, nln); // get memory
            view2D<double>     overlap_ln(nln, nln); // get memory
            for(int i01 = 0; i01 < 2; ++i01) {
                auto const & input_lmn = i01 ?  overlap_lmn :  hamiltonian_lmn;
                auto const   label_inp = i01 ? "overlap"    : "hamiltonian"   ;
                auto const & result_ln = i01 ?  overlap_ln  :  hamiltonian_ln ;
                scattering_test::emm_average(result_ln.data(), input_lmn.data(), (int)numax, nn);
                if (echo > 4) {
                    for(int iln = 0; iln < nln; ++iln) {
                        printf("# %s emm-averaged%3i %s ", label, iln, label_inp);
                        for(int jln = 0; jln < nln; ++jln) {
                            printf(" %11.6f", true_norm[iln]*true_norm[jln]*result_ln[iln][jln]);
                        }   printf("\n");
                    }   printf("\n");
                } // echo
            } // i01

            if (0) { // Warning: can only produce the same eigenenergies if potentials are converged:
                     //             Y00*r*full_potential[ts][00](r) == potential[ts](r)
                if (echo > 1) printf("\n\n# %s perform a diagonalization of the pseudo Hamiltonian\n\n", label);
                double const V_rmax = full_potential[SMT][00][rg[SMT]->n - 1]*Y00;
                auto Vsmt = std::vector<double>(rg[SMT]->n, 0);
                { // scope: prepare a smooth local potential which goes to zero at Rmax
                    for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                        Vsmt[ir] = full_potential[SMT][00][ir]*Y00 - V_rmax;
                    } // ir
                } // scope

                if (echo > 1) printf("\n# %s %s eigenstate_analysis\n\n", label, __func__);
                scattering_test::eigenstate_analysis // find the eigenstates of the spherical Hamiltonian
                  (*rg[SMT], Vsmt.data(), sigma, (int)numax + 1, nn, numax, hamiltonian_ln.data(), overlap_ln.data(), 384, V_rmax, label, echo);
            } else if (echo > 0) printf("\n# eigenstate_analysis deactivated for now! %s %s:%i\n\n", __func__, __FILE__, __LINE__);

            if (0) { // Warning: can only produce the same eigenenergies if potentials are converged:
                     //             Y00*r*full_potential[ts][00](r) == potential[ts](r)
                double* rV[TRU_AND_SMT]; std::vector<double> rV_vec[TRU_AND_SMT];
                for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                    int const nr = rg[ts]->n;
                    rV_vec[ts] = std::vector<double>(nr); rV[ts] = rV_vec[ts].data();
                    product(rV[ts], nr, full_potential[ts][00], rg[ts]->r, Y00);
                } // ts
                if (echo > 1) printf("\n# %s %s logarithmic_derivative\n\n", label, __func__);
                scattering_test::logarithmic_derivative // scan the logarithmic derivatives
                  (rg, rV, sigma, (int)numax + 1, nn, numax, hamiltonian_ln.data(), overlap_ln.data(), logder_energy_range, label, echo);
            } else if (echo > 0) printf("\n# logarithmic_derivative deactivated for now! %s %s:%i\n\n", __func__, __FILE__, __LINE__);

        } // scope

        set(hamiltonian, nSHO, 0.0); // clear
        set(overlap,     nSHO, 0.0); // clear
        // Now transform _lmn quantities to Cartesian representations using sho_unitary
        transform_SHO(hamiltonian.data(), hamiltonian.stride(), hamiltonian_lmn.data(), hamiltonian_lmn.stride(), false);
        transform_SHO(    overlap.data(),     overlap.stride(),     overlap_lmn.data(),     overlap_lmn.stride(), false);
        // Mind that this transform is unitary and assumes square-normalized SHO-projectors
        // ... which requires proper normalization factors f(i)*f(j) to be multiplied in, see sho_projection::sho_prefactor

        if (echo > 1) {
            printf("\n# %s SHO-transformed Hamiltonian elements (%s-order) in %s:\n",
                       label, sho_tools::SHO_order2string(sho_tools::order_Ezyx).c_str(), _eV);
            for(int iSHO = 0; iSHO < nSHO; ++iSHO) {
                printf("# %s hamiltonian elements for iSHO=%3d  ", label, iSHO); // ToDo: show the nx,ny,nz quantum numbers
                for(int jSHO = 0; jSHO < nSHO; ++jSHO) {
                    printf(" %7.3f", hamiltonian[iSHO][jSHO]*eV); // true_norm cannot be used here
                    // since hamiltonian is given in the Cartesian representation!
                } // jSHO
                printf("\n");
            } // iSHO
        } // echo

    } // update_matrix_elements


    void check_spherical_matrix_elements(int const echo, view3D<double> & aHSm) const {
        if (echo < 1) return; // this function only plots information to the logg
        // check if the emm-averaged Hamiltonian elements produce the same scattering properties
        // as the spherical part of the full potential
        int const nln = sho_tools::nSHO_radial(numax);

        view3D<double> potential_ln(TRU_AND_SMT, nln, nln); // emm1-emm2-degenerate, no emm0-dependency
        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            int const nr = rg[ts]->n;
            std::vector<double> wave_pot_r2dr(nr); // temporary product of partial wave, potential and metric
            for(int iln = 0; iln < nln; ++iln) {
                auto const wave_i = partial_wave[iln].wave[ts];
                // potential is defined as r*V(r), so we only need r*dr to get r^2*dr as integration weights
                product(wave_pot_r2dr.data(), nr, wave_i, potential[ts].data(), rg[ts]->rdr);
                for(int jln = 0; jln < nln; ++jln) {
                    auto const wave_j = partial_wave[jln].wave[ts];
                    potential_ln(ts,iln,jln) = dot_product(nr, wave_pot_r2dr.data(), wave_j);
                } // jln
            } // iln
        } // ts: true and smooth


        // view2D<double> hamiltonian_ln(nln, nln), overlap_ln(nln, nln); // get memory
        auto hamiltonian_ln = aHSm[0], overlap_ln = aHSm[1];
        { // scope
            for(int iln = 0; iln < nln; ++iln) {
                for(int jln = 0; jln < nln; ++jln) {
                    hamiltonian_ln[iln][jln] = ( kinetic_energy(TRU,iln,jln)
                                               - kinetic_energy(SMT,iln,jln) )
                                             + ( potential_ln(TRU,iln,jln)
                                               - potential_ln(SMT,iln,jln) );
                    int const ell = 0;
                    overlap_ln[iln][jln] = ( charge_deficit(ell,TRU,iln,jln)
                                           - charge_deficit(ell,SMT,iln,jln) );
                } // jln
            } // iln
        } // scope

        if (1) { // show atomic matrix elements
            std::vector<int8_t> ells(nln, -1);
            {
                for(int ell = 0; ell <= numax; ++ell) {
                    for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                        int const iln = sho_tools::ln_index(numax, ell, nrn);
                        ells[iln] = ell;
                    } // nrn
                } // ell
            }

            if (echo > 4) {
                printf("\n");
                for(int i01 = 0; i01 < 2; ++i01) {
                    auto const & input_ln = i01 ?  overlap_ln :  hamiltonian_ln;
                    auto const   label_qn = i01 ? "overlap"   : "hamiltonian" ;
                    for(int iln = 0; iln < nln; ++iln) {
                        printf("# %s spherical%3i %s ", label, iln, label_qn);
                        for(int jln = 0; jln < nln; ++jln) {
                            printf(" %g", (ells[iln] == ells[jln]) ? true_norm[iln]*input_ln[iln][jln]*true_norm[jln] : 0);
                        }   printf("\n");
                    }   printf("\n");
                } // i01
            } // echo
        } // 1

        if (1) {
            double const V_rmax = potential[SMT][rg[SMT]->n - 1]*rg[SMT]->rinv[rg[SMT]->n - 1];
            auto Vsmt = std::vector<double>(rg[SMT]->n, 0);
            { // scope: prepare a smooth local potential which goes to zero at Rmax
                for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                    Vsmt[ir] = potential[SMT][ir]*rg[SMT]->rinv[ir] - V_rmax;
                } // ir
            } // scope

            if (echo > 1) printf("\n# %s %s eigenstate_analysis\n\n", label, __func__);
            scattering_test::eigenstate_analysis // find the eigenstates of the spherical Hamiltonian
              (*rg[SMT], Vsmt.data(), sigma, (int)numax + 1, nn, numax, hamiltonian_ln.data(), overlap_ln.data(), 384, V_rmax, label, echo);
        } else if (echo > 0) printf("\n# eigenstate_analysis deactivated for now! %s %s:%i\n\n", __func__, __FILE__, __LINE__);

        if (1) {
            if (echo > 1) printf("\n# %s %s logarithmic_derivative\n\n", label, __func__);
            double const *rV[TRU_AND_SMT] = {potential[TRU].data(), potential[SMT].data()};
            scattering_test::logarithmic_derivative // scan the logarithmic derivatives
              (rg, rV, sigma, (int)numax + 1, nn, numax, hamiltonian_ln.data(), overlap_ln.data(), logder_energy_range, label, echo);
        } else if (echo > 0) printf("\n# logarithmic_derivative deactivated for now! %s %s:%i\n\n", __func__, __FILE__, __LINE__);

    } // check_spherical_matrix_elements

    
    void create_isolated_density_matrix(view2D<double> & density_matrix, 
          sho_tools::SHO_order_t & order, 
          view3D<double> const & aHSm, // atomic emm-degenrate matrix elements (h0s1, nln, nln)
          int const echo=9) const {
        int const nSHO = sho_tools::nSHO(numax);
        int const nln = sho_tools::nSHO_radial(numax);
        order = sho_tools::order_lmn;
        if (echo > -1) printf("# %s %s for an %d x %d density matrix in %s_order\n", label, __func__, nSHO, nSHO,
                                  sho_tools::SHO_order2string(order).c_str());
        std::vector<double> rwave(rg[SMT]->n, 0.0), ff(rg[SMT]->n);
        for(int ics = 0; ics < ncorestates; ++ics) {
            auto const occ = -core_state[ics].occupation;
            int const enn = core_state[ics].enn;
            int const ell = core_state[ics].ell;
            if (occ > 0) {
                if (ell > numax) {
                    warn("Unable to represent %c-states with numax=%d", ellchar[ell], numax);
                } else {
                    double const energy = core_state[ics].energy;
                    int const iln_off = sho_tools::ln_index(numax, ell, 0);
                    double const *const aHm = aHSm[0].data() + iln_off*nln + iln_off;
                    double const *const aSm = aHSm[1].data() + iln_off*nln + iln_off;
                    if (echo > -1) printf("# %s %s found occ %g for %i%c-state at %g %s\n", 
                                              label, __func__, occ, enn,ellchar[ell], energy*eV, _eV);
                    // solve for the pseudo state at the energy of the eigenstate of the true spherical potential
                    scattering_test::generalized_node_count_SMT( // ToDo: list above in the includes
                        *rg[SMT], potential[SMT].data(), ell, energy, rwave.data(), ff.data(), -1,
                        view2D<double>(projectors[iln_off], projectors.stride()), nn[ell],
                        aHm, aSm, aHSm.stride(), echo);
                    // now find the correct normalization of rwave:
                    double cprj[8];
                    for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                        cprj[nrn] = dot_product(rg[SMT]->n, projectors[iln_off + nrn], rwave.data(), rg[SMT]->dr);
                        if (echo > -1) printf("# %s %i%c-state has projection coeff #%i %g\n", 
                                                  label, enn,ellchar[ell], nrn, cprj[nrn]);
                    } // nrn
                    // normalization as 1 == <rwave| ( 1 + |p> charge_deficit <p| ) |rwave>
                    double const l2norm = dot_product(rg[SMT]->n, rwave.data(), rwave.data(), rg[SMT]->dr);
                    double norm{l2norm};
                    for(int irn = 0; irn < nn[ell]; ++irn) {
                        for(int jrn = 0; jrn < nn[ell]; ++jrn) {
                            norm += ( charge_deficit(0,TRU,iln_off + irn,iln_off + jrn)
                                    - charge_deficit(0,SMT,iln_off + irn,iln_off + jrn) )
                                    * cprj[irn] * cprj[jrn];
                        } // jrn
                    } // irn
                    double const f = 1./std::sqrt(norm);
                    scale(rwave.data(), rwave.size(), f);
                    scale(cprj, nn[ell], f);
                    for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                        if (echo > -1) printf("# %s %i%c-state has normalized projection coeff #%i %g\n", 
                                                  label, enn,ellchar[ell], nrn, cprj[nrn]);
                    } // nrn

                    // now set the density_matrix in order_lmn
                    for(int emmi = -ell; emmi <= ell; ++emmi) {
                        for(int irn = 0; irn < nn[ell]; ++irn) {
                            int const ilmn = sho_tools::lmn_index(numax, ell, emmi, irn);
                            for(int emmj = -ell; emmj <= ell; ++emmj) {
                                for(int jrn = 0; jrn < nn[ell]; ++jrn) {
                                    int const jlmn = sho_tools::lmn_index(numax, ell, emmj, jrn);
                                    density_matrix[ilmn][jlmn] = occ * cprj[irn] * cprj[jrn];
                                } // jrn
                            } // emm_j
                        } // irn
                    } // emm_i

                } // ell too large
            } // occupied
        } // ics

        // show the density matrix
        if (echo > 1) {
            printf("# %s %s:\n", label, __func__);
            for(int ilmn = 0; ilmn < nSHO; ++ilmn) {
                printf("# %s dm %02i ", label, ilmn);
                for(int jlmn = 0; jlmn < nSHO; ++jlmn) {
                    printf(" %7.3f", density_matrix[ilmn][jlmn]);
                } // jlmn
                printf("\n");
            } // ilmn
        } // echo

    } // create_isolated_density_matrix
    

    void update_density(float const mixing, int const echo=0) {
//         if (echo > 2) printf("\n# %s\n", __func__);
        update_core_states(mixing, echo);
        update_partial_waves(echo); // create new partial waves for the valence description
        update_charge_deficit(echo); // update quantities derived from the partial waves
        int const nln = sho_radial::nSHO_radial(numax);
        view3D<double> aHSm(2, nln, nln, 0.0);
        check_spherical_matrix_elements(echo, aHSm); // check scattering properties for emm-averaged Hamiltonian elements
        int const nSHO = sho_tools::nSHO(numax);
        view2D<double> density_matrix(nSHO, nSHO, 0.0); // get memory
        auto dm_order = sho_tools::order_Ezyx;
        if (1) create_isolated_density_matrix(density_matrix, dm_order, aHSm);
        int const lmax = std::max(ellmax, ellmax_compensator);
        int const mlm = pow2(1 + lmax);
        view3D<double> rho_tensor(mlm, nln, nln, 0.0); // get memory
        get_rho_tensor(rho_tensor, density_matrix, dm_order, echo);
        update_full_density(rho_tensor, echo);
    } // update_density

    // ==============
    // between update_density and update_potential we need to
    // export qlm_compensator, solve 3D electrostatic problem, return here with ves_multipoles
    // ==============

    void update_potential(float const mixing, double const ves_multipoles[], int const echo=0) {
        if (echo > 2) printf("\n# %s %s\n", label, __func__);
        update_full_potential(mixing, ves_multipoles, echo);
        update_matrix_elements(echo); // this line does not compile with icpc (ICC) 19.0.2.187 20190117
    } // update_potential

    template <char what>
    status_t get_smooth_spherical_quantity(double qnt[] // result array: function on an r2-grid
        , float const ar2, int const nr2 // r2-grid parameters
        , int const echo=1) const { // logg level
        char const *qnt_name   = ('c' == what) ? "core_density"     : "zero_potential";
        auto const &qnt_vector = ('c' == what) ?  core_density[SMT] :  zero_potential;
        if (echo > 8) printf("# %s call transform_to_r2_grid(%p, %.1f, %d, %s=%p, rg=%p)\n",
                        label, (void*)qnt, ar2, nr2, qnt_name, (void*)qnt_vector.data(), (void*)rg[SMT]);
        return bessel_transform::transform_to_r2_grid(qnt, ar2, nr2, qnt_vector.data(), *rg[SMT], echo);
    } // get_smooth_spherical_quantity
    
    radial_grid_t* get_smooth_radial_grid(int const echo=0) const { return rg[SMT]; }

  }; // class LiveAtom


namespace single_atom {

  status_t update(int const na, float const *Za, float const *ion,
                  radial_grid_t **rg, double *sigma_cmp,
                  double **rho, double **qlm, double **vlm, int *lmax_vlm, int *lmax_qlm, double **zero_pot) {

      static int echo = -9;
      if (echo == -9) echo = control::get("single_atom.echo", 0.); // initialize only on the 1st call to update()

      static LiveAtom **a=nullptr;

      if (nullptr == a) {
//        SimpleTimer timer(__FILE__, __LINE__, "LiveAtom-constructor");
          a = new LiveAtom*[na];
          int  const default_numax = 3;
          bool const transfer2valence = false;
          for(int ia = 0; ia < na; ++ia) {
              a[ia] = new LiveAtom(Za[ia], default_numax, transfer2valence, ion[ia], ia, echo);
          } // ia
      } // a has not been initialized

      if (na < 0) {
          for(int ia = 0; ia < -na; ++ia) {
              a[ia]->~LiveAtom(); // envoke destructor
          } // ia
          delete[] a; a = nullptr;
      } // cleanup

      for(int ia = 0; ia < na; ++ia) {

          if (nullptr != rho) {
              int const nr2 = 1 << 12; float const ar2 = 16.f; // rcut = 15.998 Bohr
              rho[ia] = new double[nr2];
              a[ia]->get_smooth_spherical_quantity<'c'>(rho[ia], ar2, nr2);
          } // get the smooth core densities

          if (nullptr != zero_pot) {
              int const nr2 = 1 << 12; float const ar2 = 16.f; // rcut = 15.998 Bohr
              zero_pot[ia] = new double[nr2];
              a[ia]->get_smooth_spherical_quantity<'V'>(zero_pot[ia], ar2, nr2);
          } // get the zero_potential

          if (nullptr != rg) rg[ia] = a[ia]->get_smooth_radial_grid(); // pointers to smooth radial grids

          if (nullptr != sigma_cmp) sigma_cmp[ia] = a[ia]->sigma_compensator; // spreads of the compensators // ToDo: use a getter function

          if (nullptr != qlm) set(qlm[ia], pow2(1 + a[ia]->ellmax_compensator), a[ia]->qlm_compensator.data()); // copy compensator multipoles

          if (nullptr != vlm) {
              a[ia]->update_potential(.5f, vlm[ia], echo); // set electrostatic multipole shifts
              a[ia]->update_density(.5f, echo);
          } // vlm

          if (nullptr != lmax_vlm) lmax_vlm[ia] = a[ia]->ellmax;
          if (nullptr != lmax_qlm) lmax_qlm[ia] = a[ia]->ellmax_compensator;

      } // ia

      return 0;
  } // update


#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_compensator_normalization(int const echo=5) {
      if (echo > 1) printf("\n# %s: %s\n", __FILE__, __func__);
      auto const rg = radial_grid::create_exponential_radial_grid(512, 2.0);
      int const nr = rg->n, lmax = 0, nlm = pow2(1 + lmax);
      std::vector<double> qlm(nlm, 0.0);
      view2D<double> cmp(1, nr);
      double maxdev{0};
      for(double sigma = 0.5; sigma < 2.25; sigma *= 1.1) {
          set(qlm.data(), nlm, 0.0); qlm[0] = 1.0;
          set(cmp, 1, 0.0); // clear
          add_or_project_compensators<0>(cmp, qlm.data(), lmax, rg, sigma, 0); // add normalized compensator
//        add_or_project_compensators<1>(cmp, qlm.data(), lmax, rg, sigma, 0); // project
//        if (echo > 0) printf("# %s: square-norm of normalized compensator with sigma = %g is %g\n", __func__, sigma, qlm[0]);
          add_or_project_compensators<3>(cmp, qlm.data(), lmax, rg, sigma, 0); // test normalization
          maxdev = std::max(maxdev, std::abs(qlm[0] - 1.0));
          if (echo > 4) printf("# %s: for sigma = %g is 1 + %.1e\n", __func__, sigma, qlm[0] - 1);
      } // sigma
      if (echo > 2) printf("# %s: largest deviation is %.1e\n", __func__, maxdev);
      return (maxdev > 1e-15);
  } // test_compensator_normalization

  int test_LiveAtom(int const echo=9) {
    int const numax = control::get("single_atom.test.numax", 3.); // default 3: ssppdf
    bool const t2v = (control::get("single_atom.test.transfer.to.valence", 0.) > 0); //
    if (echo > 0) printf("\n# %s: new struct LiveAtom has size %ld Byte\n\n", __FILE__, sizeof(LiveAtom));
//     for(int Z = 0; Z <= 109; ++Z) { // all elements
//     for(int Z = 109; Z >= 0; --Z) { // all elements backwards
//        if (echo > 1) printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
    {   double const Z = control::get("single_atom.test.Z", 29.); // default copper
        double const ion = control::get("single_atom.test.ion", 0.); // default neutral
        if (echo > 1) printf("\n# Z = %d\n", Z);
        LiveAtom a(Z, numax, t2v, ion, -1, echo); // envoke constructor
    } // Z
    return 0;
  } // test_LiveAtom

  status_t all_tests(int const echo) {
    status_t stat{0};
    stat += test_compensator_normalization(echo);
    stat += test_LiveAtom(echo);
    return stat;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace single_atom
