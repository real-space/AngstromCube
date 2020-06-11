#include <cstdio> // printf
#include <cmath> // std::sqrt, std::abs
#include <cassert> // assert
#include <algorithm> // std::max

#include "single_atom.hxx"

#include "radial_grid.h" // radial_grid_t
#include "radial_grid.hxx" // ::create_default_radial_grid, ::destroy_radial_grid, ::dot_product
#include "radial_eigensolver.hxx" // ::shooting_method
#include "radial_potential.hxx" // ::Hartree_potential
#include "angular_grid.hxx" // ::transform, ::Lebedev_grid_size
#include "radial_integrator.hxx" // ::integrate_outwards
#include "exchange_correlation.hxx" // ::lda_PZ81_kernel
#include "inline_tools.hxx" // align<nbits>
#include "sho_unitary.hxx" // ::Unitary_SHO_Transform<real_t>
#include "solid_harmonics.hxx" // ::lm_index, ::Y00, ::Y00inv
#include "atom_core.hxx" // ::initial_density, ::rad_pot, ::nl_index
#include "sho_tools.hxx" // ::lnm_index, ::nSHO, ??? some more, ::nSHO_radial
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t, emm_Degenerate, spin_QN_t, spin_Degenerate
#include "energy_level.hxx" // TRU, SMT, TRU_AND_SMT, core_level_t, valence_level_t
#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // pow2, pow3, set, scale, product, add_product, intpow
#include "simple_math.hxx" // invert
#include "simple_timer.hxx" // SimpleTimer
#include "bessel_transform.hxx" // ::transform_to_r2_grid
#include "scattering_test.hxx" // ::eigenstate_analysis, ::emm_average
#include "linear_algebra.hxx" // ::linear_solve, ::generalized_eigenvalues
#include "data_view.hxx" // view2D<T>, view3D<T>
#include "lossful_compression.hxx" // print_compressed
#include "control.hxx" // ::get
#include "sigma_config.hxx" // ::get, element_t

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
      std::vector<double> eigs(n);
      auto const info = linear_algebra::generalized_eigenvalues(n, A.data(), A.stride(),
                                                                   B.data(), B.stride(), eigs.data());
      if (lowest) *lowest = eigs[0];
      return info;
  } // minimize_curvature


  int constexpr ELLMAX=7;
  char const ellchar[] = "spdfghijklmno";
//  char const c0s1v2_char[] = "csv?";  // 0:core, 1:semicore, 2:valence   (((unused)))
  char const c0s1v2_name[][10] = {"core", "semicore", "valence"};  // 0:core, 1:semicore, 2:valence

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
        auto const sig2inv = .5/(sigma*sigma);
        if (echo > 0) printf("# sigma = %g\n", sigma);
        std::vector<double> rl(nr), rlgauss(nr);
        for(int ell = 0; ell <= lmax; ++ell) { // serial!
            double norm{0};
            for(int ir = 0; ir < nr; ++ir) {
                auto const r = rg->r[ir];
                if (0 == ell) {
                    rl[ir] = 1; // start with r^0
                    rlgauss[ir] = std::exp(-sig2inv*r*r);
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

    
    template <typename real_t>
    status_t Lagrange_derivatives(unsigned const n, real_t const y[], double const x[], double const at_x0
                        , double *value_x0, double *deriv1st, double *deriv2nd) { // zeroth, first and second derivative at x0
        double d0{0}, d1{0}, d2{0};
        //
        //  L(x) = sum_j y_j                                                           prod_{m!=j}             (x - x_m)/(x_j - x_m)
        //
        //  L'(x) = sum_j y_j sum_{k!=j} 1/(x_j - x_k)                                 prod_{m!=j, m!=k}       (x - x_m)/(x_j - x_m)
        //
        //  L''(x) = sum_j y_j sum_{k!=j} 1/(x_j - x_k) sum_{l!=j, l!=k} 1/(x_j - x_l) prod_{m!=j, m!=k, m!=l} (x - x_m)/(x_j - x_m)
        
        for (int j = 0; j < n; ++j) {
            double d0j{1}, d1j{0}, d2j{0};
            for (int k = 0; k < n; ++k) {
                if (k != j) {
                    if (std::abs(x[j] - x[k]) < 1e-9*std::max(std::abs(x[j]), std::abs(x[j]))) return 1; // status instable
                    double d1l{1}, d2l{0};
                    for (int l = 0; l < n; ++l) {
                        if (l != j && l != k) {
                            
                            double d2m{1};
                            for (int m = 0; m < n; ++m) {
                                if (m != j && m != k && m != l) {
                                    d2m *= (at_x0 - x[m])/(x[j] - x[m]);
                                } // exclude
                            } // m
                            
                            d1l *= (at_x0 - x[l])/(x[j] - x[l]);
                            d2l += d2m/(x[j] - x[l]);
                        } // exclude
                    } // l
                    
                    d0j *= (at_x0 - x[k])/(x[j] - x[k]);
                    d1j += d1l/(x[j] - x[k]);
                    d2j += d2l/(x[j] - x[k]);
                } // exclude
            } // k
            d0 += d0j * y[j];
            d1 += d1j * y[j];
            d2 += d2j * y[j];
        } // j
        
        if(value_x0) *value_x0 = d0;
        if(deriv1st) *deriv1st = d1;
        if(deriv2nd) *deriv2nd = d2;
        
        return 0; // success
    } // Lagrange_derivatives
   
    
  class LiveAtom {
  public:
      // ToDo: separate everything which is energy-parameter-set dependent and group it into a class valence_t (or some better name)

      // general config
      int32_t atom_id; // global atom identifyer
      double Z_core; // number of protons in the core
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
      int8_t numax; // limit of the SHO projector quantum numbers
      double sigma, sigma_inv; // spread of the SHO projectors and its inverse
      uint8_t nn[1+ELLMAX+2]; // number of projectors and partial waves used in each ell-channel
      std::vector<valence_level_t> partial_wave;
      // the following quantities are energy-parameter-set dependent and spin-resolved (nspins=1 or =2)
      view2D<double> hamiltonian, overlap; // matrices [nSHO][>=nSHO]
      view3D<double> kinetic_energy; // tensor [TRU_AND_SMT][nln][nln]
      view4D<double> charge_deficit; // tensor [1 + ellmax_compensator][TRU_AND_SMT][nln][nln]
      view2D<double> projectors; // [nln][nr_smt] r*projectors, depend only on sigma and numax
      view3D<double> partial_waves[TRU_AND_SMT]; // matrix [wave0_or_wKin1][nln][nr], valence states point into this
      view3D<double> true_core_waves; // matrix [wave0_or_wKin1][nln][nr], core states point into this
      std::vector<double> true_norm; // vector[nln] for display of partial wave results
      // end of energy-parameter-set dependent members
      std::vector<double> zero_potential; // PAW potential shape correction, potentially energy-parameter-set dependent

      view2D<double> aug_density; // augmented density, core + valence + compensation, (1+ellmax)^2 radial functions
      int ncorestates; // for emm-Degenerate representations, 20 (or 32 with spin-orbit) core states are maximum
      int nspins; // 1 or 2 (order 0z) or 4 (order 0zxy)

      // spin-resolved members
      double csv_charge[3];
      std::vector<core_level_t> core_state; // 20 core states are the usual max., 32 core states are enough if spin-orbit-interaction is on
      std::vector<double> core_density[TRU_AND_SMT]; // spherical core density*4pi, no Y00 factor
      float mix_spherical_valence_density; // 1.f: use the spherical valence density only, 0.f: use the valence density from partial waves, mixtures possible.
      std::vector<double> spherical_valence_density[TRU_AND_SMT]; // spherical valence density*4pi, no Y00 factor, in use?
      view2D<double> full_density[TRU_AND_SMT]; // total density, core + valence, (1+ellmax)^2 radial functions
      view2D<double> full_potential[TRU_AND_SMT]; // (1+ellmax)^2 radial functions
      std::vector<double> potential[TRU_AND_SMT]; // spherical potential r*V(r), no Y00 factor, used for the generation of partial waves

      double core_charge_deficit; // in units of electrons
      double spherical_valence_charge_deficit; // in units of electrons

      double logder_energy_range[3]; // [start, increment, stop]
      char   partial_wave_char[32]; // [iln]
      double partial_wave_energy_split[8]; // [ell]

      view2D<double> unitary_Ezyx_lmn; // unitary sho transformation matrix [order_Ezyx][order_lmn], stride=nSHO(numax)

      bool gaunt_init;
      std::vector<gaunt_entry_t> gaunt;
      // ToDo: check if all of these 4 lists are used or can be replace but some sho_tool::???_index()
      std::vector<int16_t> ln_index_list;
      std::vector<int16_t> lm_index_list;
      std::vector<int16_t> lmn_begin;
      std::vector<int16_t> lmn_end;

  public:

    // constructor method:
    LiveAtom(double const Z_protons
            , int const nu_max=3
            , bool const atomic_valence_density=false
            , double const ionization=0
            , int const global_atom_id=-1
            , int const echo=0 // logg level for this constructor method only
            ) : atom_id{global_atom_id} 
              , Z_core{Z_protons}
              , mix_spherical_valence_density{atomic_valence_density ? 1.f : 0.f}
              , gaunt_init{false}
        { // constructor

        if (atom_id >= 0) { sprintf(label, "a#%d", atom_id); } else { label[0]=0; }
        if (echo > 0) printf("\n\n#\n# %s LiveAtom with %g protons, ionization=%g\n", label, Z_core, ionization);

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

        double custom_occ[32];
        set(custom_occ, 32, 0.0);
        set(nn, 1+ELLMAX+2, uint8_t(0)); // clear
        char const custom_config = *(control::get("single_atom.from.sigma.config", "0"));
        int const nn_limiter = control::get("single_atom.nn.limit", 2);
        
        if ('1' == custom_config) {
            if (echo > 0) printf("# %s get PAW configuration data for Z=%g\n", label, Z_core);
            auto const ec = sigma_config::get(Z_core, echo);
            if (echo > 0) printf("# %s got PAW configuration data for Z=%g: rcut=%g sigma=%g %s\n", label, ec.Z, ec.rcut*Ang, ec.sigma*Ang, _Ang);

            if (ec.Z != Z_core) warn("%s number of protons adjusted from %g to %g", label, Z_core, ec.Z);
            Z_core = ec.Z;
            sigma = std::abs(ec.sigma);
            r_cut = std::abs(ec.rcut);
            // get numax from nn[]
            numax = -1;
            for(int ell = 0; ell < 8; ++ell) {
                if (ec.nn[ell] > 0) {
                    nn[ell] = std::min(int(ec.nn[ell]), nn_limiter); // take a smaller numer of partial waves
                    int const numaxmin = 2*nn[ell] - 2 + ell; // this is the minimum numax that we need to get ec.nn[ell]
                    numax = std::max(int(numax), numaxmin);
                } // nn > 0
            } // ell
            for(int inl = 0; inl < 32; ++inl) {
                custom_occ[inl] = ec.occ[inl][0] + ec.occ[inl][1]; // spin polarization is neglected so far
            } // inl

        } else {
            // defaults:
            r_cut = 2.0; // Bohr
            sigma = 0.61; // Bohr, spread for projectors (Cu)
            numax = nu_max; // 3; // 3:up to f-projectors
            // get nn[] from numax
            for(int ell = 0; ell <= ELLMAX; ++ell) {
                nn[ell] = std::max(0, (numax + 2 - ell)/2); // suggest
                nn[ell] = std::min(int(nn[ell]), nn_limiter); // take a smaller numer of partial waves
            } // ell

        } // initialize from sigma_config
        
        if (echo > 0) printf("# %s projectors and partial waves are expanded up to numax = %d\n", label,  numax);
        ellmax = 2*numax; // could be smaller than 2*numax
        if (echo > 0) printf("# %s radial density and potentials are expanded up to lmax = %d\n", label, ellmax);
        ellmax_compensator = std::min(4, int(ellmax));
        if (echo > 0) printf("# %s compensation charges are expanded up to lmax = %d\n", label, ellmax_compensator);

        sigma_compensator = r_cut/std::sqrt(20.); // Bohr
        sigma_inv = 1./sigma;
        
        if (echo > 0) {
            printf("# %s numbers of projectors ", label);
            for(int ell = 0; ell <= ELLMAX; ++ell) {
                printf(" %d", nn[ell]);
            }   printf("\n");
        } // echo


        int const nlm = pow2(1 + ellmax);
        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            int const nr = (TRU == ts)? nrt : nrs;
            // spherically symmetric quantities:
            core_density[ts] = std::vector<double>(nr, 0.0); // get memory
            potential[ts]    = std::vector<double>(nr, 0.0); // get memory
            spherical_valence_density[ts] = std::vector<double>(nr, 0.0); // get memory
            // quantities with lm-resolution:
            int const mr = align<2>(nr);
            full_density[ts]   = view2D<double>(nlm, mr, 0.0); // get memory
            full_potential[ts] = view2D<double>(nlm, mr, 0.0); // get memory
        } // true and smooth

        auto const load_stat = atom_core::read_Zeff_from_file(potential[TRU].data(), *rg[TRU], Z_core, "pot/Zeff", -1, echo);
        if (load_stat) error("loading of potential file failed for Z=%g", Z_core);
          

//         // show the loaded Zeff(r)
//         if (echo > 0) {
//            printf("\n## loaded Z_eff(r) function:\n");
//            for(int ir = 0; ir < rg[TRU]->n; ++ir) {
//                printf("%.15g %.15g\n", rg[TRU]->r[ir], -potential[TRU][ir]);
//            } // ir
//         } // echo

        std::vector<int8_t> as_valence(96, -1);
        int8_t enn_core_ell[12] = {0,0,0,0, 0,0,0,0, 0,0,0,0}; // enn-QN of the highest occupied core level

        int const inl_core_hole = control::get("core.hole.index", -1.);;
        double const core_hole_charge = std::min(std::max(0.0, control::get("core.hole.charge", 1.)), 1.0);
        double core_hole_charge_used{0};
        

        set(csv_charge, 3, 0.0); // clear numbers of electrons for core, semicore and valence, respectively

        double const core_valence_separation  = control::get("core.valence.separation", -2.0);
        double       core_semicore_separation = control::get("core.semicore.separation",    core_valence_separation);
        double const semi_valence_separation  = control::get("semicore.valence.separation", core_valence_separation);
        if (core_semicore_separation > semi_valence_separation) {
            warn("%s core.semicore.separation=%g may not be higher than semicore.valence.separation=%g, correct", 
                 label, core_semicore_separation, semi_valence_separation);
            core_semicore_separation = semi_valence_separation; // correct --> no semicore states possible
        } // warning
        ncorestates = 20; // maximum number
        true_core_waves = view3D<double>(2, ncorestates, nrt, 0.0); // get memory for the true radial wave functions and kinetic waves
        int constexpr SRA = 1;
        core_state = std::vector<core_level_t>(ncorestates);
        {
            for(int ics = 0; ics < ncorestates; ++ics) {
                core_state[ics].c0s1v2 = -1; // -1:undefined
            } // ics
            
            int ics = 0, highest_occupied_core_state_index = -1; 
            double n_electrons = Z_core - ionization; // init number of electrons to be distributed
            for(int nq_aux = 0; nq_aux < 8; ++nq_aux) { // auxiliary quantum number, allows Z up to 120
                int enn = (nq_aux + 1)/2; // init principal quantum number n
                for(int ell = nq_aux/2; ell >= 0; --ell) { // angular momentum character l
                    ++enn; // update principal quantum number
                    for(int jj = 2*ell; jj >= 2*ell; jj -= 2) { // total angular momentum j
                        auto &cs = core_state[ics]; // abbreviate
                        cs.wave[TRU] = true_core_waves(0,ics); // the true radial function
                        cs.wKin[TRU] = true_core_waves(1,ics); // the kinetic energy wave
                        
                        std::vector<double> r2rho(nrt, 0.0);
                        double E = atom_core::guess_energy(Z_core, enn); // init with a guess
                        // solve the eigenvalue problem with a spherical potential
                        radial_eigensolver::shooting_method(SRA, *rg[TRU], potential[TRU].data(),
                                                    enn, ell, E, cs.wave[TRU], r2rho.data());
                        cs.energy = E; // store eigenenergy of the spherical potential

                        int const inl = atom_core::nl_index(enn, ell);
                        int c0s1v2_auto{-1}; // -1:undefined
                        if (E > semi_valence_separation) {
                            c0s1v2_auto = 2; // mark as valence state
//                          printf("# as_valence[nl_index(enn=%d, ell=%d) = %d] = %d\n", enn, ell, inl, ics);
                        } else // move to the valence band
                        if (E > core_semicore_separation) {
                            c0s1v2_auto = 1; // mark as semicore state
                        } else { // move to the semicore band
                            c0s1v2_auto = 0; // mark as core state
                        } // stay in the core
                        assert(c0s1v2_auto >= 0);

                        cs.nrn[TRU] = enn - ell - 1; // true number of radial nodes
                        cs.enn = enn;
                        cs.ell = ell;
                        cs.emm = emm_Degenerate;
                        cs.spin = spin_Degenerate;
                        double const max_occ = 2*(jj + 1);
                        double const core_hole = core_hole_charge*(inl_core_hole == inl);
                        double const occ_no_core_hole = std::min(std::max(0., n_electrons), max_occ);
                        double const occ_auto = std::min(std::max(0., occ_no_core_hole - core_hole), max_occ) ;
                        double const real_core_hole_charge = occ_no_core_hole - occ_auto;
                        if (real_core_hole_charge > 0) {
                            core_hole_charge_used = real_core_hole_charge;
                            if (echo > 1) printf("# %s use a core.hole.index=%i (%d%c) missing core.hole.charge=%g electrons\n", 
                                                    label, inl, enn, ellchar[ell], real_core_hole_charge);
                        } // core hole active
                        
                        int const c0s1v2_cust = 2*(custom_occ[inl] > 0); // in custom_config, negative occupation numbers for core electrons, positive for valence
                        cs.c0s1v2 = custom_config ? c0s1v2_cust : c0s1v2_auto;

                        if (2 == cs.c0s1v2) as_valence[inl] = ics; // mark as good for the valence band, store the core state index

                        double const occ = custom_config ? std::abs(custom_occ[inl]) : occ_auto;
                        cs.occupation = occ;
                        if (occ > 0) {
                          
                            if ((c0s1v2_cust != c0s1v2_auto) && (echo > 4)) printf("# %s energy selectors suggest that %d%c is a %s-state but custom configuration says %s, use %s\n",
                                    label, enn, ellchar[ell], c0s1v2_name[c0s1v2_auto], c0s1v2_name[c0s1v2_cust], c0s1v2_name[cs.c0s1v2]);
                            // ToDo: find out if it would be better to look at the q_out (charge outside rcut) for each state instead of its energy.
                          
                            highest_occupied_core_state_index = ics; // store the index of the highest occupied core state
                            if (echo > 0) printf("# %s %-9s %2d%c%6.1f E=%16.6f %s\n",
                                                    label, c0s1v2_name[cs.c0s1v2], enn, ellchar[ell], occ, cs.energy*eV,_eV);
                            if (as_valence[inl] < 0) {
                                enn_core_ell[ell] = std::max(enn, int(enn_core_ell[ell])); // find the largest enn-quantum number of the occupied core states
                            } // not as valence
                            csv_charge[cs.c0s1v2] += occ;

                            // can we move this to the core state routine?
                            double const has_norm = dot_product(rg[TRU]->n, r2rho.data(), rg[TRU]->dr);
                            if (has_norm > 0) {
                                double const norm = occ/has_norm;
                                auto const density = (2 == cs.c0s1v2) ? spherical_valence_density[TRU].data() : core_density[TRU].data();
                                add_product(density, rg[TRU]->n, r2rho.data(), norm);
                            } else {
                                warn("%s %i%c-state cannot be normalized! integral=%g", label, enn, ellchar[ell], has_norm);
                            } // can be normalized
                            
                        } // occupied

                        n_electrons -= occ; // subtract as many electrons as have been assigned to this orbital
                        ++ics;
                    } // jj
                } // ell
            } // nq_aux
            ncorestates = highest_occupied_core_state_index + 1; // correct the number of core states to those occupied
            core_state.resize(ncorestates); // destruct unused core states
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
        
        scale(core_density[TRU].data(), rg[TRU]->n, rg[TRU]->rinv); // initial_density produces r^2*rho --> reduce to r*rho
        scale(core_density[TRU].data(), rg[TRU]->n, rg[TRU]->rinv); // initial_density produces   r*rho --> reduce to   rho
        if (echo > 2) printf("# %s initial core density has %g electrons\n", 
                                label, dot_product(rg[TRU]->n, core_density[TRU].data(), rg[TRU]->r2dr));

        scale(spherical_valence_density[TRU].data(), rg[TRU]->n, rg[TRU]->rinv); // initial_density produces r^2*rho --> reduce to r*rho
        scale(spherical_valence_density[TRU].data(), rg[TRU]->n, rg[TRU]->rinv); // initial_density produces   r*rho --> reduce to   rho
        if (echo > 2) printf("# %s initial valence density has %g electrons\n", 
                                label, dot_product(rg[TRU]->n, spherical_valence_density[TRU].data(), rg[TRU]->r2dr));

        if (echo > 5) printf("# %s enn_core_ell  %i %i %i %i\n", label, enn_core_ell[0], enn_core_ell[1], enn_core_ell[2], enn_core_ell[3]);

        int const nln = sho_tools::nSHO_radial(numax); // == (numax*(numax + 4) + 4)/4
        partial_waves[TRU] = view3D<double>(2, nln, nrt, 0.0); // get memory for the true   radial wave function and kinetic wave
        partial_waves[SMT] = view3D<double>(2, nln, nrs, 0.0); // get memory for the smooth radial wave function and kinetic wave
        
        set(partial_wave_char, 32, '\0'); 
        double constexpr energy_derivative = -8.0;
        double const excited_energy = control::get("single_atom.partial.wave.energy", 1.0); // 1.0 Hartree higher, special function for -8.0:energy derivative
        set(partial_wave_energy_split, 8, excited_energy);
        
        double const energy_parameter = control::get("single_atom.partial.wave.energy.parameter", -9.9e9);
        bool const use_energy_parameter = (energy_parameter > -9e9);
        if (use_energy_parameter) {
            if (echo > 5) printf("# %s use energy parameter %g %s for all ell-channels\n", label, energy_parameter*eV, _eV);
        } // use_energy_parameter

        partial_wave = std::vector<valence_level_t>(nln);
        {
//          if (echo > 0) printf("# valence "); // no new line, compact list follows
            for(int ell = 0; ell <= numax; ++ell) {
//              for(int nrn = 0; nrn < nn[ell]; ++nrn) { // smooth number or radial nodes
                for(int nrn = 0; nrn <= (numax - ell)/2; ++nrn) { // smooth number or radial nodes

                    int const iln = sho_tools::ln_index(numax, ell, nrn);
                    assert(iln < nln);
                    auto &vs = partial_wave[iln]; // abbreviate
                    int const enn = std::max(ell + 1, int(enn_core_ell[ell]) + 1) + nrn;
//                  if (echo > 0) printf(" %d%c", enn, ellchar[ell]);
                    vs.nrn[TRU] = enn - ell - 1; // true number of radial nodes
                    vs.nrn[SMT] = nrn; // number of radial nodes in the smooth wave
                    vs.occupation = 0;
                    vs.enn = enn;
                    vs.ell = ell;
                    vs.emm = emm_Degenerate;
                    vs.spin = spin_Degenerate;
                    vs.c0s1v2 = 2; // 2:valence state
                    

                    vs.wave[SMT] = partial_waves[SMT](0,iln); // the smooth radial function
                    vs.wave[TRU] = partial_waves[TRU](0,iln); // the true radial function
                    vs.wKin[SMT] = partial_waves[SMT](1,iln); // the smooth kinetic energy
                    vs.wKin[TRU] = partial_waves[TRU](1,iln); // the true kinetic energy

                    int const inl = atom_core::nl_index(enn, ell);
                    double occ{custom_occ[inl]};
                    if (nrn < nn[ell]) {
                        double E = std::max(atom_core::guess_energy(Z_core, enn), core_valence_separation);
                        if (nrn > 0) E = std::max(E, partial_wave[iln - 1].energy); // higher than previous energy
                        radial_eigensolver::shooting_method(SRA, *rg[TRU], potential[TRU].data(), enn, ell, E, vs.wave[TRU]);
                        vs.energy = E;


                        {
                            int const ics = as_valence[inl]; // index of the corresponding spherical state
    //                      printf("# as_valence[nl_index(enn=%d, ell=%d) = %d] = %d\n", enn, ell, inl, ics);
                            if (ics >= 0) { // atomic eigenstate was marked as valence
                                occ = core_state[ics].occupation; //
                                vs.occupation = occ;
                                if (occ > 0 && echo > 3) printf("# %s transfer %.1f electrons from %d%c-core state #%d"
                                        " to valence state #%d\n", label, occ, enn, ellchar[ell], ics, iln);
                            } // ics
                        }
                        if (echo > 0) printf("# %s valence %2d%c%6.1f E=%16.6f %s\n", label, enn, ellchar[ell], vs.occupation, E*eV,_eV);

                        partial_wave_char[iln] = '0' + enn; // eigenstate: '1', '2', '3', ...
                        
                        if (use_energy_parameter) {
                            vs.energy = energy_parameter;
                            partial_wave_char[iln] = 'e';
                        } else
                        if (nrn > 0) {
                            partial_wave_char[iln] = (energy_derivative == partial_wave_energy_split[ell]) ? 'd' : '*';
                        } else
                        if (occ <= 0) {
                            partial_wave_char[iln] = 'p'; // use base energy of previous ell-channel (polarization)
                        }

                    } else { // nrn < nn[ell]
                        partial_wave_char[iln] = '_'; // 'd':energy derivative, 'e': energy parameter, '1', '2', ... for eigenstates, '_':not used
                    }
                    

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
            assert(0 == int(stat));
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


//         {   // construct initial smooth spherical potential
//             set(potential[SMT].data(), rg[SMT]->n, potential[TRU].data() + nr_diff); // copy the tail of the spherical part of r*V_tru(r)
//             if (echo > 2) printf("\n# %s construct initial smooth spherical potential as parabola\n", label);
//             pseudize_function(potential[SMT].data(), rg[SMT], ir_cut[SMT], 2, 1); // replace by a parabola
//         }
        pseudize_local_potential<1>(potential[SMT].data(), potential[TRU].data(), echo); // <1> indicates that these arrays hold r*V(r) functions
        
        {   // construct an initial valence density
            spherical_valence_charge_deficit = pseudize_spherical_density(
                spherical_valence_density[SMT].data(), 
                spherical_valence_density[TRU].data(), "spherical valence", echo);
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
            ", r^2*rho_tot_tru(r), r^2*rho_tot_smt(r)"
            ", r^2*rho_cor_tru(r), r^2*rho_cor_smt(r)"
            ", r^2*rho_val_tru(r), r^2*rho_val_smt(r)"
            ", zero_potential(r) in Ha"
            ":\n", label);
            for(int irs = 0; irs < rg[SMT]->n; irs += 1) {
                int const irt = irs + nr_diff;
                auto const r = rg[SMT]->r[irs], r2 = pow2(r);
                printf("%g %g %g %g %g %g %g %g %g %g\n", r
//                         , -full_potential[TRU][00][irt]*Y00*r // for comparison, should be the same as Z_eff(r)
                        , -potential[TRU][irt] // Z_eff(r)
                        , -potential[SMT][irs] // \tilde Z_eff(r)
                        , (core_density[TRU][irt] + spherical_valence_density[TRU][irt])*r2*Y00*Y00
                        , (core_density[SMT][irs] + spherical_valence_density[SMT][irs])*r2*Y00*Y00
                        , core_density[TRU][irt]*r2*Y00*Y00
                        , core_density[SMT][irs]*r2*Y00*Y00
                        , spherical_valence_density[TRU][irt]*r2*Y00*Y00
                        , spherical_valence_density[SMT][irs]*r2*Y00*Y00
                        , zero_potential[irs]*Y00 // not converted to eV
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
      gaunt_init = (0 == int(stat));
      return stat;
    } // initialize_Gaunt

    void show_state_analysis(int const echo, radial_grid_t const *rg, double const wave[],
            int const enn, int const ell, float const occ, double const energy, int8_t const c0s1v2, int const ir_cut=-1) const {
        if (echo < 1) return; // this function only prints to the logg
        
        double q{0}, qr{0}, qr2{0}, qrm1{0}, qout{0};
        for(int ir = 0; ir < rg->n; ++ir) {
            double const rho_wf = pow2(wave[ir]);
            double const dV = rg->r2dr[ir], r = rg->r[ir], r_inv_dV = rg->rdr[ir];
            q    += rho_wf*dV; // charge
            qr   += rho_wf*r*dV; // for <r>
            qr2  += rho_wf*r*r*dV; // for variance
            qrm1 += rho_wf*r_inv_dV; // Coulomb integral without -Z
            qout += rho_wf*dV*(ir >= ir_cut);
        } // ir
        double const qinv = (q > 0) ? 1./q : 0;

        printf("# %s %-9s %2d%c%6.1f E=%16.6f %s  <r>=%g rms=%g %s <r^-1>=%g %s q_out=%.3g e\n", label,
               c0s1v2_name[c0s1v2], enn, ellchar[ell], std::abs(occ), energy*eV,_eV, 
               qr*qinv*Ang, std::sqrt(std::max(0., qr2*qinv))*Ang,_Ang, qrm1*qinv*eV,_eV, qout*qinv);
    } // show_state_analysis

    double pseudize_spherical_density(double smooth_density[], double const true_density[]
                                    , char const *quantity="core", int const echo=0) const {
        int const nrs = rg[SMT]->n;
        set(smooth_density, nrs, true_density + nr_diff); // copy the tail of the true density into the smooth density

        auto const stat = pseudize_function(smooth_density, rg[SMT], ir_cut[SMT], 3); // 3: use r^0, r^2 and r^4
        // alternatively, pseudize_function(smooth_density, rg[SMT], ir_cut[SMT], 3, 2); // 3, 2: use r^2, r^4 and r^6
        if (stat) warn("%s Matching procedure for the smooth %s density failed! info = %d", label, quantity, int(stat));

        if (echo > 8) { // plot the densities
            printf("\n## %s radius, {smooth, true} for %s density:\n", label, quantity);
            for(int ir = 0; ir < nrs; ir += 2) { // plot only every second entry
                printf("%g %g %g\n", rg[SMT]->r[ir], smooth_density[ir], true_density[ir + nr_diff]);
            }   printf("\n\n");
        } // plot
        
        // report integrals
        auto const tru_charge = dot_product(rg[TRU]->n, rg[TRU]->r2dr, true_density);
        auto const smt_charge = dot_product(rg[SMT]->n, rg[SMT]->r2dr, smooth_density);
        double const charge_deficit = tru_charge - smt_charge;
        if (echo > 1) printf("# %s true and smooth %s density have %g and %g electrons\n", label, quantity, tru_charge, smt_charge);

        return charge_deficit;
    } // pseudize_spherical_density
    
    void update_core_states(float const mixing, int const echo=0) {
        if (echo > 1) printf("\n# %s %s Z=%g\n", label, __func__, Z_core);
        // core states are feeling the spherical part of the hamiltonian only
        int const nr = rg[TRU]->n;
        std::vector<double> r2rho(nr);
        std::vector<double> new_r2core_density(nr, 0.0);
        std::vector<double> new_r2valence_density(nr, 0.0);
        double nelectrons[] = {0, 0}; // core, valence
#ifdef DEVEL
        if (echo > 7) {
            printf("# %s %s: solve for eigenstates of the full radially symmetric potential\n", label, __func__);
            printf("\n## %s %s: r, -Zeff(r)\n", label, __func__);
            print_compressed(rg[TRU]->r, potential[TRU].data(), rg[TRU]->n);
        } // echo
#endif
        for(int ics = 0; ics < ncorestates; ++ics) {
            auto & cs = core_state[ics]; // abbreviate
            int constexpr SRA = 1; // 1:scalar relativistic approximation
            radial_eigensolver::shooting_method(SRA, *rg[TRU], potential[TRU].data(), cs.enn, cs.ell, cs.energy, cs.wave[TRU], r2rho.data());
            auto const norm = dot_product(nr, r2rho.data(), rg[TRU]->dr);
            auto const norm_factor = (norm > 0)? 1./std::sqrt(norm) : 0;
            auto const scal = pow2(norm_factor)*std::abs(cs.occupation); // scaling factor for the density contribution of this state
            // transform r*wave(r) as produced by the radial_eigensolver to wave(r)
            // and normalize the core level wave function to one
            scale(cs.wave[TRU], nr, rg[TRU]->rinv, norm_factor);
            // create wKin for the computation of the kinetic energy density
            product(cs.wKin[TRU], nr, potential[TRU].data(), cs.wave[TRU], -1.); // start as wKin = -r*V(r)*wave(r)
            add_product(cs.wKin[TRU], nr, rg[TRU]->r, cs.wave[TRU], cs.energy); // now wKin = r*(E - V(r))*wave

            if (scal > 0) {
                if (0 == cs.c0s1v2) { // core state
                    add_product(new_r2core_density.data(), nr, r2rho.data(), scal);
                    nelectrons[0] += std::max(0.0, std::abs(cs.occupation));
                } else if (2 == cs.c0s1v2) { // valence state
                    add_product(new_r2valence_density.data(), nr, r2rho.data(), scal);
                    nelectrons[1] += std::max(0.0, std::abs(cs.occupation));
                } else error("%s spherical semicore density not implemented!", label); 
            } // scal > 0
            show_state_analysis(echo, rg[TRU], cs.wave[TRU], cs.enn, cs.ell, cs.occupation, cs.energy, cs.c0s1v2, ir_cut[TRU]);
        } // ics

        // report integrals
        double old_charge[2], new_charge[2];
        old_charge[0] = dot_product(nr, rg[TRU]->r2dr, core_density[TRU].data());
        new_charge[0] = dot_product(nr, rg[TRU]->dr,  new_r2core_density.data());
        if (echo > 2) printf("# %s previous core density has %g electrons, expected %g\n"
                             "# %s new core density has %g electrons\n", label, 
                             old_charge[0], nelectrons[0], label, new_charge[0]);
        
        // spherical valence density is not mixed, changes are 100% accepted
        old_charge[1] = dot_product(nr, rg[TRU]->r2dr, spherical_valence_density[TRU].data());
        new_charge[1] = dot_product(nr, rg[TRU]->dr,  new_r2valence_density.data());
        if (echo > 4) printf("# %s previous spherical valence density has %g electrons, expected %g\n"
                             "# %s new spherical valence density has %g electrons\n", label, 
                             old_charge[1], nelectrons[1], label, new_charge[1]);
        product(spherical_valence_density[TRU].data(), nr, new_r2valence_density.data(), rg[TRU]->rinv);
        scale  (spherical_valence_density[TRU].data(), nr, rg[TRU]->rinv);
        // check again
        if (mix_spherical_valence_density != 0) {
            auto const new_q = dot_product(nr, rg[TRU]->r2dr,  spherical_valence_density[TRU].data());
            if (echo > 4) printf("# %s new spherical valence density has %g electrons\n", label, new_q);
        } // use spherical_valence_density
        
        double mix_new = mixing, mix_old = 1 - mix_new;
        // rescale the mixing coefficients such that the desired number of core electrons comes out
        auto const mixed_charge = mix_old*old_charge[0] + mix_new*new_charge[0];
        if (mixed_charge != 0) {
            auto const rescale = nelectrons[0]/mixed_charge;
            mix_old *= rescale;
            mix_new *= rescale;
        } // rescale

        double core_density_change{0}, core_density_change2{0}, core_nuclear_energy{0};
        for(int ir = 0; ir < nr; ++ir) {
            double const r2inv = pow2(rg[TRU]->rinv[ir]);
            auto const new_rho = new_r2core_density[ir]*r2inv; // *r^{-2}
            spherical_valence_density[TRU][ir] = new_r2valence_density[ir]*r2inv; // *r^{-2}
            core_density_change  += std::abs(new_rho - core_density[TRU][ir])*rg[TRU]->r2dr[ir];
            core_density_change2 +=     pow2(new_rho - core_density[TRU][ir])*rg[TRU]->r2dr[ir];
            core_nuclear_energy  +=         (new_rho - core_density[TRU][ir])*rg[TRU]->rdr[ir]; // Coulomb integral change
            core_density[TRU][ir] = mix_new*new_rho + mix_old*core_density[TRU][ir];
        } // ir
        core_nuclear_energy *= -Z_core;
        if (echo > 0) printf("# %s core density change %g e (rms %g e) energy change %g %s\n", label,
            core_density_change, std::sqrt(std::max(0.0, core_density_change2)), core_nuclear_energy*eV,_eV);

        core_charge_deficit = pseudize_spherical_density(core_density[SMT].data(), 
                              core_density[TRU].data(), "core", echo - 1); 
        spherical_valence_charge_deficit = pseudize_spherical_density(spherical_valence_density[SMT].data(),
                              spherical_valence_density[TRU].data(), "spherical valence", echo - 1); 
    } // update_core_states

    double update_sigma(int const echo=0) {
        if (echo > 1) printf("\n# %s %s Z=%g\n", label, __func__, Z_core);
        
        double const sigma_old = sigma;

        int const nln = sho_tools::nSHO_radial(numax); // == (numax*(numax + 4) + 4)/4
        view2D<double> tphi(nln, align<2>(rg[TRU]->n), 0.0); // get memory for r*phi_tru
        view2D<double> sphi(nln, align<2>(rg[SMT]->n), 0.0); // get memory for r*phi_smt
        view2D<double> rprj(nln, align<2>(rg[SMT]->n), 0.0); // get memory for r*projector
        view2D<double> skin(nln, align<2>(rg[SMT]->n), 0.0); // get memory for r*T*phi_smt
        view2D<double> prj_sho(nln, align<2>(rg[SMT]->n), 0.0); // get memory for r*projector
        
        double occ_ell[12]; set(occ_ell, 12, 0.0);
        double total_occ{0};
        
        for(int ell = 0; ell <= numax; ++ell) {
            for(int nrn = 0; nrn < nn[ell]; ++nrn) { // smooth number or radial nodes
                int const iln = sho_tools::ln_index(numax, ell, nrn);
                auto & vs = partial_wave[iln]; // abbreviate ('vs' stands for valence state)

                int nnodes{0}; int constexpr SRA = 1;
                // integrate outwards homogeneously, WARNING: energy derivative not supported
                // (^T + V - E) phi(r) == 0
                radial_integrator::shoot(SRA, *rg[TRU], potential[TRU].data(), ell, vs.energy, nnodes, tphi[iln]);
                
                // copy the true tail of the TRU wave function into the SMT wave function
                set(sphi[iln], rg[SMT]->n, tphi[iln] + nr_diff);
                
                // pseudize by matching a polynomial r^ell*(c_0 + c_2 r^2 + c_4 r^4 + c_6 r^6)
                double coeff[4]; // matching coefficients
                auto const stat = pseudize_function(sphi[iln], rg[SMT], ir_cut[SMT], 4, ell + 1, coeff);
                assert(0 == stat);
                
                // construct a preliminary projector according to the Bloechl scheme:
                //  p(r) = (^T + ~V - E) ~phi(r)
                std::vector<double> rkin(rg[SMT]->n); // get memory for the kinetic wave
                double ckin[3]; // coefficients of the kinetic wave
                // ^T = .5*( ell(ell+1)/r^2 - d^2/dr^2 ) when acting onto r*wave
                for(int k = 1; k < 4; ++k) {
                    ckin[k - 1] = 0.5*( ell*(ell + 1) - (ell + 2*k)*(ell + 1 + 2*k) )*coeff[k];
                } // k
                
                if (echo > 19) printf("\n## %s classical method for %c%i: r, smooth r*wave, smooth r*Twave, true r*wave, projector:\n", label, ellchar[ell], nrn);
                // expand kinetic wave up to r_cut
                for(int ir = 0; ir <= ir_cut[SMT]; ++ir) {
                    double const r = rg[SMT]->r[ir], r2 = pow2(r), rl1 = intpow(r, ell + 1);
                    rkin[ir] = rl1*(ckin[0] + r2*(ckin[1] + r2*ckin[2])); // expand polynomial
                    // construct the preliminary projector functions accoding to Bloechls scheme
                    rprj[iln][ir] = rkin[ir] + (potential[SMT][ir]*rg[SMT]->rinv[ir] - vs.energy)*sphi[iln][ir];
                    if (echo > 19) printf("%g %g %g %g %g\n", r, sphi[iln][ir], rkin[ir], tphi[iln][ir + nr_diff], rprj[iln][ir]*rg[SMT]->rinv[ir]);
                } // ir
                if (echo > 19) printf("\n## %s classical method for %c%i: r, smooth r*wave, smooth r*Twave, true r*wave:\n", label, ellchar[ell], nrn);
                // beyond the cutoff radius, ^T sphi should math ^T tphi
                for(int ir = ir_cut[SMT]; ir < rg[SMT]->n; ++ir) {
                    int const ir_tru = ir + nr_diff;
                    rkin[ir] = (vs.energy - potential[TRU][ir_tru]*rg[SMT]->rinv[ir])*tphi[iln][ir_tru];
                    if (echo > 19) printf("%g %g %g %g\n", rg[SMT]->r[ir], sphi[iln][ir], rkin[ir], tphi[iln][ir + nr_diff]);
                } // ir
                if (echo > 19) printf("\n\n");
                
                set(skin[iln], rg[SMT]->n, rkin.data()); // store for later usage (sugestion of a potential shape)
                
                // check that the two functions match smoothly at r_cut, seems ok
                //
                // this technique produces projectors that go to zero smoothly at r_cut, however, they feature high frequency components
                //
                // now expand the projectors in SHO projectors and minimize the deviation by optimizing sigma.
                // the total deviation is weighted with the ell-channel-summed occupation numbers, occ_ell.
                occ_ell[ell] += vs.occupation;
            } // nrn
            total_occ += occ_ell[ell];
        } // ell

        double sig = 0.75*sigma_old; // initialize, start lower than the suggested sigma

        int const scan_sigma = control::get("single_atom.scan.sigma", 199.);
        std::vector<double> best_projector_coeff(nln, 0.0);
        double best_weighted_quality{-1}, sigma_opt{-1};
        int best_weighted_quality_at_iter{-1};
        // to optimize efficiently we need to define a derivative w.r.t. sigma and use a Newton method or bisection to find its zero.
        for(int iter = 0; iter < scan_sigma; ++iter) {
            int const nr = ir_cut[SMT];

            // suggest a new sigma (named sig here)
            sig *= 1.01;
            
            // compute the value of the missmatch
            
            scattering_test::expand_sho_projectors(prj_sho.data(), prj_sho.stride(), *rg[SMT], sig, numax, 0, echo/2);
            
            std::vector<double> projector_coeff(nln, 0.0);
            
            double weighted_quality{0};
            for(int ell = 0; ell <= numax; ++ell) {
                double quality_ell{0};
                for(int nrn = 0; nrn < nn[ell]; ++nrn) { // smooth number or radial nodes
                    int const iln = sho_tools::ln_index(numax, ell, nrn);
                    double const denom_sho = dot_product(nr, prj_sho[iln], prj_sho[iln], rg[SMT]->r2dr); // should be close to 1.0
                    double const denom_num = dot_product(nr, rprj[iln], rprj[iln], rg[SMT]->dr); // norm^2 of numerically given projectors
//                  if (echo > 1) printf("# %s in iteration %i norm of %c%i projectors: sho %g classical %g\n", label, iter, ellchar[ell], nrn, denom_sho, denom_num);
                    if (denom_sho*denom_num > 0) {
                        double const inner = dot_product(nr, rprj[iln], prj_sho[iln], rg[SMT]->rdr);
                        double const quality = pow2(inner) / (denom_sho*denom_num);
                        if (echo > 13) printf("# %s in iteration %i quality for %c%i with sigma= %g %s is %g\n", label, iter, ellchar[ell], nrn, sig*Ang, _Ang, quality);
                        quality_ell += quality/nn[ell];
                        projector_coeff[iln] = inner / std::sqrt(denom_sho); // if nn[ell] > 1, we want to know the ratio of c1/c0
                    } else {
                        if (echo > 1) printf("# %s in iteration %i failed to normalize %c%i projectors: sho %g classical %g\n", label, iter, ellchar[ell], nrn, denom_sho, denom_num);
                    }
                } // nrn
                if (echo > 11) printf("# %s in iteration %i quality for %c-channel with sigma= %g %s is %g\n", label, iter, ellchar[ell], sig*Ang, _Ang, quality_ell);
                weighted_quality += occ_ell[ell]*quality_ell;
            } // ell
            if (echo > 9) printf("# %s in iteration %i weighted quality with sigma= %g %s is %g\n", label, iter, sig*Ang, _Ang, weighted_quality);
            
            if (echo > 9) printf("\n# %s in iteration %i suggest new sigma= %g %s\n", label, iter, sig*Ang, _Ang);
            if (weighted_quality > best_weighted_quality) {
                best_weighted_quality = weighted_quality;
                best_weighted_quality_at_iter = iter;
                set(best_projector_coeff.data(), nln, projector_coeff.data());
                sigma_opt = sig; // store
            } // store the best result

        } // iterate for optimization
        // the upper limit for the weighted quality is the number of valence electrons
        if (echo > 5) printf("\n# %s optimized sigma= %g %s with quality %g of max. %g, %.1f %%\n\n", label, sigma_opt*Ang, _Ang, 
                                  best_weighted_quality, total_occ, best_weighted_quality*100/std::max(1., total_occ));
        if (sigma_opt == sig || 0 == best_weighted_quality_at_iter) {
            warn("%s optimal sigma is at the end of the analyzed range!", label);
        } // border
        
        
        int const suggest_vloc = int(control::get("single_atom.suggest.local.potential", -1.));
        if (suggest_vloc > -1) {
            // if we dictate the shape of the projectors
            scattering_test::expand_sho_projectors(prj_sho.data(), prj_sho.stride(), *rg[SMT], sigma_opt, numax, 1, echo/2); // r*projector
            // and we want the shape of the true partial wave given by sphi[ell nrn=0]
            int const ell = suggest_vloc;
            assert(ell <= numax);
            int const iln = sho_tools::ln_index(numax, ell, 0); // nrn=0
            // then the local potential can be found from the inversion:
            //  ~V(r) = E + ( ~p(r) - ^T ~phi(r) )/~phi(r)  
            // where we need to respect the duality: < ~p | ~phi > = 1 which fixes the scaling
            double const duality = 1./dot_product(rg[SMT]->n, prj_sho[iln], sphi[iln], rg[SMT]->dr);
            if (echo > 0) printf("\n## %s suggested local potential for ell=%i (%c) and currently used potential (Bohr, Ha, Ha):\n", label, ell, ellchar[ell]);
            for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                double const denom = (std::abs(sphi[iln][ir]) < 1e-6) ? 1 : sphi[iln][ir];
                double const pot = partial_wave[iln].energy + ( prj_sho[iln][ir]*duality - skin[iln][ir] ) / denom;
                if (echo > 0) printf("%g %g %g\n", rg[SMT]->r[ir],  pot, potential[SMT][ir]*rg[SMT]->rinv[ir]);
            } // ir
            if (echo > 0) printf("\n\n");
            error("terminate after task single_atom.suggest.local.potential=%i", suggest_vloc);
        } // suggest smooth local potential shape
        
        if (echo > 3) printf("# %s optimize sigma from %g to %g %s\n", label, sigma_old*Ang, sigma_opt*Ang, _Ang);
        return sigma_opt; // use the optimized sigma from here on
        
    } // update_sigma
        
    template <int rpow>
    status_t pseudize_local_potential(double V_smt[], double const V_tru[], int const echo=0, double const df=1) { // df=display_factor
        if (echo > 1) printf("\n# %s %s Z=%g\n", label, __func__, Z_core);
        status_t stat(0);
        double const r_cut = rg[TRU]->r[ir_cut[TRU]];
        auto const method = control::get("single_atom.local.potential.method", "parabola");
        if ('p' == *method) { // construct initial smooth spherical potential
            set(V_smt, rg[SMT]->n, V_tru + nr_diff); // copy the tail of the spherical part of V_tru(r) or r*V_tru(r)
            if (echo > 2) printf("\n# %s construct initial smooth spherical potential as parabola\n", label);
            stat = pseudize_function(V_smt, rg[SMT], ir_cut[SMT], 2, rpow); // replace by a parabola
            if (echo > 5) printf("# %s match local potential to parabola at R_cut = %g %s, V_tru(R_cut) = %g %s\n",
                           label, r_cut*Ang, _Ang, V_tru[ir_cut[TRU]]*(rpow ? 1./r_cut : 1)*df*eV, _eV);
        } else {
            assert(rpow == rpow*rpow); // rpow either 0 or 1
            // construct a Lagrange polynomial of controllable order to fit r*V_true(r) around r_cut
            int const ir0 = ir_cut[TRU];
            double const x0 = rg[TRU]->r[ir0];
            double xi[32], yi[32];
            yi[0] = V_tru[ir0]*(rpow ? 1 : x0); // r_cut*V(r_cut)
            xi[0] = rg[TRU]->r[ir0] - x0;  // expand around r_cut
            for (int order = 1; order < 16; ++order) {
                {
                    int const ir = ir0 + order; assert(ir < rg[TRU]->n);
                    int const i = 2*order - 1;
                    yi[i] = V_tru[ir]*(rpow ? 1 : rg[TRU]->r[ir]);
                    xi[i] = rg[TRU]->r[ir] - x0;
                }
                {
                    int const ir = ir0 - order; assert(ir >= 0);
                    int const i = 2*order;
                    yi[i] = V_tru[ir]*(rpow ? 1 : rg[TRU]->r[ir]);
                    xi[i] = rg[TRU]->r[ir] - x0;
                }
            } // order
#ifdef DEVEL
            if (echo > 22) {
                printf("\n## Lagrange-N, reference-value, Lagrange(rcut), dLagrange/dr, d^2Lagrange/dr^2, status:\n");
                for (int order = 1; order < 16; ++order) {
                    unsigned const Lagrange_order = 2*order + 1;
                    {
                        double d0{0}, d1{0}, d2{0};
                        auto const stat = Lagrange_derivatives(Lagrange_order, yi, xi, 0, &d0, &d1, &d2);
                        printf("%d %.15f %.15f %.15f %.15f %i\n", Lagrange_order, V_tru[ir0], d0, d1, d2, int(stat));
                    }
                } // order
            } // echo
            
            if (echo > 27 && 1 == rpow) {
                for (int order = 1; order < 16; ++order) {
                    unsigned const Lagrange_order = 2*order + 1;
                    printf("\n\n## r, r*V_true(r), Lagrange_fit(N=%d), dLagrange/dr, d^2Lagrange/dr^2, status:\n", Lagrange_order);
                    for(int ir = 0; ir < rg[TRU]->n; ++ir) {
                        double d0{0}, d1{0}, d2{0};
                        auto const stat = Lagrange_derivatives(Lagrange_order, yi, xi, rg[TRU]->r[ir] - x0, &d0, &d1, &d2);
                        printf("%g %g %g %g %g %i\n", rg[TRU]->r[ir], V_tru[ir], d0, d1, d2, int(stat));
                    } // ir
                    printf("\n\n");
                } // order
            } // echo
#endif        
            // Fit a V_s*sin(r*k_s) + r*V_0 to r*V(r) at r_cut, fulfil three equations:
            //  (r*V(r))  |r=r_cut  = d0 =  V_s*sin(r_cut*k_s) + r_cut*V_0
            //  (r*V(r))' |r=r_rcut = d1 =  V_s*cos(r_cut*k_s)*k_s  +  V_0
            //  (r*V(r))''|r=r_rcut = d2 = -V_s*sin(r_cut*k_s)*k_s^2
            //
            //  with d0 < 0, d1 > 0, d2 < 0 and r_cut*k_s constrained to (2*pi, 2.5*pi)
            unsigned const Lagrange_order = 1 + 2*std::min(std::max(1, int(control::get("single_atom.lagrange.derivative", 7.))), 15);
            double d0{0}, d1{0}, d2{0};
            stat = Lagrange_derivatives(Lagrange_order, yi, xi, 0, &d0, &d1, &d2);
            if (echo > 7) printf("# %s use %d points, value=%g derivative=%g second=%g status=%i\n", label, Lagrange_order, yi[0], d0, d1, d2, int(stat));
            
            if (d1 <= 0) warn("positive potential slope for sinc-fit expected but found %g", d1);
            if (d2 >  0) warn("negative potential curvature for sinc-fit expected but found %g", d2);
            
            double k_s{2.25*constants::pi/r_cut}, k_s_prev{0}; // initial guess
            double V_s, V_0;
            int const max_iter = 999;
            int iterations_needed{0};
            for (int iter = max_iter; (std::abs(k_s - k_s_prev) > 1e-15) && (iter > 0); --iter) {
                k_s_prev = k_s;
                double const sin = std::sin(k_s*r_cut);
                double const cos = std::cos(k_s*r_cut);
                V_s = d2/(-k_s*k_s*sin);
                V_0 = d1 - V_s*k_s*cos;
                double const d0_new = V_s*sin + r_cut*V_0;
                if (echo > 27) printf("# %s iter=%i use V_s=%g k_s=%g V_0=%g d0=%g d0_new-d0=%g\n", label, iter, V_s, k_s, V_0, d0, d0 - d0_new);
                k_s += 1e-2*(d0 - d0_new); // maybe needs saturation function like atan
                ++iterations_needed;
            } // while
            if (iterations_needed >= max_iter) warn("sinc-fit did not converge!");
            if (k_s*r_cut <= 2*constants::pi || k_s*r_cut >= 2.5*constants::pi) {
                warn("sinc-fit failed!, k_s*r_cut=%g not in (2pi, 2.5pi)", k_s*r_cut);
            } // out of target range

            if (echo > 5) printf("# %s match local potential to sinc function at R_cut = %g %s, V_tru(R_cut) = %g %s\n",
                           label, r_cut*Ang, _Ang, yi[0]/r_cut*df*eV, _eV);

            if (echo > 3) printf("# %s smooth potential value of sinc-fit at origin is %g %s\n", label, (V_s*k_s + V_0)*df*eV, _eV);
            
            // now modify the smooth local potential
            for(int ir = 0; ir <= ir_cut[SMT]; ++ir) {
                double const r = rg[SMT]->r[ir];
                V_smt[ir] = (V_s*sin(r*k_s) + r*V_0)*(rpow ? 1 : rg[SMT]->rinv[ir]); // set values to the fitted sinc-function
            } // ir
            for(int ir = ir_cut[SMT]; ir < rg[SMT]->n; ++ir) {
                V_smt[ir] = V_tru[ir + nr_diff]; // copy the tail
            } // ir

        } // method
#ifdef DEVEL
        if (echo > 11) {
            printf("\n\n## %s r in %s, r*V_tru(r), r*V_smt(r) in %s%s:\n", __func__, _Ang, _Ang, _eV);
            auto const factors = df*eV*Ang;
            for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                double const r = rg[SMT]->r[ir], r_pow = rpow ? 1 : r;
                printf("%g %g %g\n", r*Ang, r_pow*V_tru[ir + nr_diff]*factors, r_pow*V_smt[ir]*factors);
            } // ir
            printf("\n\n");
        } // echo
#endif
        return stat;
    } // pseudize_local_potential

    
    void update_energy_parameters(int const echo=0) {
        if (echo > 1) printf("\n# %s %s Z=%g chars=\"%s\"\n", label, __func__, Z_core, partial_wave_char);
        
        double previous_ell_energy{0};
        for(int ell = 0; ell <= numax; ++ell) { // loop runs serial, loop-carried dependency on previous_ell_energy
            for(int nrn = 0; nrn < nn[ell]; ++nrn) { // smooth number or radial nodes
                int const iln = sho_tools::ln_index(numax, ell, nrn);
                auto & vs = partial_wave[iln]; // abbreviate ('vs' stands for valence state)
                
                char const c = partial_wave_char[iln];
                assert ('_' != c);
                if ('e' == c) {
                    // energy parameter, vs.energy has already been set earlier
                    if (echo > -1) printf("# %s for ell=%i nrn=%i the partial wave is constructed at energy parameter E= %g %s\n", label, ell, nrn, vs.energy*eV, _eV);
                } else
                if ('d' == c) {
                    // energy derivative at the energy of the lower partial wave
                    assert(nrn > 0);
                    vs.energy = partial_wave[iln - 1].energy;
                    if (echo > -1) printf("# %s for ell=%i nrn=%i use an energy derivative as partial wave at E= %g %s\n", label, ell, nrn, vs.energy*eV, _eV);

                } else
                if ('*' == c) {
                    assert(nrn > 0);
                    vs.energy = partial_wave[iln - 1].energy + partial_wave_energy_split[ell];
                    if (echo > -1) printf("# %s for ell=%i nrn=%i use a partial wave at E= %g %s\n", label, ell, nrn, vs.energy*eV, _eV);
                
                } else
                if ('p' == c) {
                    if (0 == ell) warn("%s energy parameter for ell=%i nrn=%i undetermined", label, ell, nrn);
                    vs.energy = previous_ell_energy;
                    if (echo > -1) printf("# %s for ell=%i nrn=%i use a partial wave at E= %g %s\n", label, ell, nrn, vs.energy*eV, _eV);

                } else {
                    assert(c == '0' + vs.enn);
                    // find the eigenenergy of the TRU spherical potential
                    int constexpr SRA = 1;
                    radial_eigensolver::shooting_method(SRA, *rg[TRU], potential[TRU].data(), vs.enn, ell, vs.energy, vs.wave[TRU]);
                    if (echo > -1) printf("# %s for ell=%i nrn=%i use a partial wave at E=%g %s, the %i%c-eigenvalue\n",
                                          label, ell, nrn, vs.energy*eV,_eV, vs.enn, ellchar[ell]);
                }
                if (0 == nrn) previous_ell_energy = vs.energy;
            } // nrn
        } // ell
        
    } // update_energy_parameters
    
    void update_partial_waves(int const echo=0) {
        int const optimize_sigma = int(control::get("single_atom.optimize.sigma", 0.));
        if (optimize_sigma) sigma = update_sigma(echo); // change member variable to optimized sigma
        
        if (echo > 1) printf("\n# %s %s Z=%g\n", label, __func__, Z_core);
        // the basis for valence partial waves is generated from the spherical part of the hamiltonian
//      auto const small_component = new double[rg[TRU]->n];
        int const nr = rg[TRU]->n;
        std::vector<double> r2rho(nr);

        // ToDo: projectors only depend on sigma, numax and the radial grid --> move this to the constructor
        scattering_test::expand_sho_projectors(projectors.data(), projectors.stride(), *rg[SMT], sigma, numax, 1, echo/2);

        r_match = 9*sigma;
        int const ir_match[] = {radial_grid::find_grid_index(*rg[TRU], r_match),
                                radial_grid::find_grid_index(*rg[SMT], r_match)};
        if (echo > 3) printf("# %s matching radius %g %s at radial indices %i and %i\n",
                               label, r_match*Ang, _Ang, ir_match[TRU], ir_match[SMT]);

        if (echo > 9) {
            printf("\n## %s show the local potentials (r, r*Vtru, r*Vsmt):\n", label);
            for(int ir = 1; ir < rg[SMT]->n; ++ir) {
                printf("%g %g %g\n", rg[SMT]->r[ir], potential[TRU][ir + nr_diff], potential[SMT][ir]);
            } // ir
            printf("\n\n");
        } // echo

//         double const excited_energy = control::get("single_atom.partial.wave.energy", 1.0); // 1.0 Hartree higher
//         bool const use_energy_derivative = (control::get("single_atom.partial.wave.energy.derivative", 1.) > 0);
//         double const energy_parameter  = control::get("single_atom.partial.wave.energy.parameter", -9.9e9);
//         bool const use_energy_parameter = (energy_parameter > -9e9);
//         if (use_energy_parameter) {
//             if (echo > 5) printf("# %s use energy parameter %g %s for all ell-channels\n", label, energy_parameter*eV, _eV);
//         } // use_energy_parameter
        
        for(int ell = 0; ell <= numax; ++ell) {
            int const ln_off = sho_tools::ln_index(numax, ell, 0); // offset where to start indexing emm-degenerate partial waves

            if (echo > 3) printf("\n# %s %s for ell=%i\n\n", label, __func__, ell);

            view2D<double> projectors_ell(projectors[ln_off], projectors.stride()); // sub-view
            
            int const n = nn[ell];
            for(int nrn_ = 0; nrn_ < n; ++nrn_) { int const nrn = nrn_; // smooth number or radial nodes
                int const iln = ln_off + nrn;
                auto & vs = partial_wave[iln]; // abbreviate ('vs' stands for valence state)
                int constexpr SRA = 1;

                set(vs.wave[TRU], nr, 0.0); // clear

                bool normalize = true, orthogonalize = false;
                bool const use_energy_derivative = ('d' == partial_wave_char[iln]);
#if 0                
                if (0 == nrn) {
                    // solve for a true valence eigenstate
                    if (use_energy_parameter) {
                        vs.energy = energy_parameter;
                        int nnodes{0};
                        // integrate outwards homogeneously
                        radial_integrator::shoot(SRA, *rg[TRU], potential[TRU].data(), ell, vs.energy, nnodes, vs.wave[TRU], r2rho.data());
                    } else {
                        radial_eigensolver::shooting_method(SRA, *rg[TRU], potential[TRU].data(), vs.enn, ell, vs.energy, vs.wave[TRU], r2rho.data());
                        if (echo > 7) printf("# %s %s found a true %i%c-eigenstate of the spherical potential at E=%g %s\n",
                                              label, __func__, vs.enn, ellchar[ell], vs.energy*eV,_eV);
                    } // use_energy_parameter
                } else {
                    assert(nrn > 0);

                    if (use_energy_derivative) {
                        // choose Psi_1 as energy derivative:
                        // solve inhomogeneous equation (H - E) Psi_1 = Psi_0
                        vs.energy = partial_wave[iln - 1].energy; // same energy parameter as for Psi_0
                        std::vector<double> ff(rg[TRU]->n), inh(rg[TRU]->n);
                        product(inh.data(), rg[TRU]->n, partial_wave[iln - 1].wave[TRU], rg[TRU]->r);
                        double dg;
                        if (echo > -1) printf("# %s for ell=%i use energy derivative at E= %g %s\n", label, ell, vs.energy*eV, _eV);
                        radial_integrator::integrate_outwards<SRA>(*rg[TRU], potential[TRU].data(), ell, vs.energy,
                                                                   vs.wave[TRU], ff.data(), -1, &dg, inh.data());
                        // and T Psi_1 = (E - V) Psi_1 + Psi_0 for the kinetic part later
                        normalize = false;
                        orthogonalize = true;

                    } else {
                        vs.energy = (use_energy_parameter ? energy_parameter : partial_wave[iln - 1].energy) + nrn*excited_energy;
                        if (echo > -1) printf("# %s for ell=%i nrn=%i use a second partial wave at E= %g %s\n", label, ell, nrn, vs.energy*eV, _eV);
                        int nnodes{0};
                        // integrate outwards homogeneously
                        radial_integrator::shoot(SRA, *rg[TRU], potential[TRU].data(), ell, vs.energy, nnodes, vs.wave[TRU], r2rho.data());
                    } // how to construct a second tru partial wave?
                    
                } // bound state?
#else
                if (use_energy_derivative) {
                    // choose Psi_1 as energy derivative:
                    // solve inhomogeneous equation (H - E) Psi_1 = Psi_0
                    vs.energy = partial_wave[iln - 1].energy; // same energy parameter as for Psi_0
                    std::vector<double> ff(rg[TRU]->n), inh(rg[TRU]->n);
                    product(inh.data(), rg[TRU]->n, partial_wave[iln - 1].wave[TRU], rg[TRU]->r);
                    double dg;
                    if (echo > -1) printf("# %s for ell=%i use energy derivative at E= %g %s\n", label, ell, vs.energy*eV, _eV);
                    radial_integrator::integrate_outwards<SRA>(*rg[TRU], potential[TRU].data(), ell, vs.energy,
                                                                vs.wave[TRU], ff.data(), -1, &dg, inh.data());
                    // and T Psi_1 = (E - V) Psi_1 + Psi_0 for the kinetic part later
                    normalize = false;
                    orthogonalize = true;
                } else {
                    int nnodes{0};
                    // integrate outwards homogeneously
                    radial_integrator::shoot(SRA, *rg[TRU], potential[TRU].data(), ell, vs.energy, nnodes, vs.wave[TRU], r2rho.data());
                } // use_energy_derivative
#endif

                if (nrn > 0 && orthogonalize) {
                        auto const psi0 = partial_wave[iln - 1].wave[TRU];
                        double const d = dot_product(nr, vs.wave[TRU], psi0, rg[TRU]->rdr);
                        double const d2 = dot_product(nr, psi0, psi0, rg[TRU]->r2dr);
                        double const p = -d/d2; // projection part
                        if (echo > -1) printf("# %s for ell=%i orthogonalize energy derivative with %g\n", label, ell, p);
                        add_product(vs.wave[TRU], nr, psi0, rg[TRU]->r, p);
                } // orthogonalize
                
                // scope: divide by r and normalize the partial waves
                {   double norm_factor{1};
                    if (normalize) {
                        auto const norm_wf2 = dot_product(nr, r2rho.data(), rg[TRU]->dr);
                        norm_factor = 1./std::sqrt(norm_wf2);
                    } // normalize
                    scale(vs.wave[TRU], nr, rg[TRU]->rinv, norm_factor); // transform r*wave(r) as produced by the radial_eigensolver to wave(r)
                } // scope

                // create wKin for the computation of the kinetic energy density
                product(vs.wKin[TRU], nr, potential[TRU].data(), vs.wave[TRU], -1.); // start as wKin = -r*V(r)*wave(r)
                add_product(vs.wKin[TRU], nr, rg[TRU]->r, vs.wave[TRU], vs.energy); // now wKin = r*(E - V(r))*wave(r)
                if (use_energy_derivative && nrn > 0) {
                    add_product(vs.wKin[TRU], nr, rg[TRU]->r, partial_wave[iln - 1].wave[TRU]); // now wKin = r*(E - V(r))*wave(r) + r*psi0
                } // kinetic energy wave in the case of energy derivative has an extra term

//                 auto const tru_norm = dot_product(ir_cut[TRU], r2rho.data(), rg[TRU]->dr)/norm_wf2; // integrate only up to rcut
//                 auto const work = r2rho.data();
//                 scale(work, nr, potential[TRU].data());
//                 auto const tru_Epot = dot_product(ir_cut[TRU], work, rg[TRU]->dr)/norm_wf2;
//                 auto const tru_kinetic_E = vs.energy*tru_norm - tru_Epot; // kinetic energy contribution up to r_cut

//                 if (echo > 1) printf("# valence %2d%c%6.1f E=%16.6f %s\n", vs.enn, ellchar[ell], vs.occupation, vs.energy*eV,_eV);
                show_state_analysis(echo, rg[TRU], vs.wave[TRU], vs.enn, ell, vs.occupation, vs.energy, 2, ir_cut[TRU]);

//                 int const prelim_waves = control::get("single_atom.preliminary.partial.waves", 0.); // 0:none, -1:all, e.g. 5: s and d
//                 if (prelim_waves & (1 << ell)) {
//                     if (echo > -1) printf("# %s for ell=%i create a smooth partial wave by pseudization of the true partial wave\n", label, ell);
//                 } // preliminary partial waves
                

                // idea: make this module flexible enough so it can load a potential and
                //       generate PAW data in XML format (GPAW, ABINIT) using the SHO method

                { // scope: generate smooth partial waves from projectors, revPAW scheme
                    assert(ir_match[TRU] > 0);
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
                        auto projector = (krn > HOM) ? projectors_ell[krn - 1] : nullptr;
                        std::vector<double> rhs;
                        if (use_energy_derivative && nrn > 0 && krn > HOM) {
                            rhs = std::vector<double>(stride);
                            product(rhs.data(), rg[SMT]->n, rg[SMT]->r, partial_wave[iln - 1].wave[SMT]); // r*phi_0
                            projector = rhs.data(); // use as inhomogeneiety
                            // comment: this will be 2x the same solution, for krn==1 and krn==2, code results from a quick fix
                        }

                        double dg; // derivative at end point
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
                            scale(rphi[krn], rg[SMT]->n, rg[SMT]->rinv, scal); // scale and divide by r
                            // ToDo: extrapolate lowest point?

                            if (use_energy_derivative && nrn > 0) {
                                // (T + V - E) phi_0 = linear combinartion of projectors
                                // (T + V - E) phi_1 = phi_0 
                                // --> T phi_1 = (E - V) phi_1 + phi_0
                                for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                                    Tphi[krn][ir] = (vs.energy*rg[SMT]->r[ir] - potential[SMT][ir])*rphi[krn][ir] + scal*rhs[ir];
                                } // ir
                                // ToDo: check these equations and normalization factors
                                if (echo > 1) printf("# %s generate Tphi with inhomogeneiety\n", label);
                                // seems like the tails of TRU and SMT wave and wKin are deviating slightly beyond r_match

                            } else {
                                // Tphi = (E - V)*phi + projector
                                for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                                    Tphi[krn][ir] = (vs.energy*rg[SMT]->r[ir] - potential[SMT][ir])*rphi[krn][ir] + scal*projector[ir];
                                } // ir
                            }

                            // now visually check that the matching of value and derivative of rphi is ok.
                            if (echo > 19) {
                                printf("\n## %s check matching of rphi for ell=%i nrn=%i krn=%i (r, phi_tru, phi_smt, prj, rTphi_tru, rTphi_smt):\n",
                                        label, ell, nrn, krn-1);
                                for(int ir = 1; ir < rg[SMT]->n; ++ir) {
                                    printf("%g  %g %g  %g %g\n", rg[SMT]->r[ir],
                                        vs.wave[TRU][ir + nr_diff], rphi[krn][ir], projectors_ell[krn][ir]*rg[SMT]->rinv[ir],
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

                    char const minimize_radial_curvature = 'm'; // as suggested by Morian Sonnet: minimize the radial curvature of the smooth partial wave
                    char const energy_ordering = 'e';           // as suggested by Baumeister+Tsukamoto in PASC19 proceedings
                    char const orthogonalize_second = '2';      // take the same lowest partial wave as for nn==1 and use the freefom of the second
                    char const orthogonalize_first = '1';       // ToDo
                    
                    char const method = *control::get("single_atom.partial.wave.method", "m");

                    std::vector<double> evec(n, 0.0);
                    evec[nrn] = 1.0; // this is everything that needs to be done for method=='e' 
                    if (n > 1) {
                        
                        if (method == minimize_radial_curvature) { 
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
                                    printf("# %s generalized eigenvalue problem failed info=%i\n", label, int(info));
                                } else {
                                    set(evec.data(), n, Ekin[0]);
                                    if (echo > 6) {
                                        printf("# %s lowest eigenvalue of the radial curvature is %g %s", label, lowest_eigenvalue*eV, _eV);
                                        if (echo > 8) {
                                            printf(", coefficients"); for(int krn = 0; krn < n; ++krn) printf(" %g", Ekin[0][krn]);
                                        } // high echo
                                        printf("\n");
                                    } // echo
                                } // info
                            } // scope

                        } else // method
                        if (method == orthogonalize_second) {
                         
                            if (nrn > 0) { // only the second
                                assert(1 == nrn); // we do not know how to treat a 3rd wave yet
                                
                                // minimize the matrix element <Psi_1|p_0>
                                double c[2];
                                for(int krn = 0; krn < 2; ++krn) {
                                    c[krn] = dot_product(nr, rphi[1 + krn], projectors_ell[0], rg[SMT]->rdr);
                                } // krn
                                // now find the angle phi such that cos(phi)*c[1] + sin(phi)*c[0] is zero;
#if 0
                                for (int iang = -18; iang <= 18; ++iang) {
                                    double const angle = iang * constants::pi / 18;
                                    double const ovl10 = std::cos(angle)*c[1] + std::sin(angle)*c[0];
                                    printf("# method=orthogonalize_second angle=%g\t<Psi_1|p_0>= %g\n", angle, ovl10);
                                } // iang
#endif
                                double const angle = std::atan2(-c[1], c[0]);
                                evec[0] = std::sin(angle);
                                evec[1] = std::cos(angle);
                                {
                                    double const ovl10 = evec[0]*c[0] + evec[1]*c[1];
                                    if (echo > 8) printf("# %s method=orthogonalize_second angle=%g\t<Psi_1|p_0>= %g coeffs= %g %g\n", label, angle, ovl10, evec[0], evec[1]);
                                }
                                
                            } // if nrn > 0

                        } else // method
                        if (method == orthogonalize_first) {
                         
                            if (0 == nrn) { // only the first
                                
                                // minimize the matrix element <Psi_0|p_1>
                                double c[2];
                                for(int krn = 0; krn < 2; ++krn) {
                                    c[krn] = dot_product(nr, rphi[1 + krn], projectors_ell[1], rg[SMT]->rdr);
                                } // krn
                                // now find the angle phi such that cos(phi)*c[1] + sin(phi)*c[0] is zero;
                                double const angle = std::atan2(-c[1], c[0]);
                                evec[0] = std::sin(angle);
                                evec[1] = std::cos(angle);
                                {
                                    double const ovl01 = evec[0]*c[0] + evec[1]*c[1];
                                    if (echo > -1) printf("# method=orthogonalize_first angle=%g\t<Psi_0|p_1>= %g coeffs= %g %g\n", angle, ovl01, evec[0], evec[1]);
                                }
                                
                            } // if nrn > 0

                        } // method
                    } // n > 1 more than one partial wave
                    if (method == energy_ordering) assert(1 == evec[nrn]);

                    set(vs.wave[SMT], rg[SMT]->n, 0.0);
                    set(vs.wKin[SMT], rg[SMT]->n, 0.0);
                    double sum{0}; for(int krn = 0; krn < n; ++krn) sum += evec[krn];
                    double const scal = 1.0/sum; // scaling such that the sum is 1
                    for(int krn = 0; krn < n; ++krn) {
                        auto const coeff = scal*evec[krn];
                        add_product(vs.wave[SMT], rg[SMT]->n, rphi[1 + krn], coeff);
                        add_product(vs.wKin[SMT], rg[SMT]->n, Tphi[1 + krn], coeff);
                    } // krn

                    if (echo > 19) {
                        printf("\n## %s check matching of partial waves ell=%i nrn=%i (r, phi_tru, phi_smt, rTphi_tru, rTphi_smt):\n",
                                label, ell, nrn);
                        for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                            printf("%g  %g %g  %g %g\n", rg[SMT]->r[ir],
                                vs.wave[TRU][ir + nr_diff], vs.wave[SMT][ir],
                                vs.wKin[TRU][ir + nr_diff], vs.wKin[SMT][ir]);
                        } // ir
                        printf("\n\n");
                    } // echo

//                  if (echo > 8) { printf("# fflush(stdout) in line %i\n", __LINE__ + 1); fflush(stdout); }

                } // scope


            } // nrn


            { // scope: establish dual orthgonality with [SHO] projectors
                int const nr = rg[SMT]->n; //, mr = align<2>(nr);
                // int const stride = align<2>(rg[SMT]->n);
                int const msub = (1 + numax/2); // max. size of the subspace

                if (echo > 24) { // show normalization and orthogonality of projectors
                    for(int nrn = 0; nrn < n; ++nrn) {
                        for(int krn = 0; krn < n; ++krn) {
                            printf("# %s %c-projector <#%d|#%d> = %i + %.1e sigma=%g %s\n", label, ellchar[ell], nrn, krn,
                                (nrn == krn), dot_product(nr, projectors_ell[nrn], projectors_ell[krn], rg[SMT]->dr) - (nrn == krn), sigma*Ang,_Ang);
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
                            if (ts == TRU && echo > 3) printf("# %s create orthogonalized partial %c-wave #%i with %g * wave #%i\n", label, ellchar[ell], nrn, inv[nrn][krn], krn);
                            add_product(partial_wave[iln].wave[ts], nrts, waves[0][krn], inv[nrn][krn]);
                            add_product(partial_wave[iln].wKin[ts], nrts, waves[1][krn], inv[nrn][krn]);
                        } // krn
                    } // nrn
                } // ts in {tru, smt}

                if (n > 0) { // scope: check the overlap again
                    double dev{0}; // largest deviations from unity
                    for(int nrn = 0; nrn < n; ++nrn) { // smooth number or radial nodes
                        int const iln = ln_off + nrn;
                        auto const wave = partial_wave[iln].wave[SMT];
                        for(int krn = 0; krn < n; ++krn) { // smooth number or radial nodes
                            ovl[nrn][krn] = dot_product(nr, wave, projectors_ell[krn], rg[SMT]->rdr);
                            if (echo > 7) printf("# %s smooth partial %c-wave #%d with %c-projector #%d new overlap %g\n",
                                                  label, ellchar[ell], nrn, ellchar[ell], krn, ovl[nrn][krn]);
                            dev = std::max(dev, std::abs(ovl[nrn][krn] - (nrn == krn)));
                        } // krn
                    } // nrn
                    if (echo > 2) printf("# %s after orthogonalization <partial waves|projectors> deviates max. %.1e from unity matrix\n", label, dev);
                    if (dev > 1e-12) warn("%s duality violated, deviates %g from unity", label, dev);
                } // scope

                if (echo > 15) {
                    for(int nrn = 0; nrn < n; ++nrn) { // smooth number or radial nodes
                        auto const & vs = partial_wave[ln_off + nrn];
                        printf("\n## %s show orthogonalized partial %c-waves for nrn=%i (r, phi_tru, phi_smt, rTphi_tru, rTphi_smt):\n",
                                label, ellchar[ell], nrn);
                        for(int ir = 1; ir < rg[SMT]->n; ++ir) {
                            printf("%g  %g %g  %g %g\n", rg[SMT]->r[ir],
                                vs.wave[TRU][ir + nr_diff], vs.wave[SMT][ir],
                                vs.wKin[TRU][ir + nr_diff], vs.wKin[SMT][ir]);
                        } // ir
                        printf("\n\n");
                    } // nrn
                } // echo

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
                        if (echo > 19) printf("# %s %c-channel <%d|T|%d> kinetic energy [unsymmetrized] (true) %g and (smooth) %g (diff) %g %s\n",
                          label, ellchar[ell], i, j, E_kin_tru*eV, E_kin_smt*eV, (E_kin_tru - E_kin_smt)*eV, _eV);
                    } // j
                } // i

                if (1) { // symmetrize the kinetic energy tensor
                    for(int iln = 0 + ln_off; iln < n + ln_off; ++iln) {
                        for(int jln = 0 + ln_off; jln < iln; ++jln) { // triangular loop excluding the diagonal elements
                            if (1) { // usual
                                for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
                                    symmetrize(kinetic_energy(ts,iln,jln), kinetic_energy(ts,jln,iln));
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
                        if (echo > 19) printf("# %s %c-channel <%d|T|%d> kinetic energy [symmetrized] (true) %g and (smooth) %g (diff) %g %s\n",
                          label, ellchar[ell], i, j, E_kin_tru*eV, E_kin_smt*eV, (E_kin_tru - E_kin_smt)*eV, _eV);
                    } // j
                } // i

            } // scope: establish dual orthgonality with [SHO] projectors

        } // ell

    } // update_partial_waves

    
    
    void update_charge_deficit(int const echo=0) {
        int const nln = sho_tools::nSHO_radial(numax);
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
                bool const echo_l = (echo > 4 + 4*(ell > 0));
                if (echo_l) printf("# %s charges for ell=%d, jln = 0, 1, ...\n", label, ell);
                if (ell > 0) scale(rl.data(), nr, rg[ts]->r); // create r^{\ell}
                for(int iln = 0; iln < nln; ++iln) {
                    if (echo_l) printf("# %s %s iln = %d ", label, ts?"smt":"tru", iln);
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
            
        } else if (sho_tools::order_lmn == order) {
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
                    if (std::abs(rho_ij) > 1e-9) printf("# %s LINE=%d rho_ij = %g for lm=%d iln=%d jln=%d\n", 
                                                           label, __LINE__, rho_ij*Y00inv, lm, iln, jln);
                } // jln
            } // iln
        } // lm
#endif
    } // get_rho_tensor


    void update_full_density(view3D<double> const & density_tensor, int const echo=0) { // density tensor rho_{lm iln jln}
        int const nlm = pow2(1 + ellmax);
        int const nln = sho_tools::nSHO_radial(numax);
        // view3D<double const> density_tensor(rho_tensor, nln, nln); // rho_tensor[mln][nln][nln]

        double const mix_valence_density = 1.0 - mix_spherical_valence_density; // fraction of valence density from partial waves
        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            size_t const nr = rg[ts]->n;
            assert(full_density[ts].stride() >= nr);
            set(full_density[ts], nlm, 0.0); // clear
            for(int lm = 0; lm < nlm; ++lm) {
                if (00 == lm) {
                    set(full_density[ts][lm], nr, core_density[ts].data(), Y00); // needs scaling with Y00 since core_density has a factor 4*pi
                    if (echo > 0) printf("# %s %s density has %g electrons after adding the core density\n", label,
                        (TRU == ts)?"true":"smooth", dot_product(nr, full_density[ts][lm], rg[ts]->r2dr)*Y00inv);
                    if (mix_spherical_valence_density > 0) {
                        add_product(full_density[ts][lm], nr, spherical_valence_density[ts].data(), Y00); // needs scaling with Y00 since core_density has a factor 4*pi
                        if (echo > 0) printf("# %s %s density has %g electrons after adding the spherical valence density\n", label,
                            (TRU == ts)?"true":"smooth", dot_product(nr, full_density[ts][lm], rg[ts]->r2dr)*Y00inv);
                    } // use spherical_valence_density
                } // 00 == lm
                if (mix_valence_density > 0) {
                    for(int iln = 0; iln < nln; ++iln) {
                        auto const wave_i = partial_wave[iln].wave[ts];
                        assert(nullptr != wave_i);
                        for(int jln = 0; jln < nln; ++jln) {
                            auto const wave_j = partial_wave[jln].wave[ts];
                            assert(nullptr != wave_j);
                            double const rho_ij = density_tensor(lm,iln,jln) * mix_valence_density;
//                          if (std::abs(rho_ij) > 1e-9) printf("# %s rho_ij = %g for lm=%d iln=%d jln=%d\n", label, rho_ij*Y00inv, lm, iln, jln);
                            add_product(full_density[ts][lm], nr, wave_i, wave_j, rho_ij);
                        } // jln
                    } // iln
                } // mix_valence_density
            } // lm
            if (echo > 0) printf("# %s %s density has %g electrons after adding the valence density\n",  label,
                    (TRU == ts)?"true":"smooth", dot_product(nr, full_density[ts][00], rg[ts]->r2dr)*Y00inv);

        } // true and smooth

        int const nlm_cmp = pow2(1 + ellmax_compensator);
        for(int ell = 0; ell <= ellmax_compensator; ++ell) {
            for(int emm = -ell; emm <= ell; ++emm) {
                int const lm = solid_harmonics::lm_index(ell, emm);
                double rho_lm{0};
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
        qlm_compensator[0] += Y00*(core_charge_deficit - Z_core + mix_spherical_valence_density*spherical_valence_charge_deficit);
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
        std::vector<double> vlm(nlm, 0.0);
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
                { // scope: transform also the exchange-correlation energy density
                    view2D<double> exc_lm(nlm, mr);
                    angular_grid::transform(exc_lm.data(), exc_on_grid, mr, ellmax, true);
                    if (SMT == ts) {
                        if (echo > 7) printf("# %s local smooth exchange-correlation potential at origin is %g %s\n",
                                                            label, full_potential[SMT][00][0]*Y00*eV,_eV);
                    } // SMT only
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
                add_or_project_compensators<1>(Ves, vlm.data(), ellmax_compensator, rg[SMT], sigma_compensator); // project Ves to compensators
                if (echo > 7) printf("# %s inner integral between normalized compensator and smooth Ves(r) = %g %s\n", label, vlm[0]*Y00*eV,_eV);

                // but the solution of the 3D problem found that these integrals should be ves_multipole, therefore
                if (nullptr == ves_multipole) {
                    set(vlm.data(), nlm, 0.); // no correction of the electrostatic potential heights for isolated atoms
                } else {
                    if (echo > 6) printf("# %s v_00 found %g but expected %g %s\n", label, vlm[0]*Y00*eV, ves_multipole[0]*Y00*eV,_eV);
                    scale(vlm.data(), nlm, -1.); add_product(vlm.data(), nlm, ves_multipole, 1.); // vlm := ves_multipole - vlm
                } // no ves_multipole given
            } // smooth only

            if (SMT == ts) {
                if (echo > 7) printf("# %s local smooth electrostatic potential at origin is %g %s\n", label, Ves[00][0]*Y00*eV,_eV);
            } // SMT only

            add_or_project_compensators<2>(Ves, vlm.data(), ellmax_compensator, rg[ts], sigma_compensator);

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
                if (echo > 5) {
                      auto const Enuc = -Z_core*dot_product(nr, full_density[TRU][00], rg[ts]->rdr)*Y00inv;
                      printf("# %s Coulomb energy is %.15g %s\n", label, Enuc*eV, _eV);
                } // echo
            } // ts

            add_product(full_potential[ts].data(), nlm*mr, Ves.data(), 1.0); // add the electrostatic potential, scale_factor=1.0

        } // true and smooth

        // construct the zero_potential V_bar
        std::vector<double> V_smt(rg[SMT]->n);
        set(zero_potential.data(), rg[SMT]->n, 0.0); // init zero
        
//         set(V_smt.data(), rg[SMT]->n, full_potential[TRU][00] + nr_diff); // copy the tail of the spherical part of the true potential
        auto const df = Y00*eV; assert(df > 0); // display factor
//         if (echo > 5) printf("# %s match local potential to parabola at R_cut = %g %s, V_tru(R_cut) = %g %s\n",
//                     label, rg[SMT]->r[ir_cut[SMT]]*Ang, _Ang, full_potential[TRU][00][ir_cut[TRU]]*df, _eV);
//         auto const stat = pseudize_function(V_smt.data(), rg[SMT], ir_cut[SMT], 2);
        auto const stat = pseudize_local_potential<0>(V_smt.data(), full_potential[TRU][00], echo, Y00); // <0> indicates that these arrays hold V(r) (not r*V(r))

        if (stat) {
            if (echo > 0) printf("# %s matching procedure for the local potential failed! status= %d\n", label, int(stat));
        } else {
//             if (echo > -1) printf("# local smooth zero_potential:\n");
            for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                zero_potential[ir] = V_smt[ir] - full_potential[SMT][00][ir];
//                 if (echo > -1) printf("%g %g\n", rg[SMT]->r[ir], zero_potential[ir]*Y00);
            } // ir
//             if (echo > -1) printf("\n\n");
            if (echo > 5) printf("# %s smooth potential: V_smt(0) = %g, V_smt(R_cut) = %g %s\n",
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
                // and how much V_smt deviates from V_tru outside the sphere
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
            printf("\n# WARNING: use rad_pot to construct the r*V_tru(r) [for DEBUGGING]\n\n");
            atom_core::rad_pot(potential[TRU].data(), *rg[TRU], rho4pi.data(), Z_core);
        } // scope

    } // update_full_potential

    void update_matrix_elements(int const echo=0) {
        int const nlm = pow2(1 + ellmax);
        int const mlm = pow2(1 + numax);
        int const nln = sho_tools::nSHO_radial(numax);
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
                    printf(" %g", true_norm[iln]*overlap_lmn[ilmn][jlmn]*true_norm[jln]);
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
                scattering_test::emm_average(result_ln.data(), input_lmn.data(), int(numax), nn);
                if (echo > 4) {
                    for(int iln = 0; iln < nln; ++iln) {
                        printf("# %s emm-averaged%3i %s ", label, iln, label_inp);
                        for(int jln = 0; jln < nln; ++jln) {
                            printf(" %11.6f", true_norm[iln]*result_ln[iln][jln]*true_norm[jln]);
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
                  (*rg[SMT], Vsmt.data(), sigma, int(numax + 1), nn, numax, hamiltonian_ln.data(), overlap_ln.data(), 384, V_rmax, label, echo);
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
                  (rg, rV, sigma, int(numax + 1), nn, numax, hamiltonian_ln.data(), overlap_ln.data(), logder_energy_range, label, echo);
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

        view3D<double> potential_ln(TRU_AND_SMT, nln, nln); // fully emm-degenerate
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

                    double lower{0}, upper{0};
                    if (1 == nn[ell]) {
                        int const iln = sho_tools::ln_index(numax, ell, 0);
                        lower = overlap_ln[iln][iln];
                    } else if (nn[ell] > 1) {
                        int const iln = sho_tools::ln_index(numax, ell, 0);
                        // get eigenvalues of the 2x2 matrix [a,b]
                        //                                   [c,d] with c==b (symmetric)
                        auto const a = overlap_ln[iln][iln];
                        auto const b = overlap_ln[iln][iln + 1];
                        auto const d = overlap_ln[iln + 1][iln + 1];
                        auto const split = std::sqrt(pow2(a - d) + 4*pow2(b));
                        upper = 0.5*(a + d + split);
                        lower = 0.5*(a + d - split);
                    }
                    if (lower <= -1.) { 
                        warn("%s eigenvalues of charge deficit matrix for ell=%i instable! %g and %g", label, ell, lower, upper);
                    } else if (lower < -0.9) {
                        warn("%s eigenvalues of charge deficit matrix for ell=%i critical! %g and %g", label, ell, lower, upper);
                    }
                } // ell
            }

            if (echo > 4) {
                printf("\n");
                view2D<char> ln_labels(nln, 4);
                sho_tools::construct_label_table<4>(ln_labels.data(), numax, sho_tools::order_ln);
                if (0) {
                    printf("\n# %s with true_norm scaling:\n", label);
                    for(int i01 = 0; i01 < 2; ++i01) {
                        auto const & input_ln = i01 ?  overlap_ln :  hamiltonian_ln;
                        auto const   label_qn = i01 ? "overlap"   : "hamiltonian" ;
                        for(int iln = 0; iln < nln; ++iln) {
                            printf("# %s spherical %s %s ", label, ln_labels[iln], label_qn);
                            for(int jln = 0; jln < nln; ++jln) {
                                printf(" %g", (ells[iln] == ells[jln]) ? true_norm[iln]*input_ln[iln][jln]*true_norm[jln] : 0);
                            }   printf("\n");
                        }   printf("\n");
                    } // i01
                } // never

                printf("\n# %s without true_norm scaling:\n", label);
                for(int i01 = 0; i01 < 2; ++i01) {
                    auto const & input_ln = i01 ?  overlap_ln :  hamiltonian_ln;
                    auto const   label_qn = i01 ? "overlap"   : "hamiltonian" ;
                    for(int iln = 0; iln < nln; ++iln) {
                        printf("# %s spherical %s %s ", label, ln_labels[iln], label_qn);
                        for(int jln = 0; jln < nln; ++jln) {
                            printf(" %g", (ells[iln] == ells[jln]) ? input_ln[iln][jln] : 0);
                        }   printf("\n");
                    }   printf("\n");
                } // i01

            } // echo
            
        } // 1

#if 0        
#ifdef NEVER
                    if (echo > 1) printf("\n# %s %s multiply the true_norm factor to hamiltonian_ln and overlap_ln\n", label, __func__);
                    for(int iln = 0; iln < nln; ++iln) {
                        for(int jln = 0; jln < nln; ++jln) {
                            hamiltonian_ln(iln,jln) *= true_norm[iln]*true_norm[jln];
                            overlap_ln    (iln,jln) *= true_norm[iln]*true_norm[jln];
                        } // ilm
                    } // iln
#endif
#endif
        
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
              (*rg[SMT], Vsmt.data(), sigma, int(numax + 1), nn, numax, hamiltonian_ln.data(), overlap_ln.data(), 384, V_rmax, label, echo);
        } else if (echo > 0) printf("\n# eigenstate_analysis deactivated for now! %s %s:%i\n\n", __func__, __FILE__, __LINE__);

        if (1) {
            if (echo > 1) printf("\n# %s %s logarithmic_derivative\n\n", label, __func__);
            double const *rV[TRU_AND_SMT] = {potential[TRU].data(), potential[SMT].data()};
            scattering_test::logarithmic_derivative // scan the logarithmic derivatives
              (rg, rV, sigma, int(numax + 1), nn, numax, hamiltonian_ln.data(), overlap_ln.data(), logder_energy_range, label, echo);
        } else if (echo > 0) printf("\n# logarithmic_derivative deactivated for now! %s %s:%i\n\n", __func__, __FILE__, __LINE__);

    } // check_spherical_matrix_elements

    void update_density(float const mixing, int const echo=0) {
//         if (echo > 2) printf("\n# %s\n", __func__);
        update_core_states(mixing, echo);
        update_energy_parameters(echo);
        update_partial_waves(echo); // create new partial waves for the valence description
        update_charge_deficit(echo); // update quantities derived from the partial waves
        int const nln = sho_tools::nSHO_radial(numax);
        view3D<double> aHSm(2, nln, nln, 0.0);
        check_spherical_matrix_elements(echo, aHSm); // check scattering properties for emm-averaged Hamiltonian elements
        int const nSHO = sho_tools::nSHO(numax);
        view2D<double> density_matrix(nSHO, nSHO, 0.0); // get memory
        auto dm_order = sho_tools::order_Ezyx;
        int const lmax = std::max(ellmax, ellmax_compensator);
        int const mlm = pow2(1 + lmax);
        view3D<double> rho_tensor(mlm, nln, nln, 0.0); // get memory
        get_rho_tensor(rho_tensor, density_matrix, dm_order, echo);
        update_full_density(rho_tensor, echo);
    } // update_density

    // ==============
    // between update_density and update_potential we need to
    // export qlm_compensator, solve the 3D electrostatic problem, return here with ves_multipoles
    // ==============

    void update_potential(float const mixing, double const ves_multipoles[], int const echo=0) {
        if (echo > 2) printf("\n# %s %s\n", label, __func__);
        update_full_potential(mixing, ves_multipoles, echo);
        update_matrix_elements(echo); // this line does not compile with icpc (ICC) 19.0.2.187 20190117
    } // update_potential

    status_t get_smooth_spherical_quantity(double qnt[] // result array: function on an r2-grid
        , float const ar2, int const nr2 // r2-grid parameters
        , char const what, int const echo=1) const { // logg level
        char const *qnt_name   = ('c' == what) ? "core_density"     : (('z' == what) ? "zero_potential" : "valence_density");
        auto const &qnt_vector = ('c' == what) ?  core_density[SMT] : (('z' == what) ?  zero_potential  : spherical_valence_density[SMT]);
        if (echo > 8) printf("# %s call transform_to_r2_grid(%p, %.1f, %d, %s=%p, rg=%p)\n",
                        label, (void*)qnt, ar2, nr2, qnt_name, (void*)qnt_vector.data(), (void*)rg[SMT]);
        double const minval = ('z' == what) ? -9e307 : 0.0; // zero potential may be negative, densities should not
#ifdef DEVEL
        double const Y00s = Y00*(('c' == what) ? Y00 : 1); // Y00 for zero_pot and Y00^2 for rho_core
        if (echo > 8) {
            printf("\n## %s %s before filtering:\n", label, qnt_name);
            for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                printf("%g %g\n", rg[SMT]->r[ir], qnt_vector[ir]*Y00s);
            }   printf("\n\n");
        } // echo
#endif

        // divide input by mask function
        std::vector<double> inp(rg[SMT]->n);
        double const r2cut = pow2(rg[SMT]->rmax)*1.1, r2inv = 1./r2cut;
        for(int ir = 0; ir < rg[SMT]->n; ++ir) {
            double const r2 = pow2(rg[SMT]->r[ir]);
            auto const mask = pow8(1. - pow8(r2*r2inv));
            inp[ir] = qnt_vector[ir] / mask;
        } // ir

        auto const stat = bessel_transform::transform_to_r2_grid(qnt, ar2, nr2, inp.data(), *rg[SMT], echo);

        // multiply output by mask function
        for(int ir2 = 0; ir2 < nr2; ++ir2) {
            double const r2 = ir2/ar2;
            auto const mask = (r2 < r2cut) ? pow8(1. - pow8(r2*r2inv)) : 0;
            qnt[ir2] = std::max(minval, qnt[ir2]) * mask;
        } // ir2

#ifdef DEVEL
        if (echo > 8) {
            printf("\n## %s %s  after filtering:\n", label, qnt_name);
            for(int ir2 = 0; ir2 < nr2; ++ir2) {
                printf("%g %g\n", std::sqrt(ir2/ar2), qnt[ir2]*Y00s);
            }   printf("\n\n");
        } // echo
#endif
        return stat;
    } // get_smooth_spherical_quantity
    
    radial_grid_t* get_smooth_radial_grid(int const echo=0) const { return rg[SMT]; }

    template <char Q='t'> double get_number_of_electrons() const { return csv_charge[0] + csv_charge[1] + csv_charge[2]; }
    
    int const get_numax() const { return numax; }
    
  }; // class LiveAtom

    template <> double LiveAtom::get_number_of_electrons<'c'>() const { return csv_charge[0]; } // core
    template <> double LiveAtom::get_number_of_electrons<'s'>() const { return csv_charge[1]; } // semicore
    template <> double LiveAtom::get_number_of_electrons<'v'>() const { return csv_charge[2]; } // valence

namespace single_atom {

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
#else
      return s[0] | (uint64_t(s[1]) <<  8)
                  | (uint64_t(s[2]) << 16) 
                  | (uint64_t(s[3]) << 24)
                  | (uint64_t(s[4]) << 32)
                  | (uint64_t(s[5]) << 40)
                  | (uint64_t(s[6]) << 48)
                  | (uint64_t(s[7]) << 56);
#endif      
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
#else
      return s[0] | (uint32_t(s[1]) <<  8) 
                  | (uint32_t(s[2]) << 16) 
                  | (uint32_t(s[3]) << 24);
#endif
  } // string2int
  
  // from https://hbfs.wordpress.com/2017/01/10/strings-in-c-switchcase-statements/
  inline uint64_t constexpr __mix_hash_(char const m, uint64_t const s) { return ((s << 7) + ~(s >> 3)) + ~m; }
  inline uint64_t constexpr string2hash(char const * m) { return (*m) ? __mix_hash_(*m, string2hash(m + 1)) : 0; }

  status_t test_string_switch(char const *what, int const echo=0) {
      #define str2int string2hash
      switch ( str2int(what) ) {
          case str2int("initialize"):
          case str2int("memory cleanup"):
          case str2int("lmax qlm"):
          case str2int("lmax vlm"): // does not work with uint32_t
          case str2int("sigma cmp"):
          case str2int("core densities"):
          case str2int("valence densities"):
          case str2int("projectors"):
          case str2int("qlm charges"):
          case str2int("update"): // needs a trailing ' ' if str2int==string2long
          case str2int("hamiltonian"):
          case str2int("zero potentials"):
          case str2int("radial grids"):
            if (echo > 0) printf("# %s found selector what=\"%s\".\n", __func__, what);
          break; default:
            if (echo > 0) printf("# %s unknown selector what=\"%s\"!\n", __func__, what);
            return 1; // error
      } // switch (what)
      #undef  str2int
      return 0; // no error
  } // test_string_switch

  // simplified interface compared to update
  status_t atom_update(char const *what, int const natoms,
              double *dp, int32_t *ip, float *fp, double **dpp) {

      static std::vector<LiveAtom*> a; // this function may not be templated!!
      static int echo = -9;
      if (-9 == echo) echo = control::get("single_atom.echo", 0.); // initialize only on the 1st call to update()

      assert(what); // should never be nullptr
      char const how = what[0] | 32; // first char, to lowercase
      if (echo > 4) printf("\n# %s %s what=\"%s\" --> \'%c\'\n\n", __FILE__, __func__, what, how);

      float constexpr ar2_default = 16.f;
      int   constexpr nr2_default = 1 << 12;
      int   constexpr numax_default = 3;
      float constexpr mix_rho_default = .5f;
      float constexpr mix_pot_default = .5f;
      
      int na{natoms};
      
      status_t stat(0);
      stat = test_string_switch(what, echo);

      switch (how) {
        
          case 'i': // interface usage: atom_update("initialize", natoms, Za[], numax[], ion[]=0);
          {
              double const *Za = dp; assert(nullptr != Za); // may not be nullptr as it holds the atomic core charge Z[ia]
              a.resize(na);
              bool const atomic_valence_density = (nullptr != dpp); // control
              int const echo_init = control::get("single_atom.init.echo", 0.); // log-level for the LiveAtom constructor
              for(int ia = 0; ia < a.size(); ++ia) {
                  float const ion = (fp) ? fp[ia] : 0;
                  a[ia] = new LiveAtom(Za[ia], numax_default, atomic_valence_density, ion, ia, echo_init);
                  if(ip) ip[ia] = a[ia]->get_numax(); // export numax
              } // ia
          } 
          break;

          case 'm': // interface usage: atom_update("memory cleanup", natoms);
          {
              for(int ia = 0; ia < a.size(); ++ia) {
                  a[ia]->~LiveAtom(); // envoke destructor
              } // ia
              a.clear();
              na = a.size(); // fulfill consistency check
              assert(!dp); assert(!ip); assert(!fp); assert(!dpp); // all other arguments must be nullptr (by default)
          } 
          break;

          case 'c': // interface usage: atom_update("core densities",    natoms, null, nr2=2^12, ar2=16.f, qnt=rho_c);
          case 'v': // interface usage: atom_update("valence densities", natoms, q_00, nr2=2^12, ar2=16.f, qnt=rho_v);
          case 'z': // interface usage: atom_update("zero potentials",   natoms, null, nr2=2^12, ar2=16.f, qnt=v_bar);
          {
              double *const *const qnt = dpp; assert(nullptr != qnt);
              if ('v' == how) assert(nullptr != dp);
              for(int ia = 0; ia < a.size(); ++ia) {
                  assert(nullptr != qnt[ia]);
                  int   const nr2 = ip ? ip[ia] : nr2_default;
                  float const ar2 = fp ? fp[ia] : ar2_default;
                  stat += a[ia]->get_smooth_spherical_quantity(qnt[ia], ar2, nr2, how);
                  if ('v' == how) dp[ia] = a[ia]->spherical_valence_charge_deficit;
              } // ia
              if ('v' != how && nullptr != dp) warn("please use \'%s\'-interface with nullptr as 3rd argument", what);
          } 
          break;

          case 'r': // interface usage: atom_update("radial grid", natoms, null, null, null, (double**)rg_ptr);
          {
              assert(nullptr != dpp);
              for(int ia = 0; ia < a.size(); ++ia) {
#ifdef  DEVEL
                  dpp[ia] = reinterpret_cast<double*>(a[ia]->get_smooth_radial_grid()); // pointers to smooth radial grids
#else
                  dpp[ia] = nullptr;
                  stat = -1;
                  warn("# %s %s only available with -D DEVEL", __func__, what);
#endif
              } // ia
              assert(!dp); assert(!ip); assert(!fp); // all other arguments must be nullptr (by default)
          } 
          break;

          case 's': // interface usage: atom_update("sigma compensator", natoms, sigma);
          {
              double *const sigma = dp; assert(nullptr != sigma);
              for(int ia = 0; ia < a.size(); ++ia) {
                  sigma[ia] = a[ia]->sigma_compensator; // spreads of the compensators // ToDo: use a getter function
              } // ia
              assert(!ip); assert(!fp); assert(!dpp); // all other arguments must be nullptr (by default)
          } 
          break;

          case 'p': // interface usage: atom_update("projectors", natoms, sigma);
          {
              double  *const sigma = dp; assert(nullptr != sigma);
              int32_t *const numax = ip; assert(nullptr != numax);
              for(int ia = 0; ia < a.size(); ++ia) {
                  sigma[ia] = a[ia]->sigma; // spreads of the projectors // ToDo: use a getter function
                  numax[ia] = a[ia]->numax; //  number of SHO-projectors // ToDo: use a getter function
              } // ia
              assert(!fp); assert(!dpp); // all other arguments must be nullptr (by default)
          } 
          break;

          case 'u': // interface usage: atom_update("update", natoms, null, null, mix={mix_pot, mix_rho}, vlm);
          {
              double const *const *const vlm = dpp; assert(nullptr != vlm);
              float const *const mix = fp; // nullptr == fp is okay
              float const mix_pot = mix ? mix[0] : mix_pot_default;
              float const mix_rho = mix ? mix[1] : mix_rho_default;
              for(int ia = 0; ia < a.size(); ++ia) {
                  a[ia]->update_potential(mix_pot, vlm[ia], echo); // set electrostatic multipole shifts
                  a[ia]->update_density(mix_rho, echo);
              } // ia
              assert(!dp); assert(!ip); // all other arguments must be nullptr (by default)
          } 
          break;

          case 'a': // interface usage: atom_update("atomic density matrix", natoms, null, null, null, atom_rho);
          {
              double const *const *const atom_rho = dpp; assert(nullptr != atom_rho);
              for(int ia = 0; ia < a.size(); ++ia) {
                  assert(nullptr != atom_rho[ia]);
                  // ToDo: set the atomic density matrices of LiveAtoms a[ia]
              } // ia
              warn("Atomic density matrices are ignored!");
              stat = 1;
              assert(!dp); assert(!ip); assert(!fp); // all other arguments must be nullptr (by default)
          } 
          break;
              
          case 'q': // interface usage: atom_update("q_lm charges", natoms, null, null, null, qlm);
          {
              double *const *const qlm = dpp; assert(nullptr != qlm);
              for(int ia = 0; ia < a.size(); ++ia) {
                  int const nlm = pow2(1 + a[ia]->ellmax_compensator);
                  set(qlm[ia], nlm, a[ia]->qlm_compensator.data()); // copy compensator multipoles
              } // ia
              assert(!dp); assert(!ip); assert(!fp); // all other arguments must be nullptr (by default)
          } 
          break;

          case 'l': // interface usage: atom_update("lmax qlm", natoms, dp=null, lmax);
                    // interface usage: atom_update("lmax vlm", natoms, dp= ~0 , lmax);
          {
              int32_t *const lmax = ip; assert(nullptr != lmax);
              for(int ia = 0; ia < a.size(); ++ia) {
                  lmax[ia] = dp ? a[ia]->ellmax : a[ia]->ellmax_compensator;
              } // ia
              assert(!fp); assert(!dpp); // all other arguments must be nullptr (by default)
          } 
          break;

          case 'n': // interface usage: atom_update("numax", natoms, null, numax);
          {
              int32_t *const numax = ip; assert(nullptr != numax);
              for(int ia = 0; ia < a.size(); ++ia) {
                  numax[ia] = a[ia]->get_numax();
              } // ia
              assert(!dp); assert(!fp); assert(!dpp); // all other arguments must be nullptr (by default)
          } 
          break;
          
          case 'h': // interface usage: atom_update("hamiltonian and overlap", natoms, null, nelements, null, atom_mat);
          {
              double *const *const atom_mat = dpp; assert(nullptr != atom_mat);
              for(int ia = 0; ia < a.size(); ++ia) {
                  assert(nullptr != atom_mat[ia]);
                  int const numax = a[ia]->get_numax();
                  int const ncoeff = sho_tools::nSHO(numax);
                  int const h1hs2 = ip ? ip[ia]/(ncoeff*ncoeff) : 2;
                  for(int i = 0; i < ncoeff; ++i) {
                      if (h1hs2 > 1)
                      set(&atom_mat[ia][(1*ncoeff + i)*ncoeff + 0], ncoeff, a[ia]->overlap[i]);
                      set(&atom_mat[ia][(0*ncoeff + i)*ncoeff + 0], ncoeff, a[ia]->hamiltonian[i]);
                  } // i
              } // ia
              assert(!dp); assert(!fp); // all other arguments must be nullptr (by default)
          } 
          break;
            
          default:
          {
              if (echo > 0) printf("# %s first argument \'%s\' undefined, no action!\n", __func__, what);
              stat = how;
          } 
          break;

      } // switch(how)

      if (a.size() != na) {
          if (echo > 0) printf("# %s inconsistency detected: internally %ld atoms active by n=%d\n",
          __func__, a.size(), na);
      } // inconsistency
      
      if (stat) warn("what='%s' returns status = %i", what, int(stat));
      return stat;
  } // atom_update
  

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
    bool const avd = (control::get("single_atom.test.atomic.valence.density", 0.) > 0); //
//     for(int Z = 0; Z <= 109; ++Z) { // all elements
//     for(int Z = 109; Z >= 0; --Z) { // all elements backwards
//        if (echo > 1) printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
    {   double const Z   = control::get("single_atom.test.Z", 29.); // default copper
        double const ion = control::get("single_atom.test.ion", 0.); // default neutral
        if (echo > 1) printf("\n# Z = %g\n", Z);
        LiveAtom a(Z, numax, avd, ion, -1, echo); // envoke constructor
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


#define LIVE_ATOM_AS_LIBRARY
#ifdef  LIVE_ATOM_AS_LIBRARY
  // C and Fortran interface
  extern "C" {
      #include <cstdint> // or <stdint.h> in C
      #define SINGLE_ATOM_SOURCE    
        #include "single_atom.h"
      #undef  SINGLE_ATOM_SOURCE
  } // extern "C"
#endif
