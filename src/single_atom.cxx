#include <cstdio> // printf
#include <cstdlib> // abs
#include <cmath> // sqrt
#include <cassert> // assert
#include <algorithm> // max

#include "single_atom.hxx"

#include "radial_grid.hxx" // create_default_radial_grid, destroy_radial_grid, dot_product
#include "radial_eigensolver.hxx" // shooting_method
#include "radial_potential.hxx" // Hartree_potential
#include "angular_grid.hxx" // transform, Lebedev_grid_size
#include "radial_integrator.hxx" // integrate_outwards
#include "exchange_correlation.hxx" // lda_PZ81_kernel
#include "inline_tools.hxx" // align<nbits>
#include "sho_unitary.hxx" // Unitary_SHO_Transform<real_t>
#include "solid_harmonics.hxx" // lm_index, Y00, Y00inv
#include "atom_core.hxx" // initial_density, rad_pot, nl_index
#include "sho_tools.hxx" // lnm_index, nSHO
#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t, emm_Degenerate, spin_QN_t, spin_Degenerate
#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // pow2, pow3, set, scale, product, add_product
#include "bessel_transform.hxx" // transform_to_r2_grid

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
  void dgesv_(const int*, const int*, double*, const int*, int*, double*, const int*, int*);
} // extern "C"

  status_t solve_Ax_b(double x[], double const b[], double A[], int const n, int const stride) {
      std::copy(b, b+n, x); // copy right-hand-side into x
      status_t info = 0; int const nRHSs = 1; auto const ipivot = new int[n];
      dgesv_(&n, &nRHSs, A, &stride, ipivot, x, &n, &info);
      delete [] ipivot;
      return info;
  } // solve_Ax_b


  int constexpr TRU=0, SMT=1, TRU_AND_SMT=2, TRU_ONLY=1;
  int constexpr ELLMAX=7;
  char const ellchar[] = "spdfghijklmno";
  double const Y00 = solid_harmonics::Y00; // == 1./sqrt(4*pi)
  
  // 20 core states are the usual max., 32 core states are enough if spin-orbit-interaction is on

  // ToDo: write a class with constructors and destructors to handle the memory for *wave
  template<int Pseudo> // Pseudo=1: core states, Pseudo=2: valence states
  struct energy_level {
      double* wave[Pseudo]; // for valence states points to the true and smooth partial waves
      double energy; // energy level in Hartree atomic units
      float occupation; // occupation number
      enn_QN_t enn; // main quantum_number
      ell_QN_t ell; // angular momentum quantum_number
      emm_QN_t emm; // usually emm == emm_Degenerate
      spin_QN_t spin; // usually spin == spin_Degenerate
      enn_QN_t nrn[Pseudo]; // number of radial nodes
  };

  typedef struct energy_level<TRU_ONLY> core_level_t;
  typedef struct energy_level<TRU_AND_SMT> valence_level_t;

  
  status_t pseudize_function(double sfun[], radial_grid_t const *rg, int const irc, 
              int const nmax=4, int const ell=0, double *coeff=nullptr) {
      // match a radial function with an even-order polynomial inside r[irc]
      double Amat[4*4], x[4] = {0,0,0,0}; double* bvec = x;
      set(Amat, 4*4, 0.0);
      int const nm = std::min(std::max(1, nmax), 4);
      for(int i4 = 0; i4 < nm; ++i4) {
          // use up to 4 radial indices [irc-2,irc-1,irc+0,irc+1]
          int const ir = irc + i4 - nm/2; // not fully centered around irc
          double const r = rg->r[ir];
          double rl = intpow(r, ell);
          // set up a basis of 4 functions: r^0, r^0, r^4, r^6 at up to 4 neighboring grid points around r(irc)
          for(int j4 = 0; j4 < nm; ++j4) {
              Amat[j4*4 + i4] = rl;
              rl *= pow2(r);
          } // j4
          // b is the inhomogeneus right side of the set of linear equations
          bvec[i4] = sfun[ir];
      } // i4
      status_t info = solve_Ax_b(x, bvec, Amat, nm, 4);
      set(x + nm, 4 - nm, 0.0); // clear the unused coefficients
      // replace the inner part of the function by the even-order polynomial
      for(int ir = 0; ir < irc; ++ir) {
          double const r = rg->r[ir];
          double const rr = pow2(r);
          double const rl = intpow(r, ell);
          sfun[ir] = rl*(x[0] + rr*(x[1] + rr*(x[2] + rr*x[3])));
      } // ir
      if (nullptr != coeff) set(coeff, nm, x); // export coefficients
      return info;
  } // pseudize_function
  
  
  
    template<int ADD0_or_PROJECT1>
    void add_or_project_compensators(double out[], int const lmax, radial_grid_t const *rg, double const in[]
                                   , double const sigma_compensator, int const echo=0) {
        int const nr = rg->n, mr = align<2>(nr);
        auto const sig2inv = -.5/(sigma_compensator*sigma_compensator);
        if (echo > 0) printf("# sigma = %g\n", sigma_compensator);
        auto const rlgauss = new double[nr];
        auto const rl      = new double[nr];
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
                if (echo > 8) printf("# ell=%d norm=%g ir=%d rlgauss=%g rl=%g r2dr=%g\n", ell, norm, ir, rlgauss[ir], rl[ir], rg->r2dr[ir]);
            } // ir
            if (echo > 1) printf("# ell=%d norm=%g nr=%d\n", ell, norm, nr);
            assert(norm > 0);
            auto const scal = 1./norm;
            for(int emm = -ell; emm <= ell; ++emm) {
                int const lm = solid_harmonics::lm_index(ell, emm);
                if (0 == ADD0_or_PROJECT1) {
                    add_product(&out[lm*mr], nr, rlgauss, in[lm] * scal); // add normalized compensator to augmented density
                } else if (2 == ADD0_or_PROJECT1) {
                    add_product(&out[lm*mr], nr, rl, in[lm]); // add q_{\ell m} * r^\ell to electrostatic potential
                    // ToDo: setup of normalized rlgauss is not necessary in this version
                } else if (3 == ADD0_or_PROJECT1) {
                    out[lm] = dot_product(nr, &in[lm*mr], rl, rg->r2dr); // dot_product with metric
                    // ToDo: setup of normalized rlgauss is not necessary in this version
                } else {
                    out[lm] = dot_product(nr, &in[lm*mr], rlgauss, rg->r2dr) * scal; // dot_product with metric
                } // add or project
            } // emm
        } // ell
        delete[] rl;
        delete[] rlgauss;
    } // add_or_project_compensators

    void correct_multipole_shift(double ves[], int const lmax, radial_grid_t const *rg, double const vlm[], int const echo=0) {
        int const nr = rg->n, mr = align<2>(nr);
        auto const rl = new double[nr];
        set(rl, nr, 1.); // init with r^0
        for(int ell = 0; ell <= lmax; ++ell) { // serial!
            for(int emm = -ell; emm <= ell; ++emm) {
                int const lm = solid_harmonics::lm_index(ell, emm);
                add_product(&ves[lm*mr], nr, rl, vlm[lm]); // add q_{\ell m} * r^\ell to electrostatic potential
            } // emm
            scale(rl, nr, rg->r); // construct r^ell for the next ell-iteration
        } // ell
        delete[] rl;
    } // correct_multipole_shift
  
  
  class LiveAtom {
  public:
//    double* mem; // memory to which all radial function pointers and matrices point
      
      // general config
      int32_t id; // global atom identifyer
      float Z; // number of nucleons in the core
      radial_grid_t* rg[TRU_AND_SMT]; // radial grid descriptor for the true and smooth grid:
              // SMT may point to TRU, but at least both radial grids must have the same tail
      ell_QN_t ellmax; // limit ell for full_potential and full_density
      double r_cut; // classical augmentation radius for potential and core density
      int ir_cut[TRU_AND_SMT]; // classical augmentation radius index for potential and core density
      float r_match; // radius for matching of true and smooth partial wave, usually 6--9*sigma
      ell_QN_t numax; // limit of the SHO projector quantum numbers
      uint8_t nn[1+ELLMAX]; // number of projectors and partial waves used in each ell-channel
      ell_QN_t ellmax_compensator; // limit ell for the charge deficit compensators
      double sigma_compensator; // Gaussian spread for the charge deficit compensators
      double* qlm_compensator; // coefficients for the charge deficit compensators, (1+ellmax_compensator)^2
      double* aug_density; // augmented density, core + valence + compensation, (1+ellmax)^2 radial functions
      int matrix_stride; // stride for the matrices hamiltonian and overlap, must be >= nSHO(numax)
      int ncorestates; // for emm-Degenerate representations, 20 (or 32 with spin-oribit) core states are maximum
      int nvalencestates; // for emm-Degenerate (numax*(numax + 4) + 4)/4 is a good choice;
      int nspins; // 1 or 2 or 4
      double* unitary_zyx_lmn; // unitary sho transformation matrix [Cartesian][Radial], stride=nSHO(numax)

      // spin-resolved members of LiveAtom
      core_level_t* core_state;
      valence_level_t* valence_state;
      double* core_density[TRU_AND_SMT]; // spherical core density*4pi, no Y00 factor
      double* full_density[TRU_AND_SMT]; // total density, core + valence, (1+ellmax)^2 radial functions
      double* full_potential[TRU_AND_SMT]; // (1+ellmax)^2 radial functions
      double* potential[TRU_AND_SMT]; // spherical potentials r*V(r), no Y00 factor
      double* zero_potential; // PAW potential shape correction
      double  sigma, sigma_inv; // spread of the SHO projectors and its inverse
      double *hamiltonian, *overlap; // matrices [nSHO*matrix_stride]
      double (*kinetic_energy)[TRU_AND_SMT]; // matrix [nln*nln][TRU_AND_SMT]
      double (*charge_deficit)[TRU_AND_SMT]; // matrix [(1 + ellmax_compensator)*nln*nln][TRU_AND_SMT]
      double  core_charge_deficit; // in units of electrons

      bool gaunt_init;
      std::vector<gaunt_entry_t> gaunt;
      std::vector<int16_t> ln_index_list;
      std::vector<int16_t> lm_index_list;
      std::vector<int16_t> lmn_begin;
      std::vector<int16_t> lmn_end;

      
  public:
    LiveAtom(float const Z_nucleons
            , bool const transfer2valence=true// results look slightly different in the shape of smooth potentials but matrix elements are the same
            , float const ionization=0
            , int const echo=0) : gaunt_init{false} { // constructor
        id = -1; // unset
        Z = Z_nucleons; // convert to float
        if (echo > 0) printf("# LiveAtom with %.1f nucleons, ionization=%g\n", Z, ionization);
        
        rg[TRU] = radial_grid::create_default_radial_grid(Z);
        
        if (0) { // flat copy, true and smooth quantities live on the same radial grid
            rg[SMT] = rg[TRU]; rg[SMT]->memory_owner = false; // avoid double free
        } else { // create a radial grid descriptor which has less points at the origin
            rg[SMT] = radial_grid::create_pseudo_radial_grid(*rg[TRU]);
        } // use the same number of radial grid points for true and smooth quantities
        
        int const nrt = align<2>(rg[TRU]->n),
                  nrs = align<2>(rg[SMT]->n);
        if (echo > 0) printf("# radial grid numbers are %d and %d\n", rg[TRU]->n, rg[SMT]->n);
        if (echo > 0) printf("# radial grid numbers are %d and %d (padded to align)\n", nrt, nrs);

        numax = 3; // 3:up to f-projectors
        if (echo > 0) printf("# projectors and partial waves are expanded up to numax = %d\n", numax);
        ellmax = 0; // should be 2*numax;
        if (echo > 0) printf("# radial density and potentials are expanded up to lmax = %d\n", ellmax);
        ellmax_compensator = 0;
        if (echo > 0) printf("# compensation charges are expanded up to lmax = %d\n", ellmax_compensator);
        r_cut = 2.0; // Bohr
        sigma_compensator = r_cut/std::sqrt(20.); // Bohr
//         sigma_compensator *= 3; // 3x as large as usually taken in GPAW, much smoother
        sigma = 0.5; // Bohr, spread for projectors
        r_match = 9*sigma;
        if (echo > 0) printf("# numbers of projectors ");
        for(int ell = 0; ell <= ELLMAX; ++ell) {
            nn[ell] = std::max(0, (numax + 2 - ell)/2);
            if (echo > 0) printf(" %d", nn[ell]);
        } // ell
        if (echo > 0) printf("\n");

        
        int const nlm = pow2(1 + ellmax);
        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            int const nr = (TRU == ts)? nrt : nrs;
            core_density[ts]   = new double[nr]; // get memory
            potential[ts]      = new double[nr]; // get memory
            full_density[ts]   = new double[nlm*nr]; // get memory
            full_potential[ts] = new double[nlm*nr]; // get memory
        } // true and smooth

        set(core_density[SMT], nrs, 0.0); // init
        set(core_density[TRU], nrt, 0.0); // init
        atom_core::read_Zeff_from_file(potential[TRU], *rg[TRU], Z, "pot/Zeff", -1, echo);

//         // show the loaded Zeff(r)
//         if (echo > 0) {
//            printf("\n## loaded Z_eff(r) function:\n");
//            for(int ir = 0; ir < rg[TRU]->n; ++ir) {
//                printf("%.15g %.15g\n", rg[TRU]->r[ir], -potential[TRU][ir]);
//            } // ir
//         } // echo

        int8_t as_valence[99];
        set(as_valence, 99, (int8_t)-1);
        
        int enn_core_ell[12] = {0,0,0,0, 0,0,0,0, 0,0,0,0};
        auto const r2rho = new double[nrt];
        ncorestates = 20;
        core_state = new core_level_t[ncorestates];
        {   int ics = 0, jcs = -1; float ne = Z - ionization;
            for(int m = 0; m < 8; ++m) { // auxiliary number
                int enn = (m + 1)/2;
                for(int ell = m/2; ell >= 0; --ell) { // angular momentum character
                    ++enn; // principal quantum number
                    for(int jj = 2*ell; jj >= 2*ell; jj -= 2) {
                        auto &cs = core_state[ics]; // abbreviate
                        cs.wave[TRU] = new double[nrt]; // get memory for the true radial function
                        double E = atom_core::guess_energy(Z, enn);
                        set(r2rho, rg[TRU]->n, 0.0);
                        radial_eigensolver::shooting_method(1, *rg[TRU], potential[TRU],
                                 enn, ell, E, cs.wave[TRU], r2rho);
                        cs.energy = E;
                        
                        int const inl = atom_core::nl_index(enn, ell);
                        if (E > -1.0) { // ToDo: make the limit (-1.0 Ha) between core and valence states dynamic
                            as_valence[inl] = ics; // mark as good for the valence band
//                          printf("# as_valence[nl_index(enn=%d, ell=%d) = %d] = %d\n", enn, ell, inl, ics);
                        } // move to the valence band

                        cs.nrn[TRU] = enn - ell - 1; // true number of radial nodes
                        cs.enn = enn;
                        cs.ell = ell;
                        cs.emm = emm_Degenerate;
                        float const max_occ = 2*(jj + 1);
                        cs.spin = spin_Degenerate;

                        float occ = std::min(std::max(0.f, ne), max_occ);
                        cs.occupation = occ;
                        if (occ > 0) {
                            jcs = ics;
                            if (echo > 0) printf("# %s %2d%c%6.1f E = %g %s\n", (as_valence[inl] < 0)?
                                              "core   ":"valence", enn, ellchar[ell], occ, E*eV,_eV);
                        } // occupied
                        if (as_valence[inl] < 0) {
                            enn_core_ell[ell] = std::max(enn, enn_core_ell[ell]);
                        } else {
                            if (transfer2valence) occ = 0;
                        } // not as valence
                        if (occ > 0) {
                            double const norm = occ/dot_product(rg[TRU]->n, r2rho, rg[TRU]->dr);
                            add_product(core_density[TRU], rg[TRU]->n, r2rho, norm);
                        } // occupied
                        ne -= max_occ;
                        
                        ++ics;
                    } // jj
                } // ell
            } // m
            ncorestates = jcs + 1; // correct the number of core states to those occupied
        } // core states
        
        scale(core_density[TRU], rg[TRU]->n, rg[TRU]->rinv); // initial_density produces r^2*rho --> reduce to r*rho
        scale(core_density[TRU], rg[TRU]->n, rg[TRU]->rinv); // initial_density produces r^2*rho --> reduce to   rho
        if (echo > 2) printf("# initial core density has %g electrons\n", dot_product(rg[TRU]->n, core_density[TRU], rg[TRU]->r2dr));

        if (echo > 5) { printf("# enn_core_ell  "); for(int ell = 0; ell <= numax; ++ell) printf(" %d", enn_core_ell[ell]); printf("\n"); }

        nvalencestates = (numax*(numax + 4) + 4)/4; // ToDo: new function in sho_tools::
        valence_state = new valence_level_t[nvalencestates];
        {   int iln = 0;
//          if (echo > 0) printf("# valence "); // no new line, compact list follows
            for(int ell = 0; ell <= numax; ++ell) {
                for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                    auto &vs = valence_state[iln]; // abbreviate
                    int const enn = std::max(ell + 1, enn_core_ell[ell] + 1) + nrn;
//                  if (echo > 0) printf(" %d%c", enn, ellchar[ell]);
                    
                    vs.wave[SMT] = new double[nrs]; // get memory for the smooth radial function
                    vs.wave[TRU] = new double[nrt]; // get memory for the true radial function
                    double E = atom_core::guess_energy(Z, enn);
                    radial_eigensolver::shooting_method(1, *rg[TRU], potential[TRU],
                                          enn, ell, E, vs.wave[TRU]);
                    vs.energy = E;

                    vs.nrn[TRU] = enn - ell - 1; // true number of radial nodes
                    vs.nrn[SMT] = nrn;
                    vs.occupation = 0;
                    {
                        int const inl = atom_core::nl_index(enn, ell);
                        int const ics = as_valence[inl];
//                      printf("# as_valence[nl_index(enn=%d, ell=%d) = %d] = %d\n", enn, ell, inl, ics);
                        if (ics >= 0) { // atomic eigenstate was marked as valence
                            if (transfer2valence) {
                                auto const occ = core_state[ics].occupation;
                                vs.occupation = occ;
                                core_state[ics].occupation = 0;
                                if (occ > 0) printf("# transfer %.1f electrons from %d%c-core state #%d"
                                        " to valence state #%d\n", occ, enn, ellchar[ell], ics, iln);
                            } // transfer2valence
                        } // ics
                    }
                    vs.enn = enn;
                    vs.ell = ell;
                    vs.emm = emm_Degenerate;
                    vs.spin = spin_Degenerate;
                    if (echo > 0) printf("# valence %2d%c%6.1f E = %g %s\n", enn, ellchar[ell], vs.occupation, E*eV,_eV);
                    ++iln;
                } // nrn
            } // ell
//          if (echo > 0) printf("  (%d states)\n", iln);
        } // valence states
        delete[] r2rho;

        int irc = 0; while (rg[SMT]->r[irc] < r_cut) ++irc; // find the radial index of r_cut on the smooth radial grid
        ir_cut[SMT] = irc;
        ir_cut[TRU] = irc + rg[TRU]->n - rg[SMT]->n;
        if (echo > 0) printf("# pseudize the core density at r[%d or %d] = %.6f, requested %.3f %s\n", 
                              ir_cut[SMT], ir_cut[TRU], rg[SMT]->r[irc]*Ang, r_cut*Ang, _Ang);
        assert(rg[SMT]->r[ir_cut[SMT]] == rg[TRU]->r[ir_cut[TRU]]); // should be exactly equal
        int const nlm_aug = pow2(1 + std::max(ellmax, ellmax_compensator));
        aug_density     = new double[nlm_aug*nrs]; // get memory
        int const nlm_cmp = pow2(1 + ellmax_compensator);
        qlm_compensator = new double[nlm_cmp]; // get memory
        int const nln = nvalencestates;
        charge_deficit  = new double[(1 + ellmax_compensator)*nln*nln][TRU_AND_SMT]; // get memory
        kinetic_energy  = new double[nln*nln][TRU_AND_SMT]; // get memory
        zero_potential  = new double[nrs]; // get memory
        
        set(zero_potential, nrs, 0.0); // clear
        set(kinetic_energy[0], nln*nln*2, 0.0); // clear

        int const nSHO = sho_tools::nSHO(numax);
        matrix_stride = align<2>(nSHO); // 2^<2> doubles = 32 Byte alignment
        if (echo > 0) printf("# matrix size for hamiltonian and overlap: dim = %d, stride = %d\n", nSHO, matrix_stride);
        hamiltonian = new double[nSHO*matrix_stride]; // get memory
        overlap     = new double[nSHO*matrix_stride]; // get memory

        unitary_zyx_lmn = new double[nSHO*nSHO];
        {   auto const u = new sho_unitary::Unitary_SHO_Transform<double>(numax);
            auto const stat = u->construct_dense_matrix(unitary_zyx_lmn, numax);
            assert(0 == stat);
        } // scope to fill unitary
        
        
        // return; // early return if we only want to test the occupation configuration
 

        int const mlm = pow2(1 + numax);
        ln_index_list.resize(nSHO);
        lm_index_list.resize(nSHO);
        lmn_begin.resize(mlm);
        lmn_end.resize(mlm);
        get_valence_mapping(ln_index_list.data(), lm_index_list.data(), nSHO, nln, lmn_begin.data(), lmn_end.data(), mlm);
        
        
        int const nr_diff = rg[TRU]->n - rg[SMT]->n;
        {   // construct initial smooth spherical density
            set(potential[SMT], rg[SMT]->n, potential[TRU] + nr_diff); // copy the tail of the spherical part of r*V_tru(r)
            pseudize_function(potential[SMT], rg[SMT], ir_cut[SMT], 2, 1); // replace by a parabola
        }
        
        int const maxit_scf = 1;
        float const mixing = 0.45; // mixing with .45 works well for Cu (Z=29)
        
        for(int scf = 0; scf < maxit_scf; ++scf) {
            if (echo > 1) printf("\n\n# SCF-iteration %d\n\n", scf);
//             update((scf >= maxit_scf - 333)*9); // switch full echo on in the last 3 iterations
            update_density(mixing, echo);
            double* const ves_multipoles = nullptr;
            update_potential(mixing, ves_multipoles, echo);
        } // self-consistency iterations

        // show the smooth and true potential
        if (true && (echo > 0)) {
            printf("\n## spherical parts: r, "
            "Zeff_tru(r), Zeff_smt(r)"
            ", r^2*rho_tru(r), r^2*rho_smt(r)"
            ", zero_potential(r) in a.u."
            ":\n");
            for(int ir = 0; ir < rg[SMT]->n; ir += 1) {
                auto const r = rg[SMT]->r[ir];
                printf("%g %g %g %g %g %g\n", r
//                         , -full_potential[TRU][ir + nr_diff]*Y00*r // for comparison, should be the same as Z_eff(r)
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
            printf("\n\n\n# valence partial waves:\n");
            for(int ir = 0; ir < rg[SMT]->n; ir += 1) {
                auto const r = rg[SMT]->r[ir];
                auto const f = r; // scale with r? scale with Y00?
                printf("%g ", r);
                for(int ics = 0; ics < ncorestates; ++ics) {
                    printf("   %g", core_state[ics].wave[TRU][ir + nr_diff]*f);
                } // ics
                for(int iln = 0; iln < nvalencestates; ++iln) {
                    printf("   %g %g", valence_state[iln].wave[TRU][ir + nr_diff]*f
                                     , valence_state[iln].wave[SMT][ir]*f);
                } // iln
                printf("\n");
            } // ir
            printf("\n\n\n\n");
        } // echo
        
    } // constructor
      
    ~LiveAtom() { // destructor
        for(int ics = 0; ics < ncorestates; ++ics) { delete[] core_state[ics].wave[TRU]; }
        delete[] core_state;
        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            radial_grid::destroy_radial_grid(rg[ts]);
            for(int ivs = 0; ivs < nvalencestates; ++ivs) { delete[] valence_state[ivs].wave[ts]; }
            delete[] core_density[ts];
            delete[] potential[ts];
            delete[] full_density[ts];
            delete[] full_potential[ts];
        } // tru and smt
        delete[] valence_state;
        delete[] aug_density;
        delete[] qlm_compensator;
        delete[] charge_deficit;
        delete[] zero_potential;
        delete[] kinetic_energy;
        delete[] hamiltonian;
        delete[] overlap;
        delete[] unitary_zyx_lmn;
    } // destructor

    status_t initialize_Gaunt() {
      if (gaunt_init) return 0; // success
      auto const stat = angular_grid::create_numerical_Gaunt<6>(&gaunt);
      gaunt_init = (0 == stat);
      return stat;
    } // initialize_Gaunt

    void show_state_analysis(int const echo, radial_grid_t const *rg, double const wave[], 
            int const enn, int const ell, float const occ, double const energy, char const csv, int const ir_cut=-1) {
        if (echo < 1) return;
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

        double charge_outside = 0;
        if (ir_cut >= 0) {
            for(int ir = ir_cut; ir < rg->n; ++ir) {
                double const rho_wf = pow2(wave[ir]);
                double const dV = rg->r2dr[ir];
                charge_outside += rho_wf*dV;
            } // ir
        } // ir_cut

        //  printf("# core    %2d%c %g %g %g %g %g\n", cs.enn, ellchar[cs.ell], stats[0], stats[1], stats[2], stats[3], stats[4]);
        printf("# %s %2d%c%6.1f E=%16.6f %s  <r>=%g rms=%g %s <r^-1>=%g %s q_out=%.3g e\n", 
               ('c' == csv)?"core   ":"valence", enn, ellchar[ell], occ, energy*eV,_eV, stats[2]/stats[1]*Ang, 
               std::sqrt(std::max(0., stats[3]/stats[1]))*Ang,_Ang, stats[4]/stats[1]*eV,_eV,
               charge_outside/stats[1]);
    } // show_state_analysis

    void update_core_states(float const mixing, int const echo=0) {
        // core states are feeling the spherical part of the hamiltonian only
        int const nr = rg[TRU]->n;
        auto const r2rho = new double[nr];
        auto const new_r2core_density = new double[nr];
        set(new_r2core_density, nr, 0.0);
        double nelectrons = 0;
        for(int ics = 0; ics < ncorestates; ++ics) {
            auto &cs = core_state[ics]; // abbreviate
            int constexpr SRA = 1;
            radial_eigensolver::shooting_method(SRA, *rg[TRU], potential[TRU], cs.enn, cs.ell, cs.energy, cs.wave[TRU], r2rho);
            auto const norm = dot_product(nr, r2rho, rg[TRU]->dr);
            auto const norm_factor = (norm > 0)? 1./std::sqrt(norm) : 0;
            auto const scal = pow2(norm_factor)*cs.occupation; // scaling factor for the density contribution of this state
            nelectrons += cs.occupation;
            // transform r*wave(r) as produced by the radial_eigensolver to wave(r)
            // and normalize the core level wave function to one
            scale(cs.wave[TRU], nr, rg[TRU]->rinv, norm_factor);
            add_product(new_r2core_density, nr, r2rho, scal);
            show_state_analysis(echo, rg[TRU], cs.wave[TRU], cs.enn, cs.ell, cs.occupation, cs.energy, 'c', ir_cut[TRU]);
        } // ics
        delete[] r2rho;

        // report integrals
        auto const old_core_charge = dot_product(nr, rg[TRU]->r2dr, core_density[TRU]);
        auto const new_core_charge = dot_product(nr, rg[TRU]->dr, new_r2core_density);
        if (echo > 0) printf("# expect a core density with %g electrons\n", nelectrons);
        if (echo > 0) printf("# previous core density has %g electrons\n", old_core_charge);
        if (echo > 0) printf("# new core density has %g electrons\n",      new_core_charge);
        double mix_new = mixing, mix_old = 1 - mixing;
        // can we rescale the mixing coefficients such that the desired number of core electrons comes out?
        auto const mixed_charge = mix_old*old_core_charge + mix_new*new_core_charge;
        if (mixed_charge != 0) {
            auto const rescale = nelectrons/mixed_charge;
            mix_old *= rescale;
            mix_new *= rescale;
        } // rescale

        double core_density_change = 0, core_density_change2 = 0, core_nuclear_energy = 0;
        for(int ir = 0; ir < nr; ++ir) {
            auto const new_rho = new_r2core_density[ir]*pow2(rg[TRU]->rinv[ir]); // *r^{-2}
            core_density_change  += std::abs(new_rho - core_density[TRU][ir])*rg[TRU]->r2dr[ir];
            core_density_change2 +=     pow2(new_rho - core_density[TRU][ir])*rg[TRU]->r2dr[ir];
            core_nuclear_energy  +=         (new_rho - core_density[TRU][ir])*rg[TRU]->rdr[ir]; // Coulomb integral change
            core_density[TRU][ir] = mix_new*new_rho + mix_old*core_density[TRU][ir];
        } // ir
        core_nuclear_energy *= -Z;
        if (echo > 0) printf("# core density change %g e (rms %g e) energy change %g %s\n",
            core_density_change, std::sqrt(std::max(0.0, core_density_change2)), core_nuclear_energy*eV,_eV);
        delete[] new_r2core_density;

        { // scope: pseudize the core density
            int const nrs = rg[SMT]->n;
            int const nr_diff = nr - nrs;
            // copy the tail of the true core density into the smooth core density
            set(core_density[SMT], nrs, core_density[TRU] + nr_diff);

            auto const stat = pseudize_function(core_density[SMT], rg[SMT], ir_cut[SMT], 3); // 3: use r^0, r^2 and r^4
            // alternatively, pseudize_function(core_density[SMT], rg[SMT], ir_cut[SMT], 3, 2); // 3, 2: use r^2, r^4 and r^6
            if (stat && (echo > 0)) printf("# %s Matching procedure for the smooth core density failed! info = %d\n", __func__, stat);

            if (0) { // plot the core densities
                printf("\n## radius, smooth core density, true core density:\n");
                for(int ir = 0; ir < nrs; ir += 2) {
                    printf("%g %g %g\n", rg[SMT]->r[ir], core_density[SMT][ir]
                                                       , core_density[TRU][ir + nr_diff]);
                } // ir backwards
                printf("\n\n");
            } // plot

            // report integrals
            auto const tru_core_charge = dot_product(rg[TRU]->n, rg[TRU]->r2dr, core_density[TRU]);
            auto const smt_core_charge = dot_product(rg[SMT]->n, rg[SMT]->r2dr, core_density[SMT]);
            if (echo > 0) printf("# true and smooth core density have %g and %g electrons\n", tru_core_charge, smt_core_charge);
            core_charge_deficit = tru_core_charge - smt_core_charge;
        } // scope

    } // update_core_states

    void update_valence_states(int const echo=0) {
        // the basis for valence partial waves is generated from the spherical part of the hamiltonian
//      auto const small_component = new double[rg[TRU]->n];
        int const nr = rg[TRU]->n;
        auto r2rho = new double[nr];
        int const nln = nvalencestates;
        for(int iln = 0; iln < nln; ++iln) {
            auto &vs = valence_state[iln]; // abbreviate
            int constexpr SRA = 1;
            
            // solve for a true partial wave
//          radial_integrator::integrate_outwards<SRA>(*rg[TRU], potential[TRU], vs.ell, vs.energy, vs.wave[TRU], small_component);
            set(vs.wave[TRU], nr, 0.0); // clear

            // solve for a valence eigenstate
            radial_eigensolver::shooting_method(SRA, *rg[TRU], potential[TRU], vs.enn, vs.ell, vs.energy, vs.wave[TRU], r2rho);
            // normalize the partial waves
            auto const norm_wf2 = dot_product(nr, r2rho, rg[TRU]->dr);
            auto const norm_factor = 1./std::sqrt(norm_wf2);
            scale(vs.wave[TRU], nr, rg[TRU]->rinv, norm_factor); // transform r*wave(r) as produced by the radial_eigensolver to wave(r)
            
            auto const tru_norm = dot_product(ir_cut[TRU], r2rho, rg[TRU]->dr)/norm_wf2;
            auto const work = r2rho;
            scale(work, nr, potential[TRU]);
            auto const tru_potential_E = dot_product(ir_cut[TRU], work, rg[TRU]->dr)/norm_wf2;
            auto const tru_kinetic_E = vs.energy*tru_norm - tru_potential_E; // kinetic energy contribution up to r_cut
            
//          if (echo > 1) printf("# valence %2d%c%6.1f E=%16.6f %s\n", vs.enn, ellchar[vs.ell], vs.occupation, vs.energy*eV,_eV);
            show_state_analysis(echo, rg[TRU], vs.wave[TRU], vs.enn, vs.ell, vs.occupation, vs.energy, 'v', ir_cut[TRU]);

            int const nr_diff = rg[TRU]->n - rg[SMT]->n;
            // construct a simple pseudo wave
            int const n_poly = 3; // number of even-order polynomial terms used
            double coeff[4], smt_kinetic_E = 0, smt_norm = 0;
            set(coeff, 4, 0.0);
            set(vs.wave[SMT], rg[SMT]->n, vs.wave[TRU] + nr_diff); // workaround: simply take the tail of the true wave
            auto const stat = pseudize_function(vs.wave[SMT], rg[SMT], ir_cut[SMT], n_poly, vs.ell, coeff);
            if (stat) {
                if (echo > 0) printf("# %s Matching procedure for the smooth %d%c-valence state failed! info = %d\n", 
                                              __func__, vs.enn, ellchar[vs.ell], stat);
            } else {
                if (echo > 0) printf("# %s Matching of smooth %d%c-valence state with polynomial r^%d*(%g + r^2* %g + r^4* %g + r^6* %g)\n", 
                                              __func__, vs.enn, ellchar[vs.ell], vs.ell, coeff[0], coeff[1], coeff[2], coeff[3]);
                // psi(r) = \sum_i c_i r^{\ell + 2i}
                // \hat T_radial psi = -\frac 12 1/r d^2/dr^2 (r*psi) = -\frac 12 \sum_i c_i (\ell + 2i + 1) (\ell + 2i) r^{\ell + 2i - 2}
                // for ell=0 the i=0 term cancels with the centrifugal potential \ell(\ell + 1)r^{-2}
                // then \int_0^R dr r^2 psi(r) \hat T psi(r) 
                //    = -\frac 12 \sum_i c_i \sum_j c_j \left\[ (\ell + 2j + 1) (\ell + 2j) - (\ell + 1) \ell \right\] R^{2\ell + 2i + 2j + 1}/{2\ell + 2i + 2j + 1}
                //    =
                { 
                    double const rcut = rg[SMT]->r[ir_cut[SMT]];
                    smt_kinetic_E = 0; smt_norm = 0;
                    for(int i = 0; i < n_poly; ++i) {
                        for(int j = 0; j < n_poly; ++j) {
                            smt_kinetic_E += ((vs.ell + 2*j + 1)*(vs.ell + 2*j) - (vs.ell + 1)*vs.ell)*
                                             coeff[i]*coeff[j]*intpow(rcut, 2*vs.ell + 2*i + 2*j + 1)/(2*vs.ell + 2*i + 2*j + 1.);
                            smt_norm      += coeff[i]*coeff[j]*intpow(rcut, 2*vs.ell + 2*i + 2*j + 3)/(2*vs.ell + 2*i + 2*j + 3.);
                        } // j
                    } // i
                    smt_kinetic_E *= -0.5; // prefactor for the kinetic energy in Hartree atomic units
                }
            } // stat
            
            product(work, rg[SMT]->n, vs.wave[SMT], vs.wave[SMT]);
            double const smt_norm_numerical = dot_product(ir_cut[SMT], work, rg[SMT]->r2dr);
            if (echo > 0) printf("# %s smooth %d%c-valence state true norm %g, smooth %g and %g (smooth numerically)\n", 
               __func__, vs.enn, ellchar[vs.ell], tru_norm, smt_norm, smt_norm_numerical);
            kinetic_energy[nln*iln + iln][TRU] = tru_kinetic_E; // set diagonal entry, offdiagonals still missing
            kinetic_energy[nln*iln + iln][SMT] = smt_kinetic_E; // set diagonal entry, offdiagonals still missing
            if (echo > 0) printf("# %s %d%c-valence state has kinetic energy %g and (smooth) %g %s, norm deficit = %g electrons\n", 
               __func__, vs.enn, ellchar[vs.ell], tru_kinetic_E*eV, smt_kinetic_E*eV,_eV, tru_norm - smt_norm);
           
            // ToDo: solve for partial waves at the same energy, match and establish dual orthgonality with SHO projectors
        } // iln
        delete[] r2rho;
    } // update_valence_states

    void update_charge_deficit(int const echo=0) {
        int const nln = nvalencestates;
        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            int const nr = rg[ts]->n;
            auto const rl = new double[nr];
            auto const wave_r2rl_dr = new double[nr];
            if (echo > 1) printf("\n# charges for %s partial waves\n", (TRU==ts)?"true":"smooth");
            for(int ell = 0; ell <= ellmax_compensator; ++ell) { // loop-carried dependency on rl, run forward, run serial!
                if (echo > 1) printf("# charges for ell=%d, jln = 0, 1, ...\n", ell);
                if (0 == ell) {
                    set(rl, nr, 1.0); // init as r^0
                } else {
                    scale(rl, nr, rg[ts]->r); // create r^{\ell}
                } // 0 == ell
                for(int iln = 0; iln < nln; ++iln) {
                    if (echo > 1) printf("# iln = %d ", iln);
                    auto const wave_i = valence_state[iln].wave[ts];
                    product(wave_r2rl_dr, nr, wave_i, rl, rg[ts]->r2dr); // product of three arrays
                    for(int jln = 0; jln < nln; ++jln) {
                        auto const wave_j = valence_state[jln].wave[ts];
                        auto const cd = dot_product(nr, wave_r2rl_dr, wave_j);
                        charge_deficit[(ell*nln + iln)*nln + jln][ts] = cd;
                        if (echo > 1) printf("\t%10.6f", cd);
                    } // jln
                    if (echo > 1) printf("\n");
                } // iln
                if (echo > 1) printf("\n\n");
            } // ell
            delete[] wave_r2rl_dr;
            delete[] rl;
        } // ts
    } // update_charge_deficit

    template<typename int_t>
    void get_valence_mapping(int_t ln_index_list[], int_t lm_index_list[], int const nlmn, int const nln, 
                             int_t lmn_begin[], int_t lmn_end[], int const mlm, 
                             int const echo=0) {
        for(int lm = 0; lm < mlm; ++lm) lmn_begin[lm] = -1;
        int ilmn = 0;
        for(int ell = 0; ell <= numax; ++ell) {
            int iln_enn[8]; // create a table of iln-indices
            for(int iln = 0; iln < nln; ++iln) {
                if (ell == valence_state[iln].ell) {
                    int const nrn = valence_state[iln].nrn[SMT];
                    iln_enn[nrn] = iln;
                } // ell matches
            } // iln
            for(int emm = -ell; emm <= ell; ++emm) {
                for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                    ln_index_list[ilmn] = iln_enn[nrn]; // valence state index
                    int const lm = solid_harmonics::lm_index(ell, emm);
                    lm_index_list[ilmn] = lm;
                    if (lmn_begin[lm] < 0) lmn_begin[lm] = ilmn; // store the first index of this lm
                                             lmn_end[lm] = ilmn + 1; // store the last index of this lm
                    ++ilmn;
                } // nrn
            } // emm
        } // ell
        assert(nlmn == ilmn);
        
        if (echo > 3) {
            printf("# ln_index_list ");
            for(int i = 0; i < nlmn; ++i) {
                printf("%3d", ln_index_list[i]); 
            }   printf("\n");
            printf("# lmn_begin--lmn_end "); 
            for(int i = 0; i < mlm; ++i) {
//              if (lmn_begin[i] == lmn_end[i] - 1) { 
                if (0) { 
                    printf(" %d", lmn_begin[i]); 
                } else {
                    printf(" %d--%d", lmn_begin[i], lmn_end[i] - 1); 
                }
            }   printf("\n");
        } // echo
        
    } // get_valence_mapping
    
    
    void transform_SHO(double out[], int const out_stride, 
                  double const in[], int const in_stride, 
                  bool const in_Cartesian, double const alpha=1, int const nu_max=-1) {

        int const u_stride = sho_tools::nSHO(numax);
        int const N = (nu_max > -1)? sho_tools::nSHO(nu_max) : u_stride;
        auto const tmp = new double[N*N];
        
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
                        tij += in[nC*in_stride + kC] * unitary_zyx_lmn[kC*u_stride + mR]; // *u
                    } // kC
                    tmp[nC*N + mR] = alpha*tij;
                } // mR
            } // nC
            for(int nR = 0; nR < N; ++nR) {
                for(int mR = 0; mR < N; ++mR) {
                    double tij = 0;
                    for(int kC = 0; kC < N; ++kC) {
                        tij += unitary_zyx_lmn[kC*u_stride + nR] * tmp[kC*N + mR]; // u^T*
                    } // kC
                    out[nR*out_stride + mR] = alpha*tij;
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
                        tij += unitary_zyx_lmn[nC*u_stride + kR] * in[kR*in_stride + mR]; // u*
                    } // kR
                    tmp[nC*N + mR] = alpha*tij;
                } // mR
            } // nC
            for(int nC = 0; nC < N; ++nC) {
                for(int mC = 0; mC < N; ++mC) {
                    double tij = 0;
                    for(int kR = 0; kR < N; ++kR) {
                        tij += tmp[nC*N + kR] * unitary_zyx_lmn[mC*u_stride + kR]; // *u^T
                    } // kR
                    out[nC*out_stride + mC] = alpha*tij;
                } // mC
            } // nC

        } // in_Cartesian
        
// now using BLAS            ToDo!
//         void dgemm_(const char* tA, const char* tB, const int* M, const int* N, const int* K, const double* alpha, 
//               const double* A, const int* sA, const double* B, const int* sB, const double* beta, double* C, const int* sC);
//         performs C[n][m] := beta*C[n][m] + alpha*sum_k B[n][k] * A[k][m];

//         dgemm_(&tn, &nn, &N, &N, &N, &alpha, unitary_zyx_lmn, &N, in, &in_stride, &beta, tmp, &N);
//         dgemm_(&nn, &nt, &N, &N, &N, &alpha, tmp, &N, unitary_zyx_lmn, &N, &beta, out, &out_stride);
        
        delete[] tmp;
    } // transform_SHO
    
    void get_rho_tensor(double rho_tensor[], double const density_matrix[], int const echo=0) {
        int const nSHO = sho_tools::nSHO(numax);
        int const stride = nSHO;
        assert(stride >= nSHO);
        
        initialize_Gaunt();

        int const lmax = std::max(ellmax, ellmax_compensator);
        int const nlm = pow2(1 + lmax);
        int const mlm = pow2(1 + numax);
        int const nln = nvalencestates;
        // ToDo:
        //   transform the density_matrix[izyx*stride + jzyx]
        //   into a radial_density_matrix[ilmn*stride + jlmn]
        //   using the unitary transform from left and right
        auto const radial_density_matrix = new double[nSHO*stride];
        transform_SHO(radial_density_matrix, stride, density_matrix, stride, true);

        if (0) { // debugging
            auto const check_matrix = new double[nSHO*stride];
            transform_SHO(check_matrix, stride, radial_density_matrix, stride, false);
            double d = 0; for(int i = 0; i < nSHO*stride; ++i) d = std::max(d, std::abs(check_matrix[i] - density_matrix[i]));
            printf("# %s found max deviation %.1e when backtransforming the density matrix\n\n", __func__, d);
            assert(d < 1e-9);
            delete[] check_matrix;
        } // debugging


        if (0) {
            printf("# Radial density matrix\n");
            for(int i = 0; i < nSHO; ++i) {
                for(int j = 0; j < nSHO; ++j) {
                    printf("\t%.1f", radial_density_matrix[i*stride + j]);
                } // j
                printf("\n");
            } // i
            printf("\n");
        } // plot
        
        //   Then, contract with the Gaunt tensor over m_1 and m_2
        //   rho_tensor[(lm*nln + iln)*nln + jln] = 
        //     G_{lm l_1m_1 l_2m_2} * radial_density_matrix[il_1m_1n_1*stride + jl_2m_2n_2]

        set(rho_tensor, nlm*nln*nln, 0.0); // clear
        for(auto gnt : gaunt) {
            int const lm = gnt.lm, lm1 = gnt.lm1, lm2 = gnt.lm2; auto G = gnt.G;
            // if (0 == lm) assert(std::abs(G - Y00*(lm1 == lm2)) < 1e-12); // make sure that G_00ij = delta_ij*Y00
            if (0 == lm) G = Y00*(lm1 == lm2); // make sure that G_00ij = delta_ij*Y00
            if ((lm < nlm) && (lm1 < mlm) && (lm2 < mlm)) {
                for(int ilmn = lmn_begin[lm1]; ilmn < lmn_end[lm1]; ++ilmn) {
                    int const iln = ln_index_list[ilmn];
                    for(int jlmn = lmn_begin[lm2]; jlmn < lmn_end[lm2]; ++jlmn) {
                        int const jln = ln_index_list[jlmn];
                        rho_tensor[(lm*nln + iln)*nln + jln] += G * radial_density_matrix[ilmn*stride + jlmn];

//                         auto const rho_ij = rho_tensor[(lm*nln + iln)*nln + jln];
//                         if (std::abs(rho_ij) > 1e-9)
//                             printf("# LINE=%d rho_ij = %g for lm=%d iln=%d jln=%d\n", __LINE__, rho_ij/Y00, lm, iln, jln);
                    } // jlmn
                } // ilmn
            } // limits
        } // gnt
        delete[] radial_density_matrix;

#ifdef FULL_DEBUG        
        for(int lm = 0; lm < nlm; ++lm) {
            for(int iln = 0; iln < nln; ++iln) {
                for(int jln = 0; jln < nln; ++jln) {
                    auto const rho_ij = rho_tensor[(lm*nln + iln)*nln + jln];
                    if (std::abs(rho_ij) > 1e-9)
                        printf("# LINE=%d rho_ij = %g for lm=%d iln=%d jln=%d\n", __LINE__, rho_ij/Y00, lm, iln, jln);
                } // jln
            } // iln
        } // lm
#endif
    } // get_rho_tensor
    
    
    void update_full_density(double const rho_tensor[], int const echo=0) { // density tensor rho_{lm iln jln}
        int const nlm = pow2(1 + ellmax);
        int const nln = nvalencestates;
        
        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            int const nr = rg[ts]->n, mr = align<2>(nr);
            // full_density[ts][nlm*mr]; // memory layout
            for(int lm = 0; lm < nlm; ++lm) {
                if (0 == lm) {
                    set(full_density[ts], nr, core_density[ts], Y00); // needs scaling with Y00 since core_density has a factor 4*pi
                    if (echo > 0) printf("# %s: %s density has %g electrons after adding the core density\n",  __func__, 
                     (TRU == ts)?"true":"smooth", dot_product(nr, full_density[ts], rg[ts]->r2dr)/Y00);
                } else {
                    set(&(full_density[ts][lm*mr]), nr, 0.0); // clear
                }
                for(int iln = 0; iln < nln; ++iln) {
                    auto const wave_i = valence_state[iln].wave[ts];
                    assert(nullptr != wave_i);
                    for(int jln = 0; jln < nln; ++jln) {
                        auto const wave_j = valence_state[jln].wave[ts];
                        assert(nullptr != wave_j);
                        double const rho_ij = rho_tensor[(lm*nln + iln)*nln + jln];
//                         if (std::abs(rho_ij) > 1e-9)
//                             printf("# rho_ij = %g for lm=%d iln=%d jln=%d\n", rho_ij/Y00, lm, iln, jln);
                        add_product(&(full_density[ts][lm*mr]), nr, wave_i, wave_j, rho_ij);
                    } // jln
                } // iln
            } // lm
            if (echo > 0) printf("# %s: %s density has %g electrons after adding the valence density\n",  __func__, 
                    (TRU == ts)?"true":"smooth", dot_product(nr, full_density[ts], rg[ts]->r2dr)/Y00);

        } // true and smooth

        int const nlm_cmp = pow2(1 + ellmax_compensator);
        for(int ell = 0; ell <= ellmax_compensator; ++ell) {
            for(int emm = -ell; emm <= ell; ++emm) {
                int const lm = solid_harmonics::lm_index(ell, emm);
                double rho_lm = 0;
                for(int iln = 0; iln < nln; ++iln) {
                    for(int jln = 0; jln < nln; ++jln) {
                        double const rho_ij = rho_tensor[(lm*nln + iln)*nln + jln];
                        if (std::abs(rho_ij) > 1e-9)
                            printf("# rho_ij = %g for ell=%d emm=%d iln=%d jln=%d\n", rho_ij/Y00, ell, emm, iln, jln);
                        rho_lm += rho_ij * ( charge_deficit[(ell*nln + iln)*nln + jln][TRU]
                                           - charge_deficit[(ell*nln + iln)*nln + jln][SMT] );
                    } // jln
                } // iln
                assert(lm >= 0);
                assert(lm < nlm_cmp);
                qlm_compensator[lm] = rho_lm;
            } // emm
        } // ell
        
        // account for Z protons in the nucleus and the missing charge in the smooth core density
        qlm_compensator[0] += Y00*(core_charge_deficit - Z);
        if (echo > 5) printf("# compensator monopole charge is %g\n", qlm_compensator[0]/Y00);

        int const nlm_aug = pow2(1 + std::max(ellmax, ellmax_compensator));
        // construct the augmented density
        {   int const nr = rg[SMT]->n, mr = align<2>(nr); // on the smooth grid
            set(aug_density, nlm_aug*mr, 0.0); // clear
            set(aug_density, nlm*mr, full_density[SMT]); // copy smooth full_density, need spin summation?
            add_or_project_compensators<0>(aug_density, ellmax_compensator, rg[SMT], qlm_compensator, sigma_compensator);
            double const aug_charge = dot_product(rg[SMT]->n, rg[SMT]->r2dr, aug_density); // only aug_density[0==lm]
            if (echo > 5) printf("# augmented density has %g electrons, zero expected\n", aug_charge/Y00); // this value should be small

            double const tru_charge = dot_product(rg[TRU]->n, rg[TRU]->r2dr, full_density[TRU]); // only full_density[0==lm]
            if (echo > 5) printf("# true density has %g electrons\n", tru_charge/Y00); // this value should be of the order of Z
        } // scope

    } // update_full_density

    
    void update_full_potential(float const mixing, double const ves_multipole[], int const echo=0) {
        int const nlm = pow2(1 + ellmax);
        int const npt = angular_grid::Lebedev_grid_size(ellmax);
        auto vlm = new double[nlm];
        for(int ts = SMT; ts >= TRU; --ts) { // smooth quantities first, so we can determine vlm
            int const nr = rg[ts]->n, mr = align<2>(nr);
            // full_potential[ts][nlm*mr]; // memory layout

            auto const on_grid = new double[2*npt*mr];
            // set(on_grid, npt*mr, 0.0); // clear

            // transform the lm-index into real-space 
            // using an angular grid quadrature, e.g. Lebedev-Laikov grids
            if ((echo > 6) && (SMT == ts)) printf("# local smooth density at origin %g a.u.\n", full_density[ts][0]*Y00);
            angular_grid::transform(on_grid, full_density[ts], mr, ellmax, false);
            // envoke the exchange-correlation potential (acts in place)
//          printf("# envoke the exchange-correlation on angular grid\n");
            for(int ip = 0; ip < npt*mr; ++ip) {
                double const rho = on_grid[ip];
                double vxc = 0, exc = 0;
                exc = exchange_correlation::lda_PZ81_kernel(rho, vxc);
                on_grid[ip] = vxc;
                on_grid[npt*mr + ip] = exc;
            } // ip
            // transform back to lm-index
            angular_grid::transform(full_potential[ts], on_grid, mr, ellmax, true);
            auto const exc_lm = new double[nlm*mr];
            angular_grid::transform(exc_lm, &on_grid[npt*mr], mr, ellmax, true);
            delete[] on_grid;
            if ((echo > 6) && (SMT == ts)) printf("# local smooth exchange-correlation potential at origin is %g %s\n", full_potential[ts][0]*Y00*eV,_eV);
            if (echo > -1) {
                auto const Edc00 = dot_product(nr, full_potential[ts], full_density[ts], rg[ts]->r2dr); // dot_product with diagonal metric
                printf("# double counting correction  in %s 00 channel %.12g %s\n", (TRU == ts)?"true":"smooth", Edc00*eV,_eV);
                auto const Exc00 = dot_product(nr, exc_lm, full_density[ts], rg[ts]->r2dr); // dot_product with diagonal metric
                printf("# exchange-correlation energy in %s 00 channel %.12g %s\n", (TRU == ts)?"true":"smooth", Exc00*eV,_eV);
            } // echo
            delete[] exc_lm;

            // solve electrostatics inside the spheres
            auto   const Ves = new double[nlm*mr];
            double const q_nucleus = (TRU == ts) ? -Z*Y00 : 0; // Z = number of protons in the nucleus
            auto   const rho_aug   = (TRU == ts) ? full_density[TRU] : aug_density;
            // solve electrostatics with singularity (q_nucleus) // but no outer boundary conditions (v_lm)
            radial_potential::Hartree_potential(Ves, *rg[ts], rho_aug, mr, ellmax, q_nucleus);

            if (SMT == ts) {
                add_or_project_compensators<1>(vlm, ellmax_compensator, rg[SMT], Ves, sigma_compensator); // project to compensators
                if (echo > -1) printf("# inner integral between normalized compensator and smooth Ves(r) = %g %s\n", vlm[0]*Y00*eV,_eV);
                // this seems to high by a factor sqrt(4*pi), ToDo: find out why
//                 scale(vlm, nlm, Y00);
                
                // but the solution of the 3D problem found that these integrals should be v_lm, therefore
                if (nullptr == ves_multipole) {
                    set(vlm, nlm, 0.); // no correction of the electrostatic potential heights for isolated atoms
                } else {
                    if (echo > -1) printf("# v_00 found %g but expected %g %s\n", vlm[0]*Y00*eV, ves_multipole[0]*Y00*eV,_eV);
                    scale(vlm, nlm, -1.); add_product(vlm, nlm, ves_multipole, 1.); // vlm := ves_multipole - vlm
                } // no ves_multipole given
            } // smooth only
            
            if (echo > -1) {
                if (SMT == ts) printf("# local smooth electrostatic potential at origin is %g %s\n", Ves[0]*Y00*eV,_eV);
            }

//             correct_multipole_shift(Ves, ellmax_compensator, rg[ts], vlm); // correct heights of the electrostatic potentials
            add_or_project_compensators<2>(Ves, ellmax_compensator, rg[ts], vlm, sigma_compensator); // has the same effect as correct_multipole_shift

            if (SMT == ts) {   // debug: project again to see if the correction worked out
                double v00;
                add_or_project_compensators<1>(&v00, 0, rg[SMT], Ves, sigma_compensator); // project to compensators
                if (echo > -1) printf("# after correction v_00 is %g %s\n", v00*Y00*eV,_eV);
            }

            add_product(full_potential[ts], nlm*mr, Ves, 1.0); // add the electrostatic potential, scale_factor=1.0
            if (echo > -1) {
                if (SMT == ts) printf("# local smooth electrostatic potential at origin is %g %s\n", Ves[0]*Y00*eV,_eV);
                if (TRU == ts) printf("# local true electrostatic potential*r at origin is %g (should match -Z=%.1f)\n", 
                                          Ves[1]*(rg[TRU]->r[1])*Y00, -Z);
                if (SMT == ts) {
                    printf("\n## local smooth electrostatic potential and augmented density in a.u.:\n");
                    for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                        printf("%g %g %g\n", rg[SMT]->r[ir], Ves[ir]*Y00, aug_density[ir]*Y00);
                    }   printf("\n\n");
                } // smooth
            } // echo
            delete [] Ves;
        } // true and smooth
        if (echo > 6) printf("# local smooth augmented density at origin is %g a.u.\n", aug_density[0]*Y00);

        // construct the zero_potential V_bar
        auto const V_smt = new double[rg[SMT]->n];
        int const nr_diff = rg[TRU]->n - rg[SMT]->n;
        set(V_smt, rg[SMT]->n, full_potential[TRU] + nr_diff); // copy the tail of the spherical part of the true potential
        set(zero_potential, rg[SMT]->n, 0.0); // init zero
        auto const df = Y00*eV; assert(df > 0); // display factor
        if (echo > 5) printf("# %s: match local potential to parabola at R_cut = %g %s, V_tru(R_cut) = %g %s\n", 
                    __func__, rg[SMT]->r[ir_cut[SMT]]*Ang, _Ang, full_potential[TRU][0 + ir_cut[TRU]]*df, _eV);
        auto const stat = pseudize_function(V_smt, rg[SMT], ir_cut[SMT], 2);
        if (stat) {
            if (echo > 0) printf("# %s: Matching procedure for the potential parabola failed! info = %d\n", __func__, stat);
        } else {
//             if (echo > -1) printf("# local smooth zero_potential:\n");
            for(int ir = 0; ir < rg[SMT]->n; ++ir) {
                zero_potential[ir] = V_smt[ir] - full_potential[SMT][0 + ir];
//                 if (echo > -1) printf("%g %g\n", rg[SMT]->r[ir], zero_potential[ir]*Y00);
            } // ir
//             if (echo > -1) printf("\n\n");
            if (echo > 5) printf("# %s: potential parabola: V_smt(0) = %g, V_smt(R_cut) = %g %s\n", 
                                  __func__, V_smt[0]*df, V_smt[ir_cut[SMT]]*df, _eV);
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
            if (echo > 5) printf("# %s: zero potential statistics = %g %g %g %s\n", 
                        __func__, Vint/vol*eV, r1Vint/(vol*r_cut)*eV, r2Vint/(vol*pow2(r_cut))*eV, _eV);
                // these numbers should be small since they indicate that V_bar is localized inside the sphere
                // and how much V_smt deviates from V_tru ouside the sphere
        } // pseudization successful
        if (echo > 5) printf("# %s: zero potential: V_bar(0) = %g, V_bar(R_cut) = %g, V_bar(R_max) = %g %s\n",
            __func__, zero_potential[0]*df, zero_potential[ir_cut[SMT]]*df, zero_potential[rg[SMT]->n - 1]*df, _eV);
        delete [] V_smt;

        // add spherical zero potential for SMT==ts and 0==lm
        add_product(full_potential[SMT], rg[SMT]->n, zero_potential, 1.0);

        // feed back spherical part of the true potential into the spherical true potential r*V 
        // which determines core states and true partial waves
        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            scale(potential[ts], rg[ts]->n, 1. - mixing);
            add_product(potential[ts], rg[ts]->n, full_potential[ts] + 0, rg[ts]->r, Y00*mixing);
        } // ts true and smooth

//         // test: use the spherical routines from atom_core::rad_pot(output=r*V(r), input=rho(r)*4pi)
//         auto const rho4pi = new double[rg[TRU]->n];
//         set(rho4pi, rg[TRU]->n, full_density[TRU], 1./Y00);
//         printf("\n# WARNING: use rad_pot to construct the r*V_tru(r)\n\n");
//         atom_core::rad_pot(potential[TRU], *rg[TRU], rho4pi, Z);
//         delete[] rho4pi;

    } // update_full_potential

    void update_matrix_elements(int const echo=0) {
        int const nlm = pow2(1 + ellmax);
        int const mlm = pow2(1 + numax);
        int const nln = nvalencestates;
        int const nSHO = sho_tools::nSHO(numax);
        int const nlmn = nSHO;
        
        initialize_Gaunt();
        
        // first generate the matrix elemnts in the valence basis
        //    overlap[iln*nln + jln] and potential[(lm*nln + iln)*nln + jln]
        // then, generate the matrix elements in the radial representation
        //    overlap[ilmn*nlmn + jlmn] and hamiltonian[ilmn*nlmn + jlmn]
        //    where hamiltonian[ilmn*nlmn + jlmn] = kinetic_energy_deficit_{iln jln} 
        //            + sum_lm Gaunt_{lm ilm jlm} * potential[(lm*nln + iln)*nln + jln]
        // then, transform the matrix elements into the Cartesian representation using sho_unitary
        //    overlap[iSHO*nSHO + jSHO] and hamiltonian[iSHO*nSHO + jSHO]
        
        auto const potential_ln = new double[nlm*nln*nln][TRU_AND_SMT];
        for(int ts = TRU; ts < TRU_AND_SMT; ++ts) {
            int const nr = rg[ts]->n, mr = align<2>(nr);
            auto const wave_pot_r2dr = new double[mr];
            for(int ell = 0; ell <= ellmax; ++ell) {
                for(int emm = -ell; emm <= ell; ++emm) {
                    int const lm = solid_harmonics::lm_index(ell, emm);
                    assert(lm < nlm);
                    for(int iln = 0; iln < nln; ++iln) {
                        auto const wave_i = valence_state[iln].wave[ts];
                        product(wave_pot_r2dr, nr, wave_i, &(full_potential[ts][lm*mr]), rg[ts]->r2dr);
                        for(int jln = 0; jln < nln; ++jln) {
                            auto const wave_j = valence_state[jln].wave[ts];
                            potential_ln[(lm*nln + iln)*nln + jln][ts] =
                                dot_product(nr, wave_pot_r2dr, wave_j);
                        } // jln
                    } // iln
                } // emm
            } // ell
            delete[] wave_pot_r2dr;
        } // ts: true and smooth

        auto const hamiltonian_lmn = new double[nlmn*nlmn];
        auto const overlap_lmn     = new double[nlmn*nlmn];
        set(hamiltonian_lmn, nlmn*nlmn, 0.0); // clear
        set(overlap_lmn,     nlmn*nlmn, 0.0); // clear
        for(auto gnt : gaunt) {
            int const lm = gnt.lm, lm1 = gnt.lm1, lm2 = gnt.lm2; auto G = gnt.G;
            if (0 == lm) G = Y00*(lm1 == lm2); // make sure that G_00ij = delta_ij*Y00

            if (lm1 < mlm && lm2 < mlm) {
                if (lm < nlm) {
                    for(int ilmn = lmn_begin[lm1]; ilmn < lmn_end[lm1]; ++ilmn) {
                        int const iln = ln_index_list[ilmn];
                        for(int jlmn = lmn_begin[lm2]; jlmn < lmn_end[lm2]; ++jlmn) {
                            int const jln = ln_index_list[jlmn];
                            hamiltonian_lmn[ilmn*nlmn + jlmn] +=
                              G * ( potential_ln[(lm*nln + iln)*nln + jln][TRU] 
                                  - potential_ln[(lm*nln + iln)*nln + jln][SMT] );
                        } // jlmn
                    } // ilmn
                } // lm
                if (0 == lm) {  
                    int const ell = 0;
    //              printf("# nln = %d\n", nln);
                    for(int ilmn = lmn_begin[lm1]; ilmn < lmn_end[lm1]; ++ilmn) {
                        int const iln = ln_index_list[ilmn];
                        for(int jlmn = lmn_begin[lm2]; jlmn < lmn_end[lm2]; ++jlmn) {
                            int const jln = ln_index_list[jlmn];
                            overlap_lmn[ilmn*nlmn + jlmn] += 
                              G * ( charge_deficit[(ell*nln + iln)*nln + jln][TRU]
                                  - charge_deficit[(ell*nln + iln)*nln + jln][SMT] );
                        } // jlmn
                    } // ilmn
                } // 0 == lm
            } // limits
        } // gnt
        delete[] potential_ln;

        // add the kinetic_energy deficit to the hamiltonian
        for(int ilmn = 0; ilmn < nlmn; ++ilmn) {
            int const iln = ln_index_list[ilmn];
            int const ilm = lm_index_list[ilmn];
            if (echo > 7) printf("# hamiltonian elements in %s for ilmn=%3d  ", _eV, ilmn);
            for(int jlmn = 0; jlmn < nlmn; ++jlmn) {
                int const jlm = lm_index_list[jlmn];
                if (ilm == jlm) {
                    int const jln = ln_index_list[jlmn];
                    hamiltonian_lmn[ilmn*nlmn + jlmn] += ( kinetic_energy[iln*nln + jln][TRU]
                                                         - kinetic_energy[iln*nln + jln][SMT] )*Y00;
                } // diagonal in lm, offdiagonal in nrn
                if ((echo > 7)) printf(" %g", hamiltonian_lmn[ilmn*nlmn + jlmn]*eV);
            } // jlmn
            if ((echo > 7)) printf("\n");
        } // ilmn

        set(hamiltonian, nSHO*matrix_stride, 0.0); // clear
        set(overlap,     nSHO*matrix_stride, 0.0); // clear
        // Now transform _lmn quantities to Cartesian representations using sho_unitary
        transform_SHO(hamiltonian, matrix_stride, hamiltonian_lmn, nlmn, false);
        transform_SHO(    overlap, matrix_stride,     overlap_lmn, nlmn, false);

        delete[] hamiltonian_lmn;
        delete[] overlap_lmn;
    } // update_matrix_elements
    
    void set_pure_density_matrix(double density_matrix[], float const occ_spdf[4]=nullptr, int const echo=0) {
        float occ[12] = {0,0,0,0, 0,0,0,0, 0,0,0,0}; if (occ_spdf) std::copy(occ_spdf, 4+occ_spdf, occ);
        int const nSHO = sho_tools::nSHO(numax);
        auto const radial_density_matrix = new double[nSHO*nSHO];
        set(radial_density_matrix, nSHO*nSHO, 0.0); // clear
        for(int ell = 0; ell <= numax; ++ell) {
            for(int emm = -ell; emm <= ell; ++emm) {
                for(int enn = 0; enn <= (numax - ell)/2; ++enn) {
                    int const i = sho_tools::lmn_index(numax, ell, emm, enn);
                    assert(nSHO > i);
                    if (0 == enn) radial_density_matrix[i*nSHO + i] = occ[ell]/(2*ell + 1.); // diagonal entry
                } // enn
            } // emm
        } // ell
        transform_SHO(density_matrix, nSHO, radial_density_matrix, nSHO, false);
        if (echo > 7) {
            printf("# radial density matrix\n");
            for(int i = 0; i < nSHO; ++i) {
                for(int j = 0; j < nSHO; ++j) printf("\t%.1f", radial_density_matrix[i*nSHO + j]);
                printf("\n");
            } // i
            printf("\n# Cartesian density matrix\n");
            for(int i = 0; i < nSHO; ++i) {
                for(int j = 0; j < nSHO; ++j) printf("\t%.1f", density_matrix[i*nSHO + j]);
                printf("\n");
            } // i
        } // echo
        delete[] radial_density_matrix;
    } // set_pure_density_matrix

    void update_density(float const mixing, int const echo=0) {
//         if (echo > 2) printf("\n# %s\n", __func__);
        update_core_states(mixing, echo);
        update_valence_states(echo); // create new partial waves for the valence description
        update_charge_deficit(echo); // update quantities derived from the partial waves
        int const nSHO = sho_tools::nSHO(numax);
        auto const density_matrix = new double[nSHO*nSHO];
        {
            float occ[] = {0,0,0,0};
            for(int ivs = 0; ivs < nvalencestates; ++ivs) {
                int const ell = valence_state[ivs].ell;
                if ((ell < 4) && (0 == valence_state[ivs].nrn[SMT])) {
                    occ[ell] = valence_state[ivs].occupation;
                    if ((occ[ell] > 0) && (echo > 1))
                    printf("# Set valence density matrix to be a pure %d%c-state with occupation %.3f\n", 
                        valence_state[ivs].enn, ellchar[ell], occ[ell]);
                } // matching enn-ell quantum numbers?
            } // ivs
            set_pure_density_matrix(density_matrix, occ);
        }
        int const lmax = std::max(ellmax, ellmax_compensator);
        auto const rho_tensor = new double[pow2(1 + lmax)*pow2(nvalencestates)];
        set(rho_tensor, pow2(1 + lmax)*pow2(nvalencestates), 0.0);
        get_rho_tensor(rho_tensor, density_matrix, echo);
        delete[] density_matrix;
        update_full_density(rho_tensor, echo);
        delete[] rho_tensor;
    } // update_density

    // ==============
    // between update_density and update_potential we need to
    // export qlm_compensator, solve 3D electrostatic problem, return here with ves_multipoles
    // ==============

    void update_potential(float const mixing, double const ves_multipoles[], int const echo=0) {
        if (echo > 2) printf("\n# %s\n", __func__);
        update_full_potential(mixing, ves_multipoles, echo);
        update_matrix_elements(echo); // this line does not compile with icpc (ICC) 19.0.2.187 20190117
    } // update_potential

    status_t get_smooth_core_density(double rho[], float const ar2, int const nr2, int const echo=1) {
        if (echo > 7) printf("# %s call transform_to_r2_grid(%p, %.1f, %d, core_density=%p, rg=%p)\n", 
                              __func__, (void*)rho, ar2, nr2, (void*)core_density[SMT], (void*)rg[SMT]);
        return bessel_transform::transform_to_r2_grid(rho, ar2, nr2, core_density[SMT], *rg[SMT], echo);
    } // get_smooth_core_density

    radial_grid_t* get_smooth_radial_grid(int const echo=0) {
        return rg[SMT];
    } // get_smooth_radial_grid
    
  }; // class LiveAtom


namespace single_atom {
  // this module allows to compute a full PAW calculation on a radial grid, 
  // i.e. we assume radial symmetry of density and potential and thus 
  // we can solve all equations on 1D radial grids.
  
  // Potential radial basis functions:
  //  * grid points
  //  * Bessel functions or 
  //  * radial SHO eigenfunctions

  // What is needed?
  //  * PAW construction of smooth partial waves (SHO-scheme), see ReSpawN or juRS
  //  * projection coefficients
  //  * Diagonalizer: LAPACK dsygv.f
  //  * self-consistency, e.g. by potential mixing
  // We have only up to 4 wave functions: s, p, d and f
  // We have only a monopole charge compensator (Gaussian)
  // Maybe we should write a live_atom module first 
  //   (a PAW generator prepared for core level und partial wave update)
  
  status_t update(int const na, float const Za[], float const ion[], 
                  double **rho, radial_grid_t **rg, double *sigma_cmp, double **qlm, double **vlm) {

      static LiveAtom **a=nullptr;
      
      if (nullptr == a) {
          a = new LiveAtom*[na];
          for(int ia = 0; ia < na; ++ia) {
              a[ia] = new LiveAtom(Za[ia], false, ion[ia], 9);
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
              int const nr2 = 1 << 12; float const ar2 = 16.f;
              rho[ia] = new double[nr2];
              a[ia]->get_smooth_core_density(rho[ia], ar2, nr2);
          } // get the smooth core densities
          
          if (nullptr != rg) rg[ia] = a[ia]->get_smooth_radial_grid(); // pointers to smooth radial grids
          
          if (nullptr != sigma_cmp) sigma_cmp[ia] = a[ia]->sigma_compensator; // spreads of the compensators // ToDo: use a getter function
          
          if (nullptr != qlm) set(qlm[ia], 1, a[ia]->qlm_compensator); // copy compensator multipoles
          
          if (nullptr != vlm) a[ia]->update_potential(.5f, vlm[ia], 9); // set electrostatic multipole shifts
          
      } // ia

      return 0;
  } // update
  

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_compensator_normalization(int const echo=5) {
      if (echo > 1) printf("\n# %s: %s\n", __FILE__, __func__);
      auto const rg = radial_grid::create_exponential_radial_grid(512, 2.0);
      int const nr = rg->n, lmax = 0, nlm = pow2(1 + lmax);
      auto const cmp = new double[nr];
      for(double sigma = 0.5; sigma < 2.1; sigma *= 1.1) {
          double qlm[nlm]; set(qlm, nlm, 0.); qlm[0] = 1.;
          set(cmp, nr, 0.); // clear
          add_or_project_compensators<0>(cmp, lmax, rg, qlm, sigma, 0); // add normalized compensator
//        add_or_project_compensators<1>(qlm, lmax, rg, cmp, sigma, 0); // project
//        if (echo > 0) printf("# %s: square-norm of normalized compensator with sigma = %g is %g\n", __func__, sigma, qlm[0]);
          add_or_project_compensators<3>(qlm, lmax, rg, cmp, sigma, 0); // test normalization
          if (echo > 0) printf("# %s: normalization of compensators with sigma = %g is %g\n", __func__, sigma, qlm[0]);
      } // sigma
      return 0;
  } // test_compensator_normalization

  int test(int const echo=9) {
    if (echo > 0) printf("\n# %s: new struct LiveAtom has size %ld Byte\n\n", __FILE__, sizeof(LiveAtom));
//     for(int Z = 0; Z <= 109; ++Z) { // all elements
    // for(int Z = 109; Z >= 0; --Z) { // all elements backwards
        // if (echo > 1) printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");      
//     { int const Z = 1; // 1:hydrogen
//     { int const Z = 2; // 2:helium
//     { int const Z = 29; // 29:copper
//     { int const Z = 47; // 47:silver
//     { int const Z = 79; // 79:gold
    { int const Z = 13; // 13:aluminum
        if (echo > 1) printf("\n# Z = %d\n", Z);
        LiveAtom a(Z, false, 0.f, echo); // envoke constructor
    } // Z
    return 0;
  } // test

  status_t all_tests() {
    auto status = 0;
    status += test_compensator_normalization();
    status += test();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace single_atom
