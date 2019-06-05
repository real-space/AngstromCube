#include <cstdio> // printf
#include <cstdlib> // abs
#include <cmath> // sqrt
#include <cassert> // assert
#include <algorithm> // max

#include "single_atom.hxx"

#include "radial_grid.hxx" // create_exponential_radial_grid, destroy_radial_grid
#include "radial_eigensolver.hxx" // shooting_method
#include "radial_integrator.hxx" // integrate_outwards
#include "atom_core.hxx" // dot_product, initial_density
using namespace atom_core;

#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t, emm_Degenerate, spin_QN_t, spin_Degenerate
#include "output_units.h" // eV, _eV

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

  int constexpr TRU=0, SMT=1, TRU_AND_SMT=2;
  int constexpr CORE=0, VALENCE=1;
  int constexpr ELLMAX=7;
  char const *const ellchar = "spdfghijkl";
  
//   int constexpr NUMCORESTATES=20; // 20 are ok, 32 are enough if spin-orbit-interaction is on
//   int constexpr NUMVALENCESTATES=(ELLMAX*(ELLMAX + 4) + 4)/4;

  
  template<int Pseudo>
  struct energy_level {
      double* wave[Pseudo]; // for valence states points to the true and smooth partial waves
      double energy; // energy level in Hartree atomic units
      float occupation; // occupation number
      enn_QN_t enn; // main quantum_number
      ell_QN_t ell; // angular momentum quantum_number
      emm_QN_t emm; // usually emm == emm_Degenerate
      spin_QN_t spin; // usually spin == spin_Degenerate
  };

  typedef struct energy_level<1+CORE>       core_level_t;
  typedef struct energy_level<1+VALENCE> valence_level_t;

  template<int bits> inline size_t align(size_t const n) { return (((n - 1) >> bits) + 1) << bits; }
  
  class LiveAtom {
  public:
//    double* mem; // memory to which all radial function pointers and matrices point
      
      // general config
      int32_t id; // global atom identifyer
      float Z; // number of nucleons in the core
      radial_grid_t* rg[TRU_AND_SMT]; // radial grid descriptor for the true and smooth grid, SMT may point to TRU
      ell_QN_t ellmax; // limit ell for full_potential and full_density
      double r_cut; // classical augmentation radius for potential and core density
      float r_match; // radius for matching of true and smooth partial wave, usually 6--9*sigma
      ell_QN_t numax; // limit of the SHO projector quantum numbers
      uint8_t nn[1+ELLMAX]; // number of projectors and partial waves used in each ell-channel
      ell_QN_t ellmax_compensator; // limit ell for the charge deficit compensators
      double sigma_compensator; // Gaussian spread for the charge deficit compensators
      double* rho_compensator; // coefficients for the charge deficit compensators, (1+ellmax_compensator)^2
      int matrix_stride; // stride for the matrices hamiltonian and overlap, must be >= ((1+numax)*(2+numax)*(3+numax))/6
      int8_t ncorestates; // for emm-Degenerate representations, 20 (or 32 with spin-oribit) core states are maximum
      int8_t nvalencestates; // for emm-Degenerate (numax*(numax + 4) + 4)/4 is a good choice;
      int8_t nspins; // 1 or 2 or 4

      // spin-resolved members of live_atom
      core_level_t* core_state;
      valence_level_t* valence_state;
      double* core_density[TRU_AND_SMT]; // spherical core density
      double* full_density[TRU_AND_SMT]; // total density, core + valence, (1+ellmax)^2 radial functions
      double* full_potential[TRU_AND_SMT]; // (1+ellmax)^2 radial functions
      double* potential[TRU_AND_SMT]; // spherical potential r*V(r), no Y00 factor
      double* bar_potential; // PAW potential shape correction
      double  sigma, sigma_inv; // spread of the SHO projectors and its inverse
      double* hamiltonian, *overlap; // matrices

  public:
    LiveAtom(double const Z_nucleons) {; // constructor
        int constexpr echo = 1;
        id = -1; // unset
        Z = Z_nucleons; // convert to float
        rg[TRU] = radial_grid::create_exponential_radial_grid(250*std::sqrt(Z + 9.));
//      rg[SMT] = rg[TRU]; // flat copy, alternatively, we could create a radial grid descriptor which has less points at the origin:
        rg[SMT] = radial_grid::create_pseudo_radial_grid(*rg[TRU]);
        int const nrt = align<2>(rg[TRU]->n),
                  nrs = align<2>(rg[SMT]->n);
        if (echo > 0) printf("# padded radial grid numbers are %d and %d\n", nrt, nrs);
        numax = 3; // up to f-projectors
        ellmax = 0; // should be 2*numax;
        ellmax_compensator = 0;
        r_cut = 2.0; // Bohr
        sigma = 0.5; // Bohr
        r_match = 9*sigma;
        if (echo > 0) printf("# numbers of projectors ");
        for(int ell = 0; ell <= ELLMAX; ++ell) {
            nn[ell] = std::max(0, (numax + 2 - ell)/2);
            if (echo > 0) printf(" %d", nn[ell]);
        } // ell
        if (echo > 0) printf("\n");

        int enn_core_max[10] = {0,1,2,3,4,5,6,7,8,9};
        ncorestates = 20;
        core_state = new core_level_t[ncorestates];
        {   int ics = 0, jcs = -1; float ne = Z;
            for(int m = 0; m < 8; ++m) { // auxiliary number
                int enn = (m + 1)/2;
                for(int ell = m/2; ell >= 0; --ell) { // angular momentum character
                    ++enn; // principal quantum number
                    for(int jj = 2*ell; jj >= 2*ell; jj -= 2) {
                        float const max_occ = 2*(jj + 1);
                        float const occ = std::min(std::max(0.f, ne), 2.f*(jj + 1));
                        core_state[ics].energy = -.5*(Z/enn)*(Z/enn); // hydrogen like energy levels
                        // warning: this memory is not freed
                        core_state[ics].wave[TRU] = new double[nrt]; // get memory for the (true) radial function
                        core_state[ics].enn = enn;
                        core_state[ics].ell = ell;
                        core_state[ics].emm = emm_Degenerate;
                        core_state[ics].spin = spin_Degenerate;
                        core_state[ics].occupation = occ;
                        if (max_occ == occ) {
                            if (echo > 0) printf("# core    %2d%c%6.1f \n", enn, ellchar[ell], occ);
                            enn_core_max[ell] = std::max(enn_core_max[ell], enn); 
                        }
                        if (occ > 0) jcs = ics;
                        ne -= max_occ;
                        ++ics;
                    } // jj
                } // ell
            } // m
            ncorestates = jcs + 1;
        } // core states
        
        nvalencestates = (numax*(numax + 4) + 4)/4;
        valence_state = new valence_level_t[nvalencestates];
        {   int ivs = 0;
            for(int ell = 0; ell <= numax; ++ell) {
                for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                    int const enn = nrn + enn_core_max[ell] + 1;
                    if (echo > 0) printf("# valence %2d%c\n", enn, ellchar[ell]);
                    valence_state[ivs].energy = -.5*(Z/enn)*(Z/enn); // hydrogen like energy levels
                    // warning: this memory is not freed
                    valence_state[ivs].wave[TRU] = new double[nrt]; // get memory for the true radial function
                    valence_state[ivs].wave[SMT] = new double[nrs]; // get memory for the smooth radial function
                    valence_state[ivs].occupation = 0;
                    valence_state[ivs].enn = enn;
                    valence_state[ivs].ell = ell;
                    valence_state[ivs].emm = emm_Degenerate;
                    valence_state[ivs].spin = spin_Degenerate;
                    ++ivs;
                } // nrn
            } // ell
        } // valence states
        
        int const nlm = (1 + ellmax)*(1 + ellmax);
        for(int ts = TRU; ts < TRU_AND_SMT; ts += (SMT - TRU)) {
            int const nr = (TRU == ts)? nrt : nrs;
            core_density[ts]   = new double[nr]; // get memory
            potential[ts]      = new double[nr]; // get memory
            full_density[ts]   = new double[nlm*nr]; // get memory
            full_potential[ts] = new double[nlm*nr]; // get memory
        } // true and smooth
        bar_potential = new double[nrs]; // get memory

        int const nSHO = ((1 + numax)*(2 + numax)*(3 + numax))/6; // this should be in sho_tools
        matrix_stride = align<2>(nSHO); // 2^<2> doubles = 32 Byte alignment
        hamiltonian = new double[nSHO*matrix_stride]; // get memory
        overlap     = new double[nSHO*matrix_stride]; // get memory

//      for(int ir = 0; ir < nrt; ++ir) { potential[TRU][ir] = -Z; } // unscreened hydrogen-type potential
        initial_density(core_density[TRU], *rg[TRU], Z, 0.0);
        for(int scf = 0; scf < 133; ++scf) {
            rad_pot(potential[TRU], *rg[TRU], core_density[TRU], Z);
            update_states(9);
            printf("\n");
        } // self-consistency iterations

    };
      
    ~LiveAtom() { // destructor
        // warning: cleanup does not cover all allocations
//      if (rg[SMT] != rg[TRU]) radial_grid::destroy_radial_grid(rg[SMT]); // gives memory corruption
        radial_grid::destroy_radial_grid(rg[TRU]);
        delete[] core_state;
        delete[] valence_state;
    };
    
    void get_rho_tensor(double rho_tensor[], double const density_matrix[], int const stride) {
        int const nSHO = ((1 + numax)*(2 + numax)*(3 + numax))/6; // this should be in sho_tools
        assert(stride >= nSHO);
//         int const nlm = (1 + ellmax)*(1 + ellmax);
//         int const nvs = nvalencestates;
        // ToDo:
        //   transform the density_matrix[iSHO*stride + jSHO]
        //   into a radial_density_matrix[inlm*stride + jnlm] 
        //   using the unitary transform from left and right
        //   Then, contract with the Gauntt tensor over m_1 and m_2
        //   rho_tensor[(lm*nvs + ivs)*nvs + jvs] = 
        //     G_{lm l_1m_1 l_2m_2} * density_matrix[in_1l_1m_1*stride + jn_2l_2m_2]
    } // get
    
    void update_full_density(double const rho_tensor[]) { // density tensor rho_{lm ivs jvs}
        int const nlm = (1 + ellmax)*(1 + ellmax);
        int const nvs = nvalencestates;
        for(int ts = TRU; ts < TRU_AND_SMT; ts += (SMT - TRU)) {
            int const nr = rg[ts]->n, mr = align<2>(nr);
            // full_density[ts][nlm*mr]; // memory layout
            for(int lm = 0; lm < nlm; ++lm) {
                if (0 == lm) {
                    for(int ir = 0; ir < mr; ++ir) {
                        full_density[ts][lm*mr + ir] = core_density[ts][ir]; // needs scaling with Y00 here?
                    } // ir
                } else {
                    for(int ir = 0; ir < mr; ++ir) {
                        full_density[ts][lm*mr + ir] = 0; // init clear
                    } // ir
                }
                for(int ivs = 0; ivs < nvs; ++ivs) {
                    auto const wave_i = valence_state[ivs].wave[ts];
                    for(int jvs = 0; jvs < nvs; ++jvs) {
                        auto const wave_j = valence_state[jvs].wave[ts];
                        double const rho_ij = rho_tensor[(lm*nvs + ivs)*nvs + jvs];
                        for(int ir = 0; ir < nr; ++ir) {
                            full_density[ts][lm*mr + ir] += rho_ij * wave_i[ir] * wave_j[ir];
                        } // ir
                    } // jvs
                } // ivs
            } // lm
        } // true and smooth
    } // update

    void update_full_potential(double const q_lm[]) {
//      int const nlm = (1 + ellmax)*(1 + ellmax);
        for(int ts = TRU; ts < TRU_AND_SMT; ts += (SMT - TRU)) {
//          int const nr = rg[ts]->n, mr = align<2>(nr);
            // full_potential[ts][nlm*mr]; // memory layout
            // transform the lm-index into real-space 
            // using an angular grid quadrature, e.g. Lebendev-Laikov grids
            // envoke the exchange-correlation potential there
            // transform back to lm-index
            // add zero potential for SMT==ts and 0==lm
            // add electrostatic boundary elements q_lm*r^\ell
        } // true and smooth
    } // update 
    
    
    void update_core_states(float const mixing, int echo=0) {
        int const nr = rg[TRU]->n;
        auto r2rho = new double[nr];
        auto new_r2core_density = std::vector<double>(nr, 0.0);
        for(int ics = 0; ics < ncorestates; ++ics) {
            int constexpr SRA = 1;
            radial_eigensolver::shooting_method(SRA, *rg[TRU], potential[TRU], core_state[ics].enn, core_state[ics].ell, 
                            core_state[ics].energy, core_state[ics].wave[TRU], r2rho);
            if (echo > 0) printf("# core    %2d%c%6.1f E=%16.6f %s\n", core_state[ics].enn, ellchar[core_state[ics].ell], 
                            core_state[ics].occupation, core_state[ics].energy*eV,_eV);
            auto const norm = dot_product(nr, r2rho, rg[TRU]->dr);
            auto const scal = (norm > 0)? core_state[ics].occupation/norm : 0;
            for(int ir = 0; ir < nr; ++ir) {
                new_r2core_density[ir] += scal*r2rho[ir];
            } // ir
        } // ics
        double core_density_change = 0, core_nuclear_energy = 0;
        for(int ir = 0; ir < nr; ++ir) {
            auto const new_rho = new_r2core_density[ir]*(rg[TRU]->rinv[ir]*rg[TRU]->rinv[ir]); // *r^{-2}
            core_density_change += std::abs(new_rho - core_density[TRU][ir])*rg[TRU]->r2dr[ir];
            core_nuclear_energy +=      -Z*(new_rho - core_density[TRU][ir])*rg[TRU]->rdr[ir];
            core_density[TRU][ir] = mixing*new_rho + (1. - mixing)*core_density[TRU][ir];
        } // ir
        if (echo > 0) printf("# core density change %g e, energy change %g %s\n", core_density_change, core_nuclear_energy*eV,_eV);
    } // update

    void update_valence_states(int echo=0) {
        auto small_component = new double[rg[TRU]->n];
        for(int ivs = 0; ivs < nvalencestates; ++ivs) {
            int constexpr SRA = 1;
            valence_state[ivs].energy = -.25; // random number
            radial_integrator::integrate_outwards<SRA>(*rg[TRU], potential[TRU], valence_state[ivs].ell, valence_state[ivs].energy, 
                                    valence_state[ivs].wave[TRU], small_component);
            if (echo > 0) printf("# valence %2d%c%6.1f E=%16.6f %s\n", valence_state[ivs].enn, ellchar[valence_state[ivs].ell], 
                            valence_state[ivs].occupation, valence_state[ivs].energy*eV,_eV);
        } // ivs
    } // update

    void update_states(int echo=0) {
        update_core_states(.5, echo);
//      update_valence_states(echo);
    } // update

  }; // LiveAtom


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
  

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  int test(int echo=9) {
    if (echo > 0) printf("\n%s: new struct live_atom has size %ld Byte\n\n", __FILE__, sizeof(LiveAtom));
    for(int Z = 29; Z <= 29; ++Z) {
        if (echo > 1) printf("\n# Z = %d\n", Z);      
        LiveAtom a(Z);
    }
    return 0;
  } // test

  status_t all_tests() {
    auto status = 0;
    status += test();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace single_atom
