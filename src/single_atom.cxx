#include <cstdio> // printf
#include <cstdlib> // abs
#include <cmath> // sqrt
#include <algorithm> // max

#include "single_atom.hxx"
#include "radial_grid.hxx" // create_exponential_radial_grid, destroy_radial_grid

#include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t, emm_Degenerate, spin_QN_t, spin_Degenerate

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
      ell_QN_t ellmax_density; // limit ell for the full_potential
      double r_cut; // classical augmentation radius for potential and core density
      float r_match; // radius for matching of true and smooth partial wave, usually 6--9*sigma
      ell_QN_t ellmax_potential; // limit ell for the full_potential
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
      double* full_density[TRU_AND_SMT]; // total density, core + valence, (1+ellmax_density)^2 radial functions
      double* full_potential[TRU_AND_SMT]; // Y00 components not included here, (1+ellmax_potential)^2-1 radial functions
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
        ellmax_density = 0;
        ellmax_potential = 0;
        ellmax_compensator = 0;
        numax = 3; // up to f-projectors
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
        {   int ics = 0; float ne = Z;
            for(int m = 0; m < 8; ++m) { // auxiliary number
                int enn = (m + 1)/2;
                for(int ell = m/2; ell >= 0; --ell) { // angular momentum character
                    ++enn; // principal quantum number
                    for(int jj = 2*ell; jj >= 2*ell; jj -= 2) {
                        int const max_occ = std::min(2.f*(jj + 1), ne);
                        core_state[ics].energy = -.5*(Z/enn)*(Z/enn); // hydrogen like energy levels
                        core_state[ics].wave[TRU] = new double[nrt]; // get memory for the (true) radial function
                        core_state[ics].enn = enn;
                        core_state[ics].ell = ell;
                        core_state[ics].emm = emm_Degenerate;
                        core_state[ics].spin = spin_Degenerate;
                        core_state[ics].occupation = max_occ;
                        if (max_occ > 0) {
                            if (echo > 0) printf("%4d%c%4d \n", enn, ellchar[ell], max_occ);
                            enn_core_max[ell] = std::max(enn_core_max[ell], enn); 
                        }
                        ne -= max_occ;
                        ++ics;
                    } // jj
                } // ell
            } // m
        } // core states
        
        nvalencestates = (numax*(numax + 4) + 4)/4;
        valence_state = new valence_level_t[nvalencestates];
        {   int ivs = 0;
            for(int ell = 0; ell <= numax; ++ell) {
                for(int nrn = 0; nrn < nn[ell]; ++nrn) {
                    int const enn = nrn + enn_core_max[ell] + 1;
                    if (echo > 0) printf("# valence%4d%c\n", enn, ellchar[ell]);
                    valence_state[ivs].energy = -.5*(Z/enn)*(Z/enn); // hydrogen like energy levels
                    valence_state[ivs].occupation = 0;
                    valence_state[ivs].enn = enn;
                    valence_state[ivs].ell = ell;
                    valence_state[ivs].emm = emm_Degenerate;
                    valence_state[ivs].spin = spin_Degenerate;
                    valence_state[ivs].wave[TRU] = new double[nrt]; // get memory for the true radial function
                    valence_state[ivs].wave[SMT] = new double[nrs]; // get memory for the smooth radial function
                    ++ivs;
                } // nrn
            } // ell
        } // valence states
        

        int const nlm_rho = (ellmax_density + 1)*(ellmax_density + 1);
        int const nlm_pot = ellmax_potential*(ellmax_potential + 2); // Y00 components are not stored here
        for(int ts = TRU; ts < TRU_AND_SMT; ts += (SMT - TRU)) {
            int const nr = (TRU == ts)? nrt : nrs;
            core_density[ts]   = new double[nr]; // get memory
            potential[ts]      = new double[nr]; // get memory
            full_density[ts]   = new double[nlm_rho*nr]; // get memory
            full_potential[ts] = new double[nlm_pot*nr]; // get memory
        } // true and smooth
        bar_potential = new double[nrs]; // get memory

        int const nSHO = ((1 + numax)*(2 + numax)*(3 + numax))/6;
        matrix_stride = align<2>(nSHO); // 2^<2> doubles = 32 Byte alignment
        hamiltonian = new double[nSHO*matrix_stride]; // get memory
        overlap     = new double[nSHO*matrix_stride]; // get memory

    };
      
    ~LiveAtom() {; // destructor
//      if (rg[SMT] != rg[TRU]) radial_grid::destroy_radial_grid(rg[SMT]); // gives memory corruption
        radial_grid::destroy_radial_grid(rg[TRU]);
        delete[] core_state;
        delete[] valence_state;
    };
      
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
    printf("\n%s: new struct live_atom has size %ld Byte\n\n", __FILE__, sizeof(LiveAtom));
    LiveAtom a(29);
    return 0;
  } // test

  status_t all_tests() {
    auto status = 0;
    status += test();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace single_atom
