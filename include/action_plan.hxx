#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <complex> // std::complex<real_t>

// ToDo: get rid of memWindow_t
#ifdef    HAS_TFQMRGPU

    #include "tfqmrgpu_memWindow.h" // memWindow_t

#else  // HAS_TFQMRGPU

    #include <utility> // std::pair<T>
    typedef std::pair<size_t,size_t> memWindow_t;

#endif // HAS_TFQMRGPU

#ifdef    DEBUG
    #undef DEBUG
#endif // DEBUG

#include "status.hxx" // status_t

#include "kinetic_plan.hxx" // kinetic_plan_t
#include "dyadic_plan.hxx" // dyadic_plan_t
#include "green_parallel.hxx" // ::RequestList_t

class action_plan_t {
public: // TODo: check which members could be private

    // members needed for the usage with tfQMRgpu

    // for the inner products and axpy/xpay
    std::vector<uint16_t> colindx; // [nnzbX], must be a std::vector since nnzb is derived from colindx.size()
    memWindow_t colindxwin; // column indices in GPU memory

    // for the matrix-submatrix addition/subtraction Y -= B:
    std::vector<uint32_t> subset; // [nnzbB], list of inzbX-indices at which B is also non-zero
    memWindow_t subsetwin; // subset indices in GPU memory
    memWindow_t matBwin; // data of the right hand side operator B

    // memory positions
    memWindow_t matXwin; // solution vector in GPU memory
    memWindow_t vec3win; // random number vector in GPU memory

    uint32_t nRows = 0; // number of block rows in the Green function
    uint32_t nCols = 0; // number of block columns, max 65,536 columns
    std::vector<uint32_t> rowstart; // [nRows + 1] does not need to be transferred to the GPU

    // for memory management:
    size_t gpu_mem = 0; // device memory requirement in Byte (can this be tracked?)

    // stats:
    float residuum_reached    = 3e38;
    float flops_performed     = 0.f;
    float flops_performed_all = 0.f;
    int   iterations_needed   = -99;

    // =====================================================================================
    // additional members to define the ultra-sparse PAW Hamiltonian action

    int echo = 9;

    std::vector<int64_t> global_target_indices; // [nRows]
    std::vector<int64_t> global_source_indices; // [nCols]
    double r_truncation   = 9e18; // radius beyond which the Green function is truncated, in Bohr
    float r_confinement   = 9e18; // radius beyond which the confinement potential is added, in Bohr
    float V_confinement   = 1; // potential prefactor
    std::complex<double> E_param; // energy parameter

    kinetic_plan_t kinetic[3]; // plan to execute the kinetic energy operator

    uint32_t* RowStart = nullptr; // [nRows + 1] Needs to be transfered to the GPU?
    uint32_t* rowindx  = nullptr; // [nnzbX] // allows different parallelization strategies
    // int16_t (*source_coords)[3+1] = nullptr; // [nCols][3+1] internal coordinates
    // int16_t (*target_coords)[3+1] = nullptr; // [nRows][3+1] internal coordinates
    float   (*rowCubePos)[3+1]    = nullptr; // [nRows][3+1] internal coordinates in float, could be int16_t for most applications
    float   (*colCubePos)[3+1]    = nullptr; // [nCols][3+1] internal coordinates in float, could be int16_t for most applications
    int16_t (*target_minus_source)[3+1] = nullptr; // [nnzbX][3+1] coordinate differences                                               TODO: remove target_minus_source
    double  (**Veff)[64]          = nullptr; // effective potential, data layout [4][nRows][64], 4 >= Noco^2
    // Veff could be (*Veff[4])[64], however, then we cannot pass Veff to GPU kernels but have to pass Veff[0], Veff[1], ...
    int32_t*  veff_index          = nullptr; // [nnzbX] indirection list, values -1 for non-existent indices

    double *grid_spacing_trunc = nullptr; // [3]
    double (*phase)[2][2]      = nullptr; // [3] // phase factors for the 3 directions across the boundaries, used in kinetic and dyadic phases are derived from it

    bool noncollinear_spin = false;

    dyadic_plan_t dyadic_plan; // plan to execute the dyadic potential operator

    green_parallel::RequestList_t potential_requests; // request list to exchange potential blocks
    green_parallel::RequestList_t matrices_requests;  // request list to exchange atomic matrices


public:

    action_plan_t(int const echo=0); // default constructor

    action_plan_t( // constructor
        uint32_t const ng[3] // numbers of grid points of the unit cell in with the potential is defined
      , int8_t const bc[3] // boundary conditions in {Isolated, Periodic, Vacuum, Repeat}
      , double const hg[3] // grid spacings
      , std::vector<double> const & xyzZinso // [natoms*8]
      , int const echo // =0 // log-level
      , int const Noco // =2
    ); // declaration only

    ~action_plan_t(); // destructor

}; // action_plan_t



namespace action_plan {

  struct atom_t {
      double pos[3]; // position
      double sigma; // Gaussian spread
      int32_t gid; // global identifier
      int32_t ia; // local original atom index
      int32_t iaa; // local atom copy index
      int16_t shifts[3]; // periodic image shifts
      uint8_t nc; // number of coefficients, uint8_t sufficient up to numax=9 --> 220 coefficients
      int8_t numax; // SHO basis size
      int16_t copies[3]; // periodic copy shifts
  }; // atom_t 56 Byte, only used in CPU parts

  status_t all_tests(int const echo=0); // declaration only

} // namespace action_plan
