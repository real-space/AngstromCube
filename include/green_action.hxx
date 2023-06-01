#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <vector> // std::vector<T>

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

#include "constants.hxx" // ::pi
#include "green_memory.hxx" // get_memory, free_memory

#ifdef HAS_TFQMRGPU

//  #define DEBUG
    #ifdef HAS_NO_CUDA
        #include "tfQMRgpu/include/tfqmrgpu_cudaStubs.hxx" // cuda... (dummies)
        #define devPtr const __restrict__
    #else  // HAS_NO_CUDA
        #include <cuda.h>
    #endif // HAS_NO_CUDA
    #include "tfQMRgpu/include/tfqmrgpu_memWindow.h" // memWindow_t

#else  // HAS_TFQMRGPU

    #include <utility> // std::pair<T>
    typedef std::pair<size_t,size_t> memWindow_t;
    #ifdef HAS_NO_CUDA
        typedef size_t cudaStream_t;
    #endif // HAS_NO_CUDA

#endif // HAS_TFQMRGPU

#include "green_sparse.hxx"    // ::sparse_t<,>
// #include "green_kinetic.hxx"// ::multiply, ::finite_difference_plan_t
#include "green_kinetic.hxx"   // ::kinetic_plan_t
#include "green_potential.hxx" // ::multiply
#include "green_dyadic.hxx"    // ::multiply, ::dyadic_plan_t


#ifdef    debug_printf
  #undef  debug_printf
#endif // debug_printf

#ifdef    DEBUG
  #define debug_printf(...) { std::printf(__VA_ARGS__); std::fflush(stdout); }
#else  // DEBUG
  #define debug_printf(...)
#endif // DEBUG

namespace green_action {

  struct atom_t {
      double pos[3]; // position
      double sigma; // Gaussian spread
      int32_t gid; // global identifier
      int32_t ia; // local atom index
      int16_t shifts[3]; // periodic image shifts
      uint8_t nc; // number of coefficients, uint8_t sufficient up to numax=9 --> 220 coefficients
      int8_t numax; // SHO basis size
  }; // atom_t 48 Byte, only used in CPU parts

  // Suggestion: this could replace AtomPos + AtomLmax in the long run --> ToDo
  struct atom_image_t {
      double pos[3]; // atom position
      float  oneoversqrtsigma; // Gaussian spread, 1/sqrt(sigma)
      int8_t shifts[3]; // periodic image shifts in [-127, 127]   OR   uint16_t phase_index;
      int8_t lmax; // SHO basis size >= -1
  }; // atom_image_t 32 Byte

  struct plan_t {

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

    //   // these members are now replaced by kinetic_plan_t kinetic
    //   green_sparse::sparse_t<int32_t> kinetic_plan[3];
    //   double *grid_spacing = nullptr; // [3]
    //   int16_t kinetic_nFD[4];

      green_kinetic::kinetic_plan_t kinetic[3]; // new struct replacing [kinetic_plan, grid_spacing, kinetic_nFD]

      uint32_t* RowStart = nullptr; // [nRows + 1] Needs to be transfered to the GPU?
      uint32_t* rowindx  = nullptr; // [nnzb] // allows different parallelization strategies
      int16_t (*source_coords)[3+1] = nullptr; // [nCols][3+1] internal coordinates
      int16_t (*target_coords)[3+1] = nullptr; // [nRows][3+1] internal coordinates
      float   (*rowCubePos)[3+1]    = nullptr; // [nRows]      internal coordinates in float, could be int16_t for most applications
      float   (*colCubePos)[3+1]    = nullptr; // [nCols]      internal coordinates in float, could be int16_t for most applications
      int16_t (*target_minus_source)[3+1] = nullptr; // [nnzb][3+1] coordinate differences
      double  (**Veff)[64]          = nullptr; // effective potential, data layout [4][nRows][64], 4 >= Noco^2
      // Veff could be (*Veff[4])[64], however, then we cannot pass Veff to GPU kernels but have to pass Veff[0], Veff[1], ...
      int32_t*  veff_index          = nullptr; // [nnzb] indirection list, values -1 for non-existent indices

      double *grid_spacing_trunc = nullptr; // [3]
      double (*phase)[2][2]      = nullptr; // [3]

      bool noncollinear_spin = false;

      green_dyadic::dyadic_plan_t dyadic_plan;

      plan_t() {
          debug_printf("# default constructor for %s\n", __func__);
          // please see construct_Green_function in green_function.hxx for the construction of the plan_t
      } // constructor

      ~plan_t() {
          debug_printf("# destruct %s\n", __func__);
          free_memory(RowStart);
          free_memory(rowindx);
          free_memory(source_coords);
          free_memory(target_coords);
          free_memory(target_minus_source);
          for (int mag = 0; mag < 4*(nullptr != Veff); ++mag) {
              free_memory(Veff[mag]);
          } // mag
          free_memory(Veff);
          free_memory(veff_index);
          free_memory(colCubePos);
          free_memory(rowCubePos);
          free_memory(grid_spacing_trunc);
     //   free_memory(grid_spacing);
          free_memory(phase);
     //   for (int dd = 0; dd < 3; ++dd) kinetic_plan[dd].~sparse_t();
      } // destructor

  }; // plan_t





  template <typename floating_point_t=float, int R1C2=2, int Noco=1, int n64=64>
  class action_t { // an action as used in tfQMRgpu
  public:
      typedef floating_point_t real_t;
      static int constexpr LM = Noco*n64, // number of rows per block
                           LN = LM;    // number of columns per block
      // action_t::LN is needed to support the rectangular blocks feature in tfQMRgpu
      //
      // This action is an implicit linear operator onto block-sparse structured data.
      // compatible with the core algorithm of the tfqmrgpu-2.0 library.
      // Blocks are sized [LM][LM].
      // Arithmetic according to complex<real_t>
      // with real_t either float or double
      //
      action_t(plan_t *plan)
        : p(plan), apc(nullptr), aac(nullptr)
      {
          assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));
          debug_printf("# construct %s\n", __func__);
          char* buffer{nullptr};
          take_memory(buffer);
          assert(nullptr != plan);
      } // constructor

      ~action_t() {
          debug_printf("# destruct %s\n", __func__);
          free_memory(apc);
//        free_memory(aac); // currently not used
      } // destructor

      void take_memory(char* &buffer) {
          auto const & dp = p->dyadic_plan;
          auto const natomcoeffs = dp.AtomImageStarts ? dp.AtomImageStarts[dp.nAtomImages] : 0;
          auto const n = size_t(natomcoeffs) * p->nCols;
          // ToDo: could be using GPU memory taking it from the buffer
          // ToDo: how complicated would it be to have only one set of coefficients and multiply in-place?
          apc = get_memory<real_t[R1C2][Noco][LM]>(n, p->echo, "apc");
//        aac = get_memory<real_t[R1C2][Noco][LM]>(n, p->echo, "aac"); // currently not used
          // ToDo: alternatively, we could take GPU device memory from the buffer
      } // take_memory

      void transfer(char* const buffer, cudaStream_t const streamId=0) {
          // no transfers needed since we are using managed memory
          // but we could fill matB with unit blocks here
          // assert(p->subset.size() == p->nCols);
          // clear_on_gpu<real_t[2][LM][LM]>(matB, p->nCols);
          // assert(LM == LM);
          // for (int icol = 0; icol < p->nCols; ++icol) {
          //     for (int i = 0; i < LM; ++i) {
          //        matB[icol][0][i][i] = real_t(1);
          //     } // i
          // } // icol
      } // transfer

      bool has_preconditioner() const { return false; }


      double multiply( // returns the number of flops performed
            real_t         (*const __restrict y)[R1C2][LM][LM] // result, y[nnzb][2][LM][LM]
          , real_t   const (*const __restrict x)[R1C2][LM][LM] // input,  x[nnzb][2][LM][LM]
          , uint16_t const (*const __restrict colIndex) // column indices [nnzb], warning: must be in device memory or managed memory
          , uint32_t const nnzb // number of nonzero blocks, typically colIndex.size()
          , uint32_t const nCols=1 // should match with p->nCols, number of block columns, assert(colIndex[:] < nCols)
          , unsigned const l2nX=0  // number of levels needed for binary reductions over nnzb
          , cudaStream_t const streamId=0 // CUDA stream to run on
          , bool const precondition=false
      )
        // GPU implementation of green_potential, green_kinetic and green_dyadic
      {
          assert(p);
          if (2 == Noco) assert(p->noncollinear_spin && "Also the plan needs to be created with Noco=2");
          double nops{0};

          if (p->echo > 3) std::printf("\n");
          if (p->echo > 2) std::printf("# green_action::multiply\n");

          // start with the local potential, assign y to initial values
          nops += green_potential::multiply<real_t,R1C2,Noco>(y, x, p->Veff, p->veff_index,
                      p->target_minus_source, p->grid_spacing_trunc, nnzb, p->E_param,
                      p->V_confinement, pow2(p->r_confinement), p->echo);

          // add the kinetic energy expressions
          for (int dd = 0; dd < 3; ++dd) { // loop must run serial
              nops += p->kinetic[dd].multiply<real_t,R1C2,Noco>(y, x, p->phase[dd], p->echo);
          } // dd derivative direction

          // add the non-local potential using the dyadic action of project + add
          nops += green_dyadic::multiply<real_t,R1C2,Noco>(y, apc, x, p->dyadic_plan,
                      p->rowindx, colIndex, p->rowCubePos, nnzb, p->echo);

          if (p->echo > 4) std::printf("# green_action::multiply %g Gflop\n", nops*1e-9);

          return nops;
      } // multiply

      plan_t * get_plan() { return p; }

    private: // members

      plan_t *p; // the plan is independent of real_t and R1C2

      // temporary device memory needed for dyadic operations
      real_t (*apc)[R1C2][Noco][LM]; // atom projection coefficients apc[n_all_projection_coefficients*nCols][R1C2][Noco][Noco*64]
      real_t (*aac)[R1C2][Noco][LM]; // atom   addition coefficients aac[n_all_projection_coefficients*nCols][R1C2][Noco][Noco*64]
      // (we could live with a single copy if the application of the atom-centered matrices is in-place)

  }; // class action_t


#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_construction_and_destruction(int const echo=0) {
      {   plan_t plan;
          if (echo > 4) std::printf("# %s for action_t\n", __func__);
          { action_t<float ,1,1> action(&plan); }
          { action_t<float ,2,1> action(&plan); }
          { action_t<float ,2,2> action(&plan); }
          { action_t<double,1,1> action(&plan); }
          { action_t<double,2,1> action(&plan); }
          { action_t<double,2,2> action(&plan); }
          if (echo > 5) std::printf("# Hint: to test action_t::multiply, please envoke --test green_function\n");
      } // destruct plan
      if (echo > 6) {
          std::printf("# %s sizeof(atom_t) = %ld Byte\n", __func__, sizeof(atom_t));
          std::printf("# %s sizeof(atom_image_t) = %ld Byte\n", __func__, sizeof(atom_image_t));
          std::printf("# %s sizeof(plan_t) = %ld Byte\n", __func__, sizeof(plan_t));
      } // echo
      return 0;
  } // test_construction_and_destruction

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_construction_and_destruction(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_action
