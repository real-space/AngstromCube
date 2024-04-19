#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <vector> // std::vector<T>

#ifdef HAS_TFQMRGPU

//  #define DEBUG
    #ifdef HAS_NO_CUDA
        #include "tfqmrgpu_cudaStubs.hxx" // cuda... (dummies)
        #define devPtr const __restrict__
    #else  // HAS_NO_CUDA
        #include <cuda.h>
        #define devPtr const __restrict__
    #endif // HAS_NO_CUDA
    #include "tfqmrgpu_memWindow.h" // memWindow_t
    #include "tfqmrgpu_core.hxx" // tfqmrgpu::solve<action_t>

#else  // HAS_TFQMRGPU

    #include <utility> // std::pair<T>
    typedef std::pair<size_t,size_t> memWindow_t;
    #ifdef HAS_NO_CUDA
        typedef size_t cudaStream_t;
    #endif // HAS_NO_CUDA

#endif // HAS_TFQMRGPU

#include "green_kinetic.hxx"   // ::multiply
#include "green_potential.hxx" // ::multiply
#include "green_dyadic.hxx"    // ::multiply
#include "action_plan.hxx"     // action_plan_t
#include "kinetic_plan.hxx"    // kinetic_plan_t
#include "dyadic_plan.hxx"     // dyadic_plan_t
#include "green_sparse.hxx"    // ::sparse_t<,>
#include "simple_timer.hxx"    // SimpleTimer
#include "progress_report.hxx" // ProgressReport
#include "display_units.h"     // GByte, _GByte
#include "constants.hxx"       // ::pi
#include "green_memory.hxx"    // get_memory, free_memory
#include "status.hxx"          // status_t, STATUS_TEST_NOT_INCLUDED
#include "mpi_parallel.hxx"    // ::allreduce, ::rank
#include "recorded_warnings.hxx" // warn

#ifdef    DEBUG
  #define green_debug_printf(...) { std::printf(__VA_ARGS__); std::fflush(stdout); }
#else  // DEBUG
  #define green_debug_printf(...)
#endif // DEBUG

namespace green_action {

    template <typename floating_point_t=float, unsigned R1C2=2, unsigned Noco=1, unsigned n64=64>
    class action_t { // an action as used in tfQMRgpu
    public:
        typedef floating_point_t real_t;
        static unsigned constexpr LM = Noco*n64, // number of rows per block
                                  LN = LM;    // number of columns per block
        // action_t::LN is needed to support the rectangular blocks feature in tfQMRgpu
        //
        // This action is an implicit linear operator onto block-sparse structured data.
        // compatible with the core algorithm of the tfqmrgpu-2.0 library.
        // Blocks are sized [LM][LM].
        // Arithmetic according to complex<real_t>
        // with real_t either float or double
        //

        action_t(action_plan_t *plan, int const echo=0)
          : p_(plan), apc_(nullptr) // , aac_(nullptr)
        {
            if (echo > 1) std::printf("# construct %s<%s,R1C2=%d,Noco=%d>\n", __func__, real_t_name<real_t>(), R1C2, Noco);
            assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));
            char* buffer{nullptr};
            take_memory(buffer);
            assert(nullptr != plan);
            auto & p = *p_;
            auto const nnzbX = p.colindx.size();
#ifdef    HAS_TFQMRGPU
            if (nnzbX > 0) {
                if (echo > 0) std::printf("\n# call tfqmrgpu::mem_count\n");
                // try to instanciate tfqmrgpu::solve with this action_t<real_t,R1C2,Noco,64>
                tfqmrgpu::solve(*this); // compute GPU memory requirements
                auto const me = mpi_parallel::rank();                                  // uses MPI_COMM_WORLD
                {
                    simple_stats::Stats<> m; m.add(p.gpu_mem); mpi_parallel::allreduce(m); // uses MPI_COMM_WORLD
                    if (echo > 5) std::printf("# tfqmrgpu needs [%.1f, %.1f +/- %.1f, %.1f] %s GPU memory, %.3f %s total\n",
                                m.min()*GByte, m.mean()*GByte, m.dev()*GByte, m.max()*GByte, _GByte, m.sum()*GByte, _GByte);
                    if (echo > 9) std::printf("# rank#%i try to allocate %.9f %s green_memory\n", me, p.gpu_mem*GByte, _GByte);
                }
                memory_buffer_ = get_memory<char>(p.gpu_mem, echo, "tfQMRgpu-memoryBuffer");
// #ifdef    DEBUGGPU
                if (echo > 9) std::printf("# rank#%i allocated %.9f %s memory_buffer_ at %p\n", me, p.gpu_mem*GByte, _GByte, (void*)memory_buffer_);
// #endif // DEBUGGPU
            } else {
                if (echo > 2) std::printf("# cannot call tfqmrgpu library if X has no elements!\n");
            }
#else  // HAS_TFQMRGPU
            if (echo > 3) std::printf("# memory of a Green function is %.6f %s\n", nnzbX*R1C2*pow2(64.*Noco)*sizeof(real_t)*GByte, _GByte);
#endif // HAS_TFQMRGPU
        } // constructor

        ~action_t() {
            green_debug_printf("# destruct %s\n", __func__);
            free_memory(apc_);
  //        free_memory(aac_); // currently not used
            free_memory(memory_buffer_);
        } // destructor

      void take_memory(char* &buffer) {
          auto const & dp = p_->dyadic_plan;
          auto const natomcoeffs = dp.AtomImageStarts ? dp.AtomImageStarts[dp.nAtomImages] : 0;
          auto const n = size_t(natomcoeffs) * p_->nCols;
          apc_ = get_memory<real_t[R1C2][Noco][LM]>(n, p_->echo, "apc");
//        aac_ = get_memory<real_t[R1C2][Noco][LM]>(n, p_->echo, "aac"); // currently not used
      } // take_memory

      void transfer(char* const buffer, cudaStream_t const streamId=0) {
          // no transfers needed since we are using managed memory
      } // transfer

      bool has_preconditioner() const { return false; }


      double multiply( // returns the number of flops performed
            real_t         (*const __restrict y)[R1C2][LM][LM] // result, y[nnzb][2][LM][LM]
          , real_t   const (*const __restrict x)[R1C2][LM][LM] // input,  x[nnzb][2][LM][LM]
          , uint16_t const (*const __restrict colIndex) // column indices [nnzb], warning: must be in device memory or managed memory
          , uint32_t const nnzb // number of nonzero blocks, typically colIndex.size()
          , uint32_t const nCols=1 // should match with p.nCols, number of block columns, assert(colIndex[:] < nCols)
          , unsigned const l2nX=0  // number of levels needed for binary reductions over nnzb
          , cudaStream_t const streamId=0 // CUDA stream to run on
          , bool const precondition=false
      )
        // GPU implementation of green_potential, green_kinetic and green_dyadic
      {
          assert(p_); auto const & p = *p_;
          if (2 == Noco) assert(p.noncollinear_spin && "Also the plan needs to be created with Noco=2");
          double nops{0};

          if (p.echo > 3) std::printf("\n");
          if (p.echo > 2) std::printf("# green_action::multiply\n");

          // start with the local potential, assign y to initial values
          nops += green_potential::multiply<real_t,R1C2,Noco>(y, x, p.Veff, p.veff_index,
                      p.target_minus_source, p.grid_spacing_trunc, nnzb, p.E_param,
                      p.V_confinement, pow2(p.r_confinement), p.echo);

          // add the kinetic energy expressions
          for (int dd = 0; dd < 3; ++dd) { // loop must run serial
              nops += green_kinetic::multiply<real_t,R1C2,Noco>(y, x, p.kinetic[dd], p.phase[dd], p.echo);
          } // dd derivative direction

          // add the non-local potential using the dyadic action of project + add
          nops += green_dyadic::multiply<real_t,R1C2,Noco>(y, apc_, x, p.dyadic_plan,
                      p.rowindx, colIndex, p.rowCubePos, nnzb, p.echo);

          if (p.echo > 4) std::printf("# green_action::multiply %g Gflop\n", nops*1e-9);

          return nops;
      } // multiply

      action_plan_t * get_plan() { return p_; }


    status_t solve(
          double rho[] // result: density[nblocks][4*4*4]
        , uint32_t const nblocks // should match plan.nCols
        , int const iterations=1
        , int const echo=9
    ) {
        if (echo > 1) std::printf("# action_t<%s,R1C2=%d,Noco=%d>::%s\n", real_t_name<real_t>(), R1C2, Noco, __func__);

        auto const me = mpi_parallel::rank(); // usues MPI_COMM_WORLD
// #ifdef    DEBUGGPU
        if (echo > 7) std::printf("# rank#%i action_t at %p usues memory_buffer_ at %p\n", me, (void*)this, (void*)memory_buffer_);
// #endif // DEBUGGPU

        assert(p_); auto const & p = *p_;
        uint32_t const nnzbX = p.colindx.size();
        if (echo > 3) std::printf("# memory of a Green function is %.6f %s\n", nnzbX*R1C2*pow2(64.*Noco)*sizeof(real_t)*GByte, _GByte);

        if (0 == iterations) { 
            if (echo > 2) std::printf("# requested to run no iterations --> only check the action_t constructor\n");
            return 0;
        } // 0 iterations

#ifdef    HAS_TFQMRGPU
        if (iterations > 0) {
            if (nnzbX < 1) {
                if (echo > 2) std::printf("# cannot call tfqmrgpu library if X has no elements!\n");
                return 0;
            }
            int const maxiter = iterations;
            if (echo > 0) std::printf("\n# call tfqmrgpu::solve\n\n");
            assert(nullptr != memory_buffer_);
            double time_needed{1};
            { // scope: benchmark the solver
                SimpleTimer timer(__FILE__, __LINE__, __func__, echo);

                tfqmrgpu::solve(*this, memory_buffer_, 1e-9, maxiter, 0, true);

                time_needed = timer.stop();
            } // timer
            if (echo > 0) std::printf("\n# after tfqmrgpu::solve residuum reached= %.1e iterations needed= %d\n",
                                                               p.residuum_reached,    p.iterations_needed);
            if (echo > 6) std::printf("# after tfqmrgpu::solve flop count is %.6f %s\n", p.flops_performed*1e-9, "Gflop");
            if (echo > 6) std::printf("# estimated performance is %.6f %s\n", p.flops_performed*1e-9/time_needed, "Gflop/s");
            // export solution

            auto const Green = (real_t const (*)[2][Noco*64][Noco*64])memory_buffer_;
            if (echo > 2) std::printf("# copy %d diagonal blocks of the Green function\n", p.nCols);
            if (nblocks != p.nCols) warn("Green function solution provides %d 4x4x4 blocks, but requested %d", p.nCols, nblocks);
            for (uint32_t iCol{0}; iCol < p.nCols; ++iCol) {
                auto const inz_diagonal = p.subset.at(iCol); // works since we have non-zeros in B only on the diagonal
                for (unsigned i64{0}; i64 < 64; ++i64) {
                    int constexpr imag = 1; // imaginary part
                    rho[iCol*64u + i64] = Green[inz_diagonal][imag][i64][i64];
                } // i64
            } // iCol

            return 0;
        } // iterations > 0
#endif // HAS_TFQMRGPU

        int const niterations = std::abs(iterations);
        int constexpr LM = Noco*64;
        auto x = get_memory<real_t[R1C2][LM][LM]>(nnzbX, echo, "x");
        auto y = get_memory<real_t[R1C2][LM][LM]>(nnzbX, echo, "y");
        set(x[0][0][0], nnzbX*size_t(R1C2*LM*LM), real_t(0)); // init x

        auto colIndex = get_memory<uint16_t>(nnzbX, echo, "colIndex");
        set(colIndex, nnzbX, p.colindx.data()); // copy into GPU memory

        { // scope: benchmark the action
            SimpleTimer timer(__FILE__, __LINE__, __func__, echo);
            simple_stats::Stats<> timings;
            double nflops{0};
            ProgressReport progress(__FILE__, __LINE__, 2.5, echo); // update the line every 2.5 seconds
            for (int iteration = 0; iteration < niterations; ++iteration) {
                SimpleTimer timeit(__FILE__, __LINE__, __func__, echo*0);

                nflops += multiply(y, x, colIndex, nnzbX, p.nCols);
                cudaDeviceSynchronize();

                timings.add(timeit.stop());
                std::swap(x, y);
                progress.report(iteration, niterations);
            } // iteration
            if (echo > 1) std::printf("#\n# running action.multiply needed [%g, %g +/- %g, %g] seconds per iteration\n",
                                            timings.min(), timings.mean(), timings.dev(), timings.max());
            char const fF = (sizeof(real_t) == 8) ? 'F' : 'f';
            if (echo > 1) std::printf("# %d calls of action.multiply performed %.3e %clop in %.3e seconds, i.e. %g G%clop/s\n",
                                            niterations, nflops, fF, timings.sum(), nflops/timings.sum()*1e-9, fF);
            if (echo > 1) std::printf("# fastest call of action.multiply performed %.3e %clop in %.3e seconds, i.e. %g G%clop/s\n",
                                            nflops/niterations, fF, timings.min(), nflops/(niterations*timings.min())*1e-9, fF);
        } // scope

        free_memory(colIndex);
        free_memory(y);
        free_memory(x);
        return 0;
    } // solve




    private: // members

      action_plan_t *p_ = nullptr; // the action_plan is independent of real_t and R1C2 and stores Noco as member (to check the matching)

      // temporary device memory needed for dyadic operations
      real_t (*apc_)[R1C2][Noco][LM] = nullptr; // atom projection coefficients apc[n_all_projection_coefficients*nCols][R1C2][Noco][Noco*64]
//    real_t (*aac_)[R1C2][Noco][LM] = nullptr; // atom   addition coefficients aac[n_all_projection_coefficients*nCols][R1C2][Noco][Noco*64]
      // (we can live with a single copy as the application of the atom-centered matrices is in-place)

      char* memory_buffer_ = nullptr;

  }; // class action_t


    status_t all_tests(int const echo=0); // declaration only

} // namespace green_action

#undef green_debug_printf
