// This file is part of AngstromCube under MIT License

#include <cstdio>     // std::printf, ::snprintf, FILE, std::fprintf
#include <cstdint>    // int64_t, int32_t, uint32_t, int16_t, uint16_t, int8_t, uint8_t
#include <cassert>    // assert
#include <algorithm>  // std::max, ::min
#include <vector>     // std::vector<T>
#include <complex>    // std::complex

#include "green_action.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#ifndef NO_UNIT_TESTS
  #include "simple_timer.hxx" // SimpleTimer
  #include "simple_stats.hxx" // ::Stats<>
  #include "display_units.h" // eV, _eV, Ang, _Ang
  #include "green_memory.hxx" // get_memory, free_memory, real_t_name
  #include "progress_report.hxx" // ProgressReport
  #include "green_parallel.hxx" // ::RequestList_t
  #include "green_action.hxx" // ::plan_t, ::action_t, ::atom_t
  #include "green_potential.hxx" // ::exchange
  #include "green_dyadic.hxx" // ::dyadic_plan_t
  #include "sho_tools.hxx" // ::nSHO
  #include "control.hxx" // ::get
  #include "load_balancer.hxx" // ::get
  #include "display_units.h" // GByte, _GByte
  #include "green_input.hxx" // ::load_Hamiltonian
  #include "green_function.hxx" // ::construct_Green_function
  #include "mpi_parallel.hxx" // ::init, ::finalize, ::rank

  #ifdef    HAS_TFQMRGPU

    #ifdef    HAS_NO_CUDA
        #include "tfqmrgpu_cudaStubs.hxx" // cuda... (dummies)
        #define devPtr const __restrict__
    #else  // HAS_NO_CUDA
        #include <cuda.h>
        #define devPtr const __restrict__
    #endif // HAS_NO_CUDA
    #include "tfqmrgpu_core.hxx" // tfqmrgpu::solve<action_t>

  #endif // HAS_TFQMRGPU

    #ifdef    DEBUG
        #undef DEBUG
    #endif // DEBUG

#endif // NO_UNIT_TESTS


#ifdef    HAS_BITMAP_EXPORT
  #include "bitmap.hxx" // ::write_bmp_file
#endif // HAS_BITMAP_EXPORT

// #define   GREEN_FUNCTION_SVG_EXPORT

 /*
  *  Future plan:
  *   Support also density matrix purification scheme (McWeeney filter: x^2(3-2x)
  *   or with a non-trivial overlap operator S (3xSx - 2xSxSx from Phys. Rev. B 50, 17622)
  *   maybe better with norm-conserving PAW formalism --> S == 1
  *   in particular suitable for the real Hamiltonian action (which is not supported by tfQMRgpu)
  */

namespace green_action {

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <typename real_t, int R1C2=2, int Noco=1>
  void test_action(action_plan_t & p, int const iterations=1, int const echo=9) {
      if (echo > 1) std::printf("# %s<%s,R1C2=%d,Noco=%d>\n", __func__, real_t_name<real_t>(), R1C2, Noco);
      green_action::action_t<real_t,R1C2,Noco,64> action(&p); // constructor

      uint32_t const nnzbX = p.colindx.size();
      if (echo > 3) std::printf("# memory of a Green function is %.6f %s\n", nnzbX*R1C2*pow2(64.*Noco)*sizeof(real_t)*GByte, _GByte);

      if (0 == iterations) { 
          if (echo > 2) std::printf("# requested to run no iterations --> only check the action_t constructor\n");
          return;
      } // 0 iterations

#ifdef    HAS_TFQMRGPU
      if (iterations > 0) {
          if (nnzbX < 1) {
              if (echo > 2) std::printf("# cannot call tfqmrgpu library if X has no elements!\n");
              return;
          }
          p.echo = echo - 5;
          if (echo > 0) std::printf("\n# call tfqmrgpu::mem_count\n");

          // beware, the changes only the local potential. In a non-benchmark situation use ::update_energy_parameter
          p.E_param = std::complex<double>(control::get("green_function.energy.parameter.real", 0.0),
                                           control::get("green_function.energy.parameter.imag", 0.0));

          // try to instanciate tfqmrgpu::solve with this action_t<real_t,R1C2,Noco,64>
          tfqmrgpu::solve(action); // compute GPU memory requirements

          {
              simple_stats::Stats<> mem; mem.add(p.gpu_mem); mpi_parallel::allreduce(mem); // uses MPI_COMM_WORLD
              if (echo > 5) std::printf("# tfqmrgpu needs [%.1f, %.1f +/- %.1f, %.1f] %s GPU memory, %.3f %s total\n",
                mem.min()*GByte, mem.mean()*GByte, mem.dev()*GByte, mem.max()*GByte, _GByte, mem.sum()*GByte, _GByte);
          }
          auto memory_buffer = get_memory<char>(p.gpu_mem, echo, "tfQMRgpu-memoryBuffer");
          int const maxiter = control::get("tfqmrgpu.max.iterations", 99.);
          if (echo > 0) std::printf("\n# call tfqmrgpu::solve\n\n");
          double time_needed{1};
          { // scope: benchmark the solver
              SimpleTimer timer(__FILE__, __LINE__, __func__, echo);

              tfqmrgpu::solve(action, memory_buffer, 1e-9, maxiter, 0, true);

              time_needed = timer.stop();
          } // timer
          if (echo > 0) std::printf("\n# after tfqmrgpu::solve residuum reached= %.1e iterations needed= %d\n",
                                                             p.residuum_reached,    p.iterations_needed);
          if (echo > 6) std::printf("# after tfqmrgpu::solve flop count is %.6f %s\n", p.flops_performed*1e-9, "Gflop");
          if (echo > 6) std::printf("# estimated performance is %.6f %s\n", p.flops_performed*1e-9/time_needed, "Gflop/s");
          free_memory(memory_buffer);
          return;
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

              nflops += action.multiply(y, x, colIndex, nnzbX, p.nCols);
              cudaDeviceSynchronize();

              timings.add(timeit.stop());
              std::swap(x, y);
              p.echo = 0; // mute after the 1st iteration
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
  } // test_action


  status_t test_green_action(int const echo=0) {
      bool const already_initialized = mpi_parallel::init();

      uint32_t ng[3] = {0, 0, 0}; // grid sizes
      int8_t   bc[3] = {0, 0, 0}; // boundary conditions
      double   hg[3] = {1, 1, 1}; // grid spacings
      std::vector<double> Veff(0); // local potential
      int natoms{0}; // number of atoms
      std::vector<double> xyzZinso(0); // atom info
      std::vector<std::vector<double>> AtomMatrices(0); // non-local potential

      auto const *const filename = control::get("hamiltonian.file", "Hmt.xml");
      auto stat = green_input::load_Hamiltonian(ng, bc, hg, Veff, natoms, xyzZinso, AtomMatrices, filename, echo - 5);
      if (stat) {
          warn("failed to load_Hamiltonian with status=%d", int(stat));
          if (!already_initialized) mpi_parallel::finalize();
          return stat;
      } // stat

      int const r1c2 = control::get("green_function.benchmark.complex", 1.) + 1;
      int const noco = control::get("green_function.benchmark.noco", 1.);

//    for (int ia = 0; ia < natoms; ++ia) { xyzZinso[ia*8 + 3] = 6; } // set all atoms to carbon

      action_plan_t p;
      stat += green_function::construct_Green_function(p, ng, bc, hg, xyzZinso, echo, noco);

      assert(1 == r1c2 || 2 == r1c2);
      assert(1 == noco || r1c2 == noco);
      int const fp = control::get("green_function.benchmark.floating.point.bits", 32.);
      int const iterations = control::get("green_function.benchmark.iterations", 1.);
      int const action_key = 1000*((32 == fp) ? 32 : 64) + 10*r1c2 + noco;
      int const action = control::get("green_function.benchmark.action", action_key*1.);
                      // -1: no iterations, 0:run memory initialization only, >0: iterate
      // try one of the 6 combinations (strangely, we cannot run any two of these calls after each other, ToDo: find out what's wrong here)
      switch (action) {
          case 32022: test_action<float ,2,2>(p, iterations, echo); break; // complex non-collinear
          case 64022: test_action<double,2,2>(p, iterations, echo); break; // complex non-collinear

          case 32021: test_action<float ,2,1>(p, iterations, echo); break; // complex
          case 64021: test_action<double,2,1>(p, iterations, echo); break; // complex
#ifdef    HAS_TFQMRGPU
          case 32011:                                                       // real
          case 64011: error("tfQMRgpu needs R1C2 == 2 but found green_function.benchmark.action=%d", action); break;
#else  // HAS_TFQMRGPU
          case 32011: test_action<float ,1,1>(p, iterations, echo); break; // real
          case 64011: test_action<double,1,1>(p, iterations, echo); break; // real
#endif // HAS_TFQMRGPU
          case 0: if (echo > 1) std::printf("# green_function.benchmark.action=0 --> test_action is not called!\n"); break;
          default: warn("green_function.benchmark.action must be in {32011, 32021, 32022, 64011, 64021, 64022} but found %d", action);
                   ++stat;
      } // switch action

      if (!already_initialized) mpi_parallel::finalize();
      return stat;
  } // test_green_action

  inline status_t test_construction_and_destruction(int const echo=0) {
      {   action_plan_t plan;
          if (echo > 4) std::printf("# %s for action_t\n", __func__);
#ifndef   HAS_TFQMRGPU
          { action_t<float ,1,1> action(&plan); }
          { action_t<double,1,1> action(&plan); }
#endif // HAS_TFQMRGPU
          { action_t<float ,2,1> action(&plan); }
          { action_t<float ,2,2> action(&plan); }
          { action_t<double,2,1> action(&plan); }
          { action_t<double,2,2> action(&plan); }
          if (echo > 5) std::printf("# Hint: to test action_t::multiply, please envoke --test green_function\n");
      } // destruct plan
      if (echo > 6) std::printf("# %s sizeof(plan_t) = %ld Byte\n", __func__, sizeof(action_plan_t));
      return 0;
  } // test_construction_and_destruction

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_construction_and_destruction(echo);
      stat += test_green_action(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_action
