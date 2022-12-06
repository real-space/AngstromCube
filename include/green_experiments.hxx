#pragma once

#include <cstdio> // std::printf
#include <vector> // std::vector<T>

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "simple_timer.hxx" // SimpleTimer
#include "unit_system.hxx" // eV, _eV, Kelvin, _Kelvin
#include "green_input.hxx" // ::load_Hamitonian
#include "data_view.hxx" // view2D<T>
#include "brillouin_zone.hxx" // ::get_kpoint_mesh, ::get_kpoint_path
#include "simple_stats.hxx" // ::Stats

#ifdef HAS_TFQMRGPU

    #ifdef HAS_NO_CUDA
        #include "tfQMRgpu/include/tfqmrgpu_cudaStubs.hxx" // cuda... (dummies)
        #define devPtr const __restrict__
    #else  // HAS_NO_CUDA
        #include <cuda.h>
    #endif // HAS_NO_CUDA
    #include "tfQMRgpu/include/tfqmrgpu_memWindow.h" // memWindow_t
    #include "tfQMRgpu/include/tfqmrgpu_linalg.hxx" // ...
    #include "tfQMRgpu/include/tfqmrgpu_core.hxx" // tfqmrgpu::solve<action_t>

#endif // HAS_TFQMRGPU

#include "green_action.hxx" // ::plan_t, ::action_t
#include "green_function.hxx" // ::update_phases, ::construct_Green_function
#include "control.hxx" // ::get

namespace green_experiments {

  double const GByte = 1e-9; char const *const _GByte = "GByte";

  template <typename real_t=double, int Noco=1>
  inline status_t bandstructure(
        green_action::plan_t & p
      , std::vector<std::vector<double>> const & AtomMatrices
      , double const box[3]
      , int const echo=0
  ) {
      int constexpr R1C2 = 2;

      if (echo > 1) std::printf("# %s<%s,R1C2=%d,Noco=%d>\n", __func__, real_t_name<real_t>(), R1C2, Noco);

      green_action::action_t<real_t,R1C2,Noco,64> action(&p); // constructor


      p.gpu_mem = 0;
#ifdef  HAS_TFQMRGPU
      p.echo = echo - 5;
      if (echo > 0) std::printf("\n# call tfqmrgpu::mem_count\n");
      tfqmrgpu::solve(action); // try to instanciate tfqmrgpu::solve with this action_t<real_t,R1C2,Noco,64>
      if (echo > 5) std::printf("# tfqmrgpu::solve requires %.6f GByte GPU memory\n", p.gpu_mem*1e-9);
      double constexpr prefactor = 1./constants::pi;
      int constexpr ImaginaryPart = 1;
      assert(1 == Noco);
      auto rho = get_memory<double[Noco][Noco][64]>(p.nCols, echo, "rho");
      set(rho[0][0][0], p.nCols*Noco*Noco*64, 0.0);
#endif
      auto memory_buffer = get_memory<char>(p.gpu_mem, echo, "tfQMRgpu-memoryBuffer");

      int const maxiter = control::get("tfqmrgpu.max.iterations", 99.);
      auto const E0     = control::get("green_experiments.bandstructure.energy.offset", 0.0);
      auto const dE     = control::get("green_experiments.bandstructure.energy.spacing", 0.01);
      int const nE      = control::get("green_experiments.bandstructure.energy.points", 32.);
      auto const E_imag = control::get("green_experiments.bandstructure.energy.imag", 1e-3); // room temperature
      if (E_imag <= 0) warn("imaginary part of the energy parameter is %.1e Ha", E_imag);

      if (echo > 3) std::printf("# +tfqmrgpu.max.iterations=%d\n", maxiter);
      if (echo > 1) std::printf("# %d E-points in [%g, %g] %s\n", nE, E0*eV, (E0 + (nE - 1)*dE)*eV, _eV);

      view2D<double> k_path;
      auto const nkpoints = brillouin_zone::get_kpoint_path(k_path, echo);
      if (echo > 1) std::printf("# %s %d k-points, %d E-points, Temperature %g %s = %g %s\n",
                            __func__, nkpoints, nE, E_imag*Kelvin, _Kelvin, E_imag*eV, _eV);
      auto constexpr twopi = 2*constants::pi;
      double const reci[] = {twopi/box[0], twopi/box[1], twopi/box[2]};

      std::vector<double> bandstructure(nkpoints, -9e9);

      for (int ik = 0; ik < nkpoints; ++ik) {
          double const *const k_point = k_path[ik];

          if (echo > 0) std::printf("\n## k-point %g %g %g\n", k_point[0], k_point[1], k_point[2]);
          green_function::update_phases(p, k_point, Noco, echo);

          double E_resonance{-9};
#ifdef  HAS_TFQMRGPU
          double max_resonance{-9e9};
#endif // HAS_TFQMRGPU
          for (int iE = 0; iE < nE; ++iE) {
              double const E_real = iE*dE + E0;
              std::complex<double> E_param(E_real, E_imag);

              green_function::update_energy_parameter(p, E_param, AtomMatrices, Noco, echo);

#ifdef  HAS_TFQMRGPU
              tfqmrgpu::solve(action, memory_buffer, 1e-9, maxiter, 0, true);

              // the 1st part of the memory buffer constains the result Green function
              auto const Green = (real_t const(*)[R1C2][Noco*64][Noco*64]) memory_buffer;
              // extract the density as imaginary part of the trace of the Green function
              simple_stats::Stats<> rho_stats;
              for (unsigned icol = 0; icol < p.nCols; ++icol) {
                  auto const inzb = p.subset[icol]; // index of diagonal blocks
                  for (int i64 = 0; i64 < 64; ++i64) {
                      rho[icol][0][0][i64] = prefactor * Green[inzb][ImaginaryPart][i64][i64];
                      rho_stats.add(rho[icol][0][0][i64]);
                  } // i64
              } // icol
              // ToDo: MPIallreduce rho_stats
              auto const resonance = rho_stats.mean();
              if (echo > 0) std::printf("%.6f %.9f %.1e\n", E_real*eV, resonance, rho_stats.dev());
              if (resonance > max_resonance) { E_resonance = E_real; max_resonance = resonance; }
#else  // HAS_TFQMRGPU
              if (echo > 0) std::printf("# solve for k={%9.6f,%9.6f,%9.6f}, E=(%g, %g) %s\n",
                              k_point[0], k_point[1], k_point[2], E_real*eV, E_imag*eV, _eV);
#endif // HAS_TFQMRGPU
          } // iE
          bandstructure[ik] = E_resonance;
      } // ik
      free_memory(memory_buffer);

      if (echo > 3) {
          std::printf("\n## bandstructure in %s showing peak resonances and free electron energies\n", _eV);
          for (int ik = 0; ik < nkpoints; ++ik) {
              auto const *const k_point = k_path[ik];
              std::printf("%g %g", k_point[3], bandstructure[ik]*eV);
              for (int lat = -2; lat <= 2; ++lat) {
                  auto const E_free_electron = 0.5*(pow2((k_point[0] + lat)*reci[0]) 
                                                  + pow2(k_point[1]*reci[1])
                                                  + pow2(k_point[2]*reci[2]));
                  std::printf(" %g", E_free_electron*eV);
              } // lat
              std::printf("\n");
          } // ik
          std::printf("\n");
      } // echo
#ifndef HAS_TFQMRGPU
      warn("%s needs tfQMRgpu", __func__);
      return -1;
#endif // HAS_TFQMRGPU
      return 0;
  } // bandstructure

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_Green_bandstructure(int const echo=0) {
      uint32_t ng[3] = {0, 0, 0}; // grid sizes
      int8_t   bc[3] = {0, 0, 0}; // boundary conditions
      double   hg[3] = {1, 1, 1}; // grid spacings
      std::vector<double> Veff; // local potential
      int natoms{0}; // number of atoms
      std::vector<double> xyzZinso; // atom info
      std::vector<std::vector<double>> AtomMatrices; // non-local potential

      auto const *const filename = control::get("green_experiments.hamiltonian.file", "Hmt.xml");
      auto const load_stat = green_input::load_Hamiltonian(ng, bc, hg, Veff, natoms, xyzZinso, AtomMatrices, filename, echo - 5);
      if (load_stat) {
          warn("failed to load_Hamiltonian with status=%d", int(load_stat));
          return load_stat;
      } // load_stat

      green_action::plan_t p;
      auto const plan_stat = green_function::construct_Green_function(p, ng, bc, hg, Veff, xyzZinso, AtomMatrices, echo, nullptr, 1);
      if (plan_stat) {
          warn("failed to construct_Green_function with status=%d", int(plan_stat));
          return plan_stat;
      } // plan_stat

      double const box[] = {ng[0]*hg[0], ng[1]*hg[1], ng[2]*hg[2]};
      return bandstructure(p, AtomMatrices, box, echo);
  } // test_Green_bandstructure

  inline status_t all_tests(int const echo=0) {
      int const which = control::get("green_experiments.select.test", -1.);
      status_t stat(0);
      if (which & 0x1) stat += test_Green_bandstructure(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_experiments
