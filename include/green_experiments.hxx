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

#ifdef    HAS_TFQMRGPU

    #ifdef    HAS_NO_CUDA
        #include "tfQMRgpu/include/tfqmrgpu_cudaStubs.hxx" // cuda... (dummies)
        #define devPtr const __restrict__
    #else  // HAS_NO_CUDA
        #include <cuda.h>
    #endif // HAS_NO_CUDA
    #include "tfQMRgpu/include/tfqmrgpu_memWindow.h" // memWindow_t
    #include "tfQMRgpu/include/tfqmrgpu_linalg.hxx" // ...
    #include "tfQMRgpu/include/tfqmrgpu_core.hxx" // tfqmrgpu::solve<action_t>

#endif // HAS_TFQMRGPU

#include "green_memory.hxx" // get_memory, free_memory
#include "green_action.hxx" // ::plan_t, ::action_t
#include "green_function.hxx" // ::update_phases, ::construct_Green_function
#include "control.hxx" // ::get

namespace green_experiments {

  double const GByte = 1e-9; char const *const _GByte = "GByte";

  template <typename real_t=double, int Noco=1>
  inline status_t bandstructure(
        green_action::plan_t & p
      , std::vector<std::vector<double>> const & AtomMatrices
      , uint32_t const ng[3] // grid points
      , double const hg[3] // grid spacings
      , int const echo=0
  ) {
      int constexpr R1C2 = 2;

      if (echo > 1) std::printf("# %s<%s,R1C2=%d,Noco=%d>\n", __func__, real_t_name<real_t>(), R1C2, Noco);

      green_action::action_t<real_t,R1C2,Noco,64> action(&p); // constructor


      p.gpu_mem = 0;
#ifdef    HAS_TFQMRGPU
      p.echo = echo - 5;
      if (echo > 0) std::printf("\n# call tfqmrgpu::mem_count\n");
      tfqmrgpu::solve(action); // try to instanciate tfqmrgpu::solve with this action_t<real_t,R1C2,Noco,64>
      if (echo > 5) std::printf("# tfqmrgpu::solve requires %.6f GByte GPU memory\n", p.gpu_mem*1e-9);
      double constexpr prefactor = 1./constants::pi;
      int constexpr ImaginaryPart = 1;
      assert(1 == Noco);
      auto rho = get_memory<double[Noco][Noco][64]>(p.nCols, echo, "rho");
      set(rho[0][0][0], p.nCols*Noco*Noco*64, 0.0);
      int const maxiter = control::get("tfqmrgpu.max.iterations", 99.);
      if (echo > 3) std::printf("# +tfqmrgpu.max.iterations=%d\n", maxiter);
#endif // HAS_TFQMRGPU
      auto memory_buffer = get_memory<char>(p.gpu_mem, echo, "tfQMRgpu-memoryBuffer");

      auto const E0     = control::get("green_experiments.bandstructure.energy.offset", 0.0);
      auto const dE     = control::get("green_experiments.bandstructure.energy.spacing", 0.01);
      int const nE      = control::get("green_experiments.bandstructure.energy.points", 1.);
      auto const E_imag = control::get("green_experiments.bandstructure.energy.imag", 1e-3); // 1e-3==room temperature
      if (E_imag <= 0) warn("imaginary part of the energy parameter is %.1e Ha", E_imag);

      if (echo > 1) std::printf("# %d E-points in [%g, %g] %s\n", nE, E0*eV, (E0 + (nE - 1)*dE)*eV, _eV);

      view2D<double> k_path;
      auto const nkpoints = brillouin_zone::get_kpoint_path(k_path, echo);
      if (echo > 1) std::printf("# %s %d k-points, %d E-points, Temperature %g %s = %g %s\n",
                            __func__, nkpoints, nE, E_imag*Kelvin, _Kelvin, E_imag*eV, _eV);

      std::vector<double> bandstructure(nkpoints, -9e9);

      for (int ik = 0; ik < nkpoints; ++ik) {
          double const *const k_point = k_path[ik];

          if (echo > 0) std::printf("\n## k-point %g %g %g\n", k_point[0], k_point[1], k_point[2]);
          green_function::update_phases(p, k_point, Noco, echo);

          double E_resonance{-9};
#ifdef    HAS_TFQMRGPU
          double max_resonance{-9e9};
#endif // HAS_TFQMRGPU
          for (int iE = 0; iE < nE; ++iE) {
              double const E_real = iE*dE + E0;
              std::complex<double> E_param(E_real, E_imag);

              green_function::update_energy_parameter(p, E_param, AtomMatrices, Noco, 1.0, echo);

#ifdef    HAS_TFQMRGPU
              if (maxiter >= 0) {
                  tfqmrgpu::solve(action, memory_buffer, 1e-9, maxiter, 0, true);
              } else {
                  if(echo > 6) std::printf("# skip tfqmrgpu::solve due to maxiter=%d\n", maxiter);
              }

              // the 1st part of the memory buffer constains the result Green function
              auto const Green = (real_t const(*)[R1C2][Noco*64][Noco*64]) memory_buffer;
              // extract the density as imaginary part of the trace of the Green function
              simple_stats::Stats<> rho_stats;
              for (unsigned icol = 0; icol < p.nCols; ++icol) {
                  auto const inzb = p.subset[icol]; // index of a diagonal block
                  for (int i64 = 0; i64 < 64; ++i64) {
                      rho[icol][0][0][i64] = prefactor * Green[inzb][ImaginaryPart][i64][i64];
                      rho_stats.add(rho[icol][0][0][i64]);
                  } // i64
              } // icol
              // ToDo: MPIallreduce rho_stats
              auto const resonance = rho_stats.mean();
              if (echo > 0) std::printf("%.6f %.9f %.1e\n", E_real*eV, resonance, rho_stats.dev());
              if (resonance > max_resonance) { E_resonance = E_real; max_resonance = resonance; }

              auto const pGp = green_dyadic::get_projection_coefficients<real_t,R1C2,Noco>(Green, p.dyadic_plan, p.rowindx, p.rowCubePos, p.colCubePos, echo);
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
          auto constexpr twopi = 2*constants::pi;
          double const box[] = {ng[0]*hg[0], ng[1]*hg[1], ng[2]*hg[2]};
          double const reci[] = {twopi/box[0], twopi/box[1], twopi/box[2]};
          for (int ik = 0; ik < nkpoints; ++ik) {
              auto const *const k_point = k_path[ik];
              std::printf("%g %g", k_point[3], bandstructure[ik]*eV);
              for (int lat = -2; lat <= 2; ++lat) {
                  auto const E_free_electron = 0.5*(pow2((k_point[0] + lat)*reci[0])
                                                  + pow2( k_point[1]       *reci[1])
                                                  + pow2( k_point[2]       *reci[2]));
                  std::printf(" %g", E_free_electron*eV);
              } // lat
              std::printf("\n");
          } // ik
          std::printf("\n");
      } // echo
#ifdef    HAS_TFQMRGPU
      free_memory(rho);
#else  // HAS_TFQMRGPU
      warn("%s needs tfQMRgpu", __func__);
      return -1;
#endif // HAS_TFQMRGPU
      return 0;
  } // bandstructure





  template <typename real_t=double, int R1C2=1, int Noco=1>
  inline status_t eigensolver(
        green_action::plan_t & p
      , std::vector<std::vector<double>> const & AtomMatrices
      , uint32_t const ng[3] // grid points
      , double const hg[3] // grid spacings
      , int const echo=0
  ) {
      if (echo > 1) std::printf("# %s<%s,R1C2=%d,Noco=%d>\n", __func__, real_t_name<real_t>(), R1C2, Noco);
      for (int d = 0; d < 3; ++d) assert(0 == (ng[d] & 0x3)); // all grid numbers must be a multiple of 4

      // This module is part of the Green function code, hence its name prefix green_
      //
      // The idea of this module is to test the implementation of the action
      // against an eigensolver.
      // To do this, we need to create the action_t without truncation
      // in a cell with periodic or isolated BCs (or combinations of those)
      // and set the truncation radius to much larger than the cell extent.
      // Further, we need two action_t instances, the Hamiltonian and the overlap.

      // Then, the number of bands replaces the total number of RHS points.
      // Due to the data layout, the number of bands is always a multiple of 64.
      // i.e. we can accommodate at least 128 electrons.

      // The two action_t operators, Hmt and Ovl, are passed to a Davidson solver
      // which needs start wave functions as it is an iterative algorithm.

      // We need to implement GPU kernels to create the inner product of two sets
      // of wave functions in order to create the matrix elements.

      // The sets of wave functions would preferably be in a data layout like
      //        real_t Psi[nbands][nz*ny*nx][R1C2]; or Psi[nz*ny*nx][nbands][R1C2]
      // However, to use the GPU-action_t, the need to stay in the data layout
      //        real_t Psi[nnzb][R1C2][4*4*4][64];
      // with nnzb == (nbands/64) * (nz/4)*(ny/4)*(nx/4)
      // for a cell with nx*ny*nz grid points in total.
      // The split dimensions make it less straightforward to apply a rotation
      // to the wave function. Probably a custom code is necessary here, too.

      // R1C2 can be either 1 (   real symmetric generalized eigenvalue problem)
      //                 or 2 (complex Hermitian generalized eigenvalue problem)

      // For the initialization, we can use a wrapper that mimiques the data layout.
      // Anyway, we need an initialization of always at least 64 states.
      // Maybe plane waves could help if not enough atomic states are available.


      green_action::action_t<real_t,R1C2,Noco,64> action(&p); // constructor

      int const nbands = control::get("green_experiments.bandstructure.bands", 64.0);
      int const nb = (nbands - 1)/64 + 1;
      if (echo > 1) std::printf("# number of bands %d --> %d = %d * 64\n", nbands, nb*64, nb);
      assert(nbands <= nb*64);

      view2D<double> k_path;
      auto const nkpoints = brillouin_zone::get_kpoint_path(k_path, echo);
      if (echo > 1) std::printf("# %s %d k-points, %d bands\n", __func__, nkpoints, nb*64);

      std::vector<std::vector<double>> bandstructure(nkpoints, std::vector<double>(nbands, -9e9)); // result array

      int const ng4[] = {int(ng[0] >> 2), int(ng[1] >> 2), int(ng[2] >> 2)};
      auto const nblocks = (ng4[2]) * size_t(ng4[1]) * size_t(ng4[0]);
      if (echo > 1) std::printf("# %s data structure has %d x %d x %d = %ld blocks\n", __func__, ng4[2], ng4[1], ng4[0], nblocks);
      size_t const nnzb = nblocks * nb;
      auto  psi = get_memory<real_t[R1C2][Noco*4*4*4][Noco*64]>(nnzb, echo, "waves");

      { // create start wave functions
          auto constexpr twopi = 2*constants::pi;
          double const box[] = {ng[0]*hg[0], ng[1]*hg[1], ng[2]*hg[2]}; // in Bohr
          double const reci[] = {twopi/box[0], twopi/box[1], twopi/box[2]}; // reciprocal lattice vectors in Bohr^-1
          double const recV = reci[0]*reci[1]*reci[2]; // volume of a reciprocal lattice point in Bohr^-3
          // sphere of plane waves: V = 4*constants::pi/3 * radius^3 == nb*64 * recV
          double const radius = 1.01*std::cbrt(nb*64*recV*3/(4*constants::pi)); // in Bohr^-1
          auto const E_cut = pow2(radius); // in Rydberg
          int const npw[] = {int(radius/reci[0]), int(radius/reci[1]), int(radius/reci[2])};
          if (echo > 1) std::printf("# start waves are plane waves cutoff energy is %g Rydberg\n", E_cut);
          auto const E_pw_max = pow2(npw[0]*reci[0]) + pow2(npw[1]*reci[1]) + pow2(npw[2]*reci[2]); // in Rydberg
          if (echo > 1) std::printf("# plane wave box corner energy is %g Rydberg\n", E_pw_max);
          auto const max_npw = (2*npw[2] + 1)*(2*npw[1] + 1)*(2*npw[0] + 1);
          if (echo > 1) std::printf("# check a plane wave box of [-%d,%d] x [-%d,%d] x [-%d,%d] = %.3f k\n", npw[0], npw[0], npw[1], npw[1], npw[2], npw[2], max_npw*.001);
          assert(nb*64 <= max_npw);
          int ipw{0}, jpw{0}, mpw{0};
          bool constexpr create_now = false;
          auto kvs = get_memory<double[4]>(max_npw, echo, "plane wave vectors");
          double kvec[3];
          for (int kz = -npw[2]; kz <= npw[2]; ++kz) {     kvec[2] = kz*reci[2];
            for (int ky = -npw[1]; ky <= npw[1]; ++ky) {   kvec[1] = ky*reci[1];
              for (int kx = -npw[0]; kx <= npw[0]; ++kx) { kvec[0] = kx*reci[0];
                  auto const E_pw = pow2(kvec[0]) + pow2(kvec[1]) + pow2(kvec[2]); // in Rydberg
                  if (E_pw <= E_cut) {
                      ++jpw; // count a plane wave inside the sphere
                      if (ipw < nb*64) {
                        double const kv[] = {kvec[0]*hg[0], kvec[1]*hg[1], kvec[2]*hg[2], E_pw};
                        if (create_now) {
                          int const ib = ipw >> 6;   // divide 64
                          int const jb = ipw & 0x3f; // modulo 64
                          assert(ib*64 + jb == ipw);
                          // create the plane wave with wave vector kvec
                          for (int iz = 0; iz < ng[2]; ++iz) {      int const iz4 = iz >> 2, i4z = iz & 0x3;
                            for (int iy = 0; iy < ng[1]; ++iy) {    int const iy4 = iy >> 2, i4y = iy & 0x3;
                              for (int ix = 0; ix < ng[0]; ++ix) {  int const ix4 = ix >> 2, i4x = ix & 0x3;
                                  auto const arg = ix*kv[0] + iy*kv[1] + iz*kv[2];
                                  int const jzyx = (i4z*4 + i4y)*4 + i4x;
                                  int const izyxb = ((iz4*ng4[1] + iy4)*ng4[0] + ix4)*nb + ib;
                                  if (2 == R1C2)
                                  psi[izyxb][1][jzyx][jb] = std::sin(arg);
                                  psi[izyxb][0][jzyx][jb] = std::cos(arg);
                              } // ix
                            } // iy
                          } // iz
                        } // create_now
                        for (int d = 0; d < 4; ++d) {
                            kvs[ipw][d] = kv[d];
                        } // d
                        ++ipw; // count a plan wave that has been stored in psi
                      } // accept plane wave
                  } // plane wave inside sphere
                  ++mpw;
              } // kx
            } // ky
          } // kz
          assert(max_npw == mpw);
          if (echo > 1) std::printf("# start waves are %d of max %d (sphere) or %d (box) plane waves\n", ipw, jpw, mpw);
          auto const npws = ipw;
          assert(npws >= nb*64 && "Need to have enough plane waves to fill all bands");
          if (!create_now) {
              for (int ib = 0; ib < nb; ++ib) {
                for (int iz = 0; iz < ng[2]; ++iz) {      int const iz4 = iz >> 2, i4z = iz & 0x3;
                  for (int iy = 0; iy < ng[1]; ++iy) {    int const iy4 = iy >> 2, i4y = iy & 0x3;
                    for (int ix = 0; ix < ng[0]; ++ix) {  int const ix4 = ix >> 2, i4x = ix & 0x3;
                        int const izyxb = ((iz4*ng4[1] + iy4)*ng4[0] + ix4)*nb + ib;
                        int const jzyx = (i4z*4 + i4y)*4 + i4x;
                        for (int jb = 0; jb < 64; ++jb) {
                            int const ipw = ib*64 + jb;
                            auto const *const kv = kvs[ipw];
                            auto const arg = ix*kv[0] + iy*kv[1] + iz*kv[2];
                            if (2 == R1C2)
                            psi[izyxb][1][jzyx][jb] = std::sin(arg);
                            psi[izyxb][0][jzyx][jb] = std::cos(arg);
                        } // jb
                    } // ix
                  } // iy
                } // iz
              } // ib
          } // not create_now
          free_memory(kvs); // list of plane waves
          return 0;

      } // scope: start waves

      auto Hpsi = get_memory<real_t[R1C2][Noco*4*4*4][Noco*64]>(nnzb, echo, "H * waves");
      auto Spsi = get_memory<real_t[R1C2][Noco*4*4*4][Noco*64]>(nnzb, echo, "S * waves");

      for (int ik = 0; ik < nkpoints; ++ik) {
          double const *const k_point = k_path[ik];

          if (echo > 0) std::printf("\n## k-point %g %g %g\n", k_point[0], k_point[1], k_point[2]);
          green_function::update_phases(p, k_point, Noco, echo);

          green_function::update_energy_parameter(p,  0.0, AtomMatrices, Noco, 1.0, echo); // for H
          green_function::update_energy_parameter(p, -1.0, AtomMatrices, Noco, 0.0, echo); // for S

          if (echo > 0) {
              std::printf("# solve for k={%9.6f,%9.6f,%9.6f}, spectrum(%s) ", k_point[0], k_point[1], k_point[2], _eV);
              printf_vector(" %g", bandstructure[ik], "\n", eV);
          } // echo
      } // ik

      free_memory(Spsi); free_memory(Hpsi);
      free_memory(psi);

      return 0;
  } // eigensolver




#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_experiment(int const echo=0, char const what='g') {
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

      int constexpr Noco=1;
      green_action::plan_t p;
      auto const plan_stat = green_function::construct_Green_function(p, ng, bc, hg, Veff, xyzZinso, AtomMatrices, echo, nullptr, Noco);
      if (plan_stat) {
          warn("failed to construct_Green_function with status=%d", int(plan_stat));
          return plan_stat;
      } // plan_stat

      if ('g'  == what) {
          return bandstructure(p, AtomMatrices, ng, hg, echo);
      } else {
          return eigensolver(p, AtomMatrices, ng, hg, echo);
      }

  } // test_experiment

  inline status_t all_tests(int const echo=0) {
      int const which = control::get("green_experiments.select.test", -1.);
      status_t stat(0);
      if (which & 0x1) stat += test_experiment(echo, 'g'); // Green function spectrum
      if (which & 0x2) stat += test_experiment(echo, 'e'); // eigensolver spectrum
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_experiments
