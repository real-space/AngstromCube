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
#include "linear_algebra.hxx" // ::eigenvalues

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

      if (echo > 1) std::printf("\n# %s<%s,R1C2=%d,Noco=%d>\n", __func__, real_t_name<real_t>(), R1C2, Noco);

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
  inline size_t inner_products(
        double Hmatrix[][R1C2]
      , double Smatrix[][R1C2]
      , real_t const  psi[][R1C2][Noco*4*4*4][Noco*64]
      , real_t const Hpsi[][R1C2][Noco*4*4*4][Noco*64]
      , real_t const Spsi[][R1C2][Noco*4*4*4][Noco*64]
      , int const ng4[3]
      , int const nb
      , double const dV=1.0 // volume element
  ) {
      int constexpr Real = 0, Imag = R1C2 - 1;
      int const nbands = nb*64;

      // integrate over the real space grid
      for (int ib = 0; ib < nb; ++ib) {
          for (int jb = 0; jb < nb; ++jb) {
              for (int ib64 = 0; ib64 < 64; ++ib64) {
                  int const iband = ib*64 + ib64;
                  for (int jb64 = 0; jb64 < 64; ++jb64) {
                      int const jband = jb*64 + jb64;
                      double H_re{0}, H_im{0}, S_re{0}, S_im{0};
                      for (int kz4 = 0; kz4 < ng4[2]; ++kz4) {
                      for (int ky4 = 0; ky4 < ng4[1]; ++ky4) {
                      for (int kx4 = 0; kx4 < ng4[0]; ++kx4) {
                              int const kzyxb = ((kz4*ng4[1] + ky4)*ng4[0] + kx4)*nb + ib;
                              for (int k4z = 0; k4z < 4; ++k4z) {
                              for (int k4y = 0; k4y < 4; ++k4y) {
                              for (int k4x = 0; k4x < 4; ++k4x) {
                                    int const kzyx = (k4z*4 + k4y)*4 + k4x;

                                    auto const psi_re  =  psi[kzyxb][Real][kzyx][ib64],
                                               psi_im  =  psi[kzyxb][Imag][kzyx][ib64];

                                    auto const Hpsi_re = Hpsi[kzyxb][Real][kzyx][jb64],
                                               Hpsi_im = Hpsi[kzyxb][Imag][kzyx][jb64];

                                    auto const Spsi_re = Spsi[kzyxb][Real][kzyx][jb64],
                                               Spsi_im = Spsi[kzyxb][Imag][kzyx][jb64];

                                    H_re += psi_re * Hpsi_re; // 2 flop
                                    S_re += psi_re * Spsi_re; // 2 flop
                                    if (Imag) {
                                        H_re += psi_im * Hpsi_im;                     // 2 flop
                                        S_re += psi_im * Spsi_im;                     // 2 flop
                                        H_im += psi_re * Hpsi_im - psi_im * Hpsi_re;  // 4 flop
                                        S_im += psi_re * Spsi_im - psi_im * Spsi_re;  // 4 flop
                                    } // is complex
                              }}} // k4x k4y k4z
                      }}} // kx4 ky4 kz4
                      Hmatrix[iband*nbands + jband][Real] = H_re*dV; // 1 flop
                      Smatrix[iband*nbands + jband][Real] = S_re*dV; // 1 flop
                      if (Imag) {
                          Hmatrix[iband*nbands + jband][Imag] = H_im*dV; // 1 flop
                          Smatrix[iband*nbands + jband][Imag] = S_im*dV; // 1 flop
                      } // complex
                  } // jb64
              } // ib64
          } // jb
      } // ib

      return ng4[2]*4ul*ng4[1]*4ul*ng4[0]*4ul * 4ul * pow2(R1C2*1ul*nbands); // returns the number of floating point operations
  } // inner_products

  template <typename real_t=double, int R1C2=1, int Noco=1>
  inline status_t eigensolver(
        green_action::plan_t & pH
      , green_action::plan_t & pS
      , std::vector<std::vector<double>> const & AtomMatrices
      , uint32_t const ng[3] // grid points
      , double const hg[3] // grid spacings
      , int const nb=1 // number of bands == 64*nb
      , int const echo=0
  ) {
      if (echo > 1) std::printf("\n# %s<%s,R1C2=%d,Noco=%d>\n", __func__, real_t_name<real_t>(), R1C2, Noco);
      assert(1 == Noco && "Not prepared for Noco");
      for (int d = 0; d < 3; ++d) assert(0 == (ng[d] & 0x3)); // all grid numbers must be a multiple of 4

      if (echo > 1) std::printf("# number of bands %d = %d * 64\n", nb*64, nb);

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

      view2D<double> k_path;
      auto const nkpoints = brillouin_zone::get_kpoint_path(k_path, echo);
      if (echo > 1) std::printf("# %s %d k-points, %d bands\n", __func__, nkpoints, nb*64);

      std::vector<std::vector<double>> bandstructure(nkpoints, std::vector<double>(nb*64, -9e9)); // result array

      int const ng4[] = {int(ng[0] >> 2), int(ng[1] >> 2), int(ng[2] >> 2)};
      auto const nblocks = (ng4[2]) * size_t(ng4[1]) * size_t(ng4[0]);
      if (echo > 1) std::printf("# %s cell grid has %d x %d x %d = %ld blocks\n", __func__, ng4[2], ng4[1], ng4[0], nblocks);
      size_t const nnzb = nblocks * nb;

      // create a trivial list
      assert(nb <= (1ul << 16) && "too many bands for using uint16_t");
      auto colIndex = get_memory<uint16_t>(nnzb, echo, "colIndex");
      for (int iblock = 0; iblock < nblocks; ++iblock) {
          for (int ib = 0; ib < nb; ++ib) {
              int const inzb = iblock*nb + ib;
              colIndex[inzb] = ib;
          } // ib
      } // iblock

      if (1) { // scope: up to now, pH and pS are the same plan. We need to modify pS now
          // delete the kinetic energy lists
          for (int d = 0; d < 3; ++d) {
              pS.kinetic_plan[d] = green_sparse::sparse_t<int32_t>(); // standard constructor
              if (echo > 0) std::printf("# %s modified pS.kinetic_plan[%c].nRows() = %d\n", __func__, 'x'+d, pS.kinetic_plan[d].nRows());
          } // d
          for (int mag = 0; mag < 4; ++mag) {
              free_memory(pS.Veff[mag]);
              pS.Veff[0] = get_memory<double[64]>(1, echo, "unity instead of potential");
              for (int izyx = 0; izyx < 64; ++izyx) {
                  pS.Veff[0][0][izyx] = double(0 == mag); // replace local potential by unity operation
              } // izyx
          } // mag
          set(pS.veff_index, nnzb, 0); // overwrite veff_index by all zeros.
          pS.nCols = nb;
          pH.nCols = nb;
          pS.echo = echo - 5;
          pH.echo = echo - 5;
      } // scope

      green_action::action_t<real_t,R1C2,Noco,64> action_H(&pH); // constructor
      green_action::action_t<real_t,R1C2,Noco,64> action_S(&pS); // constructor
      green_function::update_energy_parameter(pH,  0.0, AtomMatrices, Noco, 1.0, echo); // prepare for H: A = (1*H - (0)*S)
      green_function::update_energy_parameter(pS, -1.0, AtomMatrices, Noco, 0.0, echo); // prepare for S: A = (0*H - (-1)*S)

      auto  psi = get_memory<real_t[R1C2][Noco*4*4*4][Noco*64]>(nnzb, echo, "waves");

      { // create start wave functions
          auto constexpr twopi = 2*constants::pi;
          double const box[] = {ng[0]*hg[0], ng[1]*hg[1], ng[2]*hg[2]}; // in Bohr
          double const reci[] = {twopi/box[0], twopi/box[1], twopi/box[2]}; // reciprocal lattice vectors in Bohr^-1
          double const recV = reci[0]*reci[1]*reci[2]; // volume of a reciprocal lattice point in Bohr^-3
          // sphere of plane waves: V = 4*constants::pi/3 * radius^3 == nb*64 * recV
          double const radius = 1.03*std::cbrt(nb*64*recV*3/(4*constants::pi)); // in Bohr^-1
          auto const E_cut = pow2(radius); // in Rydberg
          int const npw[] = {int(radius/reci[0]), int(radius/reci[1]), int(radius/reci[2])};
          if (echo > 1) std::printf("# start waves are plane waves cutoff energy is %g Rydberg\n", E_cut);
          auto const E_pw_max = pow2(npw[0]*reci[0]) + pow2(npw[1]*reci[1]) + pow2(npw[2]*reci[2]); // in Rydberg
          if (echo > 1) std::printf("# plane wave box corner energy is %g Rydberg\n", E_pw_max);
          auto const max_npw = (2*npw[2] + 1)*(2*npw[1] + 1)*(2*npw[0] + 1);
          if (echo > 1) std::printf("# check a plane wave box of [-%d,%d] x [-%d,%d] x [-%d,%d] = %.3f k\n", npw[0], npw[0], npw[1], npw[1], npw[2], npw[2], max_npw*.001);
          assert(nb*64 <= max_npw);
          uint32_t const stride = (((max_npw - 1) >> 2) + 1) << 2; // 2: align to 4 doubles
          int ipw{0}, jpw{0}, mpw{0};
          auto kvs = get_memory<double>(3*stride, echo, "plane wave vectors");
          // selection process
          double kvec[3];
          for (int kz = -npw[2]; kz <= npw[2]; ++kz) {     kvec[2] = kz*reci[2];
            for (int ky = -npw[1]; ky <= npw[1]; ++ky) {   kvec[1] = ky*reci[1];
              for (int kx = -npw[0]; kx <= npw[0]; ++kx) { kvec[0] = kx*reci[0];
                  auto const E_pw = pow2(kvec[0]) + pow2(kvec[1]) + pow2(kvec[2]); // in Rydberg
                  if (E_pw <= E_cut) {
                      ++jpw; // count a plane wave inside the sphere
                      if (ipw < nb*64) {
                        for (int d = 0; d < 3; ++d) {
                            kvs[d*stride + ipw] = kvec[d]*hg[d];
                        } // d
                        ++ipw; // count a plane wave that should be stored in psi
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
          {
              int constexpr Real = 0, Imag = R1C2 - 1;
              for (int ib = 0; ib < nb; ++ib) {
                for (int iz4 = 0; iz4 < ng4[2]; ++iz4) {
                for (int iy4 = 0; iy4 < ng4[1]; ++iy4) {
                for (int ix4 = 0; ix4 < ng4[0]; ++ix4) {
                        int const izyxb = ((iz4*ng4[1] + iy4)*ng4[0] + ix4)*nb + ib;
                        for (int i4z = 0; i4z < 4; ++i4z) {     int const iz = iz4*4 + i4z;
                        for (int i4y = 0; i4y < 4; ++i4y) {     int const iy = iy4*4 + i4y;
                        for (int i4x = 0; i4x < 4; ++i4x) {     int const ix = ix4*4 + i4x;
                              int const jzyx = (i4z*4 + i4y)*4 + i4x;
                              for (int jb = 0; jb < 64; ++jb) { // threadIdx.x
                                int const ipw = ib*64 + jb;
                                auto const arg = ix*kvs[0*stride + ipw] + iy*kvs[1*stride + ipw] + iz*kvs[2*stride + ipw];
                                if (Imag)
                                psi[izyxb][Imag][jzyx][jb] = std::sin(arg);
                                psi[izyxb][Real][jzyx][jb] = std::cos(arg);
                              } // jb
                            }}} // i4x i4y i4z
                    }}} // ix4 iy4 iz4
              } // ib
          }
          free_memory(kvs);

      } // scope: start waves

      auto Hpsi = get_memory<real_t[R1C2][Noco*4*4*4][Noco*64]>(nnzb, echo, "H * waves");
      auto Spsi = get_memory<real_t[R1C2][Noco*4*4*4][Noco*64]>(nnzb, echo, "S * waves");

      int const nbands = nb*64;
      auto Hmat = get_memory<double[R1C2]>(pow2(nbands), echo, "subspace Hamiltonian");
      auto Smat = get_memory<double[R1C2]>(pow2(nbands), echo, "subspace Overlap op");
      double dVol = hg[2]*hg[1]*hg[0]; // volumen element of the real space grid

      simple_stats::Stats<> Gflop_count;
      for (int ik = 0; ik < nkpoints; ++ik) {
          double const *const k_point = k_path[ik];

          if (echo > 0) std::printf("\n## k-point %g %g %g\n", k_point[0], k_point[1], k_point[2]);
          green_function::update_phases(pH, k_point, Noco, echo);
          green_function::update_phases(pS, k_point, Noco, echo);
          double nops{0};

          nops += action_H.multiply(Hpsi,  psi, colIndex, nnzb, nb);
          nops += action_S.multiply(Spsi,  psi, colIndex, nnzb, nb);

          // create inner products <psi_i|Hpsi_j> and <psi_i|Spsi_j>
          nops += inner_products<real_t,R1C2,Noco>(Hmat, Smat,
                                  psi, Hpsi, Spsi, ng4, nb, dVol);

          status_t stat(0);
          if (2 == R1C2) {
              // Hermitian generalized eigenvalue problem
              stat = linear_algebra::eigenvalues(bandstructure[ik].data(), nbands,
                                              (std::complex<double>*)Hmat, nbands,
                                              (std::complex<double>*)Smat, nbands);
          } else {
              // real symmetric generalized eigenvalue problem
              stat = linear_algebra::eigenvalues(bandstructure[ik].data(), nbands,
                                                            (double*)Hmat, nbands,
                                                            (double*)Smat, nbands);
          } // real or complex
          if (0 != stat) {
              warn("failed to diagonalize for k-point #%i", ik);
          } else {
              // ToDo: rotate 1st half of bands and generate the 2nd half from gradients
              //        gradient: phi_i = (H - E_i*S) psi_i
          }


          if (echo > 0) {
              std::printf("# solve for k={%9.6f,%9.6f,%9.6f}, spectrum(%s) ", k_point[0], k_point[1], k_point[2], _eV);
              printf_vector(" %g", bandstructure[ik], "\n", eV);
          } // echo

          Gflop_count.add(1e-9*nops);
      } // ik

      free_memory(Smat); free_memory(Hmat);
      free_memory(Spsi); free_memory(Hpsi);
      free_memory(psi);
      free_memory(colIndex);

      if (echo > 0) {
          auto const & st = Gflop_count;
          std::printf("\n# %s operations [%g, %g +/- %g, %g] Gflop per k-point, %g Gflop in total\n",
                          __func__, st.min(), st.mean(), st.dev(), st.max(), st.sum());
      } // echo

      return 0;
  } // eigensolver




#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_experiment(int const echo=0, char const how='g') {
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

      if ('g' != how) {
          double const huge = 9*std::max(std::max(ng[0]*hg[0], ng[1]*hg[1]), ng[2]*hg[2]);
          char string[32]; std::snprintf(string, 32, "%g", huge);
          control::set("green_function.truncation.radius", string);
      } // how

      int constexpr Noco=1;
      green_action::plan_t p;
      auto const plan_stat = green_function::construct_Green_function(p, ng, bc, hg, Veff, xyzZinso, AtomMatrices, echo, nullptr, Noco);
      if (plan_stat) {
          warn("failed to construct_Green_function with status=%d", int(plan_stat));
          return plan_stat;
      } // plan_stat

      if ('g' == how) {
          // compute the bandstructure using the Green function method
          return bandstructure<double,Noco>(p, AtomMatrices, ng, hg, echo);
      } else {
          // compute a bandstructure using a wave function method
          green_action::plan_t pS; // plan for the overlap operator
          int const echo_pS = echo*control::get("green_experiments.overlap.echo", 0.0);
          auto const plan_stat = green_function::construct_Green_function(pS, ng, bc, hg, Veff, xyzZinso, AtomMatrices, echo_pS, nullptr, Noco);
          if (plan_stat) {
              warn("failed to construct_Green_function with status=%d for the overlap operator", int(plan_stat));
              return plan_stat;
          } // plan_stat
          return (1 == control::get("green_experiments.eigen.real", 0.0)) ?
               eigensolver<double,1,Noco>(p, pS, AtomMatrices, ng, hg, p.nCols, echo):
               eigensolver<double,2,Noco>(p, pS, AtomMatrices, ng, hg, p.nCols, echo);
      } // how

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
