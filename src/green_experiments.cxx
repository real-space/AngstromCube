#include <cstdio> // std::printf
#include <vector> // std::vector<T>

#include "green_experiments.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "simple_timer.hxx" // SimpleTimer
#include "unit_system.hxx" // eV, _eV, Kelvin, _Kelvin
#include "green_input.hxx" // ::load_Hamitonian
#include "data_view.hxx" // view2D<T>
#include "brillouin_zone.hxx" // ::get_kpoint_path
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
#include "green_function.hxx" // ::update_phases, ::construct_Green_function, ::update_energy_parameter
#include "control.hxx" // ::get
#include "linear_algebra.hxx" // ::eigenvalues

namespace green_experiments {

  template <typename real_t=double, int Noco=1>
  status_t bandstructure(
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

              green_function::update_energy_parameter(p, E_param, AtomMatrices, hg[2]*hg[1]*hg[0], Noco, 1.0, echo);

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
#endif // HAS_TFQMRGPU
      return 0;
  } // bandstructure



  template <typename real_t=double, int R1C2=1, int Noco=1>
  size_t inner_products(
        double Hmatrix[][R1C2]
      , double Smatrix[][R1C2]
      , real_t const  psi[][R1C2][Noco*4*4*4][Noco*64]
      , real_t const Hpsi[][R1C2][Noco*4*4*4][Noco*64]
      , real_t const Spsi[][R1C2][Noco*4*4*4][Noco*64]
      , int const nblocks
      , int const nb
      , double const dV=1.0 // volume element
      , int const echo=0 // verbosity level
  )
    // Hmat[i][j] = <psi_i|Hpsi_j> and Smat[i][j] = <psi_i|Spsi_j>
  {
      int constexpr Real = 0, Imag = R1C2 - 1;
      int const nbands = nb*64;

      // integrate over the real space grid
      for (int ib = 0; ib < nb; ++ib) {
          for (int jb = 0; jb < nb; ++jb) {
              for (int ib64 = 0; ib64 < 64; ++ib64) {
                  int const iband = ib*64 + ib64;
                  for (int jb64 = 0; jb64 < 64; ++jb64) {
                      int const jband = jb*64 + jb64;
                      int const ij = iband*nbands + jband;
                      double H_re{0}, H_im{0}, S_re{0}, S_im{0};
                      for (int kzyx4 = 0; kzyx4 < nblocks; ++kzyx4) {
                              int const izyxb = kzyx4*nb + ib;
                              int const jzyxb = kzyx4*nb + jb;
                              for (int kzyx = 0; kzyx < 64; ++kzyx) {

                                    double const psi_re  =  psi[izyxb][Real][kzyx][ib64],
                                                 psi_im  = -psi[izyxb][Imag][kzyx][ib64];

                                    double const Hpsi_re = Hpsi[jzyxb][Real][kzyx][jb64],
                                                 Hpsi_im = Hpsi[jzyxb][Imag][kzyx][jb64];

                                    double const Spsi_re = Spsi[jzyxb][Real][kzyx][jb64],
                                                 Spsi_im = Spsi[jzyxb][Imag][kzyx][jb64];

                                    H_re += psi_re * Hpsi_re; // 2 flop
                                    S_re += psi_re * Spsi_re; // 2 flop
                                    if (Imag) {
                                        H_re -= psi_im * Hpsi_im;                     // 2 flop
                                        S_re -= psi_im * Spsi_im;                     // 2 flop
                                        H_im += psi_re * Hpsi_im + psi_im * Hpsi_re;  // 4 flop
                                        S_im += psi_re * Spsi_im + psi_im * Spsi_re;  // 4 flop
                                    } // is complex
                              } // kzyx
                      } // kzyx4
                      Hmatrix[ij][Real] = H_re*dV; // 1 flop
                      Smatrix[ij][Real] = S_re*dV; // 1 flop
                      if (Imag) {
                          Hmatrix[ij][Imag] = H_im*dV; // 1 flop
                          Smatrix[ij][Imag] = S_im*dV; // 1 flop
                      } // complex
                  } // jb64
              } // ib64
          } // jb
      } // ib

      for (int hs = 0; hs < 2; ++hs) { // check hermitian property
          auto const *const mat = hs ? Hmatrix : Smatrix;
          double dev[] = {0, 0, 0};
          for (int iband = 0; iband < nbands; ++iband) {
              for (int jband = 0; jband < iband; ++jband) {
                  dev[0] += std::abs(mat[iband*nbands + jband][Real] - mat[jband*nbands + iband][Real]);
                  dev[1] += std::abs(mat[iband*nbands + jband][Imag] + mat[jband*nbands + iband][Imag]);
              } // jband triangular loop
              dev[2] += std::abs(mat[iband*nbands + iband][Imag]);
          } // iband
          if (Imag) {
              if (echo > 4) std::printf("# %cmat deviation from hermitian off-diag (%.1e, %.1e), diagonal %.1e\n", hs?'H':'S', dev[0], dev[1], dev[2]);
          } else {
              if (echo > 4) std::printf("# %cmat deviation from symmetric %.2e\n", hs?'H':'S', dev[0]);
          }
      } // check hermitian property

      return nblocks*64ul * 4ul * pow2(R1C2*1ul*nbands); // returns the number of floating point operations
  } // inner_products


  template <typename real_t=double, int R1C2=1, int Noco=1>
  size_t rotate_waves(
        real_t       Rpsi[][R1C2][Noco*4*4*4][Noco*64]
      , real_t const  psi[][R1C2][Noco*4*4*4][Noco*64]
      , double const Rmat[][R1C2] // data layout [(nb*64)(nb*64)][R1C2]
      , int const nblocks
      , int const nb
  )
    // Rpsi_i = sum_j Rmat[i][j] * psi_j
  {
      int constexpr Real = 0, Imag = R1C2 - 1;
      int const nbands = nb*64;

      // sum over jband
      for (int ib = 0; ib < nb; ++ib) {
              for (int ib64 = 0; ib64 < 64; ++ib64) {
                      int const iband = ib*64 + ib64;
                      for (int kzyx4 = 0; kzyx4 < nblocks; ++kzyx4) {
                              int const izyxb = kzyx4*nb + ib;
                              for (int kzyx = 0; kzyx < 64; ++kzyx) {

                                double re{0}, im{0};
                                for (int jb = 0; jb < nb; ++jb) {
                                    int const jzyxb = kzyx4*nb + jb;
                                    for (int jb64 = 0; jb64 < 64; ++jb64) {
                                        int const jband = jb*64 + jb64;

                                        double const psi_re  = psi[jzyxb][Real][kzyx][jb64],
                                                     psi_im  = psi[jzyxb][Imag][kzyx][jb64];

                                        auto const Rmat_re = Rmat[iband*nbands + jband][Real],
                                                   Rmat_im = Rmat[iband*nbands + jband][Imag];

                                    re += Rmat_re * psi_re; // 2 flop
                                    if (Imag) {
                                        re += Rmat_im * psi_im; // 2 flop
                                        im += Rmat_re * psi_im  // 2 flop
                                            - Rmat_im * psi_re; // 2 flop
                                    } // is complex
                                  } // jb64
                                } // jb
                                if (Imag)
                                Rpsi[izyxb][Imag][kzyx][ib64] = im;
                                Rpsi[izyxb][Real][kzyx][ib64] = re;

                              } // k4x k4y k4z
                      } // kx4 ky4 kz4
              } // ib64
      } // ib

      return nblocks*64ul * 2ul * pow2(R1C2*1ul*nbands); // returns the number of floating point operations
  } // rotate_waves

  template <typename real_t=double, int R1C2=1, int Noco=1>
  size_t gradient_waves(
        real_t        psi[][R1C2][Noco*4*4*4][Noco*64]
      , real_t const Hpsi[][R1C2][Noco*4*4*4][Noco*64]
      , real_t const Spsi[][R1C2][Noco*4*4*4][Noco*64]
      , double const Eval[] // eigenvalues
      , int const nblocks
      , int const nb
      , float min_max_res[] // side result
      , int const echo=0
  )
    // psi_j = Hpsi_i - E_i * Spsi_i, j=i+nbhalf
  {
      int const nbands = nb*64; // number of all bands
      min_max_res[0] = 9e9;
      min_max_res[1] = 0.0;

      std::vector<double> res_norm2(nbands, 0.0);
      std::vector<double> psi_norm2(nbands, 0.0);
      for (int iband = 0; iband < nbands; ++iband) {
          int const ib   = iband >> 6; // divide 64
          int const ib64 = iband & 63; // modulo 64

          double norm2{0};
          for (int kzyx4 = 0; kzyx4 < nblocks; ++kzyx4) {
              int const izyxb = kzyx4*nb + ib;
              for (int kzyx = 0; kzyx < 64; ++kzyx) {
                  for (int reim = 0; reim < R1C2; ++reim) {
                      // create the residual vector
                      auto const new_psi = Hpsi[izyxb][reim][kzyx][ib64]
                           - Eval[iband] * Spsi[izyxb][reim][kzyx][ib64]; // 2 flop
                      norm2 += pow2(new_psi); // 2 flop
                  } // reim
              }// kzyx
          } // kzyx4

          if (norm2 <= 0) {
              error("Residual zero for iband=%i", iband);
          } else {
              min_max_res[0] = std::min(min_max_res[0], float(norm2));
              min_max_res[1] = std::max(min_max_res[1], float(norm2));
              res_norm2[iband] = norm2;
          }
      } // iband

      auto threshold2 = std::sqrt(min_max_res[0]*min_max_res[1]); // geometric mean
      int newbands{0};
      if (threshold2 > 1e-30) {
      // filter norms
      int iteration{0};
      newbands = nbands;
      while (newbands*4 > nbands && iteration < 999) { // reduce to less than 1/4 of all bands
          ++iteration;
          threshold2 *= 1.5;
          int isrc{0};
          for (int iband = 0; iband < nbands; ++iband) {
              isrc += (res_norm2[iband] > threshold2);
          } // iband
          newbands = isrc;
          if (0 == (iteration & 0xf) && echo > 3) std::printf("# %s iteration=%i threshold^2=%.1e\n", __func__, iteration, threshold2);
      } // while
      if (echo > 4) std::printf("# %d new bands have been selected with residuals in [%.1e, %.1e]\n", newbands, std::sqrt(threshold2), std::sqrt(min_max_res[1]));
      std::vector<int> i_index(nbands, -1);
      std::vector<int> j_index(nbands, -1);
      int isrc{0}, itrg{0};
      for (int iband = 0; iband < nbands; ++iband) {
          if (res_norm2[iband] > threshold2) {
              i_index[isrc] = iband;
              ++isrc;
          } else {
              j_index[itrg] = iband;
              ++itrg;
          }
      } // iband
      auto const oldbands = itrg;
      newbands = isrc;
      assert(newbands*2 <= nbands && "bisection failed");
      assert(oldbands + newbands == nbands);

      for (int isrc = 0; isrc < std::min(1, newbands); ++isrc) {
          // rescale psi_j
          int const iband = i_index[isrc];
          int const jband = j_index[oldbands - 1 - isrc]; // replace the upper wave functions by residual waves with a large norm
          assert(iband > -1); assert(jband > -1);

          int const jb   = jband >> 6; // divide 64
          int const jb64 = jband & 63; // modulo 64
          int const ib   = iband >> 6; // divide 64
          int const ib64 = iband & 63; // modulo 64

          real_t const f = 1./Eval[iband];
          // double const f = 1./std::sqrt(res_norm2[iband]);
          for (int kzyx4 = 0; kzyx4 < nblocks; ++kzyx4) {
              int const izyxb = kzyx4*nb + ib;
              int const jzyxb = kzyx4*nb + jb;
              for (int kzyx = 0; kzyx < 64; ++kzyx) {
                  for (int reim = 0; reim < R1C2; ++reim) {
                      auto const new_psi = Hpsi[izyxb][reim][kzyx][ib64];
                      // auto const new_psi = Hpsi[izyxb][reim][kzyx][ib64];
                           // - Eval[iband] * Spsi[izyxb][reim][kzyx][ib64]; // 2 flop
                      psi[jzyxb][reim][kzyx][jb64] = f*new_psi; // 1 flop
                  } // reim
              } // kzyx
          } // kzyx4
      } // iband
      } // threshold2 > 1e-30

      for (int mm = 0; mm < 2; ++mm) {
          min_max_res[mm] = std::sqrt(min_max_res[mm]); // export the residuals, not their squares
      } // mm
      return nblocks*64ul * R1C2*(4ul*nbands + 3ul*newbands); // returns the number of floating point operations
  } // gradient_waves


  template <typename real_t=double, int R1C2=1, int Noco=1>
  status_t eigensolver(
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
      int const nbands = nb*64;
      if (echo > 1) std::printf("# %s %d k-points, %d bands\n", __func__, nkpoints, nbands);

      std::vector<double> Sval(nbands, 1.0); // eigenvalues of the overlap operator
      std::vector<double> Eval(nbands, 0.0); // eigenvalues
      std::vector<std::vector<double>> bandstructure(nkpoints, Eval); // result array

      int const ng4[] = {int(ng[0] >> 2), int(ng[1] >> 2), int(ng[2] >> 2)};
      auto const nblocks = size_t(ng4[2]) * size_t(ng4[1]) * size_t(ng4[0]);
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
          pS.echo = echo - 15;
          pH.echo = echo - 15;
      } // scope

      // construct two different action operators
      green_action::action_t<real_t,R1C2,Noco,64> action_H(&pH); // constructor
      green_action::action_t<real_t,R1C2,Noco,64> action_S(&pS); // constructor
      double const dVol = hg[2]*hg[1]*hg[0]; // volume element of the real space grid
      green_function::update_energy_parameter(pH,  0.0, AtomMatrices, dVol, Noco, 1.0, echo); // prepare for H: A = (1*H - (0)*S)
      green_function::update_energy_parameter(pS, -1.0, AtomMatrices, dVol, Noco, 0.0, echo); // prepare for S: A = (0*H - (-1)*S)

      auto  psi = get_memory<real_t[R1C2][Noco*4*4*4][Noco*64]>(nnzb, echo, "waves");

      int constexpr Real = 0, Imag = R1C2 - 1;
      if (nb < nblocks) { // scope: create start wave functions
          // iterative solver has not yet achieved to converge, use the explicit solver with nb == nblocks
          // ToDo: delete failed code
          //
          auto constexpr twopi = 2*constants::pi;
          double const box[] = {ng[0]*hg[0], ng[1]*hg[1], ng[2]*hg[2]}; // in Bohr
          double const reci[] = {twopi/box[0], twopi/box[1], twopi/box[2]}; // reciprocal lattice vectors in Bohr^-1
          auto const recV = reci[0]*reci[1]*reci[2]; // volume of a reciprocal lattice point in Bohr^-3
          // sphere of plane waves: V = 4*constants::pi/3 * radius^3 == nb*64 * recV
          auto const radius = 2.06/R1C2*std::cbrt(nb*64*recV*3/(4*constants::pi)); // in Bohr^-1
          auto const E_cut = pow2(radius); // in Rydberg
          int const npw[] = {int(radius/reci[0]), int(radius/reci[1]), int(radius/reci[2])};
          if (echo > 1) std::printf("# start waves are plane waves cutoff energy is %g Rydberg\n", E_cut);
          auto const E_pw_max = pow2(npw[0]*reci[0]) + pow2(npw[1]*reci[1]) + pow2(npw[2]*reci[2]); // in Rydberg
          if (echo > 1) std::printf("# plane wave box corner energy is %g Rydberg\n", E_pw_max);
          auto const max_npw = (R1C2*npw[2] + 1)*(R1C2*npw[1] + 1)*(R1C2*npw[0] + 1);
          if (echo > 1) std::printf("# check a plane wave box of [-%d,%d] x [-%d,%d] x [-%d,%d] = %.3f k\n", npw[0], npw[0], npw[1], npw[1], npw[2], npw[2], max_npw*.001);
          assert(nb*64 <= max_npw);
          uint32_t const stride = (((max_npw - 1) >> 2) + 1) << 2; // 2: align to 4 doubles
          int ipw{0}, jpw{0}, mpw{0};
          auto kvs = get_memory<double>(3*stride, echo, "plane wave vectors");
          // selection process
          double kvec[3];
          for (int kz = -npw[2]*Imag; kz <= npw[2]; ++kz) {     kvec[2] = kz*reci[2];
            for (int ky = -npw[1]*Imag; ky <= npw[1]; ++ky) {   kvec[1] = ky*reci[1];
              for (int kx = -npw[0]*Imag; kx <= npw[0]; ++kx) { kvec[0] = kx*reci[0];
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
              double const f = 1./std::sqrt(dVol*(ng4[2]*4.*ng4[1]*4.*ng4[0]*4.));
              for (int ib = 0; ib < nb; ++ib) {
                for (int iz4 = 0; iz4 < ng4[2]; ++iz4) {
                for (int iy4 = 0; iy4 < ng4[1]; ++iy4) {
                for (int ix4 = 0; ix4 < ng4[0]; ++ix4) {
                        int const izyxb = ((iz4*ng4[1] + iy4)*ng4[0] + ix4)*nb + ib;
                        for (int i4z = 0; i4z < 4; ++i4z) {     int const iz = iz4*4 + i4z;
                        for (int i4y = 0; i4y < 4; ++i4y) {     int const iy = iy4*4 + i4y;
                        for (int i4x = 0; i4x < 4; ++i4x) {     int const ix = ix4*4 + i4x;
                              int const jzyx = (i4z*4 + i4y)*4 + i4x;
                              for (int ib64 = 0; ib64 < 64; ++ib64) { // threadIdx.x
                                int const ipw = ib*64 + ib64;
                                auto const arg = (ix + .5)*kvs[0*stride + ipw]
                                               + (iy + .5)*kvs[1*stride + ipw]
                                               + (iz + .5)*kvs[2*stride + ipw];
                                if (Imag) {
                                    psi[izyxb][Imag][jzyx][ib64] = f*std::cos(arg);
                                    psi[izyxb][Real][jzyx][ib64] = f*std::sin(arg); // if we treat real wave functions and isolated BCs,
                                              // the sine-solution is the eigenstate of the potential-free particle in a box problem
                                } else {
                                    psi[izyxb][Real][jzyx][ib64] = f*std::cos(arg);
                                }
                              } // ib64
                            }}} // i4x i4y i4z
                    }}} // ix4 iy4 iz4
              } // ib
          }
          free_memory(kvs);

      } else { // scope: start waves
          assert(nb == nblocks); // there are as many bands as real-space grid points
          set(psi[0][0][0], nblocks*nb*R1C2*64ul*64ul, real_t(0)); // clear
          for (int iband = 0; iband < nbands; ++iband) {
              int const ib   = iband >> 6; // divide 64
              int const ib64 = iband & 63; // modulo 64
              psi[ib*nb + ib][Real][ib64][ib64] = 1;
          } // iband
      } // start waves

      auto Hpsi = get_memory<real_t[R1C2][Noco*4*4*4][Noco*64]>(nnzb, echo, "H * waves");
      auto Spsi = get_memory<real_t[R1C2][Noco*4*4*4][Noco*64]>(nnzb, echo, "S * waves");
      auto tpsi = get_memory<real_t[R1C2][Noco*4*4*4][Noco*64]>(nnzb, echo, "temp waves");

      auto Hmat = get_memory<double[R1C2]>(pow2(nbands), echo, "subspace Hamiltonian");
      auto Smat = get_memory<double[R1C2]>(pow2(nbands), echo, "subspace Overlap op");
      auto  mat = get_memory<double[R1C2]>(pow2(nbands), echo, "matrix copy");

      int const maxiter = control::get("green_experiments.eigen.maxiter", (nb == nblocks) ? 1. : 9.);

      simple_stats::Stats<> Gflop_count;
      simple_stats::Stats<> Wtime_count;
      for (int ik = 0; ik < nkpoints; ++ik) {
          SimpleTimer timer(__FILE__, __LINE__, __func__, 0);
          double const *const k_point = k_path[ik];

          if (echo > 3) std::printf("\n## k-point %g %g %g\n", k_point[0], k_point[1], k_point[2]);
          green_function::update_phases(pH, k_point, Noco, echo);
          green_function::update_phases(pS, k_point, Noco, echo/2); // less output since it will print the same as above
          double nops{0};

        int it{0}, lastiter{maxiter - 1};
        for (; it <= lastiter && lastiter >= 0; ++it) { // Davidson iterations

          if (echo > 5) std::printf("# start Davidson iteration #%i\n", it);

          nops += action_H.multiply(Hpsi, psi, colIndex, nnzb, nb);
          nops += action_S.multiply(Spsi, psi, colIndex, nnzb, nb);

          // create inner products <psi_i|Hpsi_j> and <psi_i|Spsi_j>
          nops += inner_products<real_t,R1C2,Noco>(Hmat, Smat,
                                  psi, Hpsi, Spsi, nblocks, nb, dVol, echo);

          set(mat[0], pow2(nbands)*R1C2, Smat[0]); // deep copy of the overlap operator
          // we need a deep copy here because the eigenvalue solver changes the matrix
          status_t stat(0);
          if (2 == R1C2) { // Hermitian generalized eigenvalue problem
              stat = linear_algebra::eigenvalues(Sval.data(), nbands,
                                  (std::complex<double>*)mat, nbands);
          } else {    // real-symmetric generalized eigenvalue problem
              stat = linear_algebra::eigenvalues(Sval.data(), nbands,
                                                (double*)mat, nbands);
          } // real or complex
          if (0 != stat) { // standard eigenvalue problem failed
              warn("failed to diagonalize the overlap for k-point #%i in Davidson iteration #%i", ik, it);
              it = maxiter; // stop
          } else { // standard eigenvalue problem failed
              if (echo > 7) std::printf("# in Davidson iteration #%i overlap eigenvalues:  %g %g %g %g ...\n",
                                                        it, Sval[0], Sval[1], Sval[2], Sval[3]);
              if (Sval[0] > 0.0) {
                // overlap matrix is stable
                if (2 == R1C2) { // Hermitian generalized eigenvalue problem
                    stat = linear_algebra::eigenvalues(Eval.data(), nbands,
                                       (std::complex<double>*)Hmat, nbands,
                                       (std::complex<double>*)Smat, nbands);
                } else {    // real-symmetric generalized eigenvalue problem
                    stat = linear_algebra::eigenvalues(Eval.data(), nbands,
                                                     (double*)Hmat, nbands,
                                                     (double*)Smat, nbands);
                } // real or complex
                if (0 != stat) { // generalized eigenvalue problem failed
                    warn("failed to diagonalize for k-point #%i in Davidson iteration #%i", ik, it);
                    lastiter = it;
                    if (echo > 5) std::printf("# failed to diagonalize in Davidson iteration #%i, set to last iteration\n", it);
                } else { // generalized eigenvalue problem failed
                    if (echo > 6) std::printf("# in Davidson iteration #%i energy eigenvalues:  %g %g %g %g ... %s\n",
                                                            it, Eval[0]*eV, Eval[1]*eV, Eval[2]*eV, Eval[3]*eV, _eV);
                    set(bandstructure[ik].data(), nbands, Eval.data()); // copy
                    // ToDo: rotate 1st half of bands and generate the 2nd half from gradients
                    //        gradient: phi_i = (H - E_i*S) psi_i
                    if (echo > 5) std::printf("# rotate_waves in Davidson iteration #%i\n", it);
                    nops += rotate_waves(tpsi, psi, Hmat, nblocks, nb); // Spsi is a dummy here for a new version of psi
                    std::swap(psi, tpsi); // pointer swap

                    if (Sval[0] > .01) {
                        if (it > 0) {
                          if (echo > 9) std::printf("# gradient_waves in Davidson iteration #%i\n", it);
                          float min_max_res[2];

                          nops += action_H.multiply(Hpsi, psi, colIndex, nnzb, nb);
                          nops += action_S.multiply(Spsi, psi, colIndex, nnzb, nb);
                          nops += gradient_waves(psi, Hpsi, Spsi, Eval.data(), nblocks, nb, min_max_res, echo);
                          if (echo > 5) std::printf("# gradient_waves in iteration #%i has residual norms in [%.1e, %.1e]\n",
                                                    it, min_max_res[0], min_max_res[1]);
                        }
                    } else {
                        if (echo > 5) std::printf("# overlap becomes instable in Davidson iteration #%i\n", it);
                        lastiter = it; // exit
                    } // augment search space by gradients
                } // generalized eigenvalue problem failed

              } else { // overlap matrix is stable
                  if (echo > 7) std::printf("# in Davidson iteration #%i overlap eigenvalues are instable: %g, exit\n", it, Sval[0]);

                  if (echo > 99) {
                    auto const Omat = Smat; char const M = 'S';
                    if (echo > 11) {
                      std::printf("# Matrix %c in the %d x %d subspace of iteration #%i\n", M, nbands, nbands, it);
                      char const *const format[] =  {" %.3f", ",%.2f"}; // {" %.1e", ",%.1e"};
                      for (int iband = 0; iband < nbands; ++iband) {
                          std::printf("#%6i  ", iband);
                          for (int jband = 0; jband < nbands; ++jband) {
                              std::printf(format[0], Omat[iband*nbands + jband][0]);
                              if (R1C2 > 1)
                              std::printf(format[1], Omat[iband*nbands + jband][R1C2 - 1]);
                          } // jband
                          std::printf("\n");
                      } // iband
                    } // echo
                      std::printf("\n# Diagonal elements of Matrix %c in iteration #%i:  ", M, it);
                      for (int iband = 0; iband < nbands; ++iband) {
                          std::printf(" %.6f,%g", Omat[iband*nbands + iband][0], (R1C2 > 1)*Omat[iband*nbands + iband][R1C2 - 1]);
                      } // iband
                      std::printf("\n\n");
                  } // echo

                  lastiter = -maxiter; // stop
              } // overlap matrix is stable

          } // standard eigenvalue problem failed

        } // it Davidson iterations

          if (echo > 3) std::printf("# Davidson method ran %d of max %d iterations\n", it, maxiter);

          if (echo > 0) {
              std::printf("# solve for k={%9.6f,%9.6f,%9.6f}, spectrum(%s) ", k_point[0], k_point[1], k_point[2], _eV);
              printf_vector(" %g", bandstructure[ik], "\n", eV);
          } // echo

          Gflop_count.add(1e-9*nops);
          Wtime_count.add(timer.stop());
      } // ik

      free_memory(Smat); free_memory(Hmat);
      free_memory(Spsi); free_memory(Hpsi);
      free_memory(psi);  free_memory(tpsi);
      free_memory(colIndex);

      if (echo > 0) {
          {   auto const & st = Gflop_count;
              std::printf("\n# %s operations [%g, %g +/- %g, %g] Gflop per k-point, %g Gflop in total\n",
                            __func__, st.min(), st.mean(), st.dev(), st.max(),     st.sum()); }
          std::printf("\n# %s operations %s Gflop per k-point, %g Gflop in total\n", __func__, Gflop_count.interval().c_str(), Gflop_count.sum());
          {   auto const & st = Wtime_count;
              std::printf("# %s needed [%g, %g +/- %g, %g] seconds per k-point, %g seconds in total\n",
                            __func__, st.min(), st.mean(), st.dev(), st.max(),     st.sum()); }
      } // echo

      return 0;
  } // eigensolver




#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_experiment(int const echo=0, char const how='g') {
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
          int const r1c2 = control::get("green_experiments.eigen.real", 0.0);
          if (32 == control::get("green_experiments.eigen.floating.point.bits", 64.0)) {
              return (1 == r1c2) ? eigensolver<float ,1,Noco>(p, pS, AtomMatrices, ng, hg, p.nCols, echo):
                                   eigensolver<float ,2,Noco>(p, pS, AtomMatrices, ng, hg, p.nCols, echo);
          } // single precision
          return (1 == r1c2) ?     eigensolver<double,1,Noco>(p, pS, AtomMatrices, ng, hg, p.nCols, echo):
                                   eigensolver<double,2,Noco>(p, pS, AtomMatrices, ng, hg, p.nCols, echo);
      } // how

  } // test_experiment

  status_t all_tests(int const echo) {
      int const which = control::get("green_experiments.select.test", -1.);
      status_t stat(0);
      if (which & 0x1) stat += test_experiment(echo, 'g'); // Green function spectrum
      if (which & 0x2) stat += test_experiment(echo, 'e'); // eigensolver spectrum
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_experiments
