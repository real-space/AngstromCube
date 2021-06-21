#pragma once

#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cmath> // std::sqrt
#include <algorithm> // std::max
#include <utility> // std::swap //, std::move
#include <vector> // std::vector<T>
#include <cstdio> // std::printf

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "simple_timer.hxx" // SimpleTimer

#include "constants.hxx" // ::pi
#include "green_memory.hxx" // get_memory, free_memory
#include "green_kinetic.hxx" // ::finite_difference_plan_t
#include "simple_stats.hxx" // ::Stats<>

#ifdef HAS_TFQMRGPU
    #include "tfQMRgpu/tfqmrgpu_core.hxx" // solve<action_t>
    #include "tfQMRgpu/tfqmrgpu_memWindow.h" // memWindow_t
#else
    #include <utility> // std::pair<T>
    typedef std::pair<size_t,size_t> memWindow_t;
    typedef size_t cudaStream_t;
#endif // HAS_TFQMRGPU

namespace green_action {

  typedef struct {
      double pos[3];
      double sigma;
      int32_t gid;
      int32_t ia; // local atom index
      int16_t shifts[3];
      uint8_t nc; // number of coefficients
      int8_t numax;
  } atom_t;


  struct plan_t {

      // for the inner products and axpy/xpay
      std::vector<uint16_t> colindx; // [nnzbX]
      memWindow_t colindxwin; // column indices in GPU memory

      // for the matrix-matrix subtraction Y -= B:
      std::vector<uint32_t> subset; // [nnzbB], list of inzbX-indices where B is also non-zero
      memWindow_t matBwin; // data of the right hand side operator B
      memWindow_t subsetwin; // subset indices in GPU memory

      uint32_t nRows, nCols; // number of block rows and block columns, max 65,536 columns
      std::vector<uint32_t> rowstart; // [nRows + 1] does not need to be transfered to the GPU

      // for memory management:
      size_t gpu_mem; // device memory requirement in Byte

      // stats:
      double residuum_reached;
      double flops_performed;
      double flops_performed_all;
      int iterations_needed;

      // memory positions
      memWindow_t matXwin; // solution vector in GPU memory
      memWindow_t vec3win; // random number vector in GPU memory
      
      
      // new members to define the action
      std::vector<int64_t> global_target_indices; // [nRows]
      std::vector<int64_t> global_source_indices; // [nCols]
      double r_truncation = 9e18; // radius beyond which the Green function is truncated
      double r_Vconfinement = 9e18; // radius beyond which the confinement potential is added
      double Vconfinement = 0; // potential_prefactor in Hartree

      green_kinetic::finite_difference_plan_t fd_plan[3];
      uint32_t natom_images = 0;
      uint32_t* ApcStart = nullptr; // [natom_images + 1]
      uint32_t* RowStart = nullptr; // [nRows + 1] Needs to be transfered to the GPU?
      uint32_t* rowindx  = nullptr; // [nnzbX] // allows different parallelization strategies
      int16_t (*source_coords)[4] = nullptr; // [nCols][4] internal coordinates
      int16_t (*target_coords)[4] = nullptr; // [nRows][4] internal coordinates
      double  (*Veff)[64] = nullptr; // effective potential
      uint32_t* veff_index = nullptr; // [nRows] indirection list
      uint32_t natoms = 0;
      double (**atom_mat)[2] = nullptr; // [natoms][nc*nc][2] atomic matrices
      atom_t* atom_data = nullptr; // [natom_images]
      double *grid_spacing = nullptr;

      plan_t() {
          std::printf("# construct %s\n", __func__); std::fflush(stdout);
      } // constructor

      ~plan_t() {
          std::printf("# destruct %s\n", __func__); std::fflush(stdout);
          free_memory(ApcStart);
          free_memory(RowStart);
          free_memory(rowindx);
          free_memory(source_coords);
          free_memory(target_coords);
          free_memory(Veff);
          free_memory(veff_index);
          free_memory(atom_mat);
          free_memory(atom_data);
          for (int d = 0; d < 3; ++d) {
              fd_plan[d].~finite_difference_plan_t();
          } // d
      } // destructor

  }; // plan_t





  template <typename floating_point_t=float, unsigned block_size=64>
  class action_t {
  public:
      typedef floating_point_t real_t;
      static int constexpr LM = block_size;
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
          std::printf("# construct %s\n", __func__); std::fflush(stdout);
          char* buffer{nullptr};
          take_memory(buffer);
      } // constructor

      ~action_t() {
          std::printf("# destruct %s\n", __func__); std::fflush(stdout);
          free_memory(apc);
          free_memory(aac);
      } // destructor

      void take_memory(char* &buffer) {
          assert(p->ApcStart);
          auto const nac = p->ApcStart[p->natom_images];
          auto const n = size_t(nac) * size_t(p->nCols);
          // ToDo: could be using GPU memory taking it from the buffer
          // ToDo: how complicated would it be to have only one set of coefficients and multiply in-place?
          apc = get_memory<real_t[2][LM]>(n);
          aac = get_memory<real_t[2][LM]>(n);
      } // take_memory

      void transfer(char* const buffer, cudaStream_t const streamId=0) {
          // no transfers needed since we are using managed memory
          // but we could fill matB with unit blocks here
          // assert(p->subset.size() == p->nCols);
          // clear_on_gpu<real_t[2][LM][LM]>(matB, p->nCols);
          // for (int icol = 0; icol < p->nCols; ++icol) {
          //     for (int i = 0; i < LM; ++i) {
          //        matB[icol][0][i][i] = real_t(1);
          //     } // i
          // } // icol
      } // transfer

      bool has_preconditioner() const { return false; }

      double multiply( // returns the number of flops performed
            real_t         (*const __restrict y)[2][LM][LM] // result, y[nnzbY][2][LM][LM]
          , real_t   const (*const __restrict x)[2][LM][LM] // input,  x[nnzbX][2][LM][LM]
          , uint16_t const (*const __restrict colIndex) // column indices, uint16_t allows up to 65,536 block columns
          , uint32_t const nnzbY // == nnzbX
          , uint32_t const nCols=1 // should match with p->nCols
          , unsigned const l2nX=0
          , cudaStream_t const streamId=0
          , bool const precondition=false
      ) {

      // Action of the SHO-PAW Hamiltonian H onto a trial Green function G:
      // indicies named according to H_{ij} G_{jk} = HG_{ik}
      //                          or A_{ij} X_{jk} =  Y_{ik}
      //
      //      projection:
      //      apc.set(0); // clear
      //      for (jRow < nRows) // reduction
      //        for (RowStart[jRow] <= jnz < Rowstart[jRow + 1]) // reduction
      //          for (jai < nai) // parallel (with data-reuse on Green function element)
      //             for (jprj < nSHO(jai))
      //                for (k64 < 64) // vector parallel
      //                  for (ri < 2) // vector parallel
      //                    for (j64 < 64) // reduction
      //                    {
      //                        apc[apc_start[jai] + jprj][ColIndex[jnz]][ri][k64] +=
      //                        sho_projector[jprj][j64] *
      //                        G[jnz][ri][j64][k64];
      //                    }
          for (int iai = 0; iai < p->natom_images; ++iai) {
              // project at position p->atom_data[iai].pos
              // into apc[(p->ApcStart[iai]:p->ApcStart[iai + 1])*p->nCols][2][64]
          } // iai

      //
      //      meanwhile: kinetic energy including the local potential
      //      HG = finite_difference_kinetic_energy(G);
      //      for (d < 3) // serial
      //        for (istick < nstick[d]) // parallel
      //          while(jnz >= 0, jnz from index_list[d][istick])
      //            for (j4 < 4) // serial (data-reuse)
      //              for (-8 <= imj <= 8) // reduction
      //                for (j16 < 16) // parallel
      //                for (ri < 2) // vector parallel
      //                  for (k64 < 64) // vector parallel
      //                  {
      //                      HG[inz][ri][i4,j16][k64] +=
      //                       G[jnz][ri][j4,j16][k64] * T_FD[|imj|]
      //                  }

          // clear y
          for (size_t inz = 0; inz < nnzbY; ++inz) {
              for (int cij = 0; cij < 2*LM*LM; ++cij) {
                  y[inz][0][0][cij] = 0;
              } // cij
          } // inz

          // kinetic energy
          for (int dd = 0; dd < 3; ++dd) { // derivative direction
//            simple_stats::Stats<> ts;
              auto const & fd = p->fd_plan[dd];
              for (uint32_t il = 0; il < fd.size(); ++il) {
//                SimpleTimer list_timer(__FILE__, __LINE__, "FD-list", 0);
                  //
                  // Warning: this implementation of the finite-difference stencil
                  // is only correct for LM==1, where it performs the lowest order
                  // kinetic energy stencil -0.5*[1 -2 1], but it features the same 
                  // memory traffic profile as an 8th order FD-stencil for LM=64=4*4*4.
                  //
                  auto const *const list = fd.list(il); // we will load at least 2 indices from list
                  int ilist{0}; // counter for index list
                  int inext = list[ilist++]; // load next index of x

                  real_t block[4][2][LM][LM]; // 4 temporary blocks

                  set(block[1][0][0], 2*LM*LM, real_t(0)); // initialize zero, a non-existing block

                  // initially load one block in advance
                  assert(inext > -1 && "the 1st index must be valid!");
                  set(block[2][0][0], 2*LM*LM, x[inext][0][0]); // initial load of an existing block

                  // main loop
                  while (inext > -1) {
                      auto const ihere = inext; // central index of x == index of y
                      inext = list[ilist++]; // get next index
                      // now ilist == 2 in the first iteration
                      if (inext > -1) {
                          set(block[(ilist + 1) & 0x3][0][0], 2*LM*LM, x[inext][0][0]);
                      } else {
                          set(block[(ilist + 1) & 0x3][0][0], 2*LM*LM, real_t(0)); // non-existing block
                      } // load

                      // compute
                      set(        y[ihere][0][0], 2*LM*LM, block[ ilist      % 0x3][0][0]); // coefficient -0.5*c_{0} == 1
                      add_product(y[ihere][0][0], 2*LM*LM, block[(ilist - 1) % 0x3][0][0], real_t(-0.5)); // coefficient -0.5*c_{+/-1} == -0.5
                      add_product(y[ihere][0][0], 2*LM*LM, block[(ilist + 1) % 0x3][0][0], real_t(-0.5)); // coefficient -0.5*c_{+/-1} == -0.5
                      
                  } // while ip > -1
//                ts.add(list_timer.stop());
//                ts.add(ilist - 1); // how many blocks have been loaded and computed?
              } // il
//            std::printf("# derivative in %c-direction: %g +/- %g in [%g, %g] seconds\n", 'x'+dd, ts.avg(), ts.var(), ts.min(), ts.max());
//            std::printf("# derivative in %c-direction: %g +/- %g in [%g, %g] indices\n", 'x'+dd, ts.avg(), ts.var(), ts.min(), ts.max());
              // ===== synchronize to avoid race conditions on y ========
          } // dd

      //
      //      when projection is done: multiply atomic matrices to projection coefficients
      //      aac.set(0); // clear
      //      for (iai < nai) // parallel
      //        for (iRHS < nRHSs) // parallel
      //          for (iprj < nSHO(iai)) // parallel
      //            for (jprj < nSHO(iai)) // reduction
      //              for (kprj < nSHO(iai)) // parallel
      //                for (k64 < 64) // vector parallel
      //                {
      //                    aac[apc_start[iai] + iprj][iRHS][complex][k64] +=
      //                    atom_matrix[ia[iai]][iprj][jprj][complex] *
      //                    apc[apc_start[iai] + jprj][iRHS][complex][k64];
      //                 }
      //                 // with indirection list ia[iai]
      //
      //      after both, addition:
      //      for (iRow < nRows) // parallel
      //        for (RowStart[iRow] <= inz < Rowstart[iRow + 1]) // parallel
      //          for (iai < nai) // reduction
      //            for (iprj < nSHO(iai)) // reduction
      //              for (i64 < 64)
      //                for (ri < 2) // vector parallel
      //                  for (k64 < 64) // vector parallel
      //                  {
      //                      HG[inz][ri][i64][k64] += 
      //                      aac[apc_start[iai] + iprj][ColIndex[inz]][ri][k64] *
      //                      sho_projector[iprj][i64];
      //                  }
      //


      //
      //    furthermore, we have to multiply the local potential Veff[veff_index[iRow]],
      //    the confinement potential and apply the truncation mask
      //    which can be computed on the fly using source_coords[icol][0:2] and target_coords[iRow][0:2]
      //    here, we need 3 numbers: inner_radius^2, truncation_radius^2, potential_prefactor and potential_power
      //    then, if we determined d2 as the distance between any point in the source block from any other
      //    point in the target block, we add the confinement potential
      //        (d2 > inner_radius2)*potential_prefactor*(d2 - inner_radius2)^potential_power
      //    where potential_power should be a template argument and then apply the mask
      //        HG[...] *= (d2 < truncation_radius2)
      //    to make the truncation sphere perfectly round. inner_radius2 < truncation_radius2 assumed.
          float const truncation_radius2  = pow2(p->r_truncation);
          float const inner_radius2       = pow2(p->r_Vconfinement);
          float const potential_prefactor = p->Vconfinement;

          float const hg[3] = {float(p->grid_spacing[0]), float(p->grid_spacing[1]), float(p->grid_spacing[2])};
          
          simple_stats::Stats<> stats_inner, stats_outer, stats_conf, stats_Vconf, stats_d2; 

          for (uint32_t iRow = 0; iRow < p->nRows; ++iRow) { // parallel
              real_t V[LM]; // buffer, probably best using shared memory of the SMx
              set(V, LM, p->Veff[p->veff_index[iRow]]); // load potential values through indirection list
              auto const *const target_coords = p->target_coords[iRow];
              
              for (auto inz = p->RowStart[iRow]; inz < p->RowStart[iRow + 1]; ++inz) { // parallel
                  // apply the local effective potential to all elements in this row
                  for (int i = 0; i < LM; ++i) {
                      add_product(y[inz][0][i], LM, x[inz][0][i], V); // real
                      add_product(y[inz][1][i], LM, x[inz][1][i], V); // imag
                  } // i

                  // apply i,j-dependent potentials like the confinement and truncation mask

                  int const iCol = p->colindx[inz];
                  auto const *const source_coords = p->source_coords[iCol];
                  int const block_coord_diff[3] = {
                      source_coords[0] - target_coords[0],
                      source_coords[1] - target_coords[1],
                      source_coords[2] - target_coords[2]};
                  int constexpr n4 = (64 == LM)? 4 : ((8 == LM) ? 2 : 0);
                  for (int i4z = 0; i4z < n4; ++i4z) {
                  for (int i4y = 0; i4y < n4; ++i4y) {
                  for (int i4x = 0; i4x < n4; ++i4x) {
                      int const i64 = (i4z*n4 + i4y)*n4 + i4x;
                      float vec[3];
                      for (int j4z = 0; j4z < n4; ++j4z) { vec[2] = (block_coord_diff[2]*n4 + i4z - j4z)*hg[2];
                      for (int j4y = 0; j4y < n4; ++j4y) { vec[1] = (block_coord_diff[1]*n4 + i4y - j4y)*hg[1];
                      for (int j4x = 0; j4x < n4; ++j4x) { vec[0] = (block_coord_diff[0]*n4 + i4x - j4x)*hg[0];
                          int const j64 = (j4z*n4 + j4y)*n4 + j4x;
                          float const d2 = pow2(vec[0]) + pow2(vec[1]) + pow2(vec[2]);

                          stats_d2.add(d2);
//                           if (d2 > 71.18745) std::printf("# d2= %g iCol=%i (%i %i %i)*4 + (%i %i %i)  iRow=%i (%i %i %i)*4 + (%i %i %i)\n", d2,
//                               iCol, source_coords[0], source_coords[1], source_coords[2], i4x, i4y, i4z,
//                               iRow, target_coords[0], target_coords[1], target_coords[2], j4x, j4y, j4z);

                          // confinement potential
                          if (d2 >= inner_radius2) {
                              // truncation mask (hard)
                              if (d2 >= truncation_radius2) { // mask
                                  for (int c = 0; c < 2; ++c) { // real and imag
                                      y[inz][c][i64][j64] = 0;
                                  } // c
                                  stats_outer.add(d2);
                              } else {
                                  real_t const V_confinement = potential_prefactor * pow4(d2 - inner_radius2);
                                  for (int c = 0; c < 2; ++c) { // real and imag
                                      y[inz][c][i64][j64] += x[inz][c][i64][j64] * V_confinement;
                                  } // c
                                  stats_Vconf.add(V_confinement);
                                  stats_conf.add(d2);
                              }
                          } else { // inner_radius
                              stats_inner.add(d2);
                          }
                          // we can also define a soft mask: 
                          //    for d2 < truncation_radius2: 
                          //        y *= pow4(d2/truncation_radius2 - 1)
                          //    else
                          //        y = 0

                      }}} // j4xyz
                  }}} // i4xyz

              } // inz
          } // iRow

          { // scope: display stats
              std::printf("# stats V_conf %g +/- %g %s\n", stats_Vconf.avg()*eV, stats_Vconf.var()*eV, _eV);
              // how many grid points do we expect?
              double const f = 4*constants::pi/(3.*hg[0]*hg[1]*hg[2]) * p->nCols*LM;
              double const Vi = pow3(p->r_Vconfinement)*f;
              double const Vo = pow3(p->r_truncation)*f;
              std::printf("# expect inner %g conf %g grid points\n", Vi, Vo - Vi);
              std::printf("# stats  inner %g conf %g outer %g grid points\n", 
                      stats_inner.num(), stats_Vconf.num(), stats_outer.num());

              std::printf("# stats       distance^2 %g [%g, %g] Bohr^2\n", stats_d2.avg(),    stats_d2.min(),    stats_d2.max());
              std::printf("# stats inner distance^2 %g [%g, %g] Bohr^2\n", stats_inner.avg(), stats_inner.min(), stats_inner.max());
              std::printf("# stats conf  distance^2 %g [%g, %g] Bohr^2\n", stats_conf.avg(),  stats_conf.min(),  stats_conf.max());
              std::printf("# stats outer distance^2 %g [%g, %g] Bohr^2\n", stats_outer.avg(), stats_outer.min(), stats_outer.max());
          } // scope

          return 0; // no flops performed so far
      } // multiply

      plan_t* get_plan() { return p; }

    private: // members

      plan_t* p; // the plan is independent of real_t

      // temporary device memory needed for non-local operations 
      // (we could live with a single copy if the application of the atom-centered matrices is in-place)
      real_t (*apc)[2][LM]; // atomic projection coefficients apc[n_all_projection_coefficients][nCols][2][64]
      real_t (*aac)[2][LM]; // atomic addition   coefficients aac[n_all_projection_coefficients][nCols][2][64]

  }; // class action_t

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_Green_action(int const echo=0) {
      return STATUS_TEST_NOT_INCLUDED;
  } // test_Green_action

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_Green_action(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_action
