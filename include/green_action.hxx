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
#include "green_sparse.hxx"    // ::sparse_t<,>
#include "green_kinetic.hxx"   // ::multiply, ::finite_difference_plan_t
#include "green_potential.hxx" // ::multiply
#include "green_dyadic.hxx"    // ::multiply
#include "simple_stats.hxx" // ::Stats<>

#ifdef HAS_TFQMRGPU
    #include "tfQMRgpu/tfqmrgpu_core.hxx" // tfqmrgpu::solve<action_t>
    #include "tfQMRgpu/tfqmrgpu_memWindow.h" // memWindow_t
#else  // HAS_TFQMRGPU
    #include <utility> // std::pair<T>
    typedef std::pair<size_t,size_t> memWindow_t;
    #ifdef HAS_NO_CUDA
        typedef size_t cudaStream_t;
    #endif // HAS_NO_CUDA
#endif // HAS_TFQMRGPU

// #define debug_printf(...) { std::printf(__VA_ARGS__); std::fflush(stdout); }
#define debug_printf(...)

namespace green_action {

  typedef struct {
      double pos[3]; // position
      double sigma; // Gaussian spread
      int32_t gid; // global identifier
      int32_t ia; // local atom index
      int16_t shifts[3]; // periodic image shifts
      uint8_t nc; // number of coefficients, uint8_t sufficient up to numax=9 --> 220 coefficients
      int8_t numax; // SHO basis size
  } atom_t; // 48 Byte

  typedef struct {
      double pos[3]; // position
      float  oneoversqrtsigma; // Gaussian spread, 1/sqrt(sigma)
      int8_t shifts[3]; // periodic image shifts in [-127, 127]   OR   uint16_t phase_index;
      int8_t lmax; // SHO basis size
  } atom_image_t; // 32 Byte

  
  struct plan_t {

      // members needed for the usage with tfQMRgpu
    
      // for the inner products and axpy/xpay
      std::vector<uint16_t> colindx; // [nnzb], must be a std::vector since nnzb is derived from colindx.size()
      memWindow_t colindxwin; // column indices in GPU memory

      // for the matrix-matrix subtraction Y -= B:
      std::vector<uint32_t> subset; // [nnzbB], list of inzbX-indices at which B is also non-zero
      memWindow_t matBwin; // data of the right hand side operator B
      memWindow_t subsetwin; // subset indices in GPU memory

      // memory positions
      memWindow_t matXwin; // solution vector in GPU memory
      memWindow_t vec3win; // random number vector in GPU memory

      uint32_t nRows = 0; // number of block rows in the Green function
      uint32_t nCols = 0; // number of block columns, max 65,536 columns
      std::vector<uint32_t> rowstart; // [nRows + 1] does not need to be transfered to the GPU

      // for memory management:
      size_t gpu_mem = 0; // device memory requirement in Byte (can this be tracked?)

      // stats:
      float residuum_reached    = 3e38;
      float flops_performed     = 0.f;
      float flops_performed_all = 0.f;
      int   iterations_needed   = -99;

      // =====================================================================================

      // additional members to define the action
      std::vector<int64_t> global_target_indices; // [nRows]
      std::vector<int64_t> global_source_indices; // [nCols]
      double r_truncation   = 9e18; // radius beyond which the Green function is truncated, in Bohr
      double r_Vconfinement = 9e18; // radius beyond which the confinement potential is added, in Bohr
      double Vconfinement   = 0; // potential prefactor, in Hartree

      green_kinetic::finite_difference_plan_t fd_plan[3];
      uint32_t* RowStart = nullptr; // [nRows + 1] Needs to be transfered to the GPU?
      uint32_t* rowindx  = nullptr; // [nnzb] // allows different parallelization strategies
      int16_t (*source_coords)[3+1] = nullptr; // [nCols][3+1] internal coordinates
      int16_t (*target_coords)[3+1] = nullptr; // [nRows][3+1] internal coordinates
      float   (*CubePos)[3+1]       = nullptr; // [nRows]      internal coordinates in float
      int16_t (*target_minus_source)[3+1] = nullptr; // [nnzb][3+1] coordinate differences
      double  (*Veff)[64]  = nullptr; // effective potential, data layout [nRows*Noco*Noco][64]
      int32_t*  veff_index = nullptr; // [nnzb] indirection list

      uint32_t number_of_contributing_atoms = 0;
      double **AtomMatrices = nullptr; // [number_of_contributing_atoms][2*nc^2] atomic matrices, nc number of SHO coefficients of this atom

      uint32_t natom_images  = 0;
      uint32_t* ApcStart     = nullptr; // [natom_images + 1]
      atom_t* atom_data      = nullptr; // [natom_images]
      double (*AtomPos)[3+1] = nullptr; // [natom_images]
      int8_t* AtomLmax       = nullptr; // [natom_images]

      double *grid_spacing       = nullptr; // [3+1]
      double *grid_spacing_trunc = nullptr; // [3]
      bool noncollinear_spin = false;

      green_sparse::sparse_t<> * sparse_SHOprj = nullptr; // [nCols]
      green_sparse::sparse_t<>   sparse_SHOadd;
//    green_sparse::sparse_t<uint16_t> sparse_Green;

      plan_t() {
          debug_printf("# default constructor for %s\n", __func__);
      } // constructor

      ~plan_t() {
          debug_printf("# destruct %s\n", __func__);
          free_memory(ApcStart);
          free_memory(RowStart);
          free_memory(rowindx);
          free_memory(source_coords);
          free_memory(target_coords);
          free_memory(target_minus_source);
          free_memory(Veff);
          free_memory(veff_index);
          if (AtomMatrices) for (int iac = 0; iac < number_of_contributing_atoms; ++iac) free_memory(AtomMatrices[iac]);
          free_memory(AtomMatrices);
          free_memory(atom_data);
          free_memory(grid_spacing);
          free_memory(grid_spacing_trunc);
          for (int dd = 0; dd < 3; ++dd) { // derivative direction
              fd_plan[dd].~finite_difference_plan_t();
          } // dd
          free_memory(AtomPos);
          free_memory(AtomLmax);
          free_memory(CubePos);
          if (sparse_SHOprj) for (int icol = 0; icol < nCols; ++icol) sparse_SHOprj[icol].~sparse_t<>();
          free_memory(sparse_SHOprj);
      } // destructor

  }; // plan_t





  template <typename floating_point_t=float, int R1C2=2, int Noco=1, int n64=64>
  class action_t { // an action as used in tfQMRgpu
  public:
      typedef floating_point_t real_t;
      static int constexpr LM = Noco*n64;
      //
      // This action is an implicit linear operator onto block-sparse structured data.
      // compatible with the core algorithm of the tfqmrgpu-2.0 library.
      // Blocks are sized [LM][LM].
      // Arithmetic according to complex<real_t> 
      // with real_t either float or double
      //
      action_t(plan_t const *plan) 
        : p(plan), apc(nullptr), aac(nullptr)
      {
          assert(1 == Noco && (1 == R1C2 || 2 == R1C2) || 2 == Noco && 2 == R1C2);
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
          auto const natomcoeffs = p->ApcStart ? p->ApcStart[p->natom_images] : 0;
          auto const n = size_t(natomcoeffs) * p->nCols;
          // ToDo: could be using GPU memory taking it from the buffer
          // ToDo: how complicated would it be to have only one set of coefficients and multiply in-place?
          int const echo = 5;
          apc = get_memory<real_t[R1C2][Noco][LM]>(n, echo, "apc");
//        aac = get_memory<real_t[R1C2][Noco][LM]>(n, echo, "aac"); // currently not used
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

      double toy_multiply( // returns the number of flops performed
            real_t         (*const __restrict y)[R1C2][LM][LM] // result, y[nnzb][2][LM][LM]
          , real_t   const (*const __restrict x)[R1C2][LM][LM] // input,  x[nnzb][2][LM][LM]
          , uint16_t const (*const __restrict colIndex) // column indices, uint16_t allows up to 65,536 block columns
          , uint32_t const nnzb // == nnzb, number of nonzero blocks, typically colIndex.size()
          , uint32_t const nCols=1 // should match with p->nCols, number of block columns, assert(colIndex[:] < nCols)
          , unsigned const l2nX=0  // number of levels needed for binary reduction over nnzb
          , cudaStream_t const streamId=0 // CUDA stream to run on
          , bool const precondition=false
      )
        // Toy CPU implementation of green_kinetic and green_potential
      {

          if (echo > 1) { std::printf("# %s<%s,R1C2=%d,Noco=%d>\n", __func__, real_t_name<real_t>(), R1C2, Noco); std::fflush(stdout); }

      // Action of the SHO-PAW Hamiltonian H onto a trial Green function G:
      // indicies named according to H_{ij} G_{jk} = HG_{ik}
      //                          or A_{ij} X_{jk} =  Y_{ik}
      //
      //      projection:
      //      apc.set(0); // clear, WARNING: dimenision of apc in comments are not up-to-date
      //      for (jRow < nRows) // reduction
      //        for (RowStart[jRow] <= jnz < Rowstart[jRow + 1]) // reduction
      //          for (jai < nai) // parallel (with data-reuse on Green function element)
      //             for (jprj < nSHO(jai))
      //                for (k64 < 64) // vector parallel
      //                  for (ri < R1C2) // vector parallel
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
          for (size_t inz = 0; inz < nnzb; ++inz) {
              for (int cij = 0; cij < R1C2*LM*LM; ++cij) {
                  y[inz][0][0][cij] = 0;
              } // cij
          } // inz

          if (echo > 1) { std::printf("# %s<%s,R1C2=%d,Noco=%d>\n", __func__, real_t_name<real_t>(), R1C2, Noco); std::fflush(stdout); }
          
          // kinetic energy
          auto const max_block_index = p->colindx.size();
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

                  real_t block[4][R1C2][LM][LM]; // 4 temporary blocks

                  set(block[1][0][0], R1C2*LM*LM, real_t(0)); // initialize zero, a non-existing block

                  int ilist{0}; // counter for index list
                  int inext = list[ilist++]; // load next index of x
                  // initially load one block in advance
                  assert(inext > -1 && "the 1st index must be valid!");
                  assert(inext < max_block_index);
                  set(block[2][0][0], R1C2*LM*LM, x[inext][0][0]); // initial load of an existing block

                  // main loop
                  while (inext > -1) {
                      auto const ihere = inext; // central index of x == index of y
                      inext = list[ilist++]; // get next index
                      // now ilist == 2 in the first iteration
                      if (inext > -1) {
                          assert(inext < max_block_index);
                          set(block[(ilist + 1) & 0x3][0][0], R1C2*LM*LM, x[inext][0][0]); // load of an existing block
                      } else {
                          set(block[(ilist + 1) & 0x3][0][0], R1C2*LM*LM, real_t(0)); // non-existing block
                      } // load

                      assert(ihere > -1 && "fatal error!");
                      // compute
                      add_product(y[ihere][0][0], R1C2*LM*LM, block[(ilist    ) % 0x3][0][0], real_t( 1.0)); // coefficient -0.5*c_{0}    == 1
                      add_product(y[ihere][0][0], R1C2*LM*LM, block[(ilist - 1) % 0x3][0][0], real_t(-0.5)); // coefficient -0.5*c_{+/-1} == -0.5
                      add_product(y[ihere][0][0], R1C2*LM*LM, block[(ilist + 1) % 0x3][0][0], real_t(-0.5)); // coefficient -0.5*c_{+/-1} == -0.5

                  } // while inext > -1
//                ts.add(list_timer.stop());
//                ts.add(ilist - 1); // how many blocks have been loaded and computed?
              } // il
//            if (echo > 6) std::printf("# derivative in %c-direction: %g +/- %g in [%g, %g] seconds\n", 'x'+dd, ts.mean(), ts.dev(), ts.min(), ts.max());
              // ===== synchronize to avoid race conditions on y ========
          } // dd
          
          if (echo > 1) { std::printf("# %s<%s,R1C2=%d,Noco=%d>\n", __func__, real_t_name<real_t>(), R1C2, Noco); std::fflush(stdout); }

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
      //                    AtomMatricesrix[ia[iai]][iprj][jprj][complex] *
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

          double const *const hg = p->grid_spacing_trunc; // abbreviation

          simple_stats::Stats<> stats_inner, stats_outer, stats_conf, stats_Vconf, stats_d2; 

          for (uint32_t iRow = 0; iRow < p->nRows; ++iRow) { // parallel
              auto const *const target_coords = p->target_coords[iRow];

              for (auto inz = p->RowStart[iRow]; inz < p->RowStart[iRow + 1]; ++inz) { // parallel
                
                  real_t V[LM]; // buffer, probably best using shared memory of the SMx
                  if (p->veff_index[inz] >= 0) {
                      set(V, LM, p->Veff[p->veff_index[inz]]); // load potential values through indirection list
                  } else { set(V, LM, real_t(0)); }

                  // apply the local effective potential to all elements in this row
                  for (int i = 0; i < LM; ++i) {
                      for (int c = 0; c < R1C2; ++c) add_product(y[inz][c][i], LM, x[inz][c][i], V); // real [and imag]
                  } // i

                  // apply i,j-dependent potentials like the confinement and truncation mask

                  int const iCol = p->colindx[inz];
                  auto const *const source_coords = p->source_coords[iCol];
                  int const block_coord_diff[] = {
                      int(source_coords[0]) - int(target_coords[0]),
                      int(source_coords[1]) - int(target_coords[1]),
                      int(source_coords[2]) - int(target_coords[2])};
                  int constexpr n4 = (LM >= 64) ? 4 : 2;
                  for (int i4z = 0; i4z < n4; ++i4z) {
                  for (int i4y = 0; i4y < n4; ++i4y) {
                  for (int i4x = 0; i4x < n4; ++i4x) {
                      int const i64 = (i4z*n4 + i4y)*n4 + i4x;
                      double vec[3];
                      for (int j4z = 0; j4z < n4; ++j4z) { vec[2] = (block_coord_diff[2]*n4 + i4z - j4z)*hg[2];
                      for (int j4y = 0; j4y < n4; ++j4y) { vec[1] = (block_coord_diff[1]*n4 + i4y - j4y)*hg[1];
                      for (int j4x = 0; j4x < n4; ++j4x) { vec[0] = (block_coord_diff[0]*n4 + i4x - j4x)*hg[0];
                          int const j64 = (j4z*n4 + j4y)*n4 + j4x;
                          // distance^2
                          auto const d2 = pow2(vec[0]) + pow2(vec[1]) + pow2(vec[2]);

                          stats_d2.add(d2);
//                           if (d2 > 71.18745) std::printf("# d2= %g iCol=%i (%i %i %i)*4 + (%i %i %i)  iRow=%i (%i %i %i)*4 + (%i %i %i)\n", d2,
//                               iCol, source_coords[0], source_coords[1], source_coords[2], i4x, i4y, i4z,
//                               iRow, target_coords[0], target_coords[1], target_coords[2], j4x, j4y, j4z);

                          // confinement potential
                          if (d2 >= inner_radius2) {
                              // truncation mask (hard)
                              if (d2 >= truncation_radius2) { // mask
                                  for (int c = 0; c < R1C2; ++c) { // real and imag
                                      y[inz][c][i64][j64] = real_t(0);
                                  } // c
                                  stats_outer.add(d2);
                              } else {
                                  real_t const V_confinement = potential_prefactor * pow4(d2 - inner_radius2); // V_confinement ~ distance^8
                                  for (int c = 0; c < R1C2; ++c) { // real and imag
                                      y[inz][c][i64][j64] += x[inz][c][i64][j64] * V_confinement;
                                  } // c
                                  stats_Vconf.add(V_confinement);
                                  stats_conf.add(d2);
                              }
                          } else { // inner_radius2
                              stats_inner.add(d2);
                          } // inner_radius2
                          // we can also define a soft mask: 
                          //    for d2 < truncation_radius2: 
                          //        y *= pow4(d2/truncation_radius2 - 1)
                          //    else
                          //        y = 0

                      }}} // j4xyz
                  }}} // i4xyz

              } // inz
          } // iRow

          if (echo > 1) {
              std::printf("# stats V_conf %g +/- %g %s\n", stats_Vconf.mean()*eV, stats_Vconf.dev()*eV, _eV);
              // how many grid points do we expect?
              double const f = 4.*constants::pi/(3.*hg[0]*hg[1]*hg[2]) * p->nCols*LM;
              double const Vi = pow3(p->r_Vconfinement)*f,
                           Vo = pow3(p->r_truncation)*f;
              std::printf("# expect inner %g conf %g grid points\n", Vi, Vo - Vi);
              std::printf("# stats  inner %g conf %g outer %g grid points\n", 
                             stats_inner.num(), stats_Vconf.num(), stats_outer.num());
              std::printf("# stats       distance^2 %g [%g, %g] Bohr^2\n", stats_d2.mean(),    stats_d2.min(),    stats_d2.max());
              std::printf("# stats inner distance^2 %g [%g, %g] Bohr^2\n", stats_inner.mean(), stats_inner.min(), stats_inner.max());
              std::printf("# stats conf  distance^2 %g [%g, %g] Bohr^2\n", stats_conf.mean(),  stats_conf.min(),  stats_conf.max());
              std::printf("# stats outer distance^2 %g [%g, %g] Bohr^2\n", stats_outer.mean(), stats_outer.min(), stats_outer.max());
          } // echo

          return 0; // no flops performed so far
      } // toy_multiply

      
      
      
      double multiply( // returns the number of flops performed
            real_t         (*const __restrict y)[R1C2][LM][LM] // result, y[nnzb][2][LM][LM]
          , real_t   const (*const __restrict x)[R1C2][LM][LM] // input,  x[nnzb][2][LM][LM]
          , uint16_t const (*const __restrict colIndex) // column indices [nnzb], warning: must be in device memory
          , uint32_t const nnzb // number of nonzero blocks, typically colIndex.size()
          , uint32_t const nCols=1 // should match with p->nCols, number of block columns, assert(colIndex[:] < nCols)
          , unsigned const l2nX=0  // number of levels needed for binary reduction over nnzb
          , cudaStream_t const streamId=0 // CUDA stream to run on
          , bool const precondition=false
      )
        // GPU implementation of green_potential, green_kinetic and green_dyadic
      {
          // start with the potential, assign y to initial values
          green_potential::multiply<real_t,R1C2,Noco>(y, x, p->Veff,
              p->veff_index, p->target_minus_source, p->grid_spacing_trunc, nnzb);

          // add the kinetic energy expressions
          green_kinetic::multiply<real_t,R1C2,Noco>(y, x, p->fd_plan,
              p->grid_spacing, 4, nnzb);
          
          // add the non-local potential using the dyadic action of project + add
          green_dyadic::multiply<real_t,R1C2,Noco>(y, apc, x, p->AtomPos, p->AtomLmax, p->ApcStart, p->natom_images,
              p->sparse_SHOprj, p->AtomMatrices, p->sparse_SHOadd, p->rowindx, colIndex, p->CubePos, p->grid_spacing, nnzb, p->nCols);

          return 0; // no flops performed so far
      } // multiply

      plan_t const* get_plan() { return p; }

    private: // members

      plan_t const *p; // the plan is independent of real_t

      int echo = 9;

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
          { action_t<float ,1,1> action(&plan); }
          { action_t<float ,2,1> action(&plan); }
          { action_t<float ,2,2> action(&plan); }
          { action_t<double,1,1> action(&plan); }
          { action_t<double,2,1> action(&plan); }
          { action_t<double,2,2> action(&plan); }
      } // destruct plan
      std::printf("# %s sizeof(atom_t) = %ld Byte\n", __func__, sizeof(atom_t));
      std::printf("# %s sizeof(atom_image_t) = %ld Byte\n", __func__, sizeof(atom_image_t));
      std::printf("# %s sizeof(plan_t) = %ld Byte\n", __func__, sizeof(plan_t));
      return 0;
  } // test_construction_and_destruction

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_construction_and_destruction(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_action
