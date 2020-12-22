#pragma once

#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cmath> // std::sqrt
#include <algorithm> // std::max
#include <utility> // std::swap //, std::move
#include <vector> // std::vector<T>
#include <cstdio> // printf

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "recorded_warnings.hxx" // warn
#include "simple_math.hxx" // ::random
#include "print_tools.hxx" // printf_vector
#include "data_view.hxx" // view4D<T>, view2D<T>
#include "display_units.h" // Ang, _Ang
#include "inline_math.hxx" // pow2
#include "control.hxx" // ::get
#include "simple_stats.hxx" // ::Stats
#include "sho_tools.hxx" // ::nSHO
#include "global_coordinates.hxx" // ::get
#include "constants.hxx" // ::pi

#include "simple_timer.hxx" // SimpleTimer

#ifdef HAS_TFQMRGPU
  #include "tfqmrgpu_core.hxx" // solve<action_t>
  #include "tfqmrgpu_memWindow.h" // memWindow_t
#else
  #include <utility> // std::pair<T>
  typedef std::pair<size_t,size_t> memWindow_t;
  typedef size_t cudaStream_t;
#endif // HAS_TFQMRGPU

namespace green_function {

  double const GByte = 1e-9; char const *const _GByte = "GByte";

  template <typename T>
  T* get_memory(size_t const size=1) {
#ifdef DEBUG
      printf("# managed memory: %lu x %.3f kByte = \t%.6f %s\n", size, 1e-3*sizeof(T), size*sizeof(T)*GByte, _GByte);
#endif
      T* d{nullptr};
#ifdef HAS_CUDA
      CCheck(cudaMallocManaged(&d, size*sizeof(T)));
#else
      d = new T[size];
#endif     
      return d;
  } // get_memory

  template <typename T>
  void free_memory(T* &d) {
      if (nullptr != d) {
#ifdef HAS_CUDA
          CCheck(cudaFree(d));
#else
          delete[] d;
#endif
      } // d
      d = nullptr;
  } // free_memory



#if 0
  template <typename T>
  class block4 {
  public:
      block4(T const initial_value=T(0)) {
          for(int i = 0; i < 64; ++i) { _data[0][0][i] = initial_value; }
      } // constructor

      // access operators (i2,i1,i0) with i0,i1,i2 in [0, 3]
      T const & operator () (int i2, int i1, int i0) const { return _data[i2 & 0x3][i1 & 0x3][i0 & 0x3]; }
      T       & operator () (int i2, int i1, int i0)       { return _data[i2 & 0x3][i1 & 0x3][i0 & 0x3]; }
      // access operators [i] with i in [0, 63]
      T const & operator [] (int i64) const { return _data[0][0][i64 & 0x3f]; }
      T       & operator [] (int i64)       { return _data[0][0][i64 & 0x3f]; }

  private:
      T _data[4][4][4];
  }; // class
#endif  

  struct TruncationMask4 {
      // this is a struct that tells us if source and target grid point are too distant or not.
      // the data is stored in a 3D bitarray[8][8][8] --> 64 Byte

      TruncationMask4(
            int bx, int by, int bz // block coordinate differences, block size 4x4x4
          , double const h[3] // grid spacings
          , double const rc2 // truncation radius squared
      ) {                    
          for(int i8x8 = 0; i8x8 < 8*8; ++i8x8) {
              _data[0][i8x8] = 0; // clear, ToDo: faster using memset
          } // i8x8
          for(int iz = -3; iz < 4; ++iz) {
              double const z = (bz*4 + iz)*h[2], z2 = z*z - rc2;
              // if z2 > 0, skip the following loop
              for(int iy = -3; iy < 4; ++iy) {
                  double const y = (by*4 + iy)*h[1], y2z2 = y*y + z2;
                  // if y2z2 > 0, skip the following loop
                  int i8{0};
                  for(int ix = -3; ix < 4; ++ix) {
                      double const x = (bx*4 + ix)*h[0], r2 = x*x + y2z2;
                      i8 |= (int(r2 < 0.0) << (ix & 0x7)); // if inside rc^2, switch bit on
                  } // ix
                  _data[iz & 0x7][iy & 0x7] = i8; // store
              } // iy
          } // iz
      } // constructor

      // access operators (ix, iy, iz) with ix,iy,iz in [-3, 3]
      bool const operator () (int const ix, int const iy, int const iz) const
          { return (_data[iz & 0x7][iy & 0x7] >> (ix & 0x7)) & 0x1; }

  private:
      int8_t _data[8][8]; // 8x8x8 bit
  }; // class


  typedef struct {
      double pos[3];
      double sigma;
      int32_t gid;
      int32_t ia; // local atom index
      int16_t shifts[3];
      uint8_t nc; // number of coefficients
      int8_t numax;
  } atom_t;
  

  template <typename uint_t, typename int_t> inline
  size_t index3D(uint_t const n[3], int_t const i[3]) {
      // usual 3D indexing
      return size_t(i[2]*n[1] + i[1])*n[0] + i[0];
  } // index3D


  class finite_difference_plan_t {
  private:
      uint32_t *prefix; // in managed memory
      int32_t *fd_list; // in managed memory
      uint32_t n_lists;

  public:
      finite_difference_plan_t() : prefix(nullptr), fd_list(nullptr), n_lists(0) {}

      finite_difference_plan_t(
            int const dd // direction of derivative
          , uint16_t const num_target_coords[3]
          , uint32_t const RowStart[]
          , uint16_t const ColIndex[]
          , view3D<int32_t> const & iRow_of_coords // (Z,Y,X)
          , std::vector<bool> const sparsity_pattern[]
          , unsigned const nRHSs=1
          , int const echo=0
      ) {
          int constexpr X=0, Y=1, Z=2;
          // prepare the finite-difference sequence lists
          char const direction = 'x' + dd;
          assert(X == dd || Y == dd || Z == dd); 
          int num[3];
          set(num, 3, num_target_coords);
          int const num_dd = num[dd];
          num[dd] = 1; // replace number of target blocks in derivative direction
          if (echo > 0) printf("# FD lists for the %c-direction %d %d %d\n", direction, num[X], num[Y], num[Z]);
          simple_stats::Stats<> length_stats;
          std::vector<std::vector<int32_t>> list;
          size_t const max_lists = nRHSs*size_t(num[Z])*num[Y]*num[X];
          list.resize(max_lists);
          int ilist{0};
          for(int iRHS = 0; iRHS < nRHSs; ++iRHS) {
//                if (echo > 0) printf("# FD list for RHS #%i\n", iRHS);
              auto const & sparsity_RHS = sparsity_pattern[iRHS];
              for(int iz = 0; iz < num[Z]; ++iz) { //  
              for(int iy = 0; iy < num[Y]; ++iy) { //   only 2 of these 3 loops have a range > 1
              for(int ix = 0; ix < num[X]; ++ix) { // 
                  int idx[3] = {ix, iy, iz};
                  for(int id = 0; id < num_dd; ++id) { // loop over direction to derive
                      idx[dd] = id; // replace index in the derivate direction
//                           if (echo > 0) printf("# FD list for RHS #%i test coordinates %i %i %i\n",
//                                                   iRHS, idx[X], idx[Y], idx[Z]);
                      auto const idx3 = index3D(num_target_coords, idx);
                      if (sparsity_RHS[idx3]) {
                          auto const iRow = iRow_of_coords(idx[Z], idx[Y], idx[X]);
                          assert(iRow >= 0);

                          int32_t inz_found{-1};
                          for(auto inz = RowStart[iRow]; inz < RowStart[iRow + 1]; ++inz) {
                              if (ColIndex[inz] == iRHS) {
                                  inz_found = inz; // store where it was found
                                  inz = RowStart[iRow + 1]; // stop search loop
                              } // found
                          } // search
                          assert(inz_found >= 0); // fails at inconsistency between sparsity_pattern and the BSR tables

                          assert(ilist < max_lists);
                          list[ilist].push_back(inz_found);
                      } // sparsity pattern
                  } // id
                  int const list_length = list[ilist].size();
                  if (list_length > 0) {
                      length_stats.add(list_length);
//                           if (echo > 0) printf("# FD list of length %d for the %c-direction %i %i %i\n",
//                                                   list_length, direction, idx[X], idx[Y], idx[Z]);
                      // add end-of-sequence markers
                      list[ilist].push_back(-1);
                      
                      ++ilist; // create a new list index
                  } // list_length > 0
              }}} // ixyz
          } // iRHS
          n_lists = ilist;
          if (echo > 0) printf("# %d FD lists for the %c-direction (%.2f %%), length %.3f +/- %.3f, min %g max %g\n",
                                n_lists, direction, n_lists/(max_lists*.01),
                                length_stats.avg(), length_stats.var(), length_stats.min(), length_stats.max());

          // store in managed memory
          prefix = get_memory<uint32_t>(n_lists + 1); // create in GPU memory
          prefix[0] = 0;
          for(int ilist = 0; ilist < n_lists; ++ilist) {
              int const n = list[ilist].size();
              prefix[ilist + 1] = prefix[ilist] + n;
          } // ilist
          size_t const ntotal = prefix[n_lists];
          if (echo > 0) printf("# FD lists for the %c-direction require %d uint32_t, i.e. %.3f kByte\n",
                                  direction, ntotal, ntotal*sizeof(uint32_t)*1e-3);
          fd_list = get_memory<int32_t>(ntotal);
          for(int ilist = 0; ilist < n_lists; ++ilist) {
              int const n = list[ilist].size();
              set(&fd_list[prefix[ilist]], n, list[ilist].data()); // copy into GPU memory
          } // ilist
 
      } // constructor

      finite_difference_plan_t& operator= (finite_difference_plan_t && rhs) {
          // printf("# finite_difference_plan_t& operator= (finite_difference_plan_t && rhs);\n");
          std::swap(fd_list, rhs.fd_list);
          std::swap(prefix , rhs.prefix);
          std::swap(n_lists, rhs.n_lists);
          return *this;
      } // move assignment
      
      finite_difference_plan_t(finite_difference_plan_t && rhs) = delete;
//       {   printf("# finite_difference_plan_t(finite_difference_plan_t && rhs);\n");
//           *this = std::move(rhs);
//       } // move constructor

      finite_difference_plan_t(finite_difference_plan_t const & rhs) = delete; // copy constructor

      finite_difference_plan_t& operator= (finite_difference_plan_t const & rhs) = delete; // move assignment

      ~finite_difference_plan_t() {
          // printf("# destruct %s, pointers= %p and %p\n", __func__, fd_list, prefix); std::fflush(stdout);
          free_memory(fd_list);
          free_memory(prefix);
      } // destructor

      uint32_t size() const { return n_lists; } // number of lists
      int32_t const * list(uint32_t const i) const { return fd_list + prefix[i]; }
      
  }; // class finite_difference_plan_t



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

      finite_difference_plan_t fd_plan[3];
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
          printf("# construct %s\n", __func__); std::fflush(stdout);
      } // constructor

      ~plan_t() {
          printf("# destruct %s\n", __func__); std::fflush(stdout);
          free_memory(ApcStart);
          free_memory(RowStart);
          free_memory(rowindx);
          free_memory(source_coords);
          free_memory(target_coords);
          free_memory(Veff);
          free_memory(veff_index);
          free_memory(atom_mat);
          free_memory(atom_data);
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
          printf("# construct %s\n", __func__); std::fflush(stdout);
          char* buffer{nullptr};
          take_memory(buffer);
      } // constructor

      ~action_t() {
          printf("# destruct %s\n", __func__); std::fflush(stdout);
          free_memory(apc);
          free_memory(aac);
      } // destructor

      void take_memory(char* &buffer) {
          assert(p->ApcStart);
          auto const nac = p->ApcStart[p->natom_images];
          auto const n = size_t(nac) * size_t(p->nCols);
          // ToDo: could be using GPU memory taking it from the buffer
          apc = get_memory<real_t[2][LM]>(n);
          aac = get_memory<real_t[2][LM]>(n);
      } // take_memory

      void transfer(char* const buffer, cudaStream_t const streamId=0) {
          // no transfers needed since we are using managed memory
          // but we could fill matB with unit blocks here
          // assert(p->subset.size() == p->nCols);
          // clear_on_gpu<real_t[2][LM][LM]>(matB, p->nCols);
          // for(int icol = 0; icol < p->nCols; ++icol) {
          //     for(int i = 0; i < LM; ++i) {
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
      // indicies named after H_{ij} G_{jk} = HG_{ik}
      //
      //      projection:
      //      apc.set(0); // clear
      //      for(jRow < nRows) // reduction
      //        for(RowStart[jRow] <= jnz < Rowstart[jRow + 1]) // reduction
      //          for(jai < nai) // parallel (with data-reuse on Green function element)
      //             for(jprj < nSHO(jai))
      //                for(k64 < 64) // vector parallel
      //                  for(ri < 2) // vector parallel
      //                    for(j64 < 64) // reduction
      //                    {
      //                        apc[apc_start[jai] + jprj][ColIndex[jnz]][ri][k64] +=
      //                        sho_projector[jprj][j64] *
      //                        G[jnz][ri][j64][k64];
      //                    }
          for(int iai = 0; iai < p->natom_images; ++iai) {
              // project at position p->atom_data[iai].pos
              // into apc[(p->ApcStart[iai]:p->ApcStart[iai + 1])*p->nCols][2][64]
          } // iai
          
      //
      //      meanwhile: kinetic energy including the local potential
      //      HG = finite_difference_kinetic_energy(G);
      //      for(d < 3) // serial
      //        for(istick < nstick[d]) // parallel
      //          while(jnz >= 0, jnz from index_list[d][istick])
      //            for(j4 < 4) // serial (data-reuse)
      //              for(-8 <= imj <= 8) // reduction
      //                for(j16 < 16) // parallel
      //                for(ri < 2) // vector parallel
      //                  for(k64 < 64) // vector parallel
      //                  {
      //                      HG[inz][ri][i4,j16][k64] +=
      //                       G[jnz][ri][j4,j16][k64] * T_FD[|imj|]
      //                  }
          
          // clear y
          for(size_t inz = 0; inz < nnzbY; ++inz) {
              for(int cij = 0; cij < 2*LM*LM; ++cij) {
                  y[inz][0][0][cij] = 0;
              } // cij
          } // inz
          
          // kinetic energy
          for(int dd = 0; dd < 3; ++dd) { // derivative direction
//            simple_stats::Stats<> ts;
              auto const & fd = p->fd_plan[dd];
              for(uint32_t il = 0; il < fd.size(); ++il) {
//                SimpleTimer list_timer(__FILE__, __LINE__, "FD-list", 0);
                  //
                  // Warning: this implementation of the finite-difference stencil
                  // is only correct for LM==1, where it performs the lowest order
                  // kinetic energy stencil -0.5*[1 -2 1], but is features the same 
                  // memory traffic profile as an 8th order FD-stencil for LM=64.
                  //
                  auto const *const list = fd.list(il); // we will load at least 2 indices from list
                  int ilist{0}; // counter for index list
                  int inext = list[ilist++]; // load next index of x

                  real_t block[4][2][LM][LM]; // 4 temporary blocks
                  
                  set(block[1][0][0], 2*LM*LM, real_t(0)); // initialize non-existing block

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
//            printf("# derivative in %c-direction: %g +/- %g in [%g, %g] seconds\n", 'x'+dd, ts.avg(), ts.var(), ts.min(), ts.max());
//            printf("# derivative in %c-direction: %g +/- %g in [%g, %g] indices\n", 'x'+dd, ts.avg(), ts.var(), ts.min(), ts.max());
              // ===== synchronize to avoid race conditions on y ========
          } // dd

      //
      //      when projection is done: multiply atomic matrices to projection coefficients
      //      aac.set(0); // clear
      //      for(iai < nai) // parallel
      //        for(iRHS < nRHSs) // parallel
      //          for(iprj < nSHO(iai)) // parallel
      //            for(jprj < nSHO(iai)) // reduction
      //              for(kprj < nSHO(iai)) // parallel
      //                for(k64 < 64) // vector parallel
      //                {
      //                    aac[apc_start[iai] + iprj][iRHS][complex][k64] +=
      //                    atom_matrix[ia[iai]][iprj][jprj][complex] *
      //                    apc[apc_start[iai] + jprj][iRHS][complex][k64];
      //                 }
      //                 // with indirection list ia[iai]
      //
      //      after both, addition:
      //      for(iRow < nRows) // parallel
      //        for(RowStart[iRow] <= inz < Rowstart[iRow + 1]) // parallel
      //          for(iai < nai) // reduction
      //            for(iprj < nSHO(iai)) // reduction
      //              for(i64 < 64)
      //                for(ri < 2) // vector parallel
      //                  for(k64 < 64) // vector parallel
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
          float const truncation_radius2 = pow2(p->r_truncation);
          float const inner_radius2      = pow2(p->r_Vconfinement);
          float const potential_prefactor = p->Vconfinement;

          float const hg[3] = {float(p->grid_spacing[0]), float(p->grid_spacing[1]), float(p->grid_spacing[2])};
          
          simple_stats::Stats<> stats_inner, stats_outer, stats_conf, stats_Vconf, stats_d2; 
          
          for(uint32_t iRow = 0; iRow < p->nRows; ++iRow) { // parallel
              real_t V[LM]; // buffer, probably best using shared memory of the SMx
              set(V, LM, p->Veff[p->veff_index[iRow]]); // load potential values through indirection list
              auto const *const target_coords = p->target_coords[iRow];
              
              for(auto inz = p->RowStart[iRow]; inz < p->RowStart[iRow + 1]; ++inz) { // parallel
                  // apply the local effective potential to all elements in this row
                  for(int i = 0; i < LM; ++i) {
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
                  for(int i4z = 0; i4z < n4; ++i4z) {
                  for(int i4y = 0; i4y < n4; ++i4y) {
                  for(int i4x = 0; i4x < n4; ++i4x) {
                      int const i64 = (i4z*n4 + i4y)*n4 + i4x;
                      float vec[3];
                      for(int j4z = 0; j4z < n4; ++j4z) { vec[2] = (block_coord_diff[2]*n4 + i4z - j4z)*hg[2];
                      for(int j4y = 0; j4y < n4; ++j4y) { vec[1] = (block_coord_diff[1]*n4 + i4y - j4y)*hg[1];
                      for(int j4x = 0; j4x < n4; ++j4x) { vec[0] = (block_coord_diff[0]*n4 + i4x - j4x)*hg[0];
                          int const j64 = (j4z*n4 + j4y)*n4 + j4x;
                          float const d2 = pow2(vec[0]) + pow2(vec[1]) + pow2(vec[2]);

                          stats_d2.add(d2);
//                           if (d2 > 71.18745) printf("# d2= %g iCol=%i (%i %i %i)*4 + (%i %i %i)  iRow=%i (%i %i %i)*4 + (%i %i %i)\n", d2,
//                               iCol, source_coords[0], source_coords[1], source_coords[2], i4x, i4y, i4z,
//                               iRow, target_coords[0], target_coords[1], target_coords[2], j4x, j4y, j4z);

                          // confinement potential
                          if (d2 >= inner_radius2) {
                              // truncation mask (hard)
                              if (d2 >= truncation_radius2) { // mask
                                  for(int c = 0; c < 2; ++c) { // real and imag
                                      y[inz][c][i64][j64] = 0;
                                  } // c
                                  stats_outer.add(d2);
                              } else {
                                  real_t const V_confinement = potential_prefactor * pow4(d2 - inner_radius2);
                                  for(int c = 0; c < 2; ++c) { // real and imag
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

          printf("# stats V_conf %g +/- %g %s\n", stats_Vconf.avg()*eV, stats_Vconf.var()*eV, _eV);
          // how many grid points do we expect?
          double const f = 4*constants::pi/(3.*hg[0]*hg[1]*hg[2]) * p->nCols*LM;
          double const Vi = pow3(p->r_Vconfinement)*f;
          double const Vo = pow3(p->r_truncation)*f;
          printf("# expect inner %g conf %g grid points\n", Vi, Vo - Vi);
          printf("# stats  inner %g conf %g outer %g grid points\n", 
                  stats_inner.num(), stats_Vconf.num(), stats_outer.num());
          
          printf("# stats       distance^2 %g [%g, %g] Bohr^2\n",
                  stats_d2.avg(), stats_d2.min(), stats_d2.max());
          printf("# stats inner distance^2 %g [%g, %g] Bohr^2\n",
                  stats_inner.avg(), stats_inner.min(), stats_inner.max());
          printf("# stats conf  distance^2 %g [%g, %g] Bohr^2\n",
                  stats_conf.avg(), stats_conf.min(), stats_conf.max());
          printf("# stats outer distance^2 %g [%g, %g] Bohr^2\n",
                  stats_outer.avg(), stats_outer.min(), stats_outer.max());

          return 0; // no flops performed so far
      } // multiply

      plan_t* get_plan() { return p; }

    private: // members

      plan_t* p; // the plan is independent of real_t

      // temporary device memory needed for non-local operations
      real_t (*apc)[2][LM]; // atomic projection coefficients apc[n_all_projection_coefficients][nCols][2][64]
      real_t (*aac)[2][LM]; // atomic addition   coefficients aac[n_all_projection_coefficients][nCols][2][64]

  }; // class action_t



  inline status_t construct_Green_function(
        int const ng[3]
      , double const hg[3]
      , std::vector<double> const & Veff // [ng[2]*ng[1]*ng[0]]
      , std::vector<double> const & xyzZinso // [natoms*8]
      , std::vector<std::vector<double>> const & atom_mat // atomic hamiltonian and overlap matrix
      , int const echo=0
      , std::complex<double> const *energy_parameter=nullptr
  ) {
      int constexpr X=0, Y=1, Z=2;
      if (echo > 0) printf("\n#\n# %s(%i, %i, %i)\n#\n\n", __func__, ng[X], ng[Y], ng[Z]);

      std::complex<double> E_param(energy_parameter ? *energy_parameter : 0);

      auto const p = new plan_t(); // create a plan how to apply the SHO-PAW Hamiltonian to a block-sparse truncated Green function

      int32_t n_original_Veff_blocks[3] = {0, 0, 0};
      for(int d = 0; d < 3; ++d) {
          assert((0 == (ng[d] & 0x3)) && "All grid dimensions must be a multiple of 4!");
          n_original_Veff_blocks[d] = (ng[d] >> 2); // divided by 4
      } // d

      if (echo > 3) { printf("# n_original_Veff_blocks "); printf_vector(" %d", n_original_Veff_blocks, 3); }

      assert(Veff.size() == ng[Z]*ng[Y]*ng[X]);

//    view3D<block4<double>> Vtot_blocks(n_original_Veff_blocks[Z], n_original_Veff_blocks[Y], n_original_Veff_blocks[X]);
      size_t const n_all_Veff_blocks = n_original_Veff_blocks[Z]*n_original_Veff_blocks[Y]*n_original_Veff_blocks[X];

//    auto Vtot_gpu = get_memory<double>(n_all_Veff_blocks*64); // does not really work flawlessly
//    view4D<double> Vtot(Vtot_gpu, n_original_Veff_blocks[Y], n_original_Veff_blocks[X], 64); // wrap
//    view4D<double> Vtot(n_original_Veff_blocks[Z], n_original_Veff_blocks[Y], n_original_Veff_blocks[X], 64); // get memory
      p->Veff = get_memory<double[64]>(n_all_Veff_blocks); // in managed memory
      { // scope: reorder Veff into block-structured p->Veff
          for(int ibz = 0; ibz < n_original_Veff_blocks[Z]; ++ibz) {
          for(int iby = 0; iby < n_original_Veff_blocks[Y]; ++iby) {
          for(int ibx = 0; ibx < n_original_Veff_blocks[X]; ++ibx) {
              int const ibxyz[3] = {ibx, iby, ibz};
              auto const Veff_index = index3D(n_original_Veff_blocks, ibxyz);
//            auto const Vtot_xyz = Vtot(ibz,iby,ibx); // get a view1D, i.e. double*
              auto const Vtot_xyz = p->Veff[Veff_index];
              for(int i4z = 0; i4z < 4; ++i4z) {
              for(int i4y = 0; i4y < 4; ++i4y) {
              for(int i4x = 0; i4x < 4; ++i4x) {
                  int const i64 = (i4z*4 + i4y)*4 + i4x;
                  size_t const izyx = ((ibz*4 + i4z)
                              *ng[Y] + (iby*4 + i4y)) 
                              *ng[X] + (ibx*4 + i4x);
                  assert(izyx < ng[Z]*ng[Y]*ng[X]);
                  Vtot_xyz[i64] = Veff[izyx];
              }}} // i4
          }}} // xyz
      } // scope

      // Cartesian cell parameters for the unit cell in which the potential is defined
      double const cell[3] = {ng[X]*hg[X], ng[Y]*hg[Y], ng[Z]*hg[Z]};

      // assume periodic boundary conditions and an infinite host crystal,
      // so there is no need to consider k-points
      
      // determine the largest and smallest indices of target blocks
      // given a max distance r_trunc between source blocks and target blocks

      double const r_block_circumscribing_sphere = 1.5*std::sqrt(pow2(hg[X]) + pow2(hg[Y]) + pow2(hg[Z]));
      if (echo > 0) printf("# circumscribing radius= %g %s\n", r_block_circumscribing_sphere*Ang, _Ang);

      // assume that the source blocks lie compact in space
      int const nRHSs = n_original_Veff_blocks[Z] * n_original_Veff_blocks[Y] * n_original_Veff_blocks[X];
      if (echo > 0) printf("# total number of source blocks is %d\n", nRHSs);
      view2D<uint32_t> global_source_coords(nRHSs, 4, 0); // unsigned since they must be in [0, 2^21) anyway
      p->global_source_indices.resize(nRHSs); // [nRHSs]
      p->nCols = nRHSs;
      double center_of_mass_RHS[3] = {0, 0, 0};
      double center_of_RHSs[3]     = {0, 0, 0};
      int32_t min_source_coords[3] = {0, 0, 0}; // global coordinates
      int32_t max_source_coords[3] = {0, 0, 0}; // global coordinates
      int32_t internal_global_offset[3] = {0, 0, 0};
      double max_distance_from_comass{0};
      double max_distance_from_center{0};
      { // scope: determine min, max, center
          {   int iRHS{0};
              for(int ibz = 0; ibz < n_original_Veff_blocks[Z]; ++ibz) {
              for(int iby = 0; iby < n_original_Veff_blocks[Y]; ++iby) {
              for(int ibx = 0; ibx < n_original_Veff_blocks[X]; ++ibx) {
                  global_source_coords(iRHS,X) = ibx;
                  global_source_coords(iRHS,Y) = iby;
                  global_source_coords(iRHS,Z) = ibz;
                  p->global_source_indices[iRHS] = global_coordinates::get(ibx, iby, ibz);
                  for(int d = 0; d < 3; ++d) {
                      int32_t const rhs_coord = global_source_coords(iRHS,d);
                      center_of_mass_RHS[d] += (rhs_coord*4 + 1.5)*hg[d];
                      min_source_coords[d] = std::min(min_source_coords[d], rhs_coord);
                      max_source_coords[d] = std::max(max_source_coords[d], rhs_coord);
                  } // d
                  ++iRHS;
              }}} // xyz
              assert(nRHSs == iRHS);
          } // iRHS
          
          for(int d = 0; d < 3; ++d) {
              center_of_mass_RHS[d] /= std::max(1, nRHSs);
              auto const middle2 = min_source_coords[d] + max_source_coords[d];
              internal_global_offset[d] = middle2/2;
              center_of_RHSs[d] = ((middle2*0.5)*4 + 1.5)*hg[d];
          } // d
          
          if (echo > 0) printf("# internal and global coordinates differ by %d %d %d\n",
              internal_global_offset[X], internal_global_offset[Y], internal_global_offset[Z]);

          p->source_coords = get_memory<int16_t[4]>(nRHSs); // internal coordindates
          { // scope: compute also the largest distance from the center or center of mass
              double max_d2m{0}, max_d2c{0};
              for(int iRHS = 0; iRHS < nRHSs; ++iRHS) {
                  double d2m{0}, d2c{0};
                  p->source_coords[iRHS][3] = 0; // not used
                  for(int d = 0; d < 3; ++d) {
                      p->source_coords[iRHS][d] = global_source_coords(iRHS,d) - internal_global_offset[d];
                      d2m += pow2((global_source_coords(iRHS,d)*4 + 1.5)*hg[d] - center_of_mass_RHS[d]);
                      d2c += pow2((global_source_coords(iRHS,d)*4 + 1.5)*hg[d] - center_of_RHSs[d]);
                  } // d
                  max_d2c = std::max(max_d2c, d2c);
                  max_d2m = std::max(max_d2m, d2m);
              } // iRHS
              max_distance_from_center = std::sqrt(max_d2c);
              max_distance_from_comass = std::sqrt(max_d2m);
          } // scope

      } // scope
      if (echo > 0) printf("# center of mass of RHS blocks is %g %g %g %s\n",
          center_of_mass_RHS[X]*Ang, center_of_mass_RHS[Y]*Ang, center_of_mass_RHS[Z]*Ang, _Ang);
      if (echo > 0) printf("# center of of RHS blocks is %g %g %g %s\n",
          center_of_RHSs[X]*Ang, center_of_RHSs[Y]*Ang, center_of_RHSs[Z]*Ang, _Ang);
      if (echo > 0) printf("# largest distance of RHS blocks from center of mass is %g, from center is %g %s\n",
                              max_distance_from_comass*Ang, max_distance_from_center*Ang, _Ang);

      
      // truncation radius
      double const r_trunc = control::get("green.function.truncation.radius", 10.);
      if (echo > 0) printf("# green.function.truncation.radius=%g %s\n", r_trunc*Ang, _Ang);
      p->r_truncation = std::max(0., r_trunc);
      p->r_Vconfinement = std::min(std::max(0., r_trunc - 2.0), p->r_truncation);

// example:
//    truncation radius in Cu (fcc)
//    lattice constant = 3.522 Ang
//    volume per atom = alat^3 / 4 == 10.922 Ang^3 / atom
//    1000 atoms inside the truncation cluster
//    truncation sphere volume 10922 Ang^3
//    truncation radius = cbrt(3*V/(4*pi)) = 13.764 Ang = 26 Bohr
//    8000 atoms --> 52 Bohr
//    64000 atoms --> 104 Bohr


      // count the number of green function elements for each target block
      size_t nnz{0}; // number of non-zero BSR entries

      uint16_t num_target_coords[3] = {0, 0, 0};
      int16_t  min_target_coords[3] = {0, 0, 0}; // internal coordinates
      int16_t  max_target_coords[3] = {0, 0, 0}; // internal coordinates
      { // scope: create the truncated Green function block-sparsity pattern
          auto const rtrunc       = std::max(0., r_trunc);
          auto const rtrunc_plus  =              rtrunc + 2*r_block_circumscribing_sphere;
          auto const rtrunc_minus = std::max(0., rtrunc - 2*r_block_circumscribing_sphere);
          if (echo > 0) printf("# truncation radius %g, search within %g %s\n", rtrunc*Ang, rtrunc_plus*Ang, _Ang);
          if (echo > 0) printf("# blocks with center distance below %g %s are fully inside\n", rtrunc_minus*Ang, _Ang);

          int16_t itr[3]; // range [-32768, 32767] should be enough
          for(int d = 0; d < 3; ++d) {
              auto const itrunc = std::floor(rtrunc_plus/(4*hg[d]));
              assert(itrunc < 32768 && "target coordinate type is int16_t!");
              itr[d] = int16_t(itrunc);
              assert(itr[d] >= 0);
              min_target_coords[d] = min_source_coords[d] - internal_global_offset[d] - itr[d];
              max_target_coords[d] = max_source_coords[d] - internal_global_offset[d] + itr[d];
              num_target_coords[d] = max_target_coords[d] + 1 - min_target_coords[d];
          } // d
          auto const product_target_blocks = size_t(num_target_coords[Z])*
                                             size_t(num_target_coords[Y])*
                                             size_t(num_target_coords[X]);
          if (echo > 0) printf("# all targets within (%i, %i, %i) and (%i, %i, %i) --> %d x %d x %d = %.3f k\n",
              min_target_coords[X], min_target_coords[Y], min_target_coords[Z],
              max_target_coords[X], max_target_coords[Y], max_target_coords[Z], 
              num_target_coords[X], num_target_coords[Y], num_target_coords[Z], product_target_blocks*.001);
          std::vector<std::vector<uint16_t>> column_indices(product_target_blocks);

          double const r2trunc       = pow2(rtrunc),
                       r2trunc_plus  = pow2(rtrunc_plus),
                       r2trunc_minus = pow2(rtrunc_minus);

          std::vector<int32_t> tag_diagonal(product_target_blocks, -1);
          assert(nRHSs < 65536 && "the integer type of ColIndex is uint16_t!");
          std::vector<std::vector<bool>> sparsity_pattern(nRHSs);

          for(uint16_t iRHS = 0; iRHS < nRHSs; ++iRHS) {
              sparsity_pattern[iRHS] = std::vector<bool>(product_target_blocks, false);
              auto & sparsity_RHS = sparsity_pattern[iRHS];
              auto const *const source_coords = p->source_coords[iRHS]; // internal coordinates
              simple_stats::Stats<> stats[3];
              uint32_t hist[9] = {0,0,0,0, 0,0,0,0, 0}; // distribution of nci
              simple_stats::Stats<> stats_d2[9];
              for(int16_t bz = -itr[Z]; bz <= itr[Z]; ++bz) {
              for(int16_t by = -itr[Y]; by <= itr[Y]; ++by) {
              for(int16_t bx = -itr[X]; bx <= itr[X]; ++bx) {
                  int16_t const bxyz[3] = {bx, by, bz}; // difference vector
                  int16_t target_coords[3]; // internal target coordinates
                  for(int d = 0; d < 3; ++d) {
                      target_coords[d] = source_coords[d] + bxyz[d];
                      assert(target_coords[d] >= min_target_coords[d]);
                      assert(target_coords[d] <= max_target_coords[d]);
                  } // d
                  double const d2 = pow2(bx*4*hg[X]) + pow2(by*4*hg[Y]) + pow2(bz*4*hg[Z]);
                  uint8_t nci{0}; // number of corners inside
                  if (d2 < r2trunc_plus) { // potentially inside, check all 8 corner cases
//                    if (d2 < r2trunc_minus) { nci = 8; } else // skip the 8-corners test for inner blocks -> some speedup
                      { // scope: 8 corner test
                          // i = i4 - j4 --> i in [-3, 3], test only the 8 combinations of {-3, 3}
                          for(int iz = -3; iz < 4; iz += 6) {
                          for(int iy = -3; iy < 4; iy += 6) {
                          for(int ix = -3; ix < 4; ix += 6) {
                              double const d2c = pow2((bx*4 + ix)*hg[X]) + pow2((by*4 + iy)*hg[Y]) + pow2((bz*4 + iz)*hg[Z]);
                              nci += (d2c < r2trunc);
                          }}} // ixyz
                      } // scope

                      if (d2 < r2trunc_minus) assert(8 == nci); // for these, we could skip the 8-corners test
                  } // d2

                  if (nci > 0) { 
                      int16_t idx[3];
                      for(int d = 0; d < 3; ++d) {
                          idx[d] = target_coords[d] - min_target_coords[d];
                          assert(idx[d] >= 0);
                          assert(idx[d] < num_target_coords[d]);
                          stats[d].add(target_coords[d]);
                      } // d
                      auto const idx3 = index3D(num_target_coords, idx);
                      assert(idx3 < product_target_blocks);
                      column_indices[idx3].push_back(iRHS);
                      sparsity_RHS[idx3] = true;
                      if (0 == bx && 0 == by && 0 == bz) {
                          tag_diagonal[idx3] = iRHS;
                      } // diagonal entry
                  } // inside
                  ++hist[nci];
                  stats_d2[nci].add(d2);

              }}} // xyz
              if (echo > 7) printf("# RHS at %i %i %i reaches from (%g, %g, %g) to (%g, %g, %g)\n",
                      global_source_coords(iRHS,X), global_source_coords(iRHS,Y), global_source_coords(iRHS,Z),
                      stats[X].min(), stats[Y].min(), stats[Z].min(),
                      stats[X].max(), stats[Y].max(), stats[Z].max());
              
              if (echo > 2) {
                  if (0 == iRHS) {
                      auto const partial = hist[1] + hist[2] + hist[3] + hist[4] + hist[5] + hist[6] + hist[7];
                      printf("# RHS has %.3f k inside, %.3f k partial and %.3f k outside (of %.3f k checked blocks)\n",
                                    hist[8]*.001, partial*.001, hist[0]*.001, (hist[0] + partial + hist[8])*.001);
                      for(int nci = 0; nci <= 8; ++nci) {
                          printf("# RHS has%9.3f k cases with %d corners inside, d2 stats: %g +/- %g in [%g, %g] Bohr^2\n",
                              hist[nci]*.001, nci, stats_d2[nci].avg(), stats_d2[nci].var(), stats_d2[nci].min(), stats_d2[nci].max());
                      } // nci
                  } // RHS #0
              } // echo
          } // iRHS

          // create a histogram about the distribution of number of columns per row
          std::vector<uint32_t> hist(nRHSs + 1, 0);
          for(size_t idx3 = 0; idx3 < column_indices.size(); ++idx3) {
              auto const n = column_indices[idx3].size();
              ++hist[n];
          } // idx3

          // eval the histogram
          size_t nall{0};
          for(int n = 0; n <= nRHSs; ++n) {
              nall += hist[n];
              nnz  += hist[n]*n;
          } // n
          if (echo > 5) { printf("# histogram total=%.3f k: ", nall*.001); printf_vector(" %d", hist.data(), nRHSs + 1); }
          assert(nall == product_target_blocks); // sanity check

          p->nRows = nall - hist[0]; // the target block entries with no RHS do not create a row
          if (echo > 0) printf("# total number of Green function blocks is %.3f k, "
                               "average %.1f per source block\n", nnz*.001, nnz/double(nRHSs));
          if (echo > 0) printf("# %.3f k (%.1f %% of %.3f k) target blocks are active\n", 
              p->nRows*.001, p->nRows/(product_target_blocks*.01), product_target_blocks*.001);

          assert(nnz < (uint64_t(1) << 32) && "the integer type or RowStart is uint32_t!");

          // resize BSR tables: (Block-compressed Sparse Row format)
          if (echo > 3) { printf("# memory of Green function is %.6f %s (float, twice for double)\n",
                              nnz*2.*64.*64.*sizeof(float)*GByte, _GByte); std::fflush(stdout); }
          auto & ColIndex = p->colindx;
          ColIndex.resize(nnz);
          p->rowindx = get_memory<uint32_t>(nnz);
          auto & RowStart = p->RowStart;
          RowStart = get_memory<uint32_t>(p->nRows + 1);
          RowStart[0] = 0;
          p->veff_index = get_memory<uint32_t>(p->nRows); // indirection list for the local potential
          p->target_coords = get_memory<int16_t[4]>(p->nRows); // view2D<int16_t>(p->nRows, 4, 0);
          p->global_target_indices.resize(p->nRows);
          p->subset.resize(p->nCols); // we assume columns of the unit operator as RHS

          view3D<int32_t> iRow_of_coords(num_target_coords[Z],
                                         num_target_coords[Y],
                                         num_target_coords[X], -1);

          { // scope: fill BSR tables
              simple_stats::Stats<> st;
              uint32_t iRow{0};
              for(uint16_t z = 0; z < num_target_coords[Z]; ++z) { // serial
              for(uint16_t y = 0; y < num_target_coords[Y]; ++y) { // serial
              for(uint16_t x = 0; x < num_target_coords[X]; ++x) { // serial
                  uint16_t const idx[3] = {x, y, z};
                  auto const idx3 = index3D(num_target_coords, idx);
                  assert(idx3 < product_target_blocks);

                  auto const n = column_indices[idx3].size();
                  if (n > 0) {
                      st.add(n);
                      iRow_of_coords(idx[Z], idx[Y], idx[X]) = iRow;

                      RowStart[iRow + 1] = RowStart[iRow] + n;
                      // copy the column indices
                      set(ColIndex.data() + RowStart[iRow], n, column_indices[idx3].data());
                      set(p->rowindx      + RowStart[iRow], n, iRow);
                      // copy the target block coordinates
                      int32_t global_target_coords[3];
                      for(int d = 0; d < 3; ++d) {
                          p->target_coords[iRow][d] = idx[d] + min_target_coords[d];
                          global_target_coords[d] = p->target_coords[iRow][d] + internal_global_offset[d];
                      } // d
                      p->target_coords[iRow][3] = 0; // not used

                      p->global_target_indices[iRow] = global_coordinates::get(global_target_coords); 
                      // global_target_indices is needed to gather the local potential data from other MPI processes

                      { // scope: determine the diagonal entry (source == target)
                          auto const iCol = tag_diagonal[idx3];
                          if (iCol > -1) {
                              assert(iCol < (1ul << 16)); // number range if uint16_t
                              for(int d = 0; d < 3; ++d) {
                                  // sanity check onto internal coordinates
                                  assert(p->source_coords[iCol][d] == p->target_coords[iRow][d]);
                              } // d
                              auto inz = RowStart[iRow];
                              while(iCol != ColIndex[inz]) { ++inz; }
                              p->subset[iCol] = inz;
                          } // iCol valid
                      } // scope

                      if (1) { // scope: fill indirection table for having the local potential only defined in 1 unit cell and repeated periodically
                          int32_t mod[3];
                          for(int d = 0; d < 3; ++d) {
                              mod[d] = global_target_coords[d] % n_original_Veff_blocks[d];
                              mod[d] += (mod[d] < 0)*n_original_Veff_blocks[d];
                          } // d
                          p->veff_index[iRow] = index3D(n_original_Veff_blocks, mod);
                      } // scope
 
                      // count up the number of active rows
                      ++iRow;
                  } // n > 0
              }}} // idx
              assert(p->nRows == iRow);
              assert(nnz == RowStart[p->nRows]);
              printf("# source blocks per target block: average %.1f +/- %.1f in [%g, %g]\n", st.avg(), st.var(), st.min(), st.max());
          } // scope
          column_indices.clear(); // not needed beyond this point

          if (echo > 1) { // measure the difference in the number of target blocks of each RHS
              std::vector<uint32_t> nt(nRHSs, 0);
              // traverse the BSR structure
              for(uint32_t iRow = 0; iRow < p->nRows; ++iRow) {
                  for(auto inz = RowStart[iRow]; inz < RowStart[iRow + 1]; ++inz) {
                      auto const iCol = ColIndex[inz];
                      ++nt[iCol];
                  } // inz
              } // iRow
              // analyze nt
              simple_stats::Stats<> st;
              for(uint16_t iRHS = 0; iRHS < nRHSs; ++iRHS) {
                  st.add(nt[iRHS]);
              } // iRHS
              printf("# target blocks per source block: average %.1f +/- %.1f in [%g, %g]\n", st.avg(), st.var(), st.min(), st.max());
          } // echo

          // Green function is stored sparse 
          // as std::complex<real_t> green[nnz][64][64] 
          // or real_t green[nnz][2][64][64] for the GPU;

          for(int dd = 0; dd < 3; ++dd) { // derivate direction
              // create lists for the finite-difference derivatives
              p->fd_plan[dd] = finite_difference_plan_t(dd
                , num_target_coords
                , RowStart, ColIndex.data()
                , iRow_of_coords
                , sparsity_pattern.data()
                , nRHSs, echo);
          } // dd

          
          // transfer grid spacing into GPU memory
          p->grid_spacing = get_memory<double>(4);
          set(p->grid_spacing, 3, hg);

      } // scope







      int const natoms = atom_mat.size();
      if (echo > 2) printf("\n#\n# %s: Start atom part, %d atoms\n#\n", __func__, natoms);

      // compute which atoms will contribute, the list of natoms atoms may contain a subset of all atoms
      double max_projection_radius{0};
      for(int ia = 0; ia < natoms; ++ia) {
          auto const sigma = xyzZinso[ia*8 + 6];
          auto const projection_radius = std::max(0.0, 6*sigma);
          max_projection_radius = std::max(max_projection_radius, projection_radius);
      } // ia
      if (echo > 3) printf("# largest projection radius is %g %s\n", max_projection_radius*Ang, _Ang);


      p->ApcStart = nullptr;
      p->natom_images = 0;
      { // scope:
          SimpleTimer atom_list_timer(__FILE__, __LINE__, "Atom part");
          
          auto const radius = r_trunc + max_distance_from_center 
                            + 2*max_projection_radius + 2*r_block_circumscribing_sphere;
          int iimage[3];
          size_t nimages{1};
          for(int d = 0; d < 3; ++d) { // parallel
              iimage[d] = int(std::ceil(radius/cell[d]));
              nimages *= (iimage[d]*2 + 1);
          } // d
          auto const natom_images = natoms*nimages;
          if (echo > 3) printf("# replicate %d %d %d atom images, %.3f k total\n", 
                                  iimage[X], iimage[Y], iimage[Z], nimages*.001);

          std::vector<uint32_t> ApcStart(natom_images + 1, 0); // probably larger than needed, should be [nai + 1] later
          std::vector<atom_t> atom_data(natom_images);
          std::vector<uint8_t> atom_ncoeff(natoms, 0); // 0: atom does not contribute

          simple_stats::Stats<> nc_stats;
          uint32_t constexpr COUNT = 0; // 0:count how many blocks are really involved
          double sparse{0}, dense{0}; // stats to assess how much memory can be saved using sparse storage

          size_t iai{0};
          for(int z = -iimage[Z]; z <= iimage[Z]; ++z) { // serial
          for(int y = -iimage[Y]; y <= iimage[Y]; ++y) { // serial
          for(int x = -iimage[X]; x <= iimage[X]; ++x) { // serial
//            if (echo > 3) printf("# periodic shifts  %d %d %d\n", x, y, z);
              int const xyz[3] = {x, y, z};
              for(int ia = 0; ia < natoms; ++ia) { // loop over atoms in the unit cell, serial
                  // suggest a new atomic image position
                  double pos[3];
                  for(int d = 0; d < 3; ++d) { // parallel
                      pos[d] = xyzZinso[ia*8 + d] + xyz[d]*cell[d];
                  } // d
                  auto const atom_id = int32_t(xyzZinso[ia*8 + 4]); 
                  auto const numax =       int(xyzZinso[ia*8 + 5]);
                  auto const sigma =           xyzZinso[ia*8 + 6];
//                   if (echo > 5) printf("# image of atom #%i at %g %g %g %s\n", atom_id, pos[X]*Ang, pos[Y]*Ang, pos[Z]*Ang, _Ang);

                  double const r_projection = pow2(6*sigma); // atom-dependent, precision dependent, assume float here
                  double const r2projection = pow2(r_projection);
//                double const r2projection_plus = pow2(r_projection + r_block_circumscribing_sphere);

                  // check all target blocks if they are inside the projection radius
                  uint32_t ntb{0}; // number of target blocks
                  for(uint32_t iRow = 0; (iRow < p->nRows) && (0 == COUNT*ntb); ++iRow) {
                      auto const *const target_block = p->target_coords[iRow];
                      double d2{0};
                      for(int d = 0; d < 3; ++d) { // serial
                          double const center_of_block = (target_block[d]*4 + 1.5)*hg[d];
                          d2 += pow2(center_of_block - pos[d]);
                      } // d
//                    if (d2 < r2projection_plus) {
                      if (1) {
                          // do more precise checking
//                           if (echo > 9) printf("# target block #%i at %i %i %i gets corner check\n",
//                                           iRow, target_block[X], target_block[Y], target_block[Z]);
                          int nci{0}; // number of corners inside
                          // check 8 corners
                          for(int iz = 0; iz < 4; iz += 3) { // parallel, reduction
                          for(int iy = 0; iy < 4; iy += 3) { // parallel, reduction
                          for(int ix = 0; ix < 4; ix += 3) { // parallel, reduction
                              int const ixyz[3] = {ix, iy, iz};
                              double d2i{0};
                              for(int d = 0; d < 3; ++d) {
                                  double const grid_point = (target_block[d]*4 + ixyz[d])*hg[d];
                                  d2i += pow2(grid_point - pos[d]);
                              } // d
                              if (d2i < r2projection) {
                                  ++nci; // at least one corner of the block 
                                  // ... is inside the projection radius of this atom
                              } // inside the projection radius
                          }}} // ixyz
                          // three different cases: 0, 1...7, 8
                          if (nci > 0) {
                              // atom image contributes
                              ++ntb; // stop outer loop if COUNT==1
                              int const nCols = p->RowStart[iRow + 1] - p->RowStart[iRow];
                              sparse += nCols;
                              dense  += p->nCols; // all columns
                              
//                               if (echo > 7) printf("# target block #%i at %i %i %i is inside\n",
//                                       iRow, target_block[X], target_block[Y], target_block[Z]);
                          } else { // nci
//                               if (echo > 9) printf("# target block #%i at %i %i %i is outside\n",
//                                       iRow, target_block[X], target_block[Y], target_block[Z]);
                          } // nci

                      } else { // d2
//                           if (echo > 21) printf("# target block #%i at %i %i %i is far outside\n",
//                                       iRow, target_block[X], target_block[Y], target_block[Z]);
                      } // d2
                  } // iRow

                  if (ntb > 0) {
                      // atom image contributes, mark in the list to have more than 0 coefficients
                      auto const nc = sho_tools::nSHO(numax);
                      atom_ncoeff[ia] = nc;
                      assert(nc == atom_ncoeff[ia]); // conversion successful
                      nc_stats.add(nc);

                      // at least one target block has an intersection with the projection sphere of this atom image
                      auto & atom = atom_data[iai];
                      set(atom.pos, 3, pos);
                      atom.sigma = sigma;
                      atom.gid = atom_id;
                      atom.ia = ia;
                      set(atom.shifts, 3, xyz);
                      atom.nc = nc; 
                      atom.numax = numax;

                      
//                       if (echo > 5) printf("# image of atom #%i at %g %g %g %s contributes to %d target blocks\n",
//                                                 atom_id, pos[X]*Ang, pos[Y]*Ang, pos[Z]*Ang, _Ang, ntb);
                      ApcStart[iai + 1] = ApcStart[iai] + atom.nc;
                      ++iai;
                  } else {
//                       if (echo > 15) printf("# image of atom #%i at %g %g %g %s does not contribute\n",
//                                                 atom_id, pos[X]*Ang, pos[Y]*Ang, pos[Z]*Ang, _Ang);
                  } // ntb > 0
              } // ia
          }}} // xyz

          auto const nai = iai; // corrected number of atomic images
          if (echo > 3) printf("# %d of %d (%.2f %%) atom images have an overlap with projection spheres\n",
                                  nai, natom_images, nai/(natom_images*.01));
          auto const napc = ApcStart[nai];

          if (echo > 3 && 0 == COUNT) printf("# sparse %g and dense %g\n", sparse, dense);

          if (echo > 3) printf("# number of coefficients per image %.1f +/- %.1f in [%g, %g]\n",
                                nc_stats.avg(), nc_stats.var(), nc_stats.min(), nc_stats.max());
          
          if (echo > 3) printf("# %.3f k atomic projection coefficients, %.2f per atomic image\n", napc*.001, napc/double(nai));
          // projection coefficients for the non-local PAW operations are stored
          // as std::complex<real_t> apc[napc][nRHSs][64]
          // or real_t apc[napc][nRHSs][2][64] for the GPU
          if (echo > 3) printf("# memory of atomic projection coefficients is %.6f %s (float, twice for double)\n",
                                  napc*nRHSs*2.*64.*sizeof(float)*GByte, _GByte);
          p->natom_images = nai;
          p->ApcStart = get_memory<uint32_t>(nai + 1);
          set(p->ApcStart, nai + 1, ApcStart.data()); // copy into GPU memory

          // get all info for the atomic matrices:
          std::vector<int32_t> local_atom_index(natoms, -1); // translation table
          std::vector<int32_t> global_atom_index(natoms, -1); // translation table
          int iac{0};
          for(int ia = 0; ia < natoms; ++ia) { // serial
              if (atom_ncoeff[ia] > 0) {
                  global_atom_index[iac] = ia;
                  local_atom_index[ia] = iac;
                  ++iac;
              } // atom_contributes
          } // ia
          int const nac = iac;
          global_atom_index.resize(nac);
          
          // now store the atomic positions in GPU memory
          p->atom_data = get_memory<atom_t>(nai);
          for(int iai = 0; iai < nai; ++iai) {
              p->atom_data[iai] = atom_data[iai]; // copy
              // translate index
              int const ia = atom_data[iai].ia;
              int const iac = local_atom_index[ia];
              assert(iac > -1);
              p->atom_data[iai].ia = iac;
          } // copy into GPU memory

          // get memory for the matrices and fill
          p->atom_mat = get_memory<double(*)[2]>(nac);
          for(int iac = 0; iac < nac; ++iac) { // parallel
              int const ia = global_atom_index[iac];
              int const nc = atom_ncoeff[ia];
              assert(nc > 0); // the number of coefficients of contributing atoms must be non-zero
              p->atom_mat[iac] = get_memory<double[2]>(nc*nc);
              // fill this with matrix values
              // use MPI communication to find values in atom owner processes
              auto const hmt = atom_mat[ia].data(), ovl = hmt + nc*nc;
              for(int i = 0; i < nc; ++i) {
                  for(int j = 0; j < nc; ++j) {
                      int const ij = i*nc + j;
                      p->atom_mat[iac][ij][0] = hmt[ij] - E_param.real() * ovl[ij];
                      p->atom_mat[iac][ij][1] =         - E_param.imag() * ovl[ij];
                  } // j
              } // i
          } // iac

      } // scope

      { // scope: try out the operator
          typedef float real_t;
          int constexpr LM = 64;
          action_t<real_t,LM> action(p);
          auto const nnzbX = p->colindx.size();
          auto x = get_memory<real_t[2][LM][LM]>(nnzbX);
          auto y = get_memory<real_t[2][LM][LM]>(nnzbX);
          for(size_t i = 0; i < nnzbX*2*LM*LM; ++i) x[0][0][0][i] = 0; // init x

          // benchmark the action
          for(int iteration = 0; iteration < 1; ++iteration) {
              if (echo > 5) { printf("# iteration #%i\n", iteration); std::fflush(stdout); }
              action.multiply(y, x, p->colindx.data(), nnzbX, p->nCols);
              std::swap(x, y);
          } // iteration

      } // scope
      
      return 0;
  } // construct_Green_function
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_Green_function(int const echo=0) {
      status_t stat(0);
      warn("to test green_function. functionality, use --test xml_reading.");
      stat = STATUS_TEST_NOT_INCLUDED; // ToDo: not implemented
      return stat;
  } // test_Green_function

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_Green_function(echo); // deactivated for now
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_function
