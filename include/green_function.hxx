#pragma once

#include <cstdint>    // int64_t, int32_t, uint32_t, int16_t, uint16_t, int8_t, uint8_t
#include <cassert>    // assert
#include <cmath>      // std::sqrt, ::cbrt
#include <algorithm>  // std::max, ::min
#include <utility>    // std::swap
#include <vector>     // std::vector<T>
#include <cstdio>     // std::printf, ::snprintf
#include <complex>    // std::complex

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "simple_timer.hxx" // SimpleTimer
#include "unit_system.hxx" // eV, _eV, Ang, _Ang
#include "print_tools.hxx" // printf_vector
#include "global_coordinates.hxx" // ::get
#include "green_input.hxx" // ::load_Hamitonian

#include "green_memory.hxx" // get_memory, free_memory, real_t_name
#include "green_sparse.hxx" // ::sparse_t<,>

#ifdef HAS_TFQMRGPU

//  #define DEBUG
//  #define DEBUGGPU
    #ifdef HAS_NO_CUDA
        #include "tfQMRgpu/include/tfqmrgpu_cudaStubs.hxx" // cuda... (dummies)
        #define devPtr const __restrict__
    #else  // HAS_NO_CUDA
        #include <cuda.h>
    #endif // HAS_NO_CUDA
    #include "tfQMRgpu/include/tfqmrgpu_memWindow.h" // memWindow_t
    #include "tfQMRgpu/include/tfqmrgpu_linalg.hxx" // ...
    #include "tfQMRgpu/include/tfqmrgpu_core.hxx" // tfqmrgpu::solve<action_t>

#else  // HAS_TFQMRGPU

    #include <utility> // std::pair<T>
    typedef std::pair<size_t,size_t> memWindow_t;
    #ifdef HAS_NO_CUDA
        typedef size_t cudaStream_t;
    #endif // HAS_NO_CUDA

#endif // HAS_TFQMRGPU

#include "green_action.hxx" // ::plan_t, ::action_t, ::atom_t
#include "green_kinetic.hxx" // ::finite_difference_plan_t, index3D
#include "green_dyadic.hxx" // ::dyadic_plan_t
#include "sho_tools.hxx" // ::nSHO
#include "control.hxx" // ::get
#include "boundary_condition.hxx" // Isolated_Boundary, Periodic_Boundary
//int8_t constexpr Isolated_Boundary = 0, Periodic_Boundary = 1;

/*
 *  ToDo plan:
 *    Implement CPU version of SHOprj and SHOadd
 *    Implement periodic boundary conditions in CPU version [optional, but would speed up the following item]
 *    Ensure that tfQMRgpu can invert CPU-only
 *    Verify that the density of states of the Green function matches that found by eigenvalues (explicit solver)
 *    Implement GPU versions of the Hamiltonian
 * 
 */

 /*
  *  Future plan:
  *   Support also density matrix purification scheme (McWeeney filter: x^2(3-2x)
  *   or with a non-trivial overlap operator S (3xSx - 2xSxSx from Phys. Rev. B 50, 17622)
  *   maybe better with norm-conserving PAW formalism --> S == 1
  */

namespace green_function {

  double const GByte = 1e-9; char const *const _GByte = "GByte";

  template <typename real_t, int R1C2=2, int Noco=1>
  void try_action(green_action::plan_t & p, int const n_iterations=1, int const echo=9) {
      if (n_iterations >= 0) { // scope: try to apply the operator
          if (echo > 1) std::printf("# %s<%s,R1C2=%d,Noco=%d>\n", __func__, real_t_name<real_t>(), R1C2, Noco);
          green_action::action_t<real_t,R1C2,Noco,64> action(&p); // constructor

#ifdef  HAS_TFQMRGPU
          p.echo = echo - 5;
          if (echo > 0) std::printf("\n# call tfqmrgpu::mem_count\n");
          tfqmrgpu::solve(action); // try to compile
          if (echo > 5) std::printf("# tfqmrgpu::solve requires %.6f GByte GPU memory\n", p.gpu_mem*1e-9);
          auto memory_buffer = get_memory<char>(p.gpu_mem, echo, "tfQMRgpu-memoryBuffer");
          int const maxiter = control::get("tfqmrgpu.max.iterations", 99.);
          if (echo > 0) std::printf("\n# call tfqmrgpu::solve\n\n");
          {   SimpleTimer timer(__FILE__, __LINE__, __func__, echo);
              tfqmrgpu::solve(action, memory_buffer, 1e-9, maxiter, 0, true);
          } // timer
          if (echo > 0) std::printf("\n# after tfqmrgpu::solve residuum reached= %.1e iterations needed= %d\n",
                                                             p.residuum_reached,    p.iterations_needed);
          if (echo > 6) std::printf("# after tfqmrgpu::solve flop count is %.6f %s\n", p.flops_performed*1e-9, "Gflop");
          free_memory(memory_buffer);
          return;
#endif // HAS_TFQMRGPU

          uint32_t const nnzbX = p.colindx.size();
          int constexpr LM = Noco*64;
          auto x = get_memory<real_t[R1C2][LM][LM]>(nnzbX, echo, "x");
          auto y = get_memory<real_t[R1C2][LM][LM]>(nnzbX, echo, "y");
          for (size_t i = 0; i < nnzbX*R1C2*LM*LM; ++i) {
              x[0][0][0][i] = 0; // init x
          } // i

          auto const nnzb = p.colindx.size();
          auto colIndex = get_memory<uint16_t>(nnzb, echo, "colIndex");
          set(colIndex, nnzb, p.colindx.data()); // copy into GPU memory

          bool const toy = control::get("green_function.benchmark.toy", 0.);

          // benchmark the action
          for (int iteration = 0; iteration < n_iterations; ++iteration) {
              SimpleTimer timer(__FILE__, __LINE__);
              if (echo > 5) { std::printf("# iteration #%i\n", iteration); std::fflush(stdout); }
              if (toy) {
                  action.toy_multiply(y, x, p.colindx.data(), nnzbX, p.nCols);
              } else {
                  action.multiply(y, x, colIndex, nnzbX, p.nCols);
              }
              std::swap(x, y);
          } // iteration

          free_memory(colIndex);
          free_memory(y);
          free_memory(x);
      } // n_iterations >= 0
  } // try_action


  template <typename number_t>
  std::string vec2str(number_t const vec[3], double const f=1, char const*const sep=" ") {
      // convert a vector of 3 numbers into a string, with scaling f and separator
      char s[64]; std::snprintf(s, 63, "%g%s%g%s%g", vec[0]*f, sep, vec[1]*f, sep, vec[2]*f);
      return std::string(s);
  } // vec2str
  #define str(...) vec2str(__VA_ARGS__).c_str()

  inline status_t construct_dyadic_plan(
        green_dyadic::dyadic_plan_t & p
      , double const cell[3]
      , int8_t const boundary_condition[3]
      , double const grid_spacing[3]
      , std::vector<std::vector<double>> const & AtomMatrices
      , std::vector<double> const & xyzZinso
      , uint32_t const nRowsGreen
      , int const nrhs
      , uint32_t const *const rowStartGreen
      , uint16_t const *const colIndexGreen
      , int16_t const (*target_coords)[3+1]
      , double const r_block_circumscribing_sphere
      , double const max_distance_from_center
      , double const r_trunc
      , std::complex<double> E_param
      , int const Noco
      , int const echo=0 // verbosity
  ) {
          SimpleTimer timer(__FILE__, __LINE__, __func__);
          int constexpr X=0, Y=1, Z=2;

          p.nrhs = nrhs;

          // transfer grid spacing into managed GPU memory
          p.grid_spacing = get_memory<double>(3+1, echo, "grid_spacing");
          set(p.grid_spacing, 3, grid_spacing);
          double const r_proj = control::get("green_function.projection.radius", 6.); // in units of sigma
          p.grid_spacing[3] = r_proj; // radius in units of sigma at which the projectors stop

          int const natoms = AtomMatrices.size();
          if (echo > 2) std::printf("\n#\n# %s for %d atoms\n#\n", __func__, natoms);
          assert(xyzZinso.size() == natoms*8);

          // compute which atoms will contribute, the list of natoms atoms may contain a subset of all atoms
          double max_projection_radius{0};
          for (int ia = 0; ia < natoms; ++ia) { // loop over all atoms (can be parallel with reduction)
              auto const sigma = xyzZinso[ia*8 + 6];
              auto const projection_radius = std::max(0.0, r_proj*sigma);
              max_projection_radius = std::max(max_projection_radius, projection_radius);
          } // ia
          if (echo > 3) std::printf("# largest projection radius is %g %s\n", max_projection_radius*Ang, _Ang);

          auto const radius = r_trunc + max_distance_from_center + 2*max_projection_radius + 2*r_block_circumscribing_sphere;
          int iimage[3]; // number of replications of the unit cell to each side
          size_t nimages{1};
          for (int d = 0; d < 3; ++d) { // parallel
              iimage[d] = (Periodic_Boundary == boundary_condition[d]) * std::ceil(radius/cell[d]);
              nimages *= (iimage[d] + 1 + iimage[d]); // iimage images to the left and right
          } // d
          auto const nAtomImages = natoms*nimages;
          if (echo > 3) std::printf("# replicate %s atom images, %ld images total\n", str(iimage), nimages);

          std::vector<uint32_t> AtomImageStarts(nAtomImages + 1, 0); // probably larger than needed, should call resize(nai + 1) later
          std::vector<green_action::atom_t> atom_data(nAtomImages);
          std::vector<int8_t> atom_numax(natoms, -1); // -1: atom does not contribute

          simple_stats::Stats<> nc_stats;
          double sparse{0}, dense{0}; // stats to assess how much memory can be saved using sparse storage

          std::vector<std::vector<uint32_t>> cubes; // stores the row indices of Green function rows
          cubes.reserve(nAtomImages); // maximum (needs 24 Byte per atom image)

          size_t iai{0}; // counter for atomic images
          for (int z = -iimage[Z]; z <= iimage[Z]; ++z) { // serial
          for (int y = -iimage[Y]; y <= iimage[Y]; ++y) { // serial
          for (int x = -iimage[X]; x <= iimage[X]; ++x) { // serial
//            if (echo > 3) std::printf("# periodic shifts  %d %d %d\n", x, y, z);
              int const xyz_shift[] = {x, y, z};
              for (int ia = 0; ia < natoms; ++ia) { // loop over atoms in the unit cell, serial
                  // suggest a shifted atomic image position
                  double atom_pos[3];
                  for (int d = 0; d < 3; ++d) { // unroll
                      atom_pos[d] = xyzZinso[ia*8 + d] + xyz_shift[d]*cell[d];
                  } // d
                  auto const atom_id = int32_t(xyzZinso[ia*8 + 4]); 
                  auto const numax =       int(xyzZinso[ia*8 + 5]);
                  auto const sigma =           xyzZinso[ia*8 + 6] ;
//                   if (echo > 5) std::printf("# image of atom #%i at %s %s\n", atom_id, str(atom_pos, Ang), _Ang);

                  double const r_projection = r_proj*sigma; // atom-dependent, precision dependent, assume float here
                  double const r2projection = pow2(r_projection);
//                double const r2projection_plus = pow2(r_projection + r_block_circumscribing_sphere);

                  // check all target blocks if they are inside the projection radius
                  uint32_t ntb{0}; // number of target blocks
                  for (uint32_t icube = 0; icube < nRowsGreen; ++icube) { // loop over blocks
                      auto const *const target_block = target_coords[icube];
//                    double d2{0};
//                    for (int d = 0; d < 3; ++d) { // unroll
//                        double const center_of_block = (target_block[d]*4 + 1.5)*grid_spacing[d];
//                        d2 += pow2(center_of_block - atom_pos[d]);
//                    } // d
//                    if (d2 < r2projection_plus) {
                      if (1) {
                          // do more precise checking
//                        if (echo > 9) std::printf("# target block #%i at %s gets corner check\n", icube, str(target_block));
                          int nci{0}; // number of corners inside
                          // check 8 corners
                          for (int iz = 0; iz < 4; iz += 3) { // parallel, reduction on nci
                          for (int iy = 0; iy < 4; iy += 3) { // parallel, reduction on nci
                          for (int ix = 0; ix < 4; ix += 3) { // parallel, reduction on nci
                              int const ixyz[] = {ix, iy, iz};
                              double d2i{0}; // init distance^2 of the grid point from the center
                              for (int d = 0; d < 3; ++d) {
                                  double const grid_point = (target_block[d]*4 + ixyz[d] + 0.5)*grid_spacing[d];
                                  d2i += pow2(grid_point - atom_pos[d]);
                              } // d
                              if (d2i < r2projection) {
                                  ++nci; // at least one corner of the block is inside the projection radius of this atom
                              } // inside the projection radius
                          }}} // ix // iy // iz
                          // three different cases: 0, 1...7, 8
                          if (nci > 0) {
                              // atom image contributes
                              if (0 == ntb) {
                                  // this is the 1st cube to contribute
                                  std::vector<uint32_t> cube_list(0); // create an empty vector
                                  cubes.push_back(cube_list); // enlist the vector
                              } // 0 == ntb

                              cubes[iai].push_back(icube); // enlist
                              assert(cubes[iai][ntb] == icube); // check if enlisting worked

                              ++ntb;
                              assert(cubes[iai].size() == ntb);
                              int const nrhs = rowStartGreen[icube + 1] - rowStartGreen[icube];
                              sparse += nrhs;
                              dense  += p.nrhs; // all columns

//                            if (echo > 7) std::printf("# target block #%i at %s is inside\n", icube, str(target_block));
                          } else { // nci
//                            if (echo > 9) std::printf("# target block #%i at %s is outside\n", icube, str(target_block));
                          } // nci

                      } else { // d2 < r2projection_plus
//                        if (echo > 21) std::printf("# target block #%i at %s is far outside\n", icube, str(target_block));
                      } // d2 < r2projection_plus
                  } // icube

                  if (ntb > 0) {
                      // atom image contributes, mark in the list to have more than 0 coefficients
                      auto const nc = sho_tools::nSHO(numax); // number of coefficients for this atom
                      atom_numax[ia] = numax; // atom image does contribute, so this atom does
                      nc_stats.add(nc);

                      // at least one target block has an intersection with the projection sphere of this atom image
                      auto & atom = atom_data[iai];
                      set(atom.pos, 3, atom_pos);
                      atom.sigma = sigma;
                      atom.gid = atom_id;
                      atom.ia = ia; // local atom index
                      set(atom.shifts, 3, xyz_shift);
                      atom.nc = nc; 
                      atom.numax = numax;

//                    if (echo > 5) std::printf("# image of atom #%i at %s %s contributes to %d target blocks\n", atom_id, str(atom_pos, Ang), _Ang, ntb);
                      AtomImageStarts[iai + 1] = AtomImageStarts[iai] + nc;
                      ++iai;

                  } else {
//                    if (echo > 15) std::printf("# image of atom #%i at %s %s does not contribute\n", atom_id, str(atom_pos, Ang), _Ang);
                  } // ntb > 0
              } // ia
          }}} // x // y // z

          auto const nai = iai; // corrected number of atomic images
          if (echo > 3) std::printf("# %ld of %lu (%.2f %%) atom images have an overlap with projection spheres\n",
                                       nai, nAtomImages, nai/(nAtomImages*.01));
          auto const napc = AtomImageStarts[nai];

          if (echo > 3) std::printf("# sparse %g (%.2f %%) of dense %g\n", sparse, sparse/(std::max(1., dense)*.01), dense);

          if (nc_stats.num() > 0 && echo > 3) {
              std::printf("# number of coefficients per image in [%g, %.1f +/- %.1f, %g]\n",
                           nc_stats.min(), nc_stats.mean(), nc_stats.dev(), nc_stats.max());
          } // echo

          if (echo > 3) std::printf("# %.3f k atomic projection coefficients, %.2f per atomic image\n", napc*1e-3, napc/std::max(1., 1.*nai));
          // projection coefficients for the non-local PAW operations are stored
          // as real_t apc[napc*nrhs][R1C2][Noco][Noco*64] on the GPU
          if (echo > 3) std::printf("# memory of atomic projection coefficients is %.6f %s (float, twice for double)\n",
                                                        napc*nrhs*2.*pow2(Noco)*64.*sizeof(float)*GByte, _GByte);
          p.nAtomImages = nai;
          assert(nai == p.nAtomImages); // verify

          using ::green_sparse::sparse_t;

          // planning for the addition of sparse SHO projectors times dense coefficients operation
          auto const nnzb = rowStartGreen[nRowsGreen]; // number of non-zero blocks in the Green function
          std::vector<std::vector<uint32_t>> SHOadd(nnzb);
          // planning for the contraction of sparse Green function times sparse SHO projectors
          std::vector<std::vector<std::vector<uint32_t>>> SHOprj(nrhs);
          for (uint16_t irhs = 0; irhs < nrhs; ++irhs) {
              SHOprj[irhs].resize(nai);
          } // irhs

          for (uint32_t iai = 0; iai < p.nAtomImages; ++iai) {
              for (uint32_t itb = 0; itb < cubes[iai].size(); ++itb) {
                  auto const iRow = cubes[iai][itb];
                  for (auto inzb = rowStartGreen[iRow]; inzb < rowStartGreen[iRow + 1]; ++inzb) {
                      auto const irhs = colIndexGreen[inzb];
                      assert(irhs < nrhs);
                      SHOprj[irhs][iai].push_back(inzb);
                      SHOadd[inzb].push_back(iai);
                  } // inzb
              } // itb
          } // iai
          assert(cubes.size() == p.nAtomImages);
          cubes.resize(0); // release host memory

          p.sparse_SHOadd = sparse_t<>(SHOadd, false, "sparse_SHOadd", echo);
          // sparse_SHOadd: rows == Green function non-zero elements, cols == atom images
          if (echo > 22) {
              std::printf("# sparse_SHOadd.rowStart(%p)[%d + 1]= ", (void*)p.sparse_SHOadd.rowStart(), p.sparse_SHOadd.nRows());
              printf_vector(" %d", p.sparse_SHOadd.rowStart(), p.sparse_SHOadd.nRows() + 1);
              std::printf("# sparse_SHOadd.colIndex(%p)[%d]= ", (void*)p.sparse_SHOadd.colIndex(), p.sparse_SHOadd.nNonzeros());
              printf_vector(" %d", p.sparse_SHOadd.colIndex(), p.sparse_SHOadd.nNonzeros());
          } else if (echo > 5) {
              std::printf("# sparse_SHOadd has %d rows, %ld columns and %d nonzeros\n", p.sparse_SHOadd.nRows(), nai, p.sparse_SHOadd.nNonzeros());
          } // echo
          SHOadd.resize(0); // release host memory

          p.sparse_SHOprj = get_memory<sparse_t<>>(nrhs, echo, "sparse_SHOprj");
          size_t nops{0};
          for (uint16_t irhs = 0; irhs < nrhs; ++irhs) {
              char name[64]; std::snprintf(name, 63, "sparse_SHOprj[irhs=%i of %d]", irhs, nrhs);
              p.sparse_SHOprj[irhs] = sparse_t<>(SHOprj[irhs], false, name, echo - 2);
              // sparse_SHOprj: rows == atom images, cols == Green function non-zero elements
              nops += p.sparse_SHOprj[irhs].nNonzeros();
          } // irhs
          if (echo > 5) std::printf("# sparse_SHOprj has %d*%ld rows, %d columns and %ld nonzeros\n", nrhs,nai, nnzb, nops);
          SHOprj.resize(0); // release host memory

          p.AtomImageStarts = get_memory<uint32_t>(nai + 1, echo, "AtomImageStarts");
          set(p.AtomImageStarts, nai + 1, AtomImageStarts.data()); // copy into GPU memory

          // get all info for the atomic matrices
          std::vector<int32_t> global_atom_index(natoms); // translation table
          std::vector<int32_t>  local_atom_index(natoms, -1); // translation table
          int iac{0};
          for (int ia = 0; ia < natoms; ++ia) { // serial loop over all atoms
              if (atom_numax[ia] > -1) {
                  global_atom_index[iac] = ia;
                  local_atom_index[ia] = iac;
                  ++iac;
              } // atom contributes
          } // ia
          int const nac = iac; // number of contributing atoms
          global_atom_index.resize(nac);

          // store the atomic image positions in GPU memory
          p.AtomImageLmax  = get_memory<int8_t>(nai, echo, "AtomImageLmax"); 
          p.AtomImagePos   = get_memory<double[3+1]>(nai, echo, "AtomImagePos");
          p.AtomImagePhase = get_memory<double[4]>(nai, echo, "AtomImagePhase");
          p.AtomImageShift = get_memory<int8_t[4]>(nai, echo, "AtomImageShift");
          double const phase0[] = {1, 0, 0, 0};

          std::vector<std::vector<uint32_t>> SHOsum(nac);

          for (size_t iai = 0; iai < nai; ++iai) {
              int const ia = atom_data[iai].ia;
              int const iac = local_atom_index[ia]; assert(0 <= iac); assert(iac < nac);
              set(p.AtomImagePos[iai], 3, atom_data[iai].pos);
              set(p.AtomImageShift[iai], 3, atom_data[iai].shifts); p.AtomImageShift[iai][3] = 0;
              p.AtomImagePos[iai][3] = 1./std::sqrt(atom_data[iai].sigma);
              p.AtomImageLmax[iai] = atom_data[iai].numax; // SHO basis size
              if (echo > 9) std::printf("# atom image #%ld has lmax= %d\n", iai, p.AtomImageLmax[iai]);
              SHOsum[iac].push_back(uint32_t(iai));
              set(p.AtomImagePhase[iai], 4, phase0); // TODO construct correct phases
          } // iai, copy into GPU memory

          p.sparse_SHOsum = sparse_t<>(SHOsum, false, "SHOsum", echo);
          SHOsum.resize(0);


          // get memory for the matrices and fill
          p.AtomMatrices = get_memory<double*>(nac, echo, "AtomMatrices");
          p.AtomLmax     = get_memory<int8_t>(nai, echo, "AtomLmax"); 
          p.AtomStarts   = get_memory<uint32_t>(nac + 1, echo, "AtomStarts");
          p.AtomStarts[0] = 0; // init prefetch sum

          for (int iac = 0; iac < nac; ++iac) { // parallel
              auto const ia = global_atom_index[iac];
              int const nc = sho_tools::nSHO(atom_numax[ia]);
              assert(nc > 0); // the number of coefficients of contributing atoms must be non-zero
              p.AtomStarts[iac + 1] = p.AtomStarts[iac] + nc; // create prefect sum
              p.AtomLmax[iac] = atom_numax[ia];
              char name[64]; std::snprintf(name, 63, "AtomMatrices[iac=%d/ia=%d]", iac, ia);
              p.AtomMatrices[iac] = get_memory<double>(Noco*Noco*2*nc*nc, echo, name);
              set(p.AtomMatrices[iac], Noco*Noco*2*nc*nc, 0.0); // clear
              // fill this with matrix values
              assert(2*nc*nc <= AtomMatrices[iac].size());
              // use MPI communication to find values in atom owner processes
              auto const hmt = AtomMatrices[ia].data();
              auto const ovl = hmt + nc*nc;
              for (int i = 0; i < nc; ++i) {
                  for (int j = 0; j < nc; ++j) {
                      int const ij = i*nc + j;
                      p.AtomMatrices[iac][ij]         = hmt[ij] - E_param.real() * ovl[ij]; // real part
                      p.AtomMatrices[iac][ij + nc*nc] =         - E_param.imag() * ovl[ij]; // imag part
                  } // j
              } // i
              // TODO treat Noco components correctly: component 0 and 1 == V_upup and V_dndn
              //                                       component 2 and 3 == V_x and V_y
          } // iac
          p.nAtoms = nac;
          if (echo > 1) std::printf("# found %d contributing atoms and %ld atom images\n", nac, nai);

          p.update_flop_counts(echo); // prepare to count the number of floating point operations

          return 0;
  } // construct_dyadic_plan


  inline status_t update_phases(
        green_dyadic::dyadic_plan_t & p
      , double const k_point[3]
      , int const Noco
      , int const echo=0 // verbosity
  ) {
      std::complex<double> phase[3][2]; // for each direction: forward and backward phase
      for (int d = 0; d < 3; ++d) {
          double const arg = 2*constants::pi*k_point[d];
          phase[d][0] = std::complex<double>(std::cos(arg), std::sin(arg));
          phase[d][1] = std::conj(phase[d][0]); // should be the inverse if the k_point gets an imaginary part 
      } // d

      assert(p.AtomImagePhase && "AtomImagePhase must already be allocated");
      for (uint32_t iai = 0; iai < p.nAtomImages; ++iai) {
          std::complex<double> ph(1, 0);
          for (int d = 0; d < 3; ++d) {
              auto const shift = p.AtomImageShift[iai][d];
              ph *= (shift >= 0) ? intpow(phase[d][0], shift) : intpow(phase[d][1], -shift);
          } // d
          set(p.AtomImagePhase[iai], 4, 0.0);
          p.AtomImagePhase[iai][0] = ph.real();
          p.AtomImagePhase[iai][1] = ph.imag();
          // TODO Noco has been ignored here
      } // iai

      return 0;
  } // update_phases






  int8_t constexpr Vacuum_Boundary = 2;
  // The vacuum boundary condition is an addition to Isolated_Boundary and Periodic_Boundary from boundary_condition.hxx
  // Vacuum_Boundary means that the Green function extends from its source coordinate up to the truncation radius
  // even if this exceeds the isolated boundary at which potential values stop to be defined.

  // The wrap boundary condition is an addition to Periodic_Boundary from boundary_condition.hxx
  // Wrap_Boundary means that the truncation sphere fits into the cell, so k-points have no effect.
  // Nevertheless, it cannot be treated like Isolated_Boundary since target block coordinates may need to be wrapped.

  char const boundary_condition_name[][16] = {"isolated", "periodic", "vacuum"};
  char const boundary_condition_shortname[][8] = {"iso", "peri", "vacu"};


  inline status_t construct_Green_function(
        green_action::plan_t & p // result, create a plan how to apply the SHO-PAW Hamiltonian to a block-sparse truncated Green function
      , uint32_t const ng[3] // numbers of grid points of the unit cell in with the potential is defined
      , int8_t const boundary_condition[3] // boundary conditions
      , double const hg[3] // grid spacings
      , std::vector<double> const & Veff // [ng[2]*ng[1]*ng[0]]
      , std::vector<double> const & xyzZinso // [natoms*8]
      , std::vector<std::vector<double>> const & AtomMatrices // atomic hamiltonian and overlap matrix, [natoms][2*nsho^2]
      , int const echo=0 // log-level
      , std::complex<double> const *energy_parameter=nullptr // E in G = (H - E*S)^{-1}
      , int const Noco=2
  ) {
      int constexpr X=0, Y=1, Z=2;
      if (echo > 0) std::printf("\n#\n# %s(%s)\n#\n\n", __func__, str(ng, 1, ", "));

      p.E_param = energy_parameter ? *energy_parameter : 0;

      int32_t n_original_Veff_blocks[3] = {0, 0, 0};
      for (int d = 0; d < 3; ++d) {
          n_original_Veff_blocks[d] = (ng[d] >> 2); // divided by 4
          assert(n_original_Veff_blocks[d] > 0);
          assert(ng[d] == 4*n_original_Veff_blocks[d] && "All grid dimensions must be a multiple of 4!");
          assert(n_original_Veff_blocks[d] <= (1u << 21) && "Max grid is 2^21 blocks due to global_coordinates");
      } // d
      if (echo > 3) std::printf("# n_original_Veff_blocks %s\n", str(n_original_Veff_blocks));

      assert(Veff.size() == ng[Z]*size_t(ng[Y])*size_t(ng[X]));

      size_t const n_all_Veff_blocks = n_original_Veff_blocks[Z]*n_original_Veff_blocks[Y]*n_original_Veff_blocks[X];

      // regroup effective potential into blocks of 4x4x4
      assert(1 == Noco || 2 == Noco);
      p.noncollinear_spin = (2 == Noco);
      p.Veff = get_memory<double(*)[64]>(4, echo, "Veff");
      for (int mag = 0; mag < 4; ++mag) p.Veff[mag] = nullptr;
      for (int mag = 0; mag < Noco*Noco; ++mag) {
          p.Veff[mag] = get_memory<double[64]>(n_all_Veff_blocks, echo, "Veff[mag]"); // in managed memory
          for (size_t k = 0; k < n_all_Veff_blocks*64; ++k) {
              p.Veff[mag][0][k] = 0; // clear
          } // k
      } // mag

      { // scope: reorder Veff into block-structured p.Veff

          for (int ibz = 0; ibz < n_original_Veff_blocks[Z]; ++ibz) {
          for (int iby = 0; iby < n_original_Veff_blocks[Y]; ++iby) { // block index loops, parallel
          for (int ibx = 0; ibx < n_original_Veff_blocks[X]; ++ibx) {
              int const ibxyz[] = {ibx, iby, ibz};
              auto const Veff_index = index3D(n_original_Veff_blocks, ibxyz);
              for (int i4z = 0; i4z < 4; ++i4z) {
              for (int i4y = 0; i4y < 4; ++i4y) { // grid point index in [0, 4) loops, parallel
              for (int i4x = 0; i4x < 4; ++i4x) {
                  int const i64 = (i4z*4 + i4y)*4 + i4x;
                  size_t const izyx = (size_t(ibz*4 + i4z)
                              *ng[Y] + size_t(iby*4 + i4y)) 
                              *ng[X] + size_t(ibx*4 + i4x); // global grid point index
                  assert(izyx < size_t(ng[Z])*size_t(ng[Y])*size_t(ng[X]));
                  if (Noco > 1) {
                      p.Veff[3][Veff_index][i64] = 0; // set clear
                      p.Veff[2][Veff_index][i64] = 0; // set clear
                      p.Veff[1][Veff_index][i64] = Veff[izyx];
                  } // non-collinear
                  p.Veff[0][Veff_index][i64] = Veff[izyx]; // copy potential value
              }}} // i4
          }}} // xyz
      } // scope

      // Cartesian cell parameters for the unit cell in which the potential is defined
      double const cell[] = {ng[X]*hg[X], ng[Y]*hg[Y], ng[Z]*hg[Z]};
      double const cell_volume = cell[X]*cell[Y]*cell[Z];
      double const average_grid_spacing = std::cbrt(std::abs(hg[X]*hg[Y]*hg[Z]));
      if (echo > 1) { 
          std::printf("\n# Cell summary:\n");
          for (int d = 0; d < 3; ++d) {
              std::printf("#%8d %c-points,%7d blocks, spacing=%9.6f, cell.%c=%9.3f %s, boundary= %d\n",
                  ng[d], 'x' + d, n_original_Veff_blocks[d], hg[d]*Ang, 'x'+ d, cell[d]*Ang, _Ang, boundary_condition[d]);
          } // d
          std::printf("# ================ ============== ================== ====================== ============\n"
              "#%8.3f M points,%12.3f k, average=%9.6f, volume=%9.1f %s^3\n\n",
              ng[X]*1e-6*ng[Y]*ng[Z], n_all_Veff_blocks*1e-3, average_grid_spacing*Ang, cell_volume*pow3(Ang), _Ang);
      } // echo


      int32_t const source_cube = control::get("green_function.source.cube", 1.);
      int32_t n_source_blocks[3] = {0, 0, 0};
      for (int d = 0; d < 3; ++d) {
          n_source_blocks[d] = n_original_Veff_blocks[d];
          if (source_cube) n_source_blocks[d] = std::min(n_original_Veff_blocks[d], source_cube);
      } // d
      if (source_cube && echo > 0) std::printf("\n# limit n_source_blocks to +green_function.source.cube=%d\n", source_cube);
      if (echo > 3) std::printf("# n_source_blocks %s\n", str(n_source_blocks));

      // we assume that the source blocks lie compact and preferably close to each other      

      // determine the largest and smallest indices of target blocks
      // given a max distance r_trunc between source blocks and target blocks

      // assume that the source blocks lie compact in space
      uint32_t const nrhs = n_source_blocks[Z] * n_source_blocks[Y] * n_source_blocks[X];
      if (echo > 0) std::printf("# total number of source blocks is %d\n", nrhs);

      view2D<int32_t> global_source_coords(nrhs, 4, 0);
      p.global_source_indices.resize(nrhs, -1); // [nrhs]
      p.nCols = nrhs;
      double center_of_mass_RHS[] = {0, 0, 0};
      double center_of_RHSs[]     = {0, 0, 0};
      int32_t min_global_source_coords[] = {1 << 21, 1 << 21, 1 << 21}; // global coordinates
      int32_t max_global_source_coords[] = {0, 0, 0}; // global coordinates
      int32_t global_internal_offset[]   = {0, 0, 0};
      double max_distance_from_comass{0};
      double max_distance_from_center{0};
      { // scope: determine min, max, center
          double const by_nrhs = 1./std::max(1., double(nrhs));
          {   uint32_t irhs{0};
              for (int32_t ibz = 0; ibz < n_source_blocks[Z]; ++ibz) {
              for (int32_t iby = 0; iby < n_source_blocks[Y]; ++iby) { // block index loops, serial
              for (int32_t ibx = 0; ibx < n_source_blocks[X]; ++ibx) {
                  int32_t const ibxyz[] = {ibx + (n_original_Veff_blocks[X] - n_source_blocks[X])/2,
                                           iby + (n_original_Veff_blocks[Y] - n_source_blocks[Y])/2,
                                           ibz + (n_original_Veff_blocks[Z] - n_source_blocks[Z])/2, 0};
                  set(global_source_coords[irhs], 4, ibxyz);
                  p.global_source_indices[irhs] = global_coordinates::get(ibxyz);
                  for (int d = 0; d < 3; ++d) {
                      int32_t const rhs_coord = global_source_coords(irhs,d);
                      center_of_mass_RHS[d] += (rhs_coord*4 + 1.5)*hg[d]*by_nrhs;
                      min_global_source_coords[d] = std::min(min_global_source_coords[d], rhs_coord);
                      max_global_source_coords[d] = std::max(max_global_source_coords[d], rhs_coord);
                  } // d
                  ++irhs;
              }}} // xyz
              assert(nrhs == irhs);
          } // irhs
          if (echo > 0) std::printf("# all sources within (%s) and (%s)\n",
              str(min_global_source_coords), str(max_global_source_coords));

          for (int d = 0; d < 3; ++d) {
              auto const middle2 = min_global_source_coords[d] + max_global_source_coords[d];
              global_internal_offset[d] = middle2/2;
              center_of_RHSs[d] = ((middle2*0.5)*4 + 1.5)*hg[d];
          } // d

          if (echo > 0) std::printf("# internal and global coordinates differ by %s\n", str(global_internal_offset));

          p.source_coords = get_memory<int16_t[4]>(nrhs, echo, "source_coords"); // internal coordinates
          { // scope: compute also the largest distance from the center or center of mass
              double max_d2m{0}, max_d2c{0};
              for (uint32_t irhs = 0; irhs < nrhs; ++irhs) {
                  double d2m{0}, d2c{0};
                  for (int d = 0; d < 3; ++d) {
                      int32_t const source_coord = global_source_coords(irhs,d) - global_internal_offset[d];
                      p.source_coords[irhs][d] = source_coord; assert(source_coord == p.source_coords[irhs][d] && "safe assign");
                      d2m += pow2((global_source_coords(irhs,d)*4 + 1.5)*hg[d] - center_of_mass_RHS[d]);
                      d2c += pow2((global_source_coords(irhs,d)*4 + 1.5)*hg[d] - center_of_RHSs[d]);
                  } // d
                  p.source_coords[irhs][3] = 0; // not used
                  max_d2c = std::max(max_d2c, d2c);
                  max_d2m = std::max(max_d2m, d2m);
              } // irhs
              max_distance_from_center = std::sqrt(max_d2c);
              max_distance_from_comass = std::sqrt(max_d2m);
          } // scope

      } // scope
      if (echo > 0) std::printf("# center of mass of RHS blocks is %s %s\n", str(center_of_mass_RHS, Ang), _Ang);
      if (echo > 0) std::printf("# center of coords  RHS blocks is %s %s\n", str(center_of_RHSs    , Ang), _Ang);
      if (echo > 0) std::printf("# largest distance of RHS blocks from center of mass is %g, from center is %g %s\n",
                                                       max_distance_from_comass*Ang, max_distance_from_center*Ang, _Ang);

      // truncation radius
      double const r_trunc = control::get("green_function.truncation.radius", 10.);
      if (echo > 0) std::printf("# green_function.truncation.radius=%g %s, %.1f grid points\n", r_trunc*Ang, _Ang, r_trunc/average_grid_spacing);
      p.r_truncation  = std::max(0., r_trunc);
      // confinement potential
      p.r_confinement = std::min(std::max(0., r_trunc - 2.0), p.r_truncation);
      p.V_confinement = 1;
      if (echo > 0) std::printf("# confinement potential %g*(r/Bohr - %g)^4 %s\n", p.V_confinement*eV, p.r_confinement, _eV);
      if (echo > 0) std::printf("# V_confinement(r_truncation)= %g %s\n", p.V_confinement*eV*pow4(r_trunc - p.r_confinement), _eV);

// example:
//    truncation radius in Cu (fcc)
//    lattice constant = 3.522 Ang
//    volume per atom = alat^3 / 4 == 10.922 Ang^3 / atom
//    1000 atoms inside the truncation cluster
//    truncation sphere volume 10922 Ang^3
//    truncation radius = cbrt(3*V/(4*pi)) = 13.764 Ang = 26 Bohr
//    8000 atoms --> 52 Bohr
//    64000 atoms --> 104 Bohr

      double const r_block_circumscribing_sphere = 0.5*(4 - 1)*std::sqrt(pow2(hg[X]) + pow2(hg[Y]) + pow2(hg[Z]));
      if (echo > 0) std::printf("# circumscribing radius= %g %s\n", r_block_circumscribing_sphere*Ang, _Ang);


      // count the number of green function elements for each target block

      uint32_t num_target_coords[3] = {0, 0, 0};
      int32_t  min_target_coords[3] = {0, 0, 0}; // global coordinates
      int32_t  max_target_coords[3] = {0, 0, 0}; // global coordinates
      int32_t  min_targets[3]       = {0, 0, 0}; // global coordinates
      { // scope: create the truncated Green function block-sparsity pattern
          auto const rtrunc       = std::max(0., r_trunc);
          auto const rtrunc_plus  =              rtrunc + 2*r_block_circumscribing_sphere;
          auto const rtrunc_minus = std::max(0., rtrunc - 2*r_block_circumscribing_sphere);
          if (echo > 0) std::printf("# truncation radius %g %s, search within %g %s\n", rtrunc*Ang, _Ang, rtrunc_plus*Ang, _Ang);
          if (echo > 0 && rtrunc_minus > 0) std::printf("# blocks with center distance below %g %s are fully inside\n", rtrunc_minus*Ang, _Ang);

          uint8_t is_periodic[] = {0, 0, 0, 0}; // is_periodic > 0 if periodic BC and truncation sphere overlaps with itself
          double h[] = {hg[X], hg[Y], hg[Z]}; // customized grid spacing used in green_potential::multiply
          int32_t itr[3];
          for (int d = 0; d < 3; ++d) {

              if (r_trunc >= 0) {
                  char keyword[64]; std::snprintf(keyword, 63, "green_function.scale.grid.spacing.%c", 'x' + d); // this feature allows also truncation ellipsoids
                  auto const scale_h = control::get(keyword, 1.);
                  if (scale_h >= 0) {
                      h[d] = hg[d]*scale_h;
                      if (echo > 1 && 1 != scale_h) std::printf("# scale grid spacing in %c-direction for truncation from %g to %g %s\n", 'x'+d, hg[d]*Ang, h[d]*Ang, _Ang);
                  } // scale_h >= 0
                  if (Periodic_Boundary == boundary_condition[d]) {
                      // periodic boundary conditions
                      double const deformed_cell = h[d]*ng[d];
                      if (2*rtrunc > deformed_cell) {
                          if (h[d] > 0) {
                              warn("truncation sphere (diameter= %g %s) does not fit cell in %c-direction (%g %s)\n#          "
                                   "better use +%s=0 for cylindrical truncation", 2*rtrunc*Ang, _Ang, 'x' + d, deformed_cell*Ang, _Ang, keyword);
                              is_periodic[d] = std::min(255., std::ceil(rtrunc/deformed_cell)); // truncation sphere may overlap with its periodic images
                          } else { // h[d] > 0
                              is_periodic[d] = 1; // truncation sphere may overlap with its periodic images
                          }
                          if (echo > 1) std::printf("# boundary condition in %c-direction has up to %d images\n", 'x' + d, is_periodic[d]);
                      } else {
                          // truncation sphere does not overlap with its periodic images, we can treat it like Vacuum_Boundary
                          if (echo > 1) std::printf("# periodic boundary condition in %c-direction is wrapped\n", 'x' + d);
                      }
                  } // periodic boundary condition
              } // r_trunc >= 0

              // how many blocks around the source block do we need to check
              itr[d] = (h[d] > 0) ? std::floor(rtrunc_plus/(4*h[d])) : n_original_Veff_blocks[d]/2;
              assert(itr[d] >= 0);

              min_target_coords[d] = min_global_source_coords[d] - itr[d];
              max_target_coords[d] = max_global_source_coords[d] + itr[d];

              if (Isolated_Boundary == boundary_condition[d]) {
                  auto const n = n_original_Veff_blocks[d];
                  // limit to global coordinates in [0, n)
                  min_target_coords[d] = std::max(min_target_coords[d], 0);
                  max_target_coords[d] = std::min(max_target_coords[d], n - 1);
              }
              auto const n_target_coords = int32_t(max_target_coords[d]) + 1 - min_target_coords[d];
              assert(n_target_coords > 0);
              num_target_coords[d] = n_target_coords;

              min_targets[d] = min_target_coords[d];
              if (is_periodic[d]) {
                  // prepare for all blocks in this direction
                  min_targets[d] = 0;
                  num_target_coords[d] = n_original_Veff_blocks[d];
              } // is_periodic[d]

          } // d
          int constexpr ANY = 3;
          is_periodic[ANY] = (is_periodic[X] > 0) + (is_periodic[Y] > 0) + (is_periodic[Z] > 0);
          if (r_trunc < 0) {
              if (echo > 1) std::printf("# truncation deactivated\n");
          } else {
              if (echo > 1) std::printf("# truncation beyond %d %d %d blocks\n", itr[X], itr[Y], itr[Z]);
          }
          auto const product_target_blocks = size_t(num_target_coords[Z])*
                                             size_t(num_target_coords[Y])*
                                             size_t(num_target_coords[X]);
          if (echo > 0) std::printf("# all targets within (%s) and (%s) --> %s = %.3f k\n",
              str(min_target_coords), str(max_target_coords),
              str(num_target_coords, 1, " x "), product_target_blocks*.001);
          assert(product_target_blocks > 0);
          std::vector<std::vector<uint16_t>> column_indices(product_target_blocks);

          double const r2trunc        = pow2(rtrunc),
                       r2trunc_plus   = pow2(rtrunc_plus),
                       r2trunc_minus  = pow2(rtrunc_minus),
                       r2block_circum = pow2(r_block_circumscribing_sphere*3); // Maybe this can be reduced to *2

          std::vector<int32_t> tag_diagonal(product_target_blocks, -1);
          assert(nrhs < (1ul << 16) && "the integer type of ColIndex is uint16_t!");
          std::vector<std::vector<bool>> sparsity_pattern(nrhs); // std::vector<bool> is a memory-saving bit-array

          simple_stats::Stats<> inout[4];
          for (uint16_t irhs = 0; irhs < nrhs; ++irhs) { // serial loop if the order in column_indices is relevant
              auto & sparsity_RHS = sparsity_pattern[irhs]; // abbreviate
              sparsity_RHS.resize(product_target_blocks, false);
              auto const *const source_coords = global_source_coords[irhs]; // global source block coordinates
              simple_stats::Stats<> stats[3];
              int constexpr max_nci = 27;
              std::vector<int> hist(1 + max_nci, 0); // distribution of nci
              std::vector<simple_stats::Stats<>> stats_d2(1 + max_nci);
              size_t hit_single{0}, hit_multiple{0};
              int64_t idx3_diagonal{-1};
              
              int32_t b_first[3], b_last[3];
              for (int d = 0; d < 3; ++d) {
                  b_first[d] = std::max(-itr[d], min_target_coords[d] - source_coords[d]);
                  b_last[d]  = std::min( itr[d], max_target_coords[d] - source_coords[d]);
              } // d
              int32_t target_coords[3]; // global target block coordinates
              for (int32_t bz = b_first[Z]; bz <= b_last[Z]; ++bz) { target_coords[Z] = source_coords[Z] + bz;
                  assert(target_coords[Z] >= min_target_coords[Z] && target_coords[Z] <= max_target_coords[Z]);
              for (int32_t by = b_first[Y]; by <= b_last[Y]; ++by) { target_coords[Y] = source_coords[Y] + by;
                  assert(target_coords[Y] >= min_target_coords[Y] && target_coords[Y] <= max_target_coords[Y]);
              for (int32_t bx = b_first[X]; bx <= b_last[X]; ++bx) { target_coords[X] = source_coords[X] + bx;
                  assert(target_coords[X] >= min_target_coords[X] && target_coords[X] <= max_target_coords[X]);

//                for (int d = 0; d < 3; ++d ) assert(target_coords[d] >= min_target_coords[d]
//                                                 && target_coords[d] <= max_target_coords[d]);

                  // d2 is the distance^2 of the block centers
                  double const d2 = pow2(bx*4*h[X])
                                  + pow2(by*4*h[Y])
                                  + pow2(bz*4*h[Z]);

                  int nci{0}; // init number of corners inside
                  if (d2 < r2trunc_plus) { // potentially inside, check all 8 or 27 corner cases
//                    if (d2 < r2trunc_minus) { nci = max_nci; } else // skip the 8- or 27-corners test for inner blocks -> some speedup
                      { // scope: 8 or 27 corner test
                          int const far = (d2 > r2block_circum);
                          int const mci = far ? 8 : 27;
                          // i = i4 - j4 --> i in [-3, 3],
                          //     if two blocks are far from each other, we test only the 8 combinations of |{-3, 3}|^3
                          //     for blocks close to each other, we test all 27 combinations of |{-3, 0, 3}|^3
                          for (int iz = -3; iz <= 3; iz += 3 + 3*far) { double const d2z = pow2((bz*4 + iz)*h[Z]);
                          for (int iy = -3; iy <= 3; iy += 3 + 3*far) { double const d2y = pow2((by*4 + iy)*h[Y]) + d2z;
                          for (int ix = -3; ix <= 3; ix += 3 + 3*far) { double const d2c = pow2((bx*4 + ix)*h[X]) + d2y;
#if 0
                              if (0 == irhs && (d2c < r2trunc) && echo > 17) {
                                  std::printf("# %s: b= %i %i %i, i-j %i %i %i, d^2= %g %s\n", 
                                      __func__, bx,by,bz, ix,iy,iz, d2c, (d2c < r2trunc)?"in":"out");
                              }
#endif // 0
                              nci += (d2c < r2trunc); // add 1 if inside
                          }}} // ix iy iz
                          if (d2 < r2trunc_minus) assert(mci == nci); // for these, we could skip the 8-corners test
                          nci = (nci*27)/mci; // limit nci to [0, 27]
                      } // scope

                  } // d2 < r2trunc_plus

                  if (nci > 0) {
                      // any grid point in block (bx,by,bz) is closer than rtrunc to any grid point in the source block
                      int32_t idx[3];
                      for (int d = 0; d < 3; ++d) {
                          idx[d] = target_coords[d] - min_targets[d];
                          // TODO in the periodic case we have to make sure that we hit no element twice
                          if (is_periodic[d]) idx[d] = (target_coords[d] + 999*n_original_Veff_blocks[d]) % n_original_Veff_blocks[d];
                          assert(0 <= idx[d]); assert(idx[d] < num_target_coords[d]);
                      } // d
                      auto const idx3 = index3D(num_target_coords, idx); // flat index
                      assert(idx3 < product_target_blocks);
                      if (false == sparsity_RHS[idx3]) {
                          sparsity_RHS[idx3] = true;
                          column_indices[idx3].push_back(irhs);
                          for (int d = 0; d < 3; ++d) stats[d].add(target_coords[d]);
                          ++hit_single;
                      } else {
                          ++hit_multiple;
                      } // has not been hit yet
                      if (0 == bx && 0 == by && 0 == bz) {
                          assert(-1 == idx3_diagonal);
                          assert(-1 == tag_diagonal[idx3]);
                          tag_diagonal[idx3] = irhs;
                          idx3_diagonal = idx3;
                      } // diagonal entry
                  } // nci > 0
                  ++hist[nci];
                  stats_d2[nci].add(d2);

              }}} // xyz
              if (echo > 8) std::printf("# RHS#%i has %ld single and %ld multiple hits\n", irhs, hit_single, hit_multiple);
              assert((0 == hit_multiple) || (hit_multiple > 0 && is_periodic[ANY]));
              assert(hit_single <= product_target_blocks);
              if (echo > 7) {
                  std::printf("# RHS#%i at %s reaches from (%g, %g, %g) to (%g, %g, %g)\n", 
                                irhs, str(global_source_coords[irhs]),
                                stats[X].min(), stats[Y].min(), stats[Z].min(),
                                stats[X].max(), stats[Y].max(), stats[Z].max());
                  // here, we can also check if the center of all targets is the source coordinate
              } // echo
              if (echo > 2) {
                  int total_checked{0};
                  // list in detail
                  for (int nci = 0; nci <= max_nci; ++nci) {
                      if (hist[nci] > 0) {
                          if (echo > 7 + 10*(0 != irhs)) {
                              std::printf("# RHS#%i has%9.3f k cases with %2d corners inside, d2 stats: %g +/- %g in [%g, %g] Bohr^2\n",
                                  irhs, hist[nci]*.001, nci, stats_d2[nci].mean(), stats_d2[nci].dev(), stats_d2[nci].min(), stats_d2[nci].max());
                          } // echo
                          total_checked += hist[nci];
                      } // hist[nci] > 0
                  } // nci
                  auto const partial = total_checked - hist[0] - hist[max_nci];
                  if (echo > 6) std::printf("# RHS#%i has %.3f k inside, %.3f k partial and %.3f k outside (of %.3f k checked blocks)\n",
                                irhs, hist[max_nci]*.001, partial*.001, hist[0]*.001, total_checked*.001);
                  inout[0].add(hist[max_nci]); inout[1].add(partial); inout[2].add(hist[0]); inout[3].add(total_checked);
              } // echo
              assert(idx3_diagonal > -1 && "difference vector (0,0,0) must be hit once");
              assert(tag_diagonal[idx3_diagonal] == irhs && "diagonal inconsistent");
          } // irhs
          if (echo > 0) {
              char const inout_class[][8] = {"inside", "partial", "outside",  "checked"};
              for (int i = 0; i < 4; ++i) {
                  std::printf("# RHSs have [%7d,%9.1f +/-%5.1f,%7d ] blocks %s\n",
                      int(inout[i].min()), inout[i].mean(), inout[i].dev(), int(inout[i].max()), inout_class[i]); 
              } // i
          } // echo

          // create a histogram about the distribution of the number of columns per row
          std::vector<uint32_t> hist(1 + nrhs, 0);
          for (size_t idx3 = 0; idx3 < column_indices.size(); ++idx3) {
              auto const n = column_indices[idx3].size();
              ++hist[n];
          } // idx3

          // eval the histogram
          size_t nall{0};
          size_t nnzb{0}; // number of non-zero BSR entries
          for (int n = 0; n <= nrhs; ++n) {
              nall += hist[n];
              nnzb += hist[n]*n;
          } // n
          if (echo > 7) {
              std::printf("# histogram total= %.3f k: ", nall*.001);
              printf_vector(" %d", hist.data(), nrhs + 1);
          } // echo
          assert(nall == product_target_blocks && "sanity check");

          p.nRows = nall - hist[0]; // the target block entries with no RHS do not create a row
          if (echo > 0) std::printf("# total number of Green function blocks is %.3f k, "
                               "average %.1f per source block\n", nnzb*.001, nnzb/double(nrhs));
          if (echo > 0) std::printf("# %.3f k (%.1f %% of %.3f k) target blocks are active\n", 
              p.nRows*.001, p.nRows/(product_target_blocks*.01), product_target_blocks*.001);

          assert(nnzb < (1ull << 32) && "the integer type of RowStart is uint32_t!");

          // resize BSR tables: (Block-compressed Sparse Row format)
          if (echo > 3) { std::printf("# memory of a complex Green function is %.6f %s (float, twice for double)\n",
                              nnzb*2.*64.*64.*sizeof(float)*GByte, _GByte); std::fflush(stdout); }
          p.colindx.resize(nnzb);
          p.rowindx  = get_memory<uint32_t>(nnzb, echo, "rowindx");
          p.RowStart = get_memory<uint32_t>(p.nRows + 1, echo, "RowStart");
          p.RowStart[0] = 0;
          p.veff_index = get_memory<int32_t>(nnzb, echo, "veff_index"); // indirection list for the local potential
          set(p.veff_index, nnzb, -1); // init as non-existing
          p.target_coords = get_memory<int16_t[3+1]>(p.nRows, echo, "target_coords");
          p.CubePos       = get_memory<float  [3+1]>(p.nRows, echo, "CubePos"); // internal coordinates but in float
          p.target_minus_source = get_memory<int16_t[3+1]>(nnzb, echo, "target_minus_source");

          p.global_target_indices.resize(p.nRows);
          p.subset.resize(p.nCols); // we assume columns of the unit operator as right-hand-sides

          view3D<int32_t> iRow_of_coords(num_target_coords[Z],
                                         num_target_coords[Y],
                                         num_target_coords[X], -1); // init as non-existing

          { // scope: fill BSR tables
              size_t warn_needs_shortest{0}, periodic_image_nontrivial{0};
              simple_stats::Stats<> st;
              uint32_t iRow{0};
              for (uint32_t z = 0; z < num_target_coords[Z]; ++z) { // serial
              for (uint32_t y = 0; y < num_target_coords[Y]; ++y) { // serial
              for (uint32_t x = 0; x < num_target_coords[X]; ++x) { // serial
                  uint32_t const idx[] = {x, y, z};
                  auto const idx3 = index3D(num_target_coords, idx);
                  assert(idx3 < product_target_blocks);

                  auto const n = column_indices[idx3].size();
                  if (n > 0) {
                      st.add(n);
                      iRow_of_coords(idx[Z], idx[Y], idx[X]) = iRow; // set existing

                      p.RowStart[iRow + 1] = p.RowStart[iRow] + n;
                      // copy the column indices
                      set(p.colindx.data() + p.RowStart[iRow], n, column_indices[idx3].data());
                      // copy the target block coordinates
                      int32_t global_target_coords[3];
                      for (int d = 0; d < 3; ++d) {
                          global_target_coords[d] = idx[d] + min_targets[d];
                          auto const internal_target_coord = global_target_coords[d] - global_internal_offset[d];
                          p.target_coords[iRow][d] = internal_target_coord; assert(internal_target_coord == p.target_coords[iRow][d] && "safe assign");
                          p.CubePos[iRow][d]       = internal_target_coord; assert(internal_target_coord == p.CubePos[iRow][d]       && "safe assign");
                      } // d
                      p.target_coords[iRow][3] = 0; // not used
                      p.CubePos[iRow][3] = 0.f; // not used

                      p.global_target_indices[iRow] = global_coordinates::get(global_target_coords); 
                      // global_target_indices are needed to gather the local potential data from other MPI processes

                      { // scope: determine the diagonal entry (source == target)
                          auto const iCol = tag_diagonal[idx3];
                          if (iCol > -1) {
                              assert(iCol < (1ul << 16)); // number range of uint16_t
                              for (int d = 0; d < 3; ++d) {
                                  // sanity check on internal coordinates
                                  if (p.source_coords[iCol][d] != p.target_coords[iRow][d])
                                      error("internal coordinates mismatch for RHS#%i at source(%s), target(%s) was tagged diagonal",
                                            iCol, str(p.source_coords[iCol]), str(p.target_coords[iRow]));
                                  assert(p.source_coords[iCol][d] == p.target_coords[iRow][d]);
                              } // d
                              { // search inz such that p.colindx[inz] == iCol
                                  int32_t inz_found{-1};
                                  for (auto inz = p.RowStart[iRow]; inz < p.RowStart[iRow + 1] && -1 == inz_found; ++inz) {
                                      if (iCol == p.colindx[inz]) inz_found = inz;
                                  } // inz
                                  assert(-1 != inz_found && "iCol should be in the list");
                                  p.subset[iCol] = inz_found;
                              } // search
                          } // iCol valid
                      } // scope: determine the diagonal entry

                      int32_t veff_index{-1};
                      if (1) { // scope: fill indirection table for having the local potential only defined in 1 unit cell and repeated periodically
                          int32_t mod[3];
                          bool potential_given{true};
                          for (int d = 0; d < 3; ++d) {
                              mod[d] = global_target_coords[d];
                              if (Periodic_Boundary == boundary_condition[d]) {
                                  mod[d] = global_target_coords[d] % n_original_Veff_blocks[d];
                                  mod[d] += (mod[d] < 0)*n_original_Veff_blocks[d];
                              } else {
                                  potential_given = potential_given && (mod[d] >= 0 && mod[d] < n_original_Veff_blocks[d]);
                              }
                          } // d
                          if (potential_given) {
                              auto const iloc = index3D(n_original_Veff_blocks, mod);
                              veff_index = iloc; assert(iloc == veff_index && "safe assign");
                          } else { // potential_given
                              assert(Vacuum_Boundary == boundary_condition[X] ||
                                     Vacuum_Boundary == boundary_condition[Y] ||
                                     Vacuum_Boundary == boundary_condition[Z]);
                              veff_index = -1; // outside of the unit cell due to Vacuum_Boundary
                          } // potential_given
                      } // scope: fill indirection table

                      for (auto inz = p.RowStart[iRow]; inz < p.RowStart[iRow + 1]; ++inz) {
                          auto const iCol = p.colindx[inz];

                          int iimage[] = {0, 0, 0};
                          if (is_periodic[ANY]) {
                              int32_t dv[3];
                              for (int d = 0; d < 3; ++d) {
                                  dv[d] = int32_t(p.target_coords[iRow][d]) - p.source_coords[iCol][d];
                              } // d
                              auto const *const n = n_original_Veff_blocks;
                              // with periodic boundary conditions, we need to find the shortest distance vector accounting for periodic images of the cell
                              double d2shortest{9e37};
                              for (int iz = -is_periodic[Z]; iz <= is_periodic[Z]; ++iz) {
                              for (int iy = -is_periodic[Y]; iy <= is_periodic[Y]; ++iy) {
                              for (int ix = -is_periodic[X]; ix <= is_periodic[X]; ++ix) {
                                  double const d2 = pow2((dv[X] + ix*n[X])*4*h[X])
                                                  + pow2((dv[Y] + iy*n[Y])*4*h[Y])
                                                  + pow2((dv[Z] + iz*n[Z])*4*h[Z]);
                                  if (d2 < d2shortest) {
                                      iimage[X] = ix;  iimage[Y] = iy;  iimage[Z] = iz;
                                      d2shortest = d2;
                                  } // is shorter
                              }}} // ix iy iz
                              if (pow2(iimage[X]) + pow2(iimage[Y]) + pow2(iimage[Z]) > 0) {
                                  if (echo > 7) std::printf("# block pair target(%s) - source(%s) has shortest distance %g %s at image(%s)\n",
                                                              str(p.target_coords[iRow]), str(p.source_coords[iCol]),
                                                              std::sqrt(d2shortest)*Ang, _Ang, str(iimage));
                                  ++periodic_image_nontrivial;
                              } // |iimage|^2 > 0
                              ++warn_needs_shortest; // still, the confinement potential may not be correct in each grid point
                          } // is_periodic[ANY]

                          for (int d = 0; d < 3; ++d) {
                              auto const diff = int32_t(p.target_coords[iRow][d]) - p.source_coords[iCol][d] + n_original_Veff_blocks[d]*iimage[d];
                              p.target_minus_source[inz][d] = diff; assert(diff == p.target_minus_source[inz][d] && "safe assign");
                          } // d
                          p.target_minus_source[inz][3] = 0; // not used
                          p.rowindx[inz] = iRow;
                          p.veff_index[inz] = veff_index;
                      } // inz

                      // count up the number of active rows
                      ++iRow;
                  } // n > 0
              }}} // idx
              assert(p.nRows == iRow && "counting 2nd time");
              assert(nnzb == p.RowStart[p.nRows] && "sparse matrix consistency");
              std::printf("# source blocks per target block in [%g, %.1f +/- %.1f, %g]\n", st.min(), st.mean(), st.dev(), st.max());
              if (warn_needs_shortest > 0) {
                  warn("for %.3f k block pairs the nearest periodic images are only evaluated on the block level", warn_needs_shortest*1e-3);
              } // warn_needs_shortest
              if (periodic_image_nontrivial > 0 && echo > 1) {
                  std::printf("# for %.3f k of %.3f k block pairs the nearest periodic image is non-trivial\n", periodic_image_nontrivial*1e-3, warn_needs_shortest*1e-3);
              } // periodic_image_nontrivial
          } // scope: fill BSR tables
          column_indices.clear(); // not needed beyond this point

          if (echo > 1) { // measure the difference in the number of target blocks of each RHS
              std::vector<uint32_t> nt(nrhs, 0);
              // traverse the BSR structure
              for (uint32_t iRow = 0; iRow < p.nRows; ++iRow) {
                  for (auto inz = p.RowStart[iRow]; inz < p.RowStart[iRow + 1]; ++inz) {
                      auto const iCol = p.colindx[inz];
                      ++nt[iCol];
                  } // inz
              } // iRow
              // analyze nt
              simple_stats::Stats<> st;
              for (uint16_t irhs = 0; irhs < nrhs; ++irhs) {
                  st.add(nt[irhs]);
              } // irhs
              std::printf("# target blocks per source block in [%g, %.1f +/- %.1f, %g]\n", st.min(), st.mean(), st.dev(), st.max());
          } // echo

          // Green function is stored sparse as real_t green[nnzb][2][Noco*64][Noco*64];

          for (int dd = 0; dd < 3; ++dd) { // derivate direction

              // create lists for the finite-difference derivatives
              auto const stat = green_kinetic::finite_difference_plan(p.kinetic_plan[dd]
                , dd
                , (1 == boundary_condition[dd]) // is periodic?
                , num_target_coords
                , p.RowStart, p.colindx.data()
                , iRow_of_coords
                , sparsity_pattern.data()
                , nrhs, echo);
              if (stat && echo > 0) std::printf("# construct_kinetic_plan in %c-direction returned status= %i\n", 'x' + dd, int(stat));

          } // dd derivate direction

          // transfer grid spacing into managed GPU memory
          p.grid_spacing = get_memory<double>(3, echo, "grid_spacing");
          set(p.grid_spacing, 3, hg); // for kinetic

          p.grid_spacing_trunc = get_memory<double>(3, echo, "grid_spacing_trunc");
          set(p.grid_spacing_trunc, 3, h); // customized grid spacings used for the construction of the truncation sphere

      } // scope

      auto const stat = construct_dyadic_plan(p.dyadic_plan
                            , cell, boundary_condition, hg
                            , AtomMatrices, xyzZinso
                            , p.nRows, p.nCols, p.RowStart, p.colindx.data()
                            , p.target_coords, r_block_circumscribing_sphere
                            , max_distance_from_center, r_trunc
                            , p.E_param, Noco
                            , echo);
      if (stat && echo > 0) std::printf("# construct_dyadic_plan returned status= %i\n", int(stat));

      auto const nerr = p.dyadic_plan.consistency_check();
      if (nerr && echo > 0) std::printf("# dyadic_plan.consistency_check has %d errors\n", nerr);

      int const n_iterations = control::get("green_function.benchmark.iterations", 1.); 
                      // -1: no iterations, 0:run memory initialization only, >0: iterate
      if (n_iterations < 0) {
          if (echo > 2) std::printf("# green_function.benchmark.iterations=%d --> no benchmarks\n", n_iterations);
      } else { // n_iterations < 0
          // try one of the 6 combinations (strangely, we cannot run any two of these calls after each other, ToDo: find out what's wrong here)
          int const action = control::get("green_function.benchmark.action", 412.);
          switch (action) {
              case 422: try_action<float ,2,2>(p, n_iterations, echo); break; // complex non-collinear
              case 822: try_action<double,2,2>(p, n_iterations, echo); break; // complex non-collinear

              case 412: try_action<float ,2,1>(p, n_iterations, echo); break; // complex
              case 812: try_action<double,2,1>(p, n_iterations, echo); break; // complex
#ifndef HAS_TFQMRGPU
              case 411: try_action<float ,1,1>(p, n_iterations, echo); break; // real
              case 811: try_action<double,1,1>(p, n_iterations, echo); break; // real
#else  // HAS_TFQMRGPU
              case 411:                                                       // real
              case 811: error("tfQMRgpu needs R1C2 == 2 but found green_function.benchmark.action=%d", action); break;
#endif // HAS_TFQMRGPU
              default: warn("green_function.benchmark.action must be in {411, 412, 422, 811, 812, 822} but found %d", action);
          } // switch action
      } // n_iterations < 0

      return 0;
  } // construct_Green_function

  #undef str // === vec2str.c_str()

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_Green_function(int const echo=0) {
      uint32_t ng[3] = {0, 0, 0}; // grid sizes
      int8_t   bc[3] = {0, 0, 0}; // boundary conditions
      double   hg[3] = {1, 1, 1}; // grid spacings
      std::vector<double> Veff(0); // local potential
      int natoms{0}; // number of atoms
      std::vector<double> xyzZinso(0); // atom info
      std::vector<std::vector<double>> AtomMatrices(0); // non-local potential

      auto const *const filename = control::get("green_function.hamiltonian.file", "Hmt.xml");
      auto const stat = green_input::load_Hamiltonian(ng, bc, hg, Veff, natoms, xyzZinso, AtomMatrices, filename, echo - 5);
      if (stat) {
          warn("failed to load_Hamiltonian with status=%d", int(stat));
          return stat;
      } // stat

      green_action::plan_t p;
      return construct_Green_function(p, ng, bc, hg, Veff, xyzZinso, AtomMatrices, echo);
  } // test_Green_function

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_Green_function(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_function
