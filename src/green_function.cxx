// This file is part of AngstromCube under MIT License

#include <cstdio>     // std::printf, ::snprintf, FILE, std::fprintf
#include <cstdint>    // int64_t, int32_t, uint32_t, int16_t, uint16_t, int8_t, uint8_t
#include <cassert>    // assert
#include <cmath>      // std::sqrt, ::cbrt
#include <algorithm>  // std::max, ::min
#include <utility>    // std::swap
#include <vector>     // std::vector<T>
#include <complex>    // std::complex

#include "green_function.hxx" // ::update_energy_parameter

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "simple_timer.hxx" // SimpleTimer
#include "simple_stats.hxx" // ::Stats<>
#include "display_units.h" // eV, _eV, Ang, _Ang
#include "print_tools.hxx" // printf_vector
#include "global_coordinates.hxx" // ::get
#include "green_input.hxx" // ::load_Hamitonian

#include "action_plan.hxx" // ::atom_t
#include "kinetic_plan.hxx" // kinetic_plan_t
#include "green_memory.hxx" // get_memory, free_memory, real_t_name
#include "green_sparse.hxx" // ::sparse_t<,>
#include "progress_report.hxx" // ProgressReport

#include "green_parallel.hxx" // ::potential_exchange, ::RequestList_t

#ifndef NO_UNIT_TESTS
  #include "mpi_parallel.hxx" // ::init, ::finalize, ::rank

  #ifdef HAS_TFQMRGPU

//  #define DEBUG
//  #define DEBUGGPU
    #ifdef HAS_NO_CUDA
        #include "tfQMRgpu/include/tfqmrgpu_cudaStubs.hxx" // cuda... (dummies)
        #define devPtr const __restrict__
    #else  // HAS_NO_CUDA
        #include <cuda.h>
    #endif // HAS_NO_CUDA
    #include "tfQMRgpu/include/tfqmrgpu.h" // ...
    #include "tfQMRgpu/include/tfqmrgpu.hxx" // ...
    #include "tfQMRgpu/include/tfqmrgpu_core.hxx" // tfqmrgpu::solve<action_t>

  #endif // HAS_TFQMRGPU

#endif // NO_UNIT_TESTS

#include "green_action.hxx" // ::plan_t, ::action_t, ::atom_t
#include "green_kinetic.hxx" // ::finite_difference_plan_t, index3D
#include "green_potential.hxx" // ::exchange
#include "green_dyadic.hxx" // ::dyadic_plan_t
#include "sho_tools.hxx" // ::nSHO
#include "control.hxx" // ::get
#include "load_balancer.hxx" // ::get
#include "boundary_condition.hxx" // Isolated_Boundary, Periodic_Boundary
//int8_t constexpr Isolated_Boundary = 0, Periodic_Boundary = 1;

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

namespace green_function {

  double const GByte = 1e-9; char const *const _GByte = "GByte";

  int constexpr X=0, Y=1, Z=2;

  template <typename number_t>
  std::string vec2str(number_t const vec[3], double const f=1, char const *const sep=" ") {
      // convert a vector of 3 numbers into a string, with scaling f and 2 separators
      char s[64];
      std::snprintf(s, 64, "%g%s%g%s%g", vec[X]*f, sep, vec[Y]*f, sep, vec[Z]*f);
      return std::string(s);
  } // vec2str
  #define str(...) vec2str(__VA_ARGS__).c_str()

  // ToDo: move special BCs into headers
  int8_t constexpr Vacuum_Boundary = 2;
  // The vacuum boundary condition is an addition to Isolated_Boundary and Periodic_Boundary from boundary_condition.hxx
  // Vacuum_Boundary means that the Green function extends from its source coordinate up to the truncation radius
  // even if this exceeds the isolated boundary at which potential values stop to be defined.
  // The potential is continued as zero beyond the isolated boundary but the tail of the Green function is allowed to fade there.
  // k-points are not relevant.

  // We could think of a useful new boundary condition
  int8_t constexpr Repeat_Boundary = 3;
  // The repeat boundary allows to compute with small unit cells but with realistic truncation radii.
  // The local potential is repeated periodically (like Periodic_Boundary).
  // Sources are only relevant inside the small unit cell but again, the Green function tail exceeds the cell boundaries
  // and extends up to the truncation radius.
  // For the dyadic potential, we may not reduce over periodic atom images but create copies of the atoms.
  // k-points are not relevant.

  // The wrap boundary condition is an addition to Periodic_Boundary from boundary_condition.hxx
  int8_t constexpr Wrap_Boundary = 5;
  // Wrap_Boundary means that the truncation sphere fits into the cell, so k-points have no effect.
  // Nevertheless, it cannot be treated like Isolated_Boundary since target block coordinates may need to be wrapped.
  // However, it could be viewed as Repeat_Boundary...


    // ToDo: make it a method of action_plan_t
    status_t update_energy_parameter(
          action_plan_t & plan
        , std::complex<double> E_param
        , std::vector<std::vector<double>> const & AtomMatrices
        , double const dVol // volume element of the grid
        , double const scale_H // =1
        , int const echo // =0
        , int const Noco // =1
        , green_parallel::RequestList_t const *requests // =nullptr
    ) {
        plan.E_param = E_param;

        auto & p = plan.dyadic_plan;
        if (1 != Noco && p.nAtoms > 0) warn("not prepared for Noco=%d", Noco);

        view2D<double> output;
        if (requests) {
            size_t nc2{0};
            for (auto const & am : AtomMatrices) { nc2 = std::max(nc2, am.size()); }
            int const count = mpi_parallel::max(nc2);
            if (echo > 13) std::printf("# %s: AtomMatrices.size()= %ld, requests->window()= %ld, p.nAtoms= %d, requests->size()= %ld\n",
                                    __func__, AtomMatrices.size(),      requests->window(),      p.nAtoms,     requests->size());
            assert(p.nAtoms == requests->size());
            assert(AtomMatrices.size() == requests->window());
            view2D<double> input(requests->window(), count, 0.0);
            for (size_t iam{0}; iam < AtomMatrices.size(); ++iam) { 
                auto const & am = AtomMatrices.at(iam);
                auto const nc2 = am.size();
                set(input[iam], nc2, am.data()); // copy
            } // iam
            output = view2D<double>(requests->size(), count, 0.0);
            green_parallel::exchange(output.data(), input.data(), *requests, count, echo, "atom_mat");
        } else {
            if (echo > 0) std::printf("# skip green_function.matrices.exchange\n");
        } // requests

        auto const f = dVol; // we multiply the matrices by dVol so we can omit this factor in SHOprj
        for (int iac = 0; iac < p.nAtoms; ++iac) { // contributing atoms can be processed in parallel
            int const lmax = p.AtomLmax[iac];
            // auto const sigma = p.AtomSigma[iac];
            int const nc = sho_tools::nSHO(lmax);
            assert(nc > 0); // the number of coefficients of contributing atoms must be non-zero
            // auto sho_norm = green_dyadic::sho_normalization(lmax, sigma);
            // assert(sho_norm.size() == nc);
            // for (int i = 0; i < nc; ++i) {
            //     sho_norm[i] = std::sqrt(dVol/sho_norm[i]);
            // } // i
            double const *hmt{nullptr};
            if (requests) {
                hmt = output[iac];
            } else {
                auto const ia = p.original_atom_index[iac];
                assert(2*nc*nc <= AtomMatrices.at(ia).size());
                hmt = AtomMatrices[ia].data(); // Hamiltonian matrix elements
            }

            // fill this with matrix values
            auto const atomMatrix = p.AtomMatrices[iac]; // data layout [Noco*Noco*2*nc*nc], 2 for {real,imag}
            auto const ovl = hmt + nc*nc;          // charge deficit matrix elements
            for (int i = 0; i < nc; ++i) {
                for (int j = 0; j < nc; ++j) {
                    int const ij = i*nc + j;
                    atomMatrix[ij] = (scale_H * hmt[ij] - E_param.real() * ovl[ij])*f; // real part
                    atomMatrix[ij + nc*nc] = (          - E_param.imag() * ovl[ij])*f; // imag part
                } // j
            } // i
            // TODO treat Noco components correctly: component 0 and 1 == V_upup and V_dndn
            //                                       component 2 and 3 == V_x and V_y
        } // iac

        return 0;
    } // update_atom_matrices



    // ToDo: make it a method of action_plan_t
    status_t update_potential(
          action_plan_t & p // inout, create a plan how to apply the SHO-PAW Hamiltonian to a block-sparse truncated Green function
        , uint32_t const nb[3] // numbers of 4*4*4 grid blocks of the unit cell in with the potential is defined
        , std::vector<double> const & Veff // [nb[2]*4*nb[1]*4*nb[0]*4]
        , view3D<uint16_t> const & owner_rank // [nb[2]*nb[1]*nb[0]]
        , int const echo // =0 // verbosity
        , int const Noco // =2
    ) {
        auto const n_all_grid_points = size_t(nb[Z]*4)*size_t(nb[Y]*4)*size_t(nb[X]*4);
        if (Veff.size() == n_all_grid_points) { // scope: restructure and communicate potential
            auto const nrhs = p.nCols;
//          auto const n_all_blocks = size_t(nb[2])*size_t(nb[1])*size_t(nb[0]);
            assert(owner_rank.stride() == nb[0] && owner_rank.dim1() == nb[1]);

            auto const Vinp = new double[nrhs*Noco*Noco][64];
            // reorder Veff[ng[Z]*ng[Y]*ng[X]] into block-structured Vinp
            for (uint16_t irhs{0}; irhs < nrhs; ++irhs) {
                // assume that in this MPI rank the potential values of the right-hand-sides that are to be determined are known
                int32_t ib[3]; global_coordinates::get(ib, p.global_source_indices[irhs]);
                for (int d = 0; d < 3; ++d) { assert(ib[d] >= 0); }
                for (int i4z = 0; i4z < 4; ++i4z) { size_t const iz = ib[Z]*4 + i4z;
                for (int i4y = 0; i4y < 4; ++i4y) { size_t const iy = ib[Y]*4 + i4y;
                for (int i4x = 0; i4x < 4; ++i4x) { size_t const ix = ib[X]*4 + i4x;
                    auto const izyx = (iz*(nb[Y]*4) + iy)*(nb[X]*4) + ix; // global grid point index
                    assert(izyx < n_all_grid_points);
                    auto const i64 = (i4z*4 + i4y)*4 + i4x;
                    if (2 == Noco) {
                        Vinp[irhs*4 + 3][i64] = 0.0;  // set clear V_y
                        Vinp[irhs*4 + 2][i64] = 0.0;  // set clear V_x
                        Vinp[irhs*4 + 1][i64] = Veff[izyx]; // set V_upup
                    } // non-collinear
                    Vinp[irhs*Noco*Noco][i64] = Veff[izyx]; // copy potential value to V_dndn
                }}} // i4x i4y i4z
            } // irhs

#ifdef    HAS_NO_MPI
            auto const default_exchange = 0.; // skip the exchange since we have no parallel processes, beware that p.Veff is set to 0.0
#else  // HAS_NO_MPI
            auto const default_exchange = 1.; // do the exchange of potential blocks
#endif // HAS_NO_MPI
            if (1. == control::get("green_function.potential.exchange", default_exchange)) {
                // ToDo: move this into the constructors
                green_parallel::RequestList_t const requests(p.global_target_indices,  // requests
                                                             p.global_source_indices, // offerings
                                                             owner_rank.data(), nb, echo);
                green_parallel::potential_exchange(p.Veff, Vinp, requests, Noco, echo);
            } else {
                if (echo > 0) std::printf("# skip green_function.potential.exchange\n");
            } // needs exchange?

            delete[] Vinp;

            if (echo > 4) {
                int constexpr mag = 0; // only for the Noco=1 case
                simple_stats::Stats<> pot;
                for (int iRow = 0; iRow < p.nRows; ++iRow) {
                    auto const *const V = p.Veff[mag][iRow];
                    for (int i64 = 0; i64 < 64; ++i64) {
                        pot.add(V[i64]);
                    } // i64
                } // iRow
                std::printf("# %s effective local potential %s %s\n", __func__, pot.interval(eV).c_str(), _eV);
            } // echo

        } else {
            warn("wrong number of potential grid points found %ld expect %lld", Veff.size(), n_all_grid_points);
        }
        return 0;
    } // update_potential







#ifdef    GREEN_FUNCTION_SVG_EXPORT
  static FILE* svg;
#endif // GREEN_FUNCTION_SVG_EXPORT


  status_t construct_dyadic_plan(
        dyadic_plan_t & p
      , double const cell[3]
      , int8_t const boundary_condition[3]
      , double const grid_spacing[3]
//    , std::vector<std::vector<double>> const & AtomMatrices // [natoms][2*nSHO(numax[ia])]
      , std::vector<double> const & xyzZinso // [natoms*8]
      , uint32_t const nRowsGreen
      , uint32_t const nrhs
      , uint32_t const *const rowStartGreen
      , uint16_t const *const colIndexGreen
      , int16_t const (*internal_target_coords)[3+1]
      , int32_t const global_internal_offset[3]
      , double const r_block_circumscribing_sphere
      , double const max_distance_from_center
      , double const r_trunc
      , int const echo=0 // verbosity
      , int const Noco=1 // 1:collinear spins, 2:Non-collinear
  ) {
      SimpleTimer timer(__FILE__, __LINE__, __func__, echo);

      p.nrhs = nrhs;

      // transfer grid spacing into managed GPU memory
      p.grid_spacing = get_memory<double>(3+1, echo, "grid_spacing");
      set(p.grid_spacing, 3, grid_spacing);
      double const r_proj = control::get("green_function.projection.radius", 6.); // in units of sigma
      p.grid_spacing[3] = r_proj; // radius in units of sigma at which the projectors stop

      int32_t const natoms = xyzZinso.size()/8; // number of original atoms
      if (echo > 2) std::printf("\n#\n# %s for %d atoms\n#\n", __func__, natoms);

      // compute which atoms will contribute, the list of natoms atoms may contain a subset of all atoms
      double max_projection_radius{0};
      for (int32_t ia = 0; ia < natoms; ++ia) { // loop over all original atoms (can be parallel with reduction)
          auto const sigma = xyzZinso[ia*8 + 6];
          auto const projection_radius = std::max(0.0, r_proj*sigma);
          max_projection_radius = std::max(max_projection_radius, projection_radius);
      } // ia
      if (echo > 3) std::printf("# largest projection radius is %g %s\n", max_projection_radius*Ang, _Ang);

      auto const radius = r_trunc + max_distance_from_center + 2*max_projection_radius + 2*r_block_circumscribing_sphere;

      int iimage[3] = {0, 0, 0}; // number of images replications of the unit cell to each side
      int icopies[3] = {0, 0, 0}; // number of copied atoms
      size_t nimages{1}, ncopies{1};
      for (int d = 0; d < 3; ++d) { // parallel
          // what happens if there is a wrap_boundary but the non-local projection overlaps with two periodic ends of the sphere?
          // we may have to diffentiate between images and copies of atoms! Repeat_Boundary needs copies, Wrap_boundary needs copies
          int const nmx = std::max(0., std::ceil(radius/cell[d]));
          if (Periodic_Boundary == boundary_condition[d]) {
              iimage[d] = nmx;
          } else if (Repeat_Boundary == boundary_condition[d]) {
              icopies[d] = nmx;
          } else if (Wrap_Boundary == boundary_condition[d]) {
              icopies[d] = 1; // in WRAP the truncation sphere fits into the cell, ...
              // however, left targets of a left source and right targets of a right source could potentially see the same atom
          }
          nimages *= (2*iimage[d]  + 1); // iimage  images to the left and right
          ncopies *= (2*icopies[d] + 1); // icopies copies to the left and right
      } // d
      size_t const nAtomCopies = natoms*ncopies;
      if (ncopies > 1 && echo > 3) std::printf("# copy %d atoms %s times (%ld times) for %ld atom copies\n", natoms, str(icopies), ncopies, nAtomCopies);
      assert(ncopies >= 1); // later we divide by ncopies to retrieve the original atom index
      size_t const nAtomImages = nAtomCopies*nimages;
      if (echo > 3) std::printf("# replicate %s atom images, %ld images total\n", str(iimage), nimages);

      std::vector<uint32_t> AtomImageStarts(nAtomImages + 1, 0); // probably larger than needed, should call resize(nai + 1) later
      std::vector<action_plan::atom_t> atom_data(nAtomImages);
      std::vector<int8_t> atom_numax(nAtomCopies, -1); // -1: atom does not contribute

      simple_stats::Stats<> nc_stats;
      double sparse{0}, dense{0}; // stats to assess how much memory can be saved using sparse storage

      std::vector<std::vector<uint32_t>> cubes; // stores the row indices of Green function rows
      cubes.reserve(nAtomImages); // maximum (needs 24 Byte per atom image)

      size_t iai{0}; // counter for atomic images
{   SimpleTimer timer(strip_path(__FILE__), __LINE__, "computing distances with all atoms", echo);

      for (int z = -iimage[Z]; z <= iimage[Z]; ++z) { // serial
      for (int y = -iimage[Y]; y <= iimage[Y]; ++y) { // serial
      for (int x = -iimage[X]; x <= iimage[X]; ++x) { // serial
//            if (echo > 3) std::printf("# periodic shifts  %d %d %d\n", x, y, z);
          int const xyz_image_shift[] = {x, y, z};
          for (int ia = 0; ia < natoms; ++ia) { // loop over original atoms in the unit cell, serial
              auto const atom_id = int32_t(xyzZinso[ia*8 + 4]); // global_atom_id
              auto const numax =       int(xyzZinso[ia*8 + 5]);
              auto const sigma =           xyzZinso[ia*8 + 6] ;
              int iaa{0};
          for (int zc = -icopies[Z]; zc <= icopies[Z]; ++zc) { // serial
          for (int yc = -icopies[Y]; yc <= icopies[Y]; ++yc) { // serial
          for (int xc = -icopies[X]; xc <= icopies[X]; ++xc) { // serial
              int const xyz_copy_shift[] = {xc, yc, zc};
              // suggest a shifted atomic image position
              double atom_pos[3];
              for (int d = 0; d < 3; ++d) { // unroll
                  atom_pos[d] = xyzZinso[ia*8 + d] + (xyz_image_shift[d] + xyz_copy_shift[d] + 0.5)*cell[d];
              } // d
//            if (echo > 5) std::printf("# image of atom #%i at %s %s\n", atom_id, str(atom_pos, Ang), _Ang);

              double const r_projection = r_proj*sigma; // atom-dependent, precision dependent, assume float here
              double const r2projection = pow2(r_projection);
//            double const r2projection_plus = pow2(r_projection + r_block_circumscribing_sphere);

              // check all target blocks if they are inside the projection radius
              uint32_t ntb{0}; // number of target blocks
              for (uint32_t icube = 0; icube < nRowsGreen; ++icube) { // loop over blocks
                  auto const *const target_block = internal_target_coords[icube];
                  if (true) {
                      // do more precise checking
//                    if (echo > 9) std::printf("# target block #%i at %s gets corner check\n", icube, str(target_block));
                      int nci{0}; // number of corners inside
                      // check 8 corners
                      for (int iz = 0; iz < 4; iz += 3) { // parallel, reduction on nci
                      for (int iy = 0; iy < 4; iy += 3) { // parallel, reduction on nci
                      for (int ix = 0; ix < 4; ix += 3) { // parallel, reduction on nci
                          int const ixyz[] = {ix, iy, iz};
                          double d2i{0}; // init distance^2 of the grid point from the center
                          for (int d = 0; d < 3; ++d) {
                              double const grid_point = ((target_block[d] + global_internal_offset[d])*4 + ixyz[d] + 0.5)*grid_spacing[d];
                              d2i += pow2(grid_point - atom_pos[d]);
                          } // d
                          if (d2i < r2projection) {
                              ++nci; // at least one corner of the block is inside the projection radius of this atom
                          } // inside the projection radius
                      }}} // ix // iy // iz
                      // three different cases: 0, 1...7, 8, i.e. none, partial, full
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
                          int const nrhs_icube = rowStartGreen[icube + 1] - rowStartGreen[icube];
                          sparse += nrhs_icube;
                          dense  += p.nrhs; // all columns

//                        if (echo > 7) std::printf("# target block #%i at %s is inside\n", icube, str(target_block));
                      } else { // nci
//                        if (echo > 9) std::printf("# target block #%i at %s is outside\n", icube, str(target_block));
                      } // nci

                  } else { // d2 < r2projection_plus
                      assert(false);
//                    if (echo > 21) std::printf("# target block #%i at %s is far outside\n", icube, str(target_block));
                  } // d2 < r2projection_plus
              } // icube

              if (ntb > 0) {
                  // atom image contributes, mark in the list to have more than 0 coefficients
                  auto const nc = sho_tools::nSHO(numax); // number of coefficients for this atom
                  atom_numax[iaa] = numax; // atom image contributes, so this atom contributes
                  nc_stats.add(nc);

                  // at least one target block has an intersection with the projection sphere of this atom image
                  auto & atom = atom_data[iai];
                  set(atom.pos, 3, atom_pos);
                  atom.sigma = sigma;
                  atom.gid = atom_id;
                  atom.ia = ia;
                  atom.iaa = iaa; // local atom index, where is this used?
                  set(atom.shifts, 3, xyz_image_shift);
                  set(atom.copies, 3, xyz_copy_shift);
                  atom.nc = nc;
                  atom.numax = numax;

//                if (echo > 5) std::printf("# image of atom #%i at %s %s contributes to %d target blocks\n", atom_id, str(atom_pos, Ang), _Ang, ntb);
                  AtomImageStarts[iai + 1] = AtomImageStarts[iai] + nc;
                  ++iai;

#ifdef    GREEN_FUNCTION_SVG_EXPORT
                  if (nullptr != svg) { // these circles show the projection spheres of the atoms
                      std::fprintf(svg, "  <ellipse cx=\"%g\" cy=\"%g\" rx=\"%g\" ry=\"%g\" fill=\"none\" stroke=\"red\" />\n",
                                             atom_pos[X]/grid_spacing[X],  atom_pos[Y]/grid_spacing[Y],
                                            r_projection/grid_spacing[X], r_projection/grid_spacing[Y]);
                      std::fprintf(svg, "  <ellipse cx=\"%g\" cy=\"%g\" rx=\".1\" ry=\".1\" fill=\"none\" stroke=\"red\" />\n",
                                             atom_pos[X]/grid_spacing[X],  atom_pos[Y]/grid_spacing[Y]); // center
                      auto const Z = xyzZinso[ia*8 + 3];
                      if (Z > 0) std::fprintf(svg, "  <text x=\"%g\" y=\"%g\" style=\"font-size: 3;\">%g</text>\n",
                                             atom_pos[X]/grid_spacing[X],  atom_pos[Y]/grid_spacing[Y], Z);
                  } // nullptr != svg
#endif // GREEN_FUNCTION_SVG_EXPORT

              } else {
//                if (echo > 15) std::printf("# image of atom #%i at %s %s does not contribute\n", atom_id, str(atom_pos, Ang), _Ang);
              } // ntb > 0
              ++iaa;
      }}} // xc // yc // zc  --- copies
          } // ia           --- atoms
      }}} // x // y // z  --- images

} // timer

      auto const nai = iai; // corrected number of atomic images
      if (echo > 3) std::printf("# %ld of %lu (%.2f %%) atom images have an overlap with projection spheres\n",
                                    nai, nAtomImages, nai/std::max(nAtomImages*.01, .01));
      auto const napc = AtomImageStarts[nai]; // number of atomic projection coefficients

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
      for (unsigned irhs = 0; irhs < nrhs; ++irhs) {
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
          char name[64]; std::snprintf(name, 64, "sparse_SHOprj[irhs=%i of %d]", irhs, nrhs);
          p.sparse_SHOprj[irhs] = sparse_t<>(SHOprj[irhs], false, name, echo - 4);
          // sparse_SHOprj: rows == atom images, cols == Green function non-zero elements
          nops += p.sparse_SHOprj[irhs].nNonzeros();
      } // irhs
      if (echo > 5) std::printf("# sparse_SHOprj has %d*%ld rows, %d columns and %ld nonzeros\n", nrhs,nai, nnzb, nops);
      SHOprj.resize(0); // release host memory

      p.AtomImageStarts = get_memory<uint32_t>(nai + 1, echo, "AtomImageStarts");
      set(p.AtomImageStarts, nai + 1, AtomImageStarts.data()); // copy into GPU memory
      p.AtomImageIndex = get_memory<uint32_t>(nai, echo, "AtomImageIndex");

      // get all info for the atomic matrices
      p.global_atom_index.resize(nAtomCopies); // translation table
      p.original_atom_index.resize(nAtomCopies); // translation table
      std::vector<int64_t> local_atom_index(nAtomCopies, -1); // translation table
      size_t iac{0};
      for (int iaa = 0; iaa < nAtomCopies; ++iaa) { // serial loop over all atom copies
          if (atom_numax[iaa] > -1) {
              p.global_atom_index[iac] = iaa;
              int const ia = iaa/ncopies; // integer division
              p.original_atom_index[iac] = ia;
              local_atom_index[iaa] = iac;
              ++iac;
          } // atom contributes
      } // iaa
      auto const nac = iac; // number of contributing atoms
      p.global_atom_index.resize(nac);
      p.original_atom_index.resize(nac);
      p.global_atom_ids.resize(nac, -1);

      // store the atomic image positions in GPU memory
      p.AtomImageLmax  = get_memory<int8_t>(nai, echo, "AtomImageLmax");
      p.AtomImagePos   = get_memory<double[3+1]>(nai, echo, "AtomImagePos");
      p.AtomImagePhase = get_memory<double[4]>(nai, echo, "AtomImagePhase");
      p.AtomImageShift = get_memory<int8_t[4]>(nai, echo, "AtomImageShift");
      double const phase0[] = {1, 0, 0, 0};

      std::vector<std::vector<uint32_t>> SHOsum(nac);

      for (size_t iai = 0; iai < nai; ++iai) { // serial
          auto const iaa = atom_data[iai].iaa; // index of the atom copy
          auto const iac = local_atom_index[iaa]; assert(0 <= iac); assert(iac < nac);
          p.global_atom_ids.at(iac) = atom_data[iai].gid;
          p.AtomImageIndex[iai] = iac;
          set(p.AtomImagePos[iai], 3, atom_data[iai].pos);
          set(p.AtomImageShift[iai], 3, atom_data[iai].shifts); p.AtomImageShift[iai][3] = 0;
          p.AtomImagePos[iai][3] = 1./std::sqrt(atom_data[iai].sigma);
          p.AtomImageLmax[iai] = atom_data[iai].numax; // SHO basis size
          if (echo > 9) std::printf("# atom image #%ld has lmax= %d\n", iai, p.AtomImageLmax[iai]);
          SHOsum[iac].push_back(uint32_t(iai));
          set(p.AtomImagePhase[iai], 4, phase0); // TODO construct correct phases, phase0 is just a dummy
      } // iai, copy into GPU memory

      p.sparse_SHOsum = sparse_t<>(SHOsum, false, "sparse_SHOsum", echo);
      SHOsum.resize(0);


      // get memory for the matrices and fill
      p.AtomMatrices = get_memory<double*>(nac, echo, "AtomMatrices");
      p.AtomLmax     = get_memory<int8_t>(nac, echo, "AtomLmax");
      p.AtomStarts   = get_memory<uint32_t>(nac + 1, echo, "AtomStarts");
      p.AtomStarts[0] = 0; // init prefetch sum
      // p.AtomSigma.resize(nac);

      for (uint32_t iac = 0; iac < nac; ++iac) { // parallel
          auto const iaa = p.global_atom_index[iac];
          // p.AtomSigma[iac] = xyzZinso[ia*8 + 6];
          uint32_t const nc = sho_tools::nSHO(atom_numax[iaa]);
          assert(nc > 0); // the number of coefficients of a contributing atom copy must be non-zero
          p.AtomStarts[iac + 1] = p.AtomStarts[iac] + nc; // create prefetch sum
          p.AtomLmax[iac] = atom_numax[iaa];
          char name[64]; std::snprintf(name, 64, "AtomMatrices[iac=%d/iaa=%d]", iac, iaa);
          p.AtomMatrices[iac] = get_memory<double>(Noco*Noco*2*nc*nc, echo, name);
          set(p.AtomMatrices[iac], Noco*Noco*2*nc*nc, 0.0); // clear
      } // iac
      p.nAtoms = nac; // number of contributing atom copies

      warn("missing call to update_atom_matrices, Noco=%d", Noco);

      if (echo > 1) std::printf("# found %lu contributing atoms with %lu atom images\n", nac, nai);

      p.update_flop_counts(echo); // prepare to count the number of floating point operations

      return 0;
  } // construct_dyadic_plan


  status_t update_phases(
        action_plan_t & p
      , double const k_point[3]
      , int const echo // =0 // verbosity
      , int const Noco // =1
  ) {
      assert(p.phase && "phase[3][2][2] must already be allocated in device memory");
      green_kinetic::set_phase(p.phase, k_point, echo); // neutral (Gamma-point) phase factors

      std::complex<double> phase[3][2]; // for each direction: forward and backward phase
      for (int d = 0; d < 3; ++d) {
//        double const arg = 2*constants::pi*k_point[d];
//        phase[d][0] = std::complex<double>(std::cos(arg), std::sin(arg));
//        phase[d][1] = std::conj(phase[d][0]); // must be the inverse if the k_point gets an imaginary part
          phase[d][0] = std::complex<double>(p.phase[d][0][0], p.phase[d][0][1]);
          phase[d][1] = std::complex<double>(p.phase[d][1][0], p.phase[d][1][1]);
      } // d

      auto & dp = p.dyadic_plan;
      assert(dp.AtomImagePhase && "AtomImagePhase must already be allocated");
      for (uint32_t iai = 0; iai < dp.nAtomImages; ++iai) { // parallel
          std::complex<double> ph(1, 0);
          for (int d = 0; d < 3; ++d) {
              auto const shift = dp.AtomImageShift[iai][d];
              ph *= (shift >= 0) ? intpow(phase[d][0], shift) : intpow(phase[d][1], -shift);
          } // d
          set(dp.AtomImagePhase[iai], 4, 0.0);
          dp.AtomImagePhase[iai][0] = ph.real();
          dp.AtomImagePhase[iai][1] = ph.imag();
          assert(1 == Noco);
      } // iai

      return 0;
  } // update_phases



  std::vector<int64_t> get_right_hand_sides(
        uint32_t const nb[3] // number of blocks
      , std::vector<uint16_t> & owner_rank // result: who owns which RHS block
      , int const echo=0
  ) {
      std::vector<int64_t> global_source_indices; // result array

      auto const comm = mpi_parallel::comm();
      int const true_comm_size = mpi_parallel::size(comm);
      int const fake_comm = (true_comm_size > 1) ? 0 : control::get("green_function.fake.comm", 0.);
      auto const comm_size = (fake_comm > 0) ? fake_comm : true_comm_size;
      auto const nall = size_t(nb[Z])*size_t(nb[Y])*size_t(nb[X]);
      owner_rank.resize(nall, 0);
      if (comm_size > 1) {

          if (echo > 3) std::printf("# MPI parallelization of %.3f k right hand sides\n", nall*1e-3);
          assert(nall > 0);
          int const comm_rank = (fake_comm > 0) ? control::get("green_function.fake.rank", fake_comm - 1.)
                                                : mpi_parallel::rank(comm);
          double rank_center[4]; // rank_center[0/1/2] are the coordinates of the center of weight of the RHSs assigned to this rank
                                 // rank_center[3] is the number of tasks with nonzero weight
          load_balancer::get(comm_size, comm_rank, nb, echo, rank_center, owner_rank.data());
          if (fake_comm < 1) mpi_parallel::max(owner_rank.data(), nall); // MPI_Allreduce(MPI_MAX)
          if (echo > 9) {
              std::printf("# rank#%i owner_rank after  MPI_MAX ", comm_rank);
              printf_vector(" %i", owner_rank);
          } // echo
          auto const nrhs = size_t(rank_center[3]); // number of tasks with nonzero weight
          if (echo > 5) std::printf("# rank#%d of %d procs has %ld tasks\n", comm_rank, comm_size, nrhs);
          {
              simple_stats::Stats<> nt; // number of tasks
              nt.add(nrhs);
              mpi_parallel::allreduce(nt, comm);
              if (echo > 4) std::printf("# number of tasks per rank is in [%g, %g +/- %g, %g]\n", nt.min(), nt.mean(), nt.dev(), nt.max());
          }

          {
              global_source_indices.resize(nrhs, -1);
              uint32_t irhs{0};
              int64_t const nbX = nb[X], nbY = nb[Y];
              for (size_t iall = 0; iall < nall; ++iall) {
                  if (comm_rank == owner_rank[iall]) {
                      int32_t const iz = iall/(nbX*nbY);
                      int32_t const iy = (iall - iz*nbX*nbY)/nbX;
                      int32_t const ix = iall - ((iz*nbY) + iy)*nbX;
                      assert(iall == (iz*nbY + iy)*nbX + ix);
                      global_source_indices[irhs] = global_coordinates::get(ix, iy, iz);
                      ++irhs;
                  } // owned
              } // iall
              if (nrhs != irhs) error("rank#%i number of right hand sides=%d inconsistent with number of owned tasks=%d", comm_rank, nrhs, irhs);
              assert(nrhs == irhs && "number of right hand sides inconsistent");
          }

      } else { // comm_size > 1

#ifdef    HAS_NO_MPI
          auto const default_sources =  1.; // 1: 1x1x1 right-hand-side only (suitable default for ./a43 --test green_function)
#else  // HAS_NO_MPI
          auto const default_sources = -1.; // -1: all right-hand-sides (suitable default for ./green --test green_function)
#endif // HAS_NO_MPI
          // generate a box of source points
          double nsb[3] = {0, 0, 0}; // number of source blocks
          int32_t const source_cube = control::get(nsb, "green_function.sources", "xyz", default_sources);
          int32_t n_source_blocks[] = {int(nsb[X]), int(nsb[Y]), int(nsb[Z])};
          int32_t off[3];
          for (int d = 0; d < 3; ++d) {
              n_source_blocks[d] = (n_source_blocks[d] < 0) ? nb[d] : // "green_function.sources" negative means all
                    std::min(std::max(1, n_source_blocks[d]), int32_t(nb[d])); // clamp
              off[d] = (nb[d] - n_source_blocks[d])/2; // place at center of the unit cell
          } // d
          if (source_cube && echo > 0) std::printf("\n# use green_function.sources= %d x %d x %d\n",
                                        n_source_blocks[X], n_source_blocks[Y], n_source_blocks[Z]);

          // total number of right-hand-sides treated in this MPI process
          auto const nrhs = size_t(n_source_blocks[Z])*size_t(n_source_blocks[Y])*size_t(n_source_blocks[X]);

          if (echo > 3) std::printf("# number of RHS source blocks %s = %ld\n", str(n_source_blocks, 1, " x "), nrhs);

          global_source_indices.resize(nrhs, -1);
          uint32_t irhs{0};
          for (int32_t ibz = 0; ibz < n_source_blocks[Z]; ++ibz) {
          for (int32_t iby = 0; iby < n_source_blocks[Y]; ++iby) {
          for (int32_t ibx = 0; ibx < n_source_blocks[X]; ++ibx) {
              global_source_indices[irhs] = global_coordinates::get(ibx + off[X], iby + off[Y], ibz + off[Z]);
              ++irhs;
          }}} // xyz
          assert(nrhs == irhs);

      } // comm_size > 1

      if (echo > 5) std::printf("# total number of source blocks is %ld\n", global_source_indices.size());
      return global_source_indices;
  } // get_right_hand_sides


  status_t construct_Green_function(
        action_plan_t & p // result, create a plan how to apply the SHO-PAW Hamiltonian to a block-sparse truncated Green function
      , uint32_t const ng[3] // numbers of grid points of the unit cell in with the potential is defined
      , int8_t const boundary_condition[3] // boundary conditions in {Isolated, Periodic, Vacuum, Repeat}
      , double const hg[3] // grid spacings
//    , std::vector<double> const & Veff // [ng[2]*ng[1]*ng[0]]
      , std::vector<double> const & xyzZinso // [natoms*8]
//    , std::vector<std::vector<double>> const & AtomMatrices // atomic hamiltonian and overlap matrix, [natoms][2*nsho^2]
      , int const echo // =0 // log-level
//    , std::complex<double> const *energy_parameter // =nullptr // E in G = (H - E*S)^{-1}
      , int const Noco // =2
  ) {
      if (echo > 0) std::printf("\n#\n# %s(%s)\n#\n\n", __func__, str(ng, 1, ", "));

      p.E_param = 0;

      int8_t bc[3] = {boundary_condition[X], boundary_condition[Y], boundary_condition[Z]};
      uint32_t n_blocks[3] = {0, 0, 0};
      for (int d = 0; d < 3; ++d) {
          n_blocks[d] = (ng[d] >> 2); // divided by 4
          assert(n_blocks[d] > 0 && "Needs at least one block (=4 grid points) per direction");
          assert(ng[d] == 4*n_blocks[d] && "All grid dimensions must be a multiple of 4!");
          assert(n_blocks[d] <= (1ul << 21) && "Max grid is 2^21 blocks due to global_coordinates");
          assert(hg[d] > 0. && "Needs a positive grid spacing");
          assert(Isolated_Boundary <= bc[d] && bc[d] <= Repeat_Boundary);
      } // d
      if (echo > 3) std::printf("# n_blocks %s\n", str(n_blocks));


      auto const n_all_blocks = size_t(n_blocks[Z])*size_t(n_blocks[Y])*size_t(n_blocks[X]);

      // regroup effective potential into blocks of 4x4x4
      assert(1 == Noco || 2 == Noco);
      p.noncollinear_spin = (2 == Noco);

      // Cartesian cell parameters for the unit cell in which the potential is defined
      double const cell[] = {ng[X]*hg[X], ng[Y]*hg[Y], ng[Z]*hg[Z]};
      double const cell_volume = cell[X]*cell[Y]*cell[Z];
      double const average_grid_spacing = std::cbrt(std::abs(hg[X]*hg[Y]*hg[Z]));
      if (echo > 1) {
          std::printf("\n# Cell summary:\n");
          for (int d = 0; d < 3; ++d) {
              std::printf("# %7d %c-points, %6d blocks, spacing= %8.6f, cell.%c= %8.3f %s, boundary= %d\n",
                  ng[d], 'x' + d, n_blocks[d], hg[d]*Ang, 'x'+ d, cell[d]*Ang, _Ang, bc[d]);
          } // d
          std::printf("# ================ ============== ================== ====================== ============\n"
              "# %7.3f M points, %11.3f k, average= %8.6f, volume= %8.1f %s^3\n\n",
              ng[X]*1e-6*ng[Y]*ng[Z], n_all_blocks*1e-3, average_grid_spacing*Ang, cell_volume*pow3(Ang), _Ang);
      } // echo

      // we assume that the source blocks lie compact in space and preferably close to each other
      std::vector<uint16_t> owner_rank(0);
      p.global_source_indices = get_right_hand_sides(n_blocks, owner_rank, echo);
      // now owner_rank[] tells the MPI rank of the process responsible for a RHS block
      uint32_t const nrhs = p.global_source_indices.size();
      if (echo > 0) std::printf("# total number of source blocks is %d\n", nrhs);

      view2D<int32_t> global_source_coords(nrhs, 4, 0);
      p.nCols = nrhs;
      double center_of_mass_RHS[] = {0, 0, 0};
      double center_of_RHSs[]     = {0, 0, 0};
      // determine the largest and smallest indices of target blocks
      // given a max distance r_trunc between source blocks and target blocks
      int32_t constexpr MinMaxLim = 1 << 21; // global coordinates with int64_t can only hold up to 3*21 bits
      int32_t max_global_source_coords[] = {-MinMaxLim, -MinMaxLim, -MinMaxLim};
      int32_t min_global_source_coords[] = { MinMaxLim,  MinMaxLim,  MinMaxLim};
      int32_t global_internal_offset[]   = {0, 0, 0};
      double max_distance_from_comass{0}, max_distance_from_center{0};
      { // scope: determine min, max, center
          double const by_nrhs = 1./std::max(1., double(nrhs));
          for (uint32_t irhs = 0; irhs < nrhs; ++irhs) {
              global_coordinates::get(global_source_coords[irhs], p.global_source_indices[irhs]);
              for (int d = 0; d < 3; ++d) {
                  auto const rhs_coord = global_source_coords(irhs,d);
                  center_of_mass_RHS[d] += (rhs_coord*4 + 1.5)*hg[d]*by_nrhs;
                  min_global_source_coords[d] = std::min(min_global_source_coords[d], rhs_coord);
                  max_global_source_coords[d] = std::max(max_global_source_coords[d], rhs_coord);
              } // d
          } // irhs
          if (echo > 0) std::printf("# all sources within (%s) and (%s)\n",
              str(min_global_source_coords), str(max_global_source_coords));

          for (int d = 0; d < 3; ++d) {
              auto const middle2 = min_global_source_coords[d] + max_global_source_coords[d];
              global_internal_offset[d] = middle2/2; // integer division by 2, "lesser half" 
              center_of_RHSs[d] = ((middle2*0.5)*4 + 1.5)*hg[d];
          } // d

          if (echo > 0) std::printf("# internal and global coordinates differ by %s\n", str(global_internal_offset));

          p.source_coords = get_memory<int16_t[4]>(nrhs, echo, "source_coords"); // internal coordinates
          { // scope: fill p.source_coords and compute the largest distance from the center or center of mass
              double max_d2m{0}, max_d2c{0};
              for (uint32_t irhs = 0; irhs < nrhs; ++irhs) {
                  double d2m{0}, d2c{0};
                  for (int d = 0; d < 3; ++d) {
                      auto const source_coord = global_source_coords(irhs,d) - global_internal_offset[d];
                      auto const src_coord_16 = int16_t(source_coord); // convert to shorter integer type
                      assert(source_coord == src_coord_16 && "internal source_coords use int16_t, maybe too large");
                      p.source_coords[irhs][d] = src_coord_16;
                      auto const cube_center = (global_source_coords(irhs,d)*4 + 1.5)*hg[d];
                      d2m += pow2(cube_center - center_of_mass_RHS[d]);
                      d2c += pow2(cube_center - center_of_RHSs[d]);
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
      auto const r_trunc = control::get("green_function.truncation.radius", 10.);
      if (echo > 0) std::printf("# green_function.truncation.radius=%g %s, %.1f grid points\n", r_trunc*Ang, _Ang, r_trunc/average_grid_spacing);
      p.r_truncation  = std::max(0., r_trunc);
      // confinement potential
      p.r_confinement = std::min(std::max(0., r_trunc - 2.0), p.r_truncation);
      p.V_confinement = control::get("green_function.confinement.potential", 1.);
      if (echo > 0) std::printf("# confinement potential %g*(r/Bohr - %g)^4 %s\n", p.V_confinement*eV, p.r_confinement, _eV);
      if (echo > 0) std::printf("# V_confinement(r_truncation)= %g %s\n", p.V_confinement*eV*pow4(r_trunc - p.r_confinement), _eV);

      // count the number of green function elements for each target block

      double r_block_circumscribing_sphere{0};

      uint32_t num_target_coords[3] = {0, 0, 0}; // range of target coordinates
      int32_t  min_target_coords[3] = {0, 0, 0}; // minimum of global target coordinates
      int32_t  max_target_coords[3] = {0, 0, 0}; // maximum of global target coordinates
      { // scope: create the truncated Green function block-sparsity pattern
          auto const rtrunc = std::max(0., r_trunc);
          double scale_grid_spacing[] = {1, 1, 1};
          {
              double const def = control::get(scale_grid_spacing, "green_function.scale.grid.spacing", "xyz", 1.);
              if (1. != def) warn("using +green_function.scale.grid.spacing=%g without .x, .y or .z may be confusing", def);
          }

          double h[] = {hg[X], hg[Y], hg[Z]}; // customize grid spacing for the truncation sphere, also used in green_potential::multiply
          for (int d = 0; d < 3; ++d) { // spatial directions

              if (r_trunc >= 0) {
                  // char keyword[64]; std::snprintf(keyword, 64, "green_function.scale.grid.spacing.%c", 'x' + d);
                  // auto const scale_h = control::get(keyword, 1.);
                  // this feature allows also truncation ellipsoids (e.g. by .x != .y == .z)
                  // and cyliders (e.g. .x=1, .y=1, .z=0) and planes (e.g. .x=0, .y=0, .z=1)
                  auto const scale_h = scale_grid_spacing[d];
                  if (scale_h >= 0) {
                      h[d] = hg[d]*scale_h;
                      if (1. != scale_h) {
                          if (echo > 1) std::printf("# scale grid spacing in %c-direction for truncation from %g to %g %s\n", 'x'+d, hg[d]*Ang, h[d]*Ang, _Ang);
                      } // scale factor deviates from unity
                  } else error("cannot use +green_function.scale.grid.spacing.%c=%g with a negative number", 'x'+d, scale_h);

                  if (Periodic_Boundary == bc[d]) { // periodic boundary conditions
                      auto const deformed_cell = h[d]*ng[d];
                      if (2*rtrunc > deformed_cell) {
                          if (h[d] > 0) {
                              warn("truncation sphere (diameter= %g %s) does not fit cell in %c-direction (%g %s)" "\n#               "
                                   "better use +green_function.scale.grid.spacing.%c=0 for cylindrical truncation",
                                   2*rtrunc*Ang, _Ang, 'x' + d, deformed_cell*Ang, _Ang, 'x' + d);
                              h[d] = 0; // cylindrical or plane or no truncation (depending on the other boundaries)
                              assert(deformed_cell > 0);
                          } else { // h[d] > 0
                              assert(0 == h[d]); // truncation in this direction has been switched off manually --> no warning
                          }
                      } else {
                          bc[d] = Wrap_Boundary; // now bc[d] may differ from boundary_condition[d]
                          // truncation sphere does not overlap with its periodic images, we can treat it like a Vacuum_Boundary
                          if (echo > 1) std::printf("# periodic boundary condition in %c-direction is wrapped\n", 'x' + d);
                      }
                  } // periodic boundary condition
              } // r_trunc >= 0

          } // d

          r_block_circumscribing_sphere = 0.5*(4 - 1)*std::sqrt(pow2(h[X]) + pow2(h[Y]) + pow2(h[Z]));
          if (echo > 0) std::printf("# circumscribing radius= %g %s\n", r_block_circumscribing_sphere*Ang, _Ang);
          auto const rtrunc_plus  =              rtrunc + 2*r_block_circumscribing_sphere;
          auto const rtrunc_minus = std::max(0., rtrunc - 2*r_block_circumscribing_sphere);
          if (echo > 0) std::printf("# truncation radius %g %s, search within %g %s\n", rtrunc*Ang, _Ang, rtrunc_plus*Ang, _Ang);
          if (echo > 0 && rtrunc_minus > 0) std::printf("# blocks with center distance below %g %s are fully inside\n", rtrunc_minus*Ang, _Ang);

          int32_t itr[3];
          for (int d = 0; d < 3; ++d) { // spatial directions

              // how many blocks around each source block do we need to check
              itr[d] = (h[d] > 0) ? std::floor(rtrunc_plus/(4*h[d])) : (n_blocks[d]/2);
              assert(itr[d] >= 0);

              min_target_coords[d] = min_global_source_coords[d] - itr[d];
              max_target_coords[d] = max_global_source_coords[d] + itr[d];

              if (Isolated_Boundary == bc[d]) {
                  int32_t const nb = n_blocks[d];
                  // limit to global coordinates in [0, nb)
                  min_target_coords[d] = std::max(min_target_coords[d], 0);
                  max_target_coords[d] = std::min(max_target_coords[d], nb - 1);
              } // Isolated_Boundary

              num_target_coords[d] = std::max(0, max_target_coords[d] + 1 - min_target_coords[d]);

          } // d
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
                       r2trunc_plus   = pow2(rtrunc_plus);
#ifndef   USE_SIMPLE_RANGE_TRUNCATION
          double const r2trunc_minus  = pow2(rtrunc_minus),
                       r2block_circum = pow2(r_block_circumscribing_sphere*3); // Maybe this can be reduced to *2
#endif // USE_SIMPLE_RANGE_TRUNCATION

          std::vector<int32_t> tag_diagonal(product_target_blocks, -1);
          assert(nrhs < (1ul << 16) && "the integer type of ColIndex is uint16_t!");
          std::vector<std::vector<bool>> sparsity_pattern(nrhs); // std::vector<bool> is a memory-saving bit-array

          simple_stats::Stats<> inout[4]; // 4 classes {inside, partial, outside, checked}
          for (uint32_t irhs = 0; irhs < nrhs; ++irhs) { // must be a serial loop if the order in column_indices is relevant
              auto & sparsity_RHS = sparsity_pattern[irhs]; // abbreviate
              sparsity_RHS.resize(product_target_blocks, false);
              auto const *const source_coords = global_source_coords[irhs]; // global source block coordinates
              simple_stats::Stats<> stats[3];
              int constexpr max_nci = 27; // nci == number_of corners inside
              std::vector<uint32_t> hist(1 + max_nci, 0); // distribution of nci
              std::vector<simple_stats::Stats<>> stats_d2(1 + max_nci);
              size_t hit_single{0}, hit_multiple{0};
              int64_t idx3_diagonal{-1};

              int32_t b_first[3], b_last[3]; // box extent relative to source block
              for (int d = 0; d < 3; ++d) {
                  if (Isolated_Boundary == bc[d]) {
                      b_first[d] = std::max(-itr[d], min_target_coords[d] - source_coords[d]);
                      b_last[d]  = std::min( itr[d], max_target_coords[d] - source_coords[d]);
                  } else
                  if (Periodic_Boundary == bc[d]) {
                      b_first[d] = (1 - int32_t(n_blocks[d]))/2;
                      b_last[d]  =      int32_t(n_blocks[d]) /2;
                  } else { // boundary_condition
                      assert(Wrap_Boundary == bc[d] || Repeat_Boundary == bc[d] || Vacuum_Boundary == bc[d]);
                      b_first[d] = -itr[d]; 
                      b_last[d]  =  itr[d];
                  } // boundary_condition
              } // d
              if (echo > 7) std::printf("# RHS#%i checks target box from (%s) to (%s)\n", irhs, str(b_first), str(b_last));

              int32_t target_coords[3]; // global target block coordinates
              for (int32_t bz = b_first[Z]; bz <= b_last[Z]; ++bz) { target_coords[Z] = source_coords[Z] + bz;
                  assert(target_coords[Z] >= min_target_coords[Z] && target_coords[Z] <= max_target_coords[Z]);
              for (int32_t by = b_first[Y]; by <= b_last[Y]; ++by) { target_coords[Y] = source_coords[Y] + by;
                  assert(target_coords[Y] >= min_target_coords[Y] && target_coords[Y] <= max_target_coords[Y]);
              for (int32_t bx = b_first[X]; bx <= b_last[X]; ++bx) { target_coords[X] = source_coords[X] + bx;
                  assert(target_coords[X] >= min_target_coords[X] && target_coords[X] <= max_target_coords[X]);

//                for (int d = 0; d < 3; ++d ) {
//                    assert(target_coords[d] >= min_target_coords[d] && target_coords[d] <= max_target_coords[d]);
//                } // d

                  // d2 is the distance^2 of the block centers
                  auto const d2 = pow2(bx*4*h[X]) + pow2(by*4*h[Y]) + pow2(bz*4*h[Z]);

                  int nci{0}; // init number of corners inside
                  if (d2 < r2trunc_plus) { // potentially inside, check all 8 or 27 corner cases
#ifdef    USE_SIMPLE_RANGE_TRUNCATION
                      if (d2 <= r2trunc) { nci = max_nci; } // similar to tfqmrgpu_generate_FD_example.cxx
                                                            // skip the 8- or 27-corners test for inner blocks
#else  // USE_SIMPLE_RANGE_TRUNCATION
                      int const far = (d2 > r2block_circum); // far in {0, 1}
                      // i = i4 - j4 --> i in [-3, 3],
                      //     if two blocks are far from each other, we test only the 8 combinations of |{-3, 3}|^3
                      //     for blocks close to each other, we test all 27 combinations of |{-3, 0, 3}|^3
                      int const inc = 3 + 3*far;
                      for (int iz = -3; iz <= 3; iz += inc) { auto const d2z   = pow2((bz*4 + iz)*h[Z]);
                      for (int iy = -3; iy <= 3; iy += inc) { auto const d2yz  = pow2((by*4 + iy)*h[Y]) + d2z;
                      for (int ix = -3; ix <= 3; ix += inc) { auto const d2xyz = pow2((bx*4 + ix)*h[X]) + d2yz;
#if 0
                          if (0 == irhs && (d2xyz < r2trunc) && echo > 17) {
                              std::printf("# %s: b= %i %i %i, i-j %i %i %i, d^2= %g %s\n",
                                  __func__, bx,by,bz, ix,iy,iz, d2xyz, (d2xyz < r2trunc)?"in":"out");
                          }
#endif // 0
                          nci += (d2xyz < r2trunc); // add 1 if inside
                      }}} // ix iy iz
                      int const mci = far ? 8 : 27;
                      if (d2 < r2trunc_minus) assert(mci == nci); // for these, we could skip the 8-corners test
                      nci = (nci*27)/mci; // limit nci to [0, 27], also the 9 different far cases are cast into the bin numbers {0,3,6,10,13,16,20,23,27}
#endif // USE_SIMPLE_RANGE_TRUNCATION
                  } // d2 < r2trunc_plus

                  if (nci > 0) {
                      // at least one grid point in target block (bx,by,bz) is closer than rtrunc to a grid point in the source block
                      int32_t idx[3]; // internal coordinates
                      for (int d = 0; d < 3; ++d) {
                          idx[d] = target_coords[d] - min_target_coords[d];
                          assert(0 <= idx[d]); assert(idx[d] < num_target_coords[d]);
                      } // d
                      auto const idx3 = index3D(num_target_coords, idx); // idx3 is a flat index into the target box
                      assert(idx3 < product_target_blocks);
                      if (sparsity_RHS[idx3]) {
                          ++hit_multiple; // already set
                      } else {
                          sparsity_RHS[idx3] = true;
                          column_indices[idx3].push_back(irhs);
                          for (int d = 0; d < 3; ++d) {
                              stats[d].add(target_coords[d]);
                          } // d
                          ++hit_single;
                      } // has not been hit yet
                      if (0 == bz && 0 == by && 0 == bx) {
                          assert(-1 == idx3_diagonal);
                          assert(-1 == tag_diagonal[idx3]);
                          tag_diagonal[idx3] = irhs;
                          idx3_diagonal = idx3;
                      } // diagonal entry

                      if (nci < max_nci) { // partial hit
                          // insert these lines into green_function.svg to visualize the partially hit target blocks (2D)
//                        std::printf("  <rect width=\"4\" height=\"4\" x=\"%d\" y=\"%d\" fill=\"none\" stroke=\"blue\" />\n", target_coords[X]*4, target_coords[Y]*4); // SVG
                      } // partial

                  } // nci > 0
                  ++hist[nci];
                  stats_d2[nci].add(d2);

              }}} // bx // by // bz
              if (echo > 8) std::printf("# RHS#%i has %ld single and %ld multiple hits\n", irhs, hit_single, hit_multiple);
              assert(0 == hit_multiple); // in this version, target blocks may not be hit more than once
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
                  if (echo > 8) std::printf("# RHS#%i has %.3f k inside, %.3f k partial and %.3f k outside (of %.3f k checked blocks)\n",
                                irhs, hist[max_nci]*.001, partial*.001, hist[0]*.001, total_checked*.001);
                  inout[0].add(hist[max_nci]);  // inside
                  inout[1].add(partial);        // partial
                  inout[2].add(hist[0]);        // outside
                  inout[3].add(total_checked);  // checked
              } // echo
              assert(idx3_diagonal > -1 && "difference vector (0,0,0) must be hit once");
              assert(tag_diagonal[idx3_diagonal] == irhs && "diagonal inconsistent");
          } // irhs

          if (echo > 0) {
              char const inout_class[][8] = {"inside", "partial", "outside",  "checked"};
              for (int i = 0; i < 4; ++i) {
                  std::printf("# RHSs have [%7g,%9.1f +/-%5.1f, %7g] blocks %s\n",
                      inout[i].min(), inout[i].mean(), inout[i].dev(), inout[i].max(), inout_class[i]);
              } // i
          } // echo

          // a histogram about the distribution of the number of columns per row
          std::vector<uint32_t> hist(1 + nrhs, 0);
          for (size_t idx3 = 0; idx3 < column_indices.size(); ++idx3) {
              auto const nc = column_indices[idx3].size();
              assert(nc <= nrhs);
              ++hist[nc];
          } // idx3

          // eval the histogram
          size_t nall{0}; // checksum for histogram
          size_t nnzb{0}; // number of non-zero BSR entries in X
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
                               "average %.1f per source block\n", nnzb*.001, nnzb/std::max(nrhs*1., 1.));
          if (echo > 0) std::printf("# %.3f k (%.1f %% of %.3f k) target blocks are active\n",
              p.nRows*.001, p.nRows/(product_target_blocks*.01), product_target_blocks*.001);

#ifdef    HAS_BITMAP_EXPORT
          { // scope: export_as_bitmap, reduce over z-coordinate
              int const nx = num_target_coords[X], ny = num_target_coords[Y];
              view3D<float> image(ny, nx, 4, 0.f);
              for (uint32_t irhs = 0; irhs < nrhs; ++irhs) {
                  auto const & sparsity_RHS = sparsity_pattern[irhs]; // abbreviate
                  for (size_t idx3 = 0; idx3 < product_target_blocks; ++idx3) {
                      if (sparsity_RHS[idx3]) {
                          int32_t const ix = idx3 % nx, iy = (idx3/nx) % ny;
                          int constexpr GREEN = 1;
                          image(iy,ix,GREEN) += 1;
                      }
                  } // idx3
              } // irhs
              float maxval{0};
              for (int iy = 0; iy < ny; ++iy) {
                  for (int ix = 0; ix < nx; ++ix) {
                      for (int rgba = 0; rgba < 4; ++rgba) {
                          auto & f = image(iy,ix,rgba);
                          // apply non-linear transforms here, e.g. sqrt to pronounce the smaller values
                          f = std::sqrt(f);
                          maxval = std::max(maxval, f);
                      } // rgba
                  } // ix
              } // iy
              bitmap::write_bmp_file("green_function", image.data(), ny, nx, -1, 254.999/maxval);
          } // scope: export_as_bitmap
#endif // HAS_BITMAP_EXPORT


          assert(nnzb < (1ull << 32) && "the integer type of RowStart is uint32_t!");

          // resize BSR tables: (Block-compressed Sparse Row format)
          if (echo > 3) std::printf("# memory of a complex Green function is %.6f %s (float, twice for double)\n",
                                       nnzb*2.*64.*64.*sizeof(float)*GByte, _GByte);
          p.colindx.resize(nnzb);
          p.rowindx  = get_memory<uint32_t>(nnzb, echo, "rowindx");
          p.RowStart = get_memory<uint32_t>(p.nRows + 1, echo, "RowStart");
          p.RowStart[0] = 0;
          p.veff_index = get_memory<int32_t>(nnzb, echo, "veff_index"); // indirection list for the local potential
          set(p.veff_index, nnzb, -1); // init as non-existing
          p.target_coords = get_memory<int16_t[3+1]>(p.nRows, echo, "target_coords");
          p.rowCubePos    = get_memory<float  [3+1]>(p.nRows, echo, "rowCubePos"); // internal coordinates but in float
          p.colCubePos    = get_memory<float  [3+1]>(p.nCols, echo, "colCubePos"); // internal coordinates but in float
          p.target_minus_source = get_memory<int16_t[3+1]>(nnzb, echo, "target_minus_source");

          for (unsigned iCol = 0; iCol < p.nCols; ++iCol) {
              for (int d = 0; d < 3; ++d) {
                  p.colCubePos[iCol][d] = global_source_coords(iCol,d) - global_internal_offset[d];
              } // d
              p.colCubePos[iCol][3] = 0.f; // not used
          } // iCol

          p.global_target_indices.resize(p.nRows);
          p.subset.resize(p.nCols); // we assume columns of the unit operator as right-hand-sides

          view3D<int32_t> iRow_of_coords(num_target_coords[Z],num_target_coords[Y],num_target_coords[X], -1); // init as non-existing

          { // scope: fill BSR tables
              simple_stats::Stats<> st;
              uint32_t iRow{0}; // init as 1st index
              for (uint32_t z = 0; z < num_target_coords[Z]; ++z) { // serial
              for (uint32_t y = 0; y < num_target_coords[Y]; ++y) { // serial
              for (uint32_t x = 0; x < num_target_coords[X]; ++x) { // serial
                  uint32_t const idx[] = {x, y, z};
                  auto const idx3 = index3D(num_target_coords, idx);
                  assert(idx3 < product_target_blocks);

                  auto const ncols = column_indices[idx3].size();
                  if (ncols > 0) {
                      st.add(ncols);
                      iRow_of_coords(idx[Z], idx[Y], idx[X]) = iRow; // set existing

                      p.RowStart[iRow + 1] = p.RowStart[iRow] + ncols;
                      // copy the column indices
                      set(p.colindx.data() + p.RowStart[iRow], ncols, column_indices[idx3].data());
                      // copy the target block coordinates
                      int32_t global_target_coords[3];
                      for (int d = 0; d < 3; ++d) {
                          global_target_coords[d] = idx[d] + min_target_coords[d];
                          auto const internal_target_coord = global_target_coords[d] - global_internal_offset[d];
                          p.target_coords[iRow][d] = internal_target_coord; assert(internal_target_coord == p.target_coords[iRow][d] && "safe assign");
                          p.rowCubePos[iRow][d]    = internal_target_coord; assert(internal_target_coord == p.rowCubePos[iRow][d]    && "safe assign");
                      } // d
                      p.target_coords[iRow][3] = 0; // not used
                      p.rowCubePos[iRow][3] = 0.f; // not used

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
                              { // scope: search inz such that p.colindx[inz] == iCol
                                  int64_t inz_found{-1};
                                  for (auto inz = p.RowStart[iRow]; inz < p.RowStart[iRow + 1] && -1 == inz_found; ++inz) {
                                      if (iCol == p.colindx[inz]) inz_found = inz;
                                  } // inz
                                  assert(-1 != inz_found && "iCol should be in the list");
                                  p.subset[iCol] = inz_found;
                              } // scope
                          } // iCol valid
                      } // scope: determine the diagonal entry

                      int32_t veff_index{-1};
                      if (1) { // scope: fill indirection table for having the local potential only defined in 1 unit cell and repeated periodically
                          int32_t mod[3];
                          bool potential_given{true};
                          for (int d = 0; d < 3; ++d) {
                              mod[d] = global_target_coords[d]; // global target coordinates may be negative
                              if (Periodic_Boundary == bc[d] || Wrap_Boundary == bc[d] || Repeat_Boundary == bc[d]) {
                                  mod[d] = global_target_coords[d] % n_blocks[d];
                                  mod[d] += (mod[d] < 0)*n_blocks[d]; // cast into range [0, n_blocks[d] - 1]
                              } else {
                                  assert(Vacuum_Boundary == bc[d] || Isolated_Boundary == bc[d]);
                                  potential_given = potential_given && (mod[d] >= 0 && mod[d] < n_blocks[d]);
                                  // potential element is only given if inside the unit cell
                              } // bc
                          } // d
                          if (potential_given) {
                              auto const iloc = index3D(n_blocks, mod);
                              veff_index = iloc; assert(iloc == veff_index && "safe assign");
                          } else { // potential_given
                              assert(Vacuum_Boundary == bc[X] || Vacuum_Boundary == bc[Y] || Vacuum_Boundary == bc[Z]);
                              veff_index = -1; // outside of the unit cell due to Vacuum_Boundary
                          } // potential_given
                      } // scope: fill indirection table

                      for (auto inz = p.RowStart[iRow]; inz < p.RowStart[iRow + 1]; ++inz) {
                          auto const iCol = p.colindx[inz];
                          for (int d = 0; d < 3; ++d) {
                              auto const diff = int32_t(p.target_coords[iRow][d]) - p.source_coords[iCol][d];
                              p.target_minus_source[inz][d] = diff; assert(diff == p.target_minus_source[inz][d] && "safe assign");
                          } // d
                          p.target_minus_source[inz][3] = 0; // component #3 not used
                          p.rowindx[inz] = iRow;
                          p.veff_index[inz] = veff_index;
                      } // inz

                      ++iRow; // count up the number of active rows
                  } // n > 0
              }}} // idx
              assert(p.nRows == iRow && "counting 2nd time");
              assert(nnzb == p.RowStart[p.nRows] && "sparse matrix consistency");
              if (echo > 2) std::printf("# source blocks per target block in [%g, %.1f +/- %.1f, %g]\n", st.min(), st.mean(), st.dev(), st.max());
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

          { // scope: set up kinetic plans
              auto const keyword = "green_kinetic.range";
              int16_t const kinetic_nFD_default = control::get(keyword, 8.); // if possible use 16th order Laplace operator
              for (int dd = 0; dd < 3; ++dd) { // derivate direction
                  int16_t kinetic_nFD_dd{kinetic_nFD_default}; // suggestion for this direction
                  char keyword_dd[32]; std::snprintf(keyword_dd, 32, "%s.%c", keyword, 'x' + dd);

                  // create lists for the finite-difference derivatives
                  auto const new_stat = kinetic_plan::finite_difference_plan(p.kinetic[dd].sparse_, kinetic_nFD_dd // results
                      , dd
                      , (Periodic_Boundary == bc[dd]) // derivative direction is periodic? (not wrapped)
                      , num_target_coords
                      , p.RowStart, p.colindx.data()
                      , iRow_of_coords
                      , sparsity_pattern.data()
                      , nrhs, echo);
                  if (0 == new_stat) {
                      p.kinetic[dd].set(dd, hg[dd], nnzb, echo);
                      p.kinetic[dd].FD_range_ = control::get(keyword_dd, double(kinetic_nFD_dd));
                  } else error("failed to create new kinetic_plan_t in %c-direction", 'x' + dd);   

              } // dd derivate direction
          } // scope: set up kinetic plans

          // transfer stuff into managed GPU memory

          p.grid_spacing_trunc = get_memory<double>(3, echo, "grid_spacing_trunc");
          set(p.grid_spacing_trunc, 3, h); // customized grid spacings used for the construction of the truncation sphere

          p.phase = get_memory<double[2][2]>(3, echo, "phase");
          green_kinetic::set_phase(p.phase, nullptr, echo); // init with Gamma point

      } // scope

      p.Veff = get_memory<double(*)[64]>(4, echo, "Veff");
      for (int mag = 0; mag < 4; ++mag) { p.Veff[mag] = nullptr; }

      for (int mag = 0; mag < Noco*Noco; ++mag) {
          p.Veff[mag] = get_memory<double[64]>(p.nRows, echo, "Veff[mag]"); // in managed memory
          set(p.Veff[mag][0], p.nRows*64, 0.0);
      } // mag

      warn("potential has not been filled!", __LINE__);

#ifdef    GREEN_FUNCTION_SVG_EXPORT
      { // scope
        assert(nullptr == svg);
        svg = std::fopen("green_function.svg", "w");
        if (nullptr != svg) {
            std::fprintf(svg, "<!-- SVG code generated by %s -->\n", __FILE__); // can be viewed at https://editsvgcode.com/
            std::fprintf(svg, "<svg viewBox=\"%d %d %d %d\" xmlns=\"http://www.w3.org/2000/svg\">\n",
                min_target_coords[X]*4-20, min_target_coords[Y]*4-20, num_target_coords[X]*4+20, num_target_coords[Y]*4+20);
        } // nullptr != svg
      } // scope
#endif // GREEN_FUNCTION_SVG_EXPORT

      auto const stat = green_function::construct_dyadic_plan(p.dyadic_plan
                            , cell, bc, hg, xyzZinso
                            , p.nRows, p.nCols, p.RowStart, p.colindx.data()
                            , p.target_coords, global_internal_offset, r_block_circumscribing_sphere
                            , max_distance_from_center, r_trunc, echo, Noco);
      if (stat && echo > 0) std::printf("# construct_dyadic_plan returned status= %i\n", int(stat));

      auto const nerr = p.dyadic_plan.consistency_check();
      if (nerr && echo > 0) std::printf("# dyadic_plan.consistency_check has %d errors\n", nerr);

#ifdef    GREEN_FUNCTION_SVG_EXPORT
      if (nullptr != svg) {
            // show the cell boundaries if in range
            std::fprintf(svg, "  <rect width=\"%d\" height=\"%d\" x=\"%d\" y=\"%d\" fill=\"none\" stroke=\"black\" />\n", n_blocks[X]*4, n_blocks[Y]*4, 0, 0);
            // show a square box for each source block
            for (uint32_t irhs = 0; irhs < nrhs; ++irhs) { auto const *const v = global_source_coords[irhs];
                std::fprintf(svg, "  <rect width=\"%d\" height=\"%d\" x=\"%d\" y=\"%d\" fill=\"none\" stroke=\"grey\" />\n", 4, 4, v[X]*4, v[Y]*4);
            } // irhs
            // show a 4 truncation spheres around each source block
            auto const h = p.grid_spacing_trunc;
            if (h[X] > 0 && h[Y] > 0) {
                auto const rx = r_trunc/h[X], ry = r_trunc/h[Y];
                for (uint32_t irhs = 0; irhs < nrhs; ++irhs) { auto const *const v = global_source_coords[irhs];
                    std::fprintf(svg, "  <ellipse cx=\"%g\" cy=\"%g\" rx=\"%g\" ry=\"%g\" fill=\"none\" stroke=\"green\" />\n", v[X]*4+0.5, v[Y]*4+0.5, rx, ry);
                    std::fprintf(svg, "  <ellipse cx=\"%g\" cy=\"%g\" rx=\"%g\" ry=\"%g\" fill=\"none\" stroke=\"green\" />\n", v[X]*4+3.5, v[Y]*4+0.5, rx, ry);
                    std::fprintf(svg, "  <ellipse cx=\"%g\" cy=\"%g\" rx=\"%g\" ry=\"%g\" fill=\"none\" stroke=\"green\" />\n", v[X]*4+0.5, v[Y]*4+3.5, rx, ry);
                    std::fprintf(svg, "  <ellipse cx=\"%g\" cy=\"%g\" rx=\"%g\" ry=\"%g\" fill=\"none\" stroke=\"green\" />\n", v[X]*4+3.5, v[Y]*4+3.5, rx, ry);
                } // irhs
            } // grid spacings relevant for truncation are positive
          std::fprintf(svg, "</svg>\n");
          std::fclose(svg);
          svg = nullptr;
          if (echo > 3) std::printf("# file green_function.svg (Scalable Vector Graphics) written\n");
      } // nullptr != svg
#endif // GREEN_FUNCTION_SVG_EXPORT

      return stat;
  } // construct_Green_function

  #undef str // === vec2str.c_str()











#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
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


  status_t test_Green_function(int const echo=0) {
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
  } // test_Green_function

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_Green_function(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_function
