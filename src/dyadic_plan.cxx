// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, FILE, ::fprintf
#include <cassert> // assert

#define DEBUG

#include "dyadic_plan.hxx"

#include "status.hxx" // status_t
#include "green_memory.hxx" // get_memory, free_memory, dim3, real_t_name
#include "green_sparse.hxx" // ::sparse_t<>
#include "inline_math.hxx" // pow2, pow3
#include "sho_tools.hxx" // ::nSHO, ::n2HO, ::n1HO
#include "simple_timer.hxx" // SimpleTimer
#include "control.hxx" // ::get
#include "display_units.h" // Ang, _Ang
#include "boundary_condition.hxx" // Periodic_Boundary, Repeat_B*, Wrap_B*
#include "action_plan.hxx" // ::atom_t
#include "mpi_parallel.hxx" // ::max, MPI_COMM_WORLD


#if 0
  // Suggestion: this could replace AtomPos + AtomLmax in the long run --> ToDo
  struct atom_image_t {
      double pos[3]; // atom position
      float  oneoversqrtsigma; // Gaussian spread, 1/sqrt(sigma)
      int8_t shifts[3]; // periodic image shifts in [-127, 127]   OR   uint16_t phase_index;
      int8_t lmax; // SHO basis size >= -1
  }; // atom_image_t 32 Byte
#endif // 0


    dyadic_plan_t::dyadic_plan_t(int const echo) { // default constructor
        if (echo > 0) std::printf("# construct %s\n", __func__);
        // please see construct_dyadic_plan in green_function.hxx for the construction of the dyadic_plan_t
    } // empty and default constructor


    dyadic_plan_t::~dyadic_plan_t() { // destructor
#ifdef    DEBUG
        std::printf("# destruct %s with %d right-hand-sides\n", __func__, nrhs);
#endif // DEBUG
        free_memory(AtomImageIndex);
        free_memory(AtomImageStarts);
        free_memory(AtomStarts);
        if (AtomMatrices) for (uint32_t ia = 0; ia < nAtoms; ++ia) free_memory(AtomMatrices[ia]);
        free_memory(AtomMatrices);
        free_memory(grid_spacing);
        free_memory(AtomImagePos);
        free_memory(AtomImageLmax);
        free_memory(AtomImagePhase);
        free_memory(AtomImageShift);
        free_memory(AtomLmax);
        if (sparse_SHOprj) for (int32_t irhs = 0; irhs < nrhs; ++irhs) sparse_SHOprj[irhs].~sparse_t<>();
        free_memory(sparse_SHOprj);
    } // constructor


    dyadic_plan_t & dyadic_plan_t::operator=(dyadic_plan_t && rhs) { // move assignment
        std::swap(this->AtomStarts          , rhs.AtomStarts          );
        std::swap(this->AtomLmax            , rhs.AtomLmax            );
        std::swap(this->AtomMatrices        , rhs.AtomMatrices        );
        std::swap(this->nAtoms              , rhs.nAtoms              );
        std::swap(this->AtomImageIndex      , rhs.AtomImageIndex      );
        std::swap(this->AtomImageStarts     , rhs.AtomImageStarts     );
        std::swap(this->AtomImagePos        , rhs.AtomImagePos        );
        std::swap(this->AtomImageLmax       , rhs.AtomImageLmax       );
        std::swap(this->AtomImagePhase      , rhs.AtomImagePhase      );
        std::swap(this->AtomImageShift      , rhs.AtomImageShift      );
        std::swap(this->nAtomImages         , rhs.nAtomImages         );
        std::swap(this->grid_spacing        , rhs.grid_spacing        );
        std::swap(this->nrhs                , rhs.nrhs                );
        std::swap(this->sparse_SHOprj       , rhs.sparse_SHOprj       );
        std::swap(this->sparse_SHOadd       , rhs.sparse_SHOadd       );
        std::swap(this->sparse_SHOsum       , rhs.sparse_SHOsum       );
        std::swap(this->global_atom_ids     , rhs.global_atom_ids     );
        std::swap(this->global_atom_index   , rhs.global_atom_index   );
        std::swap(this->original_atom_index , rhs.original_atom_index );
        std::swap(this->AtomMatrices_       , rhs.AtomMatrices_       );
        this->update_flop_counts();
        return *this;
    } // move assignment


    status_t dyadic_plan_t::consistency_check() const {
        status_t stat(0);
        assert(nAtomImages >= nAtoms);
        auto const rowStart = sparse_SHOsum.rowStart();
        auto const colIndex = sparse_SHOsum.colIndex();
        stat += (sparse_SHOsum.nRows() != nAtoms);
        if (!rowStart) return stat;
        for (uint32_t ia = 0; ia < nAtoms; ++ia) {
            for (auto bsr = rowStart[ia]; bsr < rowStart[ia + 1]; ++bsr) {
                auto const iai = colIndex[bsr];
                stat += (AtomImageLmax[iai] != AtomLmax[ia]);
            } // bsr
        } // ia
        return stat; // number of errors
    } // consistency_check


    void dyadic_plan_t::update_flop_counts(int const echo) {

        flop_count_SHOgen = 0;
        flop_count_SHOadd = 0;
        {
            auto const iai_of_bsr = sparse_SHOadd.colIndex();
            auto const nnz        = sparse_SHOadd.nNonzeros();
            for (size_t bsr = 0; bsr < nnz; ++bsr) {
                int const lmax = AtomImageLmax[iai_of_bsr[bsr]];
                flop_count_SHOadd += 64 * flop_count_SHOprj_SHOadd(lmax); // 64 threads per block (real version)
                flop_count_SHOgen += 12 * flop_count_Hermite_Gauss(lmax); // 12 threads per block eval the Hermite polynomials
                // this flop count does not account for masked blocks
            } // bsr
        }
//        flop_count_SHOprj = flop_count_SHOadd; // symmetric, both missing factor R1C2 Noco^2

        flop_count_SHOmul = 0;
        for (int ia = 0; ia < nAtoms; ++ia) {
            int const lmax = AtomImageLmax[ia];
            flop_count_SHOmul += pow2(sho_tools::nSHO(lmax));
        } // ia
        flop_count_SHOmul *= 2*nrhs*64; // missing factor R1C2^2 Noco^3

        flop_count_SHOsum = 0;
        if (nAtomImages > nAtoms) {
            for (int iai = 0; iai < nAtomImages; ++iai) {
                int const lmax = AtomImageLmax[iai];
                flop_count_SHOsum += sho_tools::nSHO(lmax);
            } // iai
            flop_count_SHOsum *= 2*nrhs*64; // missing factor R1C2^2 Noco^2
        } else { /* no need to run SHOsum */ }

    } // update_flop_counts


  template <typename number_t>
  std::string vec2str(number_t const vec[3], double const f=1, char const *const sep=" ") {
      // convert a vector of 3 numbers into a string, with scaling f and 2 separators
      char s[64];
      std::snprintf(s, 64, "%g%s%g%s%g", vec[0]*f, sep, vec[1]*f, sep, vec[2]*f);
      return std::string(s);
  } // vec2str
  #define str(...) vec2str(__VA_ARGS__).c_str()

    dyadic_plan_t::dyadic_plan_t( // constructor
        double const cell[3]
      , int8_t const boundary_condition[3]
      , double const grid_spacing[3]
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
      , int const echo // =0 // verbosity
      , int const Noco // =1 // 1:collinear spins, 2:Non-collinear
      , FILE* const svg // =nullptr // for exporting Scalabe Vector Graphics
  ) {
      SimpleTimer timer(__FILE__, __LINE__, __func__, echo);
      if (echo > 0) std::printf("# construct %s\n", __func__);
      auto & p = *this;

      p.nrhs = nrhs;

      int constexpr X=0, Y=1, Z=2;



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
      p.AtomImagePhase = get_memory<double[4]>(  nai, echo, "AtomImagePhase");
      p.AtomImageShift = get_memory<int8_t[3+1]>(nai, echo, "AtomImageShift");
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


      // get GPU memory for the matrices
      p.AtomMatrices = get_memory<double*>(nac, echo, "AtomMatrices");
      p.AtomLmax     = get_memory<int8_t>(nac, echo, "AtomLmax");
      p.AtomStarts   = get_memory<uint32_t>(nac + 1, echo, "AtomStarts");
      p.AtomStarts[0] = 0; // init prefetch sum
      // p.AtomSigma.resize(nac);

      size_t nc2_max{1}; // green_parallel::exchange does not work with count==0
      for (uint32_t iac = 0; iac < nac; ++iac) { // parallel loop
          auto const iaa = p.global_atom_index[iac];
          // p.AtomSigma[iac] = xyzZinso[ia*8 + 6];
          uint32_t const nc = sho_tools::nSHO(atom_numax[iaa]);
          assert(nc > 0); // the number of coefficients of a contributing atom copy must be non-zero
          p.AtomStarts[iac + 1] = p.AtomStarts[iac] + nc; // create prefetch sum
          p.AtomLmax[iac] = atom_numax[iaa];
          char name[64]; std::snprintf(name, 64, "AtomMatrices[iac=%d/iaa=%d]", iac, iaa);
          p.AtomMatrices[iac] = get_memory<double>(Noco*Noco*2*nc*nc, echo, name);
          set(p.AtomMatrices[iac], Noco*Noco*2*nc*nc, 0.0); // clear
          nc2_max = std::max(nc2_max, size_t(2*nc*nc));
      } // iac
      p.nAtoms = nac; // number of contributing atom copies

      auto const count = Noco*Noco*mpi_parallel::max(nc2_max, MPI_COMM_WORLD);
      if (echo > 2) std::printf("# MPI data exchange for atom matrices with %ld doubles = %.3f kByte, Noco= %d\n", count, count*.008, Noco);
      p.AtomMatrices_ = view2D<double>(nac, count); // get CPU memory for atomic matrices

      if (echo > 1) std::printf("# found %lu contributing atoms with %lu atom images\n", nac, nai);

      p.update_flop_counts(echo); // prepare to count the number of floating point operations

  } // construct_dyadic_plan

#undef str
