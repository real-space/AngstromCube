// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, FILE, ::fprintf
#include <cassert> // assert

// #define DEBUG

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
#include "mpi_parallel.hxx" // ::max, MPI_COMM_WORLD
#include "data_view.hxx" // view2D<T>
#include "debug_output.hxx" // here

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
        free_memory(AtomSigma);
        if (sparse_SHOprj) for (int32_t irhs = 0; irhs < nrhs; ++irhs) sparse_SHOprj[irhs].~sparse_t<>();
        free_memory(sparse_SHOprj);
    } // constructor


    dyadic_plan_t & dyadic_plan_t::operator=(dyadic_plan_t && rhs) { // move assignment
        std::swap(this->AtomStarts          , rhs.AtomStarts          );
        std::swap(this->AtomLmax            , rhs.AtomLmax            );
        std::swap(this->AtomMatrices        , rhs.AtomMatrices        );
        std::swap(this->AtomSigma           , rhs.AtomSigma           );
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
        // std::swap(this->global_atom_index   , rhs.global_atom_index   );
        // std::swap(this->original_atom_index , rhs.original_atom_index );
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


  struct atom_t {
      double pos[3]; // position
      size_t  iaa; // atom copy index
      int32_t ia; // local original atom index
      int32_t i_images; // local atom image index
      int32_t i_copies; // local atom copy index
      int16_t shifts[3]; // periodic image shifts
      int16_t copies[3]; // periodic copy shifts
  }; // atom_t


    dyadic_plan_t::dyadic_plan_t( // constructor
        uint32_t const grid_blocks[3] // number of 4*4*4 grid blocks
      , int8_t const boundary_condition[3]
      , double const grid_spacing[3]
      , std::vector<double> const & xyzZinso // [natoms*8]
      , uint32_t const nRowsGreen // number of target blocks in the Green function
      , uint32_t const nrhs
      , uint32_t const *const rowStartGreen
      , uint16_t const *const colIndexGreen
    //   , int16_t const (*internal_target_coords)[3+1]
    //   , int32_t const global_internal_offset[3]
      , float const (*rowCubePos)[3+1]
      , double const r_block_circumscribing_sphere
      , double const max_distance_from_center
      , double const r_trunc
      , int const echo // =0 // verbosity
      , int const Noco // =1 // 1:collinear spins, 2:Non-collinear
      , FILE* const svg // =nullptr // for exporting Scalabe Vector Graphics
  ) {
      here;

      SimpleTimer timer(__FILE__, __LINE__, __func__, echo);
      if (echo > 0) std::printf("# construct %s\n", __func__);
      auto & p = *this;

      auto const me = mpi_parallel::rank(MPI_COMM_WORLD);

      p.nrhs = nrhs;

      int constexpr X=0, Y=1, Z=2;

      double const cell[] = {grid_blocks[X]*4*grid_spacing[X],
                             grid_blocks[Y]*4*grid_spacing[Y],
                             grid_blocks[Z]*4*grid_spacing[Z]};

      // transfer grid spacing into managed GPU memory
      p.grid_spacing = get_memory<double>(3+1, echo, "grid_spacing");
      set(p.grid_spacing, 3, grid_spacing);
      double const r_proj = control::get("green_function.projection.radius", 6.); // in units of sigma
      p.grid_spacing[3] = r_proj; // store radius in units of sigma at which the projectors stop

      int32_t const natoms = xyzZinso.size()/8; // number of original atoms
      if (echo > 2) std::printf("\n#\n# %s for %d atoms\n#\n", __func__, natoms);

      // compute which atoms will contribute, the list of natoms atoms may contain a subset of all atoms
      double min_sigma{9e9}, max_sigma{0};
      for (int32_t ia{0}; ia < natoms; ++ia) { // loop over all original atoms, can be parallel with reduction(max:max_projection_radius)
          double const sigma = xyzZinso[ia*8 + 6];
          max_sigma = std::max(max_sigma, sigma);
          min_sigma = std::min(min_sigma, sigma);
      } // ia
      assert(min_sigma > 0);
      if (echo > 3) std::printf("# projection radii are between %g and %g %s\n", min_sigma*r_proj*Ang, max_sigma*r_proj*Ang, _Ang);
      if (echo > 6) std::printf("# cubes have a circumscribing radius of %g %s\n", r_block_circumscribing_sphere*Ang, _Ang);

      if (min_sigma*r_proj <= r_block_circumscribing_sphere) {
          warn("a small projection radius (%g %s) could fall between cube corners (diagonal %g %s), enlarge +green_function.projection.radius=%.1f",
              min_sigma*r_proj*Ang, _Ang, r_block_circumscribing_sphere*Ang, _Ang, r_block_circumscribing_sphere/std::max(.01, min_sigma) + .1);
      }

      auto const radius = r_trunc + max_distance_from_center + 2*max_sigma*r_proj + 2*r_block_circumscribing_sphere;
      if (echo > 3) std::printf("# search radius is %g %s\n", radius*Ang, _Ang);
      auto const r2block_circumscribing_sphere = pow2(r_block_circumscribing_sphere);

      int iimages[3] = {0, 0, 0}; // number of images replications of the unit cell to each side (for Periodic_Boundary conditions)
      int icopies[3] = {0, 0, 0}; // number of copied atoms .................... (for Wrap_Boundary and Repeat_Boundary conditions)
      uint32_t nimages{1}, ncopies{1}; // init box products
      for (int d{0}; d < 3; ++d) { // parallel
          // what happens if there is a wrap_boundary but the non-local projection overlaps with two periodic ends of the sphere?
          // we may have to diffentiate between images and copies of atoms! Repeat_Boundary needs copies, Wrap_boundary needs copies
          int const nmx = std::max(0., std::ceil(radius/cell[d]));
          if (Periodic_Boundary == boundary_condition[d]) {
              iimages[d] = nmx;
          } else if (Repeat_Boundary == boundary_condition[d]) {
              icopies[d] = nmx;
          } else if (Wrap_Boundary == boundary_condition[d]) {
              icopies[d] = 1; // in WRAP the truncation sphere fits into the cell, ...
              // however, left targets of a left source and right targets of a right source could potentially see the same atom
          }
          if (iimages[d] > 0) assert(0 == icopies[d]);
          if (icopies[d] > 0) assert(0 == iimages[d]);
          nimages *= uint32_t(2*iimages[d] + 1); // iimage  images to the left and right
          ncopies *= uint32_t(2*icopies[d] + 1); // icopies copies to the left and right
      } // d
      size_t const nAtomCopies = natoms*ncopies;
      if (echo > 3) std::printf("# copy %d atoms %s times (%ld times) for %ld atom copies\n", natoms, str(icopies), ncopies, nAtomCopies);
      assert(ncopies >= 1); // later we divide by ncopies to retrieve the original atom index
      size_t const nAtomImages = nAtomCopies*nimages;
      if (echo > 3) std::printf("# replicate %s atom images, %ld images total\n", str(iimages), nimages);

      std::vector<uint32_t> AtomImageStarts(nAtomImages + 1, 0); // probably larger than needed, should call resize(nai + 1) later
      std::vector<atom_t> atom_data(nAtomImages);
      std::vector<bool> image_is_contribution(nAtomImages, false);
      std::vector<bool> copy_is_contributing(nAtomCopies, false);
      std::vector<bool> atom_is_contributing(natoms, false);

      simple_stats::Stats<> nc_stats;
      double sparse{0}, dense{0}; // stats to assess how much memory can be saved using sparse storage

      std::vector<std::vector<uint32_t>> cubes; // stores the row indices of Green function rows
      cubes.reserve(nAtomImages); // maximum (needs 24 Byte per atom image)

      here;

      size_t nci_stats[65]; set(nci_stats, 65, size_t(0));
      size_t far_outside{0};
      size_t iai{0}; // counter for relevant atomic images
{   SimpleTimer timer(strip_path(__FILE__), __LINE__, "computing distances with all atoms", echo);

      size_t i_images{0};
      for (int zi = -iimages[Z]; zi <= iimages[Z]; ++zi) { // serial
      for (int yi = -iimages[Y]; yi <= iimages[Y]; ++yi) { // serial
      for (int xi = -iimages[X]; xi <= iimages[X]; ++xi) { // serial
//            if (echo > 3) std::printf("# periodic shifts  %d %d %d\n", xi, yi, zi);
          int const xyz_image_shift[] = {xi, yi, zi};
          size_t iaa{0};
          for (int ia{0}; ia < natoms; ++ia) { // loop over original atoms in the unit cell, serial
              auto const gid = int32_t(xyzZinso[ia*8 + 4]); // global_atom_id
              auto const numax =   int(xyzZinso[ia*8 + 5]);
              auto const sigma =       xyzZinso[ia*8 + 6] ;

              double const r_projection = r_proj*sigma; // atom-dependent, precision dependent, assume float here
              double const r2projection = pow2(r_projection);
              double const r2projection_plus = pow2(r_projection + r_block_circumscribing_sphere);

              size_t i_copies{0};
          for (int zc = -icopies[Z]; zc <= icopies[Z]; ++zc) { // serial
          for (int yc = -icopies[Y]; yc <= icopies[Y]; ++yc) { // serial
          for (int xc = -icopies[X]; xc <= icopies[X]; ++xc) { // serial
              int const xyz_copy_shift[] = {xc, yc, zc};

              double atom_pos[3]; // suggest a shifted atomic image position
              for (int d{0}; d < 3; ++d) { // unroll
                  atom_pos[d] = xyzZinso[ia*8 + d] + (xyz_image_shift[d] + xyz_copy_shift[d])*cell[d];
              } // d
//            if (echo > 5) std::printf("# image of atom #%i at %s %s\n", gid, str(atom_pos, Ang), _Ang);

              // check all target blocks if they are inside the projection radius
              uint32_t ntb{0}; // number of hit target blocks for this image
              for (uint32_t icube{0}; icube < nRowsGreen; ++icube) { // loop over blocks
                  auto const *const target_block_coords = rowCubePos[icube];
                  // do we need to do precise checking?
                  double center_distance2{0};
                  for (int d{0}; d < 3; ++d) {
                      double const cube_center = (target_block_coords[d]*4.f + 2.0)*grid_spacing[d];
                      center_distance2 += pow2(cube_center - atom_pos[d]);
                  } // d
                  if (center_distance2 <= r2projection_plus) { // do more precise checking
//                    if (echo > 9) std::printf("# target block #%i at %s gets corner check with image at %s Bohr, radius= %g Bohr\n", icube, str(target_block_coords), str(atom_pos), r_projection);
                      int nci{0}; // number of corners inside
                      { // scope: check 8 corners
                        double d2xyz[3][2]; // squares of difference coordinates with the corners
                        for (int d{0}; d < 3; ++d) {
                            for (int ii{0}; ii < 2; ++ii) { // ii=0:  leftmost grid point (pos 0.5h)
                                                              // ii=1: rightmost grid point (pos 3.5h)
                                double const grid_point = (target_block_coords[d]*4.0 + ii*3 + 0.5)*grid_spacing[d];
                                d2xyz[d][ii] = pow2(grid_point - atom_pos[d]);
                            } // ii
                        } // d
                        for (int i8{0}; i8 < 8; ++i8) { // parallel over 8 corners, reduction on nci
                            auto const d2i = d2xyz[X][i8 & 1] + d2xyz[Y][(i8 & 2) >> 1] + d2xyz[Z][i8 >> 2];
                            if (d2i < r2projection) {
                                ++nci; // at least one corner of the block is inside the projection radius of this atom
                            } // inside the projection radius
                        } // i8
                      } // scope
                      // three different cases: 0, 1...7, 8, i.e. none, partial, full
                      if (nci > 0) {
                          // atom image contributes
                          if (0 == ntb) {
                              // this is the 1st cube to contribute
                              std::vector<uint32_t> cube_list(0); // create an empty vector
                              cubes.push_back(cube_list); // enlist the vector
                          } // 0 == ntb

                          cubes.at(iai).push_back(icube); // enlist
                          assert(cubes.at(iai).at(ntb) == icube); // check if enlisting worked

                          ++ntb;
                          assert(cubes[iai].size() == ntb);
                          int const nrhs_icube = rowStartGreen[icube + 1] - rowStartGreen[icube];
                          sparse += nrhs_icube;
                          dense  += p.nrhs; // all columns

//                        if (echo > 7) std::printf("# target block #%i at %s is inside\n", icube, str(target_block_coords));
                      } else { // nci
//                        if (echo > 9) std::printf("# target block #%i at %s is outside\n", icube, str(target_block_coords));
                          // assert(0 == nci);
                          assert(center_distance2 > r2block_circumscribing_sphere && "enlarge +green_function.projection.radius"); // an projection sphere fell between the corners
                      } // nci
                      ++nci_stats[nci];

                  } else { // d2 < r2projection_plus
                      ++far_outside;
//                    if (echo > 21) std::printf("# target block #%i at %s is far outside\n", icube, str(target_block_coords));
                  } // d2 < r2projection_plus
              } // icube

              if (ntb > 0) {
                  // atom image contributes, mark in the list to have more than 0 coefficients
                  auto const nc = sho_tools::nSHO(numax); // number of coefficients for this atom
                  nc_stats.add(nc);

                  atom_is_contributing.at(ia) = true;
                  image_is_contribution.at(iai) = true;
                  copy_is_contributing.at(iaa) = true;

                  // at least one target block has an intersection with the projection sphere of this atom image
                  auto & atom = atom_data[iai];
                  set(atom.pos, 3, atom_pos);
                  atom.iaa = iaa;
                  atom.ia = ia;
                  atom.i_images = i_images;
                  set(atom.shifts, 3, xyz_image_shift);
                  atom.i_copies = i_copies;
                  set(atom.copies, 3, xyz_copy_shift);

                  if (echo > 45) std::printf("# image of atom #%i at %s %s contributes to %d target blocks\n", gid, str(atom_pos, Ang), _Ang, ntb);
                  AtomImageStarts[iai + 1] = AtomImageStarts[iai] + nc;
                  ++iai;

                  if (nullptr != svg) { // these circles show the projection spheres of the atoms
                      double const h[] = {grid_spacing[X], grid_spacing[Y]};
                      std::fprintf(svg, "  <ellipse cx=\"%g\" cy=\"%g\" rx=\"%g\" ry=\"%g\" fill=\"none\" stroke=\"red\" />\n",
                                                    atom_pos[X]/h[X], atom_pos[Y]/h[Y], r_projection/h[X], r_projection/h[Y]);
                      std::fprintf(svg, "  <ellipse cx=\"%g\" cy=\"%g\" rx=\".1\" ry=\".1\" fill=\"none\" stroke=\"red\" />\n",
                                                    atom_pos[X]/h[X], atom_pos[Y]/h[Y]); // small circle for nucleus
                      auto const Z = xyzZinso[ia*8 + 3]; // atomic number
                      if (Z > 0) std::fprintf(svg, "  <text x=\"%g\" y=\"%g\" style=\"font-size: 3;\">%g</text>\n",
                                                    atom_pos[X]/h[X], atom_pos[Y]/h[Y], Z); // label
                  } // nullptr != svg

              } else {
//                if (echo > 15) std::printf("# image of atom #%i at %s %s does not contribute\n", atom_id, str(atom_pos, Ang), _Ang);
              } // ntb > 0
              ++i_copies;
      }}} // xc // yc // zc  --- copies
              assert(ncopies == i_copies);
              ++iaa;
          } // ia --- atoms
          assert(nAtomCopies == iaa);
          ++i_images;
      }}} // x // y // z --- images
      assert(nimages == i_images);
} // timer

      if (echo > 3) {   auto const *const s = nci_stats;
          std::printf("# nci_stats %.3f k  %ld %ld %ld %ld %ld %ld %ld  %ld\n", s[0]*.001, s[1],s[2],s[3],s[4],s[5],s[6],s[7], s[8]);
          std::printf("# %.6f M cube-atom-image pairs do not require corner checking\n", far_outside*1e-6);
      } // echo

      here;

      auto const nai = iai; // corrected number of atomic images
      if (echo > 3) std::printf("# %ld of %lu (%.2f %%) atom images have an overlap with projection spheres\n",
                                    nai, nAtomImages, nai/std::max(nAtomImages*.01, .01));
      auto const napc = AtomImageStarts[nai]; // number of atomic projection coefficients

      if (echo > 3) std::printf("# sparse %g (%.2f %%) of dense %g (100 %%)\n", sparse, sparse/(std::max(1., dense)*.01), dense);

      if (nc_stats.num() > 0 && echo > 3) std::printf("# number of coefficients per image in %s\n", nc_stats.interval().c_str());

      if (echo > 3) std::printf("# %.3f k atomic projection coefficients, %.2f per atomic image\n", napc*1e-3, napc/std::max(1., 1.*nai));
      // projection coefficients for the non-local PAW operations are stored
      // as real_t apc[napc*nrhs][R1C2][Noco][Noco*64] on the GPU
      if (echo > 3) std::printf("# memory of atomic projection coefficients is %.6f (float) and %.6f %s (double)\n",
          napc*nrhs*2.*pow2(Noco)*64.*sizeof(float)*GByte, napc*nrhs*2.*pow2(Noco)*64.*sizeof(double)*GByte, _GByte);
      p.nAtomImages = nai;
      assert(nai == p.nAtomImages); // verify

      here;

      // planning for the addition of sparse SHO projectors times dense coefficients operation
      auto const nnzb = rowStartGreen[nRowsGreen]; // number of non-zero blocks in the Green function
      std::vector<std::vector<uint32_t>> SHOadd(nnzb);
      // planning for the contraction of sparse Green function times sparse SHO projectors
      std::vector<std::vector<std::vector<uint32_t>>> SHOprj(nrhs);
      for (unsigned irhs{0}; irhs < nrhs; ++irhs) {
          SHOprj[irhs].resize(nai);
      } // irhs

      for (uint32_t iai{0}; iai < p.nAtomImages; ++iai) {
          for (uint32_t itb{0}; itb < cubes[iai].size(); ++itb) {
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

      here;

      p.sparse_SHOadd = green_sparse::sparse_t<>(SHOadd, false, "sparse_SHOadd", echo);
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

      here;

      p.sparse_SHOprj = get_memory<green_sparse::sparse_t<>>(nrhs, echo, "sparse_SHOprj");
      set((char*)p.sparse_SHOprj, nrhs*sizeof(green_sparse::sparse_t<>), '\0'); // clear
      size_t nops{0};
      for (uint32_t irhs{0}; irhs < nrhs; ++irhs) {
          char name[64]; std::snprintf(name, 64, "sparse_SHOprj[irhs=%i of %d]", irhs, nrhs);
          p.sparse_SHOprj[irhs] = green_sparse::sparse_t<>(SHOprj[irhs], false, name, echo - 4);
          // sparse_SHOprj: rows == atom images, cols == Green function non-zero elements
          nops += p.sparse_SHOprj[irhs].nNonzeros();
      } // irhs
      if (echo > 5) std::printf("# sparse_SHOprj has %d*%ld rows, %d columns and %ld nonzeros\n", nrhs,nai, nnzb, nops);
      SHOprj.resize(0); // release host memory

      here;

      p.AtomImageStarts = get_memory<uint32_t>(nai + 1, echo, "AtomImageStarts");
      set(p.AtomImageStarts, nai + 1, AtomImageStarts.data()); // copy into GPU memory
      p.AtomImageIndex = get_memory<uint32_t>(nai, echo, "AtomImageIndex");

      here;

      // get all info for the atomic matrices
      std::vector<int32_t> iac_of_iaa(nAtomCopies, -1);
      p.global_atom_ids.resize(nAtomCopies, -1);

      size_t iac{0};
      for (size_t iaa{0}; iaa < nAtomCopies; ++iaa) { // serial loop over all atoms
          if (copy_is_contributing[iaa]) {
              iac_of_iaa[iaa] = iac;
              auto const ia = iaa/ncopies;
              auto const gid = int32_t(xyzZinso[ia*8 + 4]); // global_atom_id
              p.global_atom_ids[iac] = gid;
              ++iac; // new contributing atom copy
          } // atom contributes
      } // iaa
      auto const nac = iac; // number of contributing atoms
      if (echo > 7) std::printf("# rank#%i has %ld contributing atoms\n", me, nac);
      p.global_atom_ids.resize(nac);
      // now p.global_atom_ids may have multiple entries if ncopies > 1, however, it is likely
      // that this coincides with a serial run anyway

      here;

      // store the atomic image positions in GPU memory
      p.AtomImageLmax  = get_memory<int8_t>(nai, echo, "AtomImageLmax");
      p.AtomImagePos   = get_memory<double[3+1]>(nai, echo, "AtomImagePos");
      p.AtomImagePhase = get_memory<double[4]>(  nai, echo, "AtomImagePhase");
      p.AtomImageShift = get_memory<int8_t[3+1]>(nai, echo, "AtomImageShift");
      double const phase0[] = {1, 0, 0, 0};

      std::vector<std::vector<uint32_t>> SHOsum(nac);

      here;

      for (size_t iai{0}; iai < nai; ++iai) { // serial
          auto const & atom = atom_data.at(iai);
          auto const ia = atom.ia;
          assert(0 <= ia); assert(ia < natoms);
          auto const iaa = atom.iaa;
          assert(0 <= iaa); assert(iaa < nAtomCopies);
          auto const iac = iac_of_iaa[iaa];
          assert(0 <= iac); assert(iac < nac);
          auto const gid = int32_t(xyzZinso[ia*8 + 4]); // global_atom_id
          auto const numax =   int(xyzZinso[ia*8 + 5]);
          auto const sigma =       xyzZinso[ia*8 + 6] ;
          p.AtomImageIndex[iai] = iac;
          set(p.AtomImagePos[iai], 3, atom.pos);
          set(p.AtomImageShift[iai], 3, atom.shifts); p.AtomImageShift[iai][3] = 0;
          p.AtomImagePos[iai][3] = 1./std::sqrt(sigma);
          p.AtomImageLmax[iai] = numax; // SHO basis size
          if (echo > 19) std::printf("# atom image #%ld of atom#%li has lmax= %d\n", iai, gid, p.AtomImageLmax[iai]);
          SHOsum[iac].push_back(uint32_t(iai));
          set(p.AtomImagePhase[iai], 4, phase0); // TODO construct correct phases, phase0 is just a dummy
      } // iai, copy into GPU memory

      p.sparse_SHOsum = green_sparse::sparse_t<>(SHOsum, false, "sparse_SHOsum", echo);
      SHOsum.resize(0);

      here;

      // get GPU memory for the matrices
      p.AtomMatrices = get_memory<double*>(nac, echo, "AtomMatrices");
      p.AtomLmax     = get_memory<int8_t>(nac, echo, "AtomLmax");
      p.AtomSigma    = get_memory<double>(nac, echo, "AtomSigma");
      p.AtomStarts   = get_memory<uint32_t>(nac + 1, echo, "AtomStarts");
      p.AtomStarts[0] = 0; // init prefetch sum

      size_t nc2_max{1}; // green_parallel::exchange does not work with count==0
      for (uint32_t iac{0}; iac < nac; ++iac) { // parallel loop
          auto const ia = p.global_atom_ids[iac]; // global atom index
          auto const numax = int(xyzZinso[ia*8 + 5]);
          p.AtomSigma[iac] =     xyzZinso[ia*8 + 6];
          auto const nc = sho_tools::nSHO(numax);
          assert(nc > 0); // the number of coefficients of a contributing atom copy must be non-zero
          p.AtomStarts[iac + 1] = p.AtomStarts[iac] + nc; // create prefetch sum
          p.AtomLmax[iac] = numax;
          char name[64]; std::snprintf(name, 64, "AtomMatrices[a#%i]", ia);
          p.AtomMatrices[iac] = get_memory<double>(Noco*Noco*2*nc*nc, echo, name);
          set(p.AtomMatrices[iac], Noco*Noco*2*nc*nc, 0.0); // clear GPU memory
          nc2_max = std::max(nc2_max, size_t(2*nc*nc));
      } // iac
      p.nAtoms = nac; // number of contributing atom copies

      here;

      auto const count = Noco*Noco*mpi_parallel::max(nc2_max, MPI_COMM_WORLD);
      if (echo > 2) std::printf("# MPI data exchange for atom matrices with %ld doubles = %.3f kByte, Noco= %d\n", count, count*.008, Noco);
      p.AtomMatrices_ = view2D<double>(nac, count, 0.0); // get CPU memory for atomic matrices

      here;

      if (echo > 1) std::printf("# found %lu contributing atoms with %lu atom images\n", nac, nai);

      p.update_flop_counts(echo); // prepare to count the number of floating point operations

      here;
  } // construct_dyadic_plan

#undef str
