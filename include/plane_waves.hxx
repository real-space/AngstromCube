#pragma once

#include <cstdio> // std::printf, std::snprintf
#include <vector> // std::vector<T>

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>
#include "real_space.hxx" // ::grid_t
#include "inline_math.hxx" // align<n>

namespace plane_waves {

  class DensityIngredients {
    public:

      view2D<std::complex<double>> psi_r; // real-space representation of the wave function psi_r[nbands][ng[2]*ng[1]*ng[0]+]
      view2D<std::complex<double>> coeff; // atom projection coefficients coeff[nbands][ncoeff+]
      std::vector<double> energies; // energies[nbands]
      std::vector<uint32_t> offset; // offset[natoms + 1], marks the separations between atoms in the coefficient vectors
      double kpoint_weight;
      int kpoint_index;
      int ng[3];  // number of grid points
      int ncoeff; // total number of coefficients
      int nbands;
      int natoms;
      char tag[32];

//    DensityIngredients() {} // default constructor

      void constructor( // constructor
            int const nG[3] // grid sizes
          , int const nBands=0
          , int const nAtoms=0
          , int const nCoeff=0
          , double k_weight=1.
          , int const kpoint_id=-1
          , int const echo=0 // log-level
//    ) : kpoint_weight(k_weight), ncoeff(nCoeff), nbands(nBands), natoms(nAtoms), kpoint_index(kpoint_id) {
      ) {
          kpoint_weight = std::max(0.0, k_weight);
          ncoeff = nCoeff;
          nbands = nBands;
          natoms = nAtoms;
          kpoint_index = kpoint_id;
          set(ng, 3, nG);
          auto const nG_all = size_t(nG[2])*size_t(nG[1])*size_t(nG[0]);
          auto const nG_aligned = align<0>(nG_all);
          int  const nC_aligned = align<0>(ncoeff);
          if (echo > 8) std::printf("# DensityIngredients allocates %.3f MByte for wave functions + %.3f kByte for coefficients\n",
                                                               nbands*16e-6*nG_aligned, nbands*16e-3*nC_aligned);
          std::complex<double> const zero(0);
          psi_r = view2D<std::complex<double>>(nbands, nG_aligned, zero);
          coeff = view2D<std::complex<double>>(nbands, nC_aligned, zero);
          std::snprintf(tag, 31, "kpoint #%i weight %g", kpoint_id, kpoint_weight);
       } // constructor

  }; // DensityIngredients

  
  status_t solve(
        int const natoms_prj // number of PAW atoms
      , view2D<double const> const & xyzZ // (natoms, 4)
      , real_space::grid_t const & g // Cartesian grid descriptor for vtot
      , double const *const vtot // total effective potential on grid
      , double const *const sigma_prj=nullptr
      , int    const *const numax_prj=nullptr
      , double *const *const atom_mat=nullptr // PAW charge-deficit and Hamiltonian correction matrices
      , int const echo=0 // log-level
      , std::vector<DensityIngredients> *export_rho=nullptr        
  ); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace plane_waves
