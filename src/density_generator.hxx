#pragma once

#include <complex> // std::norm, std::real

#include "status.hxx" // status_t
#include "grid_operators.hxx" // ::grid_operator_t
#include "complex_tools.hxx" // conjugate
#include "sho_tools.hxx" // ::nSHO
#include "sho_projection.hxx" // ::sho_prefactor
#include "data_view.hxx" // view3D<T>
#include "data_list.hxx" // data_list<T>
#include "inline_math.hxx" // scale
#include "control.hxx" // ::get
#include "display_units.h" // Ang, _Ang, eV, _eV
#include "fermi_distribution.hxx" // ::Fermi_level, ::FermiLevel_t
#include "print_tools.hxx" // print_stats

namespace density_generator {

  template <typename complex_t>
  void add_to_density(
        double rho[] // is modified
      , size_t const nzyx // number of grid points , real_space::grid_t const & g // grid descriptor
      , complex_t const wave[] // wave function on real-space grid g
      , double const weight=1 // weight (should be considerable, i.e. > 1e-16)
      , int const echo=0 // log-level
      , int const iband=-1 // log-info
      , int const ikpoint=-1 // log-info
  ) {
      // add the absolute squared of a wave function to the density
      for(size_t izyx = 0; izyx < nzyx; ++izyx) { // parallel
          rho[izyx] += weight * std::norm(wave[izyx]);
      } // izyx
  } // add_to_density

  
  template <typename complex_t>
  void add_to_density_matrices(
        double *const atom_rho[] // atomic density matrices will be modified
      , complex_t const *const atom_coeff[] // projection coefficents of the wave function
      , uint16_t const natom_coeff[] // number of coefficients per atom 
      , int const natoms // number of atoms
      , double const weight=1 // weight (should be considerable, i.e. > 1e-16)
      , int const echo=0 // log-level
      , int const iband=-1 // log-info
      , int const ikpoint=-1 // log-info
  ) {
      // construct the atomic density matrix atom_rho
      for(int ia = 0; ia < natoms; ++ia) {
          int const ncoeff = natom_coeff[ia];
          // add to the atomic density matrix
          for(int i = 0; i < ncoeff; ++i) {
              auto const c_i = conjugate(atom_coeff[ia][i]);
#ifdef DEVEL
              if (echo > 6) printf("# kpoint #%i band #%i atom #%i coeff[%i]= %.6e %g\n",
                                  ikpoint, iband, ia, i, std::real(c_i), std::imag(c_i));
#endif // DEVEL
              for(int j = 0; j < ncoeff; ++j) {
                  auto const c_j = atom_coeff[ia][j];
                  atom_rho[ia][i*ncoeff + j] += weight * std::real(c_i * c_j);
              } // j
          } // i
      } // ia
  } // add_to_density




  template <class grid_operator_t>
  status_t density(
        double rho[]                   // result density on grid
      , double *const *const atom_rho  // result atomic density matrices
      , fermi_distribution::FermiLevel_t & Fermi
      , typename grid_operator_t::complex_t const eigenfunctions[]
      , double const eigenenergies[]
      , grid_operator_t const & op
      , int const nbands=1, int const nkpoints=1
      , int const echo=0
  ) {
      // density generation for eigenstates represented on the real-space grid
      using complex_t = typename grid_operator_t::complex_t; // abbreviate
      using real_t    = decltype(std::real(complex_t(1)));

      // SimpleTimer init_function_timer(__FILE__, __LINE__, __func__, echo);
      status_t stat{0};

      if (nullptr == eigenfunctions) { warn("eigenfunctions received are nullptr"); return -1; }

      auto const g = op.get_grid();
      
      if (echo > 3) printf("# %s assume %d bands and %d k-points\n", __func__, nbands, nkpoints);
      if (echo > 3) printf("# %s assume eigenfunctions on a %d x %d x %d grid\n",
                              __func__, g('x'), g('y'), g('z'));

      int const na = op.get_natoms();
      std::vector<std::vector<real_t>> scale_factors(na);
      std::vector<uint16_t> natom_coeff(na);
      for(int ia = 0; ia < na; ++ia) {
          int const numax = op.get_numax(ia);
          natom_coeff[ia] = sho_tools::nSHO(numax);
          auto const sigma = op.get_sigma(ia);
          if (echo > 6) printf("# %s atom #%i has numax=%d and %d coefficients, sigma= %g %s\n",
                                  __func__, ia, numax, natom_coeff[ia], sigma*Ang,_Ang);
          assert(sigma > 0);
          scale_factors[ia] = sho_projection::get_sho_prefactors<real_t>(numax, sigma);
      } // ia
      data_list<complex_t> atom_coeff(natom_coeff); // get memory

      view3D<complex_t const> const psi(eigenfunctions, nbands, g.all()); // wrap
      view2D<double const>    const ene(eigenenergies,  nbands); // wrap
      std::vector<double> d_rho(g.all(), 0.0); // get memory

//       double const n_electrons = control::get("valence.electrons", 0.0);
//       double const kT = control::get("electronic.temperature", 9.765625e-4);

      double constexpr occ_threshold = 1e-16;

// #pragma omp parallel
      for(int ikpoint = 0; ikpoint < nkpoints; ++ikpoint) {
          double const kpoint_weight = 1; // depends on ikpoint
          auto const psi_k = psi[ikpoint];

          // determine the occupation numbers
          std::vector<double> occupation(nbands, 0.0);
          std::vector<double> d_occupation(nbands, 0.0);
//           double Fermi_energy{-9e99};
//           stat += fermi_distribution::Fermi_level(Fermi_energy, occupation.data(),
//                              ene[ikpoint], nbands, kT, n_electrons, 2, echo);
          Fermi.get_occupations(occupation.data(), ene[ikpoint], nbands, echo + 10, d_occupation.data());

          int ilub{nbands}; // lowest_empty_band
          for(int iband = 0; iband < nbands; ++iband) {
              double const weight_nk = occupation[iband] * kpoint_weight;
              auto const psi_nk = psi_k[iband];

              if (occupation[iband] > occ_threshold) {

                  add_to_density(rho, g.all(), psi_nk, weight_nk, echo, iband, ikpoint);

                  // project eigenstates onto the atomic projector functions
                  op.get_atom_coeffs(atom_coeff.data(), psi_nk, echo/2); // project

                  // scale the coefficients (real-space grid only)
                  for(int ia = 0; ia < na; ++ia) {
                      scale(atom_coeff[ia], natom_coeff[ia], scale_factors[ia].data());
                  } // ia

                  add_to_density_matrices(atom_rho, atom_coeff.data(),
                               natom_coeff.data(), na, weight_nk, echo, iband, ikpoint);

                  if (d_occupation[iband] > occ_threshold) {
                      double const d_weight_nk = d_occupation[iband] * kpoint_weight;

                      add_to_density(d_rho.data(), g.all(), psi_nk, d_weight_nk, echo, iband, ikpoint);

                  } // derived occupations

              } else {
                  ilub = std::min(ilub, iband);
              } // weight is considerable

          } // iband
          if (ilub < nbands) {
              if (echo > 6) printf("# %s: k-point #%i bands above band #%i at %g %s did not"
                  " contribute to the density\n", __func__, ikpoint, ilub, ene(ikpoint,ilub)*eV,_eV);
          } // some bands did not contribute to the density as they are too high above the Fermi energy
      } // ikpoint

      if (echo > 1) { printf("\n# Total valence density "); print_stats(rho, g.all(), g.dV()); }
      if (echo > 3) { printf("# Total response density"); print_stats(d_rho.data(), g.all(), g.dV()); }

#ifdef DEVEL
      if (echo > 6) {
          for(int ia = 0; ia < na; ++ia) {
              int const ncoeff = sho_tools::nSHO(op.get_numax(ia));
              double max_rho{1e-12}; for(int ij = 0; ij < pow2(ncoeff); ++ij) max_rho = std::max(max_rho, atom_rho[ia][ij]);
              printf("\n# show %d x %d density matrix for atom #%i in %s-order, normalized to maximum %.6e\n", 
                  ncoeff, ncoeff, ia, sho_tools::SHO_order2string(sho_tools::order_zyx).c_str(), max_rho);
              char labels[220*8]; sho_tools::construct_label_table(labels, op.get_numax(ia), sho_tools::order_zyx);
              for(int i = 0; i < ncoeff; ++i) {
                  printf("# %s\t", &labels[i*8]);
                  for(int j = 0; j < ncoeff; ++j) {
                      printf(" %9.6f", atom_rho[ia][i*ncoeff + j]/max_rho);
                  } // j
                  printf("\n");
              } // i
              printf("\n");
          } // ia
      } // echo
#endif // DEVEL

      return stat;
  } // density
  
#ifdef NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS
  status_t all_tests(int const echo=0); // declaration only
#endif // NO_UNIT_TESTS

} // namespace density_generator
