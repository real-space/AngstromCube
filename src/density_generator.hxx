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
#include "fermi_distribution.hxx" // ::FermiLevel_t
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
      , complex_t const atom_coeff[] // projection coefficents of the wave function, flat list
      , uint32_t const coeff_starts[] // where in the flat list do the coefficients of atom #i start?, shape[natoms + 1]
      , int const natoms // number of atoms
      , double const weight=1 // weight (should be considerable, i.e. > 1e-16)
      , int const echo=0 // log-level
      , int const iband=-1 // log-info
      , int const ikpoint=-1 // log-info
  ) {
      // construct the atomic density matrix atom_rho
      for(int ia = 0; ia < natoms; ++ia) {
          int const offset = coeff_starts[ia];
          int const ncoeff = coeff_starts[ia + 1] - offset;
          auto const a_rho = atom_rho[ia];
          // add to the atomic density matrix
          for(int i = 0; i < ncoeff; ++i) {
              auto const c_i = conjugate(atom_coeff[offset + i]);
#ifdef DEVEL
              if (echo > 16) printf("# k-point #%i band #%i atom #%i coeff[%i]= %.6e %g |c|= %g\n",
                                  ikpoint, iband, ia, i, std::real(c_i), std::imag(c_i), std::abs(c_i));
#endif // DEVEL
              for(int j = 0; j < ncoeff; ++j) {
                  auto const c_j = atom_coeff[offset + j];
                  a_rho[i*ncoeff + j] += weight * std::real(c_i * c_j);
              } // j
          } // i
      } // ia
  } // add_to_density_matrices


  template <typename uint_t>
  std::vector<uint32_t> prefetch_sum(std::vector<uint_t> const & natom_coeff) {
      int const na = natom_coeff.size();
      std::vector<uint32_t> coeff_starts(na + 1);
      coeff_starts[0] = 0u;
      for(int ia = 0; ia < na; ++ia) {
          coeff_starts[ia + 1] = coeff_starts[ia] + natom_coeff[ia];
      } // ia
      return coeff_starts;
  } // prefetch_sum

  
  template <class grid_operator_t>
  view2D<typename grid_operator_t::complex_t> atom_coefficients(
        std::vector<uint32_t> & coeff_starts
      , typename grid_operator_t::complex_t const eigenfunctions[]
      , grid_operator_t const & op
      , int const nbands=1
      , int const echo=0 // log-level
      , int const ikpoint=-1 // log-info
  ) {
      // density generation for eigenstates represented on the real-space grid
      using complex_t = typename grid_operator_t::complex_t; // abbreviate
      using real_t    = decltype(std::real(complex_t(1)));

      // SimpleTimer init_function_timer(__FILE__, __LINE__, __func__, echo);
      status_t stat(0);

      assert(eigenfunctions); // may not be nullptr

      auto const natoms = op.get_natoms();
      auto const nzyx   = op.get_degrees_of_freedom();
      if (echo > 3) printf("# %s assume %d atoms and %d bands with %ld grid elements\n", __func__, natoms, nbands, nzyx);

      // the preparation of prefactors is the the same for each kpoint
      std::vector<std::vector<real_t>> scale_factors(natoms);
      std::vector<uint16_t> natom_coeff(natoms);
      for(int ia = 0; ia < natoms; ++ia) {
          int const numax = op.get_numax(ia);
          natom_coeff[ia] = sho_tools::nSHO(numax);
          auto const sigma = op.get_sigma(ia);
          if (echo > 6) printf("# %s atom #%i has numax=%d and %d coefficients, sigma= %g %s\n",
                                  __func__, ia, numax, natom_coeff[ia], sigma*Ang,_Ang);
          assert(sigma > 0);
          scale_factors[ia] = sho_projection::get_sho_prefactors<real_t>(numax, sigma);
      } // ia

      coeff_starts = prefetch_sum(natom_coeff);
      int const n_all_coeff = coeff_starts[natoms];
      int const stride = align<0>(n_all_coeff); // ToDo: test if we get the same result if we use alignment here
      if (echo > 3) printf("# %s %d atoms have %d coefficients, stride %d\n", __func__, natoms, n_all_coeff, stride);
      if (echo > 3) printf("# %s coefficients(%d,%d of %d)\n", __func__, nbands, n_all_coeff, stride);
      view2D<complex_t> coeff(nbands, stride, complex_t(0)); // get memory

      if (echo > 3) printf("# %s assume psi(%d,%ld)\n", __func__, nbands, nzyx);
      view2D<complex_t const> const psi_k(eigenfunctions, nzyx); // wrap

      data_list<complex_t> atom_coeff(natom_coeff); // get temporary memory

      { // ikpoint
          for(int iband = 0; iband < nbands; ++iband) {

              // project grid-represented eigenstates onto atomic projector functions
              stat += op.get_atom_coeffs(atom_coeff.data(), psi_k[iband], echo/2); // project

              // scale the coefficients (real-space grid only)
              for(int ia = 0; ia < natoms; ++ia) {
                  auto const c_ia = &coeff(iband,coeff_starts[ia]);
                  product(c_ia, natom_coeff[ia], atom_coeff[ia], scale_factors[ia].data());
#ifdef DEVEL
                  for(int i = 0; i < natom_coeff[ia]; ++i) {
                      auto const c_i = c_ia[i];
                      if (echo > 16) printf("# k-point #%i band #%i atom #%i coeff[%i]= %.6e %g new\n",
                                                ikpoint, iband, ia, i, std::real(c_i), std::imag(c_i));
                  } // i
#endif // DEVEL
              } // ia

          } // iband

      } // ikpoint

      if (int(stat) && echo > 3) printf("# %s get_atom_coeffs produced status sum= %i\n", __func__, int(stat));
      return coeff;
  } // atom_coefficients


  template <typename complex_t>
  status_t density(
        double *const rho              // result density on grid
      , double *const *const atom_rho  // result atomic density matrices
      , fermi_distribution::FermiLevel_t & Fermi
      , double const eigenenergies[]     // assumed shape [nbands]
      , complex_t const eigenfunctions[] // assumed shape [nbands][g.all()]
      , complex_t const atom_coeff[]     // assumed shape [nbands][n_all_coeff]
      , uint32_t const coeff_starts[]    // shape[natoms + 1]
      , int const natoms
      , real_space::grid_t const & g
      , int const nbands=1
      , double const weight_k=1
      , int const echo=0 // log-level
      , int const ikpoint=-1 // log-info
      , double *const d_rho=nullptr             // result response density
      , double *const *const d_atom_rho=nullptr // result atomic density matrices derived
      , double charges[]=nullptr // nominal charge accumulators: 0:kpoint_denominator, 1:charge, 2:d_charge
  ) {
      // density generation for eigenstates represented on the real-space grid

      // SimpleTimer init_function_timer(__FILE__, __LINE__, __func__, echo);
      status_t stat(0);

      if (nullptr == eigenfunctions) { warn("eigenfunctions received are nullptr"); return -1; }
      
      if (echo > 3) printf("# %s assume %d bands on a %d x %d x %d grid\n", __func__, nbands, g('x'), g('y'), g('z'));

      auto const n_all_coeff = coeff_starts[natoms];
      if (echo > 3) printf("# %s assume psi(%d,%ld)\n", __func__, nbands, g.all());
      view2D<complex_t const> const psi_k(eigenfunctions, g.all()); // wrap
      if (echo > 3) printf("# %s assume atom_coeff with stride %d\n", __func__, n_all_coeff);
      view2D<complex_t const> const a_coeff(atom_coeff, n_all_coeff); // wrap

      double constexpr occ_threshold = 1e-16;
      double const kT = Fermi.get_temperature(); // for display
      int const spinfactor = Fermi.get_spinfactor(); // 1 or 2
      if (echo > 3) printf("# %s spin factor %d and temperature %g %s\n", __func__, spinfactor, kT*Kelvin, _Kelvin);

      double charge[3] = {0, 0, 0}; // {kpoint weights, charge of all bands, derivative}

      { // ikpoint
          double const weight_sk = spinfactor * weight_k;

          // determine the occupation numbers
          std::vector<double> occupation(nbands, 0.0);
          std::vector<double> d_occupation(nbands, 0.0); // derivative of occupation number w.r.t. E_Fermi
          Fermi.get_occupations(occupation.data(), eigenenergies, nbands, weight_k, echo, d_occupation.data());

          int ilub{nbands}; // index of lowest completely unoccupied band
//        double old_charge{0};
          double charge_k[3] = {1, 0, 0}; // {1, charge of bands, derivative}
          for(int iband = 0; iband < nbands; ++iband) {
              auto const psi_nk = psi_k[iband];

              if (occupation[iband] >= occ_threshold) {
                  charge_k[1] += occupation[iband]*spinfactor;
                  double const weight_nk = occupation[iband] * weight_sk;
                  if (echo > 6) printf("# %s: k-point #%i bands #%i \toccupation= %.6f d_occ= %g E= %g %s\n",
                      __func__, ikpoint, iband, occupation[iband], d_occupation[iband]*kT, eigenenergies[iband]*eV, _eV);

                  add_to_density(rho, g.all(), psi_nk, weight_nk, echo, iband, ikpoint);
#if 0 // def DEBUG
                  if (echo > 0) { 
                      printf("# valence density ");
                      auto const new_charge = print_stats(rho, g.all(), g.dV());
                      printf("# valence density of band #%i added %g electrons\n", iband, new_charge - old_charge);
                      old_charge = new_charge;
                  } // echo
#endif // DEBUG
                  add_to_density_matrices(atom_rho, a_coeff[iband],
                                  coeff_starts, natoms, weight_nk, echo, iband, ikpoint);

                  if (d_occupation[iband] >= occ_threshold) {
                      charge_k[2] += d_occupation[iband]*spinfactor;
                      double const d_weight_nk = d_occupation[iband] * weight_sk;

                      if (d_rho) {
                          add_to_density(d_rho, g.all(), psi_nk, d_weight_nk, echo, iband, ikpoint);
                      } // d_rho != nullptr

                      if (d_atom_rho) {
                          add_to_density_matrices(d_atom_rho, a_coeff[iband],
                                  coeff_starts, natoms, d_weight_nk, echo, iband, ikpoint);
                      } // d_atom_rho != nullptr

                  } // derived occupations

              } else {
                  ilub = std::min(ilub, iband);
              } // weight is considerable

          } // iband
          if (ilub < nbands) {
              if (echo > 6) printf("# %s: k-point #%i band #%i at %g %s and above did not"
                  " contribute to the density\n", __func__, ikpoint, ilub, eigenenergies[ilub]*eV,_eV);
          } // some bands did not contribute to the density as they are too high above the Fermi energy
          if (echo > 5) printf("# %s: k-point #%i has charge %g electrons and derivative %g\n",
                                  __func__, ikpoint, charge_k[1], charge_k[2]*kT); 
          add_product(charge, 3, charge_k, weight_k);
      } // ikpoint
      if (charges) add_product(charges, 3, charge, 1.0);
      if (echo > 1) { printf("\n# Total valence density "); print_stats(rho, g.all(), g.dV()); }
      if (echo > 3) { printf("# Total response density"); print_stats(d_rho, g.all(), g.dV(), "", kT); }
#if 1
      if (echo > 6) {
          for(int ia = 0; ia < natoms; ++ia) {
              int const ncoeff = coeff_starts[ia + 1] - coeff_starts[ia];
              double max_rho{1e-12}; for(int ij = 0; ij < pow2(ncoeff); ++ij) max_rho = std::max(max_rho, atom_rho[ia][ij]);
              printf("\n# show %d x %d density matrix for atom #%i in %s-order, normalized to maximum %.6e\n", 
                  ncoeff, ncoeff, ia, sho_tools::SHO_order2string(sho_tools::order_zyx).c_str(), max_rho);
              for(int i = 0; i < ncoeff; ++i) {
                  printf("# izyx=%d\t", i);
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
