#pragma once

#include <complex> // std::norm, std::real

#include "status.hxx" // status_t
#include "grid_operators.hxx" // ::grid_operator_t
#include "complex_tools.hxx" // conjugate
#include "sho_tools.hxx" // ::nSHO
#include "sho_projection.hxx" // ::sho_prefactor
#include "data_view.hxx" // view3D<T>
#include "data_list.hxx" // data_list<T>
#include "control.hxx" // ::get
#include "display_units.h" // Ang, _Ang

namespace density_generator {

  // ToDo: move somewhere, where potential_generator and this module can access it  
  inline double print_stats(double const values[], size_t const all, double const dV=1, char const prefix=' ') {
      double gmin{9e307}, gmax{-gmin}, gsum{0}, gsum2{0};
      for(size_t i = 0; i < all; ++i) {
          gmin = std::min(gmin, values[i]);
          gmax = std::max(gmax, values[i]);
          gsum  += values[i];
          gsum2 += pow2(values[i]);
      } // i
      printf("%c grid stats min %g max %g integral %g avg %g\n", prefix, gmin, gmax, gsum*dV, gsum/all);
      return gsum*dV;
  } // print_stats

  template<class grid_operator_t>
  status_t density(double rho[]        // result density on Cartesian grid
      , double *const *const atom_rho  // result atomic density matrices
      , typename grid_operator_t::complex_t const eigenfunctions[]
      , grid_operator_t const & op
      , int const nbands=1, int const nkpoints=1
      , int const echo=0) {
      
      using complex_t = typename grid_operator_t::complex_t; // abbreviate
      using real_t    = decltype(std::real(complex_t(1)));

      // SimpleTimer init_function_timer(__FILE__, __LINE__, __func__, echo);
      status_t stat{0};

      if (nullptr == eigenfunctions) { warn("eigenfunctions received are nullptr"); return -1; }

      auto const g = op.get_grid();
      
      if (echo > 3) printf("# %s assume %d bands and %d k-points\n", __func__, nbands, nkpoints);
      if (echo > 3) printf("# %s assume eigenfunctions on a %d x %d x %d Cartesian grid\n",
                              __func__, g('x'), g('y'), g('z'));

      int const na = op.get_natoms();
      std::vector<std::vector<real_t>> scale_factors(na);
      std::vector<int> ncoeff(na, 0);
      for(int ia = 0; ia < na; ++ia) {
          int const numax = op.get_numax(ia);
          ncoeff[ia] = sho_tools::nSHO(numax);
          auto const sigma = op.get_sigma(ia);
          if (echo > 6) printf("# %s atom #%i has numax=%d and %d coefficients, sigma= %g %s\n", __func__, ia, numax, ncoeff[ia], sigma*Ang,_Ang);
          assert(sigma > 0);
          scale_factors[ia] = sho_projection::get_sho_prefactors<real_t>(numax, sigma);
      } // ia
      data_list<complex_t> atom_coeff(ncoeff);

      view3D<complex_t const> const psi(eigenfunctions, nbands, g.all()); // wrap

      double const occupied_bands = control::get("devel.occupied.bands", 0.); // as long as the Fermi function is not in here

// #pragma omp parallel
      for(int ikpoint = 0; ikpoint < nkpoints; ++ikpoint) {
          double const kpoint_weight = 1; // depends on ikpoint
          auto const psi_k = psi[ikpoint];

          for(int iband = 0; iband < nbands; ++iband) {
              double const band_occupation = 2*std::min(std::max(0.0, occupied_bands - iband), 1.0); // depends on iband and ikpoint
              double const weight_nk = band_occupation * kpoint_weight;
              auto const psi_nk = psi_k[iband];

              if (weight_nk > 0) {
// #pragma omp for
                  for(size_t izyx = 0; izyx < g.all(); ++izyx) { // parallel
                      rho[izyx] += weight_nk * std::norm(psi_nk[izyx]);
                  } // izyx

                  // construct the atomic density matrix atom_rho
                  op.get_atom_coeffs(atom_coeff.data(), psi_nk, echo/2); // project again
                  for(int ia = 0; ia < na; ++ia) {
                      int const numax = op.get_numax(ia);
                      int const ncoeff = sho_tools::nSHO(numax);
                      auto const sf = scale_factors[ia].data();
                      // add to the atomic density matrix
                      for(int i = 0; i < ncoeff; ++i) {
                          auto const c_i = conjugate(atom_coeff[ia][i] * sf[i]);
#ifdef DEVEL
                          if (echo > 6) printf("# kpoint #%i band #%i atom #%i coeff[%i] = %.6e\t%g factor=%g\n",
                                              ikpoint, iband, ia, i, std::real(c_i), std::imag(c_i), sf[i]);
#endif // DEVEL
                          for(int j = 0; j < ncoeff; ++j) {
                              auto const c_j = atom_coeff[ia][j] * sf[j];
                              atom_rho[ia][i*ncoeff + j] += weight_nk * std::real(c_i * c_j);
                          } // j
                      } // i
                  } // ia
                  
              } // weight is positive
          } // iband
      } // ikpoint

      if (echo > 1) { printf("\n# Total valence density"); print_stats(rho, g.all(), g.dV()); }

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
