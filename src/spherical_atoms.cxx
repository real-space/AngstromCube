#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::floor

#include "spherical_atoms.hxx"

#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set, pow2
#include "constants.hxx" // Y00, sqrtpi
#include "solid_harmonics.hxx" // Y00
#include "real_space_grid.hxx" // grid_t, add_function
#include "radial_grid.h" // radial_grid_t
#include "chemical_symbol.h" // element_symbols
#include "single_atom.hxx" // update
#include "sho_projection.hxx" // sho_add
#include "exchange_correlation.hxx" // lda_PZ81_kernel
#include "fourier_poisson.hxx" // fourier_solve
#include "finite_difference.hxx" // Laplacian

// #define FULL_DEBUG
// #define DEBUG

namespace spherical_atoms {
  // this module makes a DFT calculation based on atoms
  // that live in a spherical potential which is found
  // by projecting the 3D potential.
  // their wave functions do not hybridize but they 
  // feel the effect of the density of neighboring atoms
  
  double constexpr Y00sq = pow2(solid_harmonics::Y00);
  
  void print_stats(double const values[], real_space_grid::grid_t<double,1> &g) {
      double gmin = 9e9, gmax = -gmin, gsum = 0, gsum2 = 0;
      for(size_t i = 0; i < g.all(); ++i) {
          gmin = std::min(gmin, values[i]);
          gmax = std::max(gmax, values[i]);
          gsum  += values[i];
          gsum2 += pow2(values[i]);
      } // i
      printf("# grid stats min %g max %g integral %g avg %g\n", gmin, gmax, gsum*g.dV(), gsum/g.all());
  } // print_stats

  status_t init(float const ion=0.f, int const echo=0) {
      status_t stat = 0;
      double constexpr Y00 = solid_harmonics::Y00;
      
      // compute the self-consistent solution of a single_atom, all states in the core
      // get the spherical core_density and bring it to the 3D grid
      // get the ell=0 compensator charge and add it to the 3D grid
      // envoke exchange_correlation and fourier_poisson
      // add XC and electrostatic potential and zero_potential contributions
      // project the total effective potential to each center using bessel_transforms
      // feed back spherical potential into single_atom
      
      int const na = 2; double const xyzZ[na][4] = {{-2,0,0, 13}, {2,0,0, 15}}; // Al-P
//       int const na = 1; double const xyzZ[na][4] = {{0,0,0, 13}}; // Al only
//       int const na = 2; double const xyzZ[na][4] = {{-2,0,0, 5}, {2,0,0, 7}}; // B-N
//       int const na = 1; double const xyzZ[na][4] = {{0,0,0, 3}}; // Li only
//       int const na = 2; double const xyzZ[na][4] = {{-2,0,0, 3}, {2,0,0, 9}}; // Li-F
      float ionization[na]; ionization[0] = ion*(na - 1); ionization[na - 1] = -ionization[0];
      
//       int const dims[] = {160 + (na-1)*32, 160, 160}; double const grid_spacing = 0.125; // very dense grid
//       int const dims[] = {80 + (na-1)*16, 80, 80}; double const grid_spacing = 0.25;
      int const dims[] = {160 + (na-1)*32, 160, 160}; double const grid_spacing = 0.25; // twice as large grid
      real_space_grid::grid_t<double,1> g(dims);
      g.set_grid_spacing(grid_spacing);
      
      double *qlm[na];
      double *vlm[na];
      float Za[na]; // list of atomic numbers
      double center[na][4]; // list of atomic centers
      if (echo > 1) printf("# %s List of Atoms: (coordinates in %s)\n", __func__,_Ang);
      for(int ia = 0; ia < na; ++ia) {
          int const iZ = (int)std::round(xyzZ[ia][3]);
          char const *El = &(element_symbols[2*iZ]); // warning, this is not a null-determined C-string
          if (echo > 4) printf("# %c%c  %16.9f%16.9f%16.9f\n", El[0],El[1], 
                         xyzZ[ia][0]*Ang, xyzZ[ia][1]*Ang, xyzZ[ia][2]*Ang);
          Za[ia] = xyzZ[ia][3]; // list of atomic numbers
          for(int d = 0; d < 3; ++d) {
              center[ia][d] = 0.5*(g.dim(d) + 1)*g.h[d] - xyzZ[ia][d];
          }   center[ia][3] = 0; // component is not used
          vlm[ia] = new double[1];
          qlm[ia] = new double[1];
      } // ia

      double *rho_core[na]; // smooth core densities on r2-grids, nr2=2^12 points, ar2=16.f
      radial_grid_t *rg[na];
      double sigma_cmp[na]; //
      stat += single_atom::update(na, Za, ionization, rho_core, rg, sigma_cmp, qlm);

      auto const rho = new double[g.all()];
      set(rho, g.all(), 0.0); // clear
      for(int ia = 0; ia < na; ++ia) {
          int const nr2 = 1 << 12; float const ar2 = 16.f;
          if (echo > 6) {
              printf("# Real-space smooth core density for atom #%d:\n", ia);
              for(int ir2 = 0; ir2 < nr2; ++ir2) {
                  double const r2 = ir2/ar2, r = std::sqrt(r2);
                  printf("%g %g\n", r, rho_core[ia][ir2]*Y00sq);
              }   printf("\n\n");
          } // echo
          double q_added;
          stat += real_space_grid::add_function(rho, g, &q_added, rho_core[ia], nr2, ar2, center[ia], 9.f, Y00sq);
          if (echo > -1) {
              printf("# after adding %g electrons smooth core density of atom #%d: ", q_added, ia);
              print_stats(rho, g);
          } // echo
//        qlm[ia][0] = -(q_added + ionization[ia])/Y00; // spherical compensator
          printf("# 00 compensator charge for atom #%d is %g\n", ia, qlm[ia][0]/Y00);
//           qlm[ia][0] = -(q_added)/Y00; // spherical compensator
//           printf("# 00 compensator charge for atom #%d is %g\n", ia, qlm[ia][0]*Y00);
      } // ia

      auto const Vxc = new double[g.all()];
      double Exc = 0, Edc = 0;
      for(size_t i = 0; i < g.all(); ++i) {
          Exc += rho[i]*exchange_correlation::lda_PZ81_kernel(rho[i], Vxc[i]);
          Edc += rho[i]*Vxc[i]; // double counting correction
      } // i
      Exc *= g.dV(); Edc *= g.dV(); // scale with volume element
      if (echo > -1) printf("# XC energy on grid %.12g %s, double counting %.12g %s\n", Exc*eV,_eV, Edc*eV,_eV);

      auto const Laplace_Ves = new double[g.all()];
      auto const         Ves = new double[g.all()];
      auto const         cmp = new double[g.all()];
      auto const        Vtot = new double[g.all()];
      set(Ves, g.all(), 0.0);
      set(cmp, g.all(), 0.0);
      if (1) { 
          for(int ia = 0; ia < na; ++ia) {
              // todo: add the compensators
              double const sigma_compensator = sigma_cmp[ia];
              double const prefactor = 1./(Y00*pow3(std::sqrt(2*constants::pi)*sigma_compensator));
              if (0) { // also works
                  float const rcut = 9.2*sigma_compensator;
                  float const ar2 = 64.f;
                  int const nr2 = (int)std::ceil(ar2*pow2(rcut));
                  if (echo > -1) printf("# use r^2-grid with r^2 = %.1f*i with %d points for compensator of atom #%d\n", ar2, nr2, ia);
                  auto const rho_cmp = new double[nr2];
                  if (echo > -1) printf("# compensator charge density of atom #%d:\n", ia);
                  double const sig2inv = -.5/pow2(sigma_compensator);
                  for(int ir2 = 0; ir2 < nr2; ++ir2) {
                      double const r2 = ir2/ar2;
                      rho_cmp[ir2] = prefactor*std::exp(sig2inv*r2);
                      rho_cmp[ir2] *= qlm[ia][0];
                      if (echo > 3) printf("%g %g\n", std::sqrt(r2), rho_cmp[ir2]);
                  }   if (echo > 3) printf("\n\n");
                  double q_added;
                  stat += real_space_grid::add_function(cmp, g, &q_added, rho_cmp, nr2, ar2, center[ia], rcut);
                  qlm[ia][0] = q_added*Y00;
                  delete[] rho_cmp;
              } else if (1) {
                  // normalizing prefactor: 4 pi int dr r^2 exp(-r2/(2 sigma^2)) = sigma^3 \sqrt{8*pi^3}, only for lm=00
                  double coeff[1];
                  set(coeff, 1, qlm[ia], prefactor);
                  stat += sho_projection::sho_add(cmp, g, coeff, 0, center[ia], sigma_compensator, echo);
              }
              if (echo > -1) {
                  // report extremal values of the density on the grid
                  printf("# after adding %g electrons compensator density for atom #%d:  ", qlm[ia][0]/Y00, ia);
                  print_stats(cmp, g);
              } // echo
          } // ia

          // add compensators cmp to rho
          add_product(rho, g.all(), cmp, 1.);
          printf("\n# augmented charge density grid stats: ");
          print_stats(rho, g);
          
          int ng[3]; double reci[3][4]; 
          for(int d = 0; d < 3; ++d) { 
              ng[d] = g.dim(d);
              set(reci[d], 4, 0.0);
              reci[d][d] = 2*constants::pi/(ng[d]*g.h[d]);
          } // d
          stat += fourier_poisson::fourier_solve(Ves, rho, ng, reci);

          // test the potential in real space, find ves_multipoles
          for(int ia = 0; ia < na; ++ia) {
              double const sigma_compensator = sigma_cmp[ia];
              double const prefactor = 1./(Y00*pow3(std::sqrt(2*constants::pi)*sigma_compensator));
              double coeff[1];
              set(coeff, 1, 0.0);
              stat += sho_projection::sho_project(coeff, 0, center[ia], sigma_compensator, Ves, g, echo);
              set(vlm[ia], 1, coeff, prefactor); // SHO-projectors are brought to the grid unnormalized, i.e. p_{00}(0) = 1.0
              printf("# potential projection for atom #%d v_00 = %g %s\n", ia, vlm[ia][0]*Y00*eV,_eV);
          } // ia
          printf("# inner product between cmp and Ves = %g %s\n", dot_product(g.all(), cmp, Ves)*g.dV()*eV,_eV);

          auto const fd = new finite_difference::finite_difference_t<double>(grid_spacing, 1, 12);
          stat += finite_difference::Laplacian(Laplace_Ves, Ves, g, *fd, -.25/constants::pi); // compute the Laplacian using high-order finite-differences
      } // if

      stat += single_atom::update(na, Za, ionization, nullptr, nullptr, nullptr, nullptr, vlm);

      set(Vtot, g.all(), Vxc); add_product(Vtot, g.all(), Ves, 1.);
      
      double* const value_pointers[] = {Ves, rho, Laplace_Ves, cmp, Vxc, Vtot};
//       values = Ves; // analyze the electrostatic potential
//       values = rho; // analyze the augmented density
//       values = Laplace_Ves; // analyze the augmented density computed as Laplacian*Ves
//       values = cmp; // analyze only the compensator density
//       values = Vxc; // analyze the xc potential
//       values = Vtot; // analyze the total potential: Vxc + Ves

      for(int iptr = 0; iptr < 1; ++iptr) { // only loop over the first 1 for electrostatics
          auto const values = value_pointers[iptr];
      
          // report extremal values of what is stored on the grid
          printf("\n# real-space grid stats: ");
          print_stats(values, g);

          for(int ia = 0; ia < na; ++ia) {
    //           int const nq = 200; float const dq = 1.f/16; // --> 199/16 = 12.4375 sqrt(Rydberg) =~= pi/(0.25 Bohr)
              float const dq = 1.f/16; int const nq = (int)(constants::pi/(grid_spacing*dq));
              auto const qc = new double[nq];
              
              printf("\n\n\n# start bessel_projection:\n\n"); // DEBUG
              stat += real_space_grid::bessel_projection(qc, nq, dq, values, g, center[ia]);
              printf("\n\n\n#   end bessel_projection.\n\n"); // DEBUG

              scale(qc, nq, pow2(solid_harmonics::Y00));
              
              auto const qcq2 = new double[nq]; 
              qcq2[0] = 0;
              for(int iq = 1; iq < nq; ++iq) {
                  qcq2[iq] = 4*constants::pi*qc[iq]/pow2(iq*dq); // cheap Poisson solver in Bessel transform
              } // iq
    
//               if (echo > 6) {
//                   printf("# Bessel coeff for atom #%d:\n", ia);
//                   for(int iq = 0; iq < nq; ++iq) {
//                       printf("%g %g %g\n", iq*dq, qc[iq], qcq2[iq]);
//                   }   printf("\n\n");
//               } // echo

              if (echo > 3) {
                  auto const rs = new double[rg[ia]->n];
                  bessel_transform::transform_s_function(rs, qc, *rg[ia], nq, dq, true); // transform back to real-space again
                  printf("# Real-space projection for atom #%d:\n", ia);
                  for(int ir = 0; ir < rg[ia]->n; ++ir) {
                      printf("%g %g\n", rg[ia]->r[ir], rs[ir]);
                  }   printf("\n\n");
                  
                  if ((values == rho) || (values == Laplace_Ves)) {
                      bessel_transform::transform_s_function(rs, qcq2, *rg[ia], nq, dq, true); // transform electrostatic solution to real-space
                      printf("# Hartree potential computed by Bessel transform for atom #%d:\n", ia);
                      for(int ir = 0; ir < rg[ia]->n; ++ir) { printf("%g %g\n", rg[ia]->r[ir], rs[ir]); } printf("\n\n");
                  } // density
                  
                  delete[] rs;
              } // echo
              delete[] qc;
          } // ia
      
      } // iptr loop for different quantities represented on the grid.
      
      
      delete[] Vxc;
      delete[] Ves;
      delete[] Vtot;
      delete[] Laplace_Ves;
      for(int ia = 0; ia < na; ++ia) {
          delete[] rho_core[ia]; // has been allocated in single_atom::update()
      } // ia

      stat += single_atom::update(-na, nullptr, nullptr); // cleanup
      return stat;
  } // init

 
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo=5) {
//       for(float ion = 1.0; ion >= -.01; ion -= 0.1) {
//           printf("\n\n\n#\n# Ionization = %g\n\n", ion);
//           init(ion, echo);
//       } // ion
//       return 0; // experiment, see ionization_result.* files
      return init(0.f, echo); // ionization of Al-P dimer by 3.0 electrons
  } // test_create_and_destroy

  status_t all_tests() {
    auto status = 0;
    status += test_create_and_destroy();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace spherical_atoms
