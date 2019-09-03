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
  
  status_t init(int const echo=0) {
      status_t stat = 0;
//       int const na = 2; double const xyzZ[na][4] = {{-2,0,0, 13}, {2,0,0, 15}}; // Al-P
      int const na = 1; double const xyzZ[na][4] = {{0,0,0, 13}}; // Al only
//       int const dims[] = {160 + (na-1)*32, 160, 160}; double const grid_spacing = 0.125;
      int const dims[] = {80 + (na-1)*16, 80, 80}; double const grid_spacing = 0.25;
      real_space_grid::grid_t<double,1> g(dims);
      g.set_grid_spacing(grid_spacing);
      float Za[na]; // list of atomic numbers
      double center[na][4]; // list of atomic center
      if (echo > 1) printf("# %s List of Atoms: (coordinates in %s)\n", __func__,_Ang);
      for(int ia = 0; ia < na; ++ia) {
          int const iZ = (int)std::round(xyzZ[ia][3]);
          char const *El = &(element_symbols[2*iZ]); // warning, this is not a null-determined C-string
          if (echo > 4) printf("# %c%c  %16.9f%16.9f%16.9f\n", El[0],El[1], 
                         xyzZ[ia][0]*Ang, xyzZ[ia][1]*Ang, xyzZ[ia][2]*Ang);
          Za[ia] = xyzZ[ia][3]; // list of atomic numbers
          for(int d = 0; d < 3; ++d) center[ia][d] = 0.5*(g.dim(d) + 1)*g.h[d] - xyzZ[ia][d];
      } // ia

      double *rho_core[na]; // smooth core densities on r2-grids, nr2=2^12 points, ar2=16.f
      radial_grid_t *rg[na];
      double sigma_cmp[na]; //
      stat += single_atom::update(Za, na, rho_core, rg, sigma_cmp);
      
      double q00[na], v00[na];

      auto const rho = g.values;
      set(rho, g.all(), 0.0); // clear
      for(int ia = 0; ia < na; ++ia) {
          int const nr2 = 1 << 12; float const ar2 = 16.f;
          if (echo > 3) {
              printf("# Real-space smooth core density for atom #%d:\n", ia);
              for(int ir2 = 0; ir2 < nr2; ++ir2) {
                  printf("%g %g\n", std::sqrt(ir2/ar2), rho_core[ia][ir2]*Y00sq);
              }   printf("\n\n");
          } // echo
          double q_added;
          stat += real_space_grid::add_function(g, &q_added, rho_core[ia], nr2, ar2, center[ia], 9.f, Y00sq);
          if (echo > -1) {
              double s = 0; for(int i = 0; i < (int)g.all(); ++i) s += g.values[i]; s *= g.dV();
              printf("# integral over rho = %g after adding %g electrons smooth core density of atom #%d\n", s, q_added, ia);
          } // echo
          q00[ia] = -q_added; // spherical compensator
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
      set(Ves, g.all(), 0.0);
      if (1) { 
          for(int ia = 0; ia < na; ++ia) {
              // todo: add the compensators
              double const sigma_compensator = sigma_cmp[ia];
              if (1) {
                  double const sigma = sigma_compensator/std::sqrt(2.); // assume exp(-r^2/(2*sigma^2))
                  double const prefactor = 1./pow3(std::sqrt(2*constants::pi)*sigma);
                  // normalizing prefactor: 4 pi int dr r^2 exp(-r2/(2 sigma^2)) = sigma^3 \sqrt{8*pi^3}
                  double coeff[1];
                  set(coeff, 1, &q00[ia], prefactor);
                  stat += sho_projection::sho_add(rho, g, coeff, 0, center[ia], sigma, echo);
              } else if (1) {
                  float const rcut = 6.5*sigma_compensator;
                  double const prefactor = 1./(pow3(sigma_compensator*constants::sqrtpi));
                  float const ar2 = 64.f; 
                  int const nr2 = (int)std::ceil(ar2*pow2(rcut));
                  if (echo > -1) printf("# use r^2-grid with r^2 = %.1f*i with %d points for compensator of atom #%d\n", ar2, nr2, ia);
                  auto const rho_cmp = new double[nr2];
                  if (echo > -1) printf("# compensator charge density of atom #%d:\n", ia);
                  for(int ir2 = 0; ir2 < nr2; ++ir2) {
                      double const r2 = ir2/ar2;
                      rho_cmp[ir2] = prefactor*std::exp(-r2/pow2(sigma_compensator));
                      rho_cmp[ir2] *= q00[ia];
                      if (echo > 3) printf("%g %g\n", std::sqrt(r2), rho_cmp[ir2]);
                  }   if (echo > 3) printf("\n\n");
                  double q_added;
                  stat += real_space_grid::add_function(g, &q_added, rho_cmp, nr2, ar2, center[ia], rcut);
                  q00[ia] = q_added;
                  delete[] rho_cmp;
              }
              if (echo > -1) {
                  double s = 0; for(size_t i = 0; i < g.all(); ++i) s += g.values[i]; s *= g.dV();
                  printf("# integral over rho = %g after adding %g compensator electrons of atom #%d\n", s, q00[ia], ia);
              } // echo
          } // ia

          // report extremal values of the density on the grid
          double gmin = 9e9, gmax = -gmin, gsum = 0, gsum2 = 0;
          for(size_t i = 0; i < g.all(); ++i) {
              gmin = std::min(gmin, g.values[i]);
              gmax = std::max(gmax, g.values[i]);
              gsum += g.values[i];
              gsum2 += pow2(g.values[i]);
          } // i
          printf("\n# real-space grid stats min %g max %g avg %g for the density\n\n", gmin, gmax, gsum/g.all());
          
          int ng[3]; double reci[3][4]; 
          for(int d = 0; d < 3; ++d) { 
              ng[d] = g.dim(d);
              set(reci[d], 4, 0.0);
              reci[d][d] = 2*constants::pi/(ng[d]*g.h[d]);
          } // d
          stat += fourier_poisson::fourier_solve(Ves, rho, ng, reci);

          // test the potential in real space
          for(int ia = 0; ia < na; ++ia) {
              double const sigma = 2.0/std::sqrt(20.); // Bohr
              double const prefactor = 1./pow3(std::sqrt(2*constants::pi)*sigma);
              double coeff[1];
              set(coeff, 1, 0.0);
              stat += sho_projection::sho_project(coeff, 0, center[ia], sigma, Ves, g, echo);
              set(&v00[ia], 1, coeff, prefactor); // SHO-projectors are brought to the grid unnormalized, i.e. p_{00}(0) = 1.0
              printf("# potential projection for atom #%d v_00 = %g %s\n", ia, v00[ia]*eV,_eV);
          } // ia
          printf("# inner product between rho_aug and Ves = %g %s\n", dot_product(g.all(), rho, Ves)*g.dV()*eV,_eV);

          auto const fd = new finite_difference::finite_difference_t<double>(grid_spacing, 1, 12);
          g.values = Ves;
          stat += finite_difference::Laplacian(Laplace_Ves, g, *fd, -.25/constants::pi); // compute the Laplacian using high-order finite-differences
      }

      g.values = Ves; // analyze the electrostatic potential
//       g.values = Laplace_Ves; // analyze the augmented density computed as Laplacian*Ves
//       g.values = rho; // analyze the augmented density
//       g.values = Vxc; // analyze the xc potential
//       g.values = Vxc; add_product(g.values, g.all(), Ves, 1.0); // analyze the total potential: Vxc + Ves

      // report extremal values of what is stored on the grid
      double gmin = 9e9, gmax = -gmin, gsum = 0, gsum2 = 0;
      for(size_t i = 0; i < g.all(); ++i) {
          gmin = std::min(gmin, g.values[i]);
          gmax = std::max(gmax, g.values[i]);
          gsum += g.values[i];
          gsum2 += pow2(g.values[i]);
      } // i
      printf("\n# real-space grid stats min %g max %g avg %g\n\n", gmin, gmax, gsum/g.all());


      for(int ia = 0; ia < na; ++ia) {
//           int const nq = 200; float const dq = 1.f/16; // --> 199/16 = 12.4375 sqrt(Rydberg) =~= pi/(0.25 Bohr)
          float const dq = 1.f/16; int const nq = (int)(constants::pi/(grid_spacing*dq));
          auto const qc = new double[nq];
          
          printf("\n\n\n# start bessel_projection:\n\n"); // DEBUG
          stat += real_space_grid::bessel_projection(qc, nq, dq, g, center[ia]);
          printf("\n\n\n#   end bessel_projection.\n\n"); // DEBUG

          scale(qc, nq, pow2(solid_harmonics::Y00));
          
//           auto const qcq2 = new double[nq]; 
//           qcq2[0] = 0;
//           for(int iq = 1; iq < nq; ++iq) {
//               qcq2[iq] = qc[iq]/pow2(iq*dq); // cheap Poisson solver in Bessel transform
//           } // iq
// 
//           if (echo > 6) {
//               printf("# Bessel coeff for atom #%d:\n", ia);
//               for(int iq = 0; iq < nq; ++iq) {
//                   printf("%g %g %g\n", iq*dq, qc[iq], qcq2[iq]);
//               }   printf("\n\n");
//           } // echo

          if (echo > 3) {
              auto const rs = new double[rg[ia]->n];
              bessel_transform::transform_s_function(rs, qc, *rg[ia], nq, dq, true); // transform back to real-space again
              printf("# Real-space projection for atom #%d:\n", ia);
              for(int ir = 0; ir < rg[ia]->n; ++ir) {
                  printf("%g %g\n", rg[ia]->r[ir], rs[ir]); // seems like we are missing some factor
              }   printf("\n\n");
              
//               bessel_transform::transform_s_function(rs, qcq2, *rg[ia], nq, dq, true); // transform electrostatic solution to real-space
//               printf("# Hartree potential computed by Bessel transform for atom #%d:\n", ia);
//               for(int ir = 0; ir < rg[ia]->n; ++ir) { printf("%g %g\n", rg[ia]->r[ir], rs[ir]); } printf("\n\n");
              
              delete[] rs;
          } // echo
          delete[] qc;
      } // ia
      
      
      // compute the self-consistent solution of a single_atom, all states in the core
      // get the spherical core_density and bring it to the 3D grid
      // get the ell=0 compensator charge and add it to the 3D grid
      // envoke exchange_correlation and fourier_poisson
      // add XC and electrostatic potential and zero_potential contributions
      // project the total effective potential to each center using bessel_transforms
      // feed back spherical potential into single_atom
      
      delete[] Vxc;
      delete[] Ves;
      delete[] Laplace_Ves;
      for(int ia = 0; ia < na; ++ia) {
          delete[] rho_core[ia]; // has been allocated in single_atom::update()
      } // ia
          
      return stat;
  } // init

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo=5) {
      return init(echo);
  } // test_create_and_destroy

  status_t all_tests() {
    auto status = 0;
    status += test_create_and_destroy();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace spherical_atoms
