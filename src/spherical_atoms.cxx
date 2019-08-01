#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::floor

#include "spherical_atoms.hxx"

#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set, pow2
#include "constants.hxx" // Y00
#include "solid_harmonics.hxx" // Y00
#include "real_space_grid.hxx" // grid_t, add_function
#include "radial_grid.h" // radial_grid_t
#include "chemical_symbol.h" // element_symbols
#include "single_atom.hxx" // update
#include "exchange_correlation.hxx" // lda_PZ81_kernel
#include "fourier_poisson.hxx" // lda_PZ81_kernel

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
      int const dims[] = {96, 80, 80};
      real_space_grid::grid_t<double,1> g(dims);
      g.set_grid_spacing(0.25);
      int const na = 2;
      double const xyzZ[na][4] = {{-2,0,0, 13}, {2,0,0, 15}}; // Al-P
      float Za[na];
      if (echo > 1) printf("# %s List of Atoms: (coordinates in %s)\n", __func__,_Ang);
      for(int ia = 0; ia < na; ++ia) {
          int const iZ = (int)std::round(xyzZ[ia][3]);
          if (echo > 4) printf("%c%c  %16.9f%16.9f%16.9f\n", element_symbols[2*iZ], element_symbols[2*iZ + 1], 
                                  xyzZ[ia][0]*Ang, xyzZ[ia][1]*Ang, xyzZ[ia][2]*Ang);
          Za[ia] = xyzZ[ia][3];
      } // ia
      
      double *rho_core[na]; // smooth core densities on r2-grids, nr2=2^11 points, ar2=16.f
      radial_grid_t *rg[na];
      stat += single_atom::update(Za, 2, rho_core, rg);
      
      auto const rho = g.values;
      set(rho, g.all(), 0.0); // clear
      for(int ia = 0; ia < na; ++ia) {
          double center[3]; for(int d = 0; d < 3; ++d) center[d] = 0.5*(g.dim(d) + 1)*g.h[d] - xyzZ[ia][d];
          stat += real_space_grid::add_function(g, rho_core[ia], 1 << 11, 16.f, center, 9.f, Y00sq);
          if (echo > -1) {
              double s = 0; for(int i = 0; i < (int)g.all(); ++i) s += g.values[i]; s *= g.dV();
              printf("# integral over rho = %g after adding smooth core density of atom #%d\n", s, ia);
          } // echo
      } // ia
     
      auto const exc = new double[g.all()];
      auto const Vxc = new double[g.all()];
      double Exc = 0, Edc = 0;
      for(size_t i = 0; i < g.all(); ++i) {
          exc[i] = exchange_correlation::lda_PZ81_kernel(rho[i], Vxc[i]);
          Exc += rho[i]*exc[i];
          Edc += rho[i]*Vxc[i]; // double counting correction
      } // i
      Exc *= g.dV();
      Edc *= g.dV();
      if (echo > -1) printf("# XC energy on grid %.12g %s, double counting %.12g %s\n", Exc*eV,_eV, Edc*eV,_eV);

      delete[] exc;
      auto const Ves = Vxc; // re-use the memory
      { 
          int ng[3]; double reci[3][4]; 
          for(int d = 0; d < 3; ++d) { 
              ng[d] = g.dim(d);
              set(reci[d], 4, 0.0);
              reci[d][d] = 2*constants::pi/(ng[d]*g.h[d]);
          } // d
          stat += fourier_poisson::fourier_solve(Ves, rho, ng, reci);
      }
      
      g.values = Ves;
      for(int ia = 0; ia < na; ++ia) {
          double center[3]; for(int d = 0; d < 3; ++d) center[d] = 0.5*(g.dim(d) + 1)*g.h[d] - xyzZ[ia][d];
          int const nq = 80; float const dq = 1.f/32;
          auto const qc = new double[nq];
          stat += real_space_grid::bessel_projection(qc, nq, dq, g, center);

          if (echo > 7) {
              printf("# Bessel coeff for atom #%d:\n", ia);
              for(int iq = 0; iq < nq; ++iq) {
                  printf("%g %g\n", iq*dq, qc[iq]);
              }   printf("\n\n");
          } // echo
          
          if (echo > -1) {
              auto const rs = new double[rg[ia]->n];
              bessel_transform::transform_s_function(rs, qc, *rg[ia], nq, dq, true); // transform back to real-space again
              printf("# Real-space projection for atom #%d:\n", ia);
              for(int ir = 0; ir < rg[ia]->n; ++ir) {
                  printf("%g %g\n", rg[ia]->r[ir], rs[ir]);
              }   printf("\n\n");
              delete[] rs;
          } // echo

      } // ia
      
      // compute the self-consistent solution of a single_atom, all states in the core
      // get the spherical core_density and bring it to the 3D grid
      // get the ell=0 compensator charge and add it to the 3D grid
      // envoke exchange_correlation and fourier_poisson
      // add XC and electrostatic potential and zero_potential contributions
      // project the total effective potential to each center using bessel_transforms
      // feed back spherical potential into single_atom
      
      return stat;
  } // init

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo=0) {
      return init(echo);
  } // test_create_and_destroy

  status_t all_tests() {
    auto status = 0;
    status += test_create_and_destroy();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace spherical_atoms
