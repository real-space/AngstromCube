#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::floor
#include <vector> // std::vector

#include "potential_generator.hxx"

#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set, pow2
#include "inline_tools.hxx" // align<n>
#include "constants.hxx" // ::sqrtpi, ::pi
#include "solid_harmonics.hxx" // ::Y00
#include "real_space_grid.hxx" // ::grid_t, ::add_function
#include "radial_grid.h" // ::radial_grid_t
#include "chemical_symbol.h" // element_symbols
#include "single_atom.hxx" // ::update
#include "sho_projection.hxx" // ::sho_add, ::sho_project
#include "exchange_correlation.hxx" // lda_PZ81_kernel
#include "boundary_condition.hxx" // ::periodic_images
#include "data_view.hxx" // view2D<T>
#include "fourier_poisson.hxx" // ::fourier_solve
#include "finite_difference.hxx" // ::Laplacian
#include "geometry_analysis.hxx" // ::read_xyz_file
#include "simple_timer.hxx" // // SimpleTimer
#include "control.hxx" // control::get
#include "lossful_compression.hxx" // RDP_lossful_compression
#include "debug_output.hxx" // dump_to_file
#include "sho_unitary.hxx" // ::Unitary_SHO_Transform<real_t>

// #define FULL_DEBUG
// #define DEBUG

namespace potential_generator {
  // this module makes a DFT calculation based on atoms
  // that live in a spherical potential which is found
  // by projecting the 3D potential.
  // their wave functions do not hybridize but they 
  // feel the effect of the density of neighboring atoms
  
  double constexpr Y00sq = pow2(solid_harmonics::Y00);
  
  double print_stats(double const values[], size_t const all, double const dV=1, char const prefix=' ') {
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

//   inline int n_grid_points(double const suggest) { return (int)align<1>((int)std::ceil(suggest)); }
  inline int even(int const any) { return (((any - 1) >> 1) + 1) << 1;}
  inline int n_grid_points(double const suggest) { return (int)even((int)std::ceil(suggest)); }
  
  inline double fold_back(double x, double const cell_extend) { 
      while(x > 0.5*cell_extend) x -= cell_extend; 
      while(x < -.5*cell_extend) x += cell_extend;
      return x;
  } // fold_back
  
  status_t init_geometry_and_grid(real_space_grid::grid_t<1> & g, double **coordinates_and_Z, 
                                  int & natoms, int bc[3], int const echo=0) {
      // SimpleTimer init_function_timer(__FILE__, __LINE__, __func__, echo);
      status_t stat{0};
      
      // compute the self-consistent solution of a single_atom, all states in the core
      // get the spherical core_density and bring it to the 3D grid
      // get the ell=0 compensator charge and add it to the 3D grid
      // envoke exchange_correlation and fourier_poisson
      // add XC and electrostatic potential and zero_potential contributions
      // project the total effective potential to each center using bessel_transforms
      // feed back spherical potential into single_atom

      natoms = 0;
      double cell[3];
      auto const geo_file = control::get("geometry.file", "atoms.xyz");
      stat += geometry_analysis::read_xyz_file(coordinates_and_Z, &natoms, geo_file, cell, bc, echo);

      // choose the box large enough not to require any periodic images
      double const h1 = 0.2378; // works for GeSbTe with alat=6.04
      int const dims[3] = {n_grid_points(cell[0]/h1), n_grid_points(cell[1]/h1), n_grid_points(cell[2]/h1)};
      if (echo > 1) printf("# use  %d x %d x %d  grid points\n", dims[0],dims[1],dims[2]);
      g = real_space_grid::grid_t<1>(dims);
      g.set_grid_spacing(cell[0]/dims[0], cell[1]/dims[1], cell[2]/dims[2]);
      if (echo > 1) printf("# use  %g %g %g  %s grid spacing\n", g.h[0]*Ang,g.h[1]*Ang,g.h[2]*Ang,_Ang);
      if (echo > 1) printf("# cell is  %g %g %g  %s\n", g.h[0]*g.dim(0)*Ang,g.h[1]*g.dim(1)*Ang,g.h[2]*g.dim(2)*Ang,_Ang);
      for(int d = 0; d < 3; ++d) {
          if (std::abs(g.h[d]*g.dim(d) - cell[d]) >= 1e-6) {
              warn("# grid in %c-direction seems inconsistent, %d * %g differs from %g %s", 
                     'x'+d, g.dim(d), g.h[d]*Ang, cell[d]*Ang, _Ang);
          }
          assert(std::abs(g.h[d]*g.dim(d) - cell[d]) < 1e-6);
          assert(std::abs(g.h[d]*g.inv_h[d] - 1) < 4e-16);
      } // d
      
      return stat;
  } // init_geometry_and_grid
  
  status_t init(float const ion=0.f, int const echo=0) {
    
      double constexpr Y00 = solid_harmonics::Y00;

      status_t stat{0};
      
      real_space_grid::grid_t<1> g;
      double *coordinates_and_Z{nullptr};
      int na{0};
      int bc[3];
      stat += init_geometry_and_grid(g, &coordinates_and_Z, na, bc, echo);

      double cell[3];
      for(int d = 0; d < 3; ++d) {
          cell[d] = g.dim(d)*g.h[d];
      } // d
     
      float ionization[na]; set(ionization, na, 0.f);
      if ((ion != 0.0) && (na > 1)) {
          if (echo > 2) printf("# %s distribute ionization of %g electrons between first and last atom\n", __func__, ion);
          ionization[0] = ion; ionization[na - 1] = -ionization[0];
      } // ionized
      
      view2D<double const> const xyzZ(coordinates_and_Z, 4); // wrap
      float Za[na]; // list of atomic numbers
      
      view2D<double> center(na, 4); // get memory for a list of atomic centers
      if (echo > 1) printf("# %s List of Atoms: (coordinates in %s)\n", __func__,_Ang);
      for(int ia = 0; ia < na; ++ia) {
          int const iZ = (int)std::round(xyzZ[ia][3]);
          char const *El = &(element_symbols[2*iZ]); // warning, this is not a null-determined C-string
          if (echo > 4) printf("# %c%c   %15.9f %15.9f %15.9f\n", El[0],El[1],
                          xyzZ[ia][0]*Ang, xyzZ[ia][1]*Ang, xyzZ[ia][2]*Ang);
          Za[ia] = (float)xyzZ[ia][3]; // list of atomic numbers
          for(int d = 0; d < 3; ++d) {
              center[ia][d] = fold_back(xyzZ[ia][d], cell[d]) + 0.5*(g.dim(d) - 1)*g.h[d]; // w.r.t. to the center of grid point (0,0,0)
          }   center[ia][3] = 0; // 4th component is not used
          if (echo > 1) printf("# relative%12.3f%16.3f%16.3f\n", center[ia][0]*g.inv_h[0],
                                       center[ia][1]*g.inv_h[1], center[ia][2]*g.inv_h[2]);
      } // ia

      float const rcut = 32; // radial grids usually and at 9.45 Bohr
      
      double *periodic_images = nullptr;
      int const n_periodic_images = boundary_condition::periodic_images(&periodic_images, cell, bc, rcut, echo);
      if (echo > 1) printf("# %s consider %d periodic images\n", __FILE__, n_periodic_images);
      
      double *rho_core[na]; // smooth core densities on r2-grids, nr2=2^12 points, ar2=16.f
      double *zero_pot[na]; // smooth zero_potential on r2-grids, nr2=2^12 points, ar2=16.f
      radial_grid_t *rg[na]; // smooth radial grid descriptors
      double sigma_cmp[na]; //
      double *qlm[na]; int lmax_qlm[na];
      double *vlm[na]; int lmax_vlm[na];
      stat += single_atom::update(na, Za, ionization, rg, sigma_cmp, nullptr, nullptr, nullptr, lmax_vlm, lmax_qlm);
      for(int ia = 0; ia < na; ++ia) {
          vlm[ia] = new double[pow2(1 + lmax_vlm[ia])];
          qlm[ia] = new double[pow2(1 + lmax_qlm[ia])];
      } // ia

      sho_unitary::Unitary_SHO_Transform<double> const unitary(9);

      std::vector<double> Laplace_Ves(g.all());
      std::vector<double>         Ves(g.all());
      std::vector<double>         cmp(g.all());
      std::vector<double>         rho(g.all());
      std::vector<double>        Vtot(g.all());
      std::vector<double>         Vxc(g.all());

  int const max_scf_iterations = control::get("potential_generator.max.scf", 3.);
  for(int scf_iteration = 0; scf_iteration < max_scf_iterations; ++scf_iteration) {
      // SimpleTimer scf_iteration_timer(__FILE__, __LINE__, "scf_iteration", echo);
      if (echo > 1) printf("\n\n#\n# %s  SCF-Iteration #%d:\n#\n\n", __FILE__, scf_iteration);

      stat += single_atom::update(na, nullptr, nullptr, nullptr, nullptr, rho_core, qlm);

      set(rho.data(), g.all(), 0.0); // clear
      for(int ia = 0; ia < na; ++ia) {
          // ToDo: these parameters are silently assumed for the r2-grid of rho_core
          int const nr2 = 1 << 12; float const ar2 = 16.f; // rcut = 15.998 Bohr
          double const r2cut = pow2(rg[ia]->rmax), r2inv = 1./r2cut;
          if (echo > 6) {
              printf("\n## Real-space smooth core density for atom #%d:\n", ia);
              for(int ir2 = 0; ir2 < nr2; ++ir2) {
                  double const r2 = ir2/ar2, r = std::sqrt(r2);
                  printf("%g %g\n", r, rho_core[ia][ir2]*Y00sq);
              }   printf("\n\n");
          } // echo

          if (1) {
              for(int ir2 = 0; ir2 < nr2; ++ir2) {
                  double const r2 = ir2/ar2;
                  rho_core[ia][ir2] *= (r2 < r2cut) ? pow8(1. - pow8(r2*r2inv)) : 0;
              } // irs
          } // mask out the high frequency oscillations that appear from Bessel transfer to the r2-grid

          if (echo > 6) {
              printf("\n## masked real-space smooth core density for atom #%d:\n", ia);
              for(int ir2 = 0; ir2 < nr2; ++ir2) {
                  double const r2 = ir2/ar2, r = std::sqrt(r2);
                  printf("%g %g\n", r, rho_core[ia][ir2]*Y00sq);
              }   printf("\n\n");
          } // echo

          double q_added = 0;
          for(int ii = 0; ii < n_periodic_images; ++ii) {
              double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, &periodic_images[4*ii], 1.0);
              double q_added_image = 0;
              stat += real_space_grid::add_function(rho.data(), g, &q_added_image, rho_core[ia], nr2, ar2, cnt, Y00sq);
//               if (echo > 7) printf("# %g electrons smooth core density of atom #%d added for image #%i\n", q_added_image, ia, ii);
              q_added += q_added_image;
          } // periodic images
          if (echo > 1) {
              printf("# after adding %g electrons smooth core density of atom #%d:", q_added, ia);
              print_stats(rho.data(), g.all(), g.dV());
          } // echo
          if (echo > 3) printf("# added smooth core charge for atom #%d is  %g electrons\n", ia, q_added);
          if (echo > 3) printf("#    00 compensator charge for atom #%d is %g electrons\n",  ia, qlm[ia][00]/Y00);
      } // ia

      { // scope: eval the XC potential
          double Exc{0}, Edc{0};
          for(size_t i = 0; i < g.all(); ++i) {
              auto const exc_i = exchange_correlation::lda_PZ81_kernel(rho[i], Vxc[i]);
              Exc += rho[i]*exc_i;
              Edc += rho[i]*Vxc[i]; // double counting correction
          } // i
          Exc *= g.dV(); Edc *= g.dV(); // scale with volume element
          if (echo > 2) printf("# exchange-correlation energy on grid %.12g %s, double counting %.12g %s\n", Exc*eV,_eV, Edc*eV,_eV);
      } // scope

      set(Ves.data(), g.all(), 0.0);
      set(cmp.data(), g.all(), 0.0);
      { // scope
          for(int ia = 0; ia < na; ++ia) {
              if (echo > 6) printf("# 00 compensator charge for atom #%d is %g\n", ia, qlm[ia][00]);
              double const sigma = sigma_cmp[ia];
              int    const ellmax = lmax_qlm[ia];
              std::vector<double> coeff(sho_tools::nSHO(ellmax), 0.0);
              double const prefactor = 1./(Y00*pow3(std::sqrt(2*constants::pi)*sigma)); // warning! only gets 00 correct
              if (0) { // old method, works only for the monopole
                  // normalizing prefactor: 4 pi int dr r^2 exp(-r2/(2 sigma^2)) = sigma^3 \sqrt{8*pi^3}, only for lm=00
                  set(coeff.data(), 1, qlm[ia], prefactor); // copy only monopole moment, ToDo
              } else {
                  if (echo > 6) printf("# theoretical value for atom #%d coeff[000] = %g\n", ia, qlm[ia][00]*prefactor);
                  stat += sho_projection::denormalize_electrostatics(coeff.data(), qlm[ia], ellmax, sigma, unitary, echo);
              }
              if (echo > 7) printf("# before SHO-adding compensators for atom #%d coeff[000] = %g\n", ia, coeff[0]);
              for(int ii = 0; ii < n_periodic_images; ++ii) {
                  double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, &periodic_images[4*ii], 1.0);
                  stat += sho_projection::sho_add(cmp.data(), g, coeff.data(), ellmax, cnt, sigma, 0);
              } // periodic images

              if (echo > 1) {
                  // report extremal values of the density on the grid
                  printf("# after adding %g electrons compensator density for atom #%d:", qlm[ia][00]/Y00, ia);
                  print_stats(cmp.data(), g.all(), g.dV());
              } // echo
              
          } // ia
          
          
          // add compensators cmp to rho
          add_product(rho.data(), g.all(), cmp.data(), 1.);
          if (echo > 1) {
              printf("\n# augmented charge density grid stats:");
              print_stats(rho.data(), g.all(), g.dV());
          } // echo

          {
              int ng[3]; double reci[3][4]; 
              for(int d = 0; d < 3; ++d) { 
                  ng[d] = g.dim(d);
                  set(reci[d], 4, 0.0);
                  reci[d][d] = 2*constants::pi/(ng[d]*g.h[d]);
              } // d
          
              // solve the Poisson equation
              stat += fourier_poisson::fourier_solve(Ves.data(), rho.data(), ng, reci);
          }

          // test the potential in real space, find ves_multipoles
          for(int ia = 0; ia < na; ++ia) {
              double const sigma = sigma_cmp[ia];
              int    const ellmax = lmax_vlm[ia];
              int const nc = sho_tools::nSHO(ellmax);
              std::vector<double> coeff(nc, 0.0);
              for(int ii = 0; ii < n_periodic_images; ++ii) {
                  std::vector<double> coeff_image(nc, 0.0);
                  double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, &periodic_images[4*ii], 1.0);
                  stat += sho_projection::sho_project(coeff_image.data(), ellmax, cnt, sigma, Ves.data(), g, 0);
                  add_product(coeff.data(), nc, coeff_image.data(), 1.0); // need phase factors? no since we only test the electrostatic
              } // periodic images
              // SHO-projectors are brought to the grid unnormalized, i.e. p_{000}(0) = 1.0 and p_{200}(0) = -.5

              if (0) { // old method, works only for the monopole
                  set(vlm[ia], pow2(1 + ellmax), 0.0);
                  double const prefactor = 1./(Y00*pow3(std::sqrt(2*constants::pi)*sigma));
                  vlm[ia][00] = prefactor*coeff[00]; // only monopole
              } else {
                  stat += sho_projection::renormalize_electrostatics(vlm[ia], coeff.data(), ellmax, sigma, unitary, echo);
              }

              if (echo > 3) {
                  printf("# potential projection for atom #%d v_00 = %.9f %s\n", ia, vlm[ia][00]*Y00*eV,_eV);
                  int const ellmax_show = std::min(ellmax, 2);
                  for(int ell = 1; ell <= ellmax_show; ++ell) {
                      printf("# potential projection for atom #%d v_%im =", ia, ell);
                      double const unitfactor = Y00 * eV * std::pow(Ang, -ell);
                      for(int emm = -ell; emm <= ell; ++emm) {
                          int const lm = sho_tools::lm_index(ell, emm);
                          printf(" %.6f", vlm[ia][lm]*unitfactor);
                      } // emm
                      printf(" %s %s^%i\n", _eV, _Ang, -ell);
                  } // ell
              } // echo
              
          } // ia
          if (echo > 3) printf("# inner product between cmp and Ves = %g %s\n", dot_product(g.all(), cmp.data(), Ves.data())*g.dV()*eV,_eV);

      } // scope

      stat += single_atom::update(na, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, vlm, nullptr, nullptr, zero_pot);

      set(Vtot.data(), g.all(), Vxc.data()); add_product(Vtot.data(), g.all(), Ves.data(), 1.);
      
      if (echo > 1) {
          printf("\n# Total effective potential (before adding zero_potentials) grid stats:");
          print_stats(Vtot.data(), g.all(), g.dV());
      } // echo
      
      // now also add the zero_potential to Vtot
                  
      for(int ia = 0; ia < na; ++ia) {
          // ToDo: these parameters are silently assumed for the r2-grid of zero_pot
          int const nr2 = 1 << 12; float const ar2 = 16.f; // rcut = 15.998 Bohr
          double const r2cut = pow2(rg[ia]->rmax), r2inv = 1./r2cut;
          if (echo > 6) {
              printf("\n## zero-potential for atom #%d:\n", ia);
              for(int ir2 = 0; ir2 < nr2; ++ir2) {
                  double const r2 = ir2/ar2, r = std::sqrt(r2);
                  printf("%g %g\n", r, zero_pot[ia][ir2]*Y00);
              }   printf("\n\n");
          } // echo

          if (1) {
              for(int ir2 = 0; ir2 < nr2; ++ir2) {
                  double const r2 = ir2/ar2;
                  zero_pot[ia][ir2] *= (r2 < r2cut) ? pow8(1. - pow8(r2*r2inv)) : 0;
              } // irs
          } // mask out the high frequency oscillations that appear from Bessel transfer to the r2-grid

          if (echo > 6) {
              printf("\n## Masked zero-potential for atom #%d:\n", ia);
              for(int ir2 = 0; ir2 < nr2; ++ir2) {
                  double const r2 = ir2/ar2, r = std::sqrt(r2);
                  printf("%g %g\n", r, zero_pot[ia][ir2]*Y00);
              }   printf("\n\n");
          } // echo
          
          for(int ii = 0; ii < n_periodic_images; ++ii) {
              double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, &periodic_images[4*ii], 1.0);
              double dummy = 0;
              stat += real_space_grid::add_function(Vtot.data(), g, &dummy, zero_pot[ia], nr2, ar2, cnt, Y00);
          } // periodic images

      } // ia
      
      if (echo > 1) {
          printf("\n# Total effective potential  (after adding zero_potentials) grid stats:");
          print_stats(Vtot.data(), g.all(), g.dV());
      } // echo
      
  } // scf_iteration
  
//   return 1; // warning! no cleanup has been run
  
//   printf("\n\n# Early exit in %s line %d\n\n", __FILE__, __LINE__); exit(__LINE__);

      { // scope: compute the Laplacian using high-order finite-differences
          int const fd_nn[3] = {12, 12, 12}; // nearest neighbors in the finite-difference approximation
          finite_difference::finite_difference_t<double> const fd(g.h, bc, fd_nn);
          {   SimpleTimer timer(__FILE__, __LINE__, "finite-difference", echo);
              stat += finite_difference::Laplacian(Laplace_Ves.data(), Ves.data(), g, fd, -.25/constants::pi);
              { // Laplace_Ves should match rho
                  double res_a{0}, res_2{0};
                  for(size_t i = 0; i < g.all(); ++i) {
                      res_a += std::abs(Laplace_Ves[i] - rho[i]);
                      res_2 +=     pow2(Laplace_Ves[i] - rho[i]);
                  } // i
                  res_a *= g.dV(); res_2 = std::sqrt(res_2*g.dV());
                  if (echo > 1) printf("# Laplace*Ves - rho: residuals abs %.2e rms %.2e\n", res_a, res_2);
              }
          } // timer
          

          if (1) {
              // Jacobi solver for the Poisson equation: try if Ves is already the solution
              // Problem A*x == f    Jacobi update:  x_{k+1} = D^{-1}(f - T*x_{k}) where D+T==A and D is diagonal
            
              // SimpleTimer timer(__FILE__, __LINE__, "Jacobi solver", echo);
              finite_difference::finite_difference_t<double> fdj(g.h, bc, fd_nn);
              fdj.scale_coefficients(-.25/constants::pi);
              double const diagonal = fdj.clear_diagonal_elements();
              double const inverse_diagonal = 1/diagonal;
              auto & Ves_copy = Ves;
              std::vector<double> jac_work(g.all());
              for(int k = 0; k < 3; ++k) {
                  double const *x_k = Ves_copy.data(); // read from x_k
                  double      *Tx_k = jac_work.data();
                  double     *x_kp1 = Ves_copy.data();
                  stat += finite_difference::Laplacian(Tx_k, x_k, g, fdj);
                  double res_a{0}, res_2{0}; // internal change between x_{k} and x_{k+1}
                  for(size_t i = 0; i < g.all(); ++i) {
                      double const old_value = x_k[i];
                      double const new_value = inverse_diagonal*(rho[i] - Tx_k[i]);
                      res_a += std::abs(new_value - old_value);
                      res_2 +=     pow2(new_value - old_value);
                      x_kp1[i] = new_value;
                  } // i
                  res_a *= g.dV(); res_2 = std::sqrt(res_2*g.dV());
                  if (echo > 1) {
                      printf("# Jacobi iteration %i: residuals abs %.2e rms %.2e\n", k, res_a, res_2);
                      if (echo > 9) { printf("# it%i Ves grid stats:", k); print_stats(x_kp1, g.all(), g.dV()); }
                  } // echo
              } // k
          } // Jacobi solver
          
          {   SimpleTimer timer(__FILE__, __LINE__, "finite-difference onto Jacobi solution", echo);
              stat += finite_difference::Laplacian(Laplace_Ves.data(), Ves.data(), g, fd, -.25/constants::pi);
              { // Laplace_Ves should match rho
                  double res_a{0}, res_2{0};
                  for(size_t i = 0; i < g.all(); ++i) {
                      res_a += std::abs(Laplace_Ves[i] - rho[i]);
                      res_2 +=     pow2(Laplace_Ves[i] - rho[i]);
                  } // i
                  res_a *= g.dV(); res_2 = std::sqrt(res_2*g.dV());
                  if (echo > 1) printf("# Laplace*Ves - rho: residuals abs %.2e rms %.2e\n", res_a, res_2);
              }
          } // timer
          
          
      } // scope

      double* const value_pointers[] = {Ves.data(), rho.data(), Laplace_Ves.data(), cmp.data(), Vxc.data(), Vtot.data()};
//       values = Ves; // analyze the electrostatic potential
//       values = rho; // analyze the augmented density
//       values = Laplace_Ves; // analyze the augmented density computed as Laplacian*Ves
//       values = cmp; // analyze only the compensator density
//       values = Vxc; // analyze the xc potential
//       values = Vtot; // analyze the total potential: Vxc + Ves

      for(int iptr = 0; iptr < 1; iptr += 1) { // only loop over the first 1 for electrostatics
          // SimpleTimer timer(__FILE__, __LINE__, "Bessel-projection-analysis", echo);
          auto const values = value_pointers[iptr];

          // report extremal values of what is stored on the grid
          if (echo > 1) { printf("\n# real-space grid stats:"); print_stats(values, g.all(), g.dV()); }

          for(int ia = 0; ia < na; ++ia) {
    //           int const nq = 200; float const dq = 1.f/16; // --> 199/16 = 12.4375 sqrt(Rydberg) =~= pi/(0.25 Bohr)
              float const dq = 1.f/16; int const nq = (int)(constants::pi/(g.smallest_grid_spacing()*dq));
              std::vector<double> qc(nq, 0.0);

//            printf("\n\n# start bessel_projection:\n"); // DEBUG
              {
                  std::vector<double> qc_image(nq, 0.0);
                  for(int ii = 0; ii < n_periodic_images; ++ii) {
                      double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, &periodic_images[4*ii], 1.0);
                      stat += real_space_grid::bessel_projection(qc_image.data(), nq, dq, values, g, cnt);
                      add_product(qc.data(), nq, qc_image.data(), 1.0);
                  } // ii
              }
//            printf("\n# end bessel_projection.\n\n"); // DEBUG

              scale(qc.data(), nq, pow2(solid_harmonics::Y00));
              
              std::vector<double> qcq2(nq, 0.0);
              for(int iq = 1; iq < nq; ++iq) {
                  qcq2[iq] = 4*constants::pi*qc[iq]/pow2(iq*dq); // cheap Poisson solver in Bessel transform
              } // iq
    
//               if (echo > 6) {
//                   printf("\n## Bessel coeff for atom #%d:\n", ia);
//                   for(int iq = 0; iq < nq; ++iq) {
//                       printf("%g %g %g\n", iq*dq, qc[iq], qcq2[iq]);
//                   }   printf("\n\n");
//               } // echo

              if (echo > 3) {
                  std::vector<double> rs(rg[ia]->n);
                  bessel_transform::transform_s_function(rs.data(), qc.data(), *rg[ia], nq, dq, true); // transform back to real-space again
                  printf("\n## Real-space projection for atom #%d:\n", ia);
                  {   auto const mask = RDP_lossful_compression(rg[ia]->r, rs.data(), rg[ia]->n);
                      for(int ir = 0; ir < rg[ia]->n; ++ir) {
                          if (mask[ir]) printf("%g %g\n", rg[ia]->r[ir], rs[ir]);
                      }   printf("\n\n");
                  } // scope: RDP-mask

                  if ((values == rho.data()) || (values == Laplace_Ves.data())) {
                      bessel_transform::transform_s_function(rs.data(), qcq2.data(), *rg[ia], nq, dq, true); // transform electrostatic solution to real-space
                      printf("\n## Hartree potential computed by Bessel transform for atom #%d:\n", ia);
                      for(int ir = 0; ir < rg[ia]->n; ++ir) {
                          printf("%g %g\n", rg[ia]->r[ir], rs[ir]); 
                      }   printf("\n\n");
                  } // density
              } // echo
          } // ia
      
      } // iptr loop for different quantities represented on the grid.
      
      for(int ia = 0; ia < na; ++ia) {
          delete[] rho_core[ia]; // has been allocated in single_atom::update()
          delete[] vlm[ia];
          delete[] qlm[ia];
      } // ia

      // generate a file which contains the full potential Vtot
      if (echo > 0) { char title[96];
          sprintf(title, "%i x %i x %i", g.dim(2), g.dim(1), g.dim(0));
          dump_to_file("vtot.dat", Vtot.size(), Vtot.data(), nullptr, 1, 1, title, echo);
      } // unless all output is suppressed

      stat += single_atom::update(-na); // cleanup
      return stat;
  } // init

 
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_init(int const echo=3) {
      float const ion = control::get("potential_generator.test.ion", 0.0);
      return init(ion, echo); // ionization of Al-P dimer by -ion electrons
  } // test_init

  status_t all_tests(int const echo) {
    auto status = 0;
    status += test_init(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace potential_generator
