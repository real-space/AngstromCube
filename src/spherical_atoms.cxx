#include <cstdio> // printf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::floor
#include <vector> // std::vector

#include "spherical_atoms.hxx"

#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set, pow2
#include "inline_tools.hxx" // align<n>
#include "constants.hxx" // Y00, sqrtpi
#include "solid_harmonics.hxx" // Y00
#include "real_space_grid.hxx" // grid_t, add_function
#include "radial_grid.h" // radial_grid_t
#include "chemical_symbol.h" // element_symbols
#include "single_atom.hxx" // update
#include "sho_projection.hxx" // sho_add, sho_prefactors
#include "exchange_correlation.hxx" // lda_PZ81_kernel
#include "boundary_condition.hxx" // periodic_images
#include "fourier_poisson.hxx" // fourier_solve
#include "finite_difference.hxx" // Laplacian
#include "geometry_analysis.hxx" // read_xyz_file
#include "simple_timer.hxx" // SimpleTimer
#include "control.hxx" // control::get
#include "lossful_compression.hxx" // RDP_lossful_compression
#include "debug_output.hxx" // dump_to_file
#include "sho_unitary.hxx" // Unitary_SHO_Transform<real_t>

// #define FULL_DEBUG
// #define DEBUG

namespace spherical_atoms {
  // this module makes a DFT calculation based on atoms
  // that live in a spherical potential which is found
  // by projecting the 3D potential.
  // their wave functions do not hybridize but they 
  // feel the effect of the density of neighboring atoms
  
  double constexpr Y00sq = pow2(solid_harmonics::Y00);
  
  double print_stats(double const values[], size_t const all, double const dV=1, char const prefix=' ') {
      double gmin = 9e9, gmax = -gmin, gsum = 0, gsum2 = 0;
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
  
  status_t init(float const ion=0.f, int const echo=0) {
      SimpleTimer init_function_timer(__FILE__, __LINE__, __func__);
      status_t stat = 0;
      double constexpr Y00 = solid_harmonics::Y00;
      
      // compute the self-consistent solution of a single_atom, all states in the core
      // get the spherical core_density and bring it to the 3D grid
      // get the ell=0 compensator charge and add it to the 3D grid
      // envoke exchange_correlation and fourier_poisson
      // add XC and electrostatic potential and zero_potential contributions
      // project the total effective potential to each center using bessel_transforms
      // feed back spherical potential into single_atom
      
      double *xyzZ;
      int na = 0;
      double cell[3];
      int bc[3];
      stat += geometry_analysis::read_xyz_file(&xyzZ, &na, "atoms.xyz", cell, bc, echo);

      float ionization[na]; set(ionization, na, 0.f);
      if ((ion != 0.0) && (na > 1)) {
          if (echo > 2) printf("# %s distribute ionization of %g electrons between first and last atom\n", __func__, ion);
          ionization[0] = ion; ionization[na - 1] = -ionization[0];
      } // ionized

      // choose the box large enough not to require any periodic images
      double const h1 = 0.2378; // works for GeSbTe with alat=6.04
      int const dims[3] = {n_grid_points(cell[0]/h1), n_grid_points(cell[1]/h1), n_grid_points(cell[2]/h1)};
//       int const dims[] = {160 + (na-1)*32, 160, 160}; double const grid_spacing = 0.125; // very dense grid
//       int const dims[] = {80 + (na-1)*16, 80, 80}; double const grid_spacing = 0.25;
//       int const dims[] = {160 + (na-1)*32, 160, 160}; double const grid_spacing = 0.25; // twice as large grid
      if (echo > 1) printf("# use  %d x %d x %d  grid points\n", dims[0],dims[1],dims[2]);
      real_space_grid::grid_t<1> g(dims);
      g.set_grid_spacing(cell[0]/dims[0], cell[1]/dims[1], cell[2]/dims[2]);
      if (echo > 1) printf("# use  %g %g %g  %s grid spacing\n", g.h[0]*Ang,g.h[1]*Ang,g.h[2]*Ang,_Ang);
      if (echo > 1) printf("# cell is  %g %g %g  %s\n", g.h[0]*g.dim(0)*Ang,g.h[1]*g.dim(1)*Ang,g.h[2]*g.dim(2)*Ang,_Ang);
      for(int d = 0; d < 3; ++d) {
          if (std::abs(g.h[d]*g.dim(d) - cell[d]) >= 1e-6) {
              printf("# grid in %c-direction seems inconsistent, %d * %g differs from %g %s\n", 
                               'x'+d, g.dim(d), g.h[d]*Ang, cell[d]*Ang, _Ang);
          }
          assert(std::abs(g.h[d]*g.dim(d) - cell[d]) < 1e-6);
          assert(std::abs(g.h[d]*g.inv_h[d] - 1) < 4e-16);
      } // d
      double const min_grid_spacing = std::min(std::min(g.h[0], g.h[1]), g.h[2]);

      
      float Za[na]; // list of atomic numbers
      auto const center = new double[na][4]; // list of atomic centers
      
      if (echo > 1) printf("# %s List of Atoms: (coordinates in %s)\n", __func__,_Ang);
      for(int ia = 0; ia < na; ++ia) {
          int const iZ = (int)std::round(xyzZ[ia*4 + 3]);
          char const *El = &(element_symbols[2*iZ]); // warning, this is not a null-determined C-string
          if (echo > 4) printf("# %c%c  %16.9f%16.9f%16.9f\n", El[0],El[1],
                          xyzZ[ia*4 + 0]*Ang, xyzZ[ia*4 + 1]*Ang, xyzZ[ia*4 + 2]*Ang);
          Za[ia] = (float)xyzZ[ia*4 + 3]; // list of atomic numbers
          for(int d = 0; d < 3; ++d) {
              center[ia][d] = fold_back(xyzZ[ia*4 + d], cell[d]) + 0.5*(g.dim(d) - 1)*g.h[d]; // w.r.t. to the center of grid point (0,0,0)
          }   center[ia][3] = 0; // 4th component is not used
          if (echo > 1) printf("# relative%12.3f%16.3f%16.3f\n", center[ia][0]*g.inv_h[0],
                                       center[ia][1]*g.inv_h[1], center[ia][2]*g.inv_h[2]);
      } // ia
      
      float const rcut = 32; // radial grids usually and at 9.45 Bohr
      
      double *periodic_images = nullptr;
      int const n_periodic_images = boundary_condition::periodic_images(&periodic_images, cell, bc, rcut, echo);
      if (echo > 1) printf("# %s consider %d periodic images\n", __FILE__, n_periodic_images);
      
      double *rho_core[na]; // smooth core densities on r2-grids, nr2=2^12 points, ar2=16.f
      radial_grid_t *rg[na]; // smooth radial grid descriptors
      double sigma_cmp[na]; //
      double *qlm[na]; int lmax_qlm[na];
      double *vlm[na]; int lmax_vlm[na];
      stat += single_atom::update(na, Za, ionization, rg, sigma_cmp, nullptr, nullptr, nullptr, lmax_vlm, lmax_qlm);
      for(int ia = 0; ia < na; ++ia) {
          vlm[ia] = new double[pow2(1 + std::max(lmax_vlm[ia], 1))];
          qlm[ia] = new double[pow2(1 + lmax_qlm[ia])];
      } // ia

      sho_unitary::Unitary_SHO_Transform<double> const unitary(9);

      std::vector<double> Laplace_Ves(g.all());
      std::vector<double>         Ves(g.all());
      std::vector<double>         cmp(g.all());
      std::vector<double>         rho(g.all());
      std::vector<double>        Vtot(g.all());
      std::vector<double>         Vxc(g.all());

  int const max_scf_iterations = control::get("spherical_atoms.max.scf", 3.);
  for(int scf_iteration = 0; scf_iteration < max_scf_iterations; ++scf_iteration) {
      SimpleTimer scf_iteration_timer(__FILE__, __LINE__, "scf_iteration");
      printf("\n\n#\n# %s  SCF-Iteration #%d:\n#\n\n", __FILE__, scf_iteration);

      stat += single_atom::update(na, Za, ionization, nullptr, nullptr, rho_core, qlm);

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
              printf("\n## masked Real-space smooth core density for atom #%d:\n", ia);
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
          if (echo > -1) {
              printf("# after adding %g electrons smooth core density of atom #%d:", q_added, ia);
              print_stats(rho.data(), g.all(), g.dV());
          } // echo
//        qlm[ia][00] = -(q_added + ionization[ia])/Y00; // spherical compensator
          printf("# 00 compensator charge for atom #%d is %g\n", ia, qlm[ia][00]/Y00);
          printf("# added smooth core charge for atom #%d is %g\n", ia, q_added);
//        qlm[ia][00] = -(q_added)/Y00; // spherical compensator
      } // ia

      double Exc = 0, Edc = 0;
      for(size_t i = 0; i < g.all(); ++i) {
          Exc += rho[i]*exchange_correlation::lda_PZ81_kernel(rho[i], Vxc[i]);
          Edc += rho[i]*Vxc[i]; // double counting correction
      } // i
      Exc *= g.dV(); Edc *= g.dV(); // scale with volume element
      if (echo > -1) printf("# exchange-correlation energy on grid %.12g %s, double counting %.12g %s\n", Exc*eV,_eV, Edc*eV,_eV);

      set(Ves.data(), g.all(), 0.0);
      set(cmp.data(), g.all(), 0.0);
      { // scope
          for(int ia = 0; ia < na; ++ia) {
              // todo: add the compensators
              double const sigma_compensator = sigma_cmp[ia];
              if (1) {
                  int const nc = sho_tools::nSHO(lmax_qlm[ia]);
                  std::vector<double> coeff(nc, 0.0);
                  // normalizing prefactor: 4 pi int dr r^2 exp(-r2/(2 sigma^2)) = sigma^3 \sqrt{8*pi^3}, only for lm=00
                  double const prefactor = 1./(Y00*pow3(std::sqrt(2*constants::pi)*sigma_compensator)); // warning! only gets 00 correct
                  set(coeff.data(), 1, qlm[ia], prefactor); // copy only monopole moment, ToDo
                  for(int ii = 0; ii < n_periodic_images; ++ii) {
                      double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, &periodic_images[4*ii], 1.0);
                      stat += sho_projection::sho_add(cmp.data(), g, coeff.data(), lmax_qlm[ia], cnt, sigma_compensator, 0);
                  } // periodic images
              } // 1
              if (echo > -1) {
                  // report extremal values of the density on the grid
                  printf("# after adding %g electrons compensator density for atom #%d:", qlm[ia][00]/Y00, ia);
                  print_stats(cmp.data(), g.all(), g.dV());
              } // echo
          } // ia
          
          // add compensators cmp to rho
          add_product(rho.data(), g.all(), cmp.data(), 1.);
          printf("\n# augmented charge density grid stats:");
          print_stats(rho.data(), g.all(), g.dV());

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
              double const sigma_compensator = sigma_cmp[ia];
              int const nc = sho_tools::nSHO(lmax_vlm[ia]);
              std::vector<double> coeff(nc, 0.0);
              for(int ii = 0; ii < n_periodic_images; ++ii) {
                  std::vector<double> coeff_image(nc, 0.0);
                  double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, &periodic_images[4*ii], 1.0);
                  stat += sho_projection::sho_project(coeff_image.data(), lmax_vlm[ia], cnt, sigma_compensator, Ves.data(), g, 0);
                  add_product(coeff.data(), nc, coeff_image.data(), 1.0); // need phase factors?
              } // periodic images
              // SHO-projectors are brought to the grid unnormalized, i.e. p_{000}(0) = 1.0 and p_{200}(0) = -.5

              std::vector<double> vEzyx(nc, 0.0);
              sho_projection::normalize_and_reorder_coefficients(vEzyx.data(), coeff.data(), lmax_vlm[ia], sigma_compensator, 1./Y00);

              // convert SHO-coefficients from order_Ezyx to order_nlm
              std::vector<double> vnlm(nc, 0.0);
              unitary.transform_vector(vnlm.data(), sho_tools::order_nlm, 
                                      vEzyx.data(), sho_tools::order_Ezyx, lmax_vlm[ia], 9);

              set(vlm[ia], pow2(1 + lmax_vlm[ia]), vnlm.data()); // neglect all nrn > 0 contributions
              
              printf("# potential projection for atom #%d v_00 = %g %s\n", ia, vlm[ia][00]*Y00*eV,_eV);
              if (lmax_vlm[ia] > 0) {
                  printf("# potential projection for atom #%d v_1 = %g %g %g %s\n", ia, 
                  vlm[ia][01]*Y00*eV, vlm[ia][02]*Y00*eV, vlm[ia][03]*Y00*eV,_eV);
              } // more than monopole

              set(vlm[ia], pow2(1 + lmax_vlm[ia]), 0.0);
              double const prefactor = 1./(Y00*pow3(std::sqrt(2*constants::pi)*sigma_compensator));
              vlm[ia][00] = prefactor*coeff[00]; // only monopole
          } // ia
          printf("# inner product between cmp and Ves = %g %s\n", dot_product(g.all(), cmp.data(), Ves.data())*g.dV()*eV,_eV);

      } // scope

      stat += single_atom::update(na, Za, ionization, nullptr, nullptr, nullptr, nullptr, vlm);

      set(Vtot.data(), g.all(), Vxc.data()); add_product(Vtot.data(), g.all(), Ves.data(), 1.);

  } // scf_iteration
  return 1; // warning! no cleanup has been run
  
//   printf("\n\n# Early exit in %s line %d\n\n", __FILE__, __LINE__); exit(__LINE__);

      { // scope: compute the Laplacian using high-order finite-differences
          int const fd_nn[3] = {12, 12, 12}; // nearest neighbors in the finite-difference approximation
          finite_difference::finite_difference_t<double> fd(g.h, bc, fd_nn);
          {   SimpleTimer timer(__FILE__, __LINE__, "finite-difference");
              stat += finite_difference::Laplacian(Laplace_Ves.data(), Ves.data(), g, fd, -.25/constants::pi);
          } // timer
      } // scope

      double* const value_pointers[] = {Ves.data(), rho.data(), Laplace_Ves.data(), cmp.data(), Vxc.data(), Vtot.data()};
//       values = Ves; // analyze the electrostatic potential
//       values = rho; // analyze the augmented density
//       values = Laplace_Ves; // analyze the augmented density computed as Laplacian*Ves
//       values = cmp; // analyze only the compensator density
//       values = Vxc; // analyze the xc potential
//       values = Vtot; // analyze the total potential: Vxc + Ves

      for(int iptr = 0; iptr < 1; ++iptr) { // only loop over the first 1 for electrostatics
          SimpleTimer timer(__FILE__, __LINE__, "Bessel-projection-analysis");
          auto const values = value_pointers[iptr];

          // report extremal values of what is stored on the grid
          printf("\n# real-space grid stats:"); print_stats(values, g.all(), g.dV());

          for(int ia = 0; ia < na; ++ia) {
    //           int const nq = 200; float const dq = 1.f/16; // --> 199/16 = 12.4375 sqrt(Rydberg) =~= pi/(0.25 Bohr)
              float const dq = 1.f/16; int const nq = (int)(constants::pi/(min_grid_spacing*dq));
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
      delete[] center;

      // generate a file which contains the full potential Vtot
      {   char title[96];
          sprintf(title, "%i x %i x %i", g.dim(2), g.dim(1), g.dim(0));
          dump_to_file("vtot.dat", Vtot.size(), Vtot.data(), nullptr, 1, 1, title, 9);
      }

      stat += single_atom::update(-na, nullptr, nullptr); // cleanup
      return stat;
  } // init

 
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo=5) {
      float const ion = control::get("spherical_atoms.test.ion", 0.0);
      return init(ion, echo); // ionization of Al-P dimer by -ion electrons
  } // test_create_and_destroy

  status_t all_tests() {
    auto status = 0;
    status += test_create_and_destroy();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace spherical_atoms
