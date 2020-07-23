#include <cstdio> // printf
#include <cmath> // std::sqrt
#include <algorithm> // std::max
#include <complex> // std::complex<real_t>
#include <vector> // std::vector<T>
#include <cassert> // assert

#include "sho_hamiltonian.hxx"

#include "sho_potential.hxx" // ::load_local_potential
#include "geometry_analysis.hxx" // ::read_xyz_file
#include "control.hxx" // ::get
#include "display_units.h" // eV, _eV, Ang, _Ang
#include "real_space.hxx" // ::grid_t
#include "sho_tools.hxx" // ::nSHO, ::n1HO, ::order_*, ::SHO_index_t, ::construct_label_table
#include "boundary_condition.hxx" // Isolated_Boundary
#include "sho_overlap.hxx" // ::
#include "data_view.hxx" // view2D<T>, view3D<T>, view4D<T>
#include "linear_algebra.hxx" // ::eigenvalues, ::generalized_eigenvalues
#include "inline_tools.hxx" // align<nbits>

// #define FULL_DEBUG
#define DEBUG


#ifdef FULL_DEBUG
    #define full_debug(print) print 
#else
    #define full_debug(print)
#endif

#ifdef DEBUG
    #define debug(print) print 
#else
    #define debug(print)
#endif

namespace sho_hamiltonian {
  // computes Hamiltonian matrix elements between to SHO basis functions
  // including a PAW non-local contribution
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_Hamiltonian(int const echo=5) {
      status_t stat = 0;
      
      auto const vtotfile = control::get("sho_potential.test.vtot.filename", "vtot.dat"); // vtot.dat can be created by potential_generator.
      int dims[] = {0, 0, 0};
      std::vector<double> vtot; // total smooth potential
      stat += sho_potential::load_local_potential(vtot, dims, vtotfile, echo);

      auto const geo_file = control::get("geometry.file", "atoms.xyz");
      double *xyzZ = nullptr;
      int natoms{0};
      double cell[3] = {0, 0, 0}; 
      int bc[3] = {-7, -7, -7};
      { // scope: read atomic positions
          stat += geometry_analysis::read_xyz_file(&xyzZ, &natoms, geo_file, cell, bc, 0);
          if (echo > 2) printf("# found %d atoms in file \"%s\" with cell=[%.3f %.3f %.3f] %s and bc=[%d %d %d]\n",
                              natoms, geo_file, cell[0]*Ang, cell[1]*Ang, cell[2]*Ang, _Ang, bc[0], bc[1], bc[2]);
      } // scope

//    for(int d = 0; d < 3; ++d) assert(bc[d] == Isolated_Boundary); // ToDo: implement periodic images

      real_space::grid_t g(dims);
      g.set_grid_spacing(cell[0]/g[0], cell[1]/g[1], cell[2]/g[2]);
      if (echo > 1) printf("# use  %g %g %g %s grid spacing\n", g.h[0]*Ang, g.h[1]*Ang, g.h[2]*Ang, _Ang);
      if (echo > 1) printf("# cell is  %g %g %g %s\n", g.h[0]*g[0]*Ang, g.h[1]*g[1]*Ang, g.h[2]*g[2]*Ang, _Ang);
      double const origin[] = {.5*(g[0] - 1)*g.h[0],
                               .5*(g[1] - 1)*g.h[1], 
                               .5*(g[2] - 1)*g.h[2]};

      auto const center = new double[natoms][4]; // list of atomic centers
      for(int ia = 0; ia < natoms; ++ia) {
          for(int d = 0; d < 3; ++d) {
              center[ia][d] = xyzZ[ia*4 + d] + origin[d]; // w.r.t. to the center of grid point (0,0,0)
          }   center[ia][3] = 0; // 4th component is not used
      } // ia
      
      auto const usual_numax = int(control::get("sho_hamiltonian.test.numax", 1.));
      auto const usual_sigma =     control::get("sho_hamiltonian.test.sigma", 2.);
      std::vector<int>    numaxs(natoms, usual_numax); // define SHO basis sizes
      std::vector<double> sigmas(natoms, usual_sigma); // define SHO basis spreads
      double const sigma_asymmetry = control::get("sho_hamiltonian.test.sigma.asymmetry", 1.0);
      if (sigma_asymmetry != 1) { sigmas[0] *= sigma_asymmetry; sigmas[natoms - 1] /= sigma_asymmetry; } // manipulate the spreads
      ++numaxs[0]; // manipulate to test variable block sizes

      std::vector<int> offset(natoms + 1, 0); // the basis atoms localized on atom#i start at offset[i]
      int total_basis_size{0}, maximum_numax{-1};
      for(int ia = 0; ia < natoms; ++ia) {
          if (echo > 0) {
              printf("# atom#%i \tZ=%g \tposition %12.6f%12.6f%12.6f  numax= %d sigma= %.3f %s\n",
                  ia, xyzZ[ia*4 + 3], xyzZ[ia*4 + 0]*Ang, xyzZ[ia*4 + 1]*Ang, xyzZ[ia*4 + 2]*Ang, numaxs[ia], sigmas[ia]*Ang, _Ang);
          } // echo
          int const atom_basis_size = sho_tools::nSHO(numaxs[ia]);
          maximum_numax = std::max(maximum_numax, numaxs[ia]);
          offset[ia] = total_basis_size; // prefix sum
          total_basis_size += atom_basis_size;
      } // ia
//    int const numax_max = maximum_numax;
      offset[natoms] = total_basis_size;
      int const nB   = total_basis_size;
      int const nBa  = align<4>(nB); // memory aligned main matrix stride
      
      typedef double psi_t; // can also be std::complex<double>
      view3D<psi_t> SHW(3, nB, nBa, psi_t(0)); // get memory for Overlap, Hamiltonian and Work array

      std::vector<std::vector<view2D<psi_t>>> S_iaja(natoms); // construct sub-views for each atom pair
      std::vector<std::vector<view2D<psi_t>>> H_iaja(natoms); // construct sub-views for each atom pair
      for(int ia = 0; ia < natoms; ++ia) {
          S_iaja[ia].resize(natoms);
          H_iaja[ia].resize(natoms);
          int const nb_ia = sho_tools::nSHO(numaxs[ia]);
          int const n1i = sho_tools::n1HO(numaxs[ia]);
          for(int ja = 0; ja < natoms; ++ja) {
              S_iaja[ia][ja] = view2D<psi_t>(&(SHW(0,offset[ia],offset[ja])), nBa); // wrapper to sub-blocks of the overlap matrix
              H_iaja[ia][ja] = view2D<psi_t>(&(SHW(1,offset[ia],offset[ja])), nBa); // wrapper to sub-blocks of the overlap matrix
              int const nb_ja = sho_tools::nSHO(numaxs[ja]);
              int const n1j = sho_tools::n1HO(numaxs[ja]);
#if 0              
              for(int ib = 0; ib < nb_ia; ++ib) {
                  for(int jb = 0; jb < nb_ja; ++jb) {
                      S_iaja[ia][ja](ib,jb) = 10*ia + ja + ib*.1 + jb*.01; // dummy values to see if all matrix entries are set
                      H_iaja[ia][ja](ib,jb) = (ia + (ib + 1)*.1)*(ja + (jb + 1)*.1); // dummy values to see if all matrix entries are set
                  } // jb
              } // ib
#endif
              int const maxmoment = 0; // ToDo: adjust to the expansion of the local potential
              view4D<double> nabla2(3, 1, n1i + 1, n1j + 1, 0.0);      //  <\chi1D_i|d/dx  d/dx|\chi1D_j>
              view4D<double> ovl1Dm(3, 1 + maxmoment, n1i, n1j, 0.0);  //  <\chi1D_i| x^moment |\chi1D_j>
              for(int d = 0; d < 3; ++d) { // spatial directions x,y,z
                  double const distance = xyzZ[ia*4 + d] - xyzZ[ja*4 + d]; // does not account for periodic images, ToDo
                  stat += sho_overlap::nabla2_matrix(nabla2[d].data(), distance, n1j + 1, n1i + 1, sigmas[ja], sigmas[ia]);
                  stat += sho_overlap::moment_tensor(ovl1Dm[d].data(), distance, n1j, n1i, sigmas[ja], sigmas[ia], maxmoment);
              } // d
              
              // construct the overlap matrix of SHO basis functions 
              double const ones[1] = {1.0}; // expansion of the identity (constant==1) into x^{m_x} y^{m_y} z^{m_z}
              // Smat(i,j) := ovl_x(ix,jx) * ovl_y(iy,jy) * ovl_z(iz,jz)
              stat += sho_potential::potential_matrix(S_iaja[ia][ja], ovl1Dm, // input (dir,m,j,i) // ToDo: turn indices!!
                           ones, 0, numaxs[ia], numaxs[ja]);

              // construct the kinetic energy contribution
              double const kinetic[1] = {0.5}; // prefactor of kinetic energy in Hartree atomic units
              stat += sho_potential::potential_matrix(H_iaja[ia][ja], nabla2, // input (dir,m,j,i) // ToDo: turn indices!!
                           kinetic, 0, numaxs[ia], numaxs[ja]);
              
              // add the contribution of the local potential
              double Vcoeff[1] = {0.0}; // prelim, ToDo: fill with real values from expansion
              int const numax_V = 0; // expansion of the local potential into x^{m_x} y^{m_y} z^{m_z} around a given expansion center
              stat += sho_potential::potential_matrix(H_iaja[ia][ja], ovl1Dm, // input (dir,m,j,i) // ToDo: turn indices!!
                           Vcoeff, numax_V, numaxs[ia], numaxs[ja]);
              
              
          } // ja
      } // ia

      //
      // now construct the Hamiltonian:
      //   -- kinetic energy,             c.f. sho_overlap::test_simple_crystal
      //   -- local potential,            c.f. sho_potential::test_local_potential_matrix_elements (method=2)
      //   -- non-local PAW Hamiltonian contributions
      // and the overlap operator:
      //   -- SHO basis function overlap, c.f. sho_overlap::test_simple_crystal
      //   -- non-local PAW charge-deficit contributions
      //

      // some 1D preparations
      
      
//       int const ncut = sho_tools::n1HO(numax_max);
//       double const normalize = 0; // 0:do not normalize, we have to deal with an overlap matrix anyway -- ToDo: clarify
//       view3D<double>  Hp(natoms, ncut, ncut, 0.0); // 1D Hermite polynomials
//       view3D<double> dHp(natoms, ncut, ncut, 0.0); // first derivative of 1D Hermite polynomials
//       for(int ia = 0; ia < natoms; ++ia) {
//           double const sigma_inv = 1./sigmas[ia];
//           prepare_centered_Hermite_polynomials(Hp[ia].data(), ncut, sigma_inv, normalize);
// 
//           for(int n = 0; n < ncut; ++n) {
//               show the Hermite polynomial coefficients for H0
//               if (echo > 3) {
//                   printf("# H[%x]: ", n);
//                   for(int m = 0; m <= n; ++m) {
//                       printf("%8.4f", Hp(ia,n,m));
//                   }   printf("\n");
//               } // echo
// 
//               construct first derivatives
//               derive_Hermite_Gauss_polynomials(dHp(ia,n), Hp(ia,n), ncut, sigma_inv);
//           } // n
//       } // ia

      
      
      
      if (echo > 7) { // display S and H
          for(int s0h1 = 0; s0h1 < 2; ++s0h1) {
              printf("\n# %s matrix:\n", s0h1 ? "Hamiltonian" : "overlap");
              for(int iB = 0; iB < nB; ++iB) {
                  printf("# row%3i ", iB);
                  for(int jB = 0; jB < nB; ++jB) {
                      printf("%6.4f", SHW(s0h1,iB,jB));
                  } // jB
                  printf("\n");
              } // iB
          } // s0h1
          printf("\n");
      } // echo

      if (nullptr != center) delete[] center;
      if (nullptr != xyzZ) delete[] xyzZ;
      return stat;
  } // test_Hamiltonian

  status_t all_tests(int const echo) {
    status_t status(0);
    status += test_Hamiltonian(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace sho_hamiltonian
