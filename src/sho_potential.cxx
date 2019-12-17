#include <cstdio> // printf
#include <cstdlib> // abs
#include <cmath> // sqrt, exp
#include <fstream> // ifstream, ofstream
#include <algorithm> // max
#include <complex> // std::complex<real_t>
#include <complex>
#include <utility> // std::pair<T1,T2>, make_pair
#include <vector> // std::vector<T>
#include <array> // std::array<T,n>
#include <cassert> // assert

#include "sho_potential.hxx"

#include "geometry_analysis.hxx" // geometry_analysis::read_xyz_file
#include "constants.hxx" // pi, sqrtpi
#include "control.hxx" // control::get
#include "display_units.h" // Ang, _Ang
#include "real_space_grid.hxx" // real_space_grid::grid_t<N>
#include "sho_tools.hxx" // sho_tools::nSHO
#include "sho_projection.hxx" // ::sho_project, ::sho_add
#include "boundary_condition.hxx" // Isolated_Boundary
#include "overlap.hxx" // overlap::generate_product_tensor, ::generate_overlap_matrix
#include "data_view.hxx" // view2D<T>

// #include "quantum_numbers.h" // enn_QN_t, ell_QN_t, emm_QN_t
// #include "display_units.h" // eV, _eV, Ang, _Ang

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


namespace sho_potential {
  // computes potential matrix elements between to SHO basis functions
  
  template<typename real_t>
  status_t generate_potential_matrix(view2D<real_t> & Vmat, view3D<real_t> const & t, real_t const Vcoeff[],
                            int const numax_p, int const numax_i, int const numax_j) {
      // use the expansion of the product of two Hermite Gauss functions into another one, factorized in 3D
      int pxyz = 0;
      for    (int pz = 0; pz <= numax_p; ++pz) {
        for  (int py = 0; py <= numax_p - pz; ++py) {
          for(int px = 0; px <= numax_p - pz - py; ++px) {

            int ixyz = 0;
            for    (int iz = 0; iz <= numax_i; ++iz) {
              for  (int iy = 0; iy <= numax_i - iz; ++iy) {
                for(int ix = 0; ix <= numax_i - iz - iy; ++ix) {

                  int jxyz = 0;
                  for    (int jz = 0; jz <= numax_j; ++jz) {
                    for  (int jy = 0; jy <= numax_j - jz; ++jy) {
                      for(int jx = 0; jx <= numax_j - jz - jy; ++jx) {

                        Vmat[ixyz][jxyz] += t(px,ix,jx) * t(py,iy,jy) * t(pz,iz,jz) * Vcoeff[pxyz];

                        ++jxyz;
                      } // jx
                    } // jy
                  } // jz
                  assert(sho_tools::nSHO(numax_j) == jxyz);

                  ++ixyz;
                } // ix
              } // iy
            } // iz
            assert(sho_tools::nSHO(numax_i) == ixyz);

            ++pxyz;
          } // px
        } // py
      } // pz
      assert(sho_tools::nSHO(numax_p) == pxyz);
    
      return 0;
  } // generate_potential_matrix

  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_potential_elements(int const echo=5) {
      status_t stat = 0;
      int dims[] = {0, 0, 0};

      std::vector<double> vtot; // total smooth potential
      auto const filename = control::get("sho_potential.test.vtot.filename", "vtot.dat");
      { // scope: read in the potential from a file
          std::ifstream infile(filename);
          int npt = 0; 
          if (!infile.is_open()) {
              if (echo > 1) printf("# %s failed to open file %s\n",  __func__, filename);
              return 1; // failure
          }
          for(int d = 2; d >= 0; --d) {
              char sep;
              infile >> sep >> dims[d];
              if (echo > 3) printf("# found dim %c with %i grid points\n", 120+d, dims[d]);
          } // d
          size_t const all = dims[0]*dims[1]*dims[2];
          vtot.reserve(all);
          size_t idx;
          double val;
          while (infile >> idx >> val) {
              assert(idx < all);
              ++npt;
              vtot.push_back(val);
          } // while
          if (echo > 3) printf("# %s use %i values from file %s\n", __func__, npt, filename);
      } // scope
      real_space_grid::grid_t<1> g(dims);

      double *xyzZ = nullptr;
      int natoms = 0;
      double cell[3] = {0, 0, 0}; 
      int bc[3] = {-7, -7, -7};
      {
          auto const geofile = "atoms.xyz";
          stat += geometry_analysis::read_xyz_file(&xyzZ, &natoms, geofile, cell, bc, 0);
          if (echo > 2) printf("# found %d atoms in file \"%s\" with cell=[%.3f %.3f %.3f] %s and bc=[%d %d %d]\n",
                              natoms, geofile, cell[0]*Ang, cell[1]*Ang, cell[2]*Ang, _Ang, bc[0], bc[1], bc[2]);
      }
      
//    for(int d = 0; d < 3; ++d) assert(bc[d] == Isolated_Boundary); // ToDo: implement periodic images

      g.set_grid_spacing(cell[0]/dims[0], cell[1]/dims[1], cell[2]/dims[2]);
      if (echo > 1) printf("# use  %g %g %g  grid spacing in %s\n", g.h[0]*Ang,g.h[1]*Ang,g.h[2]*Ang,_Ang);
      if (echo > 1) printf("# cell is  %g %g %g  in %s\n", g.h[0]*g.dim(0)*Ang,g.h[1]*g.dim(1)*Ang,g.h[2]*g.dim(2)*Ang,_Ang);
      auto const center = new double[natoms][4]; // list of atomic centers
      for(int ia = 0; ia < natoms; ++ia) {
          for(int d = 0; d < 3; ++d) {
              center[ia][d] = xyzZ[ia*4 + d] + 0.5*(g.dim(d) - 1)*g.h[d]; // w.r.t. to the center of grid point (0,0,0)
          }   center[ia][3] = 0; // 4th component is not used
      } // ia
      
      int    const usual_numax = control::get("sho_potential.test.numax", 1.);
      double const usual_sigma = control::get("sho_potential.test.sigma", 2.);
      std::vector<int>    numaxs(natoms, usual_numax); // define SHO basis size
      std::vector<double> sigmas(natoms, usual_sigma); // define SHO basis spreads
      int numax_max = 0; for(int ia = 0; ia < natoms; ++ia) numax_max = std::max(numax_max, numaxs[ia]);

      int const method = control::get("sho_potential.test.method", 1.);
      if (1 & method) { // scope:
          if (echo > 2) printf("\n# %s Method=1\n", __func__);
          // Method 1) fully numerical -- cubically scaling, expensive
          //    for each pair of atoms and basis functions,
          //    add one basis function to an empty grid,
          //    multiply the potential, project with the other basis function
          std::vector<double> basis(g.all(), 0.0);
          int const mb = sho_tools::nSHO(numax_max);
          view3D<double> Vmat(natoms, mb, mb, 0.0);
          for(int i01 = 0; i01 <= 1; ++i01) { // 0:overlap, 1:potential
              if (echo > 1) printf("\n# %s\n", i01?"potential":"overlap");
              for(int ia = 0; ia < natoms; ++ia) {
                  int const nb = sho_tools::nSHO(numaxs[ia]);
                  std::vector<double> coeff(nb, 0.0);
                  for(int ib = 0; ib < nb; ++ib) {
                      set(basis.data(), g.all(), 0.0); coeff[ib] = 1; // delta
                      sho_projection::sho_add(basis.data(), g, coeff.data(), numaxs[ia], center[ia], sigmas[ia], 0);
                      // multiply Vtot to the basis function
                      if(i01) scale(basis.data(), basis.size(), vtot.data());
                      for(int ja = 0; ja < natoms; ++ja) {
                          int const mb = sho_tools::nSHO(numaxs[ja]);
                          std::vector<double> Vcoeff(mb, 0.0);
                          sho_projection::sho_project(Vcoeff.data(), numaxs[ja], center[ja], sigmas[ja], basis.data(), g, 0);
                          set(Vmat(ja,ib), mb, Vcoeff.data()); // copy data into result array
                      } // ja
                      coeff[ib] = 0;
                  } // ib

                  if (echo > 0) {
                      for(int ja = 0; ja < natoms; ++ja) {
                          printf("# ai#%i aj#%i\n", ia, ja);
                          for(int ib = 0; ib < mb; ++ib) {
                              printf("# ai#%i b#%i %c aj#%i ", ia, ib, i01?'V':'S', ja);
                              for(int jb = 0; jb < mb; ++jb) {
                                  printf("%8.3f", Vmat(ja,ib,jb)); // show potential matrix element
                              }   printf("\n");
                          } // ib
                      } // ja
                  } // echo
              } // ia
          } // i01
          if (echo > 2) printf("\n# %s ToDo: check if method=1 depends on absolute positions!\n", __func__);
      } // scope: Method 1

      if (2 & method) { // scope:
          if (echo > 2) printf("\n# %s Method=2\n", __func__);
          // Method 2) analytical -- quadratically scaling, linear scaling with truncation, cheap
          //    for each pair of atoms, find the center of weight,
          //    expand the potential in a SHO basis with sigma_V^2 = sigma_1^2 + sigma_2^2
          //    and numax_V = numax_1 + numax_2 and determine the
          //    potential matrix elements using the tensor
          double const central_pos[] = {.5*(g.dim(0) - 1)*g.h[0],
                                        .5*(g.dim(1) - 1)*g.h[1], 
                                        .5*(g.dim(2) - 1)*g.h[2]};
          for(int ia = 0; ia < natoms; ++ia) {
              auto const sigma2i = pow2(sigmas[ia]);
              for(int ja = 0; ja < natoms; ++ja) {
                  if (echo > 1) printf("# ai#%i aj#%i\n", ia, ja);
                
                  auto const sigma2j = pow2(sigmas[ja]);
                  auto const denom = 1./(sigma2i + sigma2j);
                  double const sigma_V = std::sqrt(sigma2i*sigma2j*denom);
                  double const wi = sigma2j*denom;
                  double const wj = sigma2i*denom;
                  double cnt[3];
                  for(int d = 0; d < 3; ++d) {
                      cnt[d] = central_pos[d] + wi*xyzZ[4*ia + d] + wj*xyzZ[4*ja + d];
                  } // d
                  auto const numax_V = numaxs[ia] + numaxs[ja];
                  int const nucut = 1 + std::max(numaxs[ia], numaxs[ja]);
                  view3D<double> t(2*nucut, nucut, nucut, 0.0);
                  overlap::generate_product_tensor(t.data(), nucut, sigma_V, sigmas[ia], sigmas[ja]);
                                              
                  int const nc = sho_tools::nSHO(numax_V);
                  std::vector<double> Vcoeff(nc, 0.0);
                  sho_projection::sho_project(Vcoeff.data(), numax_V, cnt, sigma_V, vtot.data(), g, 0);

                  int const nb = sho_tools::nSHO(numaxs[ia]);
                  int const mb = sho_tools::nSHO(numaxs[ja]);
                  view2D<double> Vmat(nb, mb, 0.0); // matrix
                  
                  // use the expansion of the product of two Hermite Gauss functions into another one
                  generate_potential_matrix(Vmat, t, Vcoeff.data(), numax_V, numaxs[ia], numaxs[ja]);

                  // display matrix
                  for(int ib = 0; ib < nb; ++ib) {
                      if (echo > 0) {
                          printf("# ai#%i b#%i V aj#%i  ", ia, ib, ja);
                          for(int jb = 0; jb < mb; ++jb) {
                              printf("%8.3f", Vmat[ib][jb]);
                          }   printf("\n");
                      } // echo
                  } // ib
                  
              } // ja
          } // ia
        
          if (echo > 2) printf("\n# %s method=2 seems symmetric!\n", __func__);
      } // scope: Method 2

      if (4 & method) { // scope:
          if (echo > 2) printf("\n# %s Method=3\n", __func__);
          // Method 3) approximated -- linear scaling, more expensive?
          //    for each atom expand the potential in a local SHO basis
          //    with spread sigma_V^2 = 2*sigma_1^2 at the atomic center,
          //    expand the other orbital in the local SHO basis (may need high lmax)
          //    also using the tensor
          //    the matrix elements will not converge with the same speed w.r.t. lmax
          //    so we will require symmetrization
          int const lmax = control::get("sho_potential.test.lmax", 2*usual_numax); // converge this!

          int const nucut = 1 + lmax;
          view3D<double> t(2*nucut, nucut, nucut, 0.0);
          overlap::generate_product_tensor(t.data(), nucut); // sigmap=2, sigma0=1, sigma1=1
          
          int const nc = sho_tools::nSHO(lmax);
          std::vector<double> Vcoeff(nc, 0.0);
          for(int ia = 0; ia < natoms; ++ia) {
              double const sigma = sigmas[ia];
              set(Vcoeff.data(), nc, 0.0);
              // project the potential onto an auxiliary SHO basis centered at the position of atom ia
              sho_projection::sho_project(Vcoeff.data(), lmax, center[ia], 2*sigma, vtot.data(), g, 0);
              
              int const nb = sho_tools::nSHO(numaxs[ia]);
              view2D<double> Vaux(nb, nc, 0.0);
              // now compute local matrix elements <local basis|V|large aux. basis>
              
              generate_potential_matrix(Vaux, t, Vcoeff.data(), lmax, numaxs[ia], lmax);

              
              for(int ja = 0; ja < natoms; ++ja) {
                  if (echo > 1) printf("# ai#%i aj#%i\n", ia, ja);
                  
                  int const mucut = 1 + numaxs[ja];
                  view3D<double> ovl(3, mucut, nucut);
                  for(int d = 0; d < 3; ++d) {
                      overlap::generate_overlap_matrix(ovl[d].data(), center[ja][d] - center[ia][d],
                                                  nucut, mucut, 
                                                  sigma, sigmas[ja]);
                  } // d

                  int const mb = sho_tools::nSHO(numaxs[ja]);
                  view2D<double> Vmat(nb, mb, 0.0);

                  { // scope: matrix multiply Vaux with overlap operator
                      for(int ixyz = 0; ixyz < nb; ++ixyz) {
                        
                          int const numax_j = numaxs[ja];
                          int jxyz = 0;
                          for    (int jz = 0; jz <= numax_j; ++jz) {
                            for  (int jy = 0; jy <= numax_j - jz; ++jy) {
                              for(int jx = 0; jx <= numax_j - jz - jy; ++jx) {

                                  double Vm = 0;

                                  int const numax_k = lmax;
                                  int kxyz = 0;
                                  for    (int kz = 0; kz <= numax_k; ++kz) {
                                    for  (int ky = 0; ky <= numax_k - kz; ++ky) {
                                      for(int kx = 0; kx <= numax_k - kz - ky; ++kx) {
                                  
                                          auto const Ovl_kj = ovl(0,jx,kx) * ovl(1,jy,ky) * ovl(2,jz,kz);
                                          Vm += Vaux[ixyz][kxyz] * Ovl_kj;

                                          ++kxyz;
                                      } // kx
                                    } // ky
                                  } // kz
                                  assert(nc == kxyz);

                                  Vmat[ixyz][jxyz] = Vm;
                                  
                                  ++jxyz;
                              } // jx
                            } // jy
                          } // jz
                          assert(mb == jxyz);

                      } // ixyz
                      
                  } // scope
                  
                  // display matrix
                  for(int ib = 0; ib < nb; ++ib) {
                      if (echo > 0) {
                          printf("# ai#%i b#%i V aj#%i  ", ia, ib, ja);
                          for(int jb = 0; jb < mb; ++jb) {
                              printf("%8.3f", Vmat[ib][jb]);
                          }   printf("\n");
                      } // echo
                  } // ib

              } // ja
          } // ia
        
          if (echo > 2) printf("\n# %s method=3 seems asymmetric!\n", __func__);
      } // scope: Method 3
     
      if (nullptr != xyzZ) delete[] xyzZ;
      return stat;
  } // test_potential_elements

  status_t all_tests() {
    auto status = 0;
    status += test_potential_elements(); // expensive
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace sho_potential
