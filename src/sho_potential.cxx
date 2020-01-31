#include <cstdio> // printf
#include <cmath> // std::sqrt
#include <fstream> // std::ifstream
#include <algorithm> // std::max
#include <complex> // std::complex<real_t>
#include <utility> // std::pair<T1,T2>, std::make_pair
#include <vector> // std::vector<T>
#include <array> // std::array<T,n>
#include <cassert> // assert

#include "sho_potential.hxx"

#include "geometry_analysis.hxx" // ::read_xyz_file
#include "control.hxx" // ::get
#include "display_units.h" // eV, _eV, Ang, _Ang // ToDo
#include "real_space_grid.hxx" // ::grid_t<N>
#include "sho_tools.hxx" // ::nSHO
#include "sho_projection.hxx" // ::sho_project, ::sho_add, ::renormalize_coefficients
#include "boundary_condition.hxx" // Isolated_Boundary
#include "sho_overlap.hxx" // overlap::generate_product_tensor, 
                           // overlap::generate_overlap_matrix, 
                           // overlap::generate_potential_tensor
#include "data_view.hxx" // view2D<T>

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
  status_t generate_potential_matrix(view2D<real_t> & Vmat // result Vmat(i,j)
                            , view4D<real_t> const & t // input t(dir,p,j,i)
                            , real_t const Vcoeff[], int const numax_p // expansion of the potential in x^{nx}*y^{ny}*z^{nz}
                            , int const numax_i, int const numax_j
                            , int const dir01=1) { // 1:direction dependent input tensor, 0:isotropic
      // use the expansion of the product of two Hermite Gauss functions into another one, factorized in 3D
      // can use different tensors per direction
      int pzyx{0};
      for    (int pz = 0; pz <= numax_p;           ++pz) {
        for  (int py = 0; py <= numax_p - pz;      ++py) {
          for(int px = 0; px <= numax_p - pz - py; ++px) {

            int izyx{0};
            for    (int iz = 0; iz <= numax_i;           ++iz) {
              for  (int iy = 0; iy <= numax_i - iz;      ++iy) {
                for(int ix = 0; ix <= numax_i - iz - iy; ++ix) {

                  int jzyx{0};
                  for    (int jz = 0; jz <= numax_j;           ++jz) {  auto const tz   = t(dir01*2,pz,jz,iz);
                    for  (int jy = 0; jy <= numax_j - jz;      ++jy) {  auto const tyz  = t(dir01*1,py,jy,iy) * tz;
                      for(int jx = 0; jx <= numax_j - jz - jy; ++jx) {  auto const txyz = t(      0,px,jx,ix) * tyz;

                        Vmat[izyx][jzyx] += txyz * Vcoeff[pzyx];

                        ++jzyx;
                      } // jx
                    } // jy
                  } // jz
                  assert( sho_tools::nSHO(numax_j) == jzyx );

                  ++izyx;
                } // ix
              } // iy
            } // iz
            assert( sho_tools::nSHO(numax_i) == izyx );

            ++pzyx;
          } // px
        } // py
      } // pz
      assert( sho_tools::nSHO(numax_p) == pzyx );
    
      return 0;
  } // generate_potential_matrix

  
  template<typename real_t>
  status_t multiply_potential_matrix(view2D<real_t> & Vmat // result Vmat(i,j)
     , view3D<real_t> const & ovl  // input ovl(dir,j,i)
     , view2D<real_t> const & Vaux // input Vaux(i,k) --> ToDo: check if this is consistent!!
     , int const numax_i, int const numax_j, int const numax_k) {

      for(int izyx = 0; izyx < sho_tools::nSHO(numax_i); ++izyx) {
        
          int jzyx{0};
          for    (int jz = 0; jz <= numax_j; ++jz) {
            for  (int jy = 0; jy <= numax_j - jz; ++jy) {
              for(int jx = 0; jx <= numax_j - jz - jy; ++jx) {

                  double Vm{0};

                  int kzyx{0};
                  for    (int kz = 0; kz <= numax_k; ++kz) {              auto const tz   = ovl(2,jz,kz);
                    for  (int ky = 0; ky <= numax_k - kz; ++ky) {         auto const tyz  = ovl(1,jy,ky) * tz;
                      for(int kx = 0; kx <= numax_k - kz - ky; ++kx) {    auto const txyz = ovl(0,jx,kx) * tyz;
                  
                          Vm += Vaux[izyx][kzyx] * txyz;

                          ++kzyx;
                      } // kx
                    } // ky
                  } // kz
                  assert( sho_tools::nSHO(numax_k) == kzyx );

                  Vmat[izyx][jzyx] = Vm;
                  
                  ++jzyx;
              } // jx
            } // jy
          } // jz
          assert( sho_tools::nSHO(numax_j) == jzyx );

      } // izyx
      
      return 0;
  } // multiply_potential_matrix


  status_t normalize_coefficients(double coeff[], int const numax, double const sigma) {
      return sho_projection::renormalize_coefficients(coeff, coeff, numax, sigma);
  } // normalize_coefficients

  status_t normalize_potential_coefficients(double coeff[], int const numax, double const sigma, int const echo=0) {
      status_t stat = 0;
      int const nc = sho_tools::nSHO(numax);
      int const m = sho_tools::n1HO(numax);
      view2D<double> inv1D(m, m, 0.0); // get memory
      stat += sho_overlap::moment_normalization(inv1D.data(), m, sigma, echo);
      
      std::vector<double> c_new(nc, 0.0); // get memory

      view2D<char> labels(nc*(echo > 4), 8, '\0');
      if (echo > 4) sho_tools::construct_label_table(labels.data(), numax, sho_tools::order_zyx);
      if (echo > 4) printf("\n# %s numax=%i nc=%i sigma=%g %s\n", __func__, numax, nc, sigma*Ang, _Ang);
      int mzyx{0};
      for    (int mz = 0; mz <= numax; ++mz) {
        for  (int my = 0; my <= numax - mz; ++my) {
          for(int mx = 0; mx <= numax - mz - my; ++mx) {
            double cc{0};
            
            {
              // a specialty of the result matrix inv1D is that only matrix elements 
              //    connecting even-even or odd-odd indices are non-zero
              //    ToDo: this could be exploited in the following pattern
              int kzyx{0};
              for    (int kz = 0; kz <= numax; ++kz) {            auto const tz   = inv1D(mz,kz);
                for  (int ky = 0; ky <= numax - kz; ++ky) {       auto const tyz  = inv1D(my,ky) * tz;
                  for(int kx = 0; kx <= numax - kz - ky; ++kx) {  auto const txyz = inv1D(mx,kx) * tyz;

                    cc += txyz * coeff[kzyx];
                    
                    ++kzyx;
              }}} assert( nc == kzyx );
            }

            if (echo > 4) printf("# %s %s old%9.1e =%9.3f new%9.1e =%9.3f\n", __func__, labels[mzyx], coeff[mzyx], coeff[mzyx], cc, cc);
            c_new[mzyx] = cc; // write
            ++mzyx;
      }}} assert( nc == mzyx );
      if (echo > 4) printf("# %s\n\n", __func__);

      set(coeff, nc, c_new.data()); // copy the results back into the input array
      return stat;
  } // normalize_potential_coefficients
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_potential_elements(int const echo=5, char const *geofile="atoms.xyz") {
      status_t stat = 0;
      
//       int dims[] = {86, 86, 102}; std::vector<double> vtot(754392, 1.0); // total smooth potential, constant at 1 Hartree
      
      auto const filename = control::get("sho_potential.test.vtot.filename", "vtot.dat"); // vtot.dat was written by spherical_atoms.
      int dims[] = {0, 0, 0};
      std::vector<double> vtot; // total smooth potential
      if (1) { // scope: read in the potential from a file
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
      
      double *xyzZ = nullptr;
      int natoms = 0;
      double cell[3] = {0, 0, 0}; 
      int bc[3] = {-7, -7, -7};
      {
          stat += geometry_analysis::read_xyz_file(&xyzZ, &natoms, geofile, cell, bc, 0);
          if (echo > 2) printf("# found %d atoms in file \"%s\" with cell=[%.3f %.3f %.3f] %s and bc=[%d %d %d]\n",
                              natoms, geofile, cell[0]*Ang, cell[1]*Ang, cell[2]*Ang, _Ang, bc[0], bc[1], bc[2]);
      }
      
      
      
//    for(int d = 0; d < 3; ++d) assert(bc[d] == Isolated_Boundary); // ToDo: implement periodic images

      real_space_grid::grid_t<1> g(dims);
      g.set_grid_spacing(cell[0]/dims[0], cell[1]/dims[1], cell[2]/dims[2]);
      if (echo > 1) printf("# use  %g %g %g %s grid spacing\n", g.h[0]*Ang, g.h[1]*Ang, g.h[2]*Ang, _Ang);
      if (echo > 1) printf("# cell is  %g %g %g %s\n", g.h[0]*g.dim(0)*Ang, g.h[1]*g.dim(1)*Ang, g.h[2]*g.dim(2)*Ang, _Ang);
      double const central_pos[] = {.5*(g.dim(0) - 1)*g.h[0],
                                    .5*(g.dim(1) - 1)*g.h[1], 
                                    .5*(g.dim(2) - 1)*g.h[2]};

      auto const center = new double[natoms][4]; // list of atomic centers
      for(int ia = 0; ia < natoms; ++ia) {
          for(int d = 0; d < 3; ++d) {
              center[ia][d] = xyzZ[ia*4 + d] + central_pos[d]; // w.r.t. to the center of grid point (0,0,0)
          }   center[ia][3] = 0; // 4th component is not used
      } // ia

      if (0) { // scope: artificial linear potentials (other than constants)
          for        (int iz = 0; iz < dims[2]; ++iz) {   auto const z = iz*g.h[2] - central_pos[2];
              for    (int iy = 0; iy < dims[1]; ++iy) {   auto const y = iy*g.h[1] - central_pos[1];
                  for(int ix = 0; ix < dims[0]; ++ix) {   auto const x = ix*g.h[0] - central_pos[0];
                      int const izyx = (iz*dims[1] + iy)*dims[0] + ix;
//                    vtot[izyx] = 1.0 + x*.100 + y*.010 + z*.001;  // case 1xyz
//                    vtot[izyx] = 1.0 + z*.100 + y*.010 + x*.001;  // case 1zyx
//                    vtot[izyx] = 1.0 + y*.100 + z*.010 + x*.001;   // case 1yzx
                      //
                      // DevNotes: comparing Method1 and Method2, Al has sigma 1.65, P has sigma 1.35
                      // 
                      // found that in case 1xyz and case 1zyx on both, atom #0 and #1 (however, agreement on the atom-hopping terms)
                      // there are small deviations ~6--8% in small matrix elements of the order ~0.015 Hartree
                      // in particular (using zyx-labels) 
                      // 110 <--> 120 [y], 010 <--> 020 [y] and 011 <--> 021 [y]
                      // 
                      // found in case 1yzx similarly small deviations in
                      // 100 <--> 200 [z] and 101 <--> 201 [z] on atom #0 (Al) and
                      // 100 <--> 200 [z],    101 <--> 201 [z] and 110 <--> 210 [z] on atom#1 (P)
                      // 
//                    vtot[izyx] = 1.0 + y*.3;   // case 1y.3 // good
//                    vtot[izyx] = 1.0 + x*.3;   // case 1x.3 // good
//                    vtot[izyx] = 1.0 + z*.3;   // case 1z.3 // good
//                    vtot[izyx] = 1.0 + z*.2 + x*0.2;   // case 1zx2 // good
//                    vtot[izyx] = 1.0 + x*.2 + y*0.2;   // case 1xy2 // good
//                    vtot[izyx] = 1.0 + y*.2 + z*0.2;   // case 1yz2 // good
//                    vtot[izyx] = 1.0 + x*.2 + y*.2 + z*.2 ; // case 1xyz2 // good
//                    vtot[izyx] = 1.0 + x*.100 + y*.010 + z*.001;  // case 1xyz_new
                      // found that in case 1xyz_new
                      // there are small deviations ~.6--.8% in very small matrix elements of the order ~0.0015 Hartree
                      // 100 <--> 200 [z] on atom #1 (P) and
                      // 101 <--> 201 [z] and 110 <--> 210 on atom #0 (Al)
//                     vtot[izyx] = 1.0 + z*.100 + y*.010 + x*.001;  // case 1zyx_new 
                      // found in case 1zyx_new
                      // 001 <--> 002 [x], 011 <--> 012 [x] and 101 <--> 102 on atom #0 (Al) and the same on atom #1 (P)
                      vtot[izyx] = 1.0 + y*.100 + z*.010 + x*.001;   // case 1yzx_new
                      // found in case 1yzx_new
                      // 001 <--> 002 [x], 011 <--> 012 [x] and 101 <--> 102 on atom #0 (Al) and the same on atom #1 (P)
                      
                  } // ix
              } // iy
          } // iz
      } // scope: artificial potentials

      int    const usual_numax = control::get("sho_potential.test.numax", 1.);
      double const usual_sigma = control::get("sho_potential.test.sigma", 2.);
      std::vector<int>    numaxs(natoms, usual_numax); // define SHO basis size
      std::vector<double> sigmas(natoms, usual_sigma); // define SHO basis spreads
      sigmas[0] *= 1.1; sigmas[natoms - 1] *= 0.9; // manipulate the spreads
      int numax_max = 0; for(int ia = 0; ia < natoms; ++ia) numax_max = std::max(numax_max, numaxs[ia]);
      
      
      view3D<char> labels(1+numax_max, sho_tools::nSHO(numax_max), 8, '\0');
      for(int nu = 0; nu <= numax_max; ++nu) {
          sho_tools::construct_label_table(labels[nu].data(), nu, sho_tools::order_zyx);
      } // nu

      int const method = control::get("sho_potential.test.method", -1.); // bit-string, use method=7 to activate all
      if (1 & method) { // scope:
          if (echo > 2) printf("\n# %s Method=1\n", __func__);
          // Method 1) fully numerical -- cubically scaling, expensive
          //    for each pair of atoms and basis functions,
          //    add one basis function to an empty grid,
          //    multiply the potential, project with the other basis function
          std::vector<double> basis(g.all(), 0.0);
          std::vector<double> Vbasis(g.all(), 0.0);
          int const mxb = sho_tools::nSHO(numax_max);
          view4D<double> Vmat(natoms, natoms, mxb, mxb, 0.0);
          view4D<double> Smat(natoms, natoms, mxb, mxb, 0.0);
          view2D<double> Sdiag(natoms, mxb, 0.0);
          for(int ja = 0; ja < natoms; ++ja) {
              int const mb = sho_tools::nSHO(numaxs[ja]);
              std::vector<double> coeff(mb, 0.0);
              for(int jb = 0; jb < mb; ++jb) {
                  set(basis.data(), g.all(), 0.0); // clear
                  set(coeff.data(), mb, 0.0); coeff[jb] = 1; // delta
                  sho_projection::sho_add(basis.data(), g, coeff.data(), numaxs[ja], center[ja], sigmas[ja], 0);

                  // multiply Vtot to the basis function
                  product(Vbasis.data(), basis.size(), basis.data(), vtot.data());
                  
                  for(int ia = 0; ia < natoms; ++ia) {
                      int const nb = sho_tools::nSHO(numaxs[ia]);
                      std::vector<double> Scoeff(nb, 0.0);
                      std::vector<double> Vcoeff(nb, 0.0);
                      sho_projection::sho_project(Scoeff.data(), numaxs[ia], center[ia], sigmas[ia],  basis.data(), g, 0);
                      sho_projection::sho_project(Vcoeff.data(), numaxs[ia], center[ia], sigmas[ia], Vbasis.data(), g, 0);
                      for(int ib = 0; ib < nb; ++ib) {
                          Smat(ia,ja,ib,jb) = Scoeff[ib]; // copy result into large array
                          Vmat(ia,ja,ib,jb) = Vcoeff[ib]; // copy result into large array
                      } // ib
                  } // ia
                  Sdiag(ja,jb) = Smat(ja,ja,jb,jb); // store diagonal elements of the Overlap operator

              } // jb
          } // ja
          
          if (echo > 4) { // show the normalized operators
              for(int i01 = 0; i01 <= 1; ++i01) {
                  if (echo > 7) printf("\n# %s\n", i01?"potential (V)":"overlap (S)");
                  for(int ia = 0; ia < natoms; ++ia) {        int const nb = sho_tools::nSHO(numaxs[ia]);
                      for(int ja = 0; ja < natoms; ++ja) {    int const mb = sho_tools::nSHO(numaxs[ja]);
                          printf("# ai#%i aj#%i\n", ia, ja);
                          for(int ib = 0; ib < nb; ++ib) {
                              printf("# %c ai#%i aj#%i %s ", i01?'V':'S', ia, ja, labels(numaxs[ia],ib));
                              for(int jb = 0; jb < mb; ++jb) {
                                  double const elem = (i01 ? Vmat(ia,ja,ib,jb) : Smat(ia,ja,ib,jb)) 
                                                    / std::sqrt(Sdiag(ia,ib)*Sdiag(ja,jb));
                                  printf("%8.4f", elem); // show normalized overlap matrix element
                              } // jb
                              printf("\n");
                          } // ib
                      } // ja
                  } // ia
              } // i01
          } // echo
          
          if (echo > 2) printf("\n# %s ToDo: check if method=1 depends on absolute positions!\n", __func__);
          
      } // scope: Method 1

      if (2 & method) { // scope:
          if (echo > 2) printf("\n# %s Method=2\n", __func__);
          // Method 2) analytical -- quadratically scaling, linear scaling with truncation, cheap
          //    for each pair of atoms, find the center of weight,
          //    expand the potential in a SHO basis with sigma_V^-2 = sigma_1^-2 + sigma_2^-2
          //    and numax_V = numax_1 + numax_2 and determine the
          //    potential matrix elements using the tensor
          for(int ia = 0; ia < natoms; ++ia) {
              double const sigma_i = sigmas[ia];
              for(int ja = 0; ja < natoms; ++ja) {
                  double const sigma_j = sigmas[ja];
                  double const alpha_i = 1/pow2(sigma_i);
                  double const alpha_j = 1/pow2(sigma_j);
                  double const sigma_V = 1/std::sqrt(alpha_i + alpha_j);
                  double const wi = alpha_i*pow2(sigma_V);
                  double const wj = alpha_j*pow2(sigma_V);
                  assert( std::abs( wi + wj - 1.0 ) < 1e-12 );
                  double cnt[3]; // center of weight
                  for(int d = 0; d < 3; ++d) {
                      cnt[d] = wi*xyzZ[ia*4 + d] + wj*xyzZ[ja*4 + d];
                  } // d
                  int const numax_V = numaxs[ia] + numaxs[ja] + 1 + 3;
                  if (echo > 1) printf("# ai#%i aj#%i  center of weight: %g %g %g sigma_V: %g %s numax_V=%i\n", ia, ja, 
                                            cnt[0]*Ang, cnt[1]*Ang, cnt[2]*Ang, sigma_V*Ang, _Ang, numax_V);
                  for(int d = 0; d < 3; ++d) {
                      cnt[d] += central_pos[d];
                  } // d
                  int const nucut = sho_tools::n1HO(std::max(numaxs[ia], numaxs[ja]));
                  
                  view4D<double> t(3, 1+numax_V, nucut, nucut, 0.0);
                  for(int d = 0; d < 3; ++d) {
                      auto const distance = center[ja][d] - center[ia][d];
                      sho_overlap::moment_tensor(t[d].data(), distance, nucut, nucut, sigmas[ia], sigmas[ja], numax_V);
                  } // d
                  
                  int const nc = sho_tools::nSHO(numax_V);
                  std::vector<double> Vcoeff(nc, 0.0);
                  sho_projection::sho_project(Vcoeff.data(), numax_V, cnt, sigma_V, vtot.data(), g, 0);
                  stat += normalize_potential_coefficients(Vcoeff.data(), numax_V, sigma_V); // mute
                  // now Vcoeff is in a representation w.r.t. powers of the Cartesian coords x^{nx}*y^{ny}*z^{nz}

                  int const nb = sho_tools::nSHO(numaxs[ia]);
                  int const mb = sho_tools::nSHO(numaxs[ja]);
                  view2D<double> Vmat(nb, mb, 0.0); // matrix

                  // use the expansion of the product of two Hermite Gauss functions into another one
                  generate_potential_matrix(Vmat, t, Vcoeff.data(), numax_V, numaxs[ia], numaxs[ja]);
                  
                  // display matrix
                  for(int ib = 0; ib < nb; ++ib) {
                      if (echo > 0) {
                          printf("# V ai#%i aj#%i %s ", ia, ja, labels(numaxs[ia],ib));
                          for(int jb = 0; jb < mb; ++jb) {
                              printf("%8.4f", Vmat[ib][jb]);
                          }   printf("\n");
                      } // echo
                  } // ib
                  
              } // ja
          } // ia
        
          if (echo > 2) printf("\n# %s method=2 seems symmetric!\n", __func__);
          
      } // scope: Method 2

      if (4 & method) { // scope:
          int const lmax = control::get("sho_potential.test.lmax", 3*usual_numax); // converge this?
          if (echo > 2) printf("\n# %s Method=4 lmax=%i\n", __func__, lmax);
          // Method 4) approximated -- linear scaling, more expensive?
          //    for each atom expand the potential in a local SHO basis
          //    with spread sigma_V^2 = 2*sigma_1^2 at the atomic center,
          //    expand the other orbital in the local SHO basis (may need high lmax)
          //    also using the tensor.
          //    The matrix elements will not converge with the same speed w.r.t. lmax
          //    so we will require symmetrization

          int const nucut = sho_tools::n1HO(lmax);
          view4D<double> t(1, 2*nucut, nucut, nucut, 0.0);
          sho_overlap::generate_product_tensor(t.data(), nucut); // sigmap=2, sigma0=1, sigma1=1

          int const mxb = sho_tools::nSHO(numax_max);
          view4D<double> Smat(natoms, natoms, mxb, mxb, 0.0);
          
          int const nc = sho_tools::nSHO(lmax);
          std::vector<double> Vcoeff(nc, 0.0);
          for(int ia = 0; ia < natoms; ++ia) {
              double const sigma = sigmas[ia];
              set(Vcoeff.data(), nc, 0.0);
              double const sigma_V = std::sqrt(.5)*sigma;
              int const numax_V = lmax;
              // project the potential onto an auxiliary SHO basis centered at the position of atom ia
              sho_projection::sho_project(Vcoeff.data(), numax_V, center[ia], sigma_V, vtot.data(), g, 0);
              stat += normalize_potential_coefficients(Vcoeff.data(), numax_V, sigma_V, 9);
              // now Vcoeff is in a representation w.r.t. powers of the Cartesian coords x^{nx}*y^{ny}*z^{nz}
              
              int const nb = sho_tools::nSHO(numaxs[ia]);
              view2D<double> Vaux(nb, nc, 0.0); // Vaux(i,j)
              // now compute local matrix elements <local basis|V|large aux. basis>
              
              generate_potential_matrix(Vaux, t, Vcoeff.data(), lmax, numaxs[ia], lmax, 0);
              
              for(int ja = 0; ja < natoms; ++ja) {
                  if (echo > 1) printf("# ai#%i aj#%i\n", ia, ja);
                  
                  int const mucut = sho_tools::n1HO(numaxs[ja]);
                  view3D<double> ovl(3, mucut, nucut); // index order (dir,j,i)
                  for(int d = 0; d < 3; ++d) {
                      auto const distance = center[ja][d] - center[ia][d];
                      sho_overlap::generate_overlap_matrix(ovl[d].data(), distance, nucut, mucut, sigma, sigmas[ja]);
                  } // d

                  int const mb = sho_tools::nSHO(numaxs[ja]);
                  view2D<double> Vmat(nb, mb, 0.0);

                  // matrix multiply Vaux with overlap operator
                  stat += multiply_potential_matrix(Vmat, ovl, Vaux, numaxs[ia], numaxs[ja], lmax);

                  // display matrix
                  for(int ib = 0; ib < nb; ++ib) {
                      if (echo > 0) {
                          printf("# V ai#%i aj#%i %s ", ia, ja, labels(numaxs[ia],ib));
                          for(int jb = 0; jb < mb; ++jb) {
                              printf("%8.4f", Vmat[ib][jb]);
                          }   printf("\n");
                      } // echo
                  } // ib

                  if (1) { // extra: create overlap matrix
                      std::vector<sho_tools::SHO_index_t> idx(nb);
                      stat += sho_tools::construct_index_table<sho_tools::order_zyx>(idx.data(), numaxs[ia], echo);
                      std::vector<sho_tools::SHO_index_t> jdx(nb);
                      stat += sho_tools::construct_index_table<sho_tools::order_zyx>(jdx.data(), numaxs[ja], echo);
                      for(int ib = 0; ib < nb; ++ib) {        auto const i = idx[ib].Cartesian;
                          for(int jb = 0; jb < mb; ++jb) {    auto const j = jdx[jb].Cartesian;
                              Smat(ia,ja,ib,jb) = ovl(0,j.nx,i.nx) * ovl(1,j.ny,i.ny) * ovl(2,j.nz,i.nz);
                          } // jb
                      } // ib
                  } // echo extra
                  
              } // ja
          } // ia

          
          if (echo > 9) { // extra: display normalized overlap matrix
              for(int ia = 0; ia < natoms; ++ia) {        int const nb = sho_tools::nSHO(numaxs[ia]);
                  for(int ja = 0; ja < natoms; ++ja) {    int const mb = sho_tools::nSHO(numaxs[ja]);
                      if (echo > 1) printf("# ai#%i aj#%i\n", ia, ja);
                      for(int ib = 0; ib < nb; ++ib) {
                          printf("# S ai#%i aj#%i %s ", ia, ja, labels(numaxs[ia],ib));
                          for(int jb = 0; jb < mb; ++jb) {
//                               printf("%8.4f", Smat(ia,ja,ib,jb) / std::sqrt(Smat(ia,ia,ib,ib)*Smat(ja,ja,jb,jb)));
                              printf("%8.4f", Smat(ia,ja,ib,jb)); // un-normalized, all diagonal elements are 1.0
                          } // jb
                          printf("\n");
                      } // ib
                  } // ja
              } // ia
          } // echo extra
                  
          
          if (echo > 2) printf("\n# %s method=4 seems asymmetric!\n", __func__);
      } // scope: Method 3
     
      if (nullptr != xyzZ) delete[] xyzZ;
      return stat;
  } // test_potential_elements

  status_t all_tests(int const echo) {
    auto status = 0;
    status += test_potential_elements(echo); // expensive
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace sho_potential
