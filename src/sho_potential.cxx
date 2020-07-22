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
#include "real_space.hxx" // ::grid_t
#include "sho_tools.hxx" // ::nSHO, ::n1HO, ::order_*, ::SHO_index_t, ::construct_label_table
#include "sho_projection.hxx" // ::sho_project, ::sho_add
        // ::renormalize_coefficients
#include "boundary_condition.hxx" // Isolated_Boundary
#include "sho_overlap.hxx" // ::generate_product_tensor, ...
        // ::generate_overlap_matrix, ::generate_potential_tensor
#include "data_view.hxx" // view2D<T>, view3D<T>, view4D<T>
#include "linear_algebra.hxx" // ::inverse

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
  status_t generate_potential_matrix(view2D<real_t> & Vmat // result Vmat(i,j) = sum_m Vcoeff[m] * t(m,j,i) 
                            , view4D<real_t> const & t1D // input t1D(dir,m,j,i)
                            , real_t const Vcoeff[], int const numax_m // expansion of the potential in x^{mx}*y^{my}*z^{mz}
                            , int const numax_i, int const numax_j
                            , int const dir01=1) { // 1:direction dependent input tensor, 0:isotropic
      // use the expansion of the product of two Hermite Gauss functions into another one, factorized in 3D
      // can use different (dir01==1) tensors per direction or the same (dir01==0)
      // t(m,j,i) = t1D(0,m_x,j_x,i_x) * t1D(1,m_y,j_y,i_y) * t1D(2,m_z,j_z,i_z)

      int mzyx{0}; // contraction index
      for    (int mz = 0; mz <= numax_m;           ++mz) {
        for  (int my = 0; my <= numax_m - mz;      ++my) {
          for(int mx = 0; mx <= numax_m - mz - my; ++mx) {

            int izyx{0};
            for    (int iz = 0; iz <= numax_i;           ++iz) {
              for  (int iy = 0; iy <= numax_i - iz;      ++iy) {
                for(int ix = 0; ix <= numax_i - iz - iy; ++ix) {

                  int jzyx{0};
                  for    (int jz = 0; jz <= numax_j;           ++jz) {  auto const tz   = t1D(dir01*2,mz,jz,iz);
                    for  (int jy = 0; jy <= numax_j - jz;      ++jy) {  auto const tyz  = t1D(dir01*1,my,jy,iy) * tz;
                      for(int jx = 0; jx <= numax_j - jz - jy; ++jx) {  auto const txyz = t1D(      0,mx,jx,ix) * tyz;

                        Vmat(izyx,jzyx) += Vcoeff[mzyx] * txyz;

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

            ++mzyx;
          } // mx
        } // my
      } // mz
      assert( sho_tools::nSHO(numax_m) == mzyx );
    
      return 0;
  } // generate_potential_matrix

  
  template<typename real_t>
  status_t multiply_potential_matrix(view2D<real_t> & Vmat // result Vmat(i,j) = sum_k Vaux(i,k) * ovl(j,k)
     , view3D<real_t> const & ovl1D  // input ovl1D(dir,j,k)
     , view2D<real_t> const & Vaux // input Vaux(i,k)
     , int const numax_i, int const numax_j, int const numax_k) {

      for(int izyx = 0; izyx < sho_tools::nSHO(numax_i); ++izyx) {
        
          int jzyx{0};
          for    (int jz = 0; jz <= numax_j; ++jz) {
            for  (int jy = 0; jy <= numax_j - jz; ++jy) {
              for(int jx = 0; jx <= numax_j - jz - jy; ++jx) {

                  int kzyx{0}; // contraction index
                  for    (int kz = 0; kz <= numax_k; ++kz) {              auto const tz   = ovl1D(2,jz,kz);
                    for  (int ky = 0; ky <= numax_k - kz; ++ky) {         auto const tyz  = ovl1D(1,jy,ky) * tz;
                      for(int kx = 0; kx <= numax_k - kz - ky; ++kx) {    auto const txyz = ovl1D(0,jx,kx) * tyz;

                          Vmat(izyx,jzyx) += Vaux(izyx,kzyx) * txyz;

                          ++kzyx;
                      } // kx
                    } // ky
                  } // kz
                  assert( sho_tools::nSHO(numax_k) == kzyx );
                  
                  ++jzyx;
              } // jx
            } // jy
          } // jz
          assert( sho_tools::nSHO(numax_j) == jzyx );

      } // izyx
      
      return 0;
  } // multiply_potential_matrix


//   status_t normalize_coefficients(double coeff[], int const numax, double const sigma) {
//       return sho_projection::renormalize_coefficients(coeff, coeff, numax, sigma);
//   } // normalize_coefficients

  status_t normalize_potential_coefficients(double coeff[], int const numax, double const sigma, int const echo=0) {
      status_t stat = 0;
      int const nc = sho_tools::nSHO(numax);
      int const m = sho_tools::n1HO(numax);
      
      view2D<double> inv3D(nc, nc, 0.0); // get memory
      view2D<double> mat1D(m, m, 0.0); // get memory
      stat += sho_overlap::moment_normalization(mat1D.data(), mat1D.stride(), sigma, echo);

      int constexpr debug_check = 1;
      view2D<double> mat3D_copy(debug_check*nc, nc, 0.0); // get memory
      { // scope: set up mat3D
          view2D<double> mat3D(inv3D.data(), inv3D.stride()); // wrap
          int mzyx{0}; // moments
          for    (int mz = 0; mz <= numax; ++mz) {
            for  (int my = 0; my <= numax - mz; ++my) {
              for(int mx = 0; mx <= numax - mz - my; ++mx) {

                // mat1D(n,m) = <H(n)|x^m>
                // a specialty of the result matrix mat1D is that only matrix elements 
                //    connecting even-even or odd-odd indices are non-zero
                //    ToDo: this could be exploited in the following pattern
                int kzyx{0}; // Hermite coefficients
                for    (int kz = 0; kz <= numax; ++kz) {            auto const tz   = mat1D(kz,mz);
                  for  (int ky = 0; ky <= numax - kz; ++ky) {       auto const tyz  = mat1D(ky,my) * tz;
                    for(int kx = 0; kx <= numax - kz - ky; ++kx) {  auto const txyz = mat1D(kx,mx) * tyz;

                      // despite the name
                      mat3D(kzyx,mzyx) = txyz;
                      if (debug_check) mat3D_copy(kzyx,mzyx) = txyz;

                      ++kzyx;
                }}} assert( nc == kzyx );

                ++mzyx;
          }}} assert( nc == mzyx );

          // now invert the 3D matrix
          auto const stat = linear_algebra::inverse(nc, mat3D.data(), mat3D.stride());
          if (stat) { warn("Maybe factorization failed, status=%i", int(stat)); return stat; }
          // inverse is stored in inv3D due o pointer overlap
      } // scope

      std::vector<double> c_new(nc, 0.0); // get memory

      view2D<char> zyx_label;
      if (echo > 4) {
          printf("\n# %s numax=%i nc=%i sigma=%g %s\n", __func__, numax, nc, sigma*Ang, _Ang);
          zyx_label = view2D<char>(nc, 8);
          sho_tools::construct_label_table(zyx_label.data(), numax, sho_tools::order_zyx);
      } // echo

      // just a matrix-vector multiplication
      for(int mzyx = 0; mzyx < nc; ++mzyx) {
          double cc{0};
          for(int kzyx = 0; kzyx < nc; ++kzyx) {
              cc += inv3D(mzyx,kzyx) * coeff[kzyx];
          } // kzyx
          c_new[mzyx] = cc; // store
      } // mzyx
      
      if (debug_check) {
          // check if the input vector comes out again
          double dev{0};
          for(int kzyx = 0; kzyx < nc; ++kzyx) {
              double cc{0};
              for(int mzyx = 0; mzyx < nc; ++mzyx) {
                  cc += mat3D_copy(kzyx,mzyx) * c_new[mzyx];
              } // mzyx
              dev = std::max(dev, std::abs(cc - coeff[kzyx]));
          } // kzyx
          if (echo > 4) printf("# %s debug_check dev %.1e\n", __func__, dev);
      } // debug check

      if (echo > 4) printf("# %s\n\n", __func__);

      set(coeff, nc, c_new.data()); // copy the results back into the input array
      
      // Discussion:
      // if we, for some reason, have to reconstruct mat1D every time (although it only depends on sigma as sigma^m)
      // one could investigate if the factorization property stays
      // but probably not because of the range of indices is not a cube but a tetrahedron.
      // We could use a linear equation solver instead of matrix inverson,
      // however, the inversion only needs to be done once per (sigma,numax)
      // and also here the dependency on sigma can probably be moved out
      
      return stat;
  } // normalize_potential_coefficients
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_potential_elements(int const echo=5) {
      status_t stat = 0;
      
      auto const geo_file = control::get("geometry.file", "atoms.xyz");
      auto const vtotfile = control::get("sho_potential.test.vtot.filename", "vtot.dat"); // vtot.dat was written by spherical_atoms.
      int dims[] = {0, 0, 0};
      std::vector<double> vtot; // total smooth potential
      if (1) { // scope: read in the potential from a file
          std::ifstream infile(vtotfile);
          int npt = 0; 
          if (!infile.is_open()) {
              if (echo > 1) printf("# %s failed to open file %s\n",  __func__, vtotfile);
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
          if (echo > 3) printf("# %s use %i values from file %s\n", __func__, npt, vtotfile);
      } // scope
      
      double *xyzZ = nullptr;
      int natoms{0};
      double cell[3] = {0, 0, 0}; 
      int bc[3] = {-7, -7, -7};
      {
          stat += geometry_analysis::read_xyz_file(&xyzZ, &natoms, geo_file, cell, bc, 0);
          if (echo > 2) printf("# found %d atoms in file \"%s\" with cell=[%.3f %.3f %.3f] %s and bc=[%d %d %d]\n",
                              natoms, geo_file, cell[0]*Ang, cell[1]*Ang, cell[2]*Ang, _Ang, bc[0], bc[1], bc[2]);
      }
      
      
      
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

      auto const artificial_potential = int(control::get("sho_potential.test.artificial.potential", 0.));
      if (artificial_potential) { // scope: artificial linear potentials (other than constants)
          int const mx = (artificial_potential /   1) % 10;
          int const my = (artificial_potential /  10) % 10;
          int const mz = (artificial_potential / 100) % 10;
          if (echo > 0) printf("# artificial potential z^%i y^%i x^%i\n", mz, my, mx);
          for(int iz = 0; iz < g[2]; ++iz) {
              auto const z = iz*g.h[2] - origin[2];
              auto const zmz = intpow(z, mz);
              for(int iy = 0; iy < g[1]; ++iy) {
                  auto const y = iy*g.h[1] - origin[1];
                  auto const ymy = intpow(y, my);
                  for(int ix = 0; ix < g[0]; ++ix) {
                      auto const x = ix*g.h[0] - origin[0];
                      auto const xmx = intpow(x, mx);
                      
                      int const izyx = (iz*g[1] + iy)*g[0] + ix;
                      vtot[izyx] = xmx * ymy * zmz;

                  } // ix
              } // iy
          } // iz
      } // scope: artificial potentials

      auto const usual_numax = int(control::get("sho_potential.test.numax", 1.));
      auto const usual_sigma =     control::get("sho_potential.test.sigma", 2.);
      std::vector<int>    numaxs(natoms, usual_numax); // define SHO basis size
      std::vector<double> sigmas(natoms, usual_sigma); // define SHO basis spreads
      double const sigma_asymmetry = control::get("sho_potential.test.sigma.asymmetry", 1.0);
      if (sigma_asymmetry != 1) { sigmas[0] *= sigma_asymmetry; sigmas[natoms - 1] /= sigma_asymmetry; } // manipulate the spreads
//       --numaxs[natoms - 1]; // manipulate the basis size
      int numax_max{0};
      for(int ia = 0; ia < natoms; ++ia) {
          if (echo > 0) printf("# atom #%i Z=%g \tpos %9.3f %9.3f %9.3f  sigma=%g %s numax=%d\n", 
                ia, xyzZ[ia*4 + 3], xyzZ[ia*4 + 0]*Ang, xyzZ[ia*4 + 1]*Ang, xyzZ[ia*4 + 2]*Ang, 
                sigmas[ia]*Ang, _Ang, numaxs[ia]);
          numax_max = std::max(numax_max, numaxs[ia]);
      } // ia
      
//    auto const lmax = int(control::get("sho_potential.test.lmax", 2*numax_max)); // for Method4
      
      std::vector<view2D<char>> labels(16);
      for(int nu = 0; nu < 16; ++nu) {
          labels[nu] = view2D<char>(sho_tools::nSHO(nu), 8);
          sho_tools::construct_label_table(labels[nu].data(), nu, sho_tools::order_zyx);
      } // nu

      auto const method = int(control::get("sho_potential.test.method", -1.)); // bit-string, use method=7 or -1 to activate all
      
      int constexpr Numerical = 0, Between = 1, Onsite = 2, noReference = 3;
      view4D<double> SV_matrix[2][3]; // Numerical and Onsite (no extra Smat for Between)
      char const method_name[3][16] = {"numerical", "between", "onsite"};
      bool method_active[4] = {false, false, false, false};

      int const mxb = sho_tools::nSHO(numax_max);
      
      if (1 & method) { // scope:
          if (echo > 2) printf("\n# %s Method=1\n", __func__);
          // Method 1) fully numerical integration (expensive)
          //    for each pair of atoms and basis functions,
          //    add one basis function to an empty grid,
          //    multiply the potential, project with the other basis function
          std::vector<double>  basis(g.all(), 0.0);
          std::vector<double> Vbasis(g.all(), 0.0);
          auto & Smat = SV_matrix[0][Numerical];
          auto & Vmat = SV_matrix[1][Numerical];
          Vmat = view4D<double>(natoms, natoms, mxb, mxb, 0.0); // get memory
          Smat = view4D<double>(natoms, natoms, mxb, mxb, 0.0); // get memory
          view2D<double> Sdiag(natoms, mxb, 0.0);
          
          for(int ja = 0; ja < natoms; ++ja) {
              int const nbj = sho_tools::nSHO(numaxs[ja]);
              std::vector<double> coeff(nbj, 0.0);
              for(int jb = 0; jb < nbj; ++jb) {
                  set(coeff.data(), nbj, 0.0); coeff[jb] = 1; // Kronecker-delta
                  set(basis.data(), g.all(), 0.0); // clear
                  stat += sho_projection::sho_add(basis.data(), g, coeff.data(), numaxs[ja], center[ja], sigmas[ja], 0);

                  // multiply Vtot to the basis function
                  product(Vbasis.data(), basis.size(), basis.data(), vtot.data());

                  for(int ia = 0; ia < natoms; ++ia) {
                      int const nbi = sho_tools::nSHO(numaxs[ia]);
                      std::vector<double> Scoeff(nbi, 0.0);
                      std::vector<double> Vcoeff(nbi, 0.0);
                      stat += sho_projection::sho_project(Scoeff.data(), numaxs[ia], center[ia], sigmas[ia],  basis.data(), g, 0);
                      stat += sho_projection::sho_project(Vcoeff.data(), numaxs[ia], center[ia], sigmas[ia], Vbasis.data(), g, 0);
                      for(int ib = 0; ib < nbi; ++ib) {
                          Smat(ia,ja,ib,jb) = Scoeff[ib]; // copy result into large array
                          Vmat(ia,ja,ib,jb) = Vcoeff[ib]; // copy result into large array
                      } // ib
                      if (ia == ja) Sdiag(ja,jb) = Smat(ja,ja,jb,jb);
                  } // ia

              } // jb
          } // ja

          if (1) { // normalize with diagonal elements of the overlap matrix
                  for(int ia = 0; ia < natoms; ++ia) {        int const nbi = sho_tools::nSHO(numaxs[ia]);
                      for(int ja = 0; ja < natoms; ++ja) {    int const nbj = sho_tools::nSHO(numaxs[ja]);
                          for(int ib = 0; ib < nbi; ++ib) {
                              for(int jb = 0; jb < nbj; ++jb) {
                                  double const f = 1 / std::sqrt(Sdiag(ia,ib)*Sdiag(ja,jb));
                                  Vmat(ia,ja,ib,jb) *= f;
                                  Smat(ia,ja,ib,jb) *= f;
                              } // jb
                          } // ib
                      } // ja
                  } // ia
          } // normalize wih diagonal elements of the overlap matrix
          
          if (echo > 2) printf("\n# %s ToDo: check if method=1 depends on absolute positions!\n", __func__);
          
          method_active[Numerical] = true;
      } // scope: Method 1

      if (2 & method) { // scope:
          if (echo > 2) printf("\n# %s Method=2\n", __func__);
          // Method 2) analytical (cheap to compute)
          //    for each pair of atoms, find the center of weight,
          //    expand the potential in a SHO basis with sigma_V^-2 = sigma_1^-2 + sigma_2^-2
          //    and numax_V = numax_1 + numax_2 and determine the
          //    potential matrix elements using the tensor
          
          auto & Smat = SV_matrix[0][Between];
          auto & Vmat = SV_matrix[1][Between];
          Vmat = view4D<double>(natoms, natoms, mxb, mxb, 0.0); // get memory
          Smat = view4D<double>(natoms, natoms, mxb, mxb, 0.0); // get memory
          
          for(int ia = 0; ia < natoms; ++ia) {
              for(int ja = 0; ja < natoms; ++ja) {
                  double const alpha_i = 1/pow2(sigmas[ia]);
                  double const alpha_j = 1/pow2(sigmas[ja]);
                  double const sigma_V = 1/std::sqrt(alpha_i + alpha_j);
                  double const wi = alpha_i*pow2(sigma_V);
                  double const wj = alpha_j*pow2(sigma_V);
                  assert( std::abs( wi + wj - 1.0 ) < 1e-12 );
                  double cnt[3]; // center of weight
                  for(int d = 0; d < 3; ++d) {
                      cnt[d] = wi*xyzZ[ia*4 + d] + wj*xyzZ[ja*4 + d];
                  } // d
                  int const numax_V = numaxs[ia] + numaxs[ja];
                  if (echo > 1) printf("# ai#%i aj#%i  center of weight: %g %g %g sigma_V: %g %s numax_V=%i\n", ia, ja, 
                                          cnt[0]*Ang, cnt[1]*Ang, cnt[2]*Ang, sigma_V*Ang, _Ang, numax_V);
                  for(int d = 0; d < 3; ++d) {
                      cnt[d] += origin[d];
                  } // d
                  int const nucut_i = sho_tools::n1HO(numaxs[ia]);
                  int const nucut_j = sho_tools::n1HO(numaxs[ja]);

                  view4D<double> t(3, sho_tools::n1HO(numax_V), nucut_j, nucut_i, 0.0);
                  for(int d = 0; d < 3; ++d) {
                      auto const distance = center[ja][d] - center[ia][d];
                      stat += sho_overlap::moment_tensor(t[d].data(), distance, nucut_i, nucut_j, 
                                                                     sigmas[ia], sigmas[ja], numax_V);
                  } // d

                  std::vector<double> Vcoeff(sho_tools::nSHO(numax_V), 0.0);
                  stat += sho_projection::sho_project(Vcoeff.data(), numax_V, cnt, sigma_V, vtot.data(), g, 0); // 0:mute
                  
                  stat += normalize_potential_coefficients(Vcoeff.data(), numax_V, sigma_V, 0); // 0:mute
                  // now Vcoeff is represented w.r.t. powers of the Cartesian coords x^{nx}*y^{ny}*z^{nz}
#ifdef FULL_DEBUG
                  if (echo > 16) {
                      int mzyx{0};
                      for    (int mz = 0; mz <= numax_V; ++mz) {
                        for  (int my = 0; my <= numax_V - mz; ++my) {
                          for(int mx = 0; mx <= numax_V - mz - my; ++mx) {
                              auto const v = Vcoeff[mzyx];
                              if (ja <= ia && std::abs(v) > 5e-7)
                                  printf("# V_coeff ai#%i aj#%i %x%x%x %16.6f\n", ia, ja, mz,my,mx, v);
                              ++mzyx;
                          } // mx
                        } // my
                      } // mz
                      printf("\n");
                      assert(Vcoeff.size() == mzyx);
                  } // echo
#endif

                  // use the expansion of the product of two Hermite Gauss functions into another one
                  // Vmat(i,j) = sum_p Vcoeff[p] * t^3(p,j,i)
                  auto Vmat_iaja = Vmat(ia,ja);
                  stat += generate_potential_matrix(Vmat_iaja, t, Vcoeff.data(), numax_V, numaxs[ia], numaxs[ja]);
                  double const one[1] = {1.};
                  auto Smat_iaja = Smat(ia,ja);
                  stat += generate_potential_matrix(Smat_iaja, t, one, 0, numaxs[ia], numaxs[ja]);

              } // ja
          } // ia
        
          if (echo > 2) { printf("\n# %s method=2 seems symmetric!\n", __func__); fflush(stdout); }
          
          method_active[Between] = true;
      } // scope: Method 2

      if (4 & method) { // scope:
          if (echo > 2) printf("\n# %s Method=4\n", __func__);
          // Method 4) approximated
          //    for each atom expand the potential in a local SHO basis
          //    with spread sigma_V^2 = 2*sigma_1^2 at the atomic center,
          //    expand the other orbital in the local SHO basis (may need high lmax)
          //    also using the tensor.
          //    The matrix elements will not converge with the same speed w.r.t. lmax
          //    so we will require symmetrization
          auto & Smat = SV_matrix[0][Onsite];
          auto & Vmat = SV_matrix[1][Onsite];
          Vmat = view4D<double>(natoms, natoms, mxb, mxb, 0.0); // get memory
          Smat = view4D<double>(natoms, natoms, mxb, mxb, 0.0); // get memory

          auto const coarsest_grid_spacing = std::max(std::max(g.h[0], g.h[1]), g.h[2]);
          auto const highest_kinetic_energy = 0.5*pow2(constants::pi/coarsest_grid_spacing);
          auto const ekin = control::get("sho_potential.test.ekin", .03125); // specific for Method4
          auto const kinetic_energy = ekin * highest_kinetic_energy;

          if (echo > 3) printf("# grid spacing %g %s allows for kinetic energies up to %g %s, use %g %s (%.2f %%)\n",
              coarsest_grid_spacing*Ang, _Ang, highest_kinetic_energy*eV, _eV, kinetic_energy*eV, _eV, ekin*100);
          
          for(int ia = 0; ia < natoms; ++ia) {
              double const sigma_V = sigmas[ia]*std::sqrt(.5);
//               int const numax_V = 3*numaxs[ia]; // ToDo: is this the best choice? Better:
              // determine numax_V dynamically, depending on sigma_a and the grid spacing, (external parameter lmax is ignored)
              int const numax_V = std::floor(kinetic_energy*pow2(sigmas[ia]) - 1.5);
              if (echo > 5) printf("# atom #%i expand potential up to numax=%d with sigma=%g %s\n", ia, numax_V, sigma_V*Ang, _Ang);
              if (echo > 5) fflush(stdout);
              int const nbV = sho_tools::nSHO(numax_V);
              std::vector<double> Vcoeff(nbV, 0.0);

              // project the potential onto an auxiliary SHO basis centered at the position of atom ia
              sho_projection::sho_project(Vcoeff.data(), numax_V, center[ia], sigma_V, vtot.data(), g, 0);
              stat += normalize_potential_coefficients(Vcoeff.data(), numax_V, sigma_V, 0); // mute
              // now Vcoeff is in a representation w.r.t. powers of the Cartesian coords x^{mx}*y^{my}*z^{mz}
#if 1
              if (echo > 9) {
                  int mzyx{0};
                  for    (int mz = 0; mz <= numax_V; ++mz) {
                    for  (int my = 0; my <= numax_V - mz; ++my) {
                      for(int mx = 0; mx <= numax_V - mz - my; ++mx) {
                          auto const v = Vcoeff[mzyx];
                          if (std::abs(v) > 5e-7) printf("# V_coeff ai#%i %x%x%x %16.6f\n", ia, mz,my,mx, v);
//                        printf("#V_coeff_all %d %g ai#%i  %d %d %d\n", mz+my+mx, v, ia, mz,my,mx); // for plotting
                          ++mzyx;
                      } // mx
                    } // my
                  } // mz
                  printf("\n");
                  assert(Vcoeff.size() == mzyx);
              } // echo
#endif
      
              // determine the size of the auxiliary basis
              int const numax_k = numaxs[ia] + numax_V; // triangle rule
              if (echo > 5) printf("# atom #%i auxiliary basis up to numax=%d with sigma=%g %s\n", ia, numax_k, sigmas[ia]*Ang, _Ang);
              if (echo > 5) fflush(stdout);
              
              int const nbi = sho_tools::nSHO(numaxs[ia]);
              int const nbk = sho_tools::nSHO(numax_k);
              view2D<double> Vaux(nbi, nbk, 0.0); // get memory for Vaux(i,k)
              // now compute local matrix elements <local basis_i|V|large aux. basis_k>

              int const nucut_i = sho_tools::n1HO(numaxs[ia]);
              int const nucut_k = sho_tools::n1HO(numax_k);
              view4D<double> t(1, 1+numax_V, nucut_k, nucut_i, 0.0);
              double constexpr distance = 0.0;
              stat += sho_overlap::moment_tensor(t[0].data(), distance, nucut_i, nucut_k,
                                                 sigmas[ia], sigmas[ia], numax_V);

              // Vaux(i,k) = sum_p Vcoeff[m] * t(m,k,i)
              int constexpr isotropic = 0;
              stat += generate_potential_matrix(Vaux, t, Vcoeff.data(), numax_V, numaxs[ia], numax_k, isotropic);
#ifdef DEVEL
              if (echo > 17) {
                  printf("\n# Vaux for the atom #%i:\n", ia);
                  printf("# i-legend   ");
                  for(int izyx = 0; izyx < nbi; ++izyx) {
                      printf(" %-7s", (numaxs[ia] > 15) ? "???" : labels[numaxs[ia]][izyx]);
                  }   printf("\n");
                  for(int kzyx = 0; kzyx < nbk; ++kzyx) {
                      printf("# Vaux %s  ", (numax_k > 15) ? "???" : labels[numax_k][kzyx]);
                      for(int izyx = 0; izyx < nbi; ++izyx) {
                          printf(" %7.4f", Vaux(izyx,kzyx));
                      } // i
                      printf("\n");
                  } // k
                  printf("\n");
              } // echo
#endif

              for(int ja = 0; ja < natoms; ++ja) {
                  if (echo > 9) printf("# ai#%i aj#%i\n", ia, ja);

                  int const nucut_j = sho_tools::n1HO(numaxs[ja]);
                  view3D<double> ovl(3, nucut_j, nucut_k); // get memory, index order (dir,j,k)
                  for(int d = 0; d < 3; ++d) {
                      auto const distance = center[ja][d] - center[ia][d];
                      sho_overlap::generate_overlap_matrix(ovl[d].data(), distance, nucut_k, nucut_j, sigmas[ia], sigmas[ja]);
#ifdef DEVEL
                      if (echo > 19) {
                          printf("\n# ovl for the %c-direction:\n", 'x' + d);
                          for(int j = 0; j < nucut_j; ++j) {
                              printf("# ovl n%c=%x  ", 'x' + d, j);
                              for(int k = 0; k < nucut_k; ++k) {
                                  printf("%8.4f", ovl(d,j,k));
                              } // k
                              printf("\n");
                          } // j 
                          printf("\n");
                      } // echo
#endif
                  } // d

                  int const nbj = sho_tools::nSHO(numaxs[ja]);
                  view2D<double> Vmat_iaja(nbi, nbj, 0.0); // get memory and initialize

                  // matrix multiply Vaux with overlap operator
                  // Vmat_iaja(i,j) = sum_k Vaux(i,k) * ovl(j,k)
                  stat += multiply_potential_matrix(Vmat_iaja, ovl, Vaux, numaxs[ia], numaxs[ja], numax_k);

                  { // extra: create overlap matrix
                      std::vector<sho_tools::SHO_index_t> idx(nbi), jdx(nbj);
                      stat += sho_tools::construct_index_table<sho_tools::order_zyx>(idx.data(), numaxs[ia], echo);
                      stat += sho_tools::construct_index_table<sho_tools::order_zyx>(jdx.data(), numaxs[ja], echo);
                      for(int ib = 0; ib < nbi; ++ib) {        auto const i = idx[ib].Cartesian;
                          for(int jb = 0; jb < nbj; ++jb) {    auto const j = jdx[jb].Cartesian;
                              Smat(ia,ja,ib,jb) = ovl(0,j.nx,i.nx) * ovl(1,j.ny,i.ny) * ovl(2,j.nz,i.nz);
                              Vmat(ia,ja,ib,jb) = Vmat_iaja(ib,jb);
                          } // jb
                      } // ib
                  } // echo extra

              } // ja
          } // ia
          
          if (echo > 2) printf("\n# %s method=4 seems asymmetric!\n", __func__);
          method_active[Onsite] = true;
      } // scope: Method 4






      int ref_method[3] = {noReference, noReference, noReference}; method_active[noReference] = false;
      if (method_active[Numerical]) {
          if (method_active[Between]) ref_method[Between] = Numerical;
          if (method_active[Onsite])  ref_method[Onsite]  = Numerical;
      } else {
          if (method_active[Between] && method_active[Onsite]) ref_method[Onsite] = Between;
      }

      // now display all of these methods interleaved
      if (echo > 0) {
          auto const sv_start = int(control::get("sho_potential.test.show.potential.only", 0.)); // 0:show S and V, 1: V only
          for(int sv = sv_start; sv < 2; ++sv) {
              auto const sv_char = sv?'V':'S';
              if (echo > 7) printf("\n# %s (%c)\n", sv?"potential":"overlap", sv_char);
              
              double max_abs_dev[][2] = {{0, 0}, {0, 0}, {0, 0}};
              for(int ia = 0; ia < natoms; ++ia) {
                  int const nbi = sho_tools::nSHO(numaxs[ia]);
                  for(int ja = 0; ja < natoms; ++ja) {
                      int const nbj = sho_tools::nSHO(numaxs[ja]);
                      if (echo > 1) {
                          printf("\n# ai#%i aj#%i        ", ia, ja);
                          for(int jb = 0; jb < nbj; ++jb) {
                              printf(" %-7s", (numaxs[ja] > 15) ? "???" : labels[numaxs[ja]][jb]); // column legend
                          } // jb
                          printf("\n");
                      } // echo
                      for(int ib = 0; ib < nbi; ++ib) {
                          for(int m = 0; m < 3; ++m) { // method
                              if (method_active[m]) {
                                  printf("# %c ai#%i aj#%i %s ", sv_char, ia, ja, (numaxs[ia] > 15) ? "???" : labels[numaxs[ia]][ib]);
                                  double abs_dev{0};
                                  int const rm = ref_method[m]; // reference method
                                  for(int jb = 0; jb < nbj; ++jb) {
                                      auto const v = SV_matrix[sv][m](ia,ja,ib,jb);
                                      printf(" %7.4f", v);
                                      auto const ref = method_active[rm] ? SV_matrix[sv][rm](ia,ja,ib,jb) : 0;
                                      auto const d = v - ref;
                                      abs_dev = std::max(std::abs(d), abs_dev);
                                  } // jb
                                  printf(" %s", method_name[m]);
                                  if (noReference != rm) printf(", dev=%.1e", abs_dev);
                                  printf("\n");
                                  int const diag = (ia == ja);
                                  max_abs_dev[m][diag] = std::max(max_abs_dev[m][diag], abs_dev);
                              } // active?
                          } // method
                      } // ib
                  } // ja
              } // ia

              for(int m = 0; m < 3; ++m) { // method
                  if (method_active[m]) {
                      int const rm = ref_method[m]; // reference method
                      if (noReference != rm) {
                          printf("\n# %c largest abs deviation of '%s' to '%s' is (diag) %g and %g (off-diag), pot=%03d\n",
                              sv_char, method_name[m], method_name[rm], max_abs_dev[m][1], max_abs_dev[m][0], artificial_potential);
                      } // no reference
                  } // active
              } // method
              
          } // sv
      } // echo
      
      
      if (nullptr != xyzZ) delete[] xyzZ;
      return stat;
  } // test_potential_elements

  status_t all_tests(int const echo) {
    status_t status(0);
    status += test_potential_elements(echo); // expensive
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace sho_potential
