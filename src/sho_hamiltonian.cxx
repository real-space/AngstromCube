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
#include "sho_projection.hxx" // ::sho_project
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
  
  template<typename real_t>
  status_t kinetic_matrix(view2D<real_t> & Tmat // result Tmat(i,j) = 
                       , view4D<double> const & t1D // input t1D(dir,0,i,j) nabla^2 operator
                       , view4D<double> const & o1D // input o1D(dir,0,i,j) overlap operator
                       , int const numax_i, int const numax_j
                       , double const f=0.5) { // typical prefactor of the kinetic energy in Hartree atomic units

      int izyx{0};
      for    (int iz = 0; iz <= numax_i;           ++iz) {
        for  (int iy = 0; iy <= numax_i - iz;      ++iy) {
          for(int ix = 0; ix <= numax_i - iz - iy; ++ix) {

            int jzyx{0};
            for    (int jz = 0; jz <= numax_j;           ++jz) { auto const Tz = t1D(2,0,iz,jz), oz = o1D(2,0,iz,jz);
              for  (int jy = 0; jy <= numax_j - jz;      ++jy) { auto const Ty = t1D(1,0,iy,jy), oy = o1D(1,0,iy,jy);
                for(int jx = 0; jx <= numax_j - jz - jy; ++jx) { auto const Tx = t1D(0,0,ix,jx), ox = o1D(0,0,ix,jx);

                  Tmat(izyx,jzyx) = f * ( Tx*oy*oz + ox*Ty*oz + ox*oy*Tz );

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

      return 0;
  } // kinetic_matrix
  
  
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
      double *xyzZ_m = nullptr;
      int natoms{0};
      double cell[3] = {0, 0, 0}; 
      int bc[3] = {-7, -7, -7};
      { // scope: read atomic positions
          stat += geometry_analysis::read_xyz_file(&xyzZ_m, &natoms, geo_file, cell, bc, 0);
          if (echo > 2) printf("# found %d atoms in file \"%s\" with cell=[%.3f %.3f %.3f] %s and bc=[%d %d %d]\n",
                              natoms, geo_file, cell[0]*Ang, cell[1]*Ang, cell[2]*Ang, _Ang, bc[0], bc[1], bc[2]);
      } // scope
      view2D<double> xyzZ(xyzZ_m, 4); // wrap for simpler usage

//    for(int d = 0; d < 3; ++d) assert(bc[d] == Isolated_Boundary); // ToDo: implement periodic images

      real_space::grid_t g(dims);
      g.set_grid_spacing(cell[0]/g[0], cell[1]/g[1], cell[2]/g[2]);
      if (echo > 1) printf("# use  %g %g %g %s grid spacing\n", g.h[0]*Ang, g.h[1]*Ang, g.h[2]*Ang, _Ang);
      if (echo > 1) printf("# cell is  %g %g %g %s\n", g.h[0]*g[0]*Ang, g.h[1]*g[1]*Ang, g.h[2]*g[2]*Ang, _Ang);
      double const origin[] = {.5*(g[0] - 1)*g.h[0],
                               .5*(g[1] - 1)*g.h[1], 
                               .5*(g[2] - 1)*g.h[2]};


      auto const usual_numax = int(control::get("sho_hamiltonian.test.numax", 1.));
      auto const usual_sigma =     control::get("sho_hamiltonian.test.sigma", 2.);
      std::vector<int>    numaxs(natoms, usual_numax); // define SHO basis sizes
      std::vector<double> sigmas(natoms, usual_sigma); // define SHO basis spreads
      double const sigma_asymmetry = control::get("sho_hamiltonian.test.sigma.asymmetry", 1.0);
      if (sigma_asymmetry != 1) { sigmas[0] *= sigma_asymmetry; sigmas[natoms - 1] /= sigma_asymmetry; } // manipulate the spreads
//       ++numaxs[0]; // manipulate to test variable block sizes
      
      
      int const ncenters = pow2(natoms); // in principle, we only need a unique set of centers, even without periodic images, only ((natoms + 1)*natoms)/2
      view2D<double> center(ncenters, 8, 0.0); // list of potential expansion centers
      view3D<int> center_map(natoms, natoms, 2, -1);
      { // scope: set up list of centers
          int ic{0};
          for(int ia = 0; ia < natoms; ++ia) {
              for(int ja = 0; ja < natoms; ++ja) {
                  double const alpha_i = 1/pow2(sigmas[ia]);
                  double const alpha_j = 1/pow2(sigmas[ja]);
                  double const sigma_V = 1/std::sqrt(alpha_i + alpha_j);
                  int    const numax_V = numaxs[ia] + numaxs[ja];
                  double const wi = alpha_i*pow2(sigma_V);
                  double const wj = alpha_j*pow2(sigma_V);
                  assert( std::abs( wi + wj - 1.0 ) < 1e-12 );
                  if (echo > 1) printf("# ai#%i aj#%i \tcenter of weight\t", ia, ja);
                  for(int d = 0; d < 3; ++d) {
                      center(ic,d) = wi*xyzZ(ia,d) + wj*xyzZ(ja,d);
                      if (echo > 1) printf("%12.6f", center(ic,d)*Ang);
                      center(ic,d) += origin[d];
                  } // d
                  if (echo > 1) printf("  numax_V= %i sigma_V: %g %s \n", numax_V, sigma_V*Ang, _Ang);
                  center(ic,3) = sigma_V;
                  center(ic,4) = numax_V;
                  // debug info
                  center(ic,5) = 0;  // ToDo: store also the index of the periodic image, once active
                  center(ic,6) = ia; // ToDo: use global atom indices here
                  center(ic,7) = ja; // ToDo: use global atom indices here
                  
                  center_map(ia,ja,0) = ic;
                  center_map(ia,ja,1) = numax_V;

                  ++ic;
              } // ja
          } // ia
          assert( ncenters == ic );
      } // scope: set up list of centers

      // perform the projection of the local potential
      std::vector<std::vector<double>> Vcoeffs(ncenters);
      double const scale_potential = control::get("sho_hamiltonian.scale.potential", 1.0);
      if (1.0 != scale_potential) warn("scale potential by %g", scale_potential);
      for(int ic = 0; ic < ncenters; ++ic) {
          double const sigma_V = center(ic,3);
          int    const numax_V = center(ic,4);
          Vcoeffs[ic] = std::vector<double>(sho_tools::nSHO(numax_V), 0.0);
          stat += sho_projection::sho_project(Vcoeffs[ic].data(), numax_V, center[ic], sigma_V, vtot.data(), g, 0); // 0:mute
          
          stat += sho_potential::normalize_potential_coefficients(Vcoeffs[ic].data(), numax_V, sigma_V, 0); // 0:mute
          
          scale(Vcoeffs[ic].data(), Vcoeffs[ic].size(), scale_potential);
      } // ic
      

      std::vector<int> offset(natoms + 1, 0); // the basis atoms localized on atom#i start at offset[i]
      int total_basis_size{0}, maximum_numax{-1};
      for(int ia = 0; ia < natoms; ++ia) {
          if (echo > 0) {
              printf("# atom#%i \tZ=%g \tposition %12.6f%12.6f%12.6f  numax= %d sigma= %.3f %s\n",
                  ia, xyzZ(ia,3), xyzZ(ia,0)*Ang, xyzZ(ia,1)*Ang, xyzZ(ia,2)*Ang, numaxs[ia], sigmas[ia]*Ang, _Ang);
          } // echo
          int const atom_basis_size = sho_tools::nSHO(numaxs[ia]);
          maximum_numax = std::max(maximum_numax, numaxs[ia]);
          offset[ia] = total_basis_size; // prefix sum
          total_basis_size += atom_basis_size;
      } // ia
      int const numax_max = maximum_numax;
      if (echo > 5) printf("# largest basis size per atom is %d, numax=%d\n", sho_tools::nSHO(numax_max), numax_max);
      offset[natoms] = total_basis_size;
      int const nB   = total_basis_size;
      int const nBa  = align<4>(nB); // memory aligned main matrix stride
      
      typedef double psi_t; // can also be std::complex<double>
      view3D<psi_t> SHW(3, nB, nBa, psi_t(0)); // get memory for Overlap, Hamiltonian and Work array


      //
      // now construct the Hamiltonian:
      //   -- kinetic energy,             c.f. sho_overlap::test_simple_crystal
      //   -- local potential,            c.f. sho_potential::test_local_potential_matrix_elements (method=2)
      //   -- non-local PAW Hamiltonian contributions
      // and the overlap operator:
      //   -- SHO basis function overlap, c.f. sho_overlap::test_simple_crystal
      //   -- non-local PAW charge-deficit contributions
      //
      
      // prefactor of kinetic energy in Hartree atomic units
      double const kinetic = control::get("sho_hamiltonian.scale.kinetic", 1.0) * 0.5;
      if (0.5 != kinetic) warn("kinetic energy prefactor is %g", kinetic);
      double const ones[1] = {1.0}; // expansion of the identity (constant==1) into x^{m_x} y^{m_y} z^{m_z}

      std::vector<std::vector<view2D<psi_t>>> S_iaja(natoms); // construct sub-views for each atom pair
      std::vector<std::vector<view2D<psi_t>>> H_iaja(natoms); // construct sub-views for each atom pair
      for(int ia = 0; ia < natoms; ++ia) {
          S_iaja[ia].resize(natoms);
          H_iaja[ia].resize(natoms);
          int const n1i   = sho_tools::n1HO(numaxs[ia]);
          for(int ja = 0; ja < natoms; ++ja) {
              S_iaja[ia][ja] = view2D<psi_t>(&(SHW(0,offset[ia],offset[ja])), nBa); // wrapper to sub-blocks of the overlap matrix
              H_iaja[ia][ja] = view2D<psi_t>(&(SHW(1,offset[ia],offset[ja])), nBa); // wrapper to sub-blocks of the overlap matrix
              int const n1j   = sho_tools::n1HO(numaxs[ja]);
#if 0
              int const nb_ja = sho_tools::nSHO(numaxs[ja]);
              int const nb_ia = sho_tools::nSHO(numaxs[ia]);
              for(int ib = 0; ib < nb_ia; ++ib) {
                  for(int jb = 0; jb < nb_ja; ++jb) {
                      S_iaja[ia][ja](ib,jb) = 10*ia + ja + ib*.1 + jb*.01;           // dummy values to see if all matrix entries are set
                      H_iaja[ia][ja](ib,jb) = (ia + (ib + 1)*.1)*(ja + (jb + 1)*.1); // dummy values to see if all matrix entries are set
                  } // jb
              } // ib
#endif

              int const ic = center_map(ia,ja,0); // expansion center index
              int const numax_V = center_map(ia,ja,1); // expansion of the local potential into x^{m_x} y^{m_y} z^{m_z} around a given expansion center
              int const maxmoment = std::max(0, numax_V);

              view4D<double> nabla2(3, 1, n1i + 1, n1j + 1, 0.0);      //  <\chi1D_i|d/dx  d/dx|\chi1D_j>
              view4D<double> ovl1Dm(3, 1 + maxmoment, n1i, n1j, 0.0);  //  <\chi1D_i| x^moment |\chi1D_j>
              for(int d = 0; d < 3; ++d) { // spatial directions x,y,z
                  double const distance = xyzZ(ia,d) - xyzZ(ja,d); // does not account for periodic images, ToDo
                  stat += sho_overlap::nabla2_matrix(nabla2[d].data(), distance, n1i + 1, n1j + 1, sigmas[ia], sigmas[ja]);
                  stat += sho_overlap::moment_tensor(ovl1Dm[d].data(), distance, n1i, n1j, sigmas[ia], sigmas[ja], maxmoment);
              } // d

              // construct the overlap matrix of SHO basis functions
              // Smat(i,j) := ovl_x(ix,jx) * ovl_y(iy,jy) * ovl_z(iz,jz)
              stat += sho_potential::potential_matrix(S_iaja[ia][ja], ovl1Dm, ones, 0, numaxs[ia], numaxs[ja]);

              // construct the kinetic energy contribution
              stat += sho_hamiltonian::kinetic_matrix(H_iaja[ia][ja], nabla2, ovl1Dm, numaxs[ia], numaxs[ja], kinetic);

              // add the contribution of the local potential
              stat += sho_potential::potential_matrix(H_iaja[ia][ja], ovl1Dm, Vcoeffs[ic].data(), numax_V, numaxs[ia], numaxs[ja]);

          } // ja
      } // ia
      Vcoeffs.clear(); // free memory for potential coefficients

      // prepare for the PAW contributions: find the projection coefficient matrix P = <\chi3D_{ia ib}|\tilde p_{ja jb}>
      int const natoms_PAW = natoms; // keep it flexible, ToDo: if we turn it on, sometimes there is an exception:
      // "Incorrect checksum for freed object 0xsomeaddress: probably modified after being freed"
      auto const xyzZ_PAW = view2D<double>(xyzZ.data(), xyzZ.stride()); // duplicate view
      std::vector<int>    numax_PAW(natoms_PAW,  3); // ToDo: match with input from atomic PAW configuration
      std::vector<double> sigma_PAW(natoms_PAW, .5); // ToDo: match with input from atomic PAW configuration
      std::vector<view3D<double>> SH_PAW(natoms_PAW);
      for(int ka = 0; ka < natoms_PAW; ++ka) {
          int const nb_ka = sho_tools::nSHO(numax_PAW[ka]); // number of projectors
          SH_PAW[ka] = view3D<double>(2, nb_ka, nb_ka, 0.0); // get memory and initialize
          // ToDo: get atomic charge-deficit matrix into SH_PAW[ka][0]
          // ToDo: get atomic Hamiltonian    matrix into SH_PAW[ka][1]
      } // ka

      std::vector<std::vector<view2D<double>>> P_iaka(natoms); // potentially sparse lists, e.g. compressed row format
      std::vector<std::vector<view3D<double>>> Psh_iala(natoms); // potentially sparse lists, e.g. compressed row format
      for(int ia = 0; ia < natoms; ++ia) {
          P_iaka[ia].resize(natoms_PAW);
          Psh_iala[ia].resize(natoms_PAW);
          int const nb_ia = sho_tools::nSHO(numaxs[ia]);
          int const n1i   = sho_tools::n1HO(numaxs[ia]);
          for(int ka = 0; ka < natoms_PAW; ++ka) {
              int const nb_ka = sho_tools::nSHO(numax_PAW[ka]);
              int const n1k   = sho_tools::n1HO(numax_PAW[ka]);
              P_iaka[ia][ka] = view2D<double>(nb_ia, nb_ka, 0.0); // get memory and initialize

              view4D<double> ovl1D(3, 1, n1i, n1k, 0.0);  //  <\chi1D_i|\chi1D_k>
              for(int d = 0; d < 3; ++d) { // spatial directions x,y,z
                  double const distance = xyzZ(ia,d) - xyzZ_PAW(ka,d); // does not account for periodic images, ToDo
                  stat += sho_overlap::overlap_matrix(ovl1D(d,0,0), distance, n1i, n1k, sigmas[ia], sigma_PAW[ka]);
              } // d

              double const ones[1] = {1.0}; // expansion of the identity
              // P(ia,ka,i,j) := ovl_x(ix,kx) * ovl_y(iy,ky) * ovl_z(iz,kz)
              stat += sho_potential::potential_matrix(P_iaka[ia][ka], ovl1D, ones, 0, numaxs[ia], numax_PAW[ka]);

              
              // multiply P from left to SH_PAW (block diagonal --> ka == la)
              auto const la = ka;
              int const nb_la = nb_ka;
              Psh_iala[ia][la] = view3D<double>(2, nb_ia, nb_la, 0.0); // get memory and initialize
              for(int ib = 0; ib < nb_ia; ++ib) {
                  for(int lb = 0; lb < nb_la; ++lb) {
                      double s{0}, h{0};
                      for(int kb = 0; kb < nb_ka; ++kb) { // contract
                          s += P_iaka[ia][ka](ib,kb) * SH_PAW[ka](0,kb,lb);
                          h += P_iaka[ia][ka](ib,kb) * SH_PAW[ka](1,kb,lb);
                      } // kb
                      Psh_iala[ia][la](0,ib,lb) = s;
                      Psh_iala[ia][la](1,ib,lb) = h;
                  } // lb
              } // ib

          } // ka
      } // ia

      double const scale_h = control::get("sho_hamiltonian.scale.nonlocal.h", 1.0);
      double const scale_s = control::get("sho_hamiltonian.scale.nonlocal.s", 1.0);
      if (1 != scale_h || 1 != scale_s) warn("scale PAW contributions to H and S by %g and %g, respectively", scale_h, scale_s);

      // PAW contributions to H_{ij} = P_{ik} h_{kl} P_{jl} = Ph_{il} P_{jl}
      //                  and S_{ij} = P_{ik} s_{kl} P_{jl} = Ps_{il} P_{jl}
      for(int ja = 0; ja < natoms; ++ja) {
          int const nb_ja = sho_tools::nSHO(numaxs[ja]);
          for(int la = 0; la < natoms_PAW; ++la) { // contract
              int const nb_la = sho_tools::nSHO(numax_PAW[la]);
              for(int ia = 0; ia < natoms; ++ia) {
                  int const nb_ia = sho_tools::nSHO(numaxs[ia]);
                  for(int ib = 0; ib < nb_ia; ++ib) {
                      for(int jb = 0; jb < nb_ja; ++jb) {
                          double s{0}, h{0};
                          for(int lb = 0; lb < nb_la; ++lb) { // contract
                              s += Psh_iala[ia][la](0,ib,lb) * P_iaka[ja][la](jb,lb);
                              h += Psh_iala[ia][la](1,ib,lb) * P_iaka[ja][la](jb,lb);
                          } // lb
                          S_iaja[ia][ja](ib,jb) += s * scale_s;
                          H_iaja[ia][ja](ib,jb) += h * scale_h;
                      } // jb
                      
                  } // ib
              } // ia
          } // la
      } // ja

      
      auto const ovl_eig = int(control::get("sho_hamiltonian.test.overlap.eigvals", 0.));
      std::vector<double> eigvals(nB, 0.0);
      
      for(int s0h1 = 0; s0h1 < 2; ++s0h1) {
          if (echo > 0) printf("\n\n");
          auto const matrix_name = s0h1 ? "Hamiltonian" : "overlap";
          auto const  u = s0h1 ?  eV :  1 ;
          auto const _u = s0h1 ? _eV : "";
          if (echo > 9 - s0h1) { // display S and H
              printf("\n# %s matrix (%s):\n", matrix_name, _u);
              for(int iB = 0; iB < nB; ++iB) {
                  printf("# row%3i ", iB);
                  for(int jB = 0; jB < nB; ++jB) {
                      printf("%8.3f", SHW(s0h1,iB,jB)*u);
                  } // jB
                  printf("\n");
              } // iB
              printf("\n");
          } // echo

          status_t stat_eig(0);
          if ((0 == s0h1) && ovl_eig) {
              set(SHW(2,0), nB*nBa, SHW(0,0)); // copy overlap matrix S into work array W
              stat_eig = linear_algebra::eigenvalues(nB, SHW(2,0), nBa, eigvals.data());
          } else {
              stat_eig = linear_algebra::generalized_eigenvalues(nB, SHW(1,0), nBa, SHW(0,0), nBa, eigvals.data());
          } // ovl_eig

          if ((1 == s0h1) || ovl_eig) {
              if (stat_eig) {
                  warn("diagonalizing the %s matrix failed, status= %i", matrix_name, int(stat_eig));
                  stat += stat_eig;
              } else if (nB > 0) {
                  double const lowest_eigenvalue = eigvals[0], highest_eigenvalue = eigvals[nB - 1];
                  if (echo > 2) {
                      printf("# eigenvalues of the %s matrix: ", matrix_name);
                      int constexpr mB = 8; // show at most the 6 lowest + 2 highest eigenvalues
                      for(int iB = 0; iB < std::min(nB - 2, mB - 2); ++iB) {
                          printf(" %g", eigvals[iB]*u);
                      } // iB
                      if (nB > mB) printf(" ..."); // there are more eigenavalues than we display
                      printf(" %g %g %s\n", eigvals[nB - 2]*u, eigvals[nB - 1]*u, _u); // last two
                  } // echo
                  if (echo > 4) printf("# lowest and highest eigenvalue of the %s matrix are %g and %g %s, respectively\n", 
                                            matrix_name, lowest_eigenvalue*u, highest_eigenvalue*u, _u);
                  if (s0h1 == 0 && lowest_eigenvalue < .1) warn("overlap matrix has small eigenvalues, lowest= %g", lowest_eigenvalue);
              } // stat_eig
          } // ovl_eig

      } // s0h1
      
      if (nullptr != xyzZ_m) delete[] xyzZ_m;
      return stat;
  } // test_Hamiltonian

  status_t all_tests(int const echo) {
    status_t status(0);
    status += test_Hamiltonian(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace sho_hamiltonian
