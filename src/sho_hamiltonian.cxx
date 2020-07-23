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

      auto const center = new double[natoms][4]; // list of atomic centers
      for(int ia = 0; ia < natoms; ++ia) {
          for(int d = 0; d < 3; ++d) {
              center[ia][d] = xyzZ(ia,d) + origin[d]; // w.r.t. to the center of grid point (0,0,0)
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
                  ia, xyzZ(ia,3), xyzZ(ia,0)*Ang, xyzZ(ia,1)*Ang, xyzZ(ia,2)*Ang, numaxs[ia], sigmas[ia]*Ang, _Ang);
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


      //
      // now construct the Hamiltonian:
      //   -- kinetic energy,             c.f. sho_overlap::test_simple_crystal
      //   -- local potential,            c.f. sho_potential::test_local_potential_matrix_elements (method=2)
      //   -- non-local PAW Hamiltonian contributions
      // and the overlap operator:
      //   -- SHO basis function overlap, c.f. sho_overlap::test_simple_crystal
      //   -- non-local PAW charge-deficit contributions
      //
      
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
              int const maxmoment = 0; // ToDo: adjust to the expansion of the local potential
              view4D<double> nabla2(3, 1, n1i + 1, n1j + 1, 0.0);      //  <\chi1D_i|d/dx  d/dx|\chi1D_j>
              view4D<double> ovl1Dm(3, 1 + maxmoment, n1i, n1j, 0.0);  //  <\chi1D_i| x^moment |\chi1D_j>
              for(int d = 0; d < 3; ++d) { // spatial directions x,y,z
                  double const distance = xyzZ(ia,d) - xyzZ(ja,d); // does not account for periodic images, ToDo
                  stat += sho_overlap::nabla2_matrix(nabla2[d].data(), distance, n1i + 1, n1j + 1, sigmas[ia], sigmas[ja]);
                  stat += sho_overlap::moment_tensor(ovl1Dm[d].data(), distance, n1i, n1j, sigmas[ia], sigmas[ja], maxmoment);
              } // d
              
              // construct the overlap matrix of SHO basis functions
              double const ones[1] = {1.0}; // expansion of the identity (constant==1) into x^{m_x} y^{m_y} z^{m_z}
              // Smat(i,j) := ovl_x(ix,jx) * ovl_y(iy,jy) * ovl_z(iz,jz)
              stat += sho_potential::potential_matrix(S_iaja[ia][ja], ovl1Dm, ones, 0, numaxs[ia], numaxs[ja]);

              // construct the kinetic energy contribution
              double const kinetic[1] = {0.5}; // prefactor of kinetic energy in Hartree atomic units
              stat += sho_potential::potential_matrix(H_iaja[ia][ja], nabla2, kinetic, 0, numaxs[ia], numaxs[ja]);
              
              // add the contribution of the local potential
              double const Vcoeff[1] = {0.0}; // prelim, ToDo: fill with real values from expansion
              int const numax_V = 0; // expansion of the local potential into x^{m_x} y^{m_y} z^{m_z} around a given expansion center
              stat += sho_potential::potential_matrix(H_iaja[ia][ja], ovl1Dm, Vcoeff, numax_V, numaxs[ia], numaxs[ja]);

          } // ja
      } // ia


      // prepare for the PAW contributions: find the projection coefficient matrix P = <\chi3D_{ia ib}|\tilde p_{ja jb}>
      int const natoms_PAW = natoms; // keep it flexible, ToDo: if we turn it on, sometimes there is an exception:
      // "Incorrect checksum for freed object 0xsomeaddress: probably modified after being freed"
      auto const xyzZ_PAW = view2D<double>(xyzZ.data(), xyzZ.stride()); // duplicate view
      std::vector<int>    numax_PAW(natoms_PAW,  3); // ToDo: match with input from atomic PAW configuration
      std::vector<double> sigma_PAW(natoms_PAW, .5); // ToDo: match with input from atomic PAW configuration
      std::vector<view3D<double>> SH_PAW(natoms_PAW);
      for(int ka = 0; ka < natoms_PAW; ++ka) {
          int const np_ka = sho_tools::nSHO(numax_PAW[ka]); // number of projectors
          SH_PAW[ka] = view3D<double>(2, np_ka, np_ka, 0.0); // get memory and initialize
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
                  stat += sho_overlap::overlap_matrix(ovl1D[d].data(), distance, n1i + 1, n1k + 1, sigmas[ia], sigma_PAW[ka]);
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
                          S_iaja[ia][ja](ib,jb) += s;
                          H_iaja[ia][ja](ib,jb) += h;
                      } // jb
                      
                  } // ib
              } // ia
          } // la
      } // ja

      if (echo > 7) { // display S and H
          for(int s0h1 = 0; s0h1 < 2; ++s0h1) {
              printf("\n# %s matrix:\n", s0h1 ? "Hamiltonian" : "overlap");
              for(int iB = 0; iB < nB; ++iB) {
                  printf("# row%3i ", iB);
                  for(int jB = 0; jB < nB; ++jB) {
                      printf("%8.3f", SHW(s0h1,iB,jB));
                  } // jB
                  printf("\n");
              } // iB
          } // s0h1
          printf("\n");
      } // echo

      if (nullptr != center) delete[] center;
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
