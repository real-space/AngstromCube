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
#include "boundary_condition.hxx" // Isolated_Boundary, ::periodic_image
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
  
  template<typename complex_t, typename phase_t>
  status_t kinetic_matrix(view2D<complex_t> & Tmat // result Tmat(i,j) += 0.5*<\chi3D_i|\vec\nabla \cdot \vec\nabla|\chi3D_j>
                       , view3D<double> const & t1D // input t1D(dir,i,j)   nabla^2 operator
                       , view4D<double> const & o1D // input o1D(dir,0,i,j) overlap operator
                       , int const numax_i, int const numax_j
                       , phase_t const phase=1
                       , double const prefactor=0.5) { // typical prefactor of the kinetic energy in Hartree atomic units

      auto const phase_f = phase * prefactor;

      int izyx{0};
      for    (int iz = 0; iz <= numax_i;           ++iz) {
        for  (int iy = 0; iy <= numax_i - iz;      ++iy) {
          for(int ix = 0; ix <= numax_i - iz - iy; ++ix) {

            int jzyx{0};
            for    (int jz = 0; jz <= numax_j;           ++jz) { auto const Tz = t1D(2,iz,jz), oz = o1D(2,0,iz,jz);
              for  (int jy = 0; jy <= numax_j - jz;      ++jy) { auto const Ty = t1D(1,iy,jy), oy = o1D(1,0,iy,jy);
                for(int jx = 0; jx <= numax_j - jz - jy; ++jx) { auto const Tx = t1D(0,ix,jx), ox = o1D(0,0,ix,jx);

                  Tmat(izyx,jzyx) += phase_f * ( Tx*oy*oz + ox*Ty*oz + ox*oy*Tz );

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

  
  template<typename complex_t>
  char const * complex_name(complex_t const x) {
      auto const nByte = sizeof(std::real(x));
      auto const nReIm = sizeof(complex_t);
      if (8 == nByte) {
          return (nReIm > nByte) ? "complex<double>" : "double";
      } else if (4 == nByte) {
          return (nReIm > nByte) ? "complex<float>" : "float";
      } else if (2 == nByte) {
          return (nReIm > nByte) ? "complex<half>" : "half";
      } else if (16 == nByte) {
          return (nReIm > nByte) ? "complex<quad>" : "quad";
      } else {
          return (nReIm > nByte) ? "complex<unknown>" : "unknown";
      }
  } // complex_name
  
  
  template<typename complex_t, typename phase_t>
  status_t solve_k(int const natoms
          , view2D<double const> const & xyzZ // (natoms, 4)
          , int const numaxs[]
          , double const sigmas[]
          , int const n_periodic_images
          , view2D<double const> const & periodic_image
          , view2D<int8_t const> const & periodic_shift
          , std::vector<double> const Vcoeffs[]
          , view4D<int> const & center_map
          , int const nB, int const nBa // basis size and matrix stride
          , int const offset[]
          , int const natoms_PAW
          , view2D<double const> const & xyzZ_PAW // (natoms_PAW, 4)
          , int const numax_PAW[]
          , double const sigma_PAW[]
          , view3D<double> const hs_PAW[] // [natoms_PAW](2, nprj, nprj)
          , phase_t const Bloch_phase[3]
          , int const ikp=-1 // for debugging
          , int const echo=0) { // log-level
        
      status_t stat(0);
#ifdef DEVEL
      complex_t x{1};
      if (echo > 3) printf("\n\n# start %s<%s, phase_t=%s> nB=%d\n", 
              __func__, complex_name(x), complex_name(*Bloch_phase), nB);
#endif     

      view3D<complex_t> SHm(2, nB, nBa, complex_t(0)); // get memory for 0:Overlap S and 1:Hamiltonian matrix H

      //
      // now construct the Hamiltonian:
      //   -- kinetic energy,             c.f. sho_overlap::test_simple_crystal
      //   -- local potential,            c.f. sho_potential::test_local_potential_matrix_elements (method=2)
      //   -- non-local PAW Hamiltonian contributions
      // and the overlap operator:
      //   -- SHO basis function overlap, c.f. sho_overlap::test_simple_crystal
      //   -- non-local PAW charge-deficit contributions
      //
      
#ifdef DEVEL
      // prefactor of kinetic energy in Hartree atomic units
      double const kinetic = control::get("sho_hamiltonian.scale.kinetic", 1.0) * 0.5;
      if (0.5 != kinetic) warn("kinetic energy prefactor is %g", kinetic);
#else
      double constexpr kinetic = 0.5; // prefactor of kinetic energy in Hartree atomic units
#endif
      double const ones[1] = {1.0}; // expansion of the identity (constant==1) into x^{m_x} y^{m_y} z^{m_z}

      
      std::vector<std::vector<view2D<complex_t>>> S_iaja(natoms); // construct sub-views for each atom pair
      std::vector<std::vector<view2D<complex_t>>> H_iaja(natoms); // construct sub-views for each atom pair
      // the sub-views help to address the operators as S_iaja[ia][ja](ib,jb) although their true memory
      // layout is equivalent to view4D<double> S(ia,ib,ja,jb)
      // The vector<vector<>> construction allows for sparse operators in a later stage
      for(int ia = 0; ia < natoms; ++ia) {
          S_iaja[ia].resize(natoms);
          H_iaja[ia].resize(natoms);
          int const n1i = sho_tools::n1HO(numaxs[ia]);
          for(int ja = 0; ja < natoms; ++ja) {
              S_iaja[ia][ja] = view2D<complex_t>(&(SHm(0,offset[ia],offset[ja])), nBa); // wrapper to sub-blocks of the overlap matrix
              H_iaja[ia][ja] = view2D<complex_t>(&(SHm(1,offset[ia],offset[ja])), nBa); // wrapper to sub-blocks of the Hamiltonian matrix
              int const n1j = sho_tools::n1HO(numaxs[ja]);

              for(int ip = 0; ip < n_periodic_images; ++ip) {
                  int const ic = center_map(ia,ja,ip,0); // expansion center index
                  int const numax_V = center_map(ia,ja,ip,1); // expansion of the local potential into x^{m_x} y^{m_y} z^{m_z} around a given expansion center
                  int const maxmoment = std::max(0, numax_V);

                  phase_t phase{1};
                  view3D<double> nabla2(3, n1i + 1, n1j + 1, 0.0);         //  <\chi1D_i|d/dx  d/dx|\chi1D_j>
                  view4D<double> ovl1Dm(3, 1 + maxmoment, n1i, n1j, 0.0);  //  <\chi1D_i| x^moment |\chi1D_j>
                  for(int d = 0; d < 3; ++d) { // spatial directions x,y,z
                      double const distance = xyzZ(ia,d) - (xyzZ(ja,d) + periodic_image(ip,d));
                      phase *= std::pow(Bloch_phase[d], periodic_shift(ip,d));
                      stat += sho_overlap::nabla2_matrix(nabla2[d].data(), distance, n1i + 1, n1j + 1, sigmas[ia], sigmas[ja]);
                      stat += sho_overlap::moment_tensor(ovl1Dm[d].data(), distance, n1i, n1j, sigmas[ia], sigmas[ja], maxmoment);
                  } // d

                  // construct the overlap matrix of SHO basis functions
                  // Smat(i,j) := ovl_x(ix,jx) * ovl_y(iy,jy) * ovl_z(iz,jz)
                  stat += sho_potential::potential_matrix(S_iaja[ia][ja], ovl1Dm, ones, 0, numaxs[ia], numaxs[ja], phase);

                  // add the kinetic energy contribution
                  stat += sho_hamiltonian::kinetic_matrix(H_iaja[ia][ja], nabla2, ovl1Dm, numaxs[ia], numaxs[ja], phase, kinetic);

                  // add the contribution of the local potential
                  stat += sho_potential::potential_matrix(H_iaja[ia][ja], ovl1Dm, Vcoeffs[ic].data(), numax_V, numaxs[ia], numaxs[ja], phase);

              } // ip
          } // ja
      } // ia

      std::vector<std::vector<view2D<complex_t>>> P_iaka(natoms); // potentially sparse lists, e.g. compressed row format
      std::vector<std::vector<view3D<complex_t>>> Psh_iala(natoms); // atom-centered PAW matrices muliplied to P_iaka
      for(int ia = 0; ia < natoms; ++ia) {
          P_iaka[ia].resize(natoms_PAW);
          Psh_iala[ia].resize(natoms_PAW);
          int const nb_ia = sho_tools::nSHO(numaxs[ia]);
          int const n1i   = sho_tools::n1HO(numaxs[ia]);
          for(int ka = 0; ka < natoms_PAW; ++ka) { // does not account for periodic images, ToDo
              int const nb_ka = sho_tools::nSHO(numax_PAW[ka]);
              int const n1k   = sho_tools::n1HO(numax_PAW[ka]);
              P_iaka[ia][ka] = view2D<complex_t>(nb_ia, nb_ka, 0.0); // get memory and initialize

              view4D<double> ovl1D(3, 1, n1i, n1k, 0.0);  //  <\chi1D_i|\chi1D_k>
              for(int ip = 0; ip < n_periodic_images; ++ip) {
                  phase_t phase{1};
                  for(int d = 0; d < 3; ++d) { // spatial directions x,y,z
                      double const distance = xyzZ(ia,d) - (xyzZ_PAW(ka,d) + periodic_image(ip,d));
                      phase *= std::pow(Bloch_phase[d], periodic_shift(ip,d));
                      stat += sho_overlap::overlap_matrix(ovl1D(d,0,0), distance, n1i, n1k, sigmas[ia], sigma_PAW[ka]);
                  } // d

                  // P(ia,ka,i,k) := ovl_x(ix,kx) * ovl_y(iy,ky) * ovl_z(iz,kz)
                  stat += sho_potential::potential_matrix(P_iaka[ia][ka], ovl1D, ones, 0, numaxs[ia], numax_PAW[ka], phase);
              } // ip

              // multiply P from left to hs_PAW (block diagonal --> ka == la)
              auto const la = ka;
              int const nb_la = nb_ka;
              Psh_iala[ia][la] = view3D<complex_t>(2, nb_ia, nb_la, 0.0); // get memory and initialize
              for(int ib = 0; ib < nb_ia; ++ib) {
                  for(int lb = 0; lb < nb_la; ++lb) {
                      complex_t s{0}, h{0};
                      for(int kb = 0; kb < nb_ka; ++kb) { // contract
                          s += P_iaka[ia][ka](ib,kb) * hs_PAW[ka](1,kb,lb);
                          h += P_iaka[ia][ka](ib,kb) * hs_PAW[ka](0,kb,lb);
                      } // kb
                      Psh_iala[ia][la](0,ib,lb) = s;
                      Psh_iala[ia][la](1,ib,lb) = h;
                  } // lb
              } // ib

          } // ka
      } // ia

#ifdef DEVEL
      double const scale_h = control::get("sho_hamiltonian.scale.nonlocal.h", 1.0);
      double const scale_s = control::get("sho_hamiltonian.scale.nonlocal.s", 1.0);
      if (1 != scale_h || 1 != scale_s) warn("scale PAW contributions to H and S by %g and %g, respectively", scale_h, scale_s);
#else
      double constexpr scale_h = 1, scale_s = 1;
#endif

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
                          complex_t s{0}, h{0};
                          for(int lb = 0; lb < nb_la; ++lb) { // contract
                              s += Psh_iala[ia][la](0,ib,lb) * P_iaka[ja][la](jb,lb); // ToDo: needs a conjugation if complex?
                              h += Psh_iala[ia][la](1,ib,lb) * P_iaka[ja][la](jb,lb);
                          } // lb
                          S_iaja[ia][ja](ib,jb) += s * scale_s;
                          H_iaja[ia][ja](ib,jb) += h * scale_h;
                      } // jb
                      
                  } // ib
              } // ia
          } // la
      } // ja

      // release the memory
      Psh_iala.clear();
      P_iaka.clear();

      std::vector<double> eigvals(nB, 0.0);
      auto const ovl_eig = int(control::get("sho_hamiltonian.test.overlap.eigvals", 0.));

      for(int s0h1 = 0; s0h1 < 2; ++s0h1) { // loop must run forward and serial
          if (echo > 0) printf("\n");
          auto const matrix_name = s0h1 ? "Hamiltonian" : "overlap";
          auto const  u = s0h1 ?  eV :  1; // output unit conversion factor 
          auto const _u = s0h1 ? _eV : ""; // unit symbol

          // display S and H
          if (echo > 9 - s0h1) {
              printf("\n# %s matrix (%s) for Bloch phase", matrix_name, _u);
              for(int d = 0; d < 3; ++d) {
                  printf(" %g %g ", std::real(Bloch_phase[d]), std::imag(Bloch_phase[d]));
              } // d
              printf(":\n");
              for(int iB = 0; iB < nB; ++iB) {
                  printf("# row%3i ", iB);
                  for(int jB = 0; jB < nB; ++jB) {
                      printf("%8.3f", SHm(s0h1,iB,jB)*u);
                  } // jB
                  printf("\n");
              } // iB
              printf("\n");
          } // echo

          // diagonalize matrix
          status_t stat_eig(0);
          if (1 == s0h1) {
              stat_eig = linear_algebra::generalized_eigenvalues(nB, SHm(1,0), nBa, SHm(0,0), nBa, eigvals.data());
          } else if (ovl_eig) {
              view2D<complex_t> S_copy(nB, nBa); // get memory
              set(S_copy.data(), nB*nBa, SHm(0,0)); // copy overlap matrix S into work array W
              stat_eig = linear_algebra::eigenvalues(nB, S_copy.data(), nBa, eigvals.data());
          } // ovl_eig

          // show result
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
                      if (nB > mB) printf(" ..."); // there are more eigenvalues than we display
                      printf(" %g %g %s\n", eigvals[nB - 2]*u, eigvals[nB - 1]*u, _u); // last two
                  } // echo
                  if (echo > 4) printf("# lowest and highest eigenvalue of the %s matrix are %g and %g %s, respectively\n", 
                                            matrix_name, lowest_eigenvalue*u, highest_eigenvalue*u, _u);
                  if (s0h1 == 0 && lowest_eigenvalue < .1) warn("overlap matrix has small eigenvalues, lowest= %g", lowest_eigenvalue);
              } // stat_eig
          } // H or ovl_eig

          if (1 == s0h1) {
              if (0 == stat_eig) {
                  // ToDo: export the entire spectrum and the lowest n_occupied_bands eigenvectors
                  //       for this k-point to be later used for the density generation.
                  //       n_occupied_bands should be chosen large so that also the tail 
                  //       of the Fermi-Dirac distribution is captured. We can check that
                  //       by ensuring that the occupation numbers of the next higher state
                  //          f_FD(eigvals[n_occupied_bands] - E_Fermi) < 10^{-15}
                  //       are small.
              } // success
          } // H

      } // s0h1
      
      return stat;
  } // solve_k
  
  
  status_t solve(int const natoms // number of SHO basis centers
          , view2D<double const> const & xyzZ // (natoms, 4)
          , real_space::grid_t const & g
          , double const *const vtot
          , int const natoms_prj // =-1 number of PAW atoms
          , double const *const sigma_prj // =nullptr
          , int    const *const numax_prj // =nullptr
          , double *const *const atom_mat // =nullptr
          , int const echo) {
    
      status_t stat(0);
      
      auto const usual_numax = int(control::get("sho_hamiltonian.test.numax", 1.));
      auto const usual_sigma =     control::get("sho_hamiltonian.test.sigma", 2.);
      std::vector<int>    numaxs(natoms, usual_numax); // define SHO basis sizes
      std::vector<double> sigmas(natoms, usual_sigma); // define SHO basis spreads
      double const sigma_asymmetry = control::get("sho_hamiltonian.test.sigma.asymmetry", 1.0);
      if (sigma_asymmetry != 1) { sigmas[0] *= sigma_asymmetry; sigmas[natoms - 1] /= sigma_asymmetry; } // manipulate the spreads
      
      
      std::vector<int> offset(natoms + 1, 0); // the basis atoms localized on atom#i start at offset[i]
      int total_basis_size{0}, maximum_numax{-1};
      double maximum_sigma{0};
      for(int ia = 0; ia < natoms; ++ia) {
          if (echo > 0) {
              printf("# atom#%i \tZ=%g \tposition %12.6f%12.6f%12.6f  numax= %d sigma= %.3f %s\n",
                  ia, xyzZ(ia,3), xyzZ(ia,0)*Ang, xyzZ(ia,1)*Ang, xyzZ(ia,2)*Ang, numaxs[ia], sigmas[ia]*Ang, _Ang);
          } // echo
          int const atom_basis_size = sho_tools::nSHO(numaxs[ia]);
          maximum_numax = std::max(maximum_numax, numaxs[ia]);
          maximum_sigma = std::max(maximum_sigma, sigmas[ia]);
          offset[ia] = total_basis_size; // prefix sum
          total_basis_size += atom_basis_size;
      } // ia
      int const numax_max = maximum_numax;
      if (echo > 5) printf("# largest basis size per atom is %d, numax=%d\n", sho_tools::nSHO(numax_max), numax_max);
      offset[natoms] = total_basis_size;
      int const nB   = total_basis_size;
      int const nBa  = align<4>(nB); // memory aligned main matrix stride



      float const rcut = 9*maximum_sigma; // exp(-9^2) = 6.6e-36
      double *periodic_image_ptr{nullptr};
      int8_t *periodic_shift_ptr{nullptr};
      double const cell[] = {g[0]*g.h[0], g[1]*g.h[1], g[2]*g.h[2]};
      int const n_periodic_images = boundary_condition::periodic_images(&periodic_image_ptr, cell, g.boundary_conditions(), rcut, 0, &periodic_shift_ptr);
      if (echo > 1) printf("# %s consider %d periodic images\n", __FILE__, n_periodic_images);
      view2D<double const> periodic_image(periodic_image_ptr, 4); // wrap
      view2D<int8_t const> periodic_shift(periodic_shift_ptr, 4); // wrap

      double const origin[] = {.5*(g[0] - 1)*g.h[0],
                               .5*(g[1] - 1)*g.h[1], 
                               .5*(g[2] - 1)*g.h[2]};
      
      int ncenters = natoms*natoms*n_periodic_images;
      view2D<double> center(ncenters, 8, 0.0); // list of potential expansion centers
      view4D<int> center_map(natoms, natoms, n_periodic_images, 2, -1);
      { // scope: set up list of centers for the expansion of the local potential
          int ic{0};
          for(int ia = 0; ia < natoms; ++ia) {
              for(int ja = 0; ja < natoms; ++ja) {

                  double const alpha_i = 1/pow2(sigmas[ia]);
                  double const alpha_j = 1/pow2(sigmas[ja]);
                  double const sigma_V = 1/std::sqrt(alpha_i + alpha_j);
                  double const wi = alpha_i*pow2(sigma_V);
                  double const wj = alpha_j*pow2(sigma_V);
                  assert( std::abs( wi + wj - 1.0 ) < 1e-12 );
                  for(int ip = 0; ip < n_periodic_images; ++ip) {
                      int const numax_V = numaxs[ia] + numaxs[ja]; // depending on the distance between atom#ia and the periodic image of atom#ja, this could be lowered
                      if (echo > 7) printf("# ai#%i aj#%i \tcenter of weight\t", ia, ja);
                      for(int d = 0; d < 3; ++d) { // spatial directions x,y,z
                          center(ic,d) = wi*xyzZ(ia,d) + wj*(xyzZ(ja,d) + periodic_image(ip,d));
                          if (echo > 7) printf("%12.6f", center(ic,d)*Ang);
                          center(ic,d) += origin[d];
                      } // d
                      if (echo > 7) printf("  numax_V= %i sigma_V: %g %s \n", numax_V, sigma_V*Ang, _Ang);
                      center(ic,3) = sigma_V;
                      center(ic,4) = numax_V;
                      // debug info
                      center(ic,5) = ip; // index of the periodic image
                      center(ic,6) = ia; // ToDo: use global atom indices
                      center(ic,7) = ja; // ToDo: use global atom indices

                      center_map(ia,ja,ip,0) = ic;
                      center_map(ia,ja,ip,1) = numax_V;

                      ++ic;
                  } // ip
              } // ja
          } // ia
          ncenters = ic; // may be less that initially allocated
      } // scope: set up list of centers
      if (echo > 7) printf("# project local potential at %d sites\n", ncenters);

      // 
      // Potential expansion centers that are on the same location and have the same sigma_V
      // could be merged to reduce the projection efforts. Even if the numax_V do not match, 
      // we take the higher one and take advantage of the order_Ezyx of the coefficients 
      // after normalize_potential_coefficients.
      // 
      
      
      // perform the projection of the local potential
      std::vector<std::vector<double>> Vcoeffs(ncenters);
      double const scale_potential = control::get("sho_hamiltonian.scale.potential", 1.0);
      if (1.0 != scale_potential) warn("scale potential by %g", scale_potential);
      for(int ic = 0; ic < ncenters; ++ic) {
          double const sigma_V = center(ic,3);
          int    const numax_V = center(ic,4);
          Vcoeffs[ic] = std::vector<double>(sho_tools::nSHO(numax_V), 0.0);
          stat += sho_projection::sho_project(Vcoeffs[ic].data(), numax_V, center[ic], sigma_V, vtot, g, 0); // 0:mute
          // now Vcoeff is represented w.r.t. to Hermite polynomials H_{nx}*H_{ny}*H_{nz} and order_zyx
          
          stat += sho_potential::normalize_potential_coefficients(Vcoeffs[ic].data(), numax_V, sigma_V, 0); // 0:mute
          // now Vcoeff is represented w.r.t. powers of the Cartesian coords x^{nx}*y^{ny}*z^{nz} and order_Ezyx
          
          scale(Vcoeffs[ic].data(), Vcoeffs[ic].size(), scale_potential);
      } // ic


      // prepare for the PAW contributions: find the projection coefficient matrix P = <\chi3D_{ia ib}|\tilde p_{ka kb}>
      int const natoms_PAW = (natoms_prj < 0) ? natoms : std::min(natoms, natoms_prj); // keep it flexible
      auto const xyzZ_PAW = view2D<double const>(xyzZ.data(), xyzZ.stride()); // duplicate view
      std::vector<int32_t> numax_PAW(natoms_PAW,  3);
      std::vector<double>  sigma_PAW(natoms_PAW, .5);
      std::vector<view3D<double>> hs_PAW(natoms_PAW); // atomic matrices for charge-deficit and Hamiltonian corrections
      for(int ka = 0; ka < natoms_PAW; ++ka) {
          if (numax_prj) numax_PAW[ka] = numax_prj[ka];
          if (sigma_prj) sigma_PAW[ka] = sigma_prj[ka];
          int const nb_ka = sho_tools::nSHO(numax_PAW[ka]); // number of projectors
          if (atom_mat) {
              hs_PAW[ka] = view3D<double>(atom_mat[ka], nb_ka, nb_ka); // wrap input
          } else {
              hs_PAW[ka] = view3D<double>(2, nb_ka, nb_ka, 0.0); // get memory and initialize zero
          } // atom_mat
          // for both matrices: L2-normalization w.r.t. SHO projector functions is important
      } // ka

      // all preparations done, start k-point loop
      
      auto const nkpoints = int(control::get("sho_hamiltonian.test.kpoints", 8.));
      for(int ikp = 0; ikp < nkpoints; ++ikp) {
//        std::complex<double> Bloch_phase[3] = {1 - 2.*(ikp & 1), 1. - (ikp & 2), 1. - .5*(ikp & 4)}; // one of the 8 real k-points, Gamma and X-points
          std::complex<double> constexpr minus_one = -1;
          std::complex<double> Bloch_phase[3] = {std::pow(minus_one, ikp/(nkpoints - 1.)), 1, 1}; // first and last phases are real, dispersion in x-direction

          double Bloch_phase_real[3];
          bool needs_complex{false};
          for(int d = 0; d < 3; ++d) {
              Bloch_phase_real[d] = Bloch_phase[d].real();
              if (std::abs(Bloch_phase[d].imag()) > 2e-16) needs_complex = true;
          } // d
          if (needs_complex) {
              stat += solve_k<std::complex<double>>(
                          natoms, xyzZ, numaxs.data(), sigmas.data(),
                          n_periodic_images, periodic_image, periodic_shift,
                          Vcoeffs.data(), center_map,
                          nB, nBa, offset.data(),
                          natoms_PAW, xyzZ_PAW, numax_PAW.data(), sigma_PAW.data(), hs_PAW.data(),
                          Bloch_phase, ikp, echo);
          } else {
              stat += solve_k<double>(
                          natoms, xyzZ, numaxs.data(), sigmas.data(),
                          n_periodic_images, periodic_image, periodic_shift,
                          Vcoeffs.data(), center_map,
                          nB, nBa, offset.data(),
                          natoms_PAW, xyzZ_PAW, numax_PAW.data(), sigma_PAW.data(), hs_PAW.data(),
                          Bloch_phase_real, ikp, echo);
          } // needs_complex
          // ToDo: more cases for float and std::complex<float>, needs extension of linear_algebra module
          
      } // ikp
      
      return stat;
  } // solve

  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_Hamiltonian(int const echo=5) {
      status_t stat(0);
      
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
      view2D<double const> const xyzZ(xyzZ_m, 4); // wrap for simpler usage

      real_space::grid_t g(dims);
      g.set_boundary_conditions(bc);
      g.set_grid_spacing(cell[0]/g[0], cell[1]/g[1], cell[2]/g[2]);
      if (echo > 1) printf("# use  %g %g %g %s grid spacing\n", g.h[0]*Ang, g.h[1]*Ang, g.h[2]*Ang, _Ang);
      if (echo > 1) printf("# cell is  %g %g %g %s\n", g.h[0]*g[0]*Ang, g.h[1]*g[1]*Ang, g.h[2]*g[2]*Ang, _Ang);
      
      stat += solve(natoms, xyzZ, g, vtot.data(), 0, 0, 0, 0, echo);

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
