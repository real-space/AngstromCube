#include <cstdio> // std::printf
#include <cmath> // std::sqrt, ::abs, ::pow
#include <algorithm> // std::max, ::min
#include <complex> // std::complex<real_t>, ::imag, ::real
#include <vector> // std::vector<T>
#include <cassert> // assert
#include <set> // std::set<key>
#include <numeric> // std::iota

#include "sho_hamiltonian.hxx"

#include "sho_potential.hxx" // ::load_local_potential, ::normalize_potential_coefficients, ::potential_matrix
#include "geometry_analysis.hxx" // ::read_xyz_file, ::fold_back
#include "control.hxx" // ::get
#include "display_units.h" // eV, _eV, Ang, _Ang, Kelvin, _Kelvin
#include "real_space.hxx" // ::grid_t
#include "sho_tools.hxx" // ::nSHO, ::n1HO, ::order_*, ::construct_label_table
#include "sho_projection.hxx" // ::sho_project
#include "boundary_condition.hxx" // Isolated_Boundary, ::periodic_image
#include "sho_overlap.hxx" // ::overlap_matrix, ::nabla2_matrix, ::moment_tensor
#include "data_view.hxx" // view2D<T>, view3D<T>, view4D<T>
#include "simple_math.hxx" // ::random<real_t>, align<nBits>
#include "simple_stats.hxx" // ::Stats
#include "simple_timer.hxx" // SimpleTimer
#include "vector_math.hxx" // ::vec<N,T>
#include "dense_solver.hxx" // ::solve

namespace sho_hamiltonian {
  // computes Hamiltonian matrix elements between two SHO basis functions
  // including a non-local PAW contribution

  template <typename complex_t, typename phase_t>
  status_t kinetic_matrix(
        view2D<complex_t> & Tmat // result Tmat(i,j) += prefactor*<\chi3D_i|\vec\nabla \cdot \vec\nabla|\chi3D_j>
      , view4D<double> const & o1D // input o1D(dir,0,i,j) overlap operator in 1D
      , view3D<double> const & t1D // input t1D(dir,i,j)   nabla^2 operator in 1D
      , int const numax_i
      , int const numax_j
      , phase_t const phase=1 // typically phase_t is double or complex<double>
      , double const prefactor=0.5 // typical prefactor of the kinetic energy in Hartree atomic units
  ) {
      auto const phase_f = phase * prefactor;

      int izyx{0};
      for     (int iz = 0; iz <= numax_i;           ++iz) {
        for   (int iy = 0; iy <= numax_i - iz;      ++iy) {
          for (int ix = 0; ix <= numax_i - iz - iy; ++ix) {

            int jzyx{0};
            for     (int jz = 0; jz <= numax_j;           ++jz) { auto const Tz = t1D(2,iz,jz), oz = o1D(2,0,iz,jz);
              for   (int jy = 0; jy <= numax_j - jz;      ++jy) { auto const Ty = t1D(1,iy,jy), oy = o1D(1,0,iy,jy);
                for (int jx = 0; jx <= numax_j - jz - jy; ++jx) { auto const Tx = t1D(0,ix,jx), ox = o1D(0,0,ix,jx);

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

  template <typename complex_t, typename phase_t>
  status_t solve_k(
        int const natoms // number of SHO basis centers
      , view2D<double> const & xyzZ // (natoms, 4) positions of SHO basis centers, 4th component not used
      , int    const numaxs[] // spreads of the SHO basis
      , double const sigmas[] // cutoffs of the SHO basis
      , int const n_periodic_images // number of periodic images
      , view2D<double> const & periodic_image // periodic image coordinates
      , view2D<int8_t> const & periodic_shift // periodic image shifts
      , std::vector<double> const Vcoeffs[] // Vcoeff[ic][0..Ezyx..nSHO(numax_V)]
      , view4D<int> const & center_map // ic = center_map(ia,ja,ip,0),  numax_V = center_map(ia,ja,ip,1)
      , int const nB, int const nBa // basis size and matrix stride
      , int const offset[] // beginning of matrix blocks
      , int const natoms_PAW // number of PAW centers
      , view2D<double const> const & xyzZ_PAW // (natoms_PAW, 4) positions of PAW centers, 4th component not used
      , int    const numax_PAW[] // spreads of the SHO-type PAW projectors
      , double const sigma_PAW[] // cutoffs of the SHO-type PAW projectors
      , view3D<double> const hs_PAW[] // [natoms_PAW](2, nprj, nprj) // PAW Hamiltonian correction and charge-deficit
      , phase_t const Bloch_phase[3]
      , char const *const x_axis // display this string in front of the Hamiltonian eigenvalues
      , int const echo=0 // log-level
  ) {
      using real_t = decltype(std::real(complex_t(1))); // base type

      status_t stat(0);
#ifdef DEVEL
      if (echo > 3) std::printf("\n\n# start %s<%s, phase_t=%s> nB=%d n_periodic_images=%d\n", 
          __func__, complex_name<complex_t>(), complex_name(*Bloch_phase), nB, n_periodic_images);
#endif // DEVEL

      //
      // now construct the Hamiltonian:
      //   -- kinetic energy,             c.f. sho_overlap::test_simple_crystal
      //   -- local potential,            c.f. sho_potential::test_local_potential_matrix_elements (method=between)
      //   -- non-local PAW Hamiltonian contributions
      // and the overlap operator:
      //   -- SHO basis function overlap, c.f. sho_overlap::test_simple_crystal
      //   -- non-local PAW charge-deficit contributions
      //

      double const ones[1] = {1.0}; // expansion of the identity (constant==1) into x^{m_x} y^{m_y} z^{m_z}
      int constexpr S=1, H=0; // static indices for S:overlap matrix, H:Hamiltonian matrix

      // PAW contributions to H_{ij} = P_{ik} h_{kl} P^*_{jl} = Ph_{il} P^*_{jl}
      //                  and S_{ij} = P_{ik} s_{kl} P^*_{jl} = Ps_{il} P^*_{jl}

      // we use vectors of vectors here for potentially sparse list in the future, e.g. in the compressed row format
      std::vector<std::vector<view2D<complex_t>>> P_jala(natoms); //  <\tilde p_{ka kb}|\chi3D_{ja jb}>
      std::vector<std::vector<view3D<complex_t>>> Psh_iala(natoms); // atom-centered PAW matrices multiplied to P_jala
      for (int ia = 0; ia < natoms; ++ia) {
          P_jala[ia].resize(natoms_PAW);
          Psh_iala[ia].resize(natoms_PAW);
          int const nb_ia = sho_tools::nSHO(numaxs[ia]);
          int const n1i   = sho_tools::n1HO(numaxs[ia]);
          for (int ka = 0; ka < natoms_PAW; ++ka) {
              int const nb_ka = sho_tools::nSHO(numax_PAW[ka]);
              int const n1k   = sho_tools::n1HO(numax_PAW[ka]);
              P_jala[ia][ka] = view2D<complex_t>(nb_ia, nb_ka, 0.0); // get memory and initialize

              view4D<double> ovl1D(3, 1, n1i, n1k, 0.0);  //  <\chi1D_i|\chi1D_k>
              for (int ip = 0; ip < n_periodic_images; ++ip) { // periodic images of ka
                  phase_t phase{1};
                  for (int d = 0; d < 3; ++d) { // spatial directions x,y,z
                      auto const phase_factor = std::pow(Bloch_phase[d], int(periodic_shift(ip,d)));
                      double const distance = xyzZ(ia,d) - (xyzZ_PAW(ka,d) + periodic_image(ip,d));
                      stat += sho_overlap::overlap_matrix(ovl1D(d,0,0), distance, n1i, n1k, sigmas[ia], sigma_PAW[ka]);
                      phase *= phase_factor;
                  } // d

                  // P(ia,ka,i,k) := ovl_x(ix,kx) * ovl_y(iy,ky) * ovl_z(iz,kz)
                  stat += sho_potential::potential_matrix(P_jala[ia][ka], ovl1D, ones, 0, numaxs[ia], numax_PAW[ka], phase);
              } // ip

              // multiply P from left to hs_PAW (block diagonal --> ka == la)
              auto const la = ka;
              int const nb_la = nb_ka;
              Psh_iala[ia][la] = view3D<complex_t>(2, nb_ia, nb_la, 0.0); // get memory and initialize
              for (int ib = 0; ib < nb_ia; ++ib) {
//                   if (echo > -1) {
//                       std::printf("# ia=%i ka=%i ib=%i kb=0..%i ", ia, ka, ib, nb_ka);
//                       for (int kb = 0; kb < nb_ka; ++kb) {
//                           auto const p = P_jala[ia][ka](ib,kb);
//                           std::printf("%10.6f%9.6f", std::real(p), std::imag(p));
//                       } // kb
//                       std::printf("\n");
//                   } // echo
                  for (int lb = 0; lb < nb_la; ++lb) {
                      complex_t s(0), h(0);
                      for (int kb = 0; kb < nb_ka; ++kb) { // contract
                          auto const p = P_jala[ia][ka](ib,kb);
                          // Mind that hs_PAW matrices are ordered {0:H,1:S}
                          s += p * real_t(hs_PAW[ka](1,kb,lb));
                          h += p * real_t(hs_PAW[ka](0,kb,lb));
                      } // kb
                      Psh_iala[ia][la](S,ib,lb) = s;
                      Psh_iala[ia][la](H,ib,lb) = h;
                  } // lb
              } // ib

          } // ka
      } // ia
      // PAW projection matrix ready

      // allocate bulk memory for overlap and Hamiltonian
      view3D<complex_t> HSm(2, nB, nBa, complex_t(0)); // get memory for 1:Overlap S and 0:Hamiltonian matrix H
      
#ifdef DEVEL
      double const scale_k = control::get("hamiltonian.scale.kinetic", 1.0);
      if (1 != scale_k) warn("kinetic energy is scaled by %g", scale_k);
      real_t const scale_h = control::get("hamiltonian.scale.nonlocal.h", 1.0);
      real_t const scale_s = control::get("hamiltonian.scale.nonlocal.s", 1.0);
      if (1 != scale_h || 1 != scale_s) warn("scale PAW contributions to H and S by %g and %g, respectively", scale_h, scale_s);
#else
      real_t constexpr scale_h = 1, scale_s = 1;
      double constexpr scale_k = 1;
#endif // DEVEL
      double const kinetic = 0.5 * scale_k; // prefactor of kinetic energy in Hartree atomic units

      // we construct sub-views to the blocks for each atom pair
      std::vector<std::vector<view2D<complex_t>>> S_iaja(natoms); // <\chi3D_i| 1 + s_PAW |\chi3D_j>
      std::vector<std::vector<view2D<complex_t>>> H_iaja(natoms); // <\chi3D_i| \hat T + \tilde V(\vec r) + h_PAW |\chi3D_j>
      // the sub-views help to address the operators as S_iaja[ia][ja](ib,jb) although their true memory
      // layout is equivalent to view4D<double> S(ia,ib,ja,jb)
      // The vector<vector<>> construction allows for sparse operators in a later stage
      for (int ia = 0; ia < natoms; ++ia) {
          S_iaja[ia].resize(natoms);
          H_iaja[ia].resize(natoms);
          int const n1i   = sho_tools::n1HO(numaxs[ia]);
          int const nb_ia = sho_tools::nSHO(numaxs[ia]);
          for (int ja = 0; ja < natoms; ++ja) {
              S_iaja[ia][ja] = view2D<complex_t>(&(HSm(S,offset[ia],offset[ja])), HSm.stride()); // wrapper to sub-blocks of the overlap matrix
              H_iaja[ia][ja] = view2D<complex_t>(&(HSm(H,offset[ia],offset[ja])), HSm.stride()); // wrapper to sub-blocks of the Hamiltonian matrix
              int const n1j   = sho_tools::n1HO(numaxs[ja]);
              int const nb_ja = sho_tools::nSHO(numaxs[ja]);

              for (int ip = 0; ip < n_periodic_images; ++ip) { // periodic images of ja
                  int const ic      = center_map(ia,ja,ip,0); // expansion center index
                  int const numax_V = center_map(ia,ja,ip,1); // expansion of the local potential into x^{m_x} y^{m_y} z^{m_z} around a given expansion center
                  int const maxmoment = std::max(0, numax_V);

                  phase_t phase(1);
                  view3D<double> nabla2(3, n1i + 1, n1j + 1, 0.0);         //  <\chi1D_i|d/dx  d/dx|\chi1D_j>
                  view4D<double> ovl1Dm(3, 1 + maxmoment, n1i, n1j, 0.0);  //  <\chi1D_i| x^moment |\chi1D_j>
                  for (int d = 0; d < 3; ++d) { // spatial directions x,y,z
                      auto const phase_factor = std::pow(Bloch_phase[d], int(periodic_shift(ip,d)));
                      double const distance = xyzZ(ia,d) - (xyzZ(ja,d) + periodic_image(ip,d));
//                    std::printf("# (%g,%g)^%d = (%g,%g)\n", std::real(Bloch_phase[d]), std::imag(Bloch_phase[d]), int(periodic_shift(ip,d)), std::real(phase_factor), std::imag(phase_factor));
                      stat += sho_overlap::nabla2_matrix(nabla2[d].data(), distance, n1i + 1, n1j + 1, sigmas[ia], sigmas[ja]);
                      stat += sho_overlap::moment_tensor(ovl1Dm[d].data(), distance, n1i    , n1j    , sigmas[ia], sigmas[ja], maxmoment);
                      phase *= phase_factor;
                  } // d

                  // add the kinetic energy contribution
                  stat += sho_hamiltonian::kinetic_matrix(H_iaja[ia][ja], ovl1Dm, nabla2,                      numaxs[ia], numaxs[ja], phase, kinetic);

                  // add the contribution of the local potential
                  stat += sho_potential::potential_matrix(H_iaja[ia][ja], ovl1Dm, Vcoeffs[ic].data(), numax_V, numaxs[ia], numaxs[ja], phase);

                  // construct the overlap matrix of SHO basis functions
                  // Smat(i,j) := ovl_x(ix,jx) * ovl_y(iy,jy) * ovl_z(iz,jz)
//                std::printf("# phase = %g %g, shifts %d %d %d\n", std::real(phase), std::imag(phase), periodic_shift(ip,0), periodic_shift(ip,1), periodic_shift(ip,2));
                  stat += sho_potential::potential_matrix(S_iaja[ia][ja], ovl1Dm,               ones,       0, numaxs[ia], numaxs[ja], phase);

              } // ip


              // PAW contributions to H_{ij} = Ph_{il} P^*_{jl}
              //                  and S_{ij} = Ps_{il} P^*_{jl}
              for (int la = 0; la < natoms_PAW; ++la) { // contract
                  int const nb_la = sho_tools::nSHO(numax_PAW[la]);
                  // triple-loop matrix-matrix multiplication
                  for (int ib = 0; ib < nb_ia; ++ib) {
                      for (int jb = 0; jb < nb_ja; ++jb) {
                          complex_t s(0), h(0);
                          for (int lb = 0; lb < nb_la; ++lb) { // contract
                              auto const p = conjugate(P_jala[ja][la](jb,lb)); // needs a conjugation if complex
                              s += Psh_iala[ia][la](0,ib,lb) * p;
                              h += Psh_iala[ia][la](1,ib,lb) * p;
                          } // lb
                          S_iaja[ia][ja](ib,jb) += s * scale_s;
                          H_iaja[ia][ja](ib,jb) += h * scale_h;
                      } // jb
                  } // ib
              } // la

          } // ja
      } // ia

      Psh_iala.clear(); // release the memory, P_jala is still needed for the generation of density matrices
      S_iaja.clear(); H_iaja.clear(); // release the sub-views, matrix elements are still stored in HSm

      return dense_solver::solve(HSm, x_axis, echo);
  } // solve_k







  status_t solve(
        int const natoms // number of SHO basis centers
      , view2D<double> const & xyzZ // (natoms, 4)
      , real_space::grid_t const & g // Cartesian grid descriptor for vtot
      , double const *const vtot // total effective potential on grid
      , int const nkpoints
      , view2D<double> const & kmesh // kmesh(nkpoints, 4);
      , int const natoms_prj // =-1 number of PAW atoms
      , double const *const sigma_prj // =nullptr
      , int    const *const numax_prj // =nullptr
      , double *const *const atom_mat // =nullptr
      , int const echo // =0 log-level
  ) {
      status_t stat(0);
      SimpleTimer prepare_timer(__FILE__, __LINE__, "prepare", 0);

      int  const usual_numax = control::get("sho_hamiltonian.test.numax", 1.);
      auto const usual_sigma = control::get("sho_hamiltonian.test.sigma", .5);
      std::vector<int>    numaxs(natoms, usual_numax); // define SHO basis sizes
      std::vector<double> sigmas(natoms, usual_sigma); // define SHO basis spreads
      double const sigma_asymmetry = control::get("sho_hamiltonian.test.sigma.asymmetry", 1.0);
      if (sigma_asymmetry != 1) { sigmas[0] *= sigma_asymmetry; sigmas[natoms - 1] /= sigma_asymmetry; } // manipulate the spreads


      std::vector<int> offset(natoms + 1, 0); // the basis atoms localized on atom#i start at offset[i]
      int total_basis_size{0}, maximum_numax{-1};
      double maximum_sigma{0};
      for (int ia = 0; ia < natoms; ++ia) {
          if (echo > 0) {
              std::printf("# atom#%i \tZ=%g \tposition %12.6f%12.6f%12.6f  numax= %d sigma= %.3f %s\n",
                  ia, xyzZ(ia,3), xyzZ(ia,0)*Ang, xyzZ(ia,1)*Ang, xyzZ(ia,2)*Ang, numaxs[ia], sigmas[ia]*Ang, _Ang);
          } // echo
          int const atom_basis_size = sho_tools::nSHO(numaxs[ia]);
          maximum_numax = std::max(maximum_numax, numaxs[ia]);
          maximum_sigma = std::max(maximum_sigma, sigmas[ia]);
          offset[ia] = total_basis_size; // prefix sum
          total_basis_size += atom_basis_size;
      } // ia
      int const numax_max = maximum_numax;
      if (echo > 5) std::printf("# largest basis size per atom is %d, numax=%d\n", sho_tools::nSHO(numax_max), numax_max);
      offset[natoms] = total_basis_size;
      int const nB   = total_basis_size;
      int const nBa  = align<3>(nB); // memory aligned main matrix stride



      // prepare for the PAW contributions: find the projection coefficient matrix P = <\chi3D_{ia ib}|\tilde p_{ka kb}>
      int const natoms_PAW = (natoms_prj < 0) ? natoms : std::min(natoms, natoms_prj); // keep it flexible
      auto const xyzZ_PAW = view2D<double const>(xyzZ.data(), xyzZ.stride()); // duplicate view
      std::vector<int32_t> numax_PAW(natoms_PAW,  3);
      std::vector<double>  sigma_PAW(natoms_PAW, .5);
      double maximum_sigma_PAW{0};
      std::vector<view3D<double>> hs_PAW(natoms_PAW); // atomic matrices for charge-deficit and Hamiltonian corrections
      for (int ka = 0; ka < natoms_PAW; ++ka) {
          if (numax_prj) numax_PAW[ka] = numax_prj[ka];
          if (sigma_prj) sigma_PAW[ka] = sigma_prj[ka];
          maximum_sigma_PAW = std::max(maximum_sigma_PAW, sigma_PAW[ka]);
          int const nb_ka = sho_tools::nSHO(numax_PAW[ka]); // number of projectors
          if (atom_mat) {
              hs_PAW[ka] = view3D<double>(atom_mat[ka], nb_ka, nb_ka); // wrap input
          } else {
              hs_PAW[ka] = view3D<double>(2, nb_ka, nb_ka, 0.0); // get memory and initialize zero
          } // atom_mat
          // for both matrices: L2-normalization w.r.t. SHO projector functions is important
      } // ka


      view2D<double> periodic_image;
      view2D<int8_t> periodic_shift;
      float const rcut = 9*std::max(maximum_sigma, maximum_sigma_PAW); // exp(-9^2) = 6.6e-36
      double const cell[] = {g[0]*g.h[0], g[1]*g.h[1], g[2]*g.h[2]};
      int const n_periodic_images = boundary_condition::periodic_images(periodic_image, cell, g.boundary_conditions(), rcut, echo, &periodic_shift);
      if (echo > 1) std::printf("# %s consider %d periodic images\n", __FILE__, n_periodic_images);
      if (echo > 9) {
          for (int ip = 0; ip < n_periodic_images; ++ip) {
              std::printf("# periodic image #%d\t shifts %3d %3d %3d\n", ip, periodic_shift(ip,0), periodic_shift(ip,1), periodic_shift(ip,2));
          } // ip
      } // echo

      double const origin[] = {.5*(g[0] - 1)*g.h[0],
                               .5*(g[1] - 1)*g.h[1], 
                               .5*(g[2] - 1)*g.h[2]};
      
      int ncenters = natoms*natoms*n_periodic_images; // non-const
      view2D<double> center(ncenters, 8, 0.0); // list of potential expansion centers
      view4D<int> center_map(natoms, natoms, n_periodic_images, 2, -1);
      double sigma_V_max{0}, sigma_V_min{9e9};
      { // scope: set up list of centers for the expansion of the local potential
          int ic{0};
          for (int ia = 0; ia < natoms; ++ia) {
              for (int ja = 0; ja < natoms; ++ja) {

                  double const alpha_i = 1./pow2(sigmas[ia]);
                  double const alpha_j = 1./pow2(sigmas[ja]);
                  double const sigma_V = 1./std::sqrt(alpha_i + alpha_j);
                  double const wi = alpha_i*pow2(sigma_V);
                  double const wj = alpha_j*pow2(sigma_V);
                  assert( std::abs( wi + wj - 1.0 ) < 1e-12 );
                  sigma_V_max = std::max(sigma_V, sigma_V_max);
                  sigma_V_min = std::min(sigma_V, sigma_V_min);
                  // account for the periodic images around atom #ja
                  for (int ip = 0; ip < n_periodic_images; ++ip) {
                      int const numax_V = numaxs[ia] + numaxs[ja]; 
                      // depending on the distance between atom#ia and the periodic image of atom#ja, numax_V could be lowered
                      if (echo > 7) std::printf("# ai#%i aj#%i \tcenter of weight\t", ia, ja);
                      for (int d = 0; d < 3; ++d) { // spatial directions x,y,z
                          center(ic,d) = wi*xyzZ(ia,d) + wj*(xyzZ(ja,d) + periodic_image(ip,d));
                          // cast coordinates into [-cell/2, cell/2)
                          center(ic,d) = geometry_analysis::fold_back(center(ic,d), cell[d]); 
                          if (echo > 7) std::printf("%12.6f", center(ic,d)*Ang);
                          center(ic,d) += origin[d];
                      } // d
                      if (echo > 7) std::printf("  numax_V= %i sigma_V= %g %s \n", numax_V, sigma_V*Ang, _Ang);
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
          ncenters = ic; // may be less than initially allocated
      } // scope: set up list of centers
      if (echo > 7) std::printf("# project local potential at %d sites\n", ncenters);

      //
      // Potential expansion centers that are on the same location and have the same sigma_V
      // could be merged to reduce the projection efforts. Even if the numax_Vs do not match, 
      // we take the higher one and take advantage of the order_Ezyx of the coefficients 
      // after normalize_potential_coefficients.
      // 
      // For this analysis, define a spatial threshold how close expansion centers
      // should be to be considered mergable, furthermore, a threshold for sigma_Vs.
      // Then, we can divide the cell into sub-cells, integerizing the coordinates
      // and compare in an N^2-algorithm each sub-cell to the neighborhood (3^4=81 subcells)

      std::vector<bool> center_active(ncenters, true); // init all centers as active
      { // scope: count the number of similar, potentially identical, expansion centers
          double const threshold_space = 1e-15; // in Bohr
          double const threshold_sigma = 1e-12; // in Bohr
          auto const method = control::get("sho_hamiltonian.reduce.centers", "none"); // {none, quadratic, smart}
          if ('s' == (*method | 32)) { // "smart"
              if (echo > 5) std::printf("# sho_hamiltonian.reduce.centers='%s' --> 'smart'\n", method);

              double ab[4][2]; // offset and slope
              for (int d = 0; d < 3; ++d) {
                  // map coordinates from [-cell/2, cell/2) to [0, 65535]
                  ab[d][0] = .5*cell[d] - origin[d];
                  ab[d][1] = std::min(1./threshold_space, 65535/cell[d]); // slope
              } // d
              ab[3][0] = -(sigma_V_min - threshold_sigma);
              ab[3][1] = std::min(1./threshold_sigma, 65535/(sigma_V_max - sigma_V_min + 2*threshold_sigma));
              // map sigma_Vs from [sigma_V_min, sigma_V_max] to [0, 65535] with some saftely margins

              std::set<uint64_t> cset;
              for (int ic = 0; ic < ncenters; ++ic) {
                  uint64_t code{0};
                  double rc[4];
                  for (int d = 0; d < 4; ++d) {
                      double const reduced_coord = (center(ic,d) + ab[d][0])*ab[d][1];
                      rc[d] = reduced_coord;
                      uint16_t const i = std::floor(reduced_coord);
                      uint64_t const i64 = i;
                      code |= (i64 << 16*(3 - d)); // compact into 64bits
                  } // d
                  if (echo > 7) std::printf("# CoW #%i\tcode= %016llx \t%9.1f%9.1f%9.1f%9.1f\n", ic, code, rc[0], rc[1], rc[2], rc[3]);
                  // now enlist all codes in a set and measure its magnitude?
                  cset.insert(code);
              } // ic
              if (echo > 5) std::printf("# CoW set has only %ld of %d elements\n", cset.size(), ncenters);

              // ToDo: run over each element in the set, find which centers are in it,
              //       run over all 3^4 adjacent cells, find the centers in those,
              //       do an order-N^2-comparison to find which pairs are below the thresholds

          } else if ('q' == (*method | 32)) { // "quadratic"
              if (echo > 5) std::printf("# sho_hamiltonian.reduce.centers='%s' --> 'quadratic'\n", method);

              // the naive way: N^2-comparisons.
              std::vector<int> remap(ncenters); std::iota(remap.begin(), remap.end(), 0); // init remap with 0,1,2,3,...
              if (echo > 5) std::printf("# compare %ld pairs from %d expansion centers\n", (size_t(ncenters)*(ncenters - 1))/2, ncenters);
              typedef vector_math::vec<4,double> vec4;
              size_t nremap{0};
              for (int ic = 0; ic < ncenters; ++ic) {
                  vec4 const vec_i = center[ic];
                  for (int jc = ic + 1; jc < ncenters; ++jc) { // triangular loop
                      assert( jc > ic );
                      if (remap[jc] == jc) {
                          vec4 const vec_j = center[jc];
                          vec4 const vec_d = vec_i - vec_j;
                          if (std::abs(vec_d[3]) < threshold_sigma) {
                              if (std::abs(vec_d[2]) < threshold_space) {
                                  if (std::abs(vec_d[1]) < threshold_space) {
                                      if (std::abs(vec_d[0]) < threshold_space) {
                                          // center #ic and #jc are sufficiently close and can be unified
                                          remap[jc] = ic;
                                          center_active[jc] = false; // center #jc will not be referenced, so its projection is not needed
                                          ++nremap;
                                      } // threshold_space
                                  } // threshold_space
                              } // threshold_space
                          } // threshold_sigma
                      } // remap[jc] == jc
                  } // jc
              } // ic
              if (echo > 5) std::printf("# %ld pairs of %d expansion centers marked for remapping\n", nremap, ncenters);

              // run over center_map and apply remap to center_map(:,:,:,0)
              size_t mremap{0};
              for (int ia = 0; ia < natoms; ++ia) {
                  for (int ja = 0; ja < natoms; ++ja) {
                      for (int ip = 0; ip < n_periodic_images; ++ip) {
                          int const jc = center_map(ia,ja,ip,0);
                          int const ic = remap[jc];
                          assert( jc >= ic );
                          mremap += (ic != jc);
                          center_map(ia,ja,ip,0) = ic;
                      } // ip
                  } // ja
              } // ia
              if (echo > 5) std::printf("# %ld pairs of %d expansion centers remapped\n", mremap, ncenters);

          } else { // default method
              if (echo > 5) std::printf("# sho_hamiltonian.reduce.centers='%s' --> 'none' (default)\n", method);
              if ('n' != (*method | 32)) {
                  warn("unknown method sho_hamiltonian.reduce.centers='%s', default to 'none'", method);
              } // 'none'
              if (echo > 5) std::printf("# use the unreduced number of %d expansion centers\n", ncenters);
          } // method
      } // scope: reduce_ncenters


      // perform the projection of the local potential
      std::vector<std::vector<double>> Vcoeffs(ncenters);
      int ncenters_active{0};
      double const scale_potential = control::get("hamiltonian.scale.potential", 1.0);
      if (1.0 != scale_potential) warn("local potential is scaled by %g", scale_potential);
      for (int ic = 0; ic < ncenters; ++ic) {
          if (center_active[ic]) {
              double const sigma_V = center(ic,3);
              int    const numax_V = center(ic,4);
              int const nc = sho_tools::nSHO(numax_V);
              Vcoeffs[ic] = std::vector<double>(nc, 0.0);
              for (int ip = 0; ip < n_periodic_images; ++ip) {
                  std::vector<double> Vcoeff(nc, 0.0);
                  double cnt[3]; set(cnt, 3, center[ic]); add_product(cnt, 3, periodic_image[ip], 1.0);
                  stat += sho_projection::sho_project(Vcoeff.data(), numax_V, cnt, sigma_V, vtot, g, 0); // 0:mute
                  add_product(Vcoeffs[ic].data(), nc, Vcoeff.data(), 1.0);
              } // ip
              // now Vcoeff is represented w.r.t. to Hermite polynomials H_{nx}*H_{ny}*H_{nz} and order_zyx

              stat += sho_potential::normalize_potential_coefficients(Vcoeffs[ic].data(), numax_V, sigma_V, 0); // 0:mute
              // now Vcoeff is represented w.r.t. powers of the Cartesian coords x^{nx}*y^{ny}*z^{nz} and order_Ezyx

              scale(Vcoeffs[ic].data(), Vcoeffs[ic].size(), scale_potential);
              ++ncenters_active;
          } // center_active[ic]
      } // ic
      if (echo > 5) std::printf("# projection performed for %d of %d expansion centers\n", ncenters_active, ncenters);

      prepare_timer.stop(echo);
      // all preparations done, start k-point loop

      int const floating_point_bits = control::get("hamiltonian.floating.point.bits", 64.); // double by default
      bool const single_precision = (32 == floating_point_bits);
      simple_stats::Stats<> time_stats;
      for (int ikp = 0; ikp < nkpoints; ++ikp) {
          std::complex<double> constexpr minus_one = -1;
          std::complex<double> const Bloch_phase[3] = {std::pow(minus_one, 2*kmesh(ikp,0))
                                                     , std::pow(minus_one, 2*kmesh(ikp,1))
                                                     , std::pow(minus_one, 2*kmesh(ikp,2))};
          char x_axis[96]; std::snprintf(x_axis, 95, "# %g %g %g spectrum ", kmesh(ikp,0), kmesh(ikp,1), kmesh(ikp,2));
          double Bloch_phase_real[3];
          bool can_be_real{true};
          for (int d = 0; d < 3; ++d) {
              Bloch_phase_real[d] = Bloch_phase[d].real();
              can_be_real = can_be_real && (std::abs(Bloch_phase[d].imag()) < 2e-16);
//            std::printf("# phase-%c = %g %g\n", 'x'+d, Bloch_phase[d].real(), Bloch_phase[d].imag());
          } // d
          SimpleTimer timer(__FILE__, __LINE__, x_axis, 0);
          #define SOLVE_K_ARGS(BLOCH_PHASE) (natoms, xyzZ, numaxs.data(), sigmas.data(), \
                          n_periodic_images, periodic_image, periodic_shift, \
                          Vcoeffs.data(), center_map, \
                          nB, nBa, offset.data(), \
                          natoms_PAW, xyzZ_PAW, numax_PAW.data(), sigma_PAW.data(), hs_PAW.data(), \
                          BLOCH_PHASE, x_axis, echo)
          if (can_be_real) {
              stat += single_precision ?
                  solve_k<float>  SOLVE_K_ARGS(Bloch_phase_real):
                  solve_k<double> SOLVE_K_ARGS(Bloch_phase_real);
          } else { // replace this else by if(true) to test if real and complex version agree
              stat += single_precision ?
                  solve_k<std::complex<float>>  SOLVE_K_ARGS(Bloch_phase):
                  solve_k<std::complex<double>> SOLVE_K_ARGS(Bloch_phase);
          } // !can_be_real
          #undef SOLVE_K_ARGS
          time_stats.add(timer.stop());
          if (echo > 0) fflush(stdout);
      } // ikp

      if (echo > 3) std::printf("\n# average time per k-point is %.3f +/- %.3f min %.3f max %.3f seconds\n",
                                   time_stats.mean(), time_stats.dev(), time_stats.min(), time_stats.max());

      return stat;
  } // solve

  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_Hamiltonian(int const echo=5) {
      status_t stat(0);

      auto const vtotfile = control::get("sho_potential.test.vtot.filename", "vtot.dat"); // vtot.dat can be created by potential_generator or self_consistency
      int dims[3] = {0, 0, 0};
      std::vector<double> vtot; // total smooth potential
      stat += sho_potential::load_local_potential(vtot, dims, vtotfile, echo);

      auto const geo_file = control::get("geometry.file", "atoms.xyz");
      view2D<double> xyzZ;
      int natoms{0}; // number of atoms
      double cell[3] = {0, 0, 0}; // rectangular cell
      int8_t bc[3] = {-7, -7, -7}; // boundary conditions
      { // scope: read atomic positions
          stat += geometry_analysis::read_xyz_file(xyzZ, natoms, geo_file, cell, bc, 0);
          if (echo > 2) std::printf("# found %d atoms in file \"%s\" with cell=[%.3f %.3f %.3f] %s and bc=[%d %d %d]\n",
                              natoms, geo_file, cell[0]*Ang, cell[1]*Ang, cell[2]*Ang, _Ang, bc[0], bc[1], bc[2]);
      } // scope

      real_space::grid_t g(dims);
      g.set_boundary_conditions(bc);
      g.set_grid_spacing(cell[0]/g[0], cell[1]/g[1], cell[2]/g[2]);
      if (echo > 1) {
          std::printf("# use  %g %g %g %s grid spacing\n", g.h[0]*Ang, g.h[1]*Ang, g.h[2]*Ang, _Ang);
          std::printf("# cell is  %g %g %g %s\n", g.h[0]*g[0]*Ang, g.h[1]*g[1]*Ang, g.h[2]*g[2]*Ang, _Ang);
      } // echo
      
      int const nkpoints = control::get("hamiltonian.test.kpoints", 17.);
      auto const kpointdir = int(control::get("hamiltonian.test.kpoint.direction", 0.)) % 3;
      view2D<double> kmesh(nkpoints, 4, 0.0);
      for (int ikp = 0; ikp < nkpoints; ++ikp) {
          kmesh(ikp,kpointdir) = ikp*0.5/(nkpoints - 1.); // bandstructure e.g. in x-direction
      } // ikp

      stat += solve(natoms, xyzZ, g, vtot.data(), nkpoints, kmesh, 0, 0, 0, 0, echo);

      return stat;
  } // test_Hamiltonian
  
  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_Hamiltonian(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace sho_hamiltonian
