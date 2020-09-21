#include <cstdio> // printf
#include <cmath> // std::sqrt
#include <algorithm> // std::max
#include <complex> // std::complex<real_t>
#include <vector> // std::vector<T>
#include <cassert> // assert
#include <set> // std::set<key>
#include <numeric> // std::iota

#include "pw_hamiltonian.hxx"

#include "sho_potential.hxx" // ::load_local_potential
#include "geometry_analysis.hxx" // ::read_xyz_file, ::fold_back
#include "control.hxx" // ::get
#include "display_units.h" // eV, _eV, Ang, _Ang, Kelvin, _Kelvin
#include "real_space.hxx" // ::grid_t
#include "sho_tools.hxx" // ::nSHO, ::n1HO, ::order_*, ::SHO_index_t, ::construct_label_table
#include "sho_projection.hxx" // ::sho_project
#include "boundary_condition.hxx" // Isolated_Boundary, ::periodic_image
#include "sho_overlap.hxx" // ::
#include "data_view.hxx" // view2D<T>, view3D<T>, view4D<T>
#include "inline_tools.hxx" // align<nbits>
#include "simple_math.hxx" // ::random<real_t>
#include "vector_math.hxx" // ::vec<N,T>
#include "hermite_polynomial.hxx" // hermite_polys
#include "simple_stats.hxx" // ::Stats
#include "simple_timer.hxx" // SimpleTimer
#include "dense_solver.hxx" // ::solve
#include "unit_system.hxx" // ::energy_unit
#include "fourier_transform.hxx" // ::fft


namespace pw_hamiltonian {
  // computes Hamiltonian matrix elements for plane waves
  // including a PAW non-local contribution
  
  class PlaneWave {
    public:
      double g2; // length of the vector square == plane wave kinetic energy in Rydberg
      int16_t x, y, z;
      PlaneWave() : g2(0), x(0), y(0), z(0) {}
      PlaneWave(int const ix, int const iy, int const iz, double const len2) : g2(len2), x(ix), y(iy), z(iz) {}
    private:
  }; // class PlaneWave


  void Hermite_Gauss_projectors(double pzyx[], int const numax, double const sigma, double const gv[3]) {

        view2D<double> Hermite_Gauss(3, sho_tools::n1HO(numax));
        for(int d = 0; d < 3; ++d) {
            double const x = -gv[d]*sigma; // -x because the Fourier transform of HG_n(x) is (-1)^n HG_n(k) = HG(-k) [for sigma=1]
            hermite_polys(Hermite_Gauss[d], x, numax); // unnormalized Hermite-Gauss functions
        } // d

        { // scope: normalize Hermite_Gauss functions
            double nfactorial{1};
            for(int n = 0; n <= numax; ++n) {
                double const nrmf = std::sqrt(sigma/(constants::sqrtpi*nfactorial));
                for(int d = 0; d < 3; ++d) {
                    Hermite_Gauss(d,n) *= nrmf; // scale
                } // d
                nfactorial *= (n + 1)*0.5; // update nfactorial
            } // n
        } // scope

        // generate projection coefficient for factorized SHO projectors in zyx_order
        int lb{0};
        for(        int lz = 0; lz <= numax;           ++lz) {
            for(    int ly = 0; ly <= numax - lz;      ++ly) {
                for(int lx = 0; lx <= numax - lz - ly; ++lx) {
                    pzyx[lb] = Hermite_Gauss(0,lx) *
                               Hermite_Gauss(1,ly) *
                               Hermite_Gauss(2,lz);
                    ++lb;
                } // lx
            } // ly
        } // lz
        assert( sho_tools::nSHO(numax) == lb );

  } // Hermite_Gauss_projectors


  inline double Fourier_Gauss_factor_3D(double const sigma=1) {
      // \sqrt(2 pi sigma) ^3
      return pow3(constants::sqrt2 * constants::sqrtpi);
  } //  Fourier_Gauss_factor_3D

  template<typename complex_t, typename real_t>
  status_t solve_k(double const ecut // plane wave cutoff energy
          , double const reci[3][4] // reciprocal cell matrix
          , view4D<double> const & Vcoeff // <G|V|G'> = potential(0:1or3,iGz-jGz,iGy-jGy,iGx-jGx)
          , int const nG[3] // half the number of plane wave entries in the potential
          , double const norm_factor // normalization for the cell volume
          , int const natoms_PAW // number of PAW centers
          , view2D<double const> const & xyzZ_PAW // (natoms_PAW, 4) positions of PAW centers, 4th component not used
          , int    const numax_PAW[] // spreads of the SHO-type PAW projectors
          , double const sigma_PAW[] // cutoffs of the SHO-type PAW projectors
          , view3D<double> const hs_PAW[] // [natoms_PAW](2, nprj, nprj) // PAW Hamiltonian correction and charge-deficit
          , double const kpoint[3] // vector inside the Brillouin zone, all 3 components in [-.5, .5]
          , char const *const x_axis // display this string in front of the Hamiltonian eigenvalues
          , int & nPWs // export number of plane waves
          , int const echo=0) { // log-level

      int boxdim[3], max_PWs{1};
      for(int d = 0; d < 3; ++d) {
          boxdim[d] = std::ceil(ecut/pow2(reci[d][d]));
          max_PWs *= (boxdim[d] + 1 + boxdim[d]);
      } // d
      std::vector<PlaneWave> pw_basis(max_PWs);
      { // scope: generate set of plane waves
          int iB{0}, outside{0};
          for(int iGz = -boxdim[2]; iGz <= boxdim[2]; ++iGz) {
          for(int iGy = -boxdim[1]; iGy <= boxdim[1]; ++iGy) {
          for(int iGx = -boxdim[0]; iGx <= boxdim[0]; ++iGx) {
              double const gv[3] = {(iGx + kpoint[0])*reci[0][0],
                                    (iGy + kpoint[1])*reci[1][1],
                                    (iGz + kpoint[2])*reci[2][2]}; // ToDo: diagonal reci assumed here
              double const g2 = pow2(gv[0]) + pow2(gv[1]) + pow2(gv[2]);
              if (g2 < 2*ecut) {
                  pw_basis[iB] = PlaneWave(iGx, iGy, iGz, g2);
                  ++iB; // count
              } else { ++outside;
              } // within energy ball
          } // iGx
          } // iGy
          } // iGz
          int const number_of_PWs = iB;
          assert( number_of_PWs + outside == max_PWs );
          pw_basis.resize(number_of_PWs); // resize() or shrink_to_fit()
      } // scope: generate set of plane waves
      int const nB = pw_basis.size();
      nPWs = nB; // export the number of plane waves used for statistics
#ifdef DEVEL
      {   complex_t c(1); real_t r{1};
          if (echo > 3) printf("\n\n# start %s<%s, real_t=%s> nB=%d\n", 
                  __func__, complex_name(c), complex_name(r), nB);
      }
#endif

      // now we could sort the plane waves by their g2 length, optional
      bool constexpr sort_PWs = true;
      if (sort_PWs) {
          auto compare_lambda = [](PlaneWave const & lhs, PlaneWave const & rhs) { return lhs.g2 < rhs.g2; };
          std::sort(pw_basis.begin(), pw_basis.end(), compare_lambda);
#ifdef DEVEL
          if (echo > 6) { // show in ascending order
              printf("# %s|k+G|^2 =", x_axis);
              for(int i = 0; i < std::min(5, nPWs); ++i) {
                  printf(" %.4f", pw_basis[i].g2);
              } // i
              printf(" ...");
              for(int i = -std::min(2, nPWs); i < 0; ++i) {
                  printf(" %.3f", pw_basis[nPWs + i].g2);
              } // i
              printf(" Ry\n");
          } // echo
#endif
      } // sort_PWs

      //
      // now construct the Hamiltonian:
      //   -- kinetic energy,             c.f. sho_overlap::test_simple_crystal
      //   -- local potential,            c.f. sho_potential::test_local_potential_matrix_elements (method=2)
      //   -- non-local PAW Hamiltonian contributions
      // and the overlap operator:
      //   -- SHO basis function overlap, c.f. sho_overlap::test_simple_crystal
      //   -- non-local PAW charge-deficit contributions
      //

      int constexpr S=0, H=1; // static indices for S:overlap matrix, H:Hamiltonian matrix

      std::vector<int> offset(natoms_PAW + 1, 0); // prefetch sum
      for(int ka = 0; ka < natoms_PAW; ++ka) {
          int const nSHO = sho_tools::nSHO(numax_PAW[ka]);
          offset[ka + 1] = offset[ka] + nSHO;
      } // ka
      int const nC = offset[natoms_PAW]; // total number of PAW coefficients

      // PAW contributions to H_{ij} = P_{ik} h_{kl} P^*_{jl} = Ph_{il} P^*_{jl}
      //                  and S_{ij} = P_{ik} s_{kl} P^*_{jl} = Ps_{il} P^*_{jl}

      view2D<complex_t> P_jl(nB, nC, 0.0); //  <k+G_i|\tilde p_{la lb}>
      std::vector<double>   P2_l(nC, 0.0); // |<k+G_i|\tilde p_{la lb}>|^2
      view3D<complex_t> Psh_il(2, nB, nC, 0.0); // atom-centered PAW matrices multiplied to P_jl
      for(int jB = 0; jB < nB; ++jB) {
          auto const & i = pw_basis[jB];
          double const gv[3] = {(i.x + kpoint[0])*reci[0][0],
                                (i.y + kpoint[1])*reci[1][1],
                                (i.z + kpoint[2])*reci[2][2]}; // ToDo: diagonal reci-matrix assumed here
          for(int ka = 0; ka < natoms_PAW; ++ka) {
              int const nSHO = sho_tools::nSHO(numax_PAW[ka]);

              // phase factor related to the atomic position
              double const arg = xyzZ_PAW(ka,0)*gv[0] + xyzZ_PAW(ka,1)*gv[1] + xyzZ_PAW(ka,2)*gv[2];
              std::complex<double> const phase(std::cos(arg), std::sin(arg));

              {   std::vector<double> pzyx(nSHO);
                  Hermite_Gauss_projectors(pzyx.data(), numax_PAW[ka], sigma_PAW[ka], gv);
                  auto const fgf = Fourier_Gauss_factor_3D(sigma_PAW[ka]);

                  for(int lb = 0; lb < nSHO; ++lb) {
                      int const lC = offset[ka] + lb;
                      P_jl(jB,lC) = complex_t(phase * pzyx[lb] * fgf * norm_factor);
                      P2_l[lC] += pow2(pzyx[lb]);
                  } // lb
              }

              // multiply P from left to hs_PAW (block diagonal --> ka == la)
              auto const iB = jB;
              for(int lb = 0; lb < nSHO; ++lb) {
                  int const iC = offset[ka] + lb; // global PAW coefficient index
                  complex_t s(0), h(0);
                  for(int kb = 0; kb < nSHO; ++kb) { // contract
                      int const kC = offset[ka] + kb;
                      auto const p = P_jl(iB,kC); // ToDo: do we need to conjugate here instead of below?
                      // Mind that hs_PAW matrices are ordered {0:H,1:S}
                      s += p * real_t(hs_PAW[ka](1,kb,lb));
                      h += p * real_t(hs_PAW[ka](0,kb,lb));
                  } // kb
                  Psh_il(S,iB,iC) = s;
                  Psh_il(H,iB,iC) = h;
              } // lb
          } // ka
      } // jB

#ifdef DEVEL
      if (echo > 29) {
          for(int la = 0; la < natoms_PAW; ++la) {
              printf("# ^2-norm of the projector of atom #%i ", la);
              int const nSHO = sho_tools::nSHO(numax_PAW[la]);
              for(int lb = 0; lb < nSHO; ++lb) {
                  int const lC = offset[la] + lb;
                  printf(" %g", P2_l[lC]);
              } // lb
              printf("\n"); fflush(stdout);
          } // la
      } // echo
#endif

      // PAW projection matrix ready

      // allocate bulk memory for overlap and Hamiltonian
      int const nBa = align<4>(nB); // memory aligned main matrix stride
      view3D<complex_t> SHm(2, nB, nBa, complex_t(0)); // get memory for 0:Overlap S and 1:Hamiltonian matrix H

#ifdef DEVEL
      if (echo > 3) printf("# assume dimensions of Vcoeff(%d, %d, %d, %d)\n", 2, Vcoeff.dim2(), Vcoeff.dim1(), Vcoeff.stride());
      assert( nG[0] <= Vcoeff.stride() );
      assert( nG[1] <= Vcoeff.dim1() );
      assert( nG[2] <= Vcoeff.dim2() );
      
      double const scale_k = control::get("pw_hamiltonian.scale.kinetic", 1.0);
      double const scale_p = control::get("pw_hamiltonian.scale.potential", 1.0);
      if (1 != scale_k) warn("kinetic energy is scaled by %g", scale_k);
      if (1 != scale_p) warn("local potential is scaled by %g", scale_p);
      real_t const scale_h = control::get("pw_hamiltonian.scale.nonlocal.h", 1.0);
      real_t const scale_s = control::get("pw_hamiltonian.scale.nonlocal.s", 1.0);
      if (1 != scale_h || 1 != scale_s) warn("scale PAW contributions to H and S by %g and %g, respectively", scale_h, scale_s);
#else
      real_t constexpr scale_h = 1, scale_s = 1
      double constexpr scale_k = 1, scale_p = 1;
#endif
      double const kinetic = 0.5 * scale_k; // prefactor of kinetic energy in Hartree atomic units
      real_t const localpot = scale_p / (nG[0]*nG[1]*nG[2]);

      for(int iB = 0; iB < nB; ++iB) {      auto const & i = pw_basis[iB];

          SHm(S,iB,iB) += 1; // unity overlap matrix, plane waves are orthogonal
          SHm(H,iB,iB) += kinetic*i.g2; // kinetic energy contribution

          for(int jB = 0; jB < nB; ++jB) {  auto const & j = pw_basis[jB];

              // add the contribution of the local potential
              int const iVx = i.x - j.x, iVy = i.y - j.y, iVz = i.z - j.z;
              if ((2*std::abs(iVx) < nG[0]) && (2*std::abs(iVy) < nG[1]) && (2*std::abs(iVz) < nG[2])) {
                  int const imz = (iVz + nG[2])%nG[2], imy = (iVy + nG[1])%nG[1], imx = (iVx + nG[0])%nG[0];
                  auto const V_re = Vcoeff(0, imz, imy, imx), V_im = Vcoeff(1, imz, imy, imx);
                  SHm(H,iB,jB) += localpot*std::complex<real_t>(V_re, V_im);
              } // inside range

              // PAW contributions to H_{ij} = Ph_{il} P^*_{jl}
              //                  and S_{ij} = Ps_{il} P^*_{jl}

              // inner part of a matrix-matrix multiplication
              { // scope: add non-local contributions
                  complex_t s(0), h(0);
                  for(int lC = 0; lC < nC; ++lC) { // contract
                      auto const p = conjugate(P_jl(jB,lC)); // needs a conjugation, ToDo: here or above?
                      s += Psh_il(S,iB,lC) * p;
                      h += Psh_il(H,iB,lC) * p;
                  } // lC
                  SHm(S,iB,jB) += s * scale_s;
                  SHm(H,iB,jB) += h * scale_h;
              } // scope: add non-local contributions

          } // jB
      } // iB

      return dense_solver::solve<complex_t, real_t>(SHm, x_axis, echo);
  } // solve_k







  status_t solve(int const natoms_PAW // number of PAW atoms
          , view2D<double const> const & xyzZ // (natoms, 4)
          , real_space::grid_t const & g // Cartesian grid descriptor for vtot
          , double const *const vtot // total effective potential on grid
          , double const *const sigma_prj // =nullptr
          , int    const *const numax_prj // =nullptr
          , double *const *const atom_mat // =nullptr
          , int const echo) { // log-level

      status_t stat(0);
      SimpleTimer prepare_timer(__FILE__, __LINE__, "prepare", 0);
      
      for(int ia = 0; ia < natoms_PAW; ++ia) {
          if (echo > 0) {
              printf("# atom#%i \tZ=%g \tposition %12.6f%12.6f%12.6f %s\n",
                  ia, xyzZ(ia,3), xyzZ(ia,0)*Ang, xyzZ(ia,1)*Ang, xyzZ(ia,2)*Ang, _Ang);
          } // echo
      } // ia

      double const cell_matrix[3][4] = {{g[0]*g.h[0], 0, 0, 0}, {0, g[1]*g.h[1], 0, 0}, {0, 0, g[2]*g.h[2], 0}};
      double reci_matrix[3][4];
      auto const cell_volume = simple_math::invert(3, reci_matrix[0], 4, cell_matrix[0], 4);
      scale(reci_matrix[0], 12, 2*constants::pi); // scale by 2\pi
      if (echo > 0) printf("# cell volume is %g %s^3\n", cell_volume*pow3(Ang),_Ang);
      if (echo > 0) {
          auto constexpr sqRy = 1; auto const _sqRy = "sqRy";
          printf("# cell matrix in %s\t\tand\t\treciprocal matrix in %s:\n", _Ang, _sqRy);
          for(int d = 0; d < 3; ++d) {
              printf("# %12.6f%12.6f%12.6f    \t%12.6f%12.6f%12.6f\n",
                cell_matrix[d][0]*Ang,  cell_matrix[d][1]*Ang,  cell_matrix[d][2]*Ang,
                reci_matrix[d][0]*sqRy, reci_matrix[d][1]*sqRy, reci_matrix[d][2]*sqRy);
          } // d
          printf("\n");
      } // echo
      stat += (abs(cell_volume) > 1e-15);
      double const svol = 1./std::sqrt(cell_volume); // normalization factor for plane waves
      if (echo > 1) printf("# normalization factor for plane waves is %g a.u.\n", svol);

      char const *_ecut_u{nullptr};
      auto const ecut_u = unit_system::energy_unit(control::get("pw_hamiltonian.cutoff.energy.unit", "Ha"), &_ecut_u);
      auto const ecut = control::get("pw_hamiltonian.cutoff.energy", 11.)/ecut_u; // 11 Ha =~= 300 eV cutoff energy
      if (echo > 1) printf("# pw_hamiltonian.cutoff.energy=%.3f %s corresponds to %.3f^2 Ry or %.2f %s\n", 
                              ecut*ecut_u, _ecut_u, std::sqrt(2*ecut), ecut*eV,_eV);

//       int const nG[3] = {g[0]/2 - 1, g[1]/2 - 1, g[2]/2 - 1}; // 2*nG <= g
//       view4D<std::complex<double>> Vcoeffs(2*nG[2]+1, 2*nG[1]+1, 2*nG[0]+2, 1, 0.0);
      int const nG[3] = {g[0], g[1], g[2]}; //
      view4D<double> Vcoeffs(2,nG[2],nG[1],nG[0], 0.0);

      auto const fft_stat = fourier_transform::fft(Vcoeffs(0,0,0), Vcoeffs(1,0,0), vtot, nG, echo);
      if (0 == fft_stat) {
          if (echo > 5) printf("# used FFT to transform local potential into space of plane-wave differences\n");
      } else { // scope: Fourier transform the local potential "by hand"
          if (echo > 5) printf("# FFT failed, transform local potential manually\n");
          double const two_pi = 2*constants::pi;
          double const tpi_g[] = {two_pi/g[0], two_pi/g[1], two_pi/g[2]};
          for        (int iGz = 0; iGz < nG[2]; ++iGz) {  auto const Gz = iGz*tpi_g[2];
              for    (int iGy = 0; iGy < nG[1]; ++iGy) {  auto const Gy = iGy*tpi_g[1];
                  for(int iGx = 0; iGx < nG[0]; ++iGx) {  auto const Gx = iGx*tpi_g[0];
                      double re{0}, im{0};
                      for(int iz = 0; iz < g[2]; ++iz) {
                      for(int iy = 0; iy < g[1]; ++iy) {
                      for(int ix = 0; ix < g[0]; ++ix) {
                          double const arg = Gz*iz + Gy*iy + Gx*ix;
                          double const V = vtot[(iz*g[1] + iy)*g[0] + ix];
                          re += V * std::cos(arg);
                          im += V * std::sin(arg);
                      } // ix
                      } // iy
                      } // iz
                      Vcoeffs(0,iGz,iGy,iGx) = re; // store
                      Vcoeffs(1,iGz,iGy,iGx) = im; // store
                  } // iGx
              } // iGy
          } // iGz
          
          if (echo > 6) {
              auto const V0 = Vcoeffs(nG[2],nG[1],nG[0],0);
              printf("# 000-Fourier coefficient of the potential is %g %g %s\n", 
                                std::real(V0)*eV, std::imag(V0)*eV, _eV); 
          } // echo
      } // scope

      
      // prepare for the PAW contributions: find the projection coefficient matrix P = <\chi3D_{ia ib}|\tilde p_{ka kb}>
      std::vector<int32_t> numax_PAW(natoms_PAW,  3);
      std::vector<double>  sigma_PAW(natoms_PAW, .5);
      double maximum_sigma_PAW{0};
      std::vector<view3D<double>> hs_PAW(natoms_PAW); // atomic matrices for charge-deficit and Hamiltonian corrections
      for(int ka = 0; ka < natoms_PAW; ++ka) {
          if (numax_prj) numax_PAW[ka] = numax_prj[ka];
          if (sigma_prj) sigma_PAW[ka] = sigma_prj[ka];
          maximum_sigma_PAW = std::max(maximum_sigma_PAW, sigma_PAW[ka]);
          int const nSHO = sho_tools::nSHO(numax_PAW[ka]); // number of projectors
          if (atom_mat) {
              hs_PAW[ka] = view3D<double>(atom_mat[ka], nSHO, nSHO); // wrap input
          } else {
              hs_PAW[ka] = view3D<double>(2, nSHO, nSHO, 0.0); // get memory and initialize zero
          } // atom_mat
          // for both matrices: L2-normalization w.r.t. SHO projector functions is important
      } // ka

      
      
      prepare_timer.stop(echo);
      // all preparations done, start k-point loop

      auto const floating_point_bits = int(control::get("pw_hamiltonian.floating.point.bits", 64.)); // double by default
      auto const nkpoints = int(control::get("pw_hamiltonian.test.kpoints", 17.));
      simple_stats::Stats<double> nPW_stats, tPW_stats;
      for(int ikp = 0; ikp < nkpoints; ++ikp) {
          char x_axis[96]; std::snprintf(x_axis, 95, "# %.6f spectrum ", ikp*.5/(nkpoints - 1.));
          double const kpoint[3] = {ikp*.5/(nkpoints - 1.), 0, 0}; //
          SimpleTimer timer(__FILE__, __LINE__, x_axis, 0);

          bool can_be_real{false}; // real only with inversion symmetry
          if (can_be_real) {
              error("PW only implemented with complex");
          } else { // replace this "else" by "if(true)" to test if real and complex version agree
              int nPWs{0};
              if (32 == floating_point_bits) {
                  stat += solve_k<std::complex<float>, float>(ecut, reci_matrix, Vcoeffs, nG, svol,
                          natoms_PAW, xyzZ, numax_PAW.data(), sigma_PAW.data(), hs_PAW.data(),
                          kpoint, x_axis, nPWs, echo);
              } else {
                  stat += solve_k<std::complex<double>, double>(ecut, reci_matrix, Vcoeffs, nG, svol,
                          natoms_PAW, xyzZ, numax_PAW.data(), sigma_PAW.data(), hs_PAW.data(),
                          kpoint, x_axis, nPWs, echo);
              } // floating_point_bits
              nPW_stats.add(nPWs);
          } // !can_be_real
          if (echo > 0) fflush(stdout);
          tPW_stats.add(timer.stop());
      } // ikp
      
      if (echo > 3) printf("\n# average number of plane waves is %.3f +/- %.3f min %g max %g\n",
                      nPW_stats.avg(), nPW_stats.var(), nPW_stats.min(), nPW_stats.max());
      if (echo > 3) printf("# average time per k-point is %.3f +/- %.3f min %.3f max %.3f seconds\n",
                      tPW_stats.avg(), tPW_stats.var(), tPW_stats.min(), tPW_stats.max());

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
      view2D<double> xyzZ_noconst;
      int natoms{0};
      double cell[3] = {0, 0, 0}; 
      int bc[3] = {-7, -7, -7};
      { // scope: read atomic positions
          stat += geometry_analysis::read_xyz_file(xyzZ_noconst, natoms, geo_file, cell, bc, 0);
          if (echo > 2) printf("# found %d atoms in file \"%s\" with cell=[%.3f %.3f %.3f] %s and bc=[%d %d %d]\n",
                              natoms, geo_file, cell[0]*Ang, cell[1]*Ang, cell[2]*Ang, _Ang, bc[0], bc[1], bc[2]);
      } // scope
      view2D<double const> const xyzZ(xyzZ_noconst.data(), xyzZ_noconst.stride()); // wrap for simpler usage

      real_space::grid_t g(dims);
      g.set_boundary_conditions(bc);
      g.set_grid_spacing(cell[0]/g[0], cell[1]/g[1], cell[2]/g[2]);
      if (echo > 1) printf("# use  %g %g %g %s grid spacing\n", g.h[0]*Ang, g.h[1]*Ang, g.h[2]*Ang, _Ang);
      if (echo > 1) printf("# cell is  %g %g %g %s\n", g.h[0]*g[0]*Ang, g.h[1]*g[1]*Ang, g.h[2]*g[2]*Ang, _Ang);
      
      stat += solve(natoms, xyzZ, g, vtot.data(), 0, 0, 0, echo);

      return stat;
  } // test_Hamiltonian

  status_t test_Hermite_Gauss_normalization(int const echo=5, int const numax=3) {
      // we check if the Hermite_Gauss_projectors are L2-normalized 
      // when evaluated on a dense Cartesian mesh (in reciprocal space).
      // this is useful information as due to the Plancherel theorem,
      // this implies the right scaling of the projectors in both spaces.
      status_t stat(0);

      auto const nSHO = sho_tools::nSHO(numax);
      std::vector<double> pzyx(nSHO), p2(nSHO, 0.0);
      double constexpr dg = 0.125, d3g = pow3(dg);
      int constexpr ng = 40;
      for(int iz = -ng; iz <= ng; ++iz) {
      for(int iy = -ng; iy <= ng; ++iy) {
      for(int ix = -ng; ix <= ng; ++ix) {
          double const gv[3] = {dg*ix, dg*iy, dg*iz};
          Hermite_Gauss_projectors(pzyx.data(), numax, 1.0, gv);
          add_product(p2.data(), nSHO, pzyx.data(), pzyx.data()); // += |p^2|
      }}} // ig
      printf("\n# %s: norms ", __func__);
      for(int ip = 0; ip < nSHO; ++ip) {
          printf(" %g", p2[ip]*d3g);
      } // ip
      printf("\n");

      return stat;
  } // test_Hermite_Gauss_normalization
  
  status_t all_tests(int const echo) {
      status_t status(0);
      status += test_Hamiltonian(echo);
      status += test_Hermite_Gauss_normalization(echo);
      return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace pw_hamiltonian
