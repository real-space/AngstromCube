#include <cstdio> // printf, fflush, stdout
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
#include "sho_tools.hxx" // ::nSHO, ::n1HO, ::order_*, ::construct_label_table
#include "sho_projection.hxx" // ::sho_project
#include "boundary_condition.hxx" // Isolated_Boundary, ::periodic_image
#include "sho_overlap.hxx" // ::
#include "data_view.hxx" // view2D<T>, view3D<T>, view4D<T>
#include "data_list.hxx" // data_list<T>
#include "inline_tools.hxx" // align<nbits>
#include "simple_math.hxx" // ::random<real_t>
#include "vector_math.hxx" // ::vec<N,T>
#include "hermite_polynomial.hxx" // hermite_polys
#include "simple_stats.hxx" // ::Stats
#include "simple_timer.hxx" // SimpleTimer
#include "dense_solver.hxx" // ::solve
#include "unit_system.hxx" // ::energy_unit
#include "fourier_transform.hxx" // ::fft

#include "davidson_solver.hxx" // ::eigensolve
#include "conjugate_gradients.hxx" // ::eigensolve
#include "dense_operator.hxx" // ::dense_operator_t
#include "fermi_distribution.hxx" // ::Fermi_level

#ifdef DEVEL
    #include "print_tools.hxx" // printf_vector(fmt, vec, n [, final, scale, add])
#endif // DEVEL

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
                nfactorial *= (n + 1)*0.5; // update nfactorial * 2^-n
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


  inline double Fourier_Gauss_factor_3D() {
      // \sqrt(2 pi)^3
      return pow3(constants::sqrt2 * constants::sqrtpi);
  } //  Fourier_Gauss_factor_3D

  
  template<typename complex_t>
  status_t iterative_solve(view3D<complex_t> const & SHm, char const *x_axis="",
      int const echo=0, int const nbands=10, float const direct_ratio=2) {
      //
      // An iterative solver using the Davidson method
      //
      status_t stat(0);
      int constexpr H=1, S=0;

      int const nPW  = SHm.dim1();
      int const nPWa = SHm.stride();

      if (echo > 6) {
          printf("# %s for the lowest %d bands of a %d x %d Hamiltonian (stride %d)\n",
                    __func__, nbands, nPW, nPW, nPWa);
          fflush(stdout);
      } // echo

      // create nbands start waves
      view2D<complex_t> waves(nbands, nPW, complex_t(0)); // get memory
      
      if (nbands > nPW) {
          warn("tried to find %d bands in a basis set with %d plane waves", nbands, nPW);
          return -1;
      } // enough plane waves?

      { // scope: fill start wave functions with good, linear independent start vectors
          int const nsubspace = std::max(nbands, int(direct_ratio*nbands));
          int const nsub = std::min(nsubspace, nPW);
          view3D<complex_t> SHmat_b(2, nsub, align<2>(nsub));
          for(int ib = 0; ib < nsub; ++ib) {
              set(SHmat_b(S,ib), nsub, SHm(S,ib));
              set(SHmat_b(H,ib), nsub, SHm(H,ib));
          } // ib

          // assume that plane waves are energy ordered and solve for the nsub * nsub sub-Hamiltonian
          if (echo > 6) { printf("# %s get %d start waves from diagonalization of a %d x %d Hamiltonian (stride %d)\n",
                                  __func__, nbands, nsub, nsub, SHmat_b.stride()); fflush(stdout); }
          auto const stat_eig = dense_solver::solve(SHmat_b, "# start waves "); // mute
          if (stat_eig != 0) {
              warn("diagonalization of the %d x %d sub-Hamiltonian returned status= %i", nsub, nsub, int(stat_eig));
              return stat_eig;
          } // error?
          stat += stat_eig;
          // after the dense solver, the Hamiltonian matrix contains the eigenvectors

          if (echo > 6) { printf("# %s copy %d start waves from eigenstates of the %d x %d Hamiltonian (stride %d)\n",
                                  __func__, nbands, nsub, nsub, SHmat_b.stride()); fflush(stdout); }
          // copy eigenvectors into low PW coefficients of start waves
          for(int ib = 0; ib < nbands; ++ib) {
              set(waves[ib], nsub, SHmat_b(H,ib));
          } // ib
          
      } // scope
      

      if (echo > 6) { printf("# %s construct a dense operator for the %d x %d Hamiltonian (stride %d)\n",
                              __func__, nPW, nPW, nPWa); fflush(stdout); }
      // construct a dense-matrix operator op
      dense_operator::dense_operator_t<complex_t> const op(nPW, nPWa, SHm(H,0), SHm(S,0));

      char const method = *control::get("pw_hamiltonian.iterative.solver", "Davidson") | 32;

      // solve using the Davidson eigensolver
      std::vector<double> eigvals(nbands);
      status_t stat_slv(0);
      if ('c' == method) {
          int const nit = control::get("pw_hamiltonian.max.cg.iterations", 3.);
          if (echo > 6) { printf("# %s envoke CG solver with max. %d iterations\n", __func__, nit); fflush(stdout); }
          for(int it = 0; it < nit && (0 == stat_slv); ++it) {
              if (echo > 6) { printf("# %s envoke CG solver, outer iteration #%i\n", __func__, it); fflush(stdout); }
              stat_slv = conjugate_gradients::eigensolve(waves.data(), eigvals.data(), nbands, op, echo - 10);
              stat += stat_slv;
          } // it
      } else { // method
          int const nit = control::get("davidson_solver.max.iterations", 1.);
          if (echo > 6) { printf("# %s envoke Davidson solver with max. %d iterations\n", __func__, nit); fflush(stdout); }
          stat_slv = davidson_solver::eigensolve(waves.data(), eigvals.data(), nbands, op, echo - 10, 2.5f, nit, 1e-2);
      } // method

      if (0 == stat_slv) {
          if (echo > 2) {
              dense_solver::display_spectrum(eigvals.data(), nbands, x_axis, eV, _eV);
              if (echo > 4) printf("# lowest and highest eigenvalue of the Hamiltonian matrix is %g and %g %s, respectively\n", 
                                      eigvals[0]*eV, eigvals[nbands - 1]*eV, _eV);
              fflush(stdout);
          } // echo
      } else {
          warn("Davidson solver for the plane wave Hamiltonian failed with status= %i", int(stat_slv));
      } // stat_slv

      return stat;
  } // iterative_solve
  
  
  
  template<typename complex_t>
  status_t solve_k(double const ecut // plane wave cutoff energy
          , double const reci[3][4] // reciprocal cell matrix
          , view4D<double> const & Vcoeff // <G|V|G'> = potential(0:1or3,iGz-jGz,iGy-jGy,iGx-jGx)
          , int const nG[3] // half the number of plane wave entries in the potential
          , double const norm_factor // normalization for the cell volume
 
          , int const natoms_PAW // number of PAW centers
          , view2D<double const> const & xyzZ_PAW // (natoms_PAW, 4) positions of PAW centers, 4th component not used
          , double const grid_offset[3] // origin of atomic positions w.r.t. the local potential
          , int    const numax_PAW[] // spreads of the SHO-type PAW projectors
          , double const sigma_PAW[] // cutoffs of the SHO-type PAW projectors
          , view3D<double> const hs_PAW[] // [natoms_PAW](2, nprj, nprj) // PAW Hamiltonian correction and charge-deficit

          , double const kpoint[3] // vector inside the Brillouin zone, all 3 components in [-.5, .5]
          , char const *const x_axis // display this string in front of the Hamiltonian eigenvalues
          , int & nPWs // export number of plane waves
          , int const echo=0 // log-level
          , int const nbands=10 // number of bands (needed for iterative solver)
          , float const direct_ratio=2.f // try to get good start waves by envoking a dense solver on a (2.0*nbands)x(2.0*nbands) sub-Hamiltonian
      ) { // number of bands (needed by iterative solver)

      using real_t = decltype(std::real(complex_t(1)));
            
      int boxdim[3], max_PWs{1};
      for(int d = 0; d < 3; ++d) {
          boxdim[d] = std::ceil(ecut/pow2(reci[d][d]));
          if (echo > 17) printf("# %s reci[%d] = %.6f %.6f %.6f sqRy, kpoint=%9.6f\n", __func__, d, reci[d][0], reci[d][1], reci[d][2], kpoint[d]);
          max_PWs *= (boxdim[d] + 1 + boxdim[d]);
      } // d
      if (echo > 11) printf("# %s boxdim = %d x %d x %d, E_cut= %g Ry\n", __func__, boxdim[0], boxdim[1], boxdim[2], 2*ecut);
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
#ifdef FULL_DEBUG
              if (echo > 91) printf("# %s suggest PlaneWave(%d,%d,%d) with |k+G|^2= %g Ry inside= %d\n", __func__, iGx,iGy,iGz, g2, g2 < 2*ecut);
#endif // FULL_DEBUG
              if (g2 < 2*ecut) {
                  pw_basis[iB] = PlaneWave(iGx, iGy, iGz, g2);
                  ++iB; // count
              } else {
                  ++outside;
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
      if (echo > 6) {
          float const minutes = 1.5e-11*(2 + (sizeof(real_t) > 4))*pow3(1.*nB), seconds = 60*(minutes - int(minutes));
          printf("\n# start %s<%s> Ecut= %g %s nPW=%d (est. %d:%02d minutes)\n",
                 __func__, complex_name<complex_t>(), ecut*eV, _eV, nB, int(minutes), int(seconds));
          fflush(stdout);
      } // echo
#endif // DEVEL

      // now we could sort the plane waves by their g2 length, optional
      bool constexpr sort_PWs = true;
      if (sort_PWs) {
          auto compare_lambda = [](PlaneWave const & lhs, PlaneWave const & rhs) { return lhs.g2 < rhs.g2; };
          std::sort(pw_basis.begin(), pw_basis.end(), compare_lambda);
#ifdef DEVEL
          if (echo > 8) { // show in ascending order
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
#endif // DEVEL
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
          assert( hs_PAW[ka].stride() >= nSHO );
          assert( hs_PAW[ka].dim1()   == nSHO );
          offset[ka + 1] = offset[ka] + nSHO;
      } // ka
      int const nC = offset[natoms_PAW]; // total number of PAW coefficients

      // PAW contributions to H_{ij} = P_{ik} h_{kl} P^*_{jl} = Ph_{il} P^*_{jl}
      //                  and S_{ij} = P_{ik} s_{kl} P^*_{jl} = Ps_{il} P^*_{jl}

      view2D<complex_t> P_jl(nB, nC, 0.0); //  <k+G_i|\tilde p_{la lb}>
#ifdef DEVEL
      std::vector<double>   P2_l(nC, 0.0); // |<k+G_i|\tilde p_{la lb}>|^2
#endif // DEVEL
      view3D<complex_t> Psh_il(2, nB, nC, 0.0); // atom-centered PAW matrices multiplied to P_jl
      for(int jB = 0; jB < nB; ++jB) {
          auto const & pw = pw_basis[jB];
          double const gv[3] = {(pw.x + kpoint[0])*reci[0][0],
                                (pw.y + kpoint[1])*reci[1][1],
                                (pw.z + kpoint[2])*reci[2][2]}; // ToDo: diagonal reci-matrix assumed here
          for(int ka = 0; ka < natoms_PAW; ++ka) {
              int const nSHO = sho_tools::nSHO(numax_PAW[ka]);

              // phase factor related to the atomic position
              double const arg = -(  (grid_offset[0] + xyzZ_PAW(ka,0))*gv[0] 
                                   + (grid_offset[1] + xyzZ_PAW(ka,1))*gv[1] 
                                   + (grid_offset[2] + xyzZ_PAW(ka,2))*gv[2] );
//            printf("\n# effective position"); for(int d = 0; d < 3; ++d) printf(" %g", grid_offset[d] + xyzZ_PAW(ka,d));
              std::complex<double> const phase(std::cos(arg), std::sin(arg));

              {   std::vector<double> pzyx(nSHO);
                  Hermite_Gauss_projectors(pzyx.data(), numax_PAW[ka], sigma_PAW[ka], gv);
                  auto const fgf = Fourier_Gauss_factor_3D();

                  for(int lb = 0; lb < nSHO; ++lb) {
                      int const lC = offset[ka] + lb;
                      P_jl(jB,lC) = complex_t(phase * pzyx[lb] * fgf * norm_factor);
#ifdef DEVEL
                      P2_l[lC] += pow2(pzyx[lb]);
#endif // DEVEL
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
//                    printf("# non-local matrices atom=%i i=%i j=%i ovl= %g hmt= %g %s\n",
//                            ka, lb, kb, hs_PAW[ka](1,kb,lb), hs_PAW[ka](0,kb,lb)*eV, _eV);
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
              printf_vector(" %g", &P2_l[offset[la]], nSHO);
              fflush(stdout);
          } // la
      } // echo
#endif // DEVEL

      // PAW projection matrix ready

      // allocate bulk memory for overlap and Hamiltonian
      int const nBa = align<4>(nB); // memory aligned main matrix stride
      view3D<complex_t> SHm(2, nB, nBa, complex_t(0)); // get memory for 0:Overlap S and 1:Hamiltonian matrix H

#ifdef DEVEL
      if (echo > 9) printf("# assume dimensions of Vcoeff(%d, %d, %d, %d)\n", 2, Vcoeff.dim2(), Vcoeff.dim1(), Vcoeff.stride());
      assert( nG[0] <= Vcoeff.stride() );
      assert( nG[1] <= Vcoeff.dim1() );
      assert( nG[2] <= Vcoeff.dim2() );

      double const scale_k = control::get("hamiltonian.scale.kinetic", 1.);
      double const scale_p = control::get("hamiltonian.scale.potential", 1.);
      if (1 != scale_k) warn("kinetic energy is scaled by %g", scale_k);
      if (1 != scale_p) warn("local potential is scaled by %g", scale_p);
      real_t const scale_h = control::get("hamiltonian.scale.nonlocal.h", 1.);
      real_t const scale_s = control::get("hamiltonian.scale.nonlocal.s", 1.);
      if (1 != scale_h || 1 != scale_s) warn("scale PAW contributions to H and S by %g and %g, respectively", scale_h, scale_s);
#else  // DEVEL
      real_t constexpr scale_h = 1, scale_s = 1;
      double constexpr scale_k = 1, scale_p = 1;
#endif // DEVEL
      double const kinetic = 0.5 * scale_k; // prefactor of kinetic energy in Hartree atomic units
      real_t const localpot = scale_p / (nG[0]*nG[1]*nG[2]);

      for(int iB = 0; iB < nB; ++iB) {
          auto const & i = pw_basis[iB];

          SHm(S,iB,iB) += 1; // unity overlap matrix, plane waves are orthogonal
          SHm(H,iB,iB) += kinetic*i.g2; // kinetic energy contribution

          for(int jB = 0; jB < nB; ++jB) {
              auto const & j = pw_basis[jB];

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

      auto const nB_auto = int(control::get("pw_hamiltonian.dense.solver.below", 999.));
      char const solver = *control::get("pw_hamiltonian.solver", "auto") | 32; // expect one of {auto, both, direct, iterative}
      bool const run_solver[2] = {('i' == solver) || ('b' == solver) || (('a' == solver) && (nB >  nB_auto)),
                                  ('d' == solver) || ('b' == solver) || (('a' == solver) && (nB <= nB_auto))};
      status_t solver_stat(0);
      std::vector<double> eigenenergies(nbands, -9e9);
      if (run_solver[0]) solver_stat += iterative_solve(SHm, x_axis, echo, nbands, direct_ratio);
      if (run_solver[1]) solver_stat += dense_solver::solve(SHm, x_axis, echo, nbands, eigenenergies.data());
      // dense solver must runs second in case of "both" since it modifies the memory locations of SHm

      int const gen_density = control::get("pw_hamiltonian.density", 0.);
      if (gen_density) {

          double const n_electrons = control::get("valence.electrons", 0.);
          double const kT = control::get("electronic.temperature", 9.765625e-4);
          // determine the occupation numbers
          view2D<double> occupation(2, nbands, 0.0);
          double Fermi_energy{-9e99};
          auto const stat = fermi_distribution::Fermi_level(Fermi_energy, occupation[0],
              eigenenergies.data(), nbands, kT, n_electrons, 2, echo + 9, occupation[1]);
          solver_stat += stat;

          auto const nG_all = size_t(nG[2])*size_t(nG[1])*size_t(nG[0]);
          view3D<double> rho(nG[2], nG[1], nG[0], 0.0);
          std::vector<int> ncoeff(natoms_PAW, 0), ncoeff2(natoms_PAW, 0);
          for(int ka = 0; ka < natoms_PAW; ++ka) {
              ncoeff[ka] = sho_tools::nSHO(numax_PAW[ka]);
              ncoeff2[ka] = pow2(ncoeff[ka]);
          } // ka
          data_list<double> atom_rho(ncoeff2);
          std::complex<double> const zero(0);
          std::vector<std::complex<double>> atom_coeff(nC);
          view3D<std::complex<double>> psi_G(nG[2], nG[1], nG[0]); // Fourier space array
          view3D<std::complex<double>> psi_r(nG[2], nG[1], nG[0]); // real space array

          double const kpoint_weight = 1; // depends on ikpoint, ToDo
          for(int iband = 0; iband < nbands; ++iband) {
              double const band_occupation = 2*occupation(0,iband); // factor 2 for spin-paired
              double const weight_nk = band_occupation * kpoint_weight;
              if (weight_nk > 1e-16) {
                  if (echo > 6) { printf("# band #%i, energy= %g %s, occupation= %g, response= %g, weight= %g\n",
                      iband, eigenenergies[iband]*eV,_eV, occupation(0,iband), occupation(1,iband), weight_nk); fflush(stdout); }
                  // fill psi_G
                  set(psi_G.data(), nG_all, zero);
                  set(atom_coeff.data(), nC, zero);

                  for(int iB = 0; iB < nB; ++iB) {
                      auto const & i = pw_basis[iB];
                      // the eigenvectors are stored in the memory location of H if a direct solver method has been applied
                      int const iGx = (i.x + nG[0])%nG[0], iGy = (i.y + nG[1])%nG[1], iGz = (i.z + nG[2])%nG[2];
//                    if (echo > 66) { printf("# in Fourier space psi(x=%d,y=%d,z=%d) = eigenvector(band=%i,pw=%i)\n", iGx,iGy,iGz, iband,iB); fflush(stdout); }
                      std::complex<double> const eigenvector_coeff = SHm(H,iband,iB);
                      psi_G(iGz,iGy,iGx) = eigenvector_coeff;
                      for(int iC = 0; iC < nC; ++iC) {
                          atom_coeff[iC] += std::complex<double>(P_jl(iB,iC)) * eigenvector_coeff; 
                      } // iC
                  } // iB
                  if (echo > 6) { printf("# Fourier space array for band #%i has been filled\n", iband); fflush(stdout); }

                  auto const fft_stat = fourier_transform::fft(psi_r.data(), psi_G.data(), nG, false, echo);
                  if (echo > 1) printf("# Fourier transform for band #%i returned status= %i\n", iband, int(fft_stat));

                  // accumulate real-space density rho(r) = |psi(r)|^2
                  for(int iz = 0; iz < nG[2]; ++iz) {
                  for(int iy = 0; iy < nG[1]; ++iy) {
                  for(int ix = 0; ix < nG[0]; ++ix) {
                      rho(iz,iy,ix) += weight_nk * std::norm(psi_r(iz,iy,ix));
                  }}} // iz iy ix

                  // accumulate atomic density matrices
                  for(int ka = 0; ka < natoms_PAW; ++ka) {
                      auto const density_matrix = atom_rho[ka];
                      int const nc = ncoeff[ka], stride = nc;
                      for(int ic = 0; ic < nc; ++ic) {
                          auto const c_i = conjugate(atom_coeff[ic + offset[ka]]);
                          if (echo > 1) printf("# band #%i atom #%i c[%i]= %g %g |c|= %g\n", 
                                        iband, ka, ic, std::real(c_i), std::imag(c_i), std::abs(c_i));
                          for(int jc = 0; jc < nc; ++jc) { 
                              auto const c_j = atom_coeff[jc + offset[ka]];
                              density_matrix[ic*stride + jc] += weight_nk * std::real(c_i * c_j); // Caution: imaginary part is dropped!
                          } // jc
                      } // ic
                  } // ka

              } // weight nonzero
          } // iband

          double rho_sum{0};
          for(size_t izyx = 0; izyx < nG_all; ++izyx) {
              rho_sum += rho(0,0)[izyx]; // avoid index checking
          } // izyx
          if (echo > 0) printf("# pw_hamiltonian found a density of %g electrons, average %g\n", rho_sum, rho_sum/nG_all);

#ifdef DEVEL
          if (echo > 6) { // show atomic density matrices
              for(int ia = 0; ia < natoms_PAW; ++ia) {
                  int const nc = ncoeff[ia];
                  double max_rho{1e-12}; for(int ij = 0; ij < nc*nc; ++ij) max_rho = std::max(max_rho, atom_rho[ia][ij]);
                  printf("\n# show %d x %d density matrix for atom #%i in %s-order, normalized to maximum %.6e\n", 
                      nc, nc, ia, sho_tools::SHO_order2string(sho_tools::order_zyx).c_str(), max_rho);
                  for(int i = 0; i < nc; ++i) {
                      printf("# %3i\t", i);
                      printf_vector(" %9.6f", &atom_rho[ia][i*nc], nc, "\n", 1./max_rho);
                  } // i
                  printf("\n");
              } // ia
          } // echo
#endif // DEVEL
          
      } // gen_density
      
      return solver_stat;
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
      
      for(int ia = 0; ia < natoms_PAW; ++ia) {
          if (echo > 0) {
              printf("# atom#%i \tZ=%g \tposition %12.6f%12.6f%12.6f %s\n",
                  ia, xyzZ(ia,3), xyzZ(ia,0)*Ang, xyzZ(ia,1)*Ang, xyzZ(ia,2)*Ang, _Ang);
          } // echo
      } // ia

      double const grid_offset[3] = {0.5*(g[0] - 1)*g.h[0], // position of the cell center w.r.t. the position of grid point(0,0,0)
                                     0.5*(g[1] - 1)*g.h[1],
                                     0.5*(g[2] - 1)*g.h[2]};
      double const cell_matrix[3][4] = {{g[0]*g.h[0], 0, 0,   0},
                                        {0, g[1]*g.h[1], 0,   0},
                                        {0, 0, g[2]*g.h[2],   0}};
      double reci_matrix[3][4];
      auto const cell_volume = simple_math::invert(3, reci_matrix[0], 4, cell_matrix[0], 4);
      scale(reci_matrix[0], 3*4, 2*constants::pi); // scale by 2\pi
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
      stat += (std::abs(cell_volume) < 1e-9);
      double const svol = 1./std::sqrt(cell_volume); // normalization factor for plane waves
      if (echo > 1) printf("# normalization factor for plane waves is %g a.u.\n", svol);

      char const *_ecut_u{nullptr};
      auto const ecut_u = unit_system::energy_unit(control::get("pw_hamiltonian.cutoff.energy.unit", "Ha"), &_ecut_u);
      auto const ecut = control::get("pw_hamiltonian.cutoff.energy", 11.)/ecut_u; // 11 Ha =~= 300 eV cutoff energy
      if (echo > 1) printf("# pw_hamiltonian.cutoff.energy=%.3f %s corresponds to %.3f^2 Ry or %.2f %s\n", 
                              ecut*ecut_u, _ecut_u, std::sqrt(2*ecut), ecut*eV,_eV);

      int const nG[3] = {g[0], g[1], g[2]}; // same as the numbers of real-space grid points
      view4D<double> Vcoeffs(2, nG[2], nG[1], nG[0], 0.0);
      view3D<double> vtot_Im(   nG[2], nG[1], nG[0], 0.0);

      auto const fft_stat = fourier_transform::fft(Vcoeffs(0,0,0), Vcoeffs(1,0,0), vtot, vtot_Im(0,0), nG, true, echo);
      if (0 == fft_stat) {
          if (echo > 5) printf("# used FFT on %d x %d x %d to transform local potential into space of plane-wave differences\n", nG[0], nG[1], nG[2]);
      } else
      { // scope: Fourier transform the local potential "by hand"
          if (echo > 5) { printf("# FFT failed, transform local potential manually\n"); fflush(stdout); }
          SimpleTimer timer(__FILE__, __LINE__, "manual Fourier transform (not FFT)");
          double const two_pi = 2*constants::pi;
          double const tpi_g[] = {two_pi/g[0], two_pi/g[1], two_pi/g[2]};
          for(int iGz = 0; iGz < nG[2]; ++iGz) {  auto const Gz = iGz*tpi_g[2];
          for(int iGy = 0; iGy < nG[1]; ++iGy) {  auto const Gy = iGy*tpi_g[1];
          for(int iGx = 0; iGx < nG[0]; ++iGx) {  auto const Gx = iGx*tpi_g[0];
              double re{0}, im{0};
              for(int iz = 0; iz < g[2]; ++iz) {
              for(int iy = 0; iy < g[1]; ++iy) {
              for(int ix = 0; ix < g[0]; ++ix) {
                  double const arg = -(Gz*iz + Gy*iy + Gx*ix); // argument must be negative to match the FFT
                  double const V = vtot[(iz*g[1] + iy)*g[0] + ix];
                  re += V * std::cos(arg);
                  im += V * std::sin(arg);
              }}} // ix iy iz
              Vcoeffs(0,iGz,iGy,iGx) = re; // store
              Vcoeffs(1,iGz,iGy,iGx) = im; // store
          }}} // iGx iGy iGz

      } // scope
      if (echo > 6) {
          printf("# 000-Fourier coefficient of the potential is %g %g %s\n",
                            Vcoeffs(0,0,0,0)*eV, Vcoeffs(1,0,0,0)*eV, _eV); 
      } // echo
      
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

      double const nbands_per_atom = control::get("bands.per.atom", 10.);
      int const nbands = int(nbands_per_atom*natoms_PAW);

      // all preparations done, start k-point loop

      auto const floating_point_bits = int(control::get("hamiltonian.floating.point.bits", 64.)); // double by default
      auto const nkpoints = int(control::get("hamiltonian.test.kpoints", 17.));
      float const iterative_direct_ratio = control::get("pw_hamiltonian.iterative.solver.ratio", 2.);
      simple_stats::Stats<double> nPW_stats, tPW_stats;
      for(int ikp = 0; ikp < nkpoints; ++ikp) {
          double const kpoint[3] = {ikp*.5/std::max(1., nkpoints - 1.), 0, 0};
          char x_axis[96]; std::snprintf(x_axis, 95, "# %.6f spectrum ", kpoint[0]);
          SimpleTimer timer(__FILE__, __LINE__, x_axis, 0);

          bool can_be_real{false}; // real only with inversion symmetry
          if (can_be_real) {
              error("PW only implemented with complex");
          } else { // replace this "else" by "if(true)" to test if real and complex version agree
              int nPWs{0};
              if (32 == floating_point_bits) {
                  stat += solve_k<std::complex<float>>(ecut, reci_matrix, Vcoeffs, nG, svol,
                          natoms_PAW, xyzZ, grid_offset, numax_PAW.data(), sigma_PAW.data(), hs_PAW.data(),
                          kpoint, x_axis, nPWs, echo,
                          nbands, iterative_direct_ratio);
              } else {
                  stat += solve_k<std::complex<double>>(ecut, reci_matrix, Vcoeffs, nG, svol,
                          natoms_PAW, xyzZ, grid_offset, numax_PAW.data(), sigma_PAW.data(), hs_PAW.data(),
                          kpoint, x_axis, nPWs, echo,
                          nbands, iterative_direct_ratio);
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
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_Hamiltonian(int const echo=5) {
      status_t stat(0);

      auto const vtotfile = control::get("sho_hamiltonian.test.vtot.filename", "vtot.dat"); // vtot.dat can be created by potential_generator.
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
      printf_vector(" %g", p2.data(), nSHO, "\n", d3g);

      return stat;
  } // test_Hermite_Gauss_normalization

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_Hamiltonian(echo);
      stat += test_Hermite_Gauss_normalization(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace pw_hamiltonian
