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
#include "linear_algebra.hxx" // ::eigenvalues, ::generalized_eigenvalues
#include "inline_tools.hxx" // align<nbits>
#include "simple_math.hxx" // ::random<real_t>
#include "vector_math.hxx" // ::vec<N,T>
#include "hermite_polynomial.hxx" // hermite_polys
#include "simple_stats.hxx" // ::Stats

namespace pw_hamiltonian {
  // computes Hamiltonian matrix elements in for plane waves
  // including a PAW non-local contribution

  template<typename real_t>
  real_t Lorentzian(real_t const re, real_t const im) { return -im/(pow2(im) + pow2(re)); }

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

  template<typename complex_t>
  bool constexpr is_complex(complex_t const x) { return (sizeof(complex_t) > sizeof(std::real(x))); }

  template<typename real_t> real_t conjugate(real_t const x);
  template<> inline float  conjugate<float> (float  const x) { return x; };
  template<> inline double conjugate<double>(double const x) { return x; };
  
  template<typename real_t>
  std::complex<real_t> conjugate(std::complex<real_t> const x) { return std::conj(x); }

  
  template<typename complex_t, typename real_t> inline
  complex_t to_complex_t(std::complex<real_t> const x); // no generic implementation given
  template<> inline double to_complex_t(std::complex<double> const x) { return std::real(x); }
  template<> inline float  to_complex_t(std::complex<float>  const x) { return std::real(x); }
  template<> inline std::complex<double> to_complex_t(std::complex<double> const x) { return x; }
  template<> inline std::complex<float>  to_complex_t(std::complex<float>  const x) { return x; }
  
  class PlaneWave {
    public:
      double g2; // length of the vector
      int16_t x, y, z;
      PlaneWave() : g2(0), x(0), y(0), z(0) {}
      PlaneWave(int const ix, int const iy, int const iz, double const len2) : g2(len2), x(ix), y(iy), z(iz) {}
    private:
  }; // class PlaneWave


  void Hermite_Gauss_projectors(double pzyx[], int const numax, double const sigma, double const gv[3]) {

        view2D<double> Hermite_Gauss(3, sho_tools::n1HO(numax));
        for(int d = 0; d < 3; ++d) {
            double const x = gv[d]*sigma;
            hermite_polys(Hermite_Gauss[d], x, numax);
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
  
  

  template<typename complex_t, typename real_t>
  status_t solve_k(double const ecut // plane wave cutoff energy
          , double const reci[3][4] // reciprocal cell matrix
          , view4D<std::complex<double>> const & Vcoeff // <G|V|G'> = potential(iGz-jGz,iGy-jGy,iGx-jGx,0:3)
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

      status_t stat(0);

      int boxdim[3], max_PWs{1};
      for(int d = 0; d < 3; ++d) {
          boxdim[d] = std::ceil(ecut/pow2(reci[d][d]));
          max_PWs *= (boxdim[d] + 1 + boxdim[d]);
      } // d
      std::vector<PlaneWave> pw_basis(max_PWs);
      {
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
      } // generate set of plane waves
      int const nB = pw_basis.size();
      nPWs = nB; // export the number of plane waves used for statistics
#ifdef DEVEL
      {   complex_t c{1}; real_t r{1};
          if (echo > 3) printf("\n\n# start %s<%s, real_t=%s> nB=%d\n", 
                  __func__, complex_name(c), complex_name(r), nB);
      }
#endif
      // now we could sort the plane waves by their g2 length, optional

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
      view3D<complex_t> Psh_il(2, nB, nC, 0.0); // atom-centered PAW matrices multiplied to P_jl
      for(int jB = 0; jB < nB; ++jB) {
          auto const & i = pw_basis[jB];
          double const gv[3] = {(i.x + kpoint[0])*reci[0][0],
                                (i.y + kpoint[1])*reci[1][1],
                                (i.z + kpoint[2])*reci[2][2]}; // ToDo: diagonal reci assumed here
          for(int ka = 0; ka < natoms_PAW; ++ka) {
              int const nSHO = sho_tools::nSHO(numax_PAW[ka]);

              // phase factor related to the atomic position
              double const arg = xyzZ_PAW(ka,0)*gv[0] + xyzZ_PAW(ka,1)*gv[1] + xyzZ_PAW(ka,2)*gv[2];
              std::complex<double> const phase(std::cos(arg), std::sin(arg));

              {   std::vector<double> pzyx(nSHO);
                  Hermite_Gauss_projectors(pzyx.data(), numax_PAW[ka], sigma_PAW[ka], gv);

                  for(int lb = 0; lb < nSHO; ++lb) {
                      int const lC = offset[ka] + lb;
                      P_jl(jB,lC) = complex_t(phase * pzyx[lb]);
                  } // lb
              }

              // multiply P from left to hs_PAW (block diagonal --> ka == la)
              auto const iB = jB;
              for(int lb = 0; lb < nSHO; ++lb) {
                  int const iC = offset[ka] + lb; // global PAW coefficient index
                  complex_t s{0}, h{0};
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
      
      // PAW projection matrix ready

      // allocate bulk memory for overlap and Hamiltonian
      int const nBa = align<4>(nB); // memory aligned main matrix stride
      view3D<complex_t> SHm(2, nB, nBa, complex_t(0)); // get memory for 0:Overlap S and 1:Hamiltonian matrix H
      
#ifdef DEVEL
      if (echo > 3) printf("# assume dimensions of Vcoeff(%d, %d, %d, %d)\n",
                  2*nG[2]+1, Vcoeff.dim2(), Vcoeff.dim1(), Vcoeff.stride());
      assert( 2*nG[0] < Vcoeff.dim1() );
      assert( 2*nG[1] < Vcoeff.dim2() );
      // unfortunately we cannot check dim3
      
      // prefactor of kinetic energy in Hartree atomic units
      double const kinetic = control::get("pw_hamiltonian.scale.kinetic", 1.0) * 0.5;
      if (0.5 != kinetic) warn("kinetic energy prefactor is %g", kinetic);
      real_t const scale_h = control::get("pw_hamiltonian.scale.nonlocal.h", 1.0);
      real_t const scale_s = control::get("pw_hamiltonian.scale.nonlocal.s", 1.0);
      if (1 != scale_h || 1 != scale_s) warn("scale PAW contributions to H and S by %g and %g, respectively", scale_h, scale_s);
      real_t const localpot = control::get("pw_hamiltonian.scale.potential", 1.0);
      if (1 != localpot) warn("local potential prefactor is %g", localpot);
#else
      double constexpr kinetic = 0.5; // prefactor of kinetic energy in Hartree atomic units
      real_t constexpr scale_h = 1, scale_s = 1, localpot = 1;
#endif
      
      for(int iB = 0; iB < nB; ++iB) {      auto const & i = pw_basis[iB];

          SHm(S,iB,iB) += 1; // unity overlap matrix, plane waves are orthogonal
          SHm(H,iB,iB) += kinetic*i.g2; // kinetic energy contribution
          
          for(int jB = 0; jB < nB; ++jB) {  auto const & j = pw_basis[jB];

              // add the contribution of the local potential
              int const iVx = i.x - j.x, iVy = i.y - j.y, iVz = i.z - j.z;
              if ((std::abs(iVx) <= nG[0]) && (std::abs(iVy) <= nG[1]) && (std::abs(iVz) <= nG[2])) {
                  SHm(H,iB,jB) += localpot*Vcoeff(nG[2] + iVz, nG[1] + iVy, nG[0] + iVx, 0);
              }

              // PAW contributions to H_{ij} = Ph_{il} P^*_{jl}
              //                  and S_{ij} = Ps_{il} P^*_{jl}
              
              // inner part of a matrix-matrix multiplication
              { // scope: add non-local contributions
                  complex_t s{0}, h{0};
                  for(int lC = 0; lC < nC; ++lC) { // contract
                      auto const p = conjugate(P_jl(jB,lC)); // needs a conjugation, ToDo: here or above?
                      s += Psh_il(0,iB,lC) * p;
                      h += Psh_il(1,iB,lC) * p;
                  } // lC
                  SHm(S,iB,jB) += s * scale_s;
                  SHm(H,iB,jB) += h * scale_h;
              } // scope: add non-local contributions

          } // jB
      } // iB

      // from here on, the routine is almost identical to that in sho_hamiltonian.cxx, maybe unify
      
      std::vector<real_t> eigvals(nB, 0.0);
      auto const ovl_eig = int(control::get("pw_hamiltonian.test.overlap.eigvals", 0.));
      char const hermitian = *control::get("pw_hamiltonian.test.hermitian", "both") | 32; // 'n':none, 's':overlap, 'h':Hamiltonian, 'b':both
      real_t const E_imag = control::get("pw_hamiltonian.temperature", 9.765625e-4);
      double const f_dos = -1./constants::pi;
      
      
      double const energy_range[2] = {-1, 1}; // ToDo: external input
      
      { // scope: try a Green function approach
          auto const nE = int(control::get("pw_hamiltonian.test.green.function", 0.));
          if (nE > 0) {
              complex_t const minus1(-1);
              if (!is_complex(minus1)) {
                  warn("# Green functions can only be computed in complex versions"); return stat;
              } // is not complex

              double const dE = (energy_range[1] - energy_range[0])/nE;
              if (echo > 0) printf("\n## E_real (%s), DoS:\n", _eV);
              view2D<complex_t> ESmH(nB, nBa, complex_t(0)); // get memory
              for(int iE = 0; iE <= nE; ++iE) {
                  real_t const E_real = iE*dE + energy_range[0];
                  auto const E = to_complex_t<complex_t,real_t>(std::complex<real_t>(E_real, E_imag));
                  if (echo > 99) printf("# Green function for energy point (%g %s, %g %s)\n",
                                           std::real(E)*eV,_eV, std::imag(E)*Kelvin,_Kelvin);
                  int constexpr check = 1;
                  view2D<complex_t> ESmH_copy(check*nB, nBa, complex_t(0)); // get memory
                  view2D<complex_t> Sinv(nB, nBa, complex_t(0)); // get memory
                  // construct matrix to be inverted: E*S - H
                  for(int iB = 0; iB < nB; ++iB) {
                      set(Sinv[iB], nB, SHm(S,iB)); // copy S
                      set(ESmH[iB], nB, SHm(H,iB), minus1); // -H
                      add_product(ESmH[iB], nB, SHm(S,iB), E); // +E*S
                      if (check) set(ESmH_copy[iB], nB, ESmH[iB]); // copy
                  } // iB

                  // Green function G = (E*S - H)^{-1}, mind that with a non-orthogonal basis set, we have to multiply S before interpreting the trace
                  auto const stat_inv = linear_algebra::inverse(nB, ESmH.data(), ESmH.stride());
                  view2D<complex_t> Green(ESmH.data(), ESmH.stride()); // wrap
                  if (int(stat_inv)) warn("inversion failed for %s, E= %g %s", x_axis, E_real*eV,_eV);
                  stat += stat_inv;

#ifdef DEVEL
                  if (check) {
                      real_t devN(0), devT(0);
                      for(int iB = 0; iB < nB; ++iB) {
                          for(int jB = 0; jB < nB; ++jB) {
                              complex_t cN(0), cT(0);
                              for(int kB = 0; kB < nB; ++kB) {
                                  cN += Green(iB,kB) * ESmH_copy(kB,jB); 
                                  cT += ESmH_copy(iB,kB) * Green(kB,jB); 
                              } // kB
                              complex_t const diag = (iB == jB);
                              devN = std::max(devN, std::abs(cN - diag));
                              devT = std::max(devT, std::abs(cT - diag));
                              if (echo > 99) printf("# iB=%i jB=%i cN= %g %g, dev= %.1e,\tcT= %g %g, dev= %.1e\n", iB, jB,
                                  std::real(cN), std::imag(cN), std::abs(cN - diag), std::real(cT), std::imag(cT), std::abs(cT - diag));
                          } // jB
                      } // iB
                      if ((echo > 19) || ((echo > 0) && (devN + devT > 1e-7))) {
                          printf("# deviation of G * (ES-H) from unity is %.2e and %.2e transposed\n", devN, devT);
                      }
                  } // check
#endif

                  view2D<complex_t> GreenS(nB, nBa, complex_t(0)); // get memory
                  for(int iB = 0; iB < nB; ++iB) {
                      for(int jB = 0; jB < nB; ++jB) {
                          complex_t c(0);
                          for(int kB = 0; kB < nB; ++kB) {
                              c += Green(iB,kB) * SHm(S,kB,jB); 
                          } // kB
                          GreenS(iB,jB) = c;
                      } // jB
                  } // iB

                  // density of states
                  double density{0};
                  for(int iB = 0; iB < nB; ++iB) {
                      density += f_dos*std::imag(GreenS(iB,iB));
                  } // iB
                  if (echo > 0) printf("%g %g\n", std::real(E)*eV, density);
                  if (echo > 99) printf("# DoS for energy point (%g %s, %g %s) %g\n",
                      std::real(E)*eV,_eV, std::imag(E)*Kelvin,_Kelvin, density);
                
              } // iE
              if (echo > 0) printf("\n");
          } // green_function
      } // scope
      
      
      for(int s0h1 = 0; s0h1 < 2; ++s0h1) { // loop must run forward and serial
          if (echo > 0) printf("\n");
          auto const matrix_name = s0h1 ? "Hamiltonian" : "overlap";
          real_t const  u = s0h1 ?  eV :  1; // output unit conversion factor 
          auto   const _u = s0h1 ? _eV : ""; // unit symbol

#ifdef DEVEL
          // display S and H
          if (echo > 9 - s0h1) {
              printf("\n# %s matrix (%s) for Bloch phase", matrix_name, _u);
              for(int d = 0; d < 3; ++d) {
//                   printf(" %g+i*%g ", std::real(Bloch_phase[d]), std::imag(Bloch_phase[d]));
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

          // check if the matrix is symmetric/Hermitian
          if (('b' == hermitian) || ((s0h1 ? 'h' : 's') == hermitian)) {
              real_t diag{0}; // imaginary part of diagonal elements
              real_t offr{0}; // off-diagonal elements real part
              real_t offi{0}; // off-diagonal elements imaginary part
              for(int iB = 0; iB < nB; ++iB) {
                  diag = std::max(diag, std::abs(std::imag(SHm(s0h1,iB,iB))));
                  for(int jB = iB + 1; jB < nB; ++jB) { // triangular loop
                      auto const dev = SHm(s0h1,iB,jB) - conjugate(SHm(s0h1,jB,iB));
                      offr = std::max(offr, std::abs(std::real(dev)));
                      offi = std::max(offi, std::abs(std::imag(dev)));
                  } // jB
              } // iB
              if (echo > 1) printf("# %s deviates from hermitian: imag(diag)= %.1e  imag(off)= %.1e  real(off)= %.1e\n",
                                      matrix_name, diag, offi, offr);
          } // check hermiticity
#endif
          
          // diagonalize matrix
          status_t stat_eig(0);
          if (1 == s0h1) {
              stat_eig = linear_algebra::eigenvalues(eigvals.data(), nB, SHm(H,0), nBa, SHm(S,0), nBa);
          } else if (ovl_eig) {
              view2D<complex_t> S_copy(nB, nBa); // get memory
              set(S_copy.data(), nB*nBa, SHm(S,0)); // copy overlap matrix S into work array W
              stat_eig = linear_algebra::eigenvalues(eigvals.data(), nB, S_copy.data(), nBa);
          } // ovl_eig

          // show result
          if ((H == s0h1) || ovl_eig) {
              if (stat_eig) {
                  warn("diagonalizing the %s matrix failed, status= %i", matrix_name, int(stat_eig));
                  stat += stat_eig;
              } else if (nB > 0) {
                  double const lowest_eigenvalue = eigvals[0], highest_eigenvalue = eigvals[nB - 1];
                  if (echo > 2) {
                      printf("%s%s", x_axis, s0h1 ? "" : matrix_name);
                      int constexpr mB = 8; // show at most the 6 lowest + 2 highest eigenvalues
                      for(int iB = 0; iB < std::min(nB - 2, mB - 2); ++iB) {
                          printf(" %g", eigvals[iB]*u);
                      } // iB
                      if (nB > mB) printf(" ..."); // there are more eigenvalues than we display
                      printf(" %g %g %s\n", eigvals[nB - 2]*u, eigvals[nB - 1]*u, _u); // last two
                  } // echo
                  if (echo > 4) printf("# lowest and highest eigenvalue of the %s matrix are %g and %g %s, respectively\n", 
                                            matrix_name, lowest_eigenvalue*u, highest_eigenvalue*u, _u);
                  if (S == s0h1) {
                      if (lowest_eigenvalue <= 0) {
                          warn("overlap matrix has instable eigenvalues, lowest= %g", lowest_eigenvalue);
                      } else if (lowest_eigenvalue < .1) {
                          warn("overlap matrix has critical eigenvalues, lowest= %g", lowest_eigenvalue);
                      } // lowest_eigenvalue
#ifdef DEVEL
                  } else {
                      // experiment: show Lorentians where the eigenenergies are
                      auto const nE = int(control::get("pw_hamiltonian.test.green.lehmann", 0.));
                      if (nE > 0) {
                          double const dE = (energy_range[1] - energy_range[0])/nE;
                          // from the eigenvalues, we can plot the density of states via the Lehmann representation
                          if (echo > 0) printf("\n## E_real (%s), DoS (from Lehmann representation):\n", _eV);
                          for(int iE = 0; iE <= nE; ++iE) {
                              real_t const E_real = iE*dE + energy_range[0];
                              double density{0};
                              for(int iB = 0; iB < nB; ++iB) {
                                  density += f_dos*Lorentzian(E_real - eigvals[iB], E_imag);
                              } // density
                              if (echo > 0) printf("%g %g\n", E_real*eV, density);
                          } // iE
                          if (echo > 0) printf("\n");
                      } // nE > 0
#endif
                  } // S == s0h1
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

      double const cell_matrix[3][4] = {{g[0]*g.h[0], 0, 0, 0}, {0, g[1]*g.h[1], 0, 0}, {0, 0, g[2]*g.h[2], 0}};
      double reci_matrix[3][4];
      auto const cell_volume = simple_math::invert(3, reci_matrix[0], 4, cell_matrix[0], 4);
      if (echo > 0) printf("# cell volume is %g %s^3\n", cell_volume*pow3(Ang),_Ang);
      if (echo > 0) {
          auto const sqRy = 1;
          printf("# cell matrix in %s and reciprocal matrix in sqRy:\n", _Ang);
          for(int d = 0; d < 3; ++d) {
              printf("# %12.6f%12.6f%12.6f    \t%12.6f%12.6f%12.6f\n",
                    cell_matrix[d][0]*Ang,  cell_matrix[d][1]*Ang,  cell_matrix[d][2]*Ang,
                    reci_matrix[d][0]*sqRy, reci_matrix[d][1]*sqRy, reci_matrix[d][2]*sqRy);
          } // d
          printf("\n");
      } // echo
      stat += (abs(cell_volume) > 1e-15);
      scale(reci_matrix[0], 12, 2*constants::pi); // scale by 2\pi
      double const svol = 1./std::sqrt(cell_volume); // normalization factor for plane waves

      auto const ecut = control::get("pw_hamiltonian.cutoff.energy", 11.); // 11 Ha =~= 300 eV cutoff energy
      if (echo > 1) printf("# pw_hamiltonian.cutoff.energy=%.3f Ha corresponds %.3f^2 Ry or %.2f %s\n", 
                              ecut, std::sqrt(2*ecut), ecut*eV,_eV);

//       int const nG[3] = {g[0]/2 - 1, g[1]/2 - 1, g[2]/2 - 1}; // 2*nG <= g
      int const nG[3] = {2, 1, 1}; // 2*nG <= g
      view4D<std::complex<double>> Vcoeffs(2*nG[2]+1, 2*nG[1]+1, 2*nG[0]+2, 1, 0.0);
      { // scope: Fourier transform the local potential
          double const two_pi = 2*constants::pi;
          double const w8 = 1./(g[0]*g[1]*g[2]);
          double const tpi_g[] = {two_pi/g[0], two_pi/g[1], two_pi/g[2]};
          for        (int iGz = -nG[2]; iGz <= nG[2]; ++iGz) {  auto const Gz = iGz*tpi_g[2];
              for    (int iGy = -nG[1]; iGy <= nG[1]; ++iGy) {  auto const Gy = iGy*tpi_g[1];
                  for(int iGx = -nG[0]; iGx <= nG[0]; ++iGx) {  auto const Gx = iGx*tpi_g[0];
                      std::complex<double> t(0);
                      for(int iz = 0; iz < g[2]; ++iz) {
                      for(int iy = 0; iy < g[1]; ++iy) {
                      for(int ix = 0; ix < g[0]; ++ix) {
                          double const arg = Gz*iz + Gy*iy + Gx*ix;
                          t += vtot[(iz*g[1] + iy)*g[0] + ix] * std::complex<double>(std::cos(arg), std::sin(arg));
                      } // ix
                      } // iy
                      } // iz
                      Vcoeffs(nG[2]+iGz,nG[1]+iGy,nG[0]+iGx,0) = t*w8; // store
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

      
      

      // all preparations done, start k-point loop
      
      auto const floating_point_bits = int(control::get("pw_hamiltonian.floating.point.bits", 64.)); // double by default
      auto const nkpoints = int(control::get("pw_hamiltonian.test.kpoints", 17.));
      simple_stats::Stats<double> nPW_stats;
      for(int ikp = 0; ikp < nkpoints; ++ikp) {
          double const kpoint[3] = {ikp/(nkpoints - 1.), 1, 1}; //
          char x_axis[96]; std::snprintf(x_axis, 95, "# %.6f spectrum ", ikp*.5/(nkpoints - 1.));

          bool can_be_real{false}; // real only with inversion symmetry
          if (can_be_real) {
              error("PW only implemented with complex");
          } else { // replace this else by if(true) to test if real and complex version agree
              int nPWs{0};
              if (32 == floating_point_bits) {
                  error("ToDo: instanciate PW for complex<float>");
              } else {
                  stat += solve_k<std::complex<double>, double>(ecut, reci_matrix, Vcoeffs, nG, svol,
                          natoms_PAW, xyzZ, numax_PAW.data(), sigma_PAW.data(), hs_PAW.data(),
                          kpoint, x_axis, nPWs, echo);
              } // floating_point_bits
              nPW_stats.add(nPWs);
          } // !can_be_real
          if (echo > 0) fflush(stdout);
      } // ikp
      
      if (echo > 3) printf("# average number of plane waves is %.3f +/- %.3f min %g max %g\n",
                      nPW_stats.avg(), nPW_stats.var(), nPW_stats.min(), nPW_stats.max());

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

  status_t all_tests(int const echo) {
    status_t status(0);
    status += test_Hamiltonian(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace pw_hamiltonian
