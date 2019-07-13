#pragma once

typedef int status_t;

namespace sho_tools {
  

  // number of all 3D SHO states up to numax >= 0
  inline constexpr int nSHO(int const numax) { return ((1 + numax)*(2 + numax)*(3 + numax))/6; }

  // number of all 2D CHO (circular harmonic oscillator) states up to numax
  inline constexpr int n2HO(int const numax) { return ((1 + numax)*(2 + numax))/2; }

  // number of all 1D HO (harmonic oscillator) states up to numax
  inline constexpr int n1HO(int const numax) { return (1 + numax); }
  
  // number of all energy-degenerate 3D SHO states in the subspace nu >= 0
  inline constexpr int ndeg(int const nu) { return n2HO(nu); }

  // =============================== Radial indices ===============================
  
  // emm-degenerate
  inline constexpr int ln_index(int const numax, int const ell, int const nrn)
      { return nrn + (1 + ell*( (2*numax + 4) - ell))/4; } // ln_index
  template<int numax> inline constexpr
  int ln_index(int const ell, int const nrn) 
      { return ln_index(numax, ell, nrn); }

  // emm-resolved
  inline constexpr 
  int lnm_index(int const numax, int const ell, int const nrn, int const emm) {
      return (6*(ell + emm) +  // within each block of size (2*ell + 1)
              6*nrn*(2*ell + 1) + // for each block of size (2*ell + 1)
               (ell*( (3*((ell + numax)%2) - 1) // linear term in ell
              + ell*((3*numax + 6) - 2*ell))))/6; // quadratic and cubic terms
  } // lnm_index
  template<int numax> inline constexpr
  int lnm_index(int const ell, int const nrn, int const emm) 
      { return lnm_index(numax, ell, nrn, emm); }

  inline constexpr
  int lmn_index(int const numax, int const ell, int const emm, int const nrn) {
      return ((3*numax + 5)*ell + 3*(1 + numax)*ell*ell - 2*ell*ell*ell)/6 // found through fitting emm=0, nrn=0 values
            + emm*(1 + (numax - ell)/2) // nrn_max = (numax - ell)/2
            + nrn; } // linear contribution

  inline constexpr
  int lm_index(int const ell, int const emm) { return ell*ell + ell + emm; } // usual spherical harmonics index
      
  // =============================== Cartesian indices ===============================

  // Cartesian 3D index
  inline constexpr int zyx_index(int const numax, int const nx, int const ny, int const nz) { 
      return (nz*nz*nz - 3*(2 + numax)*nz*nz + (3*(2 + numax)*(2 + numax) - 1)*nz // contribution for nx=0, ny=0
              + ny*(2 + numax - nz)*6 - (ny*(1 + ny))*3  +  6*nx)/6;  // contribution from previous y and x
  } // zyx_index
  template<int numax, typename int_t>
  inline constexpr int zyx_index(int_t const nx, int_t const ny, int_t const nz)
      { return zyx_index(numax, nx, ny, nz); }
  template<typename int_t>
  inline constexpr int zyx_index(int const numax, int_t const nzyx[3])
      { return zyx_index(numax, nzyx[0], nzyx[1], nzyx[2]); }
  template<int numax, typename int_t>
  inline constexpr int zyx_index(int_t const nzyx[3]) 
      { return zyx_index<numax>(nzyx[0], nzyx[1], nzyx[2]); }

  // =============================== Energy-ordered indices ===============================

  // energy-ordered Cartesian 3D index
  inline constexpr int nzyx_index(int const nx, int const ny, int const nz) { 
      return ((nx+ny+nz)*(1 + nx+ny+nz)*(2 + nx+ny+nz))/6 // contribution from previous shells
               + nx + (nz*((2+ nx+ny+nz )*2-(nz + 1)))/2; // contribution from previous row and x
  } // nzyx_index

  // energy-ordered radial emm-degenerate index
  inline constexpr int nln_index (int const ell, int const nrn) {
      return ((ell + 2*nrn + 1)*(ell + 2*nrn + 1) + 2*ell)/4; } // (ell + 2*nrn)=nu, use ((nu + 1)^2)/4 as offset and add ell/2
  // energy-ordered radial 3D index
  inline constexpr int nlnm_index(int const ell, int const nrn, int const emm) {
      return ((ell + 2*nrn)*(ell + 2*nrn + 1)*(ell + 2*nrn + 2) // energy shell offset (nu*(nu+1)*(nu+2))/6
              + 3*ell*(ell - 1) // previous ells (ell*(ell - 1))/2
              + 6*(emm + ell))/6; } // linear emm-contribution

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  inline status_t test_radial_indices(int echo=4) {
      status_t nerrors = 0;
      for(int numax = 0; numax <= 9; ++numax) {
          if (echo > 6) printf("\n# %s: numax == %d\n", __func__, numax);
          int lnm = 0, ln = 0, lm = 0, lmn = 0;
          for(int ell = 0; ell <= numax; ++ell) {
              for(int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
                  int const k = ln_index(numax, ell, nrn);
                  if ((echo > 7) && (k != ln)) printf("# ln_index<%d>(ell=%d, nrn=%d) == %d %d diff=%d\n", numax, ell, nrn, ln, k, k - ln);
                  assert(k == ln);
                  nerrors += (k != ln);
                  ++ln;
                  for(int emm = -ell; emm <= ell; ++emm) {
                      int const k = lnm_index(numax, ell, nrn, emm);
                      if (echo > 8) printf("# lnm_index<%d>(ell=%d, nrn=%d, emm=%d) == %d %d diff=%d\n", numax, ell, nrn, emm, lnm, k, k - lnm);
                      assert(k == lnm);
                      nerrors += (k != lnm);
                      ++lnm;
                  } // emm
              } // nrn
              for(int emm = -ell; emm <= ell; ++emm) {
                  int const k = lm_index(ell, emm);
                  if (echo > 7) printf("# lm_index(ell=%d, emm=%d) == %d %d diff=%d\n", ell, emm, lm, k, k - lm);
                  assert(k == lm);
                  nerrors += (k != lm);
                  ++lm;
                  for(int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
                      int const k = lmn_index(numax, ell, emm, nrn);
                      if (echo > 8) printf("# lmn_index<%d>(ell=%d, emm=%d, nrn=%d) == %d %d diff=%d\n", numax, ell, emm, nrn, lmn, k, k - lmn);
                      assert(k == lmn);
                      nerrors += (k != lmn);
                      ++lmn;
                  } // nrn
              } // emm
              assert((1 + ell)*(1 + ell) == lm); // checksum

          } // ell
          if (echo > 6) printf("\n# lmn_index<%d>\n", numax);
      } // numax
      if (nerrors && echo > 1) printf("# Warning: %s found %d errors!\n", __func__, nerrors);
      return nerrors;
  } // test_radial_indices

  inline status_t test_Cartesian_indices(int echo=3) {
      status_t nerrors = 0;
      for(int numax = 0; numax <= 9; ++numax) {
          if (echo > 6) printf("\n# %s: numax == %d\n", __func__, numax);
          int zyx = 0;
          for(int nz = 0; nz <= numax; ++nz) {
              for(int ny = 0; ny <= numax - nz; ++ny) {
                  for(int nx = 0; nx <= numax - nz - ny; ++nx) {
                      int const k = zyx_index(numax, nx, ny, nz);
                      if (echo > 8) printf("# zyx_index<%d>(nx=%d, ny=%d, nz=%d) == %d %d diff=%d\n", numax, nx, ny, nz, zyx, k, k - zyx);
                      assert(k == zyx);
                      nerrors += (k != zyx);
                      ++zyx;
                  } // nx
              } // ny
          } // nz
      } // numax
      if (nerrors && echo > 1) printf("# Warning: %s found %d errors!\n", __func__, nerrors);
      return nerrors;
  } // test_Cartesian_indices

  inline status_t test_energy_ordered_indices(int echo=4) {
      status_t nerrors = 0;
      int const numax = 9;
      if (echo > 6) printf("\n# %s: numax == %d\n", __func__, numax);
      int nzyx = 0, nln = 0, nlnm = 0;
      for(int nu = 0; nu <= numax; ++nu) { // shell index
        
          // Cartesian energy ordered index
          if (echo > 7) printf("\n# nzyx_index<nu=%d>\n", nu);
          int xyz = 0;
          for(int nz = 0; nz <= nu; ++nz) { // z can also run backwards so we start from (0,0,nu) and proceed with (1,0,nu-1)
              for(int nx = 0; nx <= nu - nz; ++nx) {
                  int const ny = nu - nz - nx;
                  int const k = nzyx_index(nx, ny, nz);
                  if ((echo > 6) && (k != nzyx))
                      printf("# nzyx_index<nu=%d>(nx=%d, ny=%d, nz=%d) == %d %d diff=%d  xyz=%d %d\n", 
                             nu, nx, ny, nz, nzyx, k, k - nzyx, xyz,  nx + (nz*((2+nu)*2-(nz + 1)))/2 );
                  assert(k == nzyx);
                  nerrors += (k != nzyx);
                  ++nzyx;
                  ++xyz;
              } // ny
          } // nz
          assert(nSHO(nu) == nzyx); // checksum
          
          // radial energy ordered indices
          for(int ell = nu%2; ell <= nu; ell+=2) {
              int const nrn = (nu - ell)/2;
              int const k = nln_index(ell, nrn);
              if (echo > 9) printf("# nln_index<nu=%d>(ell=%d, nrn=%d) == %d %d\n", nu, ell, nrn, nln, k);
              assert(k == nln);
              ++nln;
              for(int emm = -ell; emm <= ell; ++emm) {
                  int const k = nlnm_index(ell, nrn, emm);
                  if (echo > 9) printf("# nlnm_index<nu=%d>(ell=%d, nrn=%d, emm=%d) == %d\n", nu, ell, nrn, emm, nlnm);
                  assert(k == nlnm);
                  nerrors += (k != nlnm);
                  ++nlnm;
              } // emm
          } // ell
          assert(nSHO(nu) == nlnm); // checksum
          
      } // nu
      if (nerrors && echo > 1) printf("# Warning: %s found %d errors!\n", __func__, nerrors);
      return nerrors;
  } // test_energy_ordered_indices

  inline status_t all_tests() {
    auto status = 0;
    status += test_radial_indices();
    status += test_Cartesian_indices();
    status += test_energy_ordered_indices();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  
  
  
} // namespace sho_tools
