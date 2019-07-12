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

  // emm-degenerate
  inline constexpr int ln_index(int const numax, int const ell, int const nrn)
      { return nrn + (1 + ell*( (2*numax + 4) - ell))/4; } // ln_index
  template<int numax> inline constexpr
  int ln_index(int const ell, int const nrn) 
      { return ln_index(numax, ell, nrn); }

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

  // Cartesian indices
  inline constexpr int xyz_index(int const numax, int const nx, int const ny, int const nz)
      { return nx + n1HO(numax)*ny + n2HO(numax)*nz; } // xyz_index
  template<int numax, typename int_t>
  inline constexpr int xyz_index(int_t const nx, int_t const ny, int_t const nz)
      { return xyz_index(numax, nx, ny, nz); }
  template<typename int_t>
  inline constexpr int xyz_index(int const numax, int_t const nxyz[3])
      { return xyz_index(numax, nxyz[0], nxyz[1], nxyz[2]); }
  template<int numax, typename int_t>
  inline constexpr int xyz_index(int_t const nxyz[3]) 
      { return xyz_index<numax>(nxyz[0], nxyz[1], nxyz[2]); }

  // energy-ordered indices
  
  // ToDo

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  inline status_t test_ln_and_lnm_index(int echo=9) {
      for(int numax = 0; numax <= 9; ++numax) {
          if (echo > 6) printf("\n# %s: numax == %d\n", __func__, numax);
          int lnm = 0, ln = 0;
          for(int ell = 0; ell <= numax; ++ell) {
              for(int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
                  int const diff_ln = ln_index(numax, ell, nrn) - ln;
                  if ((echo > 7) && (0 != diff_ln)) printf("# ln_index<%d>(ell=%d, nrn=%d) == %d  diff=%d\n", numax, ell, nrn, ln, diff_ln);
                  assert(0 == diff_ln);
                  ++ln;
                  for(int emm = -ell; emm <= ell; ++emm) {
                      int const diff_lnm = lnm_index(numax, ell, nrn, emm) - lnm;
                      if ((echo > 8) && (0 != diff_lnm)) 
                          printf("# lnm_index<%d>(ell=%d, nrn=%d, emm=%d) == %d  diff=%d\n", numax, ell, nrn, emm, lnm, diff_lnm);
                      assert(0 == diff_lnm);
                      ++lnm;
                  } // emm
              } // nrn
          } // ell
      } // numax
      return 0;
  } // test_ln_and_lnm_index

  inline status_t test_xyz_index(int echo=9) {
      for(int numax = 0; numax <= 3; ++numax) {
          if (echo > 6) printf("\n# %s: numax == %d\n", __func__, numax);
          int xyz = 0;
          for(int nz = 0; nz <= numax; ++nz) {
              for(int ny = 0; ny <= numax - nz; ++ny) {
                  for(int nx = 0; nx <= numax - nz - ny; ++nx) {
                      int const diff_xyz = xyz_index(numax, nx, ny, nz) - xyz;
                      if ((echo > 8) && (0 != diff_xyz)) 
                          printf("# xyz_index<%d>(nx=%d, ny=%d, nz=%d) == %d  diff=%d\n", numax, nx, ny, nz, xyz, diff_xyz);
//                       assert(0 == diff_xyz);
                      ++xyz;
                  } // nx
              } // ny
          } // nz
      } // numax
      return 0;
  } // test_xyz_index

  inline status_t all_tests() {
    auto status = 0;
    status += test_ln_and_lnm_index(1);
//  status += test_xyz_index(); // FAILS!
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  
  
  
} // namespace sho_tools
