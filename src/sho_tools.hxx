#pragma once

#include <cstdio> // printf

typedef int status_t;

#include "inline_math.hxx" // pow2

namespace sho_tools {

//   union _order_u { uint32_t i; char c[4]; };

  typedef enum { // different index orderings
      order_zyx     = 0x207a7978, // " zyx" Cartesian order best for triple loop, depends on numax
      order_lmn     = 0x206c6d6e, // " lmn" Radial order best for Gaunt treatment, depends on numax
      order_lnm     = 0x206c6e6d, // " lnm" Radial order best for radial basis functions, depends on numax
      order_nlnm    = 0x6e6c6e6d, // "nlnm" energy-ordered Radial
      order_ln      = 0x20206c6e, // "  ln" Radial emm-degenerate, depends on numax
      order_nln     = 0x206e6c6e, // " nln" energy-ordered Radial emm-degenerate
      order_nzyx    = 0x6e7a7978, // "nzyx" energy-ordered Cartesian
      order_unknown = 0x3f3f3f3f  // "????" error flag
  } SHO_order_t;

  // number of all 3D SHO states up to numax >= 0
  inline constexpr int nSHO(int const numax) { return ((1 + numax)*(2 + numax)*(3 + numax))/6; }

  // number of all 2D CHO (circular harmonic oscillator) states up to numax
  inline constexpr int n2HO(int const numax) { return ((1 + numax)*(2 + numax))/2; }

  // number of all 1D HO (harmonic oscillator) states up to numax
  inline constexpr int n1HO(int const numax) { return (1 + numax); }
  
  // number of all energy-degenerate 3D SHO states in the subspace nu >= 0
  inline constexpr int ndeg(int const nu) { return n2HO(nu); }

  // =============================== Radial indices ===============================

  inline constexpr int num_ln_indices(int const numax) { return (numax*(numax + 4) + 4)/4; }
  
  // emm-degenerate
  inline constexpr int ln_index(int const numax, int const ell, int const nrn)
      { return nrn + (1 + ell*( (2*numax + 4) - ell))/4; } // ln_index
  template<int numax> inline constexpr
  int ln_index(int const ell, int const nrn) { 
      return ln_index(numax, ell, nrn); }

  // emm-resolved
  inline constexpr 
  int lnm_index(int const numax, int const ell, int const nrn, int const emm) {
      return (6*(ell + emm) +  // within each block of size (2*ell + 1)
              6*nrn*(2*ell + 1) + // for each block of size (2*ell + 1)
               (ell*( (3*((ell + numax)%2) - 1) // linear term in ell
              + ell*((3*numax + 6) - 2*ell))))/6; // quadratic and cubic terms
  } // lnm_index
  template<int numax> inline constexpr
  int lnm_index(int const ell, int const nrn, int const emm) {
      return lnm_index(numax, ell, nrn, emm); }

  inline constexpr
  int lmn_index(int const numax, int const ell, int const emm, int const nrn) {
      return ((3*numax + 5)*ell + 3*(1 + numax)*ell*ell - 2*ell*ell*ell)/6 // found through fitting emm=0, nrn=0 values
            + emm*(1 + (numax - ell)/2) // nrn_max = (numax - ell)/2
            + nrn; } // linear contribution

  inline constexpr
  int lm_index(int const ell, int const emm) { 
      return ell*ell + ell + emm; } // usual spherical harmonics index
      
  // =============================== Cartesian indices ===============================

  // Cartesian 3D index
  inline constexpr int zyx_index(int const numax, int const nx, int const ny, int const nz) { 
      return (nz*nz*nz - 3*(2 + numax)*nz*nz + (3*pow2(2 + numax) - 1)*nz // contribution for nx=0, ny=0
              + ny*(2 + numax - nz)*6 - (ny*(1 + ny))*3  +  6*nx)/6;  // contribution from previous y and x
  } // zyx_index
  template<int numax, typename int_t> inline constexpr 
  int zyx_index(int_t const nx, int_t const ny, int_t const nz) {
      return zyx_index(numax, nx, ny, nz); }
  template<typename int_t> inline constexpr 
  int zyx_index(int const numax, int_t const nzyx[3]) {
      return zyx_index(numax, nzyx[0], nzyx[1], nzyx[2]); }
  template<int numax, typename int_t> inline constexpr
  int zyx_index(int_t const nzyx[3]) {
      return zyx_index<numax>(nzyx[0], nzyx[1], nzyx[2]); }

  // =============================== Energy-ordered indices ===============================

  // energy-ordered Cartesian 3D index
  inline constexpr
  int nzyx_index(int const nx, int const ny, int const nz) { 
      return ((nx+ny+nz)*(1 + nx+ny+nz)*(2 + nx+ny+nz))/6 // contribution from previous shells
               + nx + (nz*((2+ nx+ny+nz )*2-(nz + 1)))/2; // contribution from previous row and x
  } // nzyx_index

  // energy-ordered radial emm-degenerate index
  inline constexpr
  int nln_index(int const ell, int const nrn) {
      return (pow2(ell + 2*nrn + 1) + 2*ell)/4; } // (ell + 2*nrn)=nu, use ((nu + 1)^2)/4 as offset and add ell/2

  // energy-ordered radial 3D index
  inline constexpr
  int nlnm_index(int const ell, int const nrn, int const emm) {
      return ((ell + 2*nrn)*(ell + 2*nrn + 1)*(ell + 2*nrn + 2) // energy shell offset (nu*(nu+1)*(nu+2))/6
              + 3*ell*(ell - 1) // previous ells (ell*(ell - 1))/2
              + 6*(emm + ell))/6; } // linear emm-contribution

  template<typename int_t> inline constexpr
  int_t get_nu(int_t const nx, int_t const ny, int_t const nz) { return nx + ny + nz; }

  template<typename int_t> inline constexpr
  int_t get_nu(int_t const ell, int_t const nrn) { return ell + 2*nrn; }

  template<typename int_t> inline
  int get_nu(int_t const energy_ordered) {
      int nu = -1; while (energy_ordered >= nSHO(nu)) { ++nu; } return nu; }

//   template<typename int_t>
//   inline status_t zyx_translation_table(int_t table[], int const numax, bool const inverse=false) {
//       int const nsho = nSHO(numax);
//       int isho = 0;
//       for(int nz = 0; nz <= numax; ++nz) {
//           for(int ny = 0; ny <= numax - nz; ++ny) {
//               for(int nx = 0; nx <= numax - nz - ny; ++nx) {
//                   int const izyx = zyx_index(numax, nx, ny, nz);
//                   assert(izyx >= 0); assert(izyx < nsho);
//                   assert(izyx == isho);
//                   int const nzyx = nzyx_index(nx, ny, nz);
//                   assert(nzyx >= 0); assert(nzyx < nsho);
//                   if (inverse) {
//                       table[izyx] = nzyx;
//                   } else {
//                       table[nzyx] = izyx;
//                   }
//                   ++isho;
//               } // nx
//           } // ny
//       } // nz
//       return 0;
//   } // get_translation_table

  template<typename int_t> inline
  status_t construct_index_table(int_t energy_ordered[], int const numax, 
                    SHO_order_t const order, int_t *inverse=nullptr, int const echo=0) {
      // construct a table of energy ordered indices
      if (echo > 1) printf("# construct_index_table for <numax=%d> order=%c%c%c%c\n", numax, order>>24, order>>16, order>>8, order);
      if (echo > 3) printf("# ");
      int ii = 0;
      switch (order) {
        case order_zyx:
          for(int z = 0; z <= numax; ++z) {
              for(int y = 0; y <= numax - z; ++y) {
                  for(int x = 0; x <= numax - z - y; ++x) {
                      int const nzyx = nzyx_index(x, y, z);
                      energy_ordered[ii] = nzyx;
                      if (echo > 4) printf(" %d", nzyx);
                      if (nullptr != inverse) inverse[nzyx] = ii;
                      ++ii;
          }}} // x y z
          assert(nSHO(numax) == ii);
        break;
        case order_lmn:
          for(int l = 0; l <= numax; ++l) {
              for(int m = -l; m <= l; ++m) {
                  for(int n = 0; n <= (numax - l)/2; ++n) {
                      int const nlnm = nlnm_index(l, n, m);
                      energy_ordered[ii] = nlnm;
                      if (echo > 4) printf(" %d", nlnm);
                      if (nullptr != inverse) inverse[nlnm] = ii;
                      ++ii;
          }}} // l m n
          assert(nSHO(numax) == ii);
        break;
        case order_ln:
          for(int l = 0; l <= numax; ++l) {
              for(int n = 0; n <= (numax - l)/2; ++n) {
                  int const nln = nln_index(l, n);
                  energy_ordered[ii] = nln;
                  if (echo > 3) printf(" %d", nln);
                  if (nullptr != inverse) inverse[nln] = ii;
                  ++ii;
          }} // l n
        break;
        default:
            if (echo > 0) printf("# no such case implemented: order=%c%c%c%c\n", order>>24, order>>16, order>>8, order);
      } // switch order
      if (echo > 3) printf("\n\n");
      return 0; // success if 0
  } // construct_index_table

  status_t all_tests();
  
} // namespace sho_tools
