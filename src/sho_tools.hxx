#pragma once

#include <cstdio> // printf
#include <cstdint> // int64_t

typedef int status_t;

#include "inline_math.hxx" // pow2, pow3

namespace sho_tools {

  typedef enum : int64_t { // different index orderings
      order_zyx     = 0x78797a,   // "zyx"    Cartesian order best for triple loop,         depends on numax
      order_Ezyx    = 0x78797a45, // "Ezyx"   energy-ordered Cartesian,                 independent of numax
      order_lmn     = 0x6e6d6c,   // "lmn"    Radial order best for Gaunt treatment,        depends on numax
      order_lnm     = 0x6d6e6c,   // "lnm"    Radial order best for radial basis functions, depends on numax
      order_nlm     = 0x6d6c6e,   // "nlm"    Radial order best for spherical harmonics,    depends on numax
      order_Elnm    = 0x6d6e6c45, // "Elnm"   energy-ordered Radial,                    independent of numax
      order_ln      = 0x6e6c,     // "ln"     ell-ordered radial emm-degenerate,            depends on numax
      order_Enl     = 0x6c6e45,   // "Enl"    energy-ordered Radial emm-degenerate,     independent of numax
      order_unknown = 0x3f3f3f3f  // "????"   error flag
//    order_nl      = 0x6c6e,     // "nl"     enn-ordered radial emm-degenerate,            depends on numax, not implemented
  } SHO_order_t;
  
  
  inline char const * SHO_order2string(SHO_order_t const *order) { return (char const *)order; }
  inline char const * SHO_order2string(SHO_order_t const order) { auto o = order; return SHO_order2string(&o); }

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
  inline int constexpr nSHO_radial(int const numax)    { return (numax*(numax + 4) + 4)/4; } // number of different radial eigenstates
  inline constexpr int num_ln_indices(int const numax) { return (numax*(numax + 4) + 4)/4; }
  
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
  int lm_index(int const ell, int const emm) { 
      return ell*ell + ell + emm; } // usual spherical harmonics index
      
  inline constexpr 
  int nlm_index(int const numax, int const nrn, int const ell, int const emm) {
      return lm_index(ell, emm) + (nrn*(nrn*(nrn*4 - 6*(numax + 2))
                                        + 12*numax + 3*numax*numax + 11))/3;
  } // nlm_index
  template<int numax> inline constexpr
  int nlm_index(int const nrn, int const ell, int const emm) {
      return nlm_index(numax, nrn, ell, emm); }

  inline constexpr
  int lmn_index(int const numax, int const ell, int const emm, int const nrn) {
      return ((3*numax + 5)*ell + 3*(1 + numax)*ell*ell - 2*ell*ell*ell)/6 // found through fitting emm=0, nrn=0 values
            + emm*(1 + (numax - ell)/2) // nrn_max = (numax - ell)/2
            + nrn; } // linear contribution

      
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
  int Ezyx_index(int const nx, int const ny, int const nz) { 
      return ((nx+ny+nz)*(1 + nx+ny+nz)*(2 + nx+ny+nz))/6 // contribution from previous shells
               + nx + (nz*((2+ nx+ny+nz )*2-(nz + 1)))/2; // contribution from previous row and x
  } // Ezyx_index

  // energy-ordered radial emm-degenerate index
  inline constexpr
  int Enl_index(int const ell, int const nrn) {
      return (pow2(ell + 2*nrn + 1) + 2*ell)/4; } // (ell + 2*nrn)=nu, use ((nu + 1)^2)/4 as offset and add ell/2

  // energy-ordered radial 3D index
  inline constexpr
  int Elnm_index(int const ell, int const nrn, int const emm) {
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

  template<typename int_t> inline
  status_t construct_index_table(int_t energy_ordered[], int const numax, 
                    SHO_order_t const order, int_t *inverse=nullptr, int const echo=0) {
      // construct a table of energy ordered indices
      // this is needed to address e.g. the block-diagonal SHO-transformation operator
      if (echo > 1) printf("# construct_index_table for <numax=%d> order=%s\n", numax, SHO_order2string(&order));
      if (echo > 3) printf("# ");
      int ii = 0;
      switch (order) {
        
        case order_zyx:
          for(int z = 0; z <= numax; ++z) {
              for(int y = 0; y <= numax - z; ++y) {
                  for(int x = 0; x <= numax - z - y; ++x) {
                      int const eo = Ezyx_index(x, y, z);
                      energy_ordered[ii] = eo;
                      if (echo > 4) printf(" %d", eo);
                      if (inverse) inverse[eo] = ii;
                      ++ii;
          }}} // x y z
          assert(nSHO(numax) == ii);
        break;
        
        case order_lmn:
          for(int l = 0; l <= numax; ++l) {
              for(int m = -l; m <= l; ++m) {
                  for(int n = 0; n <= (numax - l)/2; ++n) {
                      int const eo = Elnm_index(l, n, m);
                      energy_ordered[ii] = eo;
                      if (echo > 4) printf(" %d", eo);
                      if (inverse) inverse[eo] = ii;
                      ++ii;
          }}} // l m n
          assert(nSHO(numax) == ii);
        break;

        case order_lnm:
          for(int l = 0; l <= numax; ++l) {
              for(int n = 0; n <= (numax - l)/2; ++n) {
                  for(int m = -l; m <= l; ++m) {
                      int const eo = Elnm_index(l, n, m);
                      energy_ordered[ii] = eo;
                      if (echo > 4) printf(" %d", eo);
                      if (inverse) inverse[eo] = ii;
                      ++ii;
          }}} // l n m
          assert(nSHO(numax) == ii);
        break;

        case order_nlm:
          for(int n = 0; n <= numax/2; ++n) {
              for(int l = 0; l <= numax - 2*n; ++l) {
                  for(int m = -l; m <= l; ++m) {
                      int const eo = Elnm_index(l, n, m);
                      energy_ordered[ii] = eo;
                      if (echo > 4) printf(" %d", eo);
                      if (inverse) inverse[eo] = ii;
                      ++ii;
          }}} // n l m
          assert(nSHO(numax) == ii);
        break;

        case order_ln: // emm-degenerate
          for(int l = 0; l <= numax; ++l) {
              for(int n = 0; n <= (numax - l)/2; ++n) {
                  int const eo = Enl_index(l, n);
                  energy_ordered[ii] = eo;
                  if (echo > 3) printf(" %d", eo);
                  if (inverse) inverse[eo] = ii;
                  ++ii;
          }} // l n
          assert(nSHO_radial(numax) == ii);
        break;

        case order_Enl: // already energy-ordered emm-degenerate
          for(int ii = 0; ii < nSHO_radial(numax); ++ii) {
              energy_ordered[ii] = ii;
              if (inverse) inverse[ii] = ii;
          } // ii
          if (echo > 3) printf(" <unity> ");
        break;

        case order_Ezyx: // already energy-ordered
        case order_Elnm: // already energy-ordered
          for(int ii = 0; ii < nSHO(numax); ++ii) {
              energy_ordered[ii] = ii;
              if (inverse) inverse[ii] = ii;
          } // ii
          if (echo > 3) printf(" <unity> ");
        break;

        default:
            if (echo > 0) printf("# no such case implemented: order=%s\n", SHO_order2string(&order));
            return order; // error
      } // switch order
      if (echo > 3) printf("\n\n");
      return 0; // success if 0
  } // construct_index_table

  status_t all_tests();
  
} // namespace sho_tools
