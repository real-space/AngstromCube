#pragma once

#include <cstdio> // std::printf, ::snprintf
#include <cstdint> // int64_t
#include <string> // std::string

#ifndef NO_UNIT_TESTS
  #include <vector> // std::vector<T>
#endif

#include "status.hxx" // status_t
#include "inline_math.hxx" // pow2, pow3

namespace sho_tools {

  typedef enum : int64_t { // different index orderings
      order_zyx     = 0x78797a,   // "zyx"    Cartesian order best for triple loop,         depends on numax
      order_Ezyx    = 0x78797a45, // "Ezyx"   Cartesian energy-ordered,                 independent of numax
      order_lmn     = 0x6e6d6c,   // "lmn"    Radial order best for Gaunt treatment,        depends on numax
      order_lnm     = 0x6d6e6c,   // "lnm"    Radial order best for radial basis functions, depends on numax
      order_nlm     = 0x6d6c6e,   // "nlm"    Radial order best for spherical harmonics,    depends on numax
      order_Elnm    = 0x6d6e6c45, // "Elnm"   Radial energy-ordered,                    independent of numax
      order_Enl     = 0x6c6e45,   // "Enl"    Radial energy-ordered, emm-degenerate,    independent of numax
      order_ln      = 0x6e6c,     // "ln"     ell-ordered radial     emm-degenerate,        depends on numax
      order_nl      = 0x6c6e,     // "nl"     enn-ordered radial     emm-degenerate,        depends on numax
      order_unknown = 0x3f3f3f3f  // "????"   error flag
  } SHO_order_t;

  inline constexpr bool is_energy_ordered(SHO_order_t const order) {
      return (order_Ezyx == order) || (order_Elnm == order) || (order_Enl == order); }

  inline constexpr bool is_Cartesian(SHO_order_t const order) {
      return (order_Ezyx == order) || (order_zyx == order); }

  inline constexpr bool is_emm_degenerate(SHO_order_t const order) {
      return (order_Enl == order) || (order_ln == order) || (order_nl == order); }

  inline std::string SHO_order2string(SHO_order_t const order) {
      auto const o = order;
      return std::string((char const *)&o);
  } // SHO_order2string

  // number of all 3D SHO states up to numax >= 0
  inline constexpr int nSHO(int const numax) { return ((1 + numax)*(2 + numax)*(3 + numax))/6; }

  // number of all 2D CHO (circular harmonic oscillator) states up to numax
  inline constexpr int n2HO(int const numax) { return ((1 + numax)*(2 + numax))/2; }

  // number of all 1D HO (harmonic oscillator) states up to numax
  inline constexpr int n1HO(int const numax) { return  (1 + numax); }


  // =============================== Radial indices ===============================

  // emm-degenerate
  inline constexpr int nSHO_radial(int const numax) { return (numax*(numax + 4) + 4)/4; } // number of different radial eigenstates

  inline constexpr
  int ln_index(int const numax, int const ell, int const nrn) {
      return nrn + (1 + ell*( (2*numax + 4) - ell))/4; } // ln_index
  template <int numax> inline constexpr
  int ln_index(int const ell, int const nrn) {
      return ln_index(numax, ell, nrn); }

  inline constexpr
  int nn_max(int const numax, int const ell) { return (numax + 2 - ell)/2; } // size of the ell-subspace, emm_Degenerate

  inline constexpr
  int nl_index(int const numax, int const nrn, int const ell) {
      return ell + (numax + 2 - nrn)*nrn; } // nl_index
  template <int numax> inline constexpr
  int nl_index(int const nrn, int const ell) {
      return nl_index(numax, nrn, ell); }

  // emm-resolved
  inline constexpr
  int lnm_index(int const numax, int const ell, int const nrn, int const emm) {
      return (6*(ell + emm) +  // within each block of size (2*ell + 1)
              6*nrn*(2*ell + 1) + // for each block of size (2*ell + 1)
               (ell*( (3*((ell + numax)%2) - 1) // linear term in ell
              + ell*((3*numax + 6) - 2*ell))))/6; // quadratic and cubic terms
  } // lnm_index
  template <int numax> inline constexpr
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
  template <int numax> inline constexpr
  int nlm_index(int const nrn, int const ell, int const emm) {
      return nlm_index(numax, nrn, ell, emm); }

  inline constexpr
  int lmn_index(int const numax, int const ell, int const emm, int const nrn) {
      return ((3*numax + 5)*ell + 3*(1 + numax)*ell*ell - 2*ell*ell*ell)/6
              + emm*(1 + (numax - ell)/2) + nrn; }


  // =============================== Cartesian indices ===============================

  // Cartesian 3D index
  inline constexpr
  int zyx_index(int const numax, int const nx, int const ny, int const nz) {
      return (nz*nz*nz - 3*(2 + numax)*nz*nz + (3*pow2(2 + numax) - 1)*nz // contribution for nx=0, ny=0
              + ny*(2 + numax - nz)*6 - (ny*(1 + ny))*3  +  6*nx)/6;  // contribution from previous y and x
  } // zyx_index
  template <int numax, typename int_t> inline constexpr
  int zyx_index(int_t const nx, int_t const ny, int_t const nz) {
      return zyx_index(numax, nx, ny, nz); }
  template <typename int_t> inline constexpr
  int zyx_index(int const numax, int_t const nzyx[3]) {
      return zyx_index(numax, nzyx[0], nzyx[1], nzyx[2]); }
  template <int numax, typename int_t> inline constexpr
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
  int Enl_index(int const nrn, int const ell) {
      return (pow2(ell + 2*nrn + 1) + 2*ell)/4; } // (ell + 2*nrn)=nu, use ((nu + 1)^2)/4 as offset and add ell/2

  // energy-ordered radial 3D index
  inline constexpr
  int Elnm_index(int const ell, int const nrn, int const emm) {
      return ((ell + 2*nrn)*(ell + 2*nrn + 1)*(ell + 2*nrn + 2) // energy shell offset (nu*(nu+1)*(nu+2))/6
              + 3*ell*(ell - 1) // previous ells (ell*(ell - 1))/2
              + 6*(emm + ell))/6; } // linear emm-contribution

  template <typename int_t> inline constexpr
  int_t get_nu(int_t const nx, int_t const ny, int_t const nz) { return nx + ny + nz; }

  template <typename int_t> inline constexpr
  int_t get_nu(int_t const ell, int_t const nrn) { return ell + 2*nrn; }

  template <typename int_t> inline
  int get_nu(int_t const energy_ordered) {
      int nu{-1}; while (energy_ordered >= nSHO(nu)) { ++nu; } return nu; }

  template <typename int_t> inline
  status_t construct_index_table(
        int_t energy_ordered[] // permutation indices for energy ordering
      , int const numax // SHO basis size
      , SHO_order_t const order // SHO input order
      , int_t *inverse=nullptr // inverse permutation
      , int const echo=0 // log-level
  ) {
      // construct a table of energy ordered indices
      // this is needed to address e.g. the block-diagonal SHO-transformation operator
      if (echo > 3) std::printf("# %s for <numax=%i> order_%s\n",
                              __func__, numax, SHO_order2string(order).c_str());
      if (echo > 4) std::printf("# ");
      int ii{0};
      switch (order) {

        case order_zyx:
          for (int z = 0; z <= numax; ++z) {
              for (int y = 0; y <= numax - z; ++y) {
                  for (int x = 0; x <= numax - z - y; ++x) {
                      assert( zyx_index(numax, x, y, z) == ii );
                      int const eo = Ezyx_index(x, y, z);
                      energy_ordered[ii] = eo;
                      if (echo > 5) std::printf(" %i", eo);
                      if (inverse) inverse[eo] = ii;
                      ++ii;
          }}} // x y z
          assert(nSHO(numax) == ii);
        break;

        case order_lmn:
          for (int l = 0; l <= numax; ++l) {
              for (int m = -l; m <= l; ++m) {
                  for (int n = 0; n <= (numax - l)/2; ++n) {
                      assert( lmn_index(numax, l, m, n) == ii );
                      int const eo = Elnm_index(l, n, m);
                      energy_ordered[ii] = eo;
                      if (echo > 5) std::printf(" %i", eo);
                      if (inverse) inverse[eo] = ii;
                      ++ii;
          }}} // l m n
          assert(nSHO(numax) == ii);
        break;

        case order_lnm:
          for (int l = 0; l <= numax; ++l) {
              for (int n = 0; n <= (numax - l)/2; ++n) {
                  for (int m = -l; m <= l; ++m) {
                      assert( lnm_index(numax, l, n, m) == ii );
                      int const eo = Elnm_index(l, n, m);
                      energy_ordered[ii] = eo;
                      if (echo > 5) std::printf(" %i", eo);
                      if (inverse) inverse[eo] = ii;
                      ++ii;
          }}} // l n m
          assert(nSHO(numax) == ii);
        break;

        case order_nlm:
          for (int n = 0; n <= numax/2; ++n) {
              for (int l = 0; l <= numax - 2*n; ++l) {
                  for (int m = -l; m <= l; ++m) {
                      assert( nlm_index(numax, n, l, m) == ii );
                      int const eo = Elnm_index(l, n, m);
                      energy_ordered[ii] = eo;
                      if (echo > 5) std::printf(" %i", eo);
                      if (inverse) inverse[eo] = ii;
                      ++ii;
          }}} // n l m
          assert(nSHO(numax) == ii);
        break;

        case order_ln: // emm-degenerate
          for (int l = 0; l <= numax; ++l) {
              for (int n = 0; n <= (numax - l)/2; ++n) {
                  assert( ln_index(numax, l, n) == ii );
                  int const eo = Enl_index(n, l);
                  energy_ordered[ii] = eo;
                  if (echo > 4) std::printf(" %i", eo);
                  if (inverse) inverse[eo] = ii;
                  ++ii;
          }} // l n
          assert(nSHO_radial(numax) == ii);
        break;

        case order_nl: // emm-degenerate
          for (int n = 0; n <= numax/2; ++n) {
              for (int l = 0; l <= numax - 2*n; ++l) {
                  assert( nl_index(numax, n, l) == ii );
                  int const eo = Enl_index(n, l);
                  energy_ordered[ii] = eo;
                  if (echo > 4) std::printf(" %i", eo);
                  if (inverse) inverse[eo] = ii;
                  ++ii;
          }} // n l
          assert(nSHO_radial(numax) == ii);
        break;

        case order_Enl: // already energy-ordered emm-degenerate
          for (int ii = 0; ii < nSHO_radial(numax); ++ii) {
              energy_ordered[ii] = ii;
              if (inverse) inverse[ii] = ii;
          } // ii
          if (echo > 4) std::printf(" <unity> ");
        break;

        case order_Ezyx: // already energy-ordered
        case order_Elnm: // already energy-ordered
          for (int ii = 0; ii < nSHO(numax); ++ii) {
              energy_ordered[ii] = ii;
              if (inverse) inverse[ii] = ii;
          } // ii
          if (echo > 4) std::printf(" <unity> ");
        break;

        default:
            if (echo > 0) std::printf("# %s: no such case implemented: order_%s\n",
                                    __func__, SHO_order2string(order).c_str());
            return order; // error
      } // switch order
      if (echo > 4) std::printf("\n\n");
      return 0; // success if 0
  } // construct_index_table

  template <typename uint_t> inline
  status_t quantum_number_table(
        uint_t idx[] // assumed stride 4
      , int const numax // size of SHO basis
      , SHO_order_t const order // SHO input order
      , int const echo=0 // log-level
  ) {
      // list all triples of quantum numbers
      if (is_Cartesian(order)) {
          int ii{0};
          for (int z = 0; z <= numax; ++z) {
              for (int y = 0; y <= numax - z; ++y) {
                  for (int x = 0; x <= numax - z - y; ++x) {
                      assert( zyx_index(numax, x, y, z) == ii );
                      int const jj = (order == order_Ezyx) ? Ezyx_index(x, y, z) : ii;
                      idx[jj*4 + 0] = x;
                      idx[jj*4 + 1] = y;
                      idx[jj*4 + 2] = z;
                      idx[jj*4 + 3] = get_nu(x, y, z);
                      ++ii;
                  } // x
              } // y
          } // z
          assert (nSHO(numax) == ii);
      } else {
          if (echo > 0) std::printf("\n# %s %s Error: only Cartesian implemented but found %s\n\n",
                                    __FILE__, __func__, SHO_order2string(order).c_str());
          return -1;
      } // is_Cartesian
      return 0;
  } // quantum_number_table

  // sho_hex is a single character the encrypts 0..63 and matches hexadecimal below 16
  // 0123456789abcdefghijklmnopqrstuvABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`??????
  // 0123456789012345678901234567890123456789012345678901234567890123456789
  // 0000000000111111111122222222223333333333444444444455555555556666666666
  inline char sho_hex(unsigned i) { return (i<10)?('0'+i):((i<32)?('W'+i):((i<64)?('!'+i):'?')); }

  template <unsigned nChar=8> inline // use char[4] for Cartesian or emm_degenerate, use char[6] or char[8] for radial indices
  status_t construct_label_table(
        char label[]
      , int const numax
      , SHO_order_t const order
      , int const echo=0
  ) {

      auto const ellchar = "spdfghijklmno?????????????????????"; // ToDo: 'i' should not be included by convention
      int ii{0};
      switch (order) {

        case order_zyx:
        case order_Ezyx: // energy-ordered
          for (int z = 0; z <= numax; ++z) {
              for (int y = 0; y <= numax - z; ++y) {
                  for (int x = 0; x <= numax - z - y; ++x) {
                      int const j = is_energy_ordered(order) ? Ezyx_index(x, y, z) : ii;
                      std::snprintf(&label[j*nChar], nChar, "%c%c%c", sho_hex(z), sho_hex(y), sho_hex(x));
                      ++ii;
          }}} // x y z
          assert(nSHO(numax) == ii);
        break;

        case order_lmn:
        case order_lnm:
        case order_nlm:
        case order_Elnm: // energy-ordered
          for (int l = 0; l <= numax; ++l) {
              for (int m = -l; m <= l; ++m) {
                  for (int n = 0; n <= (numax - l)/2; ++n) {
                      int j{ii}; assert( lmn_index(numax, l, m, n) == ii );
                      if (order_Elnm == order) j = Elnm_index(l, n, m);
                      if (order_lnm == order)  j =  lnm_index(numax, l, n, m);
                      if (order_nlm == order)  j =  nlm_index(numax, n, l, m);
                      std::snprintf(&label[j*nChar], nChar, "%i%c%i", n, ellchar[l], m);
                      ++ii;
          }}} // l m n
          assert(nSHO(numax) == ii);
        break;

        case order_ln: // ell-ordered emm-degenerate
        case order_nl: // nrn-ordered emm-degenerate
        case order_Enl: // energy-ordered emm-degenerate
          for (int l = 0; l <= numax; ++l) {
              for (int n = 0; n <= (numax - l)/2; ++n) {
                  int j{-1};
                  if (is_energy_ordered(order)) { j = Enl_index(n, l); } else
                  if (order_nl == order) { j = nl_index(numax, n, l); } else
                  if (order_ln == order) { j = ii; assert( ln_index(numax, l, n) == ii ); }
                  assert(j >= 0);
                  std::snprintf(&label[j*nChar], nChar, "%c%i", ellchar[l], n);
                  ++ii;
          }} // l n
          assert(nSHO_radial(numax) == ii);
        break;

        default:
            if (echo > 0) std::printf("# %s: no such case implemented: order_%s\n",
                                        __func__, SHO_order2string(order).c_str());
            return order; // error
      } // switch order
      return 0; // success
  } // construct_label_table

#ifdef NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_order_enum(int const echo=4) {
      SHO_order_t const ord[9] = {order_zyx, order_Ezyx, order_lmn, order_lnm,
                                  order_nlm, order_Elnm, order_ln, order_Enl, order_nl};
      for (int io = 0; io < 9; ++io) {
          SHO_order_t const oi = ord[io];
          if (echo > 3) std::printf("# %s: SHO_order_t %s\t= 0x%x\t= %10lli  %s-ordered emm-%s %s\n",
                      __func__, SHO_order2string(oi).c_str(), unsigned(oi), oi
                       , is_energy_ordered(oi)?"energy":"  loop"
                       , is_emm_degenerate(oi)?"degenerate":"resolved  "
                       , is_Cartesian(oi)?"Cartesian":"Radial");
      } // io
      return 0;
  } // test_order_enum

  inline status_t test_radial_indices(int const echo=4, int const numax_max=99) {
      status_t nerrors(0);
      for (int numax = 0; numax <= numax_max; ++numax) {
          if (echo > 6) std::printf("\n# %s: numax == %i\n", __func__, numax);
          int lnm{0}, ln{0}, lm{0}, lmn{0};
          for (int ell = 0; ell <= numax; ++ell) {
              for (int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
                  assert(ell + 2*nrn == get_nu(ell, nrn)); // test get_nu
                  int const k = ln_index(numax, ell, nrn);
                  if ((echo > 7) && (k != ln)) std::printf("# ln_index<%i>(ell=%i, nrn=%i) == %i %i diff=%i\n", numax, ell, nrn, ln, k, k - ln);
                  assert(k == ln);
                  nerrors += (k != ln);
                  ++ln;
                  for (int emm = -ell; emm <= ell; ++emm) {
                      int const k = lnm_index(numax, ell, nrn, emm);
                      if (echo > 8) std::printf("# lnm_index<%i>(ell=%i, nrn=%i, emm=%i) == %i %i diff=%i\n", numax, ell, nrn, emm, lnm, k, k - lnm);
                      assert(k == lnm);
                      nerrors += (k != lnm);
                      ++lnm;
                  } // emm
              } // nrn
              for (int emm = -ell; emm <= ell; ++emm) {
                  int const k = lm_index(ell, emm);
                  if (echo > 7) std::printf("# lm_index(ell=%i, emm=%i) == %i %i diff=%i\n", ell, emm, lm, k, k - lm);
                  assert(k == lm);
                  nerrors += (k != lm);
                  ++lm;
                  for (int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
                      int const k = lmn_index(numax, ell, emm, nrn);
                      if (echo > 8) std::printf("# lmn_index<%i>(ell=%i, emm=%i, nrn=%i) == %i %i diff=%i\n", numax, ell, emm, nrn, lmn, k, k - lmn);
                      assert(k == lmn);
                      nerrors += (k != lmn);
                      ++lmn;
                  } // nrn
              } // emm
              assert(pow2(1 + ell) == lm); // checksum
          } // ell
          assert(nSHO(numax) == lnm); // checksum
          assert(nSHO(numax) == lmn); // checksum

          // the nlm_index is nrn-first ordered, not energy-ordered
          int nlm{0}, nl{0};
          for (int nrn = 0; nrn <= numax/2; ++nrn) {
              for (int ell = 0; ell <= numax - 2*nrn; ++ell) {
                  {
                      int const k = nlm_index(numax, nrn, ell, -ell);
                      if (echo > 7) std::printf("# nlm_index<%i>(nrn=%i, ell=%i, emm=-ell) == %i %i diff=%i\n", numax, nrn, ell, nlm, k, nlm - k);
                      assert(k == nlm);
                      nerrors += (k != nlm)*(2*ell + 1);
                      nlm += (2*ell + 1); // forward one ell-shell emm=-ell...ell
                  }
                  int const k = nl_index(numax, nrn, ell); // emm-degenerate nrn-ordered
                  if (echo > 6 && ell == 0) std::printf("# nl_index<%i>(nrn=%i, ell=%i) == %i %i diff=%i\n", numax, nrn, ell, nl, k, nl - k);
                  assert(k == nl);
                  nerrors += (k != nl);
                  ++nl;
              } // ell
          } // nrn
          assert(nSHO(numax) == nlm); // checksum
          if (echo > 6) std::printf("\n# lmn_index<%i>\n", numax);
      } // numax
      if (nerrors && echo > 1) std::printf("# Warning: %s found %i errors!\n", __func__, nerrors);
      return nerrors;
  } // test_radial_indices

  inline status_t test_Cartesian_indices(int const echo=3, int const numax_max=9) {
      status_t nerrors(0);
      for (int numax = 0; numax <= numax_max; ++numax) {
          if (echo > 6) std::printf("\n# %s: numax == %i\n", __func__, numax);
          int zyx = 0;
          for (int nz = 0; nz <= numax; ++nz) {
              for (int ny = 0; ny <= numax - nz; ++ny) {
                  for (int nx = 0; nx <= numax - nz - ny; ++nx) {
                      int const k = zyx_index(numax, nx, ny, nz);
                      if (echo > 8) std::printf("# zyx_index<%i>(nx=%i, ny=%i, nz=%i) == %i %i diff=%i\n", numax, nx, ny, nz, zyx, k, k - zyx);
                      assert(k == zyx);
                      nerrors += (k != zyx);
                      ++zyx;
                  } // nx
              } // ny
          } // nz
          assert(nSHO(numax) == zyx); // checksum
      } // numax
      if (nerrors && echo > 1) std::printf("# Warning: %s found %i errors!\n", __func__, nerrors);
      return nerrors;
  } // test_Cartesian_indices

  inline status_t test_energy_ordered_indices(int const echo=4, int const numax=9) {
      if (echo > 6) std::printf("\n# %s: numax == %i\n", __func__, numax);
      status_t nerrors(0);
      int nzyx{0}, nln{0}, nlnm{0};
      for (int nu = 0; nu <= numax; ++nu) { // shell index

          // Cartesian energy-ordered index
          if (echo > 7) std::printf("\n# Ezyx_index<nu=%i>\n", nu);
          int xyz{0};
          for (int nz = 0; nz <= nu; ++nz) { // z can also run backwards so we start from (0,0,nu) and proceed with (1,0,nu-1)
              for (int nx = 0; nx <= nu - nz; ++nx) {
                  int const ny = nu - nz - nx;
                  int const k = Ezyx_index(nx, ny, nz);
                  if ((echo > 6) && (k != nzyx))
                      std::printf("# Ezyx_index<nu=%i>(nx=%i, ny=%i, nz=%i) == %i %i diff=%i  xyz=%i %i\n",
                             nu, nx, ny, nz, nzyx, k, k - nzyx, xyz,  nx + (nz*((2+nu)*2-(nz + 1)))/2 );
                  assert(k == nzyx);
                  nerrors += (k != nzyx);
                  if (get_nu(nzyx) != nu) std::printf("# get_nu(%i) = %i but expected %i\n", nzyx, get_nu(nzyx), nu);
                  assert(get_nu(nzyx) == nu);
                  ++nzyx;
                  ++xyz;
              } // ny
          } // nz
          assert(nSHO(nu) == nzyx); // checksum

          // radial emm-degenerate energy-ordered indices
          for (int ell = nu%2; ell <= nu; ell+=2) {
              int const nrn = (nu - ell)/2;
              int const k = Enl_index(nrn, ell);
              if (echo > 9) std::printf("# Enl_index<nu=%i>(nrn=%i, ell=%i) == %i %i\n", nu, nrn, ell, nln, k);
              assert(k == nln);
              ++nln;
              for (int emm = -ell; emm <= ell; ++emm) {
                  int const k = Elnm_index(ell, nrn, emm);
                  if (echo > 9) std::printf("# Elnm_index<nu=%i>(ell=%i, nrn=%i, emm=%i) == %i\n", nu, ell, nrn, emm, nlnm);
                  assert(k == nlnm);
                  nerrors += (k != nlnm);
                  assert(nu == get_nu(nlnm));
                  ++nlnm;
              } // emm
          } // ell
          assert(nSHO(nu) == nlnm); // checksum
          assert(nSHO_radial(nu) == nln); // checksum

      } // nu
      if (nerrors && echo > 1) std::printf("# Warning: %s found %i errors!\n", __func__, nerrors);
      return nerrors;
  } // test_energy_ordered_indices

  template <typename int_t>
  inline status_t test_index_table_construction(int const echo=1) {
      status_t stat(0);
      int const numax_max = 9;
      if (echo > 6) std::printf("\n# %s: numax == %i\n", __func__, numax_max);
      SHO_order_t const orders[] = {order_zyx, order_Ezyx, order_lmn, order_nlm, order_lnm, order_Elnm, order_ln, order_Enl, order_nl};
      for (int io = 0; io < 9; ++io) {
          auto const order = orders[io];
          if (echo > 6) std::printf("# %s order_%s\n", __func__, SHO_order2string(order).c_str());

          for (int numax = 0; numax <= numax_max; ++numax) {
              int const nsho = is_emm_degenerate(order) ? nSHO_radial(numax) : nSHO(numax);
              std::vector<int_t> list(nsho, -1), inv_list(nsho, -1);
              stat += construct_index_table(list.data(), numax, order, inv_list.data(), echo);
              std::vector<char> label(nsho*8);
              stat += construct_label_table(label.data(), numax, order);
              if (echo > 7) std::printf("# %s numax=%i order_%s labels:  ",
                              __func__, numax, SHO_order2string(order).c_str());
              for (int ii = 0; ii < nsho; ++ii) {
                  assert( list[inv_list[ii]] == ii ); // check that the lists are inverse of each other
                  assert( inv_list[list[ii]] == ii ); // check that the lists are inverse of each other
                  if (echo > 7) std::printf(" %s", &label[ii*8]);
              } // ii
              if (echo > 7) std::printf("\n");
          } // numax
      } // io
      if (stat && echo > 1) std::printf("# Warning: %s found %i errors!\n", __func__, stat);
      return stat;
  } // test_index_table_construction

  inline status_t test_sho_hex(int const echo=1) {
      if (echo < 1) return 0;
      std::printf("# %s: for i in 0..69: sho_hex= ", __func__);
      for (int i = 0; i < 70; ++i) {
          std::printf("%c", sho_hex(i));
      } // i
      std::printf("\n\n");
      return 0;
  } // test_sho_hex

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_radial_indices(echo);
      stat += test_Cartesian_indices(echo);
      stat += test_energy_ordered_indices(echo);
      stat += test_index_table_construction<int16_t>(echo);
      stat += test_order_enum(echo);
      stat += test_sho_hex(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace sho_tools
