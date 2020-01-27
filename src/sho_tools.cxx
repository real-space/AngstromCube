#include <cstdio> // printf
#include <cassert> // assert
#include <cstdint> // int16_t
#include <vector> // std::vector<T>

#include "sho_tools.hxx"

namespace sho_tools {
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_order_enum(int const echo=4) {
      status_t nerrors = 0;
      SHO_order_t const ord[9] = {order_zyx, order_Ezyx, order_lmn, order_lnm, 
                  order_nlm, order_Elnm, order_ln, order_Enl, order_unknown};
      for(int io = 0; io < 9; ++io) {
          SHO_order_t const oi = ord[io];
          if (echo > 0) printf("# %s: SHO_order_t %s\t= 0x%x\t= %10li\n", __func__, SHO_order2string(&oi), (unsigned)oi, oi);
      } // io
      if (nerrors && echo > 1) printf("# Warning: %s found %d errors!\n", __func__, nerrors);
      return nerrors;
  } // test_order_enum

  status_t test_radial_indices(int const echo=4) {
      status_t nerrors = 0;
      for(int numax = 0; numax <= 9; ++numax) {
          if (echo > 6) printf("\n# %s: numax == %d\n", __func__, numax);
          int lnm = 0, ln = 0, lm = 0, lmn = 0;
          for(int ell = 0; ell <= numax; ++ell) {
              for(int nrn = 0; nrn <= (numax - ell)/2; ++nrn) {
                  assert(ell + 2*nrn == get_nu(ell, nrn)); // test get_nu
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
          assert(nSHO(numax) == lnm); // checksum
          assert(nSHO(numax) == lmn); // checksum

          // the nlm_index is nrn-first ordered, not energy-ordered
          int nlm = 0;
          for(int nrn = 0; nrn <= numax/2; ++nrn) {
              for(int ell = 0; ell <= numax - 2*nrn; ++ell) {
                  int const k = nlm_index(numax, nrn, ell, -ell);
                  if (echo > 7) printf("# nlm_index<%i>(nrn=%i, ell=%d, emm=-ell) == %d %d diff=%d\n", numax, nrn, ell, nlm, k, nlm - k);
                  assert(k == nlm);
                  nerrors += (k != nlm)*(2*ell + 1);
                  nlm += (2*ell + 1); // forward one ell-shell emm=-ell...ell
              } // ell
          } // nrn
          assert(nSHO(numax) == nlm); // checksum
          if (echo > 6) printf("\n# lmn_index<%d>\n", numax);
      } // numax
      if (nerrors && echo > 1) printf("# Warning: %s found %d errors!\n", __func__, nerrors);
      return nerrors;
  } // test_radial_indices

  status_t test_Cartesian_indices(int const echo=3) {
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
          assert(nSHO(numax) == zyx); // checksum
      } // numax
      if (nerrors && echo > 1) printf("# Warning: %s found %d errors!\n", __func__, nerrors);
      return nerrors;
  } // test_Cartesian_indices

  status_t test_energy_ordered_indices(int const echo=4) {
      status_t nerrors = 0;
      int const numax = 9;
      if (echo > 6) printf("\n# %s: numax == %d\n", __func__, numax);
      int nzyx = 0, nln = 0, nlnm = 0;
      for(int nu = 0; nu <= numax; ++nu) { // shell index

          // Cartesian energy-ordered index
          if (echo > 7) printf("\n# nzyx_index<nu=%d>\n", nu);
          int xyz = 0;
          for(int nz = 0; nz <= nu; ++nz) { // z can also run backwards so we start from (0,0,nu) and proceed with (1,0,nu-1)
              for(int nx = 0; nx <= nu - nz; ++nx) {
                  int const ny = nu - nz - nx;
                  int const k = Ezyx_index(nx, ny, nz);
                  if ((echo > 6) && (k != nzyx))
                      printf("# Ezyx_index<nu=%d>(nx=%d, ny=%d, nz=%d) == %d %d diff=%d  xyz=%d %d\n", 
                             nu, nx, ny, nz, nzyx, k, k - nzyx, xyz,  nx + (nz*((2+nu)*2-(nz + 1)))/2 );
                  assert(k == nzyx);
                  nerrors += (k != nzyx);
                  if (get_nu(nzyx) != nu) printf("# get_nu(%d) = %d but expected %d\n", nzyx, get_nu(nzyx), nu);
                  assert(get_nu(nzyx) == nu);
                  ++nzyx;
                  ++xyz;
              } // ny
          } // nz
          assert(nSHO(nu) == nzyx); // checksum
          
          // radial emm-degenerate energy-ordered indices
          for(int ell = nu%2; ell <= nu; ell+=2) {
              int const nrn = (nu - ell)/2;
              int const k = Enl_index(nrn, ell);
              if (echo > 9) printf("# Enl_index<nu=%d>(nrn=%d, ell=%d) == %d %d\n", nu, nrn, ell, nln, k);
              assert(k == nln);
              ++nln;
              for(int emm = -ell; emm <= ell; ++emm) {
                  int const k = Elnm_index(ell, nrn, emm);
                  if (echo > 9) printf("# Elnm_index<nu=%d>(ell=%d, nrn=%d, emm=%d) == %d\n", nu, ell, nrn, emm, nlnm);
                  assert(k == nlnm);
                  nerrors += (k != nlnm);
                  assert(nu == get_nu(nlnm));
                  ++nlnm;
              } // emm
          } // ell
          assert(nSHO(nu) == nlnm); // checksum
          assert(nSHO_radial(nu) == nln); // checksum

      } // nu
      if (nerrors && echo > 1) printf("# Warning: %s found %d errors!\n", __func__, nerrors);
      return nerrors;
  } // test_energy_ordered_indices

  template<typename int_t>
  status_t test_index_table_construction(int const echo=1) {
      status_t stat = 0;
      int const numax_max = 9;
      if (echo > 6) printf("\n# %s: numax == %d\n", __func__, numax_max);
      SHO_order_t const orders[] = {order_zyx, order_Ezyx, order_lmn, order_nlm, order_lnm, order_Elnm, order_ln, order_Enl};
      for(int io = 0; io < 8; ++io) {
          auto const order = orders[io];
          if (echo > 6) printf("# %s order_%s\n", __func__, SHO_order2string(&order));

          for(int numax = 0; numax <= numax_max; ++numax) {
              int const nsho = is_emm_degenerate(order) ? nSHO_radial(numax) : nSHO(numax);
              std::vector<int_t> list(nsho, -1), inv_list(nsho, -1);
              stat += construct_index_table(list.data(), numax, order, inv_list.data(), echo);
              std::vector<char> label(nsho*8, '\0');
              stat += construct_label_table(label.data(), numax, order);
              if (echo > 7) printf("# %s numax=%i order_%s labels:  ", __func__, numax, SHO_order2string(&order));
              for(int ii = 0; ii < nsho; ++ii) {
                  assert( list[inv_list[ii]] == ii ); // check that the lists are inverse of each other
                  assert( inv_list[list[ii]] == ii ); // check that the lists are inverse of each other
                  if (echo > 7) printf(" %s", &label[ii*8]);
              } // ii
              if (echo > 7) printf("\n");
          } // numax
      } // io
      if (stat && (echo > 1)) printf("# Warning: %s found %d errors!\n", __func__, stat);
      return stat;
  } // test_index_table_construction
  
  status_t all_tests(int const echo) {
    auto status = 0;
    status += test_order_enum(echo);
    status += test_radial_indices(echo);
    status += test_Cartesian_indices(echo);
    status += test_energy_ordered_indices(echo);
    status += test_index_table_construction<int16_t>(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  
  
} // namespace sho_tools
