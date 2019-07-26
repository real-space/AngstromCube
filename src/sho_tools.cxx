
#include <cstdio> // printf
#include <cassert> // assert
#include <cstdint> // int16_t

#include "sho_tools.hxx"

namespace sho_tools {
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_order_enum(int const echo=4) {
      status_t nerrors = 0;
      SHO_order_t const ord[8] = {order_zyx, order_lmn, order_lnm, order_nlnm, order_nzyx, order_ln, order_nln, order_unknown};
      for(int io = 0; io < 8; ++io) {
          SHO_order_t const oi = ord[io];
          if (echo > 3) printf("# %s: SHO_order_t %c%c%c%c = 0x%8x = %10u\n", __func__, oi>>24, oi>>16, oi>>8, oi, oi, oi);
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
                  if (get_nu(nzyx) != nu) printf("# get_nu(%d) = %d but expected %d\n", nzyx, get_nu(nzyx), nu);
                  assert(get_nu(nzyx) == nu);
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
                  assert(nu == get_nu(nlnm));
                  ++nlnm;
              } // emm
          } // ell
          assert(nSHO(nu) == nlnm); // checksum
          
      } // nu
      if (nerrors && echo > 1) printf("# Warning: %s found %d errors!\n", __func__, nerrors);
      return nerrors;
  } // test_energy_ordered_indices

  template<typename int_t>
  status_t test_index_table_construction(int const echo=1) {
      status_t stat = 0;
      int const numax_max = 9;
      if (echo > 6) printf("\n# %s: numax == %d\n", __func__, numax_max);
      SHO_order_t const orders[3] = {order_zyx, order_lmn, order_ln};
      for(int io = 0; io < 3; ++io) {
          auto const order = orders[io];
          for(int numax = 0; numax <= numax_max; ++numax) {
              int const nsho = (order == order_ln)? pow2(numax + 2)/4 : nSHO(numax);
              auto const list     = new int_t[nsho];
              auto const inv_list = new int_t[nsho];
              stat += construct_index_table(list, numax, order, inv_list, echo);
              for(int ii = 0; ii < nsho; ++ii) {
                  assert(list[inv_list[ii]] == ii); // check that the lists are inverse of each other
                  assert(inv_list[list[ii]] == ii); // check that the lists are inverse of each other
              } // ii
              delete[] inv_list;
              delete[] list;
          } // numax
      } // io
      if (stat && (echo > 1)) printf("# Warning: %s found %d errors!\n", __func__, stat);
      return stat;
  } // test_index_table_construction
  
  status_t all_tests() {
    auto status = 0;
    status += test_order_enum();
    status += test_radial_indices();
    status += test_Cartesian_indices();
    status += test_energy_ordered_indices();
    status += test_index_table_construction<int16_t>();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  
  
} // namespace sho_tools
