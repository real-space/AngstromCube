#pragma once

#include <cstdio> // printf, std::sprintf
#include <cassert> // assert
// #include <cmath> // std::min, std::max

#include "status.hxx" // status_t

namespace parallel_domains {
  
  template <int D0=1>
  inline status_t decompose_grid(unsigned const ng, int const echo=0, int const min_ng_per_pe=4) {
      status_t stat{0};
      for(int npe = 1; npe <= ng + 9; ++npe) {
          for(int bc = 0; bc <= 1; ++bc) { // 0:isolated, 1:periodic

              int ng_more, np_more, ng_less, np_less, np_zero{0}, np_uppr{0};
              // how to parallelize ng grid points onto npe process elements?
              //
              if (1) {
                  // get a good load balance like this
                  ng_less = ng/npe; // floor by integer divide
//                ng_less = std::max(ng_less, min_ng_per_pe); // should come in here, ToDo
                  ng_more = ng_less + 1;
                  /**
                  *  Solve
                  *    g_< p_< + g_> p_> = G
                  *  while
                  *        p_< + p_> + p_0 = P
                  *  so
                  *    g_< (P - p_>) + g_> p_> = G
                  *  reads
                  *    p_> = (G - g_< P)/(g_> - g_<)
                  *  if p_0 == 0
                  */
                  np_more = (ng - ng_less*npe);
                  np_less = npe - np_more;

                  if (0 == bc && ng_less > 0) {
                      // for isolated BC distribute the number of process elements 
                      // with more grid points to the left and right. This is
                      // accounting for the distribution of atoms in the cell
                      // center and vacuum buffer regions at the borders.
                      np_uppr = np_more/2;
                      // with np_uppr <= np_lowr a distribution is found where
                      // the master (process element #0) has ng_more grid points
                  } // bc

              } // method
              int const np_lowr = np_more - np_uppr;
              
              int const gs[] = {ng_more, ng_less, ng_more, 0};
              int const ps[] = {np_lowr, np_less, np_uppr, np_zero};
              int gxp[3];
              char dec_str[3][16];
              for(int i = 0; i < 3; ++i) {
                  gxp[i] = gs[i] * ps[i];
                  if (gxp[i] > 0) {
                      std::sprintf(dec_str[i], "%dx%d", gs[i], ps[i]);
                  } else {
                      std::sprintf(dec_str[i], "0");
                  } // nonzero
              } // section i {lower, middle, upper}
              if (echo > 7) printf("# parallelize %d grid points as %s + %s + %s with %d process elements (BC=%s)\n", 
                          ng, dec_str[0], dec_str[1], dec_str[2], npe, bc?"periodic":"isolated");
              
#ifdef DEBUG
              // check that this distribution routine works correctly
              assert(ng == ng_more*np_lowr + ng_less*np_less + ng_more*np_uppr);
              assert(npe == np_lowr + np_less + np_uppr);
              assert(npe == np_more + np_less);
#endif
              stat += (ng != ng_more*np_lowr + ng_less*np_less + ng_more*np_uppr);
              stat += (npe != np_lowr + np_less + np_uppr);
              stat += (npe != np_more + np_less);
          } // bc
      } // npe
      return stat;
  } // decompose_grid

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_analysis(int const echo=0) {
      status_t stat{0};
      for(int ng = 1; ng <= 99; ++ng) {
          stat += decompose_grid(ng, echo);
      } // ng
      return stat;
  } // test_analysis

  inline status_t all_tests(int const echo=0) { 
      status_t stat{0};
      stat += test_analysis(echo);
      return stat; 
  } // all_tests
  
#endif // NO_UNIT_TESTS

} // namespace parallel_domains
