
#include <stdlib.h> // abs
#include <cmath> // fabs, exp
#include <cstdio> // printf
#include <algorithm> // min, max

#include "radial_grid.hxx"
#include "radial_grid.h" // radial_grid_t

#include "inline_tools.hxx" // align


namespace radial_grid {

  inline radial_grid_t* get_memory(size_t const n_aligned) {
      auto g = new radial_grid_t;
      g->r = new double[5*n_aligned];
      g->dr     = &g->r[1*n_aligned];
      g->rdr    = &g->r[2*n_aligned];
      g->r2dr   = &g->r[3*n_aligned];
      g->rinv   = &g->r[4*n_aligned];
      g->memory_owner = (nullptr != g->r);
      return g;
  } // get_memory
  
  inline void set_derived_grid_quantities(radial_grid_t& g, int const n_aligned) {
      for (int ir = 0; ir < n_aligned; ++ir) {
          auto const r = g.r[ir], dr = g.dr[ir];
          g.rdr[ir] = r*dr;
          g.r2dr[ir] = r*r*dr;
          g.rinv[ir] = (ir)? 1./r : 0;
//        printf("%g %g\n", r, dr); // DEBUG
      } // ir
  } // set_derived_grid_quantities

  radial_grid_t* create_exponential_radial_grid(int const npoints,
                   float const Rmax, float const anisotropy) { // optional args
      double const R = std::max(std::abs(Rmax)*1., .945);
      double const d = std::min(std::max(1e-4, anisotropy*1.), .1);
      int const N = std::max(std::abs(npoints), 32);

      int const N_aligned = align<2>(N); // padded to multiples of 4
      auto const g = get_memory(N_aligned);    
      
      double const a = R / (std::exp(d*(N - 1)) - 1.); // prefactor
      for (auto ir = 0; ir < N_aligned; ++ir) {
          double const edi = std::exp(d*ir);
          double const r = a*(edi - 1.);
          double const dr = a*d*edi * (ir < N);
          g->r[ir] = r;
          g->dr[ir] = dr;
      } // ir
      set_derived_grid_quantities(*g, N_aligned);
      g->n = N;
      g->rmax = g->r[g->n - 1];
      return g;
  } // create_exponential_radial_grid

  radial_grid_t* create_equidistant_radial_grid(int const npoints, float const Rmax) {
      int const N = std::max(1, npoints);
      int const N_aligned = align<2>(N); // padded to multiples of 4
      auto const g = get_memory(N_aligned);    

      double const dr = Rmax/N;
      for (auto ir = 0; ir < N_aligned; ++ir) {
          double const r = ir*dr;
          g->r[ir] = r;
          g->dr[ir] = dr;
      } // ir
      set_derived_grid_quantities(*g, N_aligned);
      g->n = N;
      g->rmax = g->r[g->n - 1];
      return g;
  } // create_equidistant_radial_grid
  
  
  radial_grid_t* create_pseudo_radial_grid(radial_grid_t const &tru, double const r_min, int const echo) {
      // find a suitable grid point to start from
      int ir = 0; while (tru.r[ir] < r_min) ++ir;
      if (echo > 3) printf("# start pseudo grid from r[%d]=%g Bohr\n", ir, tru.r[ir]);
      int const ir_offset = ir;

      auto g = new radial_grid_t;
      // offset pointers
      g->r      = tru.r    + ir_offset;
      g->dr     = tru.dr   + ir_offset;
      g->rdr    = tru.rdr  + ir_offset;
      g->r2dr   = tru.r2dr + ir_offset;
      g->rinv   = tru.rinv + ir_offset;
      g->memory_owner = false;

      g->n = tru.n - ir_offset; // reduced number of grid points
      g->rmax = tru.rmax; // both grids have the same tail
      return g;
  } // create_pseudo_radial_grid
  
  
  void destroy_radial_grid(radial_grid_t* g) {
      if (g->memory_owner) delete [] g->r;
      g->n = 0;
  } // destroy

  int find_grid_index(radial_grid_t const &g, double const radius) {
      // ToDo: if this becomes performance critical, replace by bisection algorithm
      for(int ir = 0; ir < g.n; ++ir) {
          if (radius < g.r[ir]) return ir - 1;
      } // ir
      return g.n - 2;
  } // find_grid_index
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo=9) {
      if (echo > 0) printf("\n# %s: \n", __func__);
      auto g = create_exponential_radial_grid(1 << 11);
      destroy_radial_grid(g);
      return 0;
  } // test_create_and_destroy

  status_t test_exp_grid(int const echo=3) {
      if (echo > 0) printf("\n# %s: \n", __func__);
      int const n = 1 << 11;
      auto g = create_exponential_radial_grid(n);
      if (echo > 9) {
          printf("\n## radial exponential grid (N=%d, anisotropy=%g, Rmax=%g %s): r, dr, 1/r\n", 
                    n, default_anisotropy, g->rmax*1.0, " Bohr");
          for(int ir = 0; ir < g->n; ++ir) {
              if (echo > 11 || ir < 3 || ir > g->n - 4) printf("%g %g %g\n", g->r[ir], g->dr[ir], g->rinv[ir]);
          } // ir
      } // echo
      return 0;
  } // test_exp_grid

  status_t all_tests(int const echo) {
      status_t status(0);
      status += test_create_and_destroy(echo);
      status += test_exp_grid(echo);
      return status;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace radial_grid
