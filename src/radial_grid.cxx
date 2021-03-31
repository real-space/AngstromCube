
#include <cmath> // std::abs, std::exp
#include <cstdio> // std::printf
#include <algorithm> // std::min, std::max

#include "radial_grid.hxx" // radial_grid_t

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

  inline void set_derived_grid_quantities(radial_grid_t & g, int const n_aligned) {
      for (int ir = 0; ir < n_aligned; ++ir) {
          auto const r = g.r[ir], dr = g.dr[ir];
          g.rdr[ir] = r*dr;
          g.r2dr[ir] = r*r*dr;
          g.rinv[ir] = (ir)? 1./r : 0;
//        std::printf("%g %g\n", r, dr); // DEBUG
      } // ir
  } // set_derived_grid_quantities

  char const equation_exponential[] = "r=a*(exp(d*i)-1)";
  char const equation_equidistant[] = "r=a*i/n";

  // ToDo: remove exponential from the function name
  radial_grid_t* create_exponential_radial_grid(
        int const npoints
      , float const Rmax // =default_Rmax
      , float const anisotropy // =default_anisotropy
  ) {
      int const n = std::max(std::abs(npoints), 32);
      double const R = std::max(std::abs(Rmax)*1., .945);

      int const n_aligned = align<2>(n); // padded to multiples of 4
      auto const g = get_memory(n_aligned);
#ifdef  USE_RECIPROCAL_RADIAL_GRID
      bool static warn_reciprocal{true};
      if (warn_reciprocal) {
          std::printf("# -D USE_RECIPROCAL_RADIAL_GRID\n");
          warn_reciprocal = false;
      } // warn_reciprocal
      for (auto i = 0; i < n; ++i) {
          double const rec = 1./(n*(n - i));
          g->r[i]  = R*i*rec;
          g->dr[i] = R*n*rec*n*rec;
      } // i
      for (auto ir = n; ir < n_aligned; ++ir) {
          g->r[ir]  = 0;
          g->dr[ir] = 0;
      } // ir
      g->equation = "r=a*i/(n-i)";
      double const d = 0; // anisotropy
#else
      double const d = std::min(std::max(1e-4, anisotropy*1.), .1);
      double const a = R / (std::exp(d*(n - 1)) - 1.); // prefactor
      for (auto ir = 0; ir < n_aligned; ++ir) {
          double const edi = std::exp(d*ir);
          g->r[ir]  = a*(edi - 1.);
          g->dr[ir] = a*d*edi * (ir < n);
      } // ir
      g->equation = equation_exponential;
#endif
      set_derived_grid_quantities(*g, n_aligned);
      g->n = n;
      g->rmax = g->r[g->n - 1];
      g->anisotropy = d;
      return g;
  } // create_exponential_radial_grid

  double get_prefactor(radial_grid_t const & g) {
      if (equation_exponential == g.equation) {
          return g.rmax/(std::exp(g.anisotropy*(g.n - 1)) - 1.);
      } else {
          return g.rmax/g.n;
      }
  } // get_prefactor

  radial_grid_t* create_equidistant_radial_grid(
        int const npoints
      , float const Rmax // =default_Rmax
  ) {
      int const n = std::max(1, npoints);
      int const n_aligned = align<2>(n); // padded to multiples of 4
      auto const g = get_memory(n_aligned);    

      double const dr = Rmax/n;
      for (auto ir = 0; ir < n_aligned; ++ir) {
          g->r[ir] = ir*dr;
          g->dr[ir] = dr;
      } // ir
      set_derived_grid_quantities(*g, n_aligned);
      g->n = n;
      g->rmax = g->r[g->n - 1];
      g->anisotropy = 0;
      g->equation = equation_equidistant;
      return g;
  } // create_equidistant_radial_grid

  radial_grid_t* create_pseudo_radial_grid(
        radial_grid_t const & tru
      , double const r_min // =1e-3 Bohr
      , int const echo // log-level
  ) {
      // find a suitable grid point to start from
      int ir{0}; while (tru.r[ir] < r_min) { ++ir; }
      if (echo > 3) std::printf("# start pseudo grid from r[%d]=%g Bohr\n", ir, tru.r[ir]);
      int const nr_diff = ir;

      auto g = new radial_grid_t;
      // offset pointers
      g->r      = tru.r    + nr_diff;
      g->dr     = tru.dr   + nr_diff;
      g->rdr    = tru.rdr  + nr_diff;
      g->r2dr   = tru.r2dr + nr_diff;
      g->rinv   = tru.rinv + nr_diff;
      g->memory_owner = false; // avoid double free

      g->n = tru.n - nr_diff; // reduced number of grid points
      g->rmax = tru.rmax; // both grids have the same tail
      g->anisotropy = tru.anisotropy;
      return g;
  } // create_pseudo_radial_grid

  void destroy_radial_grid(radial_grid_t* g, char const *name) {
//    std::printf("\n# %s name=%s memory_owner= %i\n\n", __func__, name, g->memory_owner);
      if (g->memory_owner) delete [] g->r;
      g->n = 0;
  } // destroy

  int find_grid_index(radial_grid_t const & g, double const radius) {
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
      if (echo > 0) std::printf("\n# %s: \n", __func__);
      auto const g = create_exponential_radial_grid(1 << 11);
      destroy_radial_grid(g);
      return 0;
  } // test_create_and_destroy

  status_t test_exp_grid(int const echo=3) {
      if (echo > 0) std::printf("\n# %s: \n", __func__);
      int const n = 1 << 11;
      auto & g = *create_exponential_radial_grid(n);
      double integ[] = {0, 0, 0};
      if (echo > 3) std::printf("\n## radial exponential grid (N=%d, anisotropy=%g, Rmax=%g %s):"
              " r, dr, 1/r, integrals over {1,r,r^2}, r, r^2/2, r^3/3\n", n, g.anisotropy, g.rmax*1.0, " Bohr");
      for(int ir = 0; ir < g.n; ++ir) {
          if (echo > 3 && (ir < 3 || ir > g.n - 4 || echo > 11)) {
              auto const r = g.r[ir];
              std::printf("%g %g %g %g %g %g %g %g %g\n", r, g.dr[ir], g.rinv[ir],
                    integ[0], integ[1], integ[2], r, 0.5*r*r, r*r*r/3);
          } // echo
          integ[0] += g.dr[ir];
          integ[1] += g.rdr[ir];
          integ[2] += g.r2dr[ir];
      } // ir
      destroy_radial_grid(&g);
      return (std::abs(integ[0] - g.rmax) > .05); // integral dr up to Rmax should match Rmax
  } // test_exp_grid

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_create_and_destroy(echo);
      stat += test_exp_grid(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace radial_grid
