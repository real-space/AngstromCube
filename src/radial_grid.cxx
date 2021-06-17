
#include <cmath> // std::abs, std::exp
#include <cstdio> // std::printf
#include <algorithm> // std::min, std::max

#include "radial_grid.hxx" // radial_grid_t

#include "inline_math.hxx" // align

namespace radial_grid {
  
  inline radial_grid_t* get_memory(size_t const nr_aligned) {
      auto g = new radial_grid_t;
      g->r = new double[5*nr_aligned];
      g->dr    = & g->r[1*nr_aligned];
      g->rdr   = & g->r[2*nr_aligned];
      g->r2dr  = & g->r[3*nr_aligned];
      g->rinv  = & g->r[4*nr_aligned];
      g->memory_owner = (nullptr != g->r);
      return g;
  } // get_memory

  inline void set_derived_grid_quantities(radial_grid_t & g, int const nr) {
      for (int ir = 0; ir < nr; ++ir) {
          auto const r = g.r[ir], dr = g.dr[ir];
          g.rdr[ir]  =   r*dr;
          g.r2dr[ir] = r*r*dr;
          g.rinv[ir] = (ir)? 1./r : 0;
//        std::printf("%g %g\n", r, dr); // DEBUG
      } // ir
  } // set_derived_grid_quantities

  radial_grid_t* create_radial_grid(
        int const npoints
      , float const rmax // [optional] largest radius
      , char const *equation // [optional] how to generate the grid
      , double const anisotropy // [optional] anisotropy parameter for exponential
  ) {

#ifdef  USE_RECIPROCAL_RADIAL_GRID
      auto const mR = 100; // multiplicator for the outer radius
      bool static warn_reciprocal{true}; // NOT thread-safe
      if (nullptr == equation) {
          equation = equation_reciprocal; // modify default behaviour
          if (warn_reciprocal) {
              std::printf("# radial_grid -D USE_RECIPROCAL_RADIAL_GRID with %gx larger nominal outer radius\n", mR*1.);
              warn_reciprocal = false; // switch off
          } // warn_reciprocal
      } // special grid equation requested?
#endif

      int const nr = std::max(std::abs(npoints), 32);
      double const R = std::max(std::abs(rmax)*1., .945);

      int const nr_aligned = align<2>(nr); // padded to multiples of 4
      auto const g = get_memory(nr_aligned);

      double & d = g->anisotropy;

      if (equation_reciprocal == equation) {
          g->equation = equation_reciprocal;
//        auto const n = nr + 4; // with i=nr-1 the outermost radius is i/(n-i)=(n-5)/5 ~=~ n/5
          auto const n = nr + mR/2; // ?
          for (int i = 0; i < nr_aligned; ++i) {
              double const rec = 1./((nr - 1)*(n - i));
              g->r[i]  = (mR*R)*i*rec;
              g->dr[i] = (mR*R)*n*rec*n*rec * (i < nr);
          } // i
          d = n; // store the real number used for the generation of the reciprocal grid in the anisotropy field

      } else if (equation_equidistant == equation) {
          g->equation = equation_equidistant;
          double const dr = rmax/nr;
          for (int ir = 0; ir < nr_aligned; ++ir) {
              g->r[ir]  = ir*dr;
              g->dr[ir] = dr * (ir < nr);
          } // ir

      } else {
          g->equation = equation_exponential;
          d = std::min(std::max(1e-4, anisotropy*1.), .1);
          double const a = R / (std::exp(d*(nr - 1)) - 1.); // prefactor
          for (int ir = 0; ir < nr_aligned; ++ir) {
              double const edi = std::exp(d*ir);
              g->r[ir]  = a*(edi - 1.);
              g->dr[ir] = a*d*edi * (ir < nr);
          } // ir

      } // switch equation

      set_derived_grid_quantities(*g, nr_aligned);
      g->n = nr;
      g->rmax = g->r[g->n - 1]; // implicit conversion to float

      return g;
  } // create_radial_grid

  double get_prefactor(radial_grid_t const & g) {
      if (equation_exponential == g.equation) {
          return g.rmax/(std::exp(g.anisotropy*(g.n - 1)) - 1.);
      } else if (equation_reciprocal == g.equation) {
          int const n = g.anisotropy, i = g.n - 1;
          return (g.r[i]*(n - i))/i;
      } else {
          return g.rmax/std::max(g.n - 1, 1);
      }
  } // get_prefactor

#if 0
  radial_grid_t* create_equidistant_radial_grid(
        int const npoints
      , float const Rmax // =default_Rmax
  ) {
      int const n = std::max(1, npoints); // is the ever called with n < 32?
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
#endif

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
      for (int ir = 0; ir < g.n; ++ir) {
          if (radius < g.r[ir]) return ir - 1;
      } // ir
      return g.n - 2;
  } // find_grid_index


#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo=9) {
      if (echo > 0) std::printf("\n# %s: \n", __func__);
      auto const g = create_radial_grid(1 << 11);
      destroy_radial_grid(g);
      return 0;
  } // test_create_and_destroy

  status_t test_exp_grid(int const echo=3) {
      if (echo > 0) std::printf("\n# %s: \n", __func__);
      int const n = 1 << 11;
      auto & g = *create_radial_grid(n);
      double integ[] = {0, 0, 0};
      if (echo > 3) std::printf("\n## radial grid (%d grid points, anisotropy=%g, up to %g %s, %s):\n"
              "## r, dr, 1/r, integrals over {1,r,r^2}, r, r^2/2, r^3/3\n", n, g.anisotropy, g.rmax*1.0, " Bohr", g.equation);
      for (int ir = 0; ir < g.n; ++ir) {
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
