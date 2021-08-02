
#include <cmath> // std::abs, std::exp
#include <cstdio> // std::printf
#include <algorithm> // std::min, std::max

#include "radial_grid.hxx" // radial_grid_t

#include "inline_math.hxx" // align

namespace radial_grid {
  
  char const* get_formula(char const equation) {
      if (equation_equidistant == equation) return "r=a*i/n";
      auto const eq_reciprocal  = "r=a*i/(n-i)";
      auto const eq_exponential = "r=a*(exp(d*i)-1)";
      if (equation_reciprocal  == equation) return eq_reciprocal;
      if (equation_exponential == equation) return eq_exponential;
#ifdef  USE_RECIPROCAL_RADIAL_GRID
      return eq_reciprocal;
#else
      return eq_exponential;
#endif
  } // get_formula

  radial_grid_t* get_memory(size_t const nr_aligned) {
      auto g = new radial_grid_t;
      g->r = new double[5*nr_aligned];
      g->dr    = & g->r[1*nr_aligned];
      g->rdr   = & g->r[2*nr_aligned];
      g->r2dr  = & g->r[3*nr_aligned];
      g->rinv  = & g->r[4*nr_aligned];
      g->memory_owner = (nullptr != g->r);
      return g;
  } // get_memory

  void set_derived_grid_quantities(radial_grid_t & g, int const nr) {
      for (int ir = 0; ir < nr; ++ir) {
          auto const r = g.r[ir], dr = g.dr[ir];
          g.rdr[ir]  =   r*dr;
          g.r2dr[ir] = r*r*dr;
          g.rinv[ir] = (ir)? 1./r : 0;
#ifdef FULL_DEBUG
          std::printf("%g %g\n", r, dr);
#endif // FULL_DEBUG
      } // ir
  } // set_derived_grid_quantities

  radial_grid_t* create_radial_grid(
        int const npoints
      , float const rmax // [optional] largest radius
      , char equation // [optional] how to generate the grid
      , double const anisotropy // [optional] anisotropy parameter for exponential
  ) {

      
      auto const mR = 128; // multiplicator for the outer radius of reciprocal grids
#ifdef  USE_RECIPROCAL_RADIAL_GRID
      bool static warn_reciprocal{true}; // NOT thread-safe
      if ('\0' == equation) {
          equation = equation_reciprocal; // modify default behaviour
          if (warn_reciprocal) {
              std::printf("# radial_grid -D USE_RECIPROCAL_RADIAL_GRID with %gx larger nominal outer radius\n", mR*1.);
              warn_reciprocal = false; // switch off
          } // warn_reciprocal
      } // special grid equation requested?
#endif // USE_RECIPROCAL_RADIAL_GRID

      int const nr = std::max(std::abs(npoints), 32);
      auto const R = std::max(std::abs(rmax)*1., .945);

      int const nr_aligned = align<2>(nr); // padded to multiples of 4
      auto const g = get_memory(nr_aligned);

      double & d = g->anisotropy;

//    std::printf("# create a mesh with R= %g Bohr, n=%d, formula=%s\n", R, nr, get_formula(equation));
      if (equation_reciprocal == equation) {
          g->equation = equation_reciprocal;
          auto const n = nr + mR/2; // with i=nr-1 the outermost radius is i/(n-i)=(n-1-mR/2)/(mR/2 + 1)
          for (int i = 0; i < nr_aligned; ++i) {
              double const rec = 1./((nr - 1)*(n - i));
              g->r[i]  = (mR*R)*i*rec;
              g->dr[i] = (mR*R)*n*rec*n*rec * (i < nr);
          } // i
          d = n; // store the real number used for the generation of the reciprocal grid in the anisotropy field

      } else if (equation_equidistant == equation) {
          g->equation = equation_equidistant;
          double const dr = R/nr;
//        std::printf("# create an equidistant mesh with dr=%g, R= %g Bohr, n=%d\n", dr, R, nr);
          for (int ir = 0; ir < nr_aligned; ++ir) {
              g->r[ir]  = ir*dr;
              g->dr[ir] = dr * (ir < nr);
          } // ir
          d = 0; // no anisotropy

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
  } // destroy_radial_grid

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
      auto const gp = create_radial_grid(1 << 11);
      destroy_radial_grid(gp);
      return 0;
  } // test_create_and_destroy

  status_t test_radial_grid_integral(int const echo=3) {
      if (echo > 0) std::printf("\n# %s: \n", __func__);
      int const n = 1 << 11;
      auto & g = *create_radial_grid(n);
      double integ[] = {0, 0, 0};
      if (echo > 3) std::printf("\n## radial grid (%d grid points, anisotropy= %g, up to %g %s, %s):\n"
              "## r, dr, 1/r, integral {1,r,r^2} dr, reference {r, r^2/2, r^3/3}\n",
              n, g.anisotropy, g.rmax*1.0, "Bohr", get_formula(g.equation));
      for (int ir = 0; ir < g.n; ++ir) {
          if (echo > 3 && (ir < 3 || ir > g.n - 4 || echo > 11)) {
              auto const r = g.r[ir];
              std::printf("%g %g %g %g %g %g %g %g %g\n", r, g.dr[ir], g.rinv[ir],
                  integ[0], integ[1], integ[2], r, r*r/2, r*r*r/3);
          } // echo
          integ[0] += g.dr[ir];
          integ[1] += g.rdr[ir];
          integ[2] += g.r2dr[ir];
      } // ir
      double const dev = std::abs(integ[0] - g.rmax);
      if (echo > 3) std::printf("# %s: integral dr from 0 to R deviates %g from R= %g %s\n",
                                   __func__, dev*1.0, g.rmax*1.0, "Bohr");
      destroy_radial_grid(&g);
      return (dev > .05); // integral dr up to rmax should match rmax
  } // test_radial_grid_integral

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_create_and_destroy(echo);
      stat += test_radial_grid_integral(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace radial_grid
