
#include <stdlib.h> // abs
#include <cmath> // fabs, exp
#include <cstdio> // printf
#include <algorithm> // min, max

#include "radial_grid.hxx"

#include "radial_grid.h" // radial_grid_t

typedef int status_t;

namespace radial_grid {

  radial_grid_t* create_exponential_radial_grid(int const npoints,
                   float const rmax, float const anisotropy) // optional args
  {
    double const R = std::max(std::abs(rmax)*1., .945);
    double const d = std::min(std::max(1e-4, anisotropy*1.), .1);
    int const N = std::max(std::abs(npoints), 32);

    int const N_aligned = (((N - 1) >> 2) + 1) << 2; // padded to 4

    auto g = new radial_grid_t;
    g->r = new double[5*N_aligned];
    g->dr     = &g->r[1*N_aligned];
    g->rdr    = &g->r[2*N_aligned];
    g->r2dr   = &g->r[3*N_aligned];
    g->rinv   = &g->r[4*N_aligned];
    g->memory_owner = (nullptr != g->r);

    double const a = R / (std::exp(d*(N - 1)) - 1.); // prefactor
    for (auto ir = 0; ir < N_aligned; ++ir) {
      double const edi = std::exp(d*ir);
      double const r = a*(edi - 1.);
      double const dr = a*d*edi * (ir < N);
      g->r[ir] = r;
      g->dr[ir] = dr;
      g->rdr[ir] = r*dr;
      g->r2dr[ir] = r*r*dr;
      g->rinv[ir] = (ir)? 1./r : 0;
//       printf("%g %g\n", r, dr); // DEBUG
    } // ir

    g->n = N;
    g->rmax = g->r[g->n - 1];
    return g;
  } // create_exponential_radial_grid

  radial_grid_t* create_pseudo_radial_grid(radial_grid_t const &tru, double const r_min, int const echo)
  {
      // find a suitable grid point to start from
      int ir = 0; while (tru.r[ir] < r_min) ++ir;
      if (echo > 3) printf("# start pseudo grid from r[%d]=%g Bohr\n", ir, tru.r[ir]);
      int const ir_offset = ir;

      auto g = new radial_grid_t;
      // offset pointers
      g->r      = tru.r    + ir_offset;
      g->dr     = tru.rdr  + ir_offset;
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

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  int test(int const echo=9) {
    printf("\n# %s: \n", __func__);
    auto g = create_exponential_radial_grid(1 << 11);
    destroy_radial_grid(g);
    return 0;
  } // test

  status_t all_tests() {
    auto status = 0;
    status += test();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS
  
} // namespace radial_grid
