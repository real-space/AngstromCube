
#include <stdlib.h> // abs
#include <cmath> // fabs, exp
#include <cstdio> // printf

#include "radial_grid.hxx"

#include "radial_grid.h" // radial_grid_t
#include "min_and_max.h" // min, max

typedef int status_t;

namespace radial_grid {

  radial_grid_t* create_exponential_radial_grid(int const npoints,
                   float const rmax, float const anisotropy) // optional args
  {
    double const R = max(fabs(rmax), .945);
    double const d = min(max(1e-4, anisotropy), 0.1);
    int const N = max(abs(npoints), 32);

    int const N_aligned = (((N - 1) >> 2) + 1) << 2; // padded to 4

    auto g = new radial_grid_t;
    g->r = new double[5*N_aligned];
    g->dr     = &g->r[1*N_aligned];
    g->rdr    = &g->r[2*N_aligned];
    g->r2dr   = &g->r[3*N_aligned];
    g->rinv   = &g->r[4*N_aligned];

    double const a = R / (exp(d*(N - 1)) - 1.); // prefactor
    for (auto ir = 0; ir < N_aligned; ++ir) {
      double const edi = exp(d*ir);
      double const r = a*(edi - 1.);
      double const dr = a*d*edi * (ir < N);
      g->r[ir] = r;
      g->dr[ir] = dr;
      g->rdr[ir] = r*dr;
      g->r2dr[ir] = r*r*dr;
      g->rinv[ir] = (ir)? 1./r : 0;
    } // ir

    g->n = N;
    g->rmax = g->r[g->n - 1];
    return g;
  } // create_exponential_radial_grid
  
  void destroy_radial_grid(radial_grid_t* g) {
    delete [] g->r;
    g->n = 0;
  } // destroy

#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS
  status_t all_tests() { printf("\nWarning: %s no unit test implemented!\n\n", __FILE__); return -1; }
#endif // NO_UNIT_TESTS
  
} // namespace radial_grid
