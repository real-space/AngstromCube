#include <cstdio> // printf
#include <cassert> // assert
#include <numeric> // std::iota
#include <cmath> // std::floor
#include <vector> // std::vector

#include "density_generator.hxx"

#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set, pow2
#include "inline_tools.hxx" // align<n>
#include "constants.hxx" // ::sqrtpi, ::pi
#include "real_space_grid.hxx" // ::grid_t, ::add_function
#include "boundary_condition.hxx" // ::periodic_images
#include "simple_timer.hxx" // SimpleTimer
#include "control.hxx" // control::get
#include "recorded_warnings.hxx" // warn, error
#include "data_view.hxx" // view3D<T>

// #define FULL_DEBUG
// #define DEBUG

namespace density_generator {
  
  double print_stats(double const values[], size_t const all, double const dV=1, char const prefix=' ') {
      double gmin{9e307}, gmax{-gmin}, gsum{0}, gsum2{0};
      for(size_t i = 0; i < all; ++i) {
          gmin = std::min(gmin, values[i]);
          gmax = std::max(gmax, values[i]);
          gsum  += values[i];
          gsum2 += pow2(values[i]);
      } // i
      printf("%c grid stats min %g max %g integral %g avg %g\n", prefix, gmin, gmax, gsum*dV, gsum/all);
      return gsum*dV;
  } // print_stats // ToDo: move somewhere, where potential_generator and this module can access it

  template<typename real_t, int const D0=1>
  status_t density(real_t const eigenfunctions[]
      , real_space_grid::grid_t<D0> const & g
      , int const nbands=1, int const nkpoints=1
      , int const echo=0) {
      // SimpleTimer init_function_timer(__FILE__, __LINE__, __func__, echo);
      status_t stat{0};

      if (nullptr == eigenfunctions) { warn("eigenfunctions received are nullptr"); return -1; }
      
      std::vector<double> rho_valence(g.all(), 0.0);
      if (echo > 3) printf("# %s assume %d bands and %d k-points\n", __func__, nbands, nkpoints);
      if (echo > 3) printf("# %s assume eigenfunctions on a %d x %d x %d Cartesian grid\n", 
                              __func__, g.dim('x'), g.dim('y'), g.dim('z'));
      
      view3D<real_t const> const psi(eigenfunctions, nbands, g.all()); // wrap
      
// #pragma omp parallel
      for(int ikpoint = 0; ikpoint < nkpoints; ++ikpoint) {
          double const kpoint_weight = 1; // depends on ikpoint
          auto const psi_k = psi[ikpoint];
          
          for(int iband = 0; iband < nbands; ++iband) {
              double const band_occupation = 1; // depends on iband and ikpoint
              double const weight_nk = band_occupation * kpoint_weight;
              auto const psi_nk = psi_k[iband];
              
              if (weight_nk > 0) {
// #pragma omp for
                  for(size_t izyx = 0; izyx < g.all(); ++izyx) { // parallel
                      rho_valence[izyx] += weight_nk * pow2(psi_nk[izyx]);
                  } // izyx
              } // weight is positive
          } // iband
      } // ikpoint

      if (echo > 1) {
          printf("\n# Total valence density  grid stats:");
          print_stats(rho_valence.data(), g.all(), g.dV());
      } // echo

      return stat;
  } // density
 
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_init(int const echo=3) {
      int const dims[] = {4, 5, 6};
      real_space_grid::grid_t<1> g(dims);
      std::vector<float> wave(g.all());
      std::iota(wave.begin(), wave.end(), 0);
      return density<float,1>(wave.data(), g, 1, 1, echo);
  } // test_init

  status_t all_tests(int const echo) {
    auto status{0};
    status += test_init(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace density_generator
