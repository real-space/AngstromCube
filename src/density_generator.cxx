#include <cstdio> // printf
#include <cassert> // assert
#include <numeric> // std::iota
#include <cmath> // std::floor
#include <vector> // std::vector<T>

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
#include "grid_operators.hxx" // ::get_atom_coeffs

// #define FULL_DEBUG
// #define DEBUG

namespace density_generator {

  // ToDo: move somewhere, where potential_generator and this module can access it  
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
  } // print_stats

  template<typename real_t, typename real_fd_t=double, int const D0=1>
  status_t density(real_t const eigenfunctions[]
      , grid_operators::grid_operator_t<real_t,real_fd_t,D0> const & op
      , int const nbands=1, int const nkpoints=1
      , int const echo=0) {
      // SimpleTimer init_function_timer(__FILE__, __LINE__, __func__, echo);
      status_t stat{0};

      if (nullptr == eigenfunctions) { warn("eigenfunctions received are nullptr"); return -1; }

      auto const g = op.get_grid();
      
      std::vector<double> rho_valence(g.all(), 0.0);
      if (echo > 3) printf("# %s assume %d bands and %d k-points\n", __func__, nbands, nkpoints);
      if (echo > 3) printf("# %s assume eigenfunctions on a %d x %d x %d Cartesian grid\n", 
                              __func__, g.dim('x'), g.dim('y'), g.dim('z'));
      
      int const na = op.get_natoms();
      std::vector<double*> atom_coeff(na, nullptr);
      std::vector<double*> atom_rho(na, nullptr); // atomic density matrix
      for(int ia = 0; ia < na; ++ia) {
          int const numax = op.get_numax(ia);
          int const ncoeff = sho_tools::nSHO(numax);
          atom_coeff[ia] = new double[ncoeff*D0];
          atom_rho[ia] = new double[ncoeff*ncoeff];
          set(atom_rho[ia], ncoeff*ncoeff, 0.0); // init clear
      } // ia
      assert(1 == D0); // vectorization is not implemented in all parts

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
                  
                  // construct the atomic density matrix atom_rho
                  op.get_atom_coeffs(atom_coeff.data(), psi_nk, echo/2); // project again
                  for(int ia = 0; ia < na; ++ia) {
                      int const numax = op.get_numax(ia);
                      int const ncoeff = sho_tools::nSHO(numax);
                      // add to the atomic density matrix, ToDo
                      for(int i = 0; i < ncoeff; ++i) {
                          auto const c_i = atom_coeff[ia][i*D0 + 0];
#ifdef DEVEL
                          if (echo > 9) printf("# kpoint #%i band #%i atom #%i coeff[%i] = %g\n", ikpoint, iband, ia, i, c_i);
#endif // DEVEL
                          for(int j = 0; j < ncoeff; ++j) {
                              auto const c_j = atom_coeff[ia][j*D0 + 0];
                              atom_rho[ia][i*ncoeff + j] += weight_nk * c_i * c_j; 
                          } // j
                      } // i
                  } // ia
                  
              } // weight is positive
          } // iband
      } // ikpoint

      if (echo > 1) {
          printf("\n# Total valence density  grid stats:");
          print_stats(rho_valence.data(), g.all(), g.dV());
      } // echo

      // memory cleanup
      for(int ia = 0; ia < na; ++ia) {
#ifdef DEVEL
          if (echo > 6) {
              int ncoeff = sho_tools::nSHO(op.get_numax(ia));
              printf("\n# show %d x %d density matrix for atom #%i", ncoeff, ncoeff, ia);
              for(int i = 0; i < ncoeff; ++i) {
                  printf("\n# ");
                  for(int j = 0; j < ncoeff; ++j) {
                      printf(" %.9g", atom_rho[ia][i*ncoeff + j]);
                  } // j
              }   printf("\n");
          } // echo
#endif // DEVEL
          delete[] atom_coeff[ia];
          delete[] atom_rho[ia];
      } // ia

      return stat;
  } // density
 
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_init(int const echo=3) {
      real_space_grid::grid_t<1> g(4, 5, 6);
      grid_operators::grid_operator_t<float,double,1> op(g, 1);
      std::vector<float> wave(g.all());
      std::iota(wave.begin(), wave.end(), 0);
      return density(wave.data(), op, 1, 1, echo);
  } // test_init

  status_t all_tests(int const echo) {
    auto status{0};
    status += test_init(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS

} // namespace density_generator
