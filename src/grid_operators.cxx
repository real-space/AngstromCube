#include <cstdio> // printf
#include <vector> // std::vector<T>
#include <cassert> // assert
 
#include "grid_operators.hxx"

#include "real_space_grid.hxx" // ::grid_t
#include "finite_difference.hxx" // ::finite_difference_t, ::Laplacian
#include "atom_image.hxx" // ::atom_image_t, ::sho_atom_t
#include "inline_math.hxx" // set
#include "sho_projection.hxx" // ::sho_project, ::sho_add
#include "sho_tools.hxx" // ::nSHO

// #include "display_units.h" // eV, _eV, Ang, _Ang

// #define FULL_DEBUG
#define DEBUG

#ifdef  DEBUG
    #include "debug_output.hxx" // dump_to_file
#endif

#ifdef FULL_DEBUG
    #define full_debug(print) print 
#else
    #define full_debug(print)
#endif

#ifdef DEBUG
    #define debug(print) print 
#else
    #define debug(print)
#endif


namespace grid_operators {
  // setup of the real-space grid-based Hamiltonian and overlap operator
  
  template<typename real_t, int D0>
  status_t _grid_operator(real_t Hpsi[] // result
                        , real_t const psi[] // input wave functions
                        , real_space_grid::grid_t<D0> const &g // 3D Cartesian grid descriptor
                        , std::vector<atom_image::sho_atom_t> const &a
                        , std::vector<atom_image::atom_image_t> const &ai
                        , int const h0s1 // index controlling which matrix of a[ia] we are multiplying, 0:Hamiltonian or 1:overlap
                        , double const *boundary_phase=nullptr // phase shifts at the boundary [optional]
                        , finite_difference::finite_difference_t<real_t> const *fd=nullptr // finite difference [optional]
                        , double const *potential=nullptr // diagonal potential operator [optional]
                        ) {
      status_t stat = 0;
      int constexpr echo = 0;

      if (fd) {
          stat += Laplacian(Hpsi, psi, g, *fd, -0.5); // -0.5: kinetic energy prefactor in Hartree units
          if (echo > 8) printf("# %s Apply Laplacian, status=%i\n", __func__, stat);
      } else {
          set(Hpsi, g.all(), real_t(0)); // clear
      } // fd

      size_t const nzyx = g.dim(2) * g.dim(1) * g.dim(0);
      if (echo > 8) printf("# %s Apply %s operator\n", __func__, potential ? "potential" : "unity");
      for(size_t izyx = 0; izyx < nzyx; ++izyx) {
          real_t const V = potential ? potential[izyx] : real_t(1); // apply potential or the unity operation of the overlap operator
          for(int i0 = 0; i0 < D0; ++i0) {
              Hpsi[izyx*D0 + i0] += V * psi[izyx*D0 + i0];
          } // i0 vectorization
      } // izyx

      int constexpr echo_sho = 0;

#if 1
      int const na = a.size();
      int const nai = ai.size();
      assert(nai == na); // this simple version works only of there is one image per atom
      for(int iai = 0; iai < nai; ++iai) { // warning: race condition on Hpsi if we parallelize this loop
          int const ia = ai[iai].index();
          assert(ia >= 0); assert(ia < na);
          auto const atom_id = ai[iai].atom_id();
          assert(atom_id == a[ia].atom_id()); // make sure the lists are linked correctly

          int const numax = a[ia].numax();
          int const ncoeff = sho_tools::nSHO(numax);
          
          std::vector<real_t>   image_coeff(ncoeff*D0, 0.0); // can we omit the initialization here?
          std::vector<real_t> V_image_coeff(ncoeff*D0, 0.0);

          stat += sho_projection::sho_project(image_coeff.data(), numax, ai[iai].get_pos(), a[ia].sigma(), psi, g, echo_sho);
          // warning: no k-dependent Bloch-factors in this version
          if (echo > 7) printf("# %s Apply non-local projection for a#%i, status=%i\n", __func__, atom_id, stat);

          { // scope: V_image_coeff = matrix * image_coeff
              int const stride = a[ia].stride();
              assert(stride >= ncoeff); // check internal consistency
              double const * const mat = a[ia].get_matrix<real_t>(h0s1);
              // DGEMM-style matrix multiplication
              for(int i = 0; i < ncoeff; ++i) {
                  for(int j = 0; j < ncoeff; ++j) {
                      auto const am = mat[i*stride + j];
                      for(int i0 = 0; i0 < D0; ++i0) {
                          V_image_coeff[i*D0 + i0] += am * image_coeff[j*D0 + i0];
                      } // i0 vectorization
                  } // j
              } // i
          } // scope
          
          if (echo > 8) { 
              printf("\n# a#%i coefficients and H*coeff (vector length=%i):\n", atom_id, D0);
              for(int i = 0; i < ncoeff; ++i) {
                  printf("# i=%i", i);
                  for(int i0 = 0; i0 < D0; ++i0) {
                      printf("\t%g %g", image_coeff[i*D0 + i0], V_image_coeff[i*D0 + i0]);
                  } // i0 vectorization
                  printf("\n");
              } // i
          } // echo

          // warning: no k-dependent inverse Bloch-factors in this version
          stat += sho_projection::sho_add(Hpsi, g, V_image_coeff.data(), numax, ai[iai].get_pos(), a[ia].sigma(), echo_sho);
          if (echo > 7) printf("# %s Apply non-local %c addition for a#%i, status=%i\n", __func__, 'H'+h0s1, atom_id, stat);
      } // iai
#else
      int const na = a.size();
      std::vector<double*> atom_coeff(na);
      std::vector<double*> V_atom_coeff(na);
      for(int ia = 0; ia < na; ++ia) {
          int const ncoeff = sho_tools::nSHO(a[ia].numax);
          atom_coeff[ia]   = new double[ncoeff*D0];
          V_atom_coeff[ia] = new double[ncoeff*D0];
          set(atom_coeff[ia], ncoeff*D0, 0.0);
      } // ia

      int const nai = ai.size();
      for(int iai = 0; iai < nai; ++iai) {
          int const ia = ai[iai].index(); 
          assert(ia >= 0); assert(ia < na);
          assert(ai[iai].atom_id == a[ia].atom_id); // make sure the lists are linked correctly
          int const numax = a[ia].numax;
          int const ncoeff = sho_tools::nSHO(numax);
          std::vector<real_t> image_coeff(ncoeff*D0, 0.0); // can we omit the initialization?

          stat += sho_projection::sho_project(image_coeff, numax, ai[iai].get_pos(), a[ia].sigma, psi, g, echo);
          
          int const ncoeff = sho_tools::nSHO(numax);
          double const Block_factor = 1.0; // ToDo: use k-dependent Bloch-factors
          add_product(atom_coeff[ia], ncoeff*D0, image_coeff, Block_factor);
          
      } // iai

      for(int ia = 0; ia < na; ++ia) {
          int const numax = a[ia].numax;
          int const stride = a[ia].stride;
          int const ncoeff = sho_tools::nSHO(numax);
          assert( nullptr != a[ia].matrix );
          
          V_atom_coeff = matrix * atom_coeff;
          
      } // ia
     
      for(int iai = 0; iai < nai; ++iai) { // warning: race condition on Hpsi if we parallelize this loop
          int const ia = ai[iai].index;
          int const numax = a[ia].numax;
          int const ncoeff = sho_tools::nSHO(numax);
          double const inv_Block_factor = 1.0; // ToDo: use k-dependent inverse Bloch-factors
          std::vector<real_t> V_image_coeff(ncoeff*D0);

          product(V_image_coeff, ncoeff*D0, V_atom_coeff[ia], inv_Block_factor);

          stat += sho_projection::sho_add(Hpsi, g, V_image_coeff, a[ia].numax, ai[iai].get_pos(), a[ia].sigma, echo);
      } // iai

      for(int ia = 0; ia < na; ++ia) { delete[] atom_coeff[ia]; delete[] V_atom_coeff[ia]; }
#endif
      return stat;
  } // _grid_operator

  // Hamiltonian
  template<typename real_t, int D0>
  status_t grid_Hamiltonian(real_t Hpsi[] // result
                          , real_t const psi[] // input wave functions
                          , real_space_grid::grid_t<D0> const &g // 3D Cartesian grid descriptor
                          , std::vector<atom_image::sho_atom_t> const &a
                          , std::vector<atom_image::atom_image_t> const &ai
                          , finite_difference::finite_difference_t<real_t> const &fd // finite difference
                          , double const potential[] // diagonal potential operator
                          , double const *boundary_phase // phase shifts at the boundary [optional]
  ) { return _grid_operator(Hpsi, psi, g, a, ai, 0, boundary_phase, &fd, potential); }

  // Overlap operator
  template<typename real_t, int D0>
  status_t grid_Overlapping(real_t Spsi[] // result
                          , real_t const psi[] // input wave functions
                          , real_space_grid::grid_t<D0> const &g // 3D Cartesian grid descriptor
                          , std::vector<atom_image::sho_atom_t> const &a
                          , std::vector<atom_image::atom_image_t> const &ai
                          , double const *boundary_phase // phase shifts at the boundary [optional]
  ) { return _grid_operator(Spsi, psi, g, a, ai, 1, boundary_phase); }
  //
  // idea: the identity operation of the Overlap operator could be implemented with a
  //       finite difference stencil that has nn[] = {0, -1, -1} (-1 means that there is no pass)
  //       because then, all three operators, Hamiltonian, Overlapping and Preconditioner share
  //       the same structure, whereas latter has a FD-stencil with all positive coefficients

  template // explicit template instantiation for real_t=double and D0=1
  status_t grid_Hamiltonian(double Hpsi[] // result
                          , double const psi[] // input wave functions
                          , real_space_grid::grid_t<1> const &g // 3D Cartesian grid descriptor
                          , std::vector<atom_image::sho_atom_t> const &a
                          , std::vector<atom_image::atom_image_t> const &ai
                          , finite_difference::finite_difference_t<double> const &fd // finite difference
                          , double const potential[] // diagonal potential operator
                          , double const *boundary_phase);

  template // explicit template instantiation for real_t=double and D0=1
  status_t grid_Overlapping(double Spsi[] // result
                          , double const psi[] // input wave functions
                          , real_space_grid::grid_t<1> const &g // 3D Cartesian grid descriptor
                          , std::vector<atom_image::sho_atom_t> const &a
                          , std::vector<atom_image::atom_image_t> const &ai
                          , double const *boundary_phase);
  
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t basic_test(int const echo=9) {
      status_t stat = 0;
      int constexpr D0 = 2; // vectorization
      int const dims[] = {12, 13, 14}; int const bc[] = {0, 0, 0}, nn[] = {8, 8, 8};
      real_space_grid::grid_t<D0> g(dims);
      std::vector<double> psi(2*g.all(), 1.0);
      std::vector<double> potential(dims[2]*dims[1]*dims[0], 0.5);
      std::vector<atom_image::sho_atom_t> a(1);
      std::vector<atom_image::atom_image_t> ai(1);
      a[0]  = atom_image::sho_atom_t(3, 0.5, 999); // numax=3, sigma=0.5, atom_id=999
      ai[0] = atom_image::atom_image_t(dims[0]*g.h[0]/2, dims[1]*g.h[1]/2, dims[2]*g.h[2]/2, 999, 0);
      finite_difference::finite_difference_t<double> kinetic(g.h, bc, nn);
      stat += grid_Hamiltonian(psi.data(), &psi[g.all()], g, a, ai, kinetic, potential.data());
      stat += grid_Overlapping(psi.data(), &psi[g.all()], g, a, ai);
      return stat;
  } // basic_test

  status_t all_tests(int const echo) {
    auto status = 0;
    status += basic_test(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace grid_operators
