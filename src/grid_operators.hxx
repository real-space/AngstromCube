#pragma once

typedef int status_t;

#include "grid_operators.hxx"

#include "atom_image.hxx" // ::atom_image_t, ::sho_atom_t
#include "real_space_grid.hxx" // ::grid_t
#include "finite_difference.hxx" // ::finite_difference_t<real_t>
#include "vector_math.hxx" // ::vec<N,T>
#include "boundary_condition.hxx" // Isolated_Boundary

namespace grid_operators {

  template<typename real_t, typename real_fd_t, int D0=1>
  status_t _grid_operator(real_t Hpsi[] // result
                        , real_t const psi[] // input wave functions
                        , real_space_grid::grid_t<D0> const &g // 3D Cartesian grid descriptor
                        , std::vector<atom_image::sho_atom_t> const &a
                        , std::vector<atom_image::atom_image_t> const &ai
                        , int const h0s1 // index controlling which matrix of a[ia] we are multiplying, 0:Hamiltonian or 1:overlap
                        , double const *boundary_phase=nullptr // phase shifts at the boundary [optional]
                        , finite_difference::finite_difference_t<real_fd_t> const *fd=nullptr // finite difference [optional]
                        , double const *potential=nullptr // diagonal potential operator [optional]
                        , int const echo=0) {
      status_t stat = 0;

      if (fd) {
          stat += Laplacian(Hpsi, psi, g, *fd); // prefactor -0.5 moved into fd-coefficients
          if (echo > 8) printf("# %s Apply Laplacian, status=%i\n", __func__, stat);
      } else {
          set(Hpsi, g.all(), real_t(0)); // clear
      } // fd

      size_t const nzyx = g.dim(2) * g.dim(1) * g.dim(0);
      if (echo > 8) printf("# %s Apply %s operator\n", __func__, potential ? "potential" : "unity");
      for(size_t izyx = 0; izyx < nzyx; ++izyx) {
          real_fd_t const V = potential ? potential[izyx] : real_fd_t(1); // apply potential or the unity operation of the overlap operator
          for(int i0 = 0; i0 < D0; ++i0) {
              Hpsi[izyx*D0 + i0] += V * psi[izyx*D0 + i0];
          } // i0 vectorization
      } // izyx

      int constexpr echo_sho = 0;

#if 1
      typedef vector_math::vec<D0,real_fd_t> temp_vec_t; // computation is performed in real_fd_t precision

      if (h0s1 >= 0) {
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
                  auto *const mat = a[ia].get_matrix<real_fd_t>(h0s1);
                  // DGEMM-style matrix multiplication
                  for(int i = 0; i < ncoeff; ++i) {
                      temp_vec_t ci(0.0);
                      for(int j = 0; j < ncoeff; ++j) {
                          auto const am = mat[i*stride + j];
                          temp_vec_t const cj = &image_coeff[j*D0];
                          ci += cj * am;
                      } // j
                      for(int i0 = 0; i0 < D0; ++i0) {
                          V_image_coeff[i*D0 + i0] = ci[i0]; // implicit from real_fd_t to read_t
                      } // i0 vectorization
                  } // i
              } // scope

              if (echo > 8) { 
                  printf("\n# a#%i coefficients and %c*coeff (vector length=%i):\n", atom_id, h0s1+'H', D0);
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
      } // h0s1 >= 0
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
  template<typename real_t, typename real_fd_t, int D0=1>
  status_t grid_Hamiltonian(real_t Hpsi[] // result
                          , real_t const psi[] // input wave functions
                          , real_space_grid::grid_t<D0> const &g // 3D Cartesian grid descriptor
                          , std::vector<atom_image::sho_atom_t> const &a
                          , std::vector<atom_image::atom_image_t> const &ai
                          , finite_difference::finite_difference_t<real_fd_t> const &fd // finite difference
                          , double const potential[] // diagonal potential operator
                          , double const *boundary_phase=nullptr // phase shifts at the boundary [optional]
                          , int const echo=0
  ) { return _grid_operator(Hpsi, psi, g, a, ai, 0, boundary_phase, &fd, potential); }

  // Overlap operator
  template<typename real_t, typename real_fd_t=double, int D0=1>
  status_t grid_Overlapping(real_t Spsi[] // result
                          , real_t const psi[] // input wave functions
                          , real_space_grid::grid_t<D0> const &g // 3D Cartesian grid descriptor
                          , std::vector<atom_image::sho_atom_t> const &a
                          , std::vector<atom_image::atom_image_t> const &ai
                          , double const *boundary_phase=nullptr // phase shifts at the boundary [optional]
                          , int const echo=0
  ) { return _grid_operator<real_t, real_fd_t, D0>(Spsi, psi, g, a, ai, 1, boundary_phase); }
  
  
  //
  // idea: the identity operation of the Overlap operator could be implemented with a
  //       finite difference stencil that has nn[] = {0, -1, -1} (-1 means that there is no pass)
  //       because then, all three operators, Hamiltonian, Overlapping and Preconditioner share
  //       the same structure, whereas the latter has a FD-stencil with all positive coefficients
  //

  template <typename real_t=double, typename real_fd_t=double, int D0=1>
  class grid_operator_t 
  {
    public:

      grid_operator_t(real_space_grid::grid_t<D0> const & rsg, int const natoms=0, int const nprecond=0)
      : grid(rsg), atoms(natoms), has_precond(nprecond > 0), has_overlap(true) { // constructor
//           int const nprecond = control::get("conjugate_gradients.precond", 1.);

          auto const & g = grid;  // abbrev.

          potential = std::vector<double>(g.dim(2)*g.dim(1)*g.dim(0), 0.0); // flat effective local potential, zero everywhere
          
          images = std::vector<atom_image::atom_image_t>(natoms); // atom_images
//           atoms [0] = atom_image::sho_atom_t(3, 0.5, 999); // numax=3, sigma=0.5, atom_id=999
//           images[0] = atom_image::atom_image_t(g.dim(0)*g.h[0]/2, g.dim(1)*g.h[1]/2, g.dim(2)*g.h[2]/2, 999, 0); 
          // image position at the center, index=0 maps into list of sho_atoms
          
          int const bc[] = {Isolated_Boundary, Isolated_Boundary, Isolated_Boundary};
          int const nn[] = {8, 8, 8}; // half-order of finite difference stencil for kinetic energy
          
          int const nn_precond[] = {nprecond, nprecond, nprecond};
          precond = finite_difference::finite_difference_t<real_t>(g.h, bc, nn_precond);
          // this simple preconditioner is a diffusion stencil
          for(int d = 0; d < 3; ++d) {
              precond.c2nd[d][1] = 1/12.;
              precond.c2nd[d][0] = 2/12.; // stencil [1/4 1/2 1/4] in all 3 directions
          } // d
          kinetic = finite_difference::finite_difference_t<real_fd_t>(g.h, bc, nn);
          kinetic.scale_coefficients(-0.5); // prefactor of the kinetic energy in Hartree atomic units

      } // constructor

      status_t Hamiltonian(real_t Hpsi[], real_t const psi[], int const echo=0) const {
          return _grid_operator(Hpsi, psi, grid, atoms, images, 0, boundary_phase.data(), &kinetic, potential.data(), echo);
      } // Hamiltonian

      status_t Overlapping(real_t Spsi[], real_t const psi[], int const echo=0) const {
          return _grid_operator<real_t, real_fd_t, D0>(Spsi, psi, grid, atoms, images, 1, boundary_phase.data(), nullptr, nullptr, echo);
      } // Overlapping

      status_t Conditioner(real_t Cpsi[], real_t const psi[], int const echo=0) const {
          return _grid_operator(Cpsi, psi, grid, atoms, images, -1, boundary_phase.data(), &precond, nullptr, echo);
      } // Pre-Conditioner
      
    
    private:
      real_space_grid::grid_t<D0> grid;
      std::vector<atom_image::sho_atom_t> atoms;
      std::vector<atom_image::atom_image_t> images;
      std::vector<double> boundary_phase; // could be real_t or real_fd_t in the future
      std::vector<double> potential;
      finite_difference::finite_difference_t<real_fd_t> kinetic;
      finite_difference::finite_difference_t<real_t> precond;
      bool has_precond;
      bool has_overlap;
    public:
      real_space_grid::grid_t<D0> const & get_grid() const { return grid; }
      bool use_precond() const { return has_precond; }
      bool use_overlap() const { return has_overlap; }
  }; // class grid_operator_t
  
  status_t all_tests(int const echo=0);

} // namespace grid_operators
