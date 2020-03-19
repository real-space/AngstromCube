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
                        , int const echo=0
                        , real_t** atomic_projection_coefficients=nullptr
                         ) {
      status_t stat = 0;

      if (Hpsi) {
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
      } // Hpsi != nullptr

      int constexpr echo_sho = 0;

      typedef vector_math::vec<D0,real_fd_t> temp_vec_t; // computation is performed in real_fd_t precision

      if (h0s1 >= 0) {
#if 0
          int const na = a.size(); // number of atoms
          int const nai = ai.size(); // number of atomic images
          assert(nai == na); // this simple version works only of there is one image per atom
          for(int iai = 0; iai < nai; ++iai) { // warning: race condition on Hpsi if we parallelize this loop
              int const ia = ai[iai].index();
              assert(ia >= 0); assert(ia < na);
              auto const atom_id = ai[iai].atom_id();
              assert(atom_id == a[ia].atom_id()); // make sure the lists are linked correctly

              int const numax = a[ia].numax();
              int const ncoeff = sho_tools::nSHO(numax);
              
              std::vector<real_t>   image_coeff(ncoeff*D0, 0.0); // can we omit the initialization here?
              std::vector<real_t> V_image_coeff(ncoeff*D0, 0.0); // operator applied to coefficients

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
                          V_image_coeff[i*D0 + i0] = ci[i0]; // implicit from real_fd_t to real_t
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
              if (Hpsi) {
                  stat += sho_projection::sho_add(Hpsi, g, V_image_coeff.data(), numax, ai[iai].get_pos(), a[ia].sigma(), echo_sho);
                  if (echo > 7) printf("# %s Apply non-local %c addition for a#%i, status=%i\n", __func__, 'H'+h0s1, atom_id, stat);
              } // Hpsi != nullptr
          } // iai
#else
          int const na = a.size(); // number of atoms
          std::vector<std::vector<real_t>> V_atom_coeff(na);
          std::vector<std::vector<real_t>>   atom_coeff(na);
          for(int ia = 0; ia < na; ++ia) {
              int const numax = a[ia].numax();
              int const ncoeff = sho_tools::nSHO(numax);
              V_atom_coeff[ia] = std::vector<real_t>(ncoeff*D0, 0.0);
                atom_coeff[ia] = std::vector<real_t>(ncoeff*D0, 0.0);
          } // ia

          int const nai = ai.size(); // number of atomic images
          bool const need_images = (nai > na); // if there are more periodic images than atoms, we need the Bloch phases
          for(int iai = 0; iai < nai; ++iai) {
              int const ia = ai[iai].index(); 
              assert(ia >= 0); assert(ia < na);
              assert(ai[iai].atom_id() == a[ia].atom_id()); // make sure the lists are linked correctly
              int const numax = a[ia].numax();
              int const ncoeff = sho_tools::nSHO(numax);
              if (need_images) {
                  std::vector<real_t> image_coeff(ncoeff*D0, 0.0); // can we omit the initialization?
              
                  stat += sho_projection::sho_project(image_coeff.data(), numax, ai[iai].get_pos(), a[ia].sigma(), psi, g, echo);
              
                  real_t const Bloch_factor = 1.0; // ToDo: use k-dependent Bloch-factors
                  add_product(atom_coeff[ia].data(), ncoeff*D0, image_coeff.data(), Bloch_factor);
              } else {
                  stat += sho_projection::sho_project(atom_coeff[ia].data(), numax, ai[iai].get_pos(), a[ia].sigma(), psi, g, echo);
              } // need_images

          } // iai

          if (atomic_projection_coefficients) {
              for(int ia = 0; ia < na; ++ia) {
                  int const numax = a[ia].numax();
                  int const ncoeff = sho_tools::nSHO(numax);
                  set(atomic_projection_coefficients[ia], ncoeff*D0, atom_coeff[ia].data()); // export
              } // ia
          } // atomic_projection_coefficients
                  
          if (Hpsi) {

              for(int ia = 0; ia < na; ++ia) {
                  int const numax = a[ia].numax();
                  int const ncoeff = sho_tools::nSHO(numax);
              
                  { // scope: V_atom_coeff = matrix * atom_coeff
                      int const stride = a[ia].stride();
                      assert(stride >= ncoeff); // check internal consistency
                      auto *const mat = a[ia].get_matrix<real_fd_t>(h0s1);
                      // DGEMM-style matrix multiplication
                      for(int i = 0; i < ncoeff; ++i) {
                          temp_vec_t ci(0.0);
                          for(int j = 0; j < ncoeff; ++j) {
                              auto const am = mat[i*stride + j];
                              
                              temp_vec_t const cj = &atom_coeff[ia][j*D0];
                              ci += cj * am;
                          } // j
                          for(int i0 = 0; i0 < D0; ++i0) {
                              V_atom_coeff[ia][i*D0 + i0] = ci[i0]; // implicit from real_fd_t to real_t
                          } // i0 vectorization
                      } // i
                  } // scope
              
              } // ia
            
              for(int iai = 0; iai < nai; ++iai) { // warning: race condition on Hpsi if we parallelize this loop
                  int const ia = ai[iai].index();
                  int const numax = a[ia].numax();
                  int const ncoeff = sho_tools::nSHO(numax);
                  auto add_coeff = V_atom_coeff[ia].data();
                  if (need_images) {
                      real_t const inv_Bloch_factor = 1.0; // ToDo: use k-dependent inverse Bloch-factors
                      std::vector<real_t> V_image_coeff(ncoeff*D0);
                      add_coeff = V_image_coeff.data();
                      scale(add_coeff, ncoeff*D0, V_atom_coeff[ia].data(), inv_Bloch_factor);
                  } // need images
                  stat += sho_projection::sho_add(Hpsi, g, add_coeff, numax, ai[iai].get_pos(), a[ia].sigma(), echo);
              } // iai
              
          } // Hpsi != nullptr

#endif
      } // h0s1 >= 0
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

  // Overlap operator
  template<typename real_t, int D0=1>
  status_t get_atom_coeffs(real_t** atom_coeffs // result
                          , real_t const psi[] // input wave functions
                          , real_space_grid::grid_t<D0> const &g // 3D Cartesian grid descriptor
                          , std::vector<atom_image::sho_atom_t> const &a
                          , std::vector<atom_image::atom_image_t> const &ai
                          , double const *boundary_phase=nullptr // phase shifts at the boundary [optional]
                          , int const echo=0
  ) { return _grid_operator<real_t, real_t, D0>(nullptr, psi, g, a, ai, 0, boundary_phase, nullptr, nullptr, 0, atom_coeffs); }
  
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
          for(int ia = 0; ia < natoms; ++ia) {
              int const global_atom_id = ia;
              int const index = ia;
              int const numax = 1; // s and p only (for test)
              double const sigma = 1.;
              atoms[ia]  = atom_image::sho_atom_t(sigma, numax, global_atom_id);
              images[ia] = atom_image::atom_image_t(0,0,0, global_atom_id, index);
          } // ia
          
          int const nn[] = {8, 8, 8}; // half-order of finite difference stencil for kinetic energy
          
          // this simple preconditioner is a diffusion stencil
          int const nn_precond[] = {nprecond, nprecond, nprecond};
          precond = finite_difference::finite_difference_t<real_t>(g.h, nn_precond);
          for(int d = 0; d < 3; ++d) {
              precond.c2nd[d][1] = 1/12.;
              precond.c2nd[d][0] = 2/12.; // stencil [1/4 1/2 1/4] in all 3 directions
          } // d
          
          // the kinetic energy operator
          kinetic = finite_difference::finite_difference_t<real_fd_t>(g.h, nn);
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

      status_t get_atom_coeffs(real_t** atom_coeffs, real_t const psi[], int const echo=0) const {
          return _grid_operator<real_t, real_fd_t, D0>(nullptr, psi, grid, atoms, images, 0, boundary_phase.data(), nullptr, nullptr, echo, atom_coeffs);
      } // get_atom_coeffs

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
      int  get_natoms()  const { return atoms.size(); }
      int  get_numax(int const ia) const { return atoms[ia].numax(); }
  }; // class grid_operator_t
  
  status_t all_tests(int const echo=0);

} // namespace grid_operators
