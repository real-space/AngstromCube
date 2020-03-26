#pragma once

#include "status.hxx" // status_t

#include "grid_operators.hxx"

#include "atom_image.hxx" // ::sho_atom_t
#include "real_space_grid.hxx" // ::grid_t
#include "finite_difference.hxx" // ::finite_difference_t<real_t>
#include "vector_math.hxx" // ::vec<N,T>
#include "boundary_condition.hxx" // ::periodic_images
#include "recorded_warnings.hxx" // error

#include "chemical_symbol.hxx" // ::get
#include "display_units.h" // Ang, _Ang

#define DEBUG

namespace grid_operators {

  template<typename real_t, typename real_fd_t, int D0=1>
  status_t _grid_operator(real_t Hpsi[] // result
                        , real_t const psi[] // input wave functions
                        , real_space_grid::grid_t<D0> const &g // 3D Cartesian grid descriptor
                        , std::vector<atom_image::sho_atom_t> const &a // atoms
                        , int const h0s1 // index controlling which matrix of a[ia] we are multiplying, 0:Hamiltonian or 1:overlap
                        , double const *boundary_phase=nullptr // phase shifts at the boundary [optional]
                        , finite_difference::finite_difference_t<real_fd_t> const *fd=nullptr // finite difference [optional]
                        , double const *potential=nullptr // diagonal potential operator [optional]
                        , int const echo=0
                        , real_t** atomic_projection_coefficients=nullptr
                         ) {
      status_t stat(0);

      if (Hpsi) {
          if (fd) {
              stat += finite_difference::Laplacian(Hpsi, psi, g, *fd); // prefactor -0.5 moved into fd-coefficients
              if (echo > 8) printf("# %s Apply Laplacian, status=%i\n", __func__, stat);
          } else {
              set(Hpsi, g.all(), real_t(0)); // clear
          } // fd

          size_t const nzyx = g[2] * g[1] * g[0];
          if (echo > 8) printf("# %s Apply %s operator\n", __func__, potential ? "potential" : "unity");
          for(size_t izyx = 0; izyx < nzyx; ++izyx) {
              real_fd_t const V = potential ? potential[izyx] : real_fd_t(1); // apply potential or the unity operation of the overlap operator
              for(int i0 = 0; i0 < D0; ++i0) {
                  Hpsi[izyx*D0 + i0] += V * psi[izyx*D0 + i0];
              } // i0 vectorization
          } // izyx
      } // Hpsi != nullptr

      int constexpr echo_sho = 0;

      if (h0s1 >= 0) {
        
          int const na = a.size(); // number of atoms
          std::vector<std::vector<real_t>> V_atom_coeff(na);
          std::vector<std::vector<real_t>>   atom_coeff(na);
          
          for(int ia = 0; ia < na; ++ia) {
              int const numax = a[ia].numax();
              int const ncoeff = sho_tools::nSHO(numax);
              atom_coeff[ia] = std::vector<real_t>(ncoeff*D0, 0.0);

              if (a[ia].nimages() > 1) {
                  for(int ii = 0; ii < a[ia].nimages(); ++ii) {
                      std::vector<real_t> image_coeff(ncoeff*D0, 0.0); // can we omit the initialization?
                  
                      stat += sho_projection::sho_project(image_coeff.data(), numax, a[ia].pos(ii), a[ia].sigma(), psi, g, echo_sho);

                      real_t const Bloch_factor = 1.0; // ToDo: use k-dependent Bloch-factors
                      add_product(atom_coeff[ia].data(), ncoeff*D0, image_coeff.data(), Bloch_factor);
                  } // ii
#ifdef DEBUG
              } else if (a[ia].nimages() < 1) {
                  error("atom #%i has no image!", a[ia].atom_id());
#endif
              } else {
                  stat += sho_projection::sho_project(atom_coeff[ia].data(), numax, a[ia].pos(), a[ia].sigma(), psi, g, echo_sho);
              } // need_images

          } // ia

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
                  V_atom_coeff[ia] = std::vector<real_t>(ncoeff*D0, 0.0);
              
                  { // scope: V_atom_coeff = matrix * atom_coeff
                      typedef vector_math::vec<D0,real_fd_t> temp_vec_t; // computation is performed in real_fd_t precision
                      
                      int const stride = a[ia].stride();
                      assert(stride >= ncoeff); // check internal consistency
                      auto *const mat = a[ia].get_matrix<real_fd_t>(h0s1);
                      // DGEMM-style matrix multiplication
                      for(int i = 0; i < ncoeff; ++i) {
                          temp_vec_t ci(0);
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
              
                  if (a[ia].nimages() > 1) {
                      for(int ii = 0; ii < a[ia].nimages(); ++ii) {
                          real_t const inv_Bloch_factor = 1.0; // ToDo: use k-dependent inverse Bloch-factors
                          std::vector<real_t> V_image_coeff(ncoeff*D0);
                          scale(V_image_coeff.data(), ncoeff*D0, V_atom_coeff[ia].data(), inv_Bloch_factor);
                      
                          stat += sho_projection::sho_add(Hpsi, g, V_image_coeff.data(), numax, a[ia].pos(ii), a[ia].sigma(), echo_sho);
                      } // ii
                  } else {
                      stat += sho_projection::sho_add(Hpsi, g, V_atom_coeff[ia].data(), numax, a[ia].pos(), a[ia].sigma(), echo_sho);
                  } // need images
                  
              } // ia
              
          } // Hpsi != nullptr

      } // h0s1 >= 0
      return stat;
  } // _grid_operator

  // Hamiltonian
  template<typename real_t, typename real_fd_t=double, int D0=1>
  status_t grid_Hamiltonian(real_t Hpsi[] // result
                          , real_t const psi[] // input wave functions
                          , real_space_grid::grid_t<D0> const &g // 3D Cartesian grid descriptor
                          , std::vector<atom_image::sho_atom_t> const &a
                          , finite_difference::finite_difference_t<real_fd_t> const &fd // finite difference
                          , double const potential[] // diagonal potential operator
                          , double const *boundary_phase=nullptr // phase shifts at the boundary [optional]
                          , int const echo=0
  ) { return _grid_operator(Hpsi, psi, g, a, 0, boundary_phase, &fd, potential); }

  // Overlap operator
  template<typename real_t, typename real_fd_t=double, int D0=1>
  status_t grid_Overlapping(real_t Spsi[] // result
                          , real_t const psi[] // input wave functions
                          , real_space_grid::grid_t<D0> const &g // 3D Cartesian grid descriptor
                          , std::vector<atom_image::sho_atom_t> const &a
                          , double const *boundary_phase=nullptr // phase shifts at the boundary [optional]
                          , int const echo=0
  ) { return _grid_operator<real_t, real_fd_t, D0>(Spsi, psi, g, a, 1, boundary_phase); }

  // Overlap operator
  template<typename real_t, int D0=1>
  status_t get_atom_coeffs(real_t** atom_coeffs // result
                          , real_t const psi[] // input wave functions
                          , real_space_grid::grid_t<D0> const &g // 3D Cartesian grid descriptor
                          , std::vector<atom_image::sho_atom_t> const &a
                          , double const *boundary_phase=nullptr // phase shifts at the boundary [optional]
                          , int const echo=0
  ) { return _grid_operator<real_t, real_t, D0>(nullptr, psi, g, a, 0, boundary_phase, nullptr, nullptr, 0, atom_coeffs); }

  //
  // idea: the identity operation of the Overlap operator could be implemented with a
  //       finite difference stencil that has nn[] = {0, -1, -1} (-1 means that there is no pass)
  //       because then, all three operators, Hamiltonian, Overlapping and Preconditioner share
  //       the same structure, whereas the latter has a FD-stencil with all positive coefficients
  //

  template <typename real_t, typename real_fd_t=double, int D0=1>
  class grid_operator_t
  {
    private:
      void _constructor(real_space_grid::grid_t<D0> const & g // real space grid descriptor
               , double const *local_potential // may be nullptr
               , int const nn_precond
               , int const nn_kinetic=8) {
        
//           int const nn_precond = control::get("conjugate_gradients.precond", 1.);

          // the kinetic energy operator
          kinetic = finite_difference::finite_difference_t<real_fd_t>(g.h, nn_kinetic);
          kinetic.scale_coefficients(-0.5); // prefactor of the kinetic energy in Hartree atomic units

          // the local effective potential
          potential = std::vector<double>(g.all(), 0.0); // init as zero everywhere
          if (nullptr != local_potential) {
              set(potential.data(), g.all(), local_potential); // copy data in
          } // local potential given

          // this simple grid-based preconditioner is a diffusion stencil
          preconditioner = finite_difference::finite_difference_t<real_t>(g.h, nn_precond);
          for(int d = 0; d < 3; ++d) {
              preconditioner.c2nd[d][1] = 1/12.;
              preconditioner.c2nd[d][0] = 2/12.; // stencil [1/4 1/2 1/4] in all 3 directions, normalized
          } // d

      } // _constructor
      
    public:

      grid_operator_t(real_space_grid::grid_t<D0> const & g // real space grid descriptor
                      , std::vector<atom_image::sho_atom_t> const & a // is moved to atoms
                      , double const *local_potential=nullptr // effective local potential
                      , int const nn_precond=1) // range of the preconditioner
          : grid(g), atoms(a), has_precond(nn_precond > 0), has_overlap(true) {
          _constructor(grid, local_potential, nn_precond);
      } // constructor with atoms

      grid_operator_t(real_space_grid::grid_t<D0> const & g // real space grid descriptor
                      , double const *local_potential=nullptr // effective local potential
                      , int const nn_precond=1) // range of the preconditioner
          : grid(g), atoms(0), has_precond(nn_precond > 0), has_overlap(true) {
          _constructor(grid, local_potential, nn_precond);
      } // constructor without atoms

      status_t Hamiltonian(real_t Hpsi[], real_t const psi[], int const echo=0) const {
          return _grid_operator(Hpsi, psi, grid, atoms, 0, boundary_phase.data(), &kinetic, potential.data(), echo);
      } // Hamiltonian

      status_t Overlapping(real_t Spsi[], real_t const psi[], int const echo=0) const {
          return _grid_operator<real_t, real_fd_t, D0>(Spsi, psi, grid, atoms, 1, boundary_phase.data(), nullptr, nullptr, echo);
      } // Overlapping

      status_t Conditioner(real_t Cpsi[], real_t const psi[], int const echo=0) const {
          return _grid_operator(Cpsi, psi, grid, atoms, -1, boundary_phase.data(), &preconditioner, nullptr, echo);
      } // Pre-Conditioner

      status_t get_atom_coeffs(real_t** atom_coeffs, real_t const psi[], int const echo=0) const {
          return _grid_operator<real_t, real_fd_t, D0>(nullptr, psi, grid, atoms, 0, boundary_phase.data(), nullptr, nullptr, echo, atom_coeffs);
      } // get_atom_coeffs

    private:
      real_space_grid::grid_t<D0> grid;
      std::vector<atom_image::sho_atom_t> atoms;
      std::vector<double> boundary_phase; // could be real_t or real_fd_t in the future
      std::vector<double> potential;
      finite_difference::finite_difference_t<real_fd_t> kinetic;
      finite_difference::finite_difference_t<real_t> preconditioner;
      bool has_precond;
      bool has_overlap;
    public:
      real_space_grid::grid_t<D0> const & get_grid() const { return grid; }
      bool use_precond() const { return has_precond; }
      bool use_overlap() const { return has_overlap; }
      int  get_natoms()  const { return atoms.size(); }
      int  get_numax(int const ia) const { return atoms[ia].numax(); }
  }; // class grid_operator_t
  
  
  
  
  // prepare the sho_atoms for grid_operator_t
  template <int D0=1> // vectorization 
  status_t list_of_atoms(std::vector<atom_image::sho_atom_t> & a
                       , double const xyzZins[] // data layout [na][8], collection of different quantities
                       , int const na // number of atoms
                       , int const stride // typically stride=8
                       , real_space_grid::grid_t<D0> const & g // actually we need cell info here, not grid, so the <D0> templating could be dropped
                       , int const echo=9
                       , double const *const *const atom_matrices=nullptr // data layout am[na][2*ncoeff[ia]^2]
                       , float const rcut=18 // sho_projection usually ends at 9*sigma
                        ) {
      status_t stat(0);

      double const cell[] = {g[0]*g.h[0], g[0]*g.h[1], g[2]*g.h[2]};
      double *periodic_image_positions{nullptr};
      int const n_periodic_images = boundary_condition::periodic_images(
            &periodic_image_positions, cell, g.boundary_conditions(), rcut, echo);
      if (echo > 1) printf("# %s consider %d periodic images\n", __FILE__, n_periodic_images);

      a.resize(na);
      for(int ia = 0; ia < na; ++ia) {
          double const *pos =             & xyzZins[ia*stride + 0];
          double const Z =                  xyzZins[ia*stride + 3];
          int32_t const atom_id = int32_t(  xyzZins[ia*stride + 4]);
          int const numax = (int)std::round(xyzZins[ia*stride + 5]);
          double const sigma =              xyzZins[ia*stride + 6];
          assert(sigma > 0);
          int const Zi = std::round(Z);
          a[ia] = atom_image::sho_atom_t(sigma, numax, atom_id, pos, Zi);
          if (n_periodic_images > 1) {
              a[ia].set_image_positions(pos, n_periodic_images, periodic_image_positions);
          } // more than one periodic image
          char Symbol[4] = {0,0,0,0};
          chemical_symbol::get(Symbol, Z);
          if (echo > 3) printf("# %s %s %g %g %g %s has %d images, sigma %g %s, numax %d (atom_id %i)\n", __func__, 
              Symbol, pos[0]*Ang, pos[1]*Ang, pos[2]*Ang, _Ang, n_periodic_images, sigma*Ang, _Ang, numax, atom_id);

          if (atom_matrices) {
              // set atomic Hamiltonian and charge deficit matrices
              int const ncoeff = sho_tools::nSHO(numax);
              for(int h0s1 = 0; h0s1 <= 1; ++h0s1) {
                  stat += a[ia].set_matrix(&atom_matrices[ia][h0s1*ncoeff*ncoeff], ncoeff, ncoeff, h0s1);
              } // 0:Hamiltonian and 1:overlap/charge deficit
          } // atom_matrices != nullptr
      } // ia
      delete[] periodic_image_positions;
      if (nullptr == atom_matrices) warn("Atom-centered matrices were not set");
      return stat;
  } // list_of_atoms
 
  status_t all_tests(int const echo=0);

} // namespace grid_operators
