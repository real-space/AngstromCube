#pragma once

#include <cstdio>  // std::printf, ::fopen, ::fprintf, ::fclose
#include <cmath>   // std::abs, ::round, ::pow
#include <cstdint> // int8_t
#include <complex> // std::real, ::imag
#include <vector>  // std::vector<T>
#include <algorithm> // std::min, ::max

#include "status.hxx" // status_t

#include "atom_image.hxx" // ::sho_atom_t
#include "real_space.hxx" // ::grid_t
#include "finite_difference.hxx" // ::stencil_t<real_t>
#include "boundary_condition.hxx" // ::periodic_images
#include "recorded_warnings.hxx" // error, warn
#include "complex_tools.hxx" // conjugate, to_complex_t
#include "inline_math.hxx" // intpow
#include "sho_projection.hxx" // ::sho_project, ::sho_add, ::sho_prefactor

#include "chemical_symbol.hxx" // ::get
#include "display_units.h" // Ang, _Ang
#include "print_tools.hxx" // printf_vector

#ifdef DEVEL
  #include "control.hxx" // ::get
#endif // DEVEL

// #define DEBUG

namespace grid_operators {

  template <typename complex_t>
  complex_t Bloch_phase(complex_t const boundary_phase[3][2], int8_t const idx[3], int const inv=0) {
      complex_t Bloch_factor(1);
      for (int d = 0; d < 3; ++d) {
          int const i01 = ((idx[d] < 0) + inv) & 0x1;
          Bloch_factor *= intpow(boundary_phase[d][i01], std::abs(idx[d]));
      } // d
      return Bloch_factor;
  } // Bloch_phase


  template <typename complex_t, typename real_fd_t=double>
  status_t _grid_operation(
        complex_t Hpsi[] // result
      , complex_t const psi[] // input wave functions
      , real_space::grid_t const & g // 3D Cartesian grid descriptor
      , std::vector<atom_image::sho_atom_t> const & a // atoms
      , int const h0s1 // index controlling which matrix of a[ia] we are multiplying, 0:Hamiltonian or 1:overlap
      , complex_t const boundary_phase[3][2] // phase shifts at the boundary
      , finite_difference::stencil_t<real_fd_t> const *kinetic=nullptr // finite difference [optional]
      , double const *potential=nullptr // diagonal potential operator [optional]
      , int const echo=0 // log-level [optional]
      , complex_t       *const *const atomic_projection_coefficients=nullptr // [optional]
      , complex_t const *const *const start_wave_coefficients=nullptr // [optional]
      , float const scale_sigmas=1 // only during addition (used for start wave functions) [optional]
  ) {
      using real_t = decltype(std::real(complex_t(1)));
      using atom_matrix_t = decltype(std::real(real_fd_t(1)));

      status_t stat(0);

      if (Hpsi) {
          if (kinetic) {
              if (psi) {
                  stat += finite_difference::apply(Hpsi, psi, g, *kinetic, 1, boundary_phase);
                  if (echo > 8) std::printf("# %s Apply Laplacian, status=%i\n", __func__, stat);
              } // psi != nullptr
          } else {
              set(Hpsi, g.all(), complex_t(0)); // clear
          } // kinetic

          if (psi) {
              size_t const nzyx = g[2] * g[1] * g[0];
              if (echo > 8) std::printf("# %s Apply %s operator\n", __func__, potential ? "potential" : "unity");
              for (size_t izyx = 0; izyx < nzyx; ++izyx) {
                  real_t const V = potential ? potential[izyx] : 1; // apply potential or the unity operation of the overlap operator
                  Hpsi[izyx] += V * psi[izyx];
              } // izyx
          } else {
              if (echo > 18) std::printf("# %s has no input function\n", __func__);
          } // psi != nullptr
      } else {
          if (echo > 18) std::printf("# %s has no output function\n", __func__);
      } // Hpsi != nullptr

      int const echo_sho = 0*(nullptr != atomic_projection_coefficients)
                         + 0*(nullptr != start_wave_coefficients);

      if (h0s1 >= 0) {

          int const natoms = a.size(); // number of atoms
          std::vector<std::vector<complex_t>> atom_coeff(natoms);

          for (int ia = 0; ia < natoms; ++ia) {
              int const numax = a[ia].numax();
              auto const sigma = a[ia].sigma();
              int const ncoeff = sho_tools::nSHO(numax);
              atom_coeff[ia] = std::vector<complex_t>(ncoeff, complex_t(0));

              if (psi) {
                  if (a[ia].nimages() > 1) {
                      for (int ii = 0; ii < a[ia].nimages(); ++ii) {
                          std::vector<complex_t> image_coeff(ncoeff, 0.0); // can we omit the initialization?

                          stat += sho_projection::sho_project(image_coeff.data(), numax, a[ia].pos(ii), sigma, psi, g, echo_sho);

                          complex_t const Bloch_factor = Bloch_phase(boundary_phase, a[ia].idx(ii));
                          add_product(atom_coeff[ia].data(), ncoeff, image_coeff.data(), Bloch_factor);
                      } // ii
#ifdef DEBUG
                  } else if (a[ia].nimages() < 1) {
                      error("atom #%i has no image!", a[ia].atom_id());
#endif // DEBUG
                  } else {
                      // Gamma point, no periodic images
                      stat += sho_projection::sho_project(atom_coeff[ia].data(), numax, a[ia].pos(), a[ia].sigma(), psi, g, echo_sho);
                  } // need_images
              } // psi != nullptr

              if (atomic_projection_coefficients) {
                  set(atomic_projection_coefficients[ia], ncoeff, atom_coeff[ia].data()); // export
              } // atomic_projection_coefficients
          } // ia

          if (Hpsi) {

              for (int ia = 0; ia < natoms; ++ia) {
                  int const numax = (start_wave_coefficients) ? 3 : a[ia].numax();
                  auto const sigma = a[ia].sigma()*scale_sigmas;
                  int const ncoeff = sho_tools::nSHO(numax);
                  std::vector<complex_t> V_atom_coeff_ia(ncoeff);

                  if (start_wave_coefficients) {
                      assert( 3 == numax ); // this option is only used for start wave functions.
#ifdef DEVEL
                      if (echo > 19) {
                          std::printf("# %s atomic addition coefficients for atom #%i are", __func__, ia);
                          printf_vector(" %g", start_wave_coefficients[ia], ncoeff);
                      } // echo
#endif // DEVEL
                      set(V_atom_coeff_ia.data(), ncoeff, start_wave_coefficients[ia]); // import
                  } else {

                      // scope: V_atom_coeff = matrix * atom_coeff
                      int const stride = a[ia].stride();
                      assert(stride >= ncoeff); // check internal consistency
                      auto *const mat = a[ia].get_matrix<atom_matrix_t>(h0s1);
                      // matrix-vector multiplication
                      for (int i = 0; i < ncoeff; ++i) {
                          complex_t ci(0);
                          for (int j = 0; j < ncoeff; ++j) {
                              auto const am = mat[i*stride + j];
                              auto const cj = atom_coeff[ia][j];
#ifdef DEVEL
//                               if (echo > 9) std::printf("# %s atomic %s matrix for atom #%i mat(%i,%i)= %g\n",
//                                                    __func__, h0s1?"overlap":"hamiltonian", ia, i, j, am);
#endif // DEVEL
                              ci += am * cj;
                          } // j
                          V_atom_coeff_ia[i] = ci;
                      } // i
                      // scope

                  } // start_wave_coefficients

                  if (a[ia].nimages() > 1) {
                      for (int ii = 0; ii < a[ia].nimages(); ++ii) {
                          complex_t const inv_Bloch_factor = Bloch_phase(boundary_phase, a[ia].idx(ii), 1);
                          std::vector<complex_t> V_image_coeff(ncoeff);
                          set(V_image_coeff.data(), ncoeff, V_atom_coeff_ia.data(), inv_Bloch_factor);

                          stat += sho_projection::sho_add(Hpsi, g, V_image_coeff.data(), numax, a[ia].pos(ii), sigma, echo_sho);
                      } // ii
                  } else {
                      // Gamma point, no periodic images
                      stat += sho_projection::sho_add(Hpsi, g, V_atom_coeff_ia.data(), numax, a[ia].pos(), sigma, echo_sho);
                  } // need images

              } // ia

          } // Hpsi != nullptr

      } // h0s1 >= 0
      return stat;
  } // _grid_operation

  //
  // idea: the identity operation of the Overlap operator could be implemented with a
  //       finite difference stencil that has nn[] = {0, -1, -1} (-1 means that there is no pass)
  //       because then, all three operators, Hamiltonian, Overlapping and Preconditioner share
  //       the same structure, whereas the latter has a FD-stencil with all positive coefficients
  //

  inline status_t set_nonlocal_potential(std::vector<atom_image::sho_atom_t> & a
                , double const *const *const atom_matrices // data layout am[natoms][2*ncoeff[ia]^2]
                , int const echo=0) {
      status_t stat(0);
      assert(atom_matrices); // may not be nullptr

#ifdef DEVEL
      double const scale_hs[] = {control::get("hamiltonian.scale.nonlocal.h", 1.),
                                 control::get("hamiltonian.scale.nonlocal.s", 1.)};
      if (1 != scale_hs[0] || 1 != scale_hs[1]) warn("scale PAW contributions to H and S by %g and %g, respectively", scale_hs[0], scale_hs[1]);
#else  // DEVEL
      double constexpr scale_hs[] = {1, 1};
#endif // DEVEL

      for (size_t ia = 0; ia < a.size(); ++ia) {
          // set atomic Hamiltonian and charge deficit matrices
          assert(atom_matrices[ia]);
          int const ncoeff = sho_tools::nSHO(a[ia].numax());
          for (int h0s1 = 0; h0s1 <= 1; ++h0s1) {
              stat += a[ia].set_matrix(&atom_matrices[ia][h0s1*ncoeff*ncoeff], ncoeff, ncoeff, h0s1, scale_hs[h0s1]);
          } // 0:Hamiltonian and 1:overlap/charge deficit
      } // ia
      return stat;
  } // set_nonlocal_potential


  // prepare the sho_atoms for grid_operator_t
  inline // ToDo: move body to grid_operators.cxx since this is not templated
  std::vector<atom_image::sho_atom_t> list_of_atoms(
        double const xyzZins[] // data layout [na][8], collection of different quantities
      , int const na // number of atoms
      , int const stride // typically stride=8
      , real_space::grid_t const & gc // we need cell info here (for periodic_images) and grid spacing (for grid_offset)
      , int const echo=9
      , double const *const *const atom_matrices=nullptr // data layout am[na][2*ncoeff[ia]^2]
      , float const rcut=18 // sho_projection usually ends at 9*sigma
  ) {
      double const cell[] = {gc[0]*gc.h[0], gc[1]*gc.h[1], gc[2]*gc.h[2]};
      double const grid_offset[] = {0.5*(gc[0] - 1)*gc.h[0], 0.5*(gc[1] - 1)*gc.h[1], 0.5*(gc[2] - 1)*gc.h[2]};
      view2D<double> image_positions;
      view2D<int8_t> image_indices;
      int const n_periodic_images = boundary_condition::periodic_images(
            image_positions, cell, gc.boundary_conditions(), rcut, echo, &image_indices);
      if (echo > 1) std::printf("# %s consider %d periodic images\n", __FILE__, n_periodic_images);

      assert(stride >= 7);
      std::vector<atom_image::sho_atom_t> a(na);
      for (int ia = 0; ia < na; ++ia) {
          double pos[3];
          for (int d = 0; d < 3; ++d) {
              pos[d] =    grid_offset[d] + xyzZins[ia*stride + d];
          } // d
          double const Z =                 xyzZins[ia*stride + 3];
          int32_t const atom_id = int32_t( xyzZins[ia*stride + 4]);
          int const numax = int(std::round(xyzZins[ia*stride + 5]));
          double const sigma =             xyzZins[ia*stride + 6];
          assert(sigma > 0);
          int const Zi = std::round(Z);
          a[ia] = atom_image::sho_atom_t(sigma, numax, atom_id, pos, Zi);
          a[ia].set_image_positions(pos, n_periodic_images, &image_positions, &image_indices);

          char Symbol[4]; chemical_symbol::get(Symbol, Z);
          double const *apos = &xyzZins[ia*stride + 0];
          if (echo > 3) std::printf("# %s %s %g %g %g %s has %d images, sigma %g %s, numax %d (atom_id %i)\n", __func__,
              Symbol, apos[0]*Ang, apos[1]*Ang, apos[2]*Ang, _Ang, n_periodic_images, sigma*Ang, _Ang, numax, atom_id);
          if (echo > 3) std::printf("# %s %s %g %g %g (rel) has %d images, sigma %g %s, numax %d (atom_id %i)\n", __func__,
              Symbol, pos[0]*gc.inv_h[0], pos[1]*gc.inv_h[1], pos[2]*gc.inv_h[2], n_periodic_images, sigma*Ang, _Ang, numax, atom_id);

      } // ia

      status_t stat(0);
      if (nullptr != atom_matrices) {
          stat += set_nonlocal_potential(a, atom_matrices, echo);
      } else {
//           warn("Atom-centered matrices were not set"); // this warning is deactivated as in the default use case,
            // we will create the an instance of the grid_operator_t before we know the atom_matrices, so nullptr is often passed but it is ok.
      } // atom_matrices != nullptr
      if (stat) warn("set_nonlocal_potential returned status sum %d", int(stat));

      return a;
  } // list_of_atoms

  inline std::vector<atom_image::sho_atom_t> empty_list_of_atoms()
      {  std::vector<atom_image::sho_atom_t> a(0); return a; }

  template <typename wave_function_t, typename real_FiniDiff_t=wave_function_t>
  class grid_operator_t
  {
    public:
      typedef wave_function_t complex_t;
      typedef real_FiniDiff_t real_fd_t;

    public:

      grid_operator_t(
            real_space::grid_t const & g // real space grid descriptor
          , std::vector<atom_image::sho_atom_t> const & a // is moved to atoms
          , double const *local_potential=nullptr // effective local potential
          , int const nn_precond=1 // range of the preconditioner
          , int const nn_kinetic=8
          , int const echo=0
      ) : grid(g), atoms(a), has_precond(nn_precond > 0), has_overlap(true) {

//        int const nn_precond = control::get("conjugate_gradients.precond", 1.);

          // the kinetic energy operator
          kinetic = finite_difference::stencil_t<real_fd_t>(g.h, nn_kinetic, -0.5);
          // -0.5: prefactor of the kinetic energy in Hartree atomic units

          // the local effective potential, ToDo: separate this because we might want to call it later
          potential = std::vector<double>(g.all(), 0.0); // init as zero everywhere
#ifdef DEVEL
          double const scale_k = control::get("hamiltonian.scale.kinetic", 1.);
          if (1 != scale_k) {
              kinetic.scale_coefficients(scale_k);
              warn("kinetic energy is scaled by %g", scale_k);
          } // scale_k != 1
#endif // DEVEL
          set_potential(local_potential, g.all(), nullptr, echo);

          set_kpoint<double>(); // set Gamma point, silent

//        std::printf("\n# here: %s %s:%d\n\n", __func__, __FILE__, __LINE__);

          // this simple grid-based preconditioner is a diffusion stencil
          preconditioner = finite_difference::stencil_t<complex_t>(g.h, std::min(1, nn_precond));
          for (int d = 0; d < 3; ++d) {
              preconditioner.c2nd[d][1] = 1/12.;
              preconditioner.c2nd[d][0] = 2/12.; // stencil [1/4 1/2 1/4] in all 3 directions, normalized
          } // d

      } // constructor with atoms

      status_t Hamiltonian(complex_t Hpsi[], complex_t const psi[], int const echo=0) const {
          return _grid_operation(Hpsi, psi, grid, atoms, 0, boundary_phase, &kinetic, potential.data(), echo);
      } // Hamiltonian

      status_t Overlapping(complex_t Spsi[], complex_t const psi[], int const echo=0) const {
          return _grid_operation<complex_t, real_fd_t>(Spsi, psi, grid, atoms, 1, boundary_phase, nullptr, nullptr, echo);
      } // Overlapping

      status_t Conditioner(complex_t Cpsi[], complex_t const psi[], int const echo=0) const {
          return _grid_operation(Cpsi, psi, grid, atoms, -1, boundary_phase, &preconditioner, nullptr, echo);
      } // Pre-Conditioner

      status_t get_atom_coeffs(complex_t *const *const atom_coeffs, complex_t const psi[], int const echo=0) const {
          return _grid_operation<complex_t, real_fd_t>(nullptr, psi, grid, atoms, 0, boundary_phase, nullptr, nullptr, echo, atom_coeffs);
      } // get_atom_coeffs

      status_t get_start_waves(complex_t psi0[], complex_t const *const *const atom_coeffs, float const scale_sigmas=1, int const echo=0) const {
          return _grid_operation<complex_t, real_fd_t>(psi0, nullptr, grid, atoms, 0, boundary_phase, nullptr, nullptr, echo, nullptr, atom_coeffs, scale_sigmas);
      } // get_start_waves

      double get_volume_element() const { return grid.dV(); }
      size_t get_degrees_of_freedom() const { return size_t(grid[2]) * size_t(grid[1]) * size_t(grid[0]); }

      template <typename real_t=double>
      void set_kpoint(real_t const kpoint[3]=nullptr, int const echo=0) {
          if (nullptr != kpoint) {
              std::complex<double> const minus1(-1);
              if (echo > 5) std::printf("# grid_operator.%s", __func__);
              for (int d = 0; d < 3; ++d) {
                  k_point[d] = kpoint[d];
                  if (is_complex<complex_t>()) {
                      std::complex<double> phase = std::pow(minus1, 2*k_point[d]);
                      boundary_phase[d][1] = to_complex_t<complex_t, double>(phase);
                  } else {
                      int const kk = std::round(2*k_point[d]);
                      assert(2*k_point[d] == kk);
                      boundary_phase[d][1] = kk ? -1 : 1;
                  } // real phase factor
                  boundary_phase[d][0] = conjugate(boundary_phase[d][1]);
                  if (echo > 5) std::printf("   %c: %g ph(%g, %g)", 'x'+d, k_point[d],
                      std::real(boundary_phase[d][1]), std::imag(boundary_phase[d][1]));
              } // d
              if (echo > 5) std::printf("\n");
          } else {
              if (echo > 6) std::printf("# grid_operator.%s resets kpoint to Gamma\n", __func__);
              set(k_point, 3, 0.0);
              set(boundary_phase[0], 3*2, complex_t(1));
          } // nullptr
      } // set_kpoint

      status_t set_potential(
            double const *local_potential=nullptr
          , size_t const ng=0
          , double const *const *const atom_matrices=nullptr
          , int const echo=0
      ) {
          status_t stat(0);
          if (echo > 0) std::printf("# %s %s\n", __FILE__, __func__);
          if (nullptr != local_potential) {
              if (ng == grid.all()) {
#ifdef DEVEL
                  double const scale_p = control::get("hamiltonian.scale.potential", 1.);
                  if (1 != scale_p) warn("local potential is scaled by %g", scale_p);
#else  // DEVEL
                  double constexpr scale_p = 1;
#endif // DEVEL
                  set(potential.data(), ng, local_potential, scale_p); // copy data in
                  if (echo > 0) std::printf("# %s %s local potential copied (%ld elements)\n", __FILE__, __func__, ng);
              } else {
                  error("expect %ld element for the local potential but got %ld\n", grid.all(), ng);
                  ++stat;
              }
          } else {
              if (echo > 0) std::printf("# %s %s no local potential given!\n", __FILE__, __func__);
              ++stat;
          } // local potential given

          if (nullptr != atom_matrices) {
              stat += set_nonlocal_potential(atoms, atom_matrices, echo);
              if (echo > 0) std::printf("# %s %s non-local matrices copied for %ld atoms\n", __FILE__, __func__, atoms.size());
          } else {
              if (echo > 0) std::printf("# %s %s no non-local matrices given!\n", __FILE__, __func__);
              ++stat;
          } // atom_matrices
          if (stat && (echo > 0)) std::printf("# %s %s returns status=%i\n", __FILE__, __func__, int(stat));
          return stat;
      } // set_potential

      void construct_dense_operator(complex_t Hmat[], complex_t Smat[], size_t const stride, int const echo=0) {
          // assume shapes Hmat[][stride], Smat[][stride]
          size_t const ndof = grid.all();
          if (echo > 1) { std::printf("\n# %s with %ld x %ld (stride %ld)\n", __func__, ndof, ndof, stride); std::fflush(stdout); }
          assert(ndof <= stride);
          std::vector<complex_t> psi(ndof);
          for (size_t dof = 0; dof < ndof; ++dof) {
              set(psi.data(), ndof, complex_t(0));
              psi[dof] = 1;
              Hamiltonian(&Hmat[dof*stride], psi.data(), echo);
              Overlapping(&Smat[dof*stride], psi.data(), echo);
          } // dof
          if (echo > 1) std::printf("# %s done\n\n", __func__);
      } // construct_dense_operator

      int write_to_file( // TODO: could be moved out of the templated class
            int const echo=0
          , char const *const fileformat="xml" // or "json"
          , double const energy_min_max_Fermi[3]=nullptr
          , char const *filename=nullptr // filename, nullptr: use a default name
          , char const *pathname="."
      ) {
#ifdef  DEVEL
          char file_name_buffer[512];
          if (nullptr == filename) {
              std::snprintf(file_name_buffer, 511, "%s/%s.%s", pathname, "Hmt", fileformat);
              filename = file_name_buffer;
          } // generate a default file name
          if (echo > 0) std::printf("# %s filename=%s\n", __func__, filename);

          auto *const f = std::fopen(filename, "w");
          if (nullptr == f) {
              if (echo > 0) std::printf("# %s Error opening file %s for writing!\n", __func__, filename);
              return __LINE__;
          } // failed to open

          if ('x' == (*fileformat | 32)) {

              // XML file format
              std::fprintf(f, "<?xml version=\"%.1f\"?>\n", 1.0);
              std::fprintf(f, "<grid_Hamiltonian version=\"%.1f\">\n", 0.);
              std::fprintf(f, "  <!-- Units: Hartree and Bohr radii. -->\n");
              if (nullptr != energy_min_max_Fermi) {
                  std::fprintf(f, "  <spectrum min=\"%.6f\" max=\"%.6f\" Fermi=\"%.6f\"/>\n",
                      energy_min_max_Fermi[0], energy_min_max_Fermi[1], energy_min_max_Fermi[2]);
              } // energy_min_max_Fermi

              std::fprintf(f, "  <sho_atoms number=\"%ld\">\n", atoms.size());
              for (int ia = 0; ia < atoms.size(); ++ia) {
                  std::fprintf(f, "    <atom gid=\"%i\">\n", atoms[ia].atom_id());
                  auto const pos = atoms[ia].pos();
                  std::fprintf(f, "      <position x=\"%.12f\" y=\"%.12f\" z=\"%.12f\"/>\n", pos[0], pos[1], pos[2]);
                  int const numax = atoms[ia].numax();
                  std::fprintf(f, "      <projectors type=\"sho\" numax=\"%d\" sigma=\"%.12f\"/>\n",
                                                                  numax, atoms[ia].sigma());
                  int const nSHO = sho_tools::nSHO(numax);
                  for (int h0s1 = 0; h0s1 < 2; ++h0s1) {
                      auto const mat = atoms[ia].template get_matrix<double>(h0s1);
                      auto const stride = atoms[ia].stride();
                      auto const tag = h0s1 ? "overlap" : "hamiltonian";
                      std::fprintf(f, "      <%s>\n", tag);
                      for (int i = 0; i < nSHO; ++i) {
                          std::fprintf(f, "        ");
                          for (int j = 0; j < nSHO; ++j) {
                              std::fprintf(f, " %.15e", mat[i*stride + j]);
                          } // j
                          std::fprintf(f, "\n");
                      } // i
                      std::fprintf(f, "      </%s>\n", tag);
                  } // h0s1
                  std::fprintf(f, "    </atom>\n");
              } // ia
              std::fprintf(f, "  </sho_atoms>\n");

              std::fprintf(f, "  <spacing x=\"%.17f\" y=\"%.17f\" z=\"%.17f\"/>\n", grid.h[0], grid.h[1], grid.h[2]);
              std::fprintf(f, "  <boundary x=\"%d\" y=\"%d\" z=\"%d\"/>\n", grid.boundary_condition(0), grid.boundary_condition(1), grid.boundary_condition(2));
              std::fprintf(f, "  <potential nx=\"%d\" ny=\"%d\" nz=\"%d\">", grid[0], grid[1], grid[2]);
              for (int izyx = 0; izyx < grid.all(); ++izyx) {
                  if (0 == (izyx & 3)) std::fprintf(f, "\n    ");
                  std::fprintf(f, " %.15f", potential[izyx]);
              } // ia
              std::fprintf(f, "\n  </potential>\n");

              // XML file tail
              std::fprintf(f, "</grid_Hamiltonian>\n");

          } else { // .xml or .json

              // JSON file format
              std::fprintf(f, "{\n"); // open grid_Hamiltonian
              std::fprintf(f, "  \"comment\": \"grid_Hamiltonian in units of Hartree and Bohr radii\"\n");

              if (nullptr != energy_min_max_Fermi) {
                  std::fprintf(f, "  ,\"spectrum\": {\"min\": %.6f, \"max\": %.6f, \"Fermi\": %.6f}\n",
                      energy_min_max_Fermi[0], energy_min_max_Fermi[1], energy_min_max_Fermi[2]);
              } // energy_min_max_Fermi

              std::fprintf(f, " ,\"sho_atoms\":\n  {\n");
              std::fprintf(f, "    \"number\": %ld\n", atoms.size());
              std::fprintf(f, "   ,\"atoms\": [\n");
              for (int ia = 0; ia < atoms.size(); ++ia) {
                  std::fprintf(f, "     %c{\"atom_id\": %i\n", ia?',':' ', atoms[ia].atom_id());
                  auto const pos = atoms[ia].pos();
                  std::fprintf(f, "      ,\"position\": [%.12f, %.12f, %.12f]\n", pos[0], pos[1], pos[2]);
                  int const numax = atoms[ia].numax();
                  std::fprintf(f, "      ,\"projectors\": {\"type\": \"sho\", \"numax\": %d, \"sigma\": %.12f}\n",
                                                                                numax, atoms[ia].sigma());
                  int const nSHO = sho_tools::nSHO(numax);
                  for (int h0s1 = 0; h0s1 < 2; ++h0s1) {
                      auto const mat = atoms[ia].template get_matrix<double>(h0s1);
                      auto const stride = atoms[ia].stride();
                      auto const tag = h0s1 ? "overlap" : "hamiltonian";
                      std::fprintf(f, "      ,\"%s\":\n", tag);
                      for (int i = 0; i < nSHO; ++i) {
                          std::fprintf(f, "        %c", i?',':'['); // open row
                          for (int j = 0; j < nSHO; ++j) {
                              if (0 == (j & 0x3)) std::fprintf(f, "\n          ");
                              std::fprintf(f, "%c%.15e", j?',':'[', mat[i*stride + j]);
                          } // j
                          std::fprintf(f, "]\n"); // close row
                      } // i
                      std::fprintf(f, "        ]\n"); // close matrix
                  } // h0s1
                  std::fprintf(f, "      }\n"); // close atom
              } // ia
              std::fprintf(f, "    ]\n"); // close atoms
              std::fprintf(f, "  }\n"); // close sho_atoms

              std::fprintf(f, " ,\"spacing\": [%.17f, %.17f, %.17f]\n", grid.h[0], grid.h[1], grid.h[2]);
              std::fprintf(f, " ,\"boundary\": [%d, %d, %d]\n", grid.boundary_condition(0), grid.boundary_condition(1), grid.boundary_condition(2));
              std::fprintf(f, " ,\"potential\": {\n"); // open potential
              std::fprintf(f, "    \"grid\": [%d, %d, %d]\n", grid[0], grid[1], grid[2]);
              std::fprintf(f, "   ,\"values\":");
              for (int izyx = 0; izyx < grid.all(); ++izyx) {
                  if (0 == (izyx & 0x3)) std::fprintf(f, "\n    ");
                  std::fprintf(f, "%c%.15f", izyx?',':'[', potential[izyx]);
              } // ia
              std::fprintf(f, "]\n"); // close values
              std::fprintf(f, "  }\n"); // close potential
              std::fprintf(f, "}\n"); // close grid_Hamiltonian

              // get the result file checked at https://jsonlint.com

          } // .xml or .json

          std::fclose(f);
          if (echo > 3) std::printf("# file %s written\n", filename);
          return 0; // 0:success
#else  // DEVEL
          return -1; // activate -D DEVEL for export function
#endif // DEVEL
      } // write_to_file

    private:

      real_space::grid_t grid;
      std::vector<atom_image::sho_atom_t> atoms;
      complex_t boundary_phase[3][2];
      double k_point[3];
      std::vector<double> potential;
//    std::vector<complex_t> scale_factors;
      finite_difference::stencil_t<real_fd_t> kinetic;
      finite_difference::stencil_t<complex_t> preconditioner;
      bool has_precond;
      bool has_overlap;

    public:

      real_space::grid_t const & get_grid() const { return grid; }
      bool use_precond() const { return has_precond; }
      bool use_overlap() const { return has_overlap; }
      int  get_natoms()  const { return atoms.size(); }
      int    get_numax(int const ia) const { return atoms[ia].numax(); }
      double get_sigma(int const ia) const { return atoms[ia].sigma(); }

  }; // class grid_operator_t

#ifdef NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t class_test(int const echo=9) {
      status_t stat(0);
      int const dims[] = {12, 13, 14};
      int const all = dims[0]*dims[1]*dims[2];
      grid_operator_t<double> op(dims, empty_list_of_atoms());
      std::vector<double> psi(all, 1.0), Hpsi(all);
      stat += op.Hamiltonian(Hpsi.data(),  psi.data(), echo);
      stat += op.Overlapping( psi.data(), Hpsi.data(), echo);
      stat += op.Conditioner(Hpsi.data(),  psi.data(), echo);
      return stat;
  } // class_test

  inline status_t class_with_atoms_test(int const echo=9) {
      status_t stat(0);
      real_space::grid_t g(36, 25, 24);
      //                          x   y   z    Z     id nu sigma dummy
      double const xyzZinso[] = {.1, .2, -4,  13.0,  767, 3, 1.5, 9e9,  // Al
                                -.1, .2,  3,  15.1,  757, 4, 1.7, 8e8}; // P*
      auto const a = list_of_atoms(xyzZinso, 2, 8, g, echo);
      grid_operator_t<double> op(g, a);
      std::vector<double> psi(g.all(), 1.0), Hpsi(g.all());
      stat += op.Hamiltonian(Hpsi.data(),  psi.data(), echo);
      stat += op.Overlapping( psi.data(), Hpsi.data(), echo);
      stat += op.Conditioner(Hpsi.data(),  psi.data(), echo);
      return stat;
  } // class_with_atoms_test

  inline status_t projector_normalization_test(int const echo=9) {
      status_t stat(0);
      real_space::grid_t g(32, 32, 32); // grid spacings = {1,1,1} by default
      int const numax = 3;
      double const sigma = 2.0;
      //                         x  y  z   Z   id   numax  sigma
      double const xyzZinso[] = {0, 0, 0, 7.2, 747, numax, sigma, 0}; // "Air-atom", 80% N, 20% O --> Z=7.2
      auto const a = list_of_atoms(xyzZinso, 1, 7, g, echo);
      int const ncoeff = sho_tools::nSHO(numax);
      std::vector<double> psi(g.all(), 0.0), a_vec(ncoeff, 0.0), p_vec(ncoeff);
      auto const a_coeff = a_vec.data(), p_coeff = p_vec.data();
      grid_operator_t<double> op(g, a);
      double const f2 = pow2(sho_projection::sho_prefactor(0, 0, 0, sigma));
      double dev{0};
      for (int ic = 0; ic < ncoeff; ++ic) {
          set(psi.data(), g.all(), 0.0); // clear
          a_coeff[ic] = 1; // set
          stat += op.get_start_waves(psi.data(), &a_coeff, 1, echo);
          stat += op.get_atom_coeffs(&p_coeff, psi.data(), echo);
          if (echo > 5) std::printf("# %s%3i  ", __func__, ic);
          for (int jc = 0; jc < ncoeff; ++jc) {
              if (ic == jc) {
                  if (echo > 5) std::printf(" %7.3f", p_coeff[jc]*f2);
              } else {
                  if (echo > 15) std::printf(" %7.3f", p_coeff[jc]*f2);
                  dev = std::max(dev, std::abs(p_coeff[jc]*f2));
              } // ic == jc
          } // jc
          if (echo > 5) std::printf("\n");
          a_coeff[ic] = 0; // reset
      } // ic
      if (echo > 2) std::printf("# %s: largest deviation is %.1e\n", __func__, dev);
      return stat + (dev > 3e-14);
  } // projector_normalization_test

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += class_test(echo);
      stat += class_with_atoms_test(echo);
      stat += projector_normalization_test(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace grid_operators

#ifdef DEBUG
  #undef DEBUG
#endif // DEBUG

#ifdef FULL_DEBUG
  #undef FULL_DEBUG
#endif // FULL_DEBUG
