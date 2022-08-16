#!/usr/bin/env bash

echo ""
echo "This script tests if each header file module.hxx has all the"
echo "include statements it needs to be compiled standalone."
echo "Linking will require a dependency tree of objects, so we skip that."
echo ""

for module in \
  exchange_correlation \
  spherical_harmonics \
  conjugate_gradients \
  potential_generator \
  hermite_polynomial \
  global_coordinates \
  radial_eigensolver \
  boundary_condition \
  fermi_distribution \
  recorded_warnings \
  finite_difference \
  radial_integrator \
  geometry_analysis \
  density_generator \
  fourier_transform \
  iterative_poisson \
  self_consistency \
  radial_potential \
  bessel_transform \
  parallel_domains \
  structure_solver \
  scattering_test \
  davidson_solver \
  chemical_symbol \
  linear_operator \
  sho_hamiltonian \
  fourier_poisson \
  solid_harmonics \
  bisection_tools \
  green_potential \
  green_function \
  poisson_solver \
  brillouin_zone \
  sho_projection \
  shift_boundary \
  linear_algebra \
  grid_operators \
  dense_operator \
  element_config \
  load_balancer \
  complex_tools \
  vector_layout \
  sho_potential \
  green_kinetic \
  pawxml_import \
  green_sparse \
  green_dyadic \
  green_action \
  simple_stats \
  mpi_parallel \
  angular_grid \
  pseudo_tools \
  simple_timer \
  sigma_config \
  dense_solver \
  json_reading \
  xml_reading \
  unit_system \
  simple_math \
  sho_overlap \
  radial_grid \
  single_atom \
  inline_math \
  sho_unitary \
  plane_wave \
  atom_image \
  real_space \
  multi_grid \
  sho_radial \
  sho_tools \
  sho_basis \
  atom_core \
  data_view \
  control \
; do

    ./standalone_module_test.sh $module

done
