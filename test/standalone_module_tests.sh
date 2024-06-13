#!/usr/bin/env bash

echo ""
echo "This script tests if each header file module.hxx has all the"
echo "include statements it needs to be compiled standalone."
echo "Linking will require a dependency tree of objects, so we skip that."
echo ""

for module in \
  angular_grid \
  atom_core \
  atom_image \
  bessel_transform \
  bisection_tools \
  bitmap \
  boundary_condition \
  brillouin_zone \
  chemical_symbol \
  complex_tools \
  conjugate_gradients \
  control \
  data_view \
  davidson_solver \
  debug_output \
  dense_operator \
  dense_solver \
  density_generator \
  element_config \
  exchange_correlation \
  fermi_distribution \
  finite_difference \
  fourier_poisson \
  fourier_transform \
  geometry_analysis \
  global_coordinates \
  green_action \
  green_dyadic \
  green_experiments \
  green_function \
  green_input \
  green_kinetic \
  green_parallel \
  green_potential \
  green_sparse \
  green_tests \
  grid_operators \
  hermite_polynomial \
  inline_math \
  iterative_poisson \
  json_reading \
  linear_algebra \
  load_balancer \
  mpi_parallel \
  multi_grid \
  parallel_domains \
  pawxml_export \
  pawxml_import \
  plane_wave \
  poisson_solver \
  potential_generator \
  progress_report \
  pseudo_tools \
  radial_eigensolver \
  radial_grid \
  radial_integrator \
  radial_potential \
  radial_r2grid \
  real_space \
  recorded_warnings \
  scattering_test \
  self_consistency \
  shift_boundary \
  sho_basis \
  sho_hamiltonian \
  sho_overlap \
  sho_potential \
  sho_projection \
  sho_radial \
  sho_tools \
  sho_unitary \
  sigma_config \
  simple_math \
  simple_stats \
  simple_timer \
  single_atom \
  solid_harmonics \
  spherical_harmonics \
  structure_solver \
  symmetry_group \
  unit_system \
  xml_reading \
; do

  ./standalone_module_test.sh $module

done
