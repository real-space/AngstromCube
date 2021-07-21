#!/usr/bin/env bash

### This script tests if each header file module.hxx has all the
### include statements it needs to be compiled standalone.
### Linking will require a dependency tree of objects, so we skip that.

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
  green_function \
  poisson_solver \
  brillouin_zone \
  sho_projection \
  shift_boundary \
  linear_algebra \
  grid_operators \
  dense_operator \
  element_config \
  green_kinetic \
  complex_tools \
  vector_layout \
  sho_potential \
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
  atom_core \
  data_view \
  control \
; do

    echo $module

    ## generate a short main function including the header
    echo "#include \"$module.hxx\" // ::all_tests"            > test_me.cxx
    echo "int main() { return int($module::all_tests(6)); }" >> test_me.cxx

    ## compile to check for missing include files
    g++ -std=c++11 \
        -I../include/ \
        -g -pedantic -Wall -O0 \
              -D HAS_NO_MKL \
              -D DEVEL \
              -D HAS_NO_MPI \
              -D NO_UNIT_TESTS \
              -c test_me.cxx

    ## cleanup
    rm -f test_me.cxx test_me.o
done
