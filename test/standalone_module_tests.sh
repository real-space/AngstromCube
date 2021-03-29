#!/usr/bin/env bash

### This script tests if each header file module.hxx has all the include statements it needs

for module in \
  recorded_warnings \
  finite_difference \
  hermite_polynomial \
  spherical_harmonics \
  structure_solver \
  self_consistency \
  conjugate_gradients \
  potential_generator \
  global_coordinates \
  radial_eigensolver \
  boundary_condition \
  fermi_distribution \
  radial_integrator \
  geometry_analysis \
  density_generator \
  fourier_transform \
  iterative_poisson \
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
  complex_tools \
  vector_layout \
  sho_potential \
  mpi_parallel \
  angular_grid \
  pseudo_tools \
  inline_tools \
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
  plane_waves \
  sho_unitary \
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
    echo "#include \"$module.hxx\" // ::all_tests"       > test_me.cxx
    echo "int main() { return $module::all_tests(6); }" >> test_me.cxx

    ## compile
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
