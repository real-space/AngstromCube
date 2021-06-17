#pragma once

#include <cmath> // std::sqrt

#include "radial_grid.h" // radial_grid_t
#include "status.hxx" // status_t

namespace radial_grid {

  double constexpr default_anisotropy = 0.01;
  float constexpr default_Rmax = 9.45;
  
  char const equation_exponential[] = "r=a*(exp(d*i)-1)";
  char const equation_equidistant[] = "r=a*i/n";
  char const equation_reciprocal[]  = "r=a*i/(n-i)";

  radial_grid_t* create_radial_grid( // returns a pointer to a new radial grid descriptor
        int const npoints // number of grid points
      , float const rmax=default_Rmax // [optional] largest radius
      , char const *equation=nullptr // [optional] how to generate the grid
      , double const anisotropy=default_anisotropy // [optional] anisotropy parameter
  ); // declaration only

  radial_grid_t* create_pseudo_radial_grid(
        radial_grid_t const & tru // radial grid for true quantities (usually goes down to 0)
      , double const r_min=1e-3 // start radius
      , int const echo=0 // log-level
  ); // declaration only

  inline radial_grid_t* create_default_radial_grid(float const Z_nucleons=0) {
      return create_radial_grid(250*std::sqrt(Z_nucleons + 9.) + .5); }

  void destroy_radial_grid(radial_grid_t* g, char const *name=""); // radial grid descriptor

  int find_grid_index(radial_grid_t const &g, double const radius);

  double get_prefactor(radial_grid_t const & g);

  status_t all_tests(int const echo=0); // declaration only
  
} // namespace radial_grid
