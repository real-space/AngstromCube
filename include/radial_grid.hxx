#pragma once

#include <cmath> // std::sqrt, ::abs

#include "radial_grid.h" // radial_grid_t
#include "status.hxx" // status_t

namespace radial_grid {

  double constexpr default_anisotropy = 0.01;
  float  constexpr default_Rmax = 9.45;

  char constexpr equation_exponential = 'e';
  char constexpr equation_equidistant = '-';
  char constexpr equation_reciprocal  = '/';

  radial_grid_t* create_radial_grid( // returns a pointer to a new radial grid descriptor
        int const npoints // number of grid points
      , float const rmax=default_Rmax // [optional] largest radius
      , char equation='\0' // [optional] how to generate the grid
      , double const anisotropy=default_anisotropy // [optional] anisotropy parameter
  ); // declaration only

  radial_grid_t* create_pseudo_radial_grid(
        radial_grid_t const & tru // radial grid for true quantities (usually goes down to 0)
      , double const r_min=1e-3 // start radius
      , int const echo=0 // log-level
  ); // declaration only

  inline radial_grid_t* create_default_radial_grid(float const Z_protons=0) {
      return create_radial_grid(250*std::sqrt(std::abs(Z_protons) + 9));
  } // create_default_radial_grid

  void destroy_radial_grid(radial_grid_t* g, char const *name=""); // declaration only

  int find_grid_index(radial_grid_t const & g, double const radius); // declaration only

  double get_prefactor(radial_grid_t const & g); // declaration only
  
  char const* get_formula(char const equation='\0'); // declaration only

  status_t all_tests(int const echo=0); // declaration only
  
} // namespace radial_grid
