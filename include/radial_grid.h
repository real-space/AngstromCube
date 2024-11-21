#pragma once
// This file is part of AngstromCube under MIT License


struct radial_grid_t {
    int   n = 0; // number of points
    float rmax = 0.f; // max radius
    double const*    r = nullptr; // r[0:n)
    double const*   dr = nullptr; // dr[0:n)
    double const*  rdr = nullptr; // r*dr[0:n)
    double const* r2dr = nullptr; // r^2*dr[0:n)
    double const* rinv = nullptr; // r^-1
    double anisotropy = 0.; // ToDo: this could be float
    bool  memory_owner = false;
    char  equation = '\0';

#ifdef    __cplusplus
  public:

//   double constexpr default_anisotropy = 0.01;
//   float  constexpr default_Rmax = 9.45;

    radial_grid_t(void) {}; // default constructor

    // radial_grid_t( // constructor
    //     int const npoints // number of grid points
    //   , int const echo=0  // verbosity
    //   , float const rmax=default_Rmax               // [optional] largest radius
    //   , char equation='\0'                          // [optional] how to generate the grid
    //   , double const anisotropy=default_anisotropy  // [optional] anisotropy parameter
    // ); // constructor

    // radial_grid_t( // pseudo grid constructor
    //     radial_grid_t const & tru // radial grid for true quantities (usually goes down to 0)
    //   , int const echo=0         // verbosity
    //   , double const r_min=1e-3 // start radius
    // ); // constructor

    // ~radial_grid_t(void); // destructor

#endif // __cplusplus

}; // radial_grid_t
