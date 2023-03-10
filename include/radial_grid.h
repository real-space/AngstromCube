#pragma once
// This file is part of AngstromCube under MIT License

// #ifdef __cplusplus
// extern "C" {
// #endif

struct radial_grid_t {
    int   n = 0; // number of points
    float rmax = 0.f; // max radius
    double const*    r = nullptr; // r[0:n)
    double const*   dr = nullptr; // dr[0:n)
    double const*  rdr = nullptr; // r*dr[0:n)
    double const* r2dr = nullptr; // r^2*dr[0:n)
    double const* rinv = nullptr; // r^-1
    double anisotropy = 0.;
    bool  memory_owner = true;
    char  equation = '\0';
}; // radial_grid_t

// #ifdef __cplusplus
// } // extern "C"
// #endif
