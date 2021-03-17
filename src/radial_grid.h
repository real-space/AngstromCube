#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float rmax = 0.f; // max radius
    int   n = 0; // number of points
    double*     r; // r[0:n)
    double*    dr; // dr[0:n)
    double*   rdr; // r*dr[0:n)
    double*  r2dr; // r^2*dr[0:n)
    double*  rinv;  // r^-1
    float anisotropy = 0.f;
    bool  memory_owner;
    char const* equation;
} radial_grid_t;

#ifdef __cplusplus
} // extern "C"
#endif
