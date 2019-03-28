#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int   n = 0; // numer of points
    float rmax = 0.; // max radius
    double*     r; // r[0:n)
    double*    dr; // dr[0:n)
    double*   rdr; // r*dr[0:n)
    double*  r2dr; // r^2*dr[0:n)
    double*  rinv;  // r^-1
} radial_grid_t;

#ifdef __cplusplus
} // extern "C"
#endif
