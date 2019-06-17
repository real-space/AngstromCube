#pragma once

typedef struct {
    double* Xlm2grid;
    double* grid2Xlm;
    double (*xyzw)[4];
    int Xlm2grid_stride;
    int grid2Xlm_stride;
    int npoints;
    int ellmax;
} angular_grid_t;
