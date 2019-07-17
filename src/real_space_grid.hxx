#pragma once

#include <cstdint> // uint32_t
#include <cstdio> // printf

typedef int status_t;

namespace real_space_grid {

  class real_space_grid_t {
  public:
      double* values;
  private:
      uint32_t dims[6]; // 0:inner dim, 1,2,3:real-space grid dimensions, 4:outer dim, 5:not used
  public:
      double h[3], inv_h[3]; // grid spacings and their inverse
      
      real_space_grid_t(int const d123[3], int const d0=1, int const d4=1, double* const v=nullptr)
       : h{1.,1.,1.}, inv_h{1.,1.,1.} {
          dims[0] = d0; // inner
          dims[1] = d123[0]; // x
          dims[2] = d123[1]; // y
          dims[3] = d123[2]; // z
          dims[4] = d4; // outer
          dims[5] = 0; // not used
          long const nnumbers = 1 * d4 * d123[2] * d123[1] * d123[0] * d0;
          if (nnumbers > 0) {
              if (nullptr == v) {
                  printf("# allocate a grid with %d * %d x %d x %d * %d doubles = %.6f GByte\n", 
                      dims[4], dims[3], dims[2], dims[1], dims[0], nnumbers*1e-9*sizeof(double));
                  values = new double[nnumbers];
              } else {
                  printf("# use %p as grid with %d * %d x %d x %d * %d doubles\n", 
                      v, dims[4], dims[3], dims[2], dims[1], dims[0]);
                  values = v;
              }
          } else {
              printf("# grid invalid: d0=%d d123={%d, %d, %d} d4=%d\n", 
                  d0, d123[0], d123[1], d123[2], d4);
              values = nullptr;
          }
      } // constructor

      ~real_space_grid_t() {
          delete [] values;
          for(int i6 = 0; i6 < 6; ++i6) dims[i6] = 0;
      } // destructor
      
      status_t set_grid_spacing(float const hx, float const hy=-1, float const hz=-1) {
          status_t stat = 0;
          float const h3[] = {hx, (hy<0)?hx:hy, (hz<0)?hx:hy};
          for(int i3 = 0; i3 < 3; ++i3) {
              h[i3] = h3[i3]; // convert to double
              if ((0.0 != h[i3]) && (h[i3] == h[i3])) {
                  inv_h[i3] = 1./h[i3]; // invert only here
              } else ++stat;
          } // i3
          return stat;
      } // set
      
      inline int dim(char const xyz) { return dims[(xyz | 32) - 119]; }
      inline int dim(int const x0y1z2) { return dims[1 + x0y1z2]; }
      inline double dV() { return h[0]*h[1]*h[2]; } // volume element
      
      inline size_t all() { return 1 * dims[4] * dims[3] * dims[2] * dims[1] * dims[0]; }
  };
  
  status_t all_tests();

} // namespace real_space_grid
