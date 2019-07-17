#pragma once

#include <cstdint> // uint32_t
#include <cstdio> // printf

typedef int status_t;

namespace real_space_grid {

  template<typename real_t, int D0> // D0: inner dimension, vector length
  class grid_t {
  public:
      real_t* values;
  private:
      uint32_t dims[4]; // 0,1,2:real-space grid dimensions, 3:outer dim
  public:
      double h[3], inv_h[3]; // grid spacings and their inverse

      grid_t(int const dim[3], int const dim_outer=1, real_t* const v=nullptr)
       : h{1,1,1}, inv_h{1,1,1} {
          dims[0] = std::max(1, dim[0]); // x
          dims[1] = std::max(1, dim[1]); // y
          dims[2] = std::max(1, dim[2]); // z
          dims[3] = std::max(1, dim_outer);
          long const nnumbers = dims[3] * dims[2] * dims[1] * dims[0] * D0;
          if (nnumbers > 0) {
              if (nullptr == v) {
                  printf("# allocate a grid with %d * %d x %d x %d * %d real_%ld = %.6f GByte\n", 
                      dims[3], dims[2], dims[1], dims[0], D0, sizeof(real_t), nnumbers*1e-9*sizeof(real_t));
                  values = new real_t[nnumbers];
              } else {
                  printf("# use %p as grid with %d * %d x %d x %d * %d real_%ld\n", 
                      v, dims[3], dims[2], dims[1], dims[0], D0, sizeof(real_t));
                  values = v;
              }
          } else {
              printf("# grid invalid: <D0=%d> dims={%d, %d, %d,  %d}\n", 
                      D0, dims[0], dims[1], dims[2], dims[3]);
              values = nullptr;
          }
      } // constructor

      ~grid_t() {
          delete [] values;
          for(int i4 = 0; i4 < 4; ++i4) dims[i4] = 0;
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
      
      inline int dim(char const xyz) const { return ('w' == (xyz | 32)) ? D0 : dims[(xyz | 32) - 120]; }
      inline int dim(int const x0y1z2) const { return dims[x0y1z2]; }
      inline double dV(bool const Cartesian=true) const { return h[0]*h[1]*h[2]; } // volume element, assuming a Cartesian grid

      inline size_t all() const { return dims[3] * dims[2] * dims[1] * dims[0] * D0; }
  };
  
  status_t all_tests();

} // namespace real_space_grid
