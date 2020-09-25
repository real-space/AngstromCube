#pragma once

#include <cstdint> // uint32_t
#include <cstdio> // printf
#include <vector> // std::vector<T>

#include "real_space.hxx" // grid_t
#include "boundary_condition.hxx" // *_Boundary
#include "recorded_warnings.hxx" // warn

#include "status.hxx" // status_t

namespace finite_difference {

  int constexpr nnArraySize = 16;
  
  template<typename real_t>
  int set_Laplacian_coefficients(real_t c[], // c[1+nn]
            int const nn=1, // returns nn on success
            double const grid_spacing=1,
            char const direction='?') {
    double const h2 = grid_spacing*grid_spacing;
    int constexpr nnMaxImplemented = 13;
    switch (std::min(nn, nnMaxImplemented)) {
        // regular cases first
      
        case 1: // use 3 points, lowest order
            c[1] = 1./(1.*h2);
            c[0] = -2./(1.*h2);
        break;

        case 2: // use 5 points
            c[2] = -1./(12.*h2);
            c[1] = 4./(3.*h2);
            c[0] = -5./(2.*h2);
        break;

        case 3: // use 7 points
            c[3] = 1./(90.*h2);
            c[2] = -3./(20.*h2);
            c[1] = 3./(2.*h2);
            c[0] = -49./(18.*h2);
        break;

        case 4: // use 9 points
            c[4] = -1./(560.*h2);
            c[3] = 8./(315.*h2);
            c[2] = -1./(5.*h2);
            c[1] = 8./(5.*h2);
            c[0] = -205./(72.*h2);
        break;

        case 5: // use 11 points
            c[5] = 1./(3150.*h2);
            c[4] = -5./(1008.*h2);
            c[3] = 5./(126.*h2);
            c[2] = -5./(21.*h2);
            c[1] = 5./(3.*h2);
            c[0] = -5269./(1800.*h2);
        break;

        case 6: // use 13 points
            c[6] = -1./(16632.*h2);
            c[5] = 2./(1925.*h2);
            c[4] = -1./(112.*h2);
            c[3] = 10./(189.*h2);
            c[2] = -15./(56.*h2);
            c[1] = 12./(7.*h2);
            c[0] = -5369./(1800.*h2);
        break;

        case 7: // use 15 points
            c[7] = 1./(84084.*h2);
            c[6] = -7./(30888.*h2);
            c[5] = 7./(3300.*h2);
            c[4] = -7./(528.*h2);
            c[3] = 7./(108.*h2);
            c[2] = -7./(24.*h2);
            c[1] = 7./(4.*h2);
            c[0] = -266681./(88200.*h2);
        break;

        case 8: // use 17 points
            c[8] = -1./(411840.*h2);
            c[7] = 16./(315315.*h2);
            c[6] = -2./(3861.*h2);
            c[5] = 112./(32175.*h2);
            c[4] = -7./(396.*h2);
            c[3] = 112./(1485.*h2);
            c[2] = -14./(45.*h2);
            c[1] = 16./(9.*h2);
            c[0] = -1077749./(352800.*h2);
        break;

        case 9: // use 19 points
            c[9] = 1./(1969110.*h2);
            c[8] = -9./(777920.*h2);
            c[7] = 9./(70070.*h2);
            c[6] = -2./(2145.*h2);
            c[5] = 18./(3575.*h2);
            c[4] = -63./(2860.*h2);
            c[3] = 14./(165.*h2);
            c[2] = -18./(55.*h2);
            c[1] = 9./(5.*h2);
            c[0] = -9778141./(3175200.*h2);
        break;

        case 10: // use 21 points
            c[10] = -1./(9237800.*h2);
            c[9] = 10./(3741309.*h2);
            c[8] = -5./(155584.*h2);
            c[7] = 30./(119119.*h2);
            c[6] = -5./(3432.*h2);
            c[5] = 24./(3575.*h2);
            c[4] = -15./(572.*h2);
            c[3] = 40./(429.*h2);
            c[2] = -15./(44.*h2);
            c[1] = 20./(11.*h2);
            c[0] = -1968329./(635040.*h2);
        break;

        case 11: // use 23 points
            c[11] = 1./(42678636.*h2);
            c[10] = -11./(17635800.*h2);
            c[9] = 11./(1360476.*h2);
            c[8] = -55./(806208.*h2);
            c[7] = 55./(129948.*h2);
            c[6] = -11./(5304.*h2);
            c[5] = 11./(1300.*h2);
            c[4] = -11./(364.*h2);
            c[3] = 55./(546.*h2);
            c[2] = -55./(156.*h2);
            c[1] = 11./(6.*h2);
            c[0] = -239437889./(76839840.*h2);
        break;

        case 12: // use 25 points
            c[12] = -1./(194699232.*h2);
            c[11] = 12./(81800719.*h2);
            c[10] = -3./(1469650.*h2);
            c[9] = 44./(2380833.*h2);
            c[8] = -33./(268736.*h2);
            c[7] = 132./(205751.*h2);
            c[6] = -11./(3978.*h2);
            c[5] = 396./(38675.*h2);
            c[4] = -99./(2912.*h2);
            c[3] = 88./(819.*h2);
            c[2] = -33./(91.*h2);
            c[1] = 24./(13.*h2);
            c[0] = -240505109./(76839840.*h2);
        break;

        case 13: // use 27 points
            c[13] = 1./(878850700.*h2);
            c[12] = -13./(374421600.*h2);
            c[11] = 13./(25169452.*h2);
            c[10] = -13./(2600150.*h2);
            c[9] = 13./(366282.*h2);
            c[8] = -143./(723520.*h2);
            c[7] = 143./(158270.*h2);
            c[6] = -143./(40698.*h2);
            c[5] = 143./(11900.*h2);
            c[4] = -143./(3808.*h2);
            c[3] = 143./(1260.*h2);
            c[2] = -13./(35.*h2);
            c[1] = 13./(7.*h2);
            c[0] = -40799043101./(12985932960.*h2);
        break;

        // treat irregular cases
        case -1:
            warn("Finite difference no pass in %c-direction", direction);
        break;

        case 0:
            warn("Finite difference in %c-direction switched off!", direction);
            c[0] = 0;
        break; 
        
        default:
            warn("Cannot treat case of %i finite difference neighbors in %c-direction", nn, direction);
            return -1; // failed, does not return nn
    } // switch min(nn, nnMaxImplemented)

    if (nn > nnMaxImplemented) {
        warn("Max implemented %i but requested %i finite difference neighbors in %c-direction", nnMaxImplemented, nn, direction);
        for(int i = nnMaxImplemented + 1; i < std::min(nn + 1, nnArraySize); ++i) c[i] = 0; // clear out the others.
        return nnMaxImplemented;
    } // larger than implemented

    return nn;
  } // set_Laplacian_coefficients

  
  template<typename real_t> // real_t may be float or double
  class stencil_t {
    public:
      real_t c2nd[3][nnArraySize]; // coefficients for the 2nd derivative
    private:
      int8_t _nn[3]; // number of FD neighbors
      
      void _constructor(double const grid_spacing[3], int const nneighbors[3], double const scale_factor) {
          for(int d = 0; d < 3; ++d) {
              for(int i = 0; i < nnArraySize; ++i) c2nd[d][i] = 0; // clear
              double const h = grid_spacing[d];
              _nn[d] = set_Laplacian_coefficients(c2nd[d], nneighbors[d], h, 'x'+d);
              if (_nn[d] < nneighbors[d]) {
                  warn("In stencil_t requested nn=%i but use nn=%i for %c-direction", 
                                  nneighbors[d], _nn[d], 'x'+d);
              }
          } // spatial direction d
          scale_coefficients(scale_factor);
      } // _constructor
      
    public:

      stencil_t(double const grid_spacing[3], int const nneighbors[3], double const scale_factor=1) {
          _constructor(grid_spacing, nneighbors, scale_factor);
      } // preferred constructor

      stencil_t(double const grid_spacing[3], int const nn=4, double const scale_factor=1) {
          int const nns[3] = {nn, nn, nn};
          _constructor(grid_spacing, nns, scale_factor);
      } // isotropic nn constructor

      stencil_t(double const h=1, int const nn=4, double const scale_factor=1) {
          double const hgs[3] = {h, h, h};
          int const nns[3] = {nn, nn, nn};
          _constructor(hgs, nns, scale_factor);
      } // isotropic constructor, default constructor

      double clear_diagonal_elements() { // modifies the coefficients c2nd[][]
          double diag{0};
          for(int d = 0; d < 3; ++d) {
              diag += c2nd[d][0];
              c2nd[d][0] = 0; // clear diagonal elements
          } // d
          return diag;
      } // clear_diagonal_elements

  private:
      void scale_coefficients(double const f[3]) {
          for(int d = 0; d < 3; ++d) {
              for(int i = 0; i < nnArraySize; ++i) {
                  c2nd[d][i] *= f[d];
              } // i
          } // d
      } // scale_coefficients
      
      void scale_coefficients(double const f) { double const f3[] = {f, f, f}; scale_coefficients(f3); }
  public:
      int8_t const * nearest_neighbors() const { return _nn; }
      int nearest_neighbors(int const d) const { assert(d >= 0); assert(d < 3); return _nn[d]; }
      real_t * Laplace_coefficients(int const d) { assert(d >= 0); assert(d < 3); return c2nd[d]; }
      
  }; // class stencil_t

  template <typename real_out_t, // result is stored in this precision 
            typename real_in_t, // input comes in this precision
            typename real_fd_t> // computations are executed in this precision
  status_t apply(real_out_t out[], real_in_t const in[], real_space::grid_t const &g, 
                     stencil_t<real_fd_t> const &fd, double const factor=1) {

      int const n16 = nnArraySize; // max number of finite difference neighbors
      std::vector<int> list[3]; // could be of type int16_t, needs assert(n < (1 << 15));
      for(int d = 0; d < 3; ++d) {
          int const n  = g[d];
          int const bc = g.boundary_condition(d);
          int const nf = fd.nearest_neighbors(d);
          int const nh = n16 + n + n16; // number including largest halos
          list[d] = std::vector<int>(nh, -1); // get memory, init as non-existing
          for(int i = 0; i < n; ++i) list[d][n16 + i] = i; // core region

          // lower boundary
          if (Periodic_Boundary == bc) { // periodic BC
              for(int i = n16 - nf; i < n16; ++i) list[d][i] = n + i - n16; // wrap around
          } else if (Mirrored_Boundary == bc) { // mirror BC
              for(int i = n16 - nf; i < n16; ++i) list[d][i] = n16 - 1 - i; // wrap around and mirror
          } // else open BC, list[:] = -1
          
          // upper boundary
          if (Periodic_Boundary == bc) { // periodic BC
              for(int i = 0; i < nf; ++i) list[d][n16 + n + i] = i; // wrap around
          } else if (Mirrored_Boundary == bc) { // mirror BC
              for(int i = 0; i < nf; ++i) list[d][n16 + n + i] = n - 1 - i; // wrap around and mirror
          } // else open BC, list[:] = -1
  
          if (0) { // DEBUG: show indirection list
              printf("# indirection list for %c  ", 120+d);
              for(int i = n16 - nf; i < n16 + n + nf; ++i) {
                  if ((n16 == i) || (n16 + n) == i) printf(" |");
                  printf(" %i", list[d][i]);
              }   printf("\n");
          } // show indirection list

      } // spatial direction d

      real_fd_t const scale_factor = factor;
      for(int z = 0; z < g('z'); ++z) {
          for(int y = 0; y < g('y'); ++y) {
              for(int x = 0; x < g('x'); ++x) {
                  int const i_zyx = (z*g('y') + y)*g('x') + x;

                  real_fd_t t(0); // init result

                  for(int d = 0; d < 3; ++d) {
                      int const nf = fd.nearest_neighbors(d);
                      int zyx[3] = {x, y, z};
                      int const i_center = zyx[d];
                      for(int jmi = -nf; jmi <= nf; ++jmi) {
                          int const j = i_center + jmi;
                          int const index = list[d][n16 + j];
                          if (index >= 0) {
                              zyx[d] = index;
                              int const zyx_prime = (zyx[2]*g('y') + zyx[1])*g('x') + zyx[0];
                              auto const coeff = fd.c2nd[d][std::abs(jmi)];
                              t += in[zyx_prime] * coeff;
                          } // index exists
                      } // jmi
                  } // d direction of the derivative

                  out[i_zyx] = t * scale_factor; // store

              } // x
          } // y
      } // z

      return 0; // success
  } // apply
  
  status_t all_tests(int const echo=0);

} // namespace finite_difference
