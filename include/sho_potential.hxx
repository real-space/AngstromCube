#pragma once
// This file is part of AngstromCube under MIT License

#include <cassert> // assert
#include <vector> // std::vector<T>

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D<T>, view4D<T>
#include "sho_tools.hxx" // ::nSHO, ::n1HO

namespace sho_potential {
  // computes potential matrix elements between two SHO basis functions

  status_t load_local_potential(
        std::vector<double> & vtot // output
      , int dims[3] // output dimensions found
      , char const *filename // input filename
      , int const echo=0 // log-level
  ); // declaration only

  status_t normalize_potential_coefficients(
        double coeff[] // coefficients[nSHO(numax)], input: in zyx_order, output in Ezyx_order
      , int const numax // SHO basis size
      , double const sigma // SHO basis spread
      , int const echo // log-level
  ); // declaration only

  template <typename complex_t, typename phase_t>
  status_t potential_matrix(
        view2D<complex_t> & Vmat // result Vmat(i,j) += sum_m Vcoeff[m] * t(m,i,j) 
      , view4D<double> const & t1D // input t1D(dir,m,i,j)
      , double const Vcoeff[] // coefficients of the SHO-expanded local potential (Ezyx_order)
      , int const numax_m // order of expansion of the potential in x^{mx}*y^{my}*z^{mz}
      , int const numax_i
      , int const numax_j
      , phase_t const phase=1 // typically phase_t is double or complex<double>
      , int const dir01=1 // 1:direction dependent input tensor, 0:isotropic
  ) {
      // use the expansion of the product of two Hermite Gauss functions into another one, factorized in 3D
      // can use different (dir01==1) tensors per direction or the same (dir01==0, isotropic)
      // t(m,j,i) = t1D(0,m_x,j_x,i_x) * t1D(1,m_y,j_y,i_y) * t1D(2,m_z,j_z,i_z)
      assert( t1D.stride()  >= sho_tools::n1HO(numax_j) );
      assert( t1D.dim1()    >= sho_tools::n1HO(numax_i) );
      assert( t1D.dim2()    >= sho_tools::n1HO(numax_m) );
      assert( Vmat.stride() >= sho_tools::nSHO(numax_j) );

      int mzyx{0}; // contraction index
      for     (int mu = 0; mu <= numax_m; ++mu) { // shell index for order_Ezyx
        for   (int mz = 0; mz <= mu;      ++mz) {
          for (int mx = 0; mx <= mu - mz; ++mx) {
            int const my = mu - mz - mx;

            auto const phase_Vcoeff_m = phase * Vcoeff[mzyx];

            int izyx{0};
            for     (int iz = 0; iz <= numax_i;           ++iz) {
              for   (int iy = 0; iy <= numax_i - iz;      ++iy) {
                for (int ix = 0; ix <= numax_i - iz - iy; ++ix) {

                  int jzyx{0};
                  for     (int jz = 0; jz <= numax_j;           ++jz) {  auto const tz   = t1D(dir01*2,mz,iz,jz);
                    for   (int jy = 0; jy <= numax_j - jz;      ++jy) {  auto const tyz  = t1D(dir01*1,my,iy,jy) * tz;
                      for (int jx = 0; jx <= numax_j - jz - jy; ++jx) {  auto const txyz = t1D(      0,mx,ix,jx) * tyz;

                        Vmat(izyx,jzyx) += phase_Vcoeff_m * txyz;

                        ++jzyx;
                      } // jx
                    } // jy
                  } // jz
                  assert( sho_tools::nSHO(numax_j) == jzyx );

                  ++izyx;
                } // ix
              } // iy
            } // iz
            assert( sho_tools::nSHO(numax_i) == izyx );

            ++mzyx;
          } // mx
        } // my
      } // mz
      assert( sho_tools::nSHO(numax_m) == mzyx );

      return 0;
  } // potential_matrix

  status_t all_tests(int const echo=0); // declaration only

} // namespace sho_potential
