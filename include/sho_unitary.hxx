#pragma once
// This file is part of AngstromCube under MIT License

#include "status.hxx" // status_t
#include "sho_tools.hxx" // ::SHO_order_t, ::order_*

namespace sho_unitary {

  class Unitary_SHO_Transform {
  public:

      Unitary_SHO_Transform(int const lmax=7, int const echo=8); // declaration only

      ~Unitary_SHO_Transform() {
          for (int nu = 0; nu <= numax_; ++nu) {
              delete[] u_[nu];
          } // nu
          delete[] u_;
      } // destructor

      double get_entry(int const nzyx, int const nlnm) const; // declaration only

      status_t construct_dense_matrix(
            double matrix[]
          , int const nu_max
          , int const matrix_stride=-1
          , sho_tools::SHO_order_t const row_order=sho_tools::order_Ezyx // energy-ordered Cartesian
          , sho_tools::SHO_order_t const col_order=sho_tools::order_Elnm // energy-ordered radial
      ) const; // declaration only

      status_t transform_vector(double out[], sho_tools::SHO_order_t const out_order
                        , double const inp[], sho_tools::SHO_order_t const inp_order
                        , int const nu_max, int const echo=0) const; // declaration only

      double test_unitarity(int const echo=7) const; // declaration only

      inline int numax() const { return numax_; }

  private: // members
      double **u_; // block diagonal matrix entries
      int numax_;  // largest ell

  }; // class Unitary_SHO_Transform

  status_t all_tests(int const echo=0); // declaration only

} // namespace sho_unitary
