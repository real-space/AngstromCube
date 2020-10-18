#pragma once

#include "status.hxx" // status_t

namespace sho_overlap {
  
  template <typename real_t>
  status_t moment_tensor(
        real_t tensor[] // tensor layout [1 + maxmoment][n1][n0]
      , double const distance
      , int const n1
      , int const n0
      , double const sigma1=1.0
      , double const sigma0=1.0
      , int const maxmoment=0
  ); // declaration only

  template <typename real_t> inline
  status_t overlap_matrix(
        real_t matrix[] // matrix layout [n1][n0]
      , double const distance
      , int const n1
      , int const n0
      , double const sigma1 // =1.0
      , double const sigma0 // =1.0
  ) {
      return moment_tensor(matrix, distance, n1, n0, sigma1, sigma0, 0);
  } // overlap_matrix

  template <typename real_t> inline // ToDo: rename 1 <--> 0
  status_t nabla2_matrix(
        real_t matrix[] // matrix layout [n1][n0]
      , double const distance
      , int const n1
      , int const n0
      , double const sigma1 // =1.0
      , double const sigma0 // =1.0
  ) {
      return moment_tensor(matrix, distance, n1, n0, sigma1, sigma0, -2);
  } // nabla2_matrix

  status_t moment_normalization(
        double matrix[] // matrix layout [m][m]
      , int const m
      , double const sigma=1.0
      , int const echo=0
  ); // declaration only

  template <typename real_t>
  status_t product_tensor(
        real_t tensor[]
      , int const n // tensor layout [2*n-1][n][n]
      , double const sigma=2.0 // 2:typical for density tensor (REALLY?)
      , double const sigma1=1.0
      , double const sigma0=1.0
  ); // declaration only

  status_t all_tests(int const echo=0);

} // namespace sho_overlap
