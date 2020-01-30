#pragma once

typedef int status_t;

namespace sho_overlap {

  template<typename real_t>
  status_t generate_overlap_matrix(real_t matrix[], // matrix layout [n1][n0]
                     double const distance,
                     int const n0, int const n1, 
                     double const sigma0=1, double const sigma1=1);

  template<typename real_t>
  status_t moment_tensor(real_t tensor[], // matrix layout [][n1][n0]
                     double const distance,
                     int const n0, int const n1,
                     double const sigma0=1, double const sigma1=1, int const maxmoment=0);

  template<typename real_t>
  status_t generate_product_tensor(real_t tensor[], int const n, // tensor layout [2*n][n][n]
                     double const sigma=2, // 2:typical for density tensor
                     double const sigma0=1, double const sigma1=1);

  status_t all_tests(int const echo=0);

} // namespace sho_overlap
