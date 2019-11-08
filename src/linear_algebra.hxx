#pragma once

extern "C" {
  #include "mkl_lapacke.h" // LAPACK_COL_MAJOR, MKL_INT
} // extern "C"

typedef int status_t;

namespace linear_algebra {

  inline status_t linear_solve(int const n, double a[], int const lda, double b[], int const ldb, int const nrhs=1) {
      auto const ipiv = new MKL_INT[n];
      status_t const info = LAPACKE_dgesv( LAPACK_COL_MAJOR, n, nrhs, a, lda, ipiv, b, ldb );
      delete[] ipiv;
      return info;
  } // linear_solve

  inline status_t eigenvalues(int const n, double a[], int const lda, double w[]) {
      status_t const info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', n, a, lda, w );
      return info;
  } // (standard_)eigenvalues

  inline status_t generalized_eigenvalues(int const n, double a[], int const lda, double b[], int const ldb, double w[]) {
//    printf("\n# call dsygv(1, 'v', 'u', %i, %p, %i, %p, %i, %p)\n\n",   n, a, lda, b, ldb, w);
      status_t const info = LAPACKE_dsygv( LAPACK_COL_MAJOR, 1, 'V', 'U', n, a, lda, b, ldb, w );
      return info;
  } // generalized_eigenvalues

  inline status_t all_tests() { return 0; }

} // namespace linear_algebra
