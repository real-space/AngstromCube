#pragma once

#ifndef	HAS_no_MKL
	#define HAS_MKL
#endif


extern "C" {
#ifdef 	HAS_MKL
  #include "mkl_lapacke.h" // LAPACK_COL_MAJOR, MKL_INT
#else
	void dgesv_(int const *n, int const *nrhs, double a[], int const *lda, 
				int ipiv[], double b[], int const *ldb, int *info);
	void dsyev_(char const *jobz, char const *uplo, int const* n, double a[], int const *lda, 
		        double w[], double work[], int const *lwork, int *info);
    void dsygv_(int const *itype, char const *jobz, char const *uplo, int const* n, double a[], int const *lda,
    			double b[], int const *ldb, double w[], double work[], int const *lwork, int *info);
#endif
} // extern "C"

typedef int status_t;

namespace linear_algebra {

  inline status_t linear_solve(int const n, double a[], int const lda, double b[], int const ldb, int const nrhs=1) {
#ifdef 	HAS_MKL
      auto const ipiv = new MKL_INT[n];
      status_t const info = LAPACKE_dgesv( LAPACK_COL_MAJOR, n, nrhs, a, lda, ipiv, b, ldb );
#else
      auto const ipiv = new int[2*n];
      int info = 0;
      dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info); // Fortran interface
#endif
      delete[] ipiv;
      return info;
  } // linear_solve

  inline status_t eigenvalues(int const n, double a[], int const lda, double w[]) {
#ifdef 	HAS_MKL
      status_t const info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', n, a, lda, w );
#else
      int info = 0; char const jobz = 'V', uplo = 'U'; int const lwork = (2*n + 2)*n;
      auto const work = new double[lwork];
      dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info); // Fortran interface
      delete[] work;
#endif
      return info;
  } // (standard_)eigenvalues

  inline status_t generalized_eigenvalues(int const n, double a[], int const lda, double b[], int const ldb, double w[]) {
//    printf("\n# call dsygv(1, 'v', 'u', %i, %p, %i, %p, %i, %p)\n\n",   n, a, lda, b, ldb, w);
#ifdef 	HAS_MKL
      status_t const info = LAPACKE_dsygv( LAPACK_COL_MAJOR, 1, 'V', 'U', n, a, lda, b, ldb, w );
#else
      int info = 0; char const jobz = 'V', uplo = 'U'; int const itype = 1, lwork = (2*n + 2)*n;
      auto const work = new double[lwork];
      dsygv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, &info); // Fortran interface
      delete[] work;
#endif      
      return info;
  } // generalized_eigenvalues

  inline status_t all_tests() { return 0; }

} // namespace linear_algebra
