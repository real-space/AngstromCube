#pragma once

#include <complex> // std::complex<real_t>

#ifndef HAS_no_MKL
	#define HAS_MKL
#endif


extern "C" {
#ifdef  HAS_MKL
  #include "mkl_lapacke.h" // LAPACK_COL_MAJOR, MKL_INT
#else
    // float64 linear_solve
	void dgesv_(int const *n, int const *nrhs, double a[], int const *lda, 
				int ipiv[], double b[], int const *ldb, int *info);
    // float64 eigenvalues
	void dsyev_(char const *jobz, char const *uplo, int const* n, double a[], int const *lda, 
		        double w[], double work[], int const *lwork, int *info);
	void zheev_(char const *jobz, char const *uplo, int const* n, std::complex<double> a[], int const *lda, 
		        double w[], std::complex<double> work[], int const *lwork, double rwork[], int *info);
    // float64 generalized_eigenvalues
    void dsygv_(int const *itype, char const *jobz, char const *uplo, int const* n, double a[], int const *lda,
    			double b[], int const *ldb, double w[], double work[], int const *lwork, int *info);
    void zhegv_(int const *itype, char const *jobz, char const *uplo, int const* n, std::complex<double> a[], int const *lda,
    			std::complex<double> b[], int const *ldb, double w[],
                std::complex<double> work[], int const *lwork, double rwork[], int *info);
    // float32 eigenvalues
	void ssyev_(char const *jobz, char const *uplo, int const* n, float a[], int const *lda, 
		        float w[], float work[], int const *lwork, int *info);
	void cheev_(char const *jobz, char const *uplo, int const* n, std::complex<float> a[], int const *lda, 
		        float w[], std::complex<float> work[], int const *lwork, float rwork[], int *info);
    // float32 generalized_eigenvalues
    void ssygv_(int const *itype, char const *jobz, char const *uplo, int const* n, float a[], int const *lda,
    			float b[], int const *ldb, float w[], float work[], int const *lwork, int *info);
    void chegv_(int const *itype, char const *jobz, char const *uplo, int const* n, std::complex<float> a[], int const *lda,
    			std::complex<float> b[], int const *ldb, float w[],
                std::complex<float> work[], int const *lwork, float rwork[], int *info);
#endif
} // extern "C"


extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int const *m, int const *n, double a[], int const *lda, int ipiv[], int *info);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int const *n, double a[], int const *lda, int ipiv[], double work[], int const *lwork, int *info);
} // extern "C"


#include "status.hxx" // status_t

namespace linear_algebra {
  
  inline status_t inverse(int const n, double a[], int const lda) {
      int info{0};
      int const lwork = n*n;
      auto const work = new double[lwork];
      auto const ipiv = new int[n];

      dgetrf_(&n, &n, a, &lda, ipiv, &info); // Fortran interface
      if (info) return info; // early return: factorization failed!
      dgetri_(&n, a, &lda, ipiv, work, &lwork, &info); // Fortran interface

      delete[] ipiv;
      delete[] work;
      return info;
  } // inverse
  
  inline status_t linear_solve(int const n, double a[], int const lda, double b[], int const ldb, int const nrhs=1) {
#ifdef  HAS_MKL
      auto const ipiv = new MKL_INT[n];
      status_t const info = LAPACKE_dgesv( LAPACK_COL_MAJOR, n, nrhs, a, lda, ipiv, b, ldb );
#else
      auto const ipiv = new int[2*n];
      int info{0};
      dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info); // Fortran interface
#endif
      delete[] ipiv;
      return info;
  } // linear_solve
  
  
  inline status_t eigenvalues(int const n, double a[], int const lda, double w[]) {
#ifdef  HAS_MKL
      status_t const info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', n, a, lda, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const lwork = (2*n + 2)*n;
      std::vector<double> work(lwork);
      dsyev_(&jobz, &uplo, &n, a, &lda, w, work.data(), &lwork, &info); // Fortran interface
#endif
      return info;
  } // (standard_)eigenvalues

  inline status_t eigenvalues(int const n, std::complex<double> a[], int const lda, double w[]) {
#ifdef  HAS_MKL
      status_t const info = LAPACKE_zheev( LAPACK_COL_MAJOR, 'V', 'U', n, a, lda, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const lwork = (2*n + 2)*n;
      std::vector<std::complex<double>> work(lwork);
      std::vector<double> rwork(3*n);
      zheev_(&jobz, &uplo, &n, a, &lda, w, work.data(), &lwork, rwork.data(), &info); // Fortran interface
#endif
      return info;
  } // (standard_)eigenvalues
  
  inline status_t generalized_eigenvalues(int const n, double a[], int const lda, double b[], int const ldb, double w[]) {
//    printf("\n# call dsygv(1, 'v', 'u', %i, %p, %i, %p, %i, %p)\n\n",   n, a, lda, b, ldb, w);
#ifdef  HAS_MKL
      status_t const info = LAPACKE_dsygv( LAPACK_COL_MAJOR, 1, 'V', 'U', n, a, lda, b, ldb, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const itype = 1, lwork = (2*n + 2)*n;
      std::vector<double> work(lwork);
      dsygv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work.data(), &lwork, &info); // Fortran interface
#endif
      return info;
  } // generalized_eigenvalues

  inline status_t generalized_eigenvalues(int const n, std::complex<double> a[], int const lda, std::complex<double> b[], int const ldb, double w[]) {
//    printf("\n# call zhegv(1, 'v', 'u', %i, %p, %i, %p, %i, %p)\n\n",   n, a, lda, b, ldb, w);
#ifdef  HAS_MKL
      status_t const info = LAPACKE_zhegv( LAPACK_COL_MAJOR, 1, 'V', 'U', n, a, lda, b, ldb, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const itype = 1, lwork = (2*n + 2)*n;
      std::vector<std::complex<double>> work(lwork);
      std::vector<double> rwork(3*n);
      zhegv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work.data(), &lwork, rwork.data(), &info); // Fortran interface
#endif
      return info;
  } // generalized_eigenvalues

  
  // float32 version: copy the above, replace double by float and change d --> s, z --> c

  
  inline status_t eigenvalues(int const n, float a[], int const lda, float w[]) {
#ifdef  HAS_MKL
      status_t const info = LAPACKE_ssyev( LAPACK_COL_MAJOR, 'V', 'U', n, a, lda, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const lwork = (2*n + 2)*n;
      std::vector<float> work(lwork);
      ssyev_(&jobz, &uplo, &n, a, &lda, w, work.data(), &lwork, &info); // Fortran interface
#endif
      return info;
  } // (standard_)eigenvalues

  inline status_t eigenvalues(int const n, std::complex<float> a[], int const lda, float w[]) {
#ifdef  HAS_MKL
      status_t const info = LAPACKE_cheev( LAPACK_COL_MAJOR, 'V', 'U', n, a, lda, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const lwork = (2*n + 2)*n;
      std::vector<std::complex<float>> work(lwork);
      std::vector<float> rwork(3*n);
      cheev_(&jobz, &uplo, &n, a, &lda, w, work.data(), &lwork, rwork.data(), &info); // Fortran interface
#endif
      return info;
  } // (standard_)eigenvalues
  
  inline status_t generalized_eigenvalues(int const n, float a[], int const lda, float b[], int const ldb, float w[]) {
//    printf("\n# call ssygv(1, 'v', 'u', %i, %p, %i, %p, %i, %p)\n\n",   n, a, lda, b, ldb, w);
#ifdef  HAS_MKL
      status_t const info = LAPACKE_ssygv( LAPACK_COL_MAJOR, 1, 'V', 'U', n, a, lda, b, ldb, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const itype = 1, lwork = (2*n + 2)*n;
      std::vector<float> work(lwork);
      ssygv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work.data(), &lwork, &info); // Fortran interface
#endif
      return info;
  } // generalized_eigenvalues

  inline status_t generalized_eigenvalues(int const n, std::complex<float> a[], int const lda, std::complex<float> b[], int const ldb, float w[]) {
//    printf("\n# call chegv(1, 'v', 'u', %i, %p, %i, %p, %i, %p)\n\n",   n, a, lda, b, ldb, w);
#ifdef  HAS_MKL
      status_t const info = LAPACKE_chegv( LAPACK_COL_MAJOR, 1, 'V', 'U', n, a, lda, b, ldb, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const itype = 1, lwork = (2*n + 2)*n;
      std::vector<std::complex<float>> work(lwork);
      std::vector<float> rwork(3*n);
      chegv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work.data(), &lwork, rwork.data(), &info); // Fortran interface
#endif
      return info;
  } // generalized_eigenvalues

  
  inline status_t all_tests(int const echo=0) { return 0; }

} // namespace linear_algebra
