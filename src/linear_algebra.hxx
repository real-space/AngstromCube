#pragma once

#include <complex> // std::complex<real_t>
#include <vector> // std::vector<T>

#ifndef HAS_no_MKL
	#define HAS_MKL
#endif


extern "C" {
#ifdef  HAS_MKL
  #include "mkl_lapacke.h" // LAPACK_COL_MAJOR, MKL_INT, ...
#else
  #define MKL_INT int

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
    
    // LU decomoposition of a general matrix
    void dgetrf_(int const *m, int const *n, double a[], int const *lda, int ipiv[], int *info);
    void sgetrf_(int const *m, int const *n, float  a[], int const *lda, int ipiv[], int *info);
    void zgetrf_(int const *m, int const *n, std::complex<double> a[], int const *lda, int ipiv[], int *info);
    void cgetrf_(int const *m, int const *n, std::complex<float>  a[], int const *lda, int ipiv[], int *info);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int const *n, double a[], int const *lda, int ipiv[], double work[], int const *lwork, int *info);
    void sgetri_(int const *n, float  a[], int const *lda, int ipiv[], float  work[], int const *lwork, int *info);
    void zgetri_(int const *n, std::complex<double> a[], int const *lda, int ipiv[], std::complex<double> work[], int const *lwork, int *info);
    void cgetri_(int const *n, std::complex<float>  a[], int const *lda, int ipiv[], std::complex<float>  work[], int const *lwork, int *info);
#endif
    
    // BLAS interface to matrix matrix multiplication
    void dgemm_(const char*, const char*, const int*, const int*, const int*, const double*,
                const double*, const int*, const double*, const int*, const double*, double*, const int*);
    void zgemm_(const char*, const char*, const int*, const int*, const int*, const std::complex<double>*,
                const std::complex<double>*, const int*, const std::complex<double>*, const int*, const std::complex<double>*, std::complex<double>*, const int*);
    void sgemm_(const char*, const char*, const int*, const int*, const int*, const float*,
                const float*, const int*, const float*, const int*, const float*, float*, const int*);
    void cgemm_(const char*, const char*, const int*, const int*, const int*, const std::complex<float>*,
                const std::complex<float>*, const int*, const std::complex<float>*, const int*, const std::complex<float>*, std::complex<float>*, const int*);
} // extern "C"



#include "status.hxx" // status_t

namespace linear_algebra {
  
  inline status_t inverse(int const n, double a[], int const lda) {
      int info{0};
      std::vector<MKL_INT> ipiv(2*n);
#ifdef  HAS_MKL
      info = LAPACKE_dgetrf( LAPACK_COL_MAJOR, n, n, a, lda, ipiv.data() );
#else
      dgetrf_(&n, &n, a, &lda, ipiv.data(), &info); // Fortran interface
#endif
      if (info) return info; // early return: factorization failed!

#ifdef  HAS_MKL
      info = LAPACKE_dgetri( LAPACK_COL_MAJOR, n, a, lda, ipiv.data() );
#else
      int const lwork = n*n;
      std::vector<double> work(lwork);
      dgetri_(&n, a, &lda, ipiv.data(), work.data(), &lwork, &info); // Fortran interface
#endif
      return info;
  } // inverse
  
  inline status_t inverse(int const n, float a[], int const lda) {
      int info{0};
      std::vector<MKL_INT> ipiv(2*n);
#ifdef  HAS_MKL
      info = LAPACKE_sgetrf( LAPACK_COL_MAJOR, n, n, a, lda, ipiv.data() );
#else
      sgetrf_(&n, &n, a, &lda, ipiv.data(), &info); // Fortran interface
#endif
      if (info) return info; // early return: factorization failed!

#ifdef  HAS_MKL
      info = LAPACKE_sgetri( LAPACK_COL_MAJOR, n, a, lda, ipiv.data() );
#else
      int const lwork = n*n;
      std::vector<float> work(lwork);
      sgetri_(&n, a, &lda, ipiv.data(), work.data(), &lwork, &info); // Fortran interface
#endif
      return info;
  } // inverse

  inline status_t inverse(int const n, std::complex<double> a[], int const lda) {
      int info{0};
      std::vector<MKL_INT> ipiv(2*n);
#ifdef  HAS_MKL
      info = LAPACKE_zgetrf( LAPACK_COL_MAJOR, n, n, (MKL_Complex16*)a, lda, ipiv.data() );
#else
      zgetrf_(&n, &n, a, &lda, ipiv.data(), &info); // Fortran interface
#endif
      if (info) return info; // early return: factorization failed!

#ifdef  HAS_MKL
      info = LAPACKE_zgetri( LAPACK_COL_MAJOR, n, (MKL_Complex16*)a, lda, ipiv.data() );
#else
      int const lwork = n*n;
      std::vector<std::complex<double>> work(lwork);
      zgetri_(&n, a, &lda, ipiv.data(), work.data(), &lwork, &info); // Fortran interface
#endif
      return info;
  } // inverse

  inline status_t inverse(int const n, std::complex<float> a[], int const lda) {
      int info{0};
      std::vector<MKL_INT> ipiv(2*n);
#ifdef  HAS_MKL
      info = LAPACKE_cgetrf( LAPACK_COL_MAJOR, n, n, (MKL_Complex8*)a, lda, ipiv.data() );
#else
      cgetrf_(&n, &n, a, &lda, ipiv.data(), &info); // Fortran interface
#endif
      if (info) return info; // early return: factorization failed!

#ifdef  HAS_MKL
      info = LAPACKE_cgetri( LAPACK_COL_MAJOR, n, (MKL_Complex8*)a, lda, ipiv.data() );
#else
      int const lwork = n*n;
      std::vector<std::complex<float>> work(lwork);
      cgetri_(&n, a, &lda, ipiv.data(), work.data(), &lwork, &info); // Fortran interface
#endif
      return info;
  } // inverse
  
  
  inline status_t linear_solve(int const n, double a[], int const lda, double b[], int const ldb, int const nrhs=1) {
      std::vector<MKL_INT> ipiv(2*n);
#ifdef  HAS_MKL
      return LAPACKE_dgesv( LAPACK_COL_MAJOR, n, nrhs, a, lda, ipiv.data(), b, ldb );
#else
      int info{0};
      dgesv_(&n, &nrhs, a, &lda, ipiv.data(), b, &ldb, &info); // Fortran interface
      return info;
#endif
  } // linear_solve
  

  inline status_t _eigenvalues(int const n, double a[], int const lda, double w[]) {
#ifdef  HAS_MKL
      return LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', n, a, lda, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const lwork = (2*n + 2)*n;
      std::vector<double> work(lwork);
      dsyev_(&jobz, &uplo, &n, a, &lda, w, work.data(), &lwork, &info); // Fortran interface
      return info;
#endif
  } // (standard_)eigenvalues

  inline status_t _eigenvalues(int const n, std::complex<double> a[], int const lda, double w[]) {
#ifdef  HAS_MKL
      return LAPACKE_zheev( LAPACK_COL_MAJOR, 'V', 'U', n, (MKL_Complex16*)a, lda, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const lwork = (2*n + 2)*n;
      std::vector<std::complex<double>> work(lwork);
      std::vector<double> rwork(3*n);
      zheev_(&jobz, &uplo, &n, a, &lda, w, work.data(), &lwork, rwork.data(), &info); // Fortran interface
      return info;
#endif
  } // (standard_)eigenvalues
  
  inline status_t _generalized_eigenvalues(int const n, double a[], int const lda, double b[], int const ldb, double w[]) {
//    printf("\n# call dsygv(1, 'v', 'u', %i, %p, %i, %p, %i, %p)\n\n",   n, a, lda, b, ldb, w);
#ifdef  HAS_MKL
      return LAPACKE_dsygv( LAPACK_COL_MAJOR, 1, 'V', 'U', n, a, lda, b, ldb, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const itype = 1, lwork = (2*n + 2)*n;
      std::vector<double> work(lwork);
      dsygv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work.data(), &lwork, &info); // Fortran interface
      return info;
#endif
  } // generalized_eigenvalues

  inline status_t _generalized_eigenvalues(int const n, std::complex<double> a[], int const lda, std::complex<double> b[], int const ldb, double w[]) {
//    printf("\n# call zhegv(1, 'v', 'u', %i, %p, %i, %p, %i, %p)\n\n",   n, a, lda, b, ldb, w);
#ifdef  HAS_MKL
      return LAPACKE_zhegv( LAPACK_COL_MAJOR, 1, 'V', 'U', n, (MKL_Complex16*)a, lda, (MKL_Complex16*)b, ldb, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const itype = 1, lwork = (2*n + 2)*n;
      std::vector<std::complex<double>> work(lwork);
      std::vector<double> rwork(3*n);
      zhegv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work.data(), &lwork, rwork.data(), &info); // Fortran interface
      return info;
#endif
  } // generalized_eigenvalues

  
  // float32 version: copy the above, replace double by float and change d --> s, z --> c

  
  inline status_t _eigenvalues(int const n, float a[], int const lda, float w[]) {
#ifdef  HAS_MKL
      return LAPACKE_ssyev( LAPACK_COL_MAJOR, 'V', 'U', n, a, lda, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const lwork = (2*n + 2)*n;
      std::vector<float> work(lwork);
      ssyev_(&jobz, &uplo, &n, a, &lda, w, work.data(), &lwork, &info); // Fortran interface
      return info;
#endif
  } // (standard_)eigenvalues

  inline status_t _eigenvalues(int const n, std::complex<float> a[], int const lda, float w[]) {
#ifdef  HAS_MKL
      return LAPACKE_cheev( LAPACK_COL_MAJOR, 'V', 'U', n, (MKL_Complex8*)a, lda, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const lwork = (2*n + 2)*n;
      std::vector<std::complex<float>> work(lwork);
      std::vector<float> rwork(3*n);
      cheev_(&jobz, &uplo, &n, a, &lda, w, work.data(), &lwork, rwork.data(), &info); // Fortran interface
      return info;
#endif
  } // (standard_)eigenvalues
  
  inline status_t _generalized_eigenvalues(int const n, float a[], int const lda, float b[], int const ldb, float w[]) {
//    printf("\n# call ssygv(1, 'v', 'u', %i, %p, %i, %p, %i, %p)\n\n",   n, a, lda, b, ldb, w);
#ifdef  HAS_MKL
      return LAPACKE_ssygv( LAPACK_COL_MAJOR, 1, 'V', 'U', n, a, lda, b, ldb, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const itype = 1, lwork = (2*n + 2)*n;
      std::vector<float> work(lwork);
      ssygv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work.data(), &lwork, &info); // Fortran interface
      return info;
#endif
  } // generalized_eigenvalues

  inline status_t _generalized_eigenvalues(int const n, std::complex<float> a[], int const lda, std::complex<float> b[], int const ldb, float w[]) {
//    printf("\n# call chegv(1, 'v', 'u', %i, %p, %i, %p, %i, %p)\n\n",   n, a, lda, b, ldb, w);
#ifdef  HAS_MKL
      return LAPACKE_chegv( LAPACK_COL_MAJOR, 1, 'V', 'U', n, (MKL_Complex8*)a, lda, (MKL_Complex8*)b, ldb, w );
#else
      int info{0}; char const jobz = 'V', uplo = 'U'; int const itype = 1, lwork = (2*n + 2)*n;
      std::vector<std::complex<float>> work(lwork);
      std::vector<float> rwork(3*n);
      chegv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work.data(), &lwork, rwork.data(), &info); // Fortran interface
      return info;
#endif
  } // generalized_eigenvalues

    // new interface: for hermitian/symmetric [generalized] eigenvalue problems
    template<typename real_t>
    inline status_t eigenvalues(real_t w[], int const n, real_t a[], int const lda,
                                                 real_t b[]=nullptr, int const ldb=0) {
        return b ? _generalized_eigenvalues(n, a, lda, b, (ldb < n)?lda:ldb, w) : _eigenvalues(n, a, lda, w); }

    template<typename real_t>
    inline status_t eigenvalues(real_t w[], int const n, std::complex<real_t> a[], int const lda,
                                                 std::complex<real_t> b[]=nullptr, int const ldb=0) {
        return b ? _generalized_eigenvalues(n, a, lda, b, (ldb < n)?lda:ldb, w) : _eigenvalues(n, a, lda, w); }

  // matrix-matrix multiplication: in the case 'n','n' --> this is c[n,m] += b[n,k]*a[k,m]
  inline status_t gemm(int const M, int const N, int const K, double c[], int const ldc
          , double const b[], int const ldb, double const a[], int const lda
          , double const alpha=1, double const beta=0, char const transa='n', char const transb='n') {
      dgemm_(&transa, &transb, &M, &N, &K, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); return 0; }

  inline status_t gemm(int const M, int const N, int const K, float c[], int const ldc
          , float const b[], int const ldb, float const a[], int const lda
          , float const alpha=1, float const beta=0, char const transa='n', char const transb='n') {
      sgemm_(&transa, &transb, &M, &N, &K, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); return 0; }

  inline status_t gemm(int const M, int const N, int const K, std::complex<double> c[], int const ldc
          , std::complex<double> const b[], int const ldb, std::complex<double> const a[], int const lda
          , std::complex<double> const alpha=1, std::complex<double> const beta=0, char const transa='n', char const transb='n') {
      zgemm_(&transa, &transb, &M, &N, &K, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); return 0; }

  inline status_t gemm(int const N, int const M, int const K, std::complex<float> c[], int const ldc
          , std::complex<float> const b[], int const ldb, std::complex<float> const a[], int const lda
          , std::complex<float> const alpha=1, std::complex<float> const beta=0, char const transa='n', char const transb='n') {
      if (0) { 
          printf("# cgemm(transa=\'%c\', transb=\'%c\', M=%d, N=%d, K=%d, alpha=(%g, %g), A[][%d], B[][%d], beta=(%g, %g), C[][%d])\n",
                          transa, transb, M, N, K, alpha.real(), alpha.imag(), lda, ldb, beta.real(), beta.imag(), ldc);
          fflush(stdout);
          if (('n' == (transa | 32)) && ('n' == (transb | 32))) {
              for(int n = 0; n < N; ++n) {
                  for(int m = 0; m < M; ++m) {
                      std::complex<float> cc(0);
                      for(int k = 0; k < K; ++k) {
                          cc += b[n*ldb + k] * a[k*lda + m];
                      } // k
                      c[n*ldc + m] = beta*c[n*ldc + m] + alpha*cc;
                  } // m
              } // n
              return 0;
          } // trans and transb 
      } // 1
      cgemm_(&transa, &transb, &M, &N, &K, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); return 0; }


  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }

} // namespace linear_algebra
