#include <cstdio> // printf
#include <cmath> // cos, abs
#include <cassert> // assert
#include <algorithm> // fill

#include "fourier_poisson.hxx"

#ifndef HAS_no_MKL
  #include "mkl_dfti.h" // Intel(c) Math Kernel Library, discrete Fourier transform interface
#endif

#ifdef HAS_FFTW
extern "C" {
  #include "fftw.h" // fftw_plan, fftw_plan_dft_r2r_3d(int n0, int n1, int n2, double *in, double *out, 
                    //    fftw_r2r_kind kind0, fftw_r2r_kind kind1, fftw_r2r_kind kind2, unsigned flags)
                    //    FFTW_REDFT00
}
#endif

#include "constants.hxx" // pi

#define call(x) { auto const error_status = (x); if (error_status) return error_status; }

namespace fourier_poisson {

  status_t fourier_solve(double x[] // (out) solution of Laplace*x == b
                        ,double const b[] // right hand side b
                        ,int const ng[3] // grid numbers
                        ,double const bravais[3][4] // unit cell extend
                        ) {

//     fftw_r2r_kind const kind = FFTW_REDFT00;
//     unsigned const flags = 0;
    size_t const ng_all = ng[0] * ng[1] * ng[2];
    auto x_F = new double[ng_all];

//     fftw_plan const plan = fftw_plan_dft_r2r_3d(ng[0], ng[1], ng[2], b, x_F, kind, kind, kind, flags);

    // transform to Fourier space
//   x_F = fft( b );

//     B[0] = 0; // charge neutrality, clear the k=[0 0 0]-component

    int const nh[3] = {ng[0]/2, ng[1]/2, ng[3]/2};

    for(int j2 = 0; j2 < ng[2]; ++j2) {
      int const k2 = j2 - (j2 > nh[2])*ng[2];
      for(int j1 = 0; j1 < ng[1]; ++j1) {
        int const k1 = j1 - (j1 > nh[1])*ng[1];
        for(int j0 = 0; j0 < ng[0]; ++j0) {
          int const k0 = j0 - (j0 > nh[0])*ng[0];

          // todo: implement non-trivial Bravais matrix
          double const kk = k0*k0 + k1*k1 + k2*k2 + 1e-12; // add some safety
            
          x_F[(j2*ng[1] + j1)*ng[0] + j0] /= kk;

          } // j0
        } // j1
      } // j2
    
    // transform to real space
//     x = fft( x_F );

      return 0;
  } // fourier_solve

  
  status_t fourier_solve_MKL(double x[] // (out) solution of Laplace*x == b
                        , double const b[] // right hand side b
                        , int const ng[3] // grid numbers
                        , double const bravais[3][4] // unit cell extend
                        ) {


      return 0;
  } // fourier_solve_MKL

  template<typename real_t>
  status_t fft_MKL(real_t out[] // (out) indexing out[(iz*ng[1] + iy)*ng[0] + ix]
                 , real_t imag[] // imaginary part of in when backward, imaginary part of out when forward
                 , real_t const in[] // (in) indexing in[(iz*ng[1] + iy)*ng[0] + ix]
                 , int const ng[3] // grid numbers
                 , char const direction='f' // in the case of backward, the representations are swapped
                  ) {
      MKL_LONG status;
      MKL_LONG const l[3] = {ng[2], ng[1], ng[0]};
      size_t const ngall = l[2]*l[1]*l[0];
      DFTI_DESCRIPTOR_HANDLE my_desc_handle;
      status = DftiCreateDescriptor(&my_desc_handle, (sizeof(real_t) > 4) ? DFTI_DOUBLE : DFTI_SINGLE, DFTI_COMPLEX, 3, l);
      status = DftiSetValue(my_desc_handle, DFTI_COMPLEX_STORAGE, DFTI_REAL_REAL);
      status = DftiSetValue(my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
      status = DftiCommitDescriptor(my_desc_handle);
      auto imag2 = new real_t[ngall];
      if ('f' == direction) { // forward
          std::fill(imag2, imag2 + ngall, 0); // assume that array in is real, so its imaginary part (here imag2) is zero
          status = DftiComputeForward (my_desc_handle, (void*)in, (void*)imag2, (void*)out, (void*)imag); // perform the forward FFT
      } else {
          // in the backtransform we are not interested in the imaginary part of the output, so the content of imag2 is ignored
          status = DftiComputeBackward(my_desc_handle, (void*)in, (void*)imag, (void*)out, (void*)imag2); // perform the forward FFT
      }
      status = DftiFreeDescriptor(&my_desc_handle); // cleanup, can be moved out
      delete [] imag2;
      return status;
  } // fft_MKL
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template<typename real_t>
  int test_fft(int echo=9) {
      printf("\n# %s: \n", __func__);
      int const ng[3] = {29, 13, 9};
      int const ngall = ng[2]*ng[1]*ng[0];
      auto rs = new real_t[ngall];
      double const pw[3] = {3./ng[0], 2./ng[1], 1./ng[2]};
      if (echo > 1) printf("# set up a single plane wave as [%g %g %g]\n", pw[0]*ng[0], pw[1]*ng[1], pw[2]*ng[2]);
      for(int z = 0; z < ng[2]; ++z) {
      for(int y = 0; y < ng[1]; ++y) {
      for(int x = 0; x < ng[0]; ++x) {
                  int const i = (z*ng[1] + y)*ng[0] + x;
                  rs[i] = std::cos(2*constants::pi*((pw[0]*x + pw[1]*y + pw[2]*z)));
      }}} // zyx
      auto ft = new real_t[2*ngall]; auto ft_imag = ft + ngall; // two arrays ft[_real] and ft_imag are adjacent in memory
      auto const status_fft = fft_MKL(ft, ft_imag, rs, ng, 'f'); // forward
      real_t maximum = 0; int at[4] = {-1,-1,-1,-1};
      for(int reim = 0; reim < 2; ++reim) {
      for(int z = 0; z < ng[2]; ++z) {
      for(int y = 0; y < ng[1]; ++y) {
      for(int x = 0; x < ng[0]; ++x) {
                  int const i = (z*ng[1] + y)*ng[0] + x;
                  auto const fta = std::abs(ft[reim*ngall + i]);
                  if (fta > maximum) { maximum = fta; at[0] = x; at[1] = y; at[2] = z; at[3] = reim; }
      }}}} // czyx
      if (echo > 5) printf("# detected peak at index [%d %d %d] %s-part, value %g\n", 
                             at[0], at[1], at[2], (at[3])?"imag":"real", maximum);
      auto rs_back = new real_t[ngall];
      auto const status_inv = fft_MKL(rs_back, ft_imag, ft, ng, 'b'); // backward
      if (echo > 8) printf("\n# back-transformed cos-wave values:\n");
      real_t const omega_inv = 1./ngall;
      double deva = 0, dev2 = 0;
      for(int z = 0; z < ng[2]; ++z) {
      for(int y = 0; y < ng[1]; ++y) {  // if (echo > 8) printf("\n");
      for(int x = 0; x < ng[0]; ++x) {
                  int const i = (z*ng[1] + y)*ng[0] + x;
                  auto const d = rs_back[i]*omega_inv - rs[i];
                  deva += std::abs(d); dev2 += d*d;
                  if (echo > 8) printf("%d %g %g %g\n", i, rs_back[i]*omega_inv, rs[i], d);
      }}} // zyx
      if (echo > 2) printf("# back-transformed cos-wave differs avg %g or %g\n", deva/ngall, dev2/ngall);
      return status_fft + status_inv;
  } // test_fft

  status_t all_tests() {
    auto status = 0;
    status += test_fft<float>();
    status += test_fft<double>();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS
  
} // namespace fourier_poisson
