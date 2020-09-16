#include <cstdio> // printf
#include <cmath> // std::cos, std::abs, std::sqrt
#include <cassert> // assert
#include <vector> // std::vector<T>
// #include <type_traits> // std::is_same
#include <complex> // std::complex<real_t>

#ifndef HAS_no_MKL
  #include "mkl_dfti.h" // Dfti* (discrete Fourier transform interface of the Intel(c) Math Kernel Library)
#endif

#ifdef HAS_FFTW
extern "C" {
  #include <fftw3.h> // fftw_plan, fftw_plan_dft_r2r_3d
}
#endif

namespace fourier_transform {
  
  template<typename real_t, bool forward=true>
  status_t fft(real_t out[] // (out) indexing out[(iz*ng[1] + iy)*ng[0] + ix]
             , real_t imag[] // imaginary part of in when backward, imaginary part of out when forward
             , real_t const in[] // (in) indexing in[(iz*ng[1] + iy)*ng[0] + ix]
             , int const ng[3] // grid numbers
             , int const echo=0) { // log level
#ifndef HAS_no_MKL
      MKL_LONG status;
      MKL_LONG const l[3] = {ng[2], ng[1], ng[0]};
      size_t const ngall = l[2]*l[1]*l[0];
      DFTI_DESCRIPTOR_HANDLE my_desc_handle;
      status = DftiCreateDescriptor(&my_desc_handle, (sizeof(real_t) > 4) ? DFTI_DOUBLE : DFTI_SINGLE, DFTI_COMPLEX, 3, l);
      status = DftiSetValue(my_desc_handle, DFTI_COMPLEX_STORAGE, DFTI_REAL_REAL);
      status = DftiSetValue(my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
      status = DftiCommitDescriptor(my_desc_handle);
      std::vector<real_t> imag2(ngall, 0.0);
      if (forward) { // forward
          status = DftiComputeForward (my_desc_handle, (void*)in, (void*)imag2.data(), (void*)out, (void*)imag); // perform the forward FFT
      } else {
          // in the backtransform we are not interested in the imaginary part of the output, so the content of imag2 is ignored
          status = DftiComputeBackward(my_desc_handle, (void*)in, (void*)imag, (void*)out, (void*)imag2.data()); // perform the forward FFT
      }
      DftiFreeDescriptor(&my_desc_handle); // cleanup, can be moved out
      if (status != 0 && echo > 0) printf("# MKL-FFT returns status=%li\n", status);
      return status;
#else // not defined HAS_no_MKL

#ifdef HAS_FFTW
//    if (std::is_same<real_t, double>::value)
      if (true) {
          size_t const ngall = size_t(ng[2]) * ng[1] * ng[0];
          std::vector<std::complex<double>> cvi(ngall), cvo(ngall); // complex arrays
          for(size_t i = 0; i < ngall; ++i) { // ToDo: OpenMP for, SIMD
              auto const im = forward ? 0 : imag[i];
              cvi[i] = std::complex<double>(in[i], im);
          } // i
          auto const plan = fftw_plan_dft_3d(ng[2], ng[1], ng[0], (fftw_complex*) cvi.data(), 
                                                                  (fftw_complex*) cvo.data(), 
                                      forward ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
          if (nullptr == plan) return __LINE__; // error
          fftw_execute(plan);
          fftw_destroy_plan(plan);
          for(size_t i = 0; i < ngall; ++i) { // ToDo: OpenMP for, SIMD
              out[i] = cvo[i].real();
              if (forward) imag[i] = cvo[i].imag(); // forward
          } // i
          return 0;
      } // real_t == double
#endif // defined HAS_FFTW

      return -1; // has no FFT library
#endif // defined HAS_no_MKL     
  } // fft

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=1) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  template<typename real_t>
  inline status_t test_fft(int const echo=6) {
      if (echo > 0) printf("\n# %s:\n", __func__);
      int const ng[3] = {29, 13, 9};
      int const ngall = ng[2]*ng[1]*ng[0];
      std::vector<real_t> rs(ngall);
      double const pw[3] = {3./ng[0], 2./ng[1], 1./ng[2]};
      if (echo > 1) printf("# %s: set up a single plane wave as [%g %g %g]\n", __func__, pw[0]*ng[0], pw[1]*ng[1], pw[2]*ng[2]);
      for(int z = 0; z < ng[2]; ++z) {
      for(int y = 0; y < ng[1]; ++y) {
      for(int x = 0; x < ng[0]; ++x) {
                  int const i = (z*ng[1] + y)*ng[0] + x;
                  rs[i] = std::cos(2*constants::pi*((pw[0]*x + pw[1]*y + pw[2]*z)));
      }}} // zyx
      std::vector<real_t> ft(2*ngall); auto const ft_imag = ft.data() + ngall; // two arrays ft[_real] and ft_imag are adjacent in memory
      auto const status_fft = fft(ft.data(), ft_imag, rs.data(), ng); // forward
      real_t maximum = 0; int at[4] = {-1,-1,-1,-1};
      for(int reim = 0; reim < 2; ++reim) {
      for(int z = 0; z < ng[2]; ++z) {
      for(int y = 0; y < ng[1]; ++y) {
      for(int x = 0; x < ng[0]; ++x) {
                  int const i = (z*ng[1] + y)*ng[0] + x;
                  auto const fta = std::abs(ft[reim*ngall + i]);
                  if (fta > maximum) { maximum = fta; at[0] = x; at[1] = y; at[2] = z; at[3] = reim; }
      }}}} // czyx
      if (echo > 5) printf("# %s: detected peak at index [%d %d %d] %s-part, value %g\n", 
                      __func__, at[0], at[1], at[2], (at[3])?"imag":"real", maximum);
      std::vector<real_t> rs_back(ngall);
      auto const status_inv = fft<real_t, false>(rs_back.data(), ft_imag, ft.data(), ng); // backward
      if (echo > 8) printf("\n# %s: back-transformed cos-wave values:\n", __func__);
      real_t const omega_inv = 1./ngall;
      double deva = 0, dev2 = 0;
      for(int z = 0; z < ng[2]; ++z) {
      for(int y = 0; y < ng[1]; ++y) {
      for(int x = 0; x < ng[0]; ++x) {
                  int const i = (z*ng[1] + y)*ng[0] + x;
                  auto const d = rs_back[i]*omega_inv - rs[i];
                  deva += std::abs(d); dev2 += d*d;
                  if (echo > 8) printf("%d %g %g %g\n", i, rs_back[i]*omega_inv, rs[i], d);
      }}} // zyx
      if (echo > 2) printf("# back-transformed cos-wave differs abs %.1e rms %.1e\n", deva/ngall, std::sqrt(dev2/ngall));
      if (echo > 1) printf("# %s: status = %i\n\n", __func__, int(status_fft) + int(status_inv));
      return int(status_fft) + int(status_inv);
  } // test_fft

  inline status_t all_tests(int const echo) {
    status_t status(0);
    status += test_fft<float>(echo);
    status += test_fft<double>(echo);
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS
  
} // namespace fourier_transform
