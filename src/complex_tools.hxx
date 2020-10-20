#pragma once

#include <cstdio> // printf
#include <type_traits> // std::true_type, std::false_type
#include <complex> // std::complex<real_t>

  template<typename T> struct is_complex_t                  : public std::false_type {};
  template<typename T> struct is_complex_t<std::complex<T>> : public std::true_type  {};
  template<typename T> constexpr bool is_complex() { return is_complex_t<T>::value; } // not needed when using C++14

  template<typename real_t> char const * complex_name(real_t const x=0);
  template<> inline char const * complex_name<double>(double const x) { return "double"; }
  template<> inline char const * complex_name<float> (float  const x) { return "float"; }
  template<> inline char const * complex_name<std::complex<double>>(std::complex<double> const x) { return "complex<double>"; }
  template<> inline char const * complex_name<std::complex<float>> (std::complex<float>  const x) { return "complex<float>"; }

  template<typename real_t> real_t conjugate(real_t const x);
  template<> inline double conjugate<double>(double const x) { return x; };
  template<> inline float  conjugate<float> (float  const x) { return x; };
  template<typename real_t>
  std::complex<real_t> inline conjugate(std::complex<real_t> const x) { return std::conj(x); }

  template<typename complex_t, typename real_t> inline
  complex_t to_complex_t(std::complex<real_t> const x); // no generic implementation given
  template<> inline double to_complex_t(std::complex<double> const x) { return std::real(x); }
  template<> inline float  to_complex_t(std::complex<float>  const x) { return std::real(x); }
  template<> inline std::complex<double> to_complex_t(std::complex<double> const x) { return x; }
  template<> inline std::complex<float>  to_complex_t(std::complex<float>  const x) { return x; }

  // use this function to extract the type of matrix elements e.g. 
  //    using double_complex_t = decltype(to_double_complex_t(complex_t(1)));
  inline std::complex<double> to_double_complex_t(std::complex<double> const x) { return x; }
  inline std::complex<double> to_double_complex_t(std::complex<float>  const x) { return std::complex<double>(x.real(), x.imag()); }
  inline double               to_double_complex_t(double               const x) { return x; }
  inline double               to_double_complex_t(float                const x) { return double(x); }

namespace complex_tools {

    inline char const * bool2string(bool const b) { return b ? "true" : "false"; }

    template<typename complex_t>
    inline status_t test_complex(bool const expect_complex, int const echo=0) {
        bool const is = is_complex<complex_t>();
        if (echo > 0) printf("# %s is_complex<%s>() = %s\n", __func__,
                                complex_name<complex_t>(), bool2string(is));
        if (echo > 1) printf("# typeof(conjugate(%s x)) = %s\n",
                                complex_name<complex_t>(),
                                complex_name<decltype(conjugate(complex_t(1)))>());
        return (expect_complex != is);
    } // test_complex

    inline status_t all_tests(int const echo=0) {
        status_t stat(0);
        stat += test_complex<std::complex<double>>(true, echo);
        stat += test_complex<std::complex<float>> (true, echo);
        stat += test_complex<double>(false, echo);
        stat += test_complex<float> (false, echo);
        return stat;
    } // all_tests

} // namespace complex_tools
