#pragma once

#include <cstdio> // printf, fflush, stdout

#include "status.hxx" // status_t

#include "linear_algebra.hxx" // ::gemm
#include "inline_math.hxx" // product
#include "complex_tools.hxx" // complex_name

#define DEBUG

namespace dense_operator {

  template<typename wave_function_t, typename DoubleComplex_t=wave_function_t>
  class dense_operator_t {
    //
    //  An operator offering the action of dense-stored 
    //  Hamiltonian and Overlap matrix onto single vectors
    //

    public:
      typedef wave_function_t complex_t;
      typedef DoubleComplex_t doublecomplex_t;

    private:
      complex_t const *Hmt, *Smt;
      double const *Cnd;
      int nB, nBa;
      
      inline status_t matrix_vector_multiplication(complex_t mvec[]
                     , complex_t const mat[], complex_t const vec[], int const echo=0) const {
          if (echo > 0) { printf("# %s<%s> gemm\n", __func__, complex_name<complex_t>()); fflush(stdout); }
          return linear_algebra::gemm(1, nB, nB, mvec, nB, vec, nB, mat, nBa); }

    public:

      dense_operator_t(int const nB, int const stride
                  , complex_t const *Hmt          // Hamiltonian
                  , complex_t const *Smt=nullptr  // Overlap matrix
                  , double    const *Cnd=nullptr) // diagonal preconditioner
        : Hmt{Hmt}, Smt{Smt}, Cnd{Cnd}, nB{nB}, nBa{stride}
      { assert( nB <= nBa ); } // constructor

      status_t Hamiltonian(complex_t Hpsi[], complex_t const psi[], int const echo=0) const {
          return matrix_vector_multiplication(Hpsi, Hmt, psi, echo); // multiply Hpsi = Hmt*psi
      } // Hamiltonian

      status_t Overlapping(complex_t Spsi[], complex_t const psi[], int const echo=0) const {
          return use_overlap() ? matrix_vector_multiplication(Spsi, Smt, psi, echo) : 0;
      } // Overlapping

      status_t Conditioner(complex_t Cpsi[], complex_t const psi[], int const echo=0) const {
          if (use_precond()) product(Cpsi, nB, Cnd, psi); 
          return 0;
      } // Pre-Conditioner

      double get_volume_element() const { return 1.0; }
      size_t get_degrees_of_freedom() const { return size_t(nB); }
      bool use_precond() const { return (nullptr != Cnd); }
      bool use_overlap() const { return (nullptr != Smt); }
  }; // class dense_operator_t


  inline status_t test_construct_and_destroy(int const echo=0) {
      double matrix[3][4];
      dense_operator_t<double> const op(3, 4, matrix[0]);
      return op.get_degrees_of_freedom() - 3;
  } // test_construct_and_destroy

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_construct_and_destroy(echo);
      return stat;
  } // all_tests

} // namespace dense_operator
