#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, ::fflush, stdout

#include "status.hxx" // status_t

#include "linear_algebra.hxx" // ::gemm
#include "inline_math.hxx" // product
#include "complex_tools.hxx" // complex_name
#include "grid_operators.hxx" // ::kpoint_t<>

// #define DEBUG

namespace dense_operator {

  template <typename wave_function_t> // choose from {float, double, complex<float>, complex<double>}
  class dense_operator_t {
    //
    //  An operator offering the action of dense-stored 
    //  Hamiltonian and Overlap matrix onto single vectors
    //

    public:
      typedef wave_function_t complex_t;
      typedef grid_operators::kpoint_t<complex_t> kpt_t; // dummy arguments

    private:
      complex_t const *Hmt, *Smt; // Hamiltonian matrix, Overlap matrix
      complex_t const *Cnd; // Preconditioning operator
      int nB, nBa;

      inline status_t matrix_vector_multiplication(complex_t mvec[]
                     , complex_t const mat[], complex_t const vec[], int const echo=0) const {
          if (echo > 19) std::printf("# %s<%s> gemm\n", __func__, complex_name<complex_t>());
          return linear_algebra::gemm(nB, 1, nB, mvec, nB, vec, nB, mat, nBa);
      } // matrix_vector_multiplication

    public:

      dense_operator_t(
            int const nB                 // dimension (number of Basis functions)
          , int const stride             // access stride (must be >= nB)
          , complex_t const *Hmt         // Hamiltonian
          , complex_t const *Smt=nullptr // Overlap matrix, optional
          , complex_t const *Cnd=nullptr // diagonal preconditioner, optional
      )
        : Hmt{Hmt}, Smt{Smt}, Cnd{Cnd}, nB{nB}, nBa{stride}
      {
          assert( nB <= nBa );
      } // constructor

      status_t Hamiltonian(complex_t Hpsi[], complex_t const psi[], kpt_t const & kp, int const echo=0) const {
          return matrix_vector_multiplication(Hpsi, Hmt, psi, echo); // multiply Hpsi = Hmt*psi
      } // Hamiltonian

      status_t Overlapping(complex_t Spsi[], complex_t const psi[], kpt_t const & kp, int const echo=0) const {
          return use_overlap() ? matrix_vector_multiplication(Spsi, Smt, psi, echo) : 0;
      } // Overlapping

      status_t Conditioner(complex_t Cpsi[], complex_t const psi[], kpt_t const & kp, int const echo=0) const {
          if (use_precond()) product(Cpsi, nB, Cnd, psi); // diagonal preconditioner
          return 0;
      } // Pre-Conditioner

      double get_volume_element() const { return 1.0; }
      size_t get_degrees_of_freedom() const { return size_t(nB); }
      bool use_precond() const { return (nullptr != Cnd); }
      bool use_overlap() const { return (nullptr != Smt); }
  }; // class dense_operator_t

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_construct_and_destroy(int const echo=0) {
      double const matrix[3][4] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}};
      dense_operator_t<double> const op(3, 4, matrix[0]);
      return op.get_degrees_of_freedom() - 3;
  } // test_construct_and_destroy

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_construct_and_destroy(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace dense_operator

#ifdef DEBUG
  #undef DEBUG
#endif
