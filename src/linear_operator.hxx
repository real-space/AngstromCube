#pragma once

#include "vector_layout.hxx"
// #define CRTP_printf(...)

#include <cstdio> // printf
#define CRTP_printf(...) printf(__VA_ARGS__)

template <typename real_t> char const * real_t_name() { return nullptr; }
template <> char const * real_t_name<double>() { return "double"; }
template <> char const * real_t_name<float> () { return "float"; }

template <class CRTP_t, typename real_t>
class LinOp {
  // LinOp is a base class for a linear operator A (of which we only need the action)
  // using the Curiously Recurring Template Pattern (CRTP)
  // methods within LinOp can use this template to access public members of derived classes
  // furthermore, all derived classes are forced to offer implementations for
  //     _apply
  // and needs a member layout_ which offers .axpby, and .inner products
  // that does not need to be virtual but can even carry the inline attribute

  public:

    inline void apply(real_t Ax[], real_t const x[]) const {
        CRTP_printf("# %s:%i %s<%s>\n", __FILE__, __LINE__, __func__, real_t_name<real_t>());
        return static_cast<CRTP_t const *>(this)->_apply(Ax, x);
    } // apply

    inline void axpby(real_t y[], real_t const x[], real_t const *a=nullptr, real_t const *b=nullptr) const {
        CRTP_printf("# %s:%i %s<%s>\n", __FILE__, __LINE__, __func__, real_t_name<real_t>());
        return static_cast<CRTP_t const *>(this)->layout_.axpby(y, x, a, b);
    } // axpby

    inline void inner(real_t d[], real_t const x[], real_t const *y=nullptr) const {
        CRTP_printf("# %s:%i %s<%s>\n", __FILE__, __LINE__, __func__, real_t_name<real_t>());
        return static_cast<CRTP_t const *>(this)->layout_.inner(d, x, y);
    } // inner

  private:
    // no private members

}; // class LinOp<CRTP_t,real_t>



namespace linear_operator {
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  // example class
  template <typename real_t>
  class DiagOp : public LinOp<DiagOp<real_t>,real_t> {
    // the CRTP design pattern makes sure that DiagOp behaves like a LinOp
    public:

      DiagOp(VecLayout<real_t> const & layout) 
        : layout_(layout) // clone
        , diagonal_(layout.ndof(), 3.0) { // values on the diagonal could in principle differ
          CRTP_printf("# %s:%i %s<%s> (constructor)\n", __FILE__, __LINE__, __func__, real_t_name<real_t>());
      } // constructor

      void _apply(real_t Dx[], real_t const x[]) const {
          CRTP_printf("# %s:%i %s<%s>\n", __FILE__, __LINE__, __func__, real_t_name<real_t>());
          auto const m = layout_.stride(); // stride can be larger than nrhs_ for alignment
          for(int i = 0; i < layout_.ndof(); ++i) {
              for(int j = 0; j < layout_.nrhs(); ++j) {
                  auto const ij = i*m + j;
                  Dx[ij] = diagonal_[i]*x[ij]; // diagonal operator: just multiply
              } // j
          } // i
      } // _apply (implementation)

      VecLayout<real_t> layout_;
    private:
      std::vector<real_t> diagonal_; // vector of diagonal elements
  }; // class DiagOp


  inline status_t all_tests(int const echo=9) {
    if (echo > 0) printf("# %s\n", __func__);

    typedef double real_t;
    std::vector<real_t> Hx(7*8, 0.0), x(7*8, 1.0), xHx(8, 0.0);
    VecLayout<real_t> layout(7, 8);
    DiagOp<real_t> hamiltonian(layout);
    hamiltonian.apply(Hx.data(), x.data());
    hamiltonian.inner(xHx.data(), x.data(), Hx.data());
    hamiltonian.axpby(x.data(), Hx.data(), xHx.data());

    for(uint i = 0; i < xHx.size(); ++i) { printf(" %g", xHx[i]); } printf("\n"); // show
    // printf("# hamiltonian has been applied %ld times\n", hamiltonian.get_apply_count());

    return 0;
  } // all_tests

// in addition to .apply() what else do we need to run a general tfqmr inverter?
// nrm2, dotp, axpy, xpay + routines for the extraction of the results

#endif

} // namespace linear_operator
