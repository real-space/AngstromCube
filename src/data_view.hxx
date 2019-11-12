#include <cstdio> // printf
#include <cassert> // assert

template<typename T>
class view2D {
public:
  
  view2D(T* const ptr, size_t const stride) : _data(ptr), _stride(stride), _memory_owner(false) { } // constructor
  view2D(size_t const n1, size_t const stride) : _stride(stride), _memory_owner(true) { _data = new T[n1*stride]; } // constructor

  ~view2D() { if(_memory_owner) delete[] _data; } // destructor

#ifdef _VIEW2D_HAS_PARENTHESIS
  T const & operator () (size_t const i1, size_t const i0) const { return _data[i1*_stride + i0]; } // (i,j)
  T       & operator () (size_t const i1, size_t const i0)       { return _data[i1*_stride + i0]; } // (i,j)

  T const & at(size_t const i1, size_t const i0) const { assert(i0 < _stride); return _data[i1*_stride + i0]; }
  T       & at(size_t const i1, size_t const i0)       { assert(i0 < _stride); return _data[i1*_stride + i0]; }
#endif

  T* operator [] (size_t const i1) const { return &_data[i1*_stride]; } // []

  T* data() const { return _data; }
  size_t stride() const { return _stride; }
  bool is_memory_owner() const { return _memory_owner; }

  // private data members
private:
  T* _data;
  size_t _stride;
  bool _memory_owner;

  // to prevent unwanted copying:
  view2D(view2D<T> const &); // delete the copy constructor
  view2D& operator = (view2D<T> const &); // delete the copy assignment constructor
  
}; // view2D

namespace data_view {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

inline int test_view2D(int const echo=9) {
    int constexpr n1 = 3, n0 = 5;
    view2D<double> a(3,8); // memory allocation
    assert(a.stride() >= n0);
    for(int i = 0; i < n1; ++i) {
        for(int j = 0; j < n0; ++j) {
            #ifdef _VIEW2D_HAS_PARENTHESIS
            a(i,j) = i + 0.1*j;
            if (echo > 0) printf("a(%i,%i) = %g\n", i, j, a(i,j));
            assert(a.at(i,j) == a[i][j]);
            #else
            a[i][j] = i + 0.1*j;
            if (echo > 0) printf("a(%i,%i) = %g\n", i, j, a[i][j]);
            #endif
        } // j
    } // i

    int ii = 1;
    if (echo > 0) printf("\n# ai = a[%i][:]\n", ii);
    auto const ai = a[ii]; // pointer into contiguous memory
    for(int j = 0; j < n0; ++j) {
        if (echo > 0) printf("ai[%i] = %g\n", j, ai[j]);
        #ifdef _VIEW2D_HAS_PARENTHESIS
        assert(a.at(ii,j) == ai[j]);
        #else
        assert(a[ii][j] == ai[j]);
        #endif
    } // j

} // test_view2D

inline int all_tests(int const echo=3) {
    int status = 0;
    status += test_view2D(echo);
} // all_tests

#endif // NO_UNIT_TESTS  

} // namespace data_view
