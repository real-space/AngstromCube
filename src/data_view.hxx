#include <cstdio> // printf
#include <cassert> // assert

#define DimUnknown 0

template<typename T>
class view2D {
public:
  
  view2D(T* const ptr, size_t const stride) 
    : _data(ptr), _n0(stride), _n1(DimUnknown) { } // constructor
  view2D(size_t const n1, size_t const stride, T const init_value=T(0)) 
    : _n0(stride), _n1(n1) { 
        auto const all = _n1*_n0;
        _data = new T[all];
        for(size_t i = 0; i < all; ++i) _data[i] = init_value; // warning! first touch here!
    } // constructor

  ~view2D() { if (is_memory_owner()) delete[] _data; } // destructor

#ifdef  _VIEW2D_HAS_PARENTHESIS
  T const & operator () (size_t const i1, size_t const i0) const { return _data[i1*_n0 + i0]; } // (i,j)
  T       & operator () (size_t const i1, size_t const i0)       { return _data[i1*_n0 + i0]; } // (i,j)

  T const & at(size_t const i1, size_t const i0) const { assert(i0 < _n0); return _data[i1*_n0 + i0]; }
  T       & at(size_t const i1, size_t const i0)       { assert(i0 < _n0); return _data[i1*_n0 + i0]; }
#endif

  T* operator [] (size_t const i1) const { return &_data[i1*_n0]; } // []

  T* data() const { return _data; }
  size_t stride() const { return _n0; }
  bool is_memory_owner() const { return (_n1 > DimUnknown); }

  // private data members
private:
  T* _data;
  size_t _n0, _n1; // _n1==0 --> unknown

  // to prevent unwanted copying:
  view2D(view2D<T> const &); // delete the copy constructor
  view2D& operator = (view2D<T> const &); // delete the copy assignment constructor
  
}; // view2D


template<typename T>
class view3D {
public:
  
  view3D(T* const ptr, size_t const n1, size_t const stride) 
    : _data(ptr), _n0(stride), _n1(n1), _n2(DimUnknown) { } // constructor
  view3D(size_t const n2, size_t const n1, size_t const stride, T const init_value=T(0)) 
    :  _n0(stride), _n1(n1), _n2(n2) {
        auto const all = _n2*_n1*_n0;
        _data = new T[all];
        for(size_t i = 0; i < all; ++i) _data[i] = init_value; // warning! first touch here!
    } // constructor

  ~view3D() { if (is_memory_owner()) delete[] _data; } // destructor

#define _VIEW3D_HAS_PARENTHESIS
#ifdef  _VIEW3D_HAS_PARENTHESIS
#define _access return _data[(i2*_n1 + i1)*_n0 + i0]
  T const & operator () (size_t const i2, size_t const i1, size_t const i0) const { _access; }
  T       & operator () (size_t const i2, size_t const i1, size_t const i0)       { _access; }

  T const & at(size_t const i2, size_t const i1, size_t const i0) const { assert(i1 < _n1); assert(i0 < _n0); _access; }
  T       & at(size_t const i2, size_t const i1, size_t const i0)       { assert(i1 < _n1); assert(i0 < _n0); _access; }
#undef _access
#endif

#ifdef  _VIEW3D_HAS_INDEXING
  T* operator [] (size_t const i2) const { return &_data[i2*_n1*_n0]; } // [] returns the inner 2 dimensions as flat array
#endif

  T* data() const { return _data; }
  size_t stride() const { return _n0; }
  bool is_memory_owner() const { return (_n2 > DimUnknown); }

  // private data members
private:
  T* _data;
  size_t _n0, _n1, _n2; // _n2==0 -->unknown

  // to prevent unwanted copying:
  view3D(view3D<T> const &); // delete the copy constructor
  view3D& operator = (view3D<T> const &); // delete the copy assignment constructor

}; // view3D

#undef DimUnknown

namespace data_view {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

inline int test_view2D(int const echo=9) {
    int constexpr n1 = 3, n0 = 5;
    if (echo > 0) printf("\n# %s(%i,%i)\n", __func__, n1, n0);
    view2D<double> a(n1,8); // memory allocation
    assert(a.stride() >= n0);
    for(int i = 0; i < n1; ++i) {
        for(int j = 0; j < n0; ++j) {
            #ifdef _VIEW2D_HAS_PARENTHESIS
            a(i,j) = i + 0.1*j;
            if (echo > 0) printf("# a2D(%i,%i) = %g\n", i, j, a(i,j));
            assert(a.at(i,j) == a[i][j]);
            #else
            a[i][j] = i + 0.1*j;
            if (echo > 0) printf("# a2D(%i,%i) = %g\n", i, j, a[i][j]);
            #endif
        } // j
    } // i

    int ii = 1;
    if (echo > 0) printf("\n# ai = a1D[%i][:]\n", ii);
    auto const ai = a[ii]; // pointer into contiguous memory
    for(int j = 0; j < n0; ++j) {
        if (echo > 0) printf("# ai[%i] = %g\n", j, ai[j]);
        #ifdef _VIEW2D_HAS_PARENTHESIS
        assert(a.at(ii,j) == ai[j]);
        #else
        assert(a[ii][j] == ai[j]);
        #endif
    } // j

    return 0;
} // test_view2D

inline int test_view3D(int const echo=9) {
    int constexpr n2 = 3, n1 = 2, n0 = 5;
    if (echo > 0) printf("\n# %s(%i,%i,%i)\n", __func__, n2, n1, n0);
    view3D<double> a(n2,n1,8); // memory allocation
    assert(a.stride() >= n0);
    for(int h = 0; h < n2; ++h) {
      for(int i = 0; i < n1; ++i) {
        for(int j = 0; j < n0; ++j) {
            a(h,i,j) = h + 0.1*i + 0.01*j;
            if (echo > 0) printf("# a3D(%i,%i,%i) = %g\n", h, i, j, a(h,i,j));
            #ifdef _VIEW3D_HAS_INDEXING
            assert(a(h,i,j) == a[h][i*a.stride() + j]);
            #endif
        } // j
      } // h
    } // i

    return 0;
} // test_view3D


inline int all_tests(int const echo=3) {
    int status = 0;
    status += test_view2D(echo);
    status += test_view3D(echo);
    return status;
} // all_tests

#endif // NO_UNIT_TESTS  

} // namespace data_view
