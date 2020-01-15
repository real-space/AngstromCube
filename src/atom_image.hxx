#pragma once

#include <cstdint> // int32_t
#include <cstdio> // printf
#include <vector> // std::vector<T>
#include "sho_tools.hxx" // ::nSHO, ::construct_index_table, ::order_zyx
#include "inline_tools.hxx" // align<N_bits>
#include "sho_projection.hxx" // ::sho_prefactor

typedef int status_t;

namespace atom_image {
  
  class atom_image_t {
    public:

      atom_image_t() {} // default constructor
      atom_image_t(double const x, double const y, double const z, 
                  int32_t const atom_id, int32_t const index=-1)
          : _atom_id(atom_id), _index(index) {
          _pos[0] = x; _pos[1] = y; _pos[2] = z;
      } // contructor

      double const * get_pos() const { return _pos; };
      double get_pos(int const d) const { assert(0 <= d); assert(d < 3); return _pos[d]; };
      int32_t index()   const { return _index; }
      int32_t atom_id() const { return _atom_id; }

    private:
      double  _pos[3];      // real space positions of the atomic image
      int32_t _atom_id{-1}; // global atom identifyer
      int32_t _index{-1};   // flexible index, e.g. for local atom counting
  }; // class atom_image_t

  
  class sho_atom_t { // an atom with SHO-projectors
    public:
      sho_atom_t() : _sigma(1.), _numax(-1), _atom_id(-1) {} // default constructor

      sho_atom_t(int const numax, double const sigma, int32_t const atom_id) 
        : _sigma(sigma), _numax(numax), _atom_id(atom_id) {
          _ncoeff = sho_tools::nSHO(_numax);
          _stride = align<2>(_ncoeff);
          _matrix64 = std::vector<double>(2*_ncoeff*_stride, 0.0);
          _matrix32 = std::vector<float> (2*_ncoeff*_stride, 0.0);
      } // constructor

      template <typename real_t> 
      inline real_t const * get_matrix(int const h0s1=0) const { assert(0); return nullptr; } // will fail
      
      status_t set_matrix(double const values[], int const ncoeff, int const stride, int const h0s1=0) {
          assert(0 == h0s1 || 1 == h0s1);
          
          std::vector<int> reorder(_ncoeff, 0);
          sho_tools::construct_index_table(reorder.data(), _numax, sho_tools::order_zyx);
          // reorder is necessary as values has order_Ezyx and _matrix has order_zyx, could be moved out

          std::vector<double> rescale(_ncoeff, 0.0);
          {   int ii = 0;
              for(int z = 0; z <= _numax; ++z) {
                  for(int y = 0; y <= _numax - z; ++y) {
                      for(int x = 0; x <= _numax - z - y; ++x) {
                          assert(reorder[ii] == sho_tools::Ezyx_index(x, y, z));
                          rescale[ii] = sho_projection::sho_prefactor(x, y, z, _sigma);
                          ++ii;
                      } // x
                  } // y
              } // z
          } // rescale because projector functions used in a fast SHO-transform are not normalized, coule be moved out

          for(int ij = 0; ij < _ncoeff*_stride; ++ij) {
              _matrix64[h0s1*_ncoeff*_stride + ij] = 0; // clear
              _matrix32[h0s1*_ncoeff*_stride + ij] = 0; // clear
          } // ij

          int const nc = std::min(ncoeff, _ncoeff); // take the lower number of coefficients
          for(int i = 0; i < nc; ++i) {
              for(int j = 0; j < nc; ++j) {
                  int const ij = (h0s1*_ncoeff + i)*_stride + j;
                  _matrix64[ij] = rescale[i] * values[reorder[i]*stride + reorder[j]] * rescale[j];
                  _matrix32[ij] = _matrix64[ij]; // convert to float
              } // j
          } // i
          return (ncoeff != _ncoeff); // report missmatch
      } // set_matrix

      int32_t atom_id() const { return _atom_id; }
      int     numax()   const { return _numax; }
      double  sigma()   const { return _sigma; }
      int     stride()  const { return _stride; }

    private:
      double  _sigma{1.};
      int32_t _numax{-1};
      int32_t _atom_id{-1};

      std::vector<double> _matrix64; // data layout matrix[Hmt0_Ovl1][ncoeff][stride], SHO-coefficient layout is order_zyx
      std::vector<float>  _matrix32; // data layout matrix[Hmt0_Ovl1][ncoeff][stride], SHO-coefficient layout is order_zyx
      int32_t _ncoeff{0};
      int32_t _stride{0};
  }; // class sho_atom_t
  
      template <> // specialization for real_t=double
      inline double const * sho_atom_t::get_matrix<double>(int const h0s1) const { 
          assert((0 == h0s1) || (1 == h0s1));
          return &_matrix64[h0s1*_ncoeff*_stride];
      } // get_matrix

      template <> // specialization for real_t=float
      inline float  const * sho_atom_t::get_matrix<float> (int const h0s1) const { 
          assert((0 == h0s1) || (1 == h0s1));
          return &_matrix32[h0s1*_ncoeff*_stride];
      } // get_matrix

  
  inline status_t all_tests() {
      return (3*8 + 2*4 != sizeof(atom_image_t)); // make sure that this struct has 32 Byte
  } // all_tests

} // namespace atom_image
