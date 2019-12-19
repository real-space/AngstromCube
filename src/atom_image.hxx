#pragma once

#include <cstdint> // int32_t
#include <cstdio> // printf
#include <vector> // std::vector<T>
#include "sho_tools.hxx" // ::nSHO
#include "inline_tools.hxx" // align<N_bits>

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

  
  class sho_atom_t {
    public:
      sho_atom_t() {} // default constructor
      sho_atom_t(int const numax, double const sigma, int32_t const atom_id) 
        : _sigma(sigma), _numax(numax), _atom_id(atom_id) {
          _ncoeff = sho_tools::nSHO(_numax);
          _stride = align<2>(_ncoeff);
          _matrix = std::vector<double>(2*_ncoeff*_stride, 1.0);
      } // constructor

      double const * const get_matrix(int const h0s1=0) const { 
          assert(0 == h0s1 || 1 == h0s1);
          return &_matrix[h0s1*_ncoeff*_stride];
      } // get_matrix

      int32_t atom_id() const { return _atom_id; }
      int     numax()   const { return _numax; }
      double  sigma()   const { return _sigma; }
      int     stride()  const { return _stride; }

    private:
      std::vector<double> _matrix; // data layout matrix[Ovl0_Hmt1][ncoeff][stride]
      double  _sigma{1.};
      int32_t _numax{-1};
      int32_t _atom_id{-1};
      int32_t _ncoeff{0};
      int32_t _stride{0};
  }; // class sho_atom_t
  
  
  inline status_t all_tests() {
      return (3*8 + 2*4 != sizeof(atom_image_t)); // make sure that this struct has 32 Byte
  } // all_tests

} // namespace atom_image
