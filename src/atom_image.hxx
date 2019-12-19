#pragma once

#include <cstdint> // int32_t
#include <cstdio> // printf
#include <vector> // std::vector<T>

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
      double get_pos(int const d) const { assert(d < 3); assert(d >= 0); return _pos[d]; };
      int32_t index()   const { return _index; }
      int32_t atom_id() const { return _atom_id; }

    private:
      double  _pos[3];
      int32_t _atom_id{-1};
      int32_t _index{-1};
  }; // class atom_image_t

  
  class sho_atom_t {
    public:
      double* matrix{nullptr}; // data layout matrix[Ovl0_Hmt1][nSHO(numax)][stride]
      double  sigma{1.};
      int     numax{-1};
      int     stride{0};
      int     atom_id{-1};
  }; // class sho_atom_t
  
  
  inline status_t all_tests() {
      return (3*8 + 2*4 != sizeof(atom_image_t)); // make sure that this struct has 32 Byte
  } // all_tests

} // namespace atom_image
