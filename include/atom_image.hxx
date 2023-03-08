#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdint> // int32_t, int8_t
#include <cstdio> // std::printf
#include <vector> // std::vector<T>

#include "sho_tools.hxx" // ::nSHO
#include "inline_math.hxx" // align<N_bits>
#include "sho_projection.hxx" // ::get_sho_prefactors

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

namespace atom_image {

  class atom_image_t {
    public:

      atom_image_t(void) {} // default constructor
      atom_image_t(double const x, double const y, double const z,
                   int32_t const atom_id=-1,
                   int const ix=-128, int const iy=-128, int const iz=-128, int const Zi=-128)
          : _atom_id(atom_id) {
          _pos[0] = x; _pos[1] = y; _pos[2] = z;
          _index[0] = ix; _index[1] = iy; _index[2] = iz; _index[3] = Zi;
      } // constructor

      double const * pos() const { return _pos; }
      double pos(int const d) const { assert(0 <= d); assert(d < 3); return _pos[d]; }
      int32_t atom_id() const { return _atom_id; }
      int8_t const * index() const { return _index; }
      int index(int const d) const { assert(0 <= d); assert(d < 4); return _index[d]; }

    private:
      double  _pos[3];      // real space positions of the atomic image
      int32_t _atom_id{-1}; // global atom identifyer
      int8_t  _index[4];    // flexible indices, e.g. for phases
  }; // class atom_image_t


  class sho_atom_t { // an atom with SHO-projectors
    public:

      sho_atom_t(void) : _sigma(1.), _numax(-1), _atom_id(-1) {} // default constructor
      sho_atom_t(
            double const sigma
          , int const numax
          , int32_t const atom_id
          , double const *pos=nullptr
          , int8_t const Zi=-128
      )
          : _sigma(sigma), _numax(numax), _atom_id(atom_id), _images(0)
      {
          assert(sigma > 0);
          _ncoeff = sho_tools::nSHO(_numax);
          _stride = align<2>(_ncoeff);
          _matrix64 = std::vector<double>(2*_ncoeff*_stride, 0.0);
          _matrix32 = std::vector<float> (2*_ncoeff*_stride, 0.f);
          if (nullptr != pos) {
              _images.resize(1); // only the one base image
              _images[0] = atom_image::atom_image_t(pos[0], pos[1], pos[2], atom_id, 0,0,0, Zi);
          } // pos
      } // constructor

      template <typename real_t>
      inline real_t const * get_matrix(int const h0s1=0) const; // provide no implementation for the general case

      status_t set_matrix(
            double const values[] // data layout values[ncoeff][stride_values], coefficients are assumed in order_zyx --> therefore, we must match ncoeff == _ncoeff
          , int const ncoeff
          , int const stride_values
          , int const h0s1=0
          , double const factor=1
      ) {
          assert(0 == h0s1 || 1 == h0s1); // 0:hamiltonian H, 1:overlap S (or I, contains charge deficit)
          assert(ncoeff <= stride_values);
          assert(_ncoeff <= _stride); // sanity check

          auto const rescale = sho_projection::get_sho_prefactors<double>(_numax, _sigma);

          for (int ij = 0; ij < _ncoeff*_stride; ++ij) {
              _matrix64[h0s1*_ncoeff*_stride + ij] = 0; // clear
              _matrix32[h0s1*_ncoeff*_stride + ij] = 0; // clear
          } // ij

          assert(ncoeff == _ncoeff && "for need order_nxyz for setting a matrix of different size!");
          int const nc = std::min(ncoeff, _ncoeff); // take the lower number of coefficients
          for (int i = 0; i < nc; ++i) {
              for (int j = 0; j < nc; ++j) {
                  int const ij = (h0s1*_ncoeff + i)*_stride + j;
                  int const ij_values = i*stride_values + j;
                  _matrix64[ij] = rescale[i] * values[ij_values] * factor * rescale[j];
                  _matrix32[ij] = _matrix64[ij]; // convert to float
              } // j
          } // i
          return (ncoeff != _ncoeff); // report mismatch
      } // set_matrix

      status_t set_image_positions(
            double const atom_position[3]
          , int const nimages=1
          , view2D<double> const *periodic_positions=nullptr
          , view2D<int8_t> const *indices=nullptr
      ) {
          if (nullptr == periodic_positions) {
              _images.resize(1);
              _images[0] = atom_image_t(atom_position[0], atom_position[1], atom_position[2], _atom_id, 0,0,0);
              return (nimages - 1); // return inconsistency if nullptr==periodic_positions && 1!=nimages
          } // nullptr
          _images.resize(nimages);
          int8_t const i000[] = {0,0,0};
          for (int ii = 0; ii < nimages; ++ii) {
              double p[3];
              for (int d = 0; d < 3; ++d) {
                  p[d] = atom_position[d] + (*periodic_positions)(ii,d);
              } // d
              int8_t const *iv = indices ? ((*indices)[ii]) : i000;
              _images[ii] = atom_image_t(p[0], p[1], p[2], _atom_id, iv[0], iv[1], iv[2]);
          } // ii
          return 0;
      } // set_image_positions

      int32_t atom_id() const { return _atom_id; }
      int     numax()   const { return _numax; }
      double  sigma()   const { return _sigma; }
      int     stride()  const { return _stride; }
      int     nimages() const { return _images.size(); }
      double const * pos(int const ii=0) const { assert(ii >= 0); assert(ii < _images.size()); return _images[ii].pos(); }
      int8_t const * idx(int const ii=0) const { assert(ii >= 0); assert(ii < _images.size()); return _images[ii].index(); }

    private:
      double  _sigma{1.};
      int32_t _numax{-1};
      int32_t _atom_id{-1};
      int32_t _ncoeff{0};
      int32_t _stride{0};
      std::vector<double> _matrix64; // data layout matrix[Hmt0_Ovl1][ncoeff][stride], SHO-coefficient layout is order_zyx
      std::vector<float>  _matrix32; // data layout matrix[Hmt0_Ovl1][ncoeff][stride], SHO-coefficient layout is order_zyx
      std::vector<atom_image_t> _images;
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

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t all_tests(int const echo=0) {
      if (echo > 2) std::printf("# %s sizeof(atom_image_t) = %ld Byte\n", __FILE__, sizeof(atom_image_t));
      return (3*8 + 1*4 + 4*1 != sizeof(atom_image_t)); // make sure that this struct has 32 Byte
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace atom_image
