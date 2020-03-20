#pragma once

#include <cstdint> // int32_t
#include <cstdio> // printf
#include <vector> // std::vector<T>
#include "sho_tools.hxx" // ::nSHO, ::construct_index_table, ::order_zyx
#include "inline_tools.hxx" // align<N_bits>
#include "sho_projection.hxx" // ::sho_prefactor

#include "status.hxx" // status_t

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
      sho_atom_t(double const sigma, int const numax, int32_t const atom_id) 
          : _sigma(sigma), _numax(numax), _atom_id(atom_id)
      {
          assert(sigma > 0);
          _ncoeff = sho_tools::nSHO(_numax);
          _stride = align<2>(_ncoeff);
          _matrix64 = std::vector<double>(2*_ncoeff*_stride, 0.0);
          _matrix32 = std::vector<float> (2*_ncoeff*_stride, 0.0);
      } // constructor

      template <typename real_t> 
      inline real_t const * get_matrix(int const h0s1=0) const; // provide no implementation for the general case

      status_t set_matrix(double const values[], // data layout values[ncoeff][stride]
                          int const ncoeff, int const stride, int const h0s1=0) {
          assert(0 == h0s1 || 1 == h0s1); // 0:hamiltonian H, 1:overlap S (or I, contains charge deficit)
          assert(ncoeff <= stride);
          
          std::vector<int> reorder(_ncoeff, 0);
          sho_tools::construct_index_table(reorder.data(), _numax, sho_tools::order_zyx);
          // reorder is necessary as values has order_Ezyx and _matrix has order_zyx, could be moved out

          std::vector<double> rescale(_ncoeff, 0.0);
          {   int ii = 0;
              for(int z = 0; z <= _numax; ++z) {
                  for(int y = 0; y <= _numax - z; ++y) {
                      for(int x = 0; x <= _numax - z - y; ++x) {
                          assert(reorder[ii] == sho_tools::Ezyx_index(x, y, z));
                          rescale[ii] = sho_projection::sho_prefactor(x, y, z, _sigma); // ToDo: can be replaced by the factorized version
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
          return (ncoeff != _ncoeff); // report mismatch
      } // set_matrix

      status_t set_image_positions(double const atom_position[3], int const nimages=1, double const *periodic_positions=nullptr) {
          if (nullptr == periodic_positions) {
              _images.resize(1);
              _images[0] = atom_image_t(atom_position[0], atom_position[1], atom_position[2], _atom_id, 0,0,0);
              return (nimages - 1); // return inconsistency if nullptr==periodic_positions && 1!=nimages
          } // nullptr
          _images.resize(nimages);
          for(int ii = 0; ii < nimages; ++ii) {
              double p[3];
              for(int d = 0; d < 3; ++d) {
                  p[d] = atom_position[d] + periodic_positions[4*ii + d];
              } // d
              _images[ii] = atom_image_t(p[0], p[1], p[2], _atom_id, 0,0,0);
          } // ii
          return 0;
      } // set_image_positions

      int32_t atom_id() const { return _atom_id; }
      int     numax()   const { return _numax; }
      double  sigma()   const { return _sigma; }
      int     stride()  const { return _stride; }
      int     nimages() const { return _images.size(); }
      double const * pos(int const ii=0) const { assert(ii >= 0); assert(ii < _images.size()); return _images[ii].pos(); }

    private:
      double  _sigma{1.};
      int32_t _numax{-1};
      int32_t _atom_id{-1};

      std::vector<double> _matrix64; // data layout matrix[Hmt0_Ovl1][ncoeff][stride], SHO-coefficient layout is order_zyx
      std::vector<float>  _matrix32; // data layout matrix[Hmt0_Ovl1][ncoeff][stride], SHO-coefficient layout is order_zyx
      int32_t _ncoeff{0};
      int32_t _stride{0};
      
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
  inline status_t all_tests(int const echo) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  inline status_t all_tests(int const echo=0) {
      return (3*8 + 4 + 4*1 != sizeof(atom_image_t)); // make sure that this struct has 32 Byte
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace atom_image
