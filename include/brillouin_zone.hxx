#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, ::snprintf
#include <algorithm> // std::min, ::max
#include <cmath> // std::round, ::sqrt

#include "status.hxx" // status_t
#include "data_view.hxx" // view2D, view4D
#include "inline_math.hxx" // set, pow2, product, is_integer
#include "print_tools.hxx" // printf_vector
#include "control.hxx" // ::get

namespace brillouin_zone {

  int constexpr WEIGHT = 3; // the weight is stored in the 3rd component
                            // while components 0,1,2 carry the k-vector

  inline int get_kpoint_mesh( // returns the number of k-points
        view2D<double> & mesh // on exit shape (nkpoints, 4)
      , unsigned const nv[3] // grid along each dimension
      , bool const complex_phase_factors=true
      , int const echo=0 // log-level
  ) {
      unsigned n[3];
      double shift[3];
      for (int d = 0; d < 3; ++d) {
          n[d] = std::max(1u, nv[d]);
          shift[d] = n[d] - 1.;
          if (!complex_phase_factors) {
              n[d] = std::min(n[d], 2u); // only Gamma and X point lead to real phase factors
              shift[d] = 0;
          } // d
      } // phase factors must be real

      auto const nfull = size_t(n[2])*size_t(n[1])*size_t(n[0]);
      view4D<double> full(n[2], n[1], n[0], 4); // get temporary memory
      double const denom[] = {.5/n[0], .5/n[1], .5/n[2]};
      double xyzw[] = {0, 0, 0, 1};
      for (int iz = 0; iz < n[2]; ++iz) {  xyzw[2] = (2*iz - shift[2])*denom[2];
      for (int iy = 0; iy < n[1]; ++iy) {  xyzw[1] = (2*iy - shift[1])*denom[1];
      for (int ix = 0; ix < n[0]; ++ix) {  xyzw[0] = (2*ix - shift[0])*denom[0];
          set(full(iz, iy, ix), 4, xyzw);
          if (echo > 18) std::printf("# k-point mesh entry %9.6f %9.6f %9.6f weight= %g\n", xyzw[0],xyzw[1],xyzw[2], xyzw[3]);
      }}} // iz iy iz

//       int const nsymmetries = 1;
//       view3D<int8_t> symmetry(nsymmetries, 4, 4, 0); // use only the 3x3 sub matrices of the 4x4 memory chunks
//       for (int d = 0; d < 3; ++d) symmetry(0,d,d) = -1; // time reversal

      int8_t const rot_mat[24*3*3] = {
         1,  0,  0,  0,  1,  0,  0,  0,  1,
        -1,  0,  0,  0, -1,  0,  0,  0,  1,
        -1,  0,  0,  0,  1,  0,  0,  0, -1,
         1,  0,  0,  0, -1,  0,  0,  0, -1,
         0,  1,  0,  1,  0,  0,  0,  0, -1,
         0, -1,  0, -1,  0,  0,  0,  0, -1,
         0, -1,  0,  1,  0,  0,  0,  0,  1,
         0,  1,  0, -1,  0,  0,  0,  0,  1,
         0,  0,  1,  0, -1,  0,  1,  0,  0,
         0,  0, -1,  0, -1,  0, -1,  0,  0,
         0,  0, -1,  0,  1,  0,  1,  0,  0,
         0,  0,  1,  0,  1,  0, -1,  0,  0,
        -1,  0,  0,  0,  0,  1,  0,  1,  0,
        -1,  0,  0,  0,  0, -1,  0, -1,  0,
         1,  0,  0,  0,  0, -1,  0,  1,  0,
         1,  0,  0,  0,  0,  1,  0, -1,  0,
         0,  0,  1,  1,  0,  0,  0,  1,  0,
         0,  0, -1, -1,  0,  0,  0,  1,  0,
         0,  0, -1,  1,  0,  0,  0, -1,  0,
         0,  0,  1, -1,  0,  0,  0, -1,  0,
         0,  1,  0,  0,  0,  1,  1,  0,  0,
         0, -1,  0,  0,  0, -1,  1,  0,  0,
         0, -1,  0,  0,  0,  1, -1,  0,  0,
         0,  1,  0,  0,  0, -1, -1,  0,  0};

      int const n24 = control::get("hamiltonian.kmesh.symmetries", 1.);
      assert(n24 < 25);

      int const has_inv = 1 + (control::get("hamiltonian.kmesh.inversion", 1.) > 0);

      for (int iz = 0; iz < n[2]; ++iz) {
      for (int iy = 0; iy < n[1]; ++iy) {
      for (int ix = 0; ix < n[0]; ++ix) {
          // apply time reversal symmetry (k and -k produce the same density)
        double const *const vec = full(iz,iy,ix);
        for (int inv = 0; inv < has_inv; ++inv) { // inversion symmetry
          double const f_inv = 1. - inv*2.;
          for (int s = 0; s < n24; ++s) {
              double xyz[] = {0, 0, 0};
              for (int i = 0; i < 3; ++i) {
                  for (int j = 0; j < 3; ++j) {
                      xyz[i] += f_inv*rot_mat[s*3*3 + i*3 + j] * vec[j];
                  } // j
              } // i
          // convert xyz to nearest integers
          int jxyz[3];
          for (int d = 0; d < 3; ++d) {
              jxyz[d] = int(std::round(n[d]*xyz[d] + 0.5*shift[d]));
              if (jxyz[d] < 0) jxyz[d] += n[d];
              assert(jxyz[d] < n[d]);
          } // d
          int const jx = jxyz[0], jy = jxyz[1], jz = jxyz[2];
          if (echo > 16) {
              std::printf("# symmetry #%i maps k-point", 0);
              printf_vector(" %9.6f", full(iz,iy,ix), 3, " to");
              printf_vector(" %9.6f", full(jz,jy,jx), 3);
          } // echo

          auto const w8 = full(jz,jy,jx,WEIGHT);
          if (w8 > 0) {
              double diff2{0};
              for (int d = 0; d < 3; ++d) {
                  diff2 += pow2(xyz[d] - full(jz,jy,jx,d));
              } // d
              if (diff2 < 1e-12) {
                  // k + (-k) == nullvector, transfer weight to (iz,iy,ix)
                  full(iz,iy,ix,WEIGHT) += w8;
                  full(jz,jy,jx,WEIGHT) -= w8;
                  if (echo > 14) {
                      std::printf("# transfer k-point weight %g from", w8);
                      printf_vector(" %9.6f", full(jz,jy,jx), 3, " to");
                      printf_vector(" %9.6f", full(iz,iy,ix), 3);
                  } // echo
              } // d2 is zero
          } // w8 > 0

          } // s
        } // inv
      }}} // iz iy iz

      // copy kpoints with positive weight into result list
      // up to here, all weights were integer, so we have to divide by the denominator nfull
      int constexpr COUNT=0, WRITE=1;
      double const wfull = 1./nfull;
      double const weighted[4] = {1, 1, 1, wfull};
      int nmesh{0};
      for (int i01 = COUNT; i01 <= WRITE; i01 += (WRITE - COUNT)) {
          double w8sum{0};
          int imesh{0};
          for (int iz = 0; iz < n[2]; ++iz) {
          for (int iy = 0; iy < n[1]; ++iy) {
          for (int ix = 0; ix < n[0]; ++ix) {
              auto const w8 = full(iz,iy,ix,WEIGHT);
              if (w8 > 0) {
                  if (WRITE == i01) product(mesh[imesh], 4, full(iz,iy,ix), weighted); // copy and apply weight denominator
                  w8sum += w8;
                  ++imesh;
              } // w8 > 0
          }}} // iz iy iz
          if (COUNT == i01) {
              nmesh = imesh;
              if (echo > 8) std::printf("# %d k-points have positive weight, weight sum = 1 + %.1e\n", nmesh, w8sum*wfull - 1);
              mesh = view2D<double>(nmesh, 4); // get memory
          } else {
              assert(nmesh == imesh && "1st and 2nd time counting did not agree");
          } // COUNT or WRITE
      } // twice

      auto const fraction = nmesh*wfull; // == nmesh/(n[X]*n[Y]*n[Z])
      if (echo > 3) std::printf("# k-point mesh with %d x %d x %d has %d points, %g %%, %g of %d\n",
                                 n[0],n[1],n[2], nmesh, fraction*100, 1/fraction, has_inv*24);
      return nmesh;
  } // get_kpoint_mesh


  inline int get_kpoint_mesh(
        view2D<double> & mesh
      , unsigned const n // isotropic
      , bool const complex_phase_factors=true
      , int const echo=0 // log-level
  ) {
      unsigned const nv[] = {n, n, n};
      return get_kpoint_mesh(mesh, nv, complex_phase_factors, echo);
  } // get_kpoint_mesh


  inline int get_kpoint_mesh(
        view2D<double> & mesh
  ) {
      unsigned nv[3];
      auto const iso = control::get("hamiltonian.kmesh", 1.); // isotropic default value
      nv[0]          = control::get("hamiltonian.kmesh.x", iso);
      nv[1]          = control::get("hamiltonian.kmesh.y", iso);
      nv[2]          = control::get("hamiltonian.kmesh.z", iso);
      int const echo = control::get("hamiltonian.kmesh.echo", 0.);
      int const cmpl = control::get("hamiltonian.kmesh.complex", 1.);
      return get_kpoint_mesh(mesh, nv, (0 != cmpl), echo);
  } // get_kpoint_mesh


  inline bool needs_complex(double const kvec[3]) {
      return !(is_integer(2*kvec[0])
            && is_integer(2*kvec[1])
            && is_integer(2*kvec[2]));
  } // needs_complex


  inline bool needs_complex(
        view2D<double> const & mesh
      , int const nmesh
  ) {
      for (int ik = 0; ik < nmesh; ++ik) {
          if (needs_complex(mesh[ik])) return true;
      } // ik
      return false;
  } // needs_complex


  inline char get_special_kpoint(double v[3], char const *const keyword, char const *const def="?") {
      v[0] = 0; v[1] = 0; v[2] = 0; // clear
      char const cin = *control::get(keyword, def);
      switch (cin | 32) { // convert to lower case
        case 'g': case '0':                           return '0';
        case 'x': case '1':               v[0] = 0.5; return 'X';
        case 'm':             v[1] = 0.5; v[0] = 0.5; return 'M';
        case 'r': v[2] = 0.5; v[1] = 0.5; v[0] = 0.5; return 'R';
        default : control::get(v, keyword);           return 'c';
            //    control::get(double v[], keyword) --> {keyword.x, keyword.y, keyword.z}
      } // switch
  } // get_special_kpoint


  inline int get_kpoint_path( // returns the number of k-points
        view2D<double> & path // on exit shape (nk + 1, 4)
      , int const echo=0 // log-level
  ) {
      double from[3], to[3];
      auto const c_from = get_special_kpoint(from, "hamiltonian.kpath.from", "Gamma");
      auto const c_to   = get_special_kpoint(to,   "hamiltonian.kpath.to",   "X");
      auto const dk = control::get("hamiltonian.kpath.spacing", 0.05);
      if (echo > 2) std::printf("# +hamiltonian.kpath.from=%c .to=%c .spacing=%g\n", c_from, c_to, dk);
      auto const length = std::sqrt(pow2(to[0] - from[0]) + pow2(to[1] - from[1]) + pow2(to[2] - from[2]));
      int const nk = std::ceil(length/std::max(dk, 1e-9));
      int const npath = nk + 1;
      path = view2D<double>(npath, 4, 0.0); // get memory
      double const by_nk = (nk > 0) ? 1./nk : 0;
      for (int ik = 0; ik <= nk; ++ik) {
          for (int d = 0; d < 3; ++d) {
              path(ik,d) = from[d] + (to[d] - from[d])*(ik*by_nk);
          } // d
          path(ik,3) = length*(ik*by_nk); // instead of a weight, this is a progress on the band path
      } // ik
      if (echo > 0) std::printf("# k-point path from %g %g %g (%c) to %g %g %g (%c) in %d steps, dk= %g\n",
                          from[0], from[1], from[2], c_from, to[0], to[1], to[2], c_to, nk, length*by_nk);
      return npath;
  } // get_kpoint_path









#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_mesh(int const echo=0) {
      status_t stat(0);
      for (int n = 0; n <= 8; ++n) {
          view2D<double> mesh;
          auto const nm = get_kpoint_mesh(mesh, n, true, echo);
          if (echo > 5) std::printf("# %s(n= %d) returns %d points\n", __func__, n, nm);
          stat += (nm < 1); // error if less than 1 (the Gamma points) is in the mesh
      } // n
      return stat;
  } // test_mesh

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_mesh(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace brillouin_zone
