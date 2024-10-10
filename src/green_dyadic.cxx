// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cmath> // std::exp
#include <vector> // std::vector<T>

#include "green_dyadic.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#ifndef   NO_UNIT_TESTS
    #include "dyadic_plan.hxx" // dyadic_plan_t
    #include "green_memory.hxx" // get_memory, free_memory, dim3, real_t_name
    #include "green_sparse.hxx" // ::sparse_t<>
    #include "inline_math.hxx" // pow2, pow3
    #include "sho_tools.hxx" // ::nSHO, ::n2HO, ::n1HO
    #include "constants.hxx" // ::sqrtpi
    #include "green_parallel.hxx" // ::rank, ::size, ::dyadic_exchange
    #include "control.hxx" // ::get
#endif // NO_UNIT_TESTS

#ifndef   HAS_NO_CUDA
    #include <cuda/std/complex> // std::complex
    #define std__complex cuda::std::complex
#else  // HAS_NO_CUDA
    #include <complex> // std::complex
    #define std__complex std::complex
#endif // HAS_NO_CUDA

namespace green_dyadic {


#ifdef    NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS


  status_t test_Hermite_polynomials_1D(int const echo=0, double const sigma=16) {
      int constexpr Lmax = 7;
      double H1D[Lmax + 1][3][4]; // 1D-Hermite polynomials times Gaussian decay function for 3 directions and grid points
      float  h1D[Lmax + 1][3][4]; // H1D in single precision
      float     xi_squared[3][4]; // (x/sigma)^2
      float  const xyzc[] = {0, 0, 0}; // cube position vector [in units of 4 grid points]
      double const hxyz[] = {1, 1, 1,   6.2832}; // grid spacing vector [in Bohr] and truncation_radius/sigma
      double       xyza[] = {0, 0, 0,   1./std::sqrt(sigma)}; // atom position vector [in Bohr] and spread^{-1/2}
      if (echo > 10) std::printf("\n## %s xi, H1D[0:%d]:", __func__, Lmax);
      double maxdev{0};
      for (int it = 0; it < 9; ++it) {
          for (int i3 = 0; i3 < 3; ++i3) {
              xyza[i3] = -4*(it*3 + i3);
              for (int i4 = 0; i4 < 4; ++i4) {
                  int const i = i3*4 + i4; // == threadIdx.x
                  Hermite_polynomials_1D(h1D, xi_squared, i, Lmax, xyza, xyzc, hxyz); // float
                  Hermite_polynomials_1D(H1D, xi_squared, i, Lmax, xyza, xyzc, hxyz); // double
//                if (echo > 10) std::printf("# %s idir=%d i4=%d ivec=%i xi^2=%g H1D= ", __func__, i3, i4, i, xi_squared[i3][i4]);
                  if (echo > 10) std::printf("\n%g  ", std::sqrt(xi_squared[i3][i4]));
                  for (int l = 0; l <= Lmax; ++l) {
                      if (echo > 10) std::printf(" %g", H1D[l][i3][i4]);
                      auto const dev = h1D[l][i3][i4] - H1D[l][i3][i4];
                      if (0 != H1D[l][i3][i4]) {
                          auto const reldev = std::abs(dev/H1D[l][i3][i4]);
                          maxdev = std::max(maxdev, reldev);
                      }
                  } // l
              } // i4
          } // i3
      } // it translations
      if (echo > 3) std::printf("\n# %s largest relative deviation between float and double is %.1e\n", __func__, maxdev);
      if (echo > 10) std::printf("\n\n\n");
      return 0;
  } // test_Hermite_polynomials_1D


    template <typename real_t, int R1C2=2, int Noco=1>
    size_t multiply( // deprecated, without atomic images
          real_t         (*const __restrict__ Ppsi)[R1C2][Noco*64][Noco*64] // result,  modified Green function blocks [nnzb][R1C2][Noco*64][Noco*64]
        , real_t         (*const __restrict__  Cpr)[R1C2][Noco]   [Noco*64] // projection coefficients     [natomcoeffs*nrhs][R1C2][Noco   ][Noco*64]
        , real_t   const (*const __restrict__  psi)[R1C2][Noco*64][Noco*64] // input, unmodified Green function blocks [nnzb][R1C2][Noco*64][Noco*64]
        , double   const (*const __restrict__ AtomPos)[3+1] // atomic positions and spread parameter
        , int8_t   const (*const __restrict__ AtomLmax) // SHO basis size [nAtomImages]
        , uint32_t const (*const __restrict__ AtomStarts) // prefetch sum over nSHO(AtomLmax[:])
        , uint32_t const natoms
        , green_sparse::sparse_t<> const (*const __restrict__ sparse_SHOprj)
        , double   const (*const          *const __restrict__ AtomMatrices)
        , green_sparse::sparse_t<> const &                    sparse_SHOadd
        , uint32_t const (*const __restrict__ RowIndexCubes) // Green functions rowIndex[nnzb]
        , uint16_t const (*const __restrict__ ColIndexCubes) // Green functions colIndex[nnzb]
        , float    const (*const __restrict__ CubePos)[3+1] // cube positions +alignment
        , double   const (*const __restrict__ hGrid) // grid spacings + projection radius
        , uint32_t const nnzb // number of blocks of the Green function
        , uint32_t const nrhs // number of block columns of the Green function
        , int const echo=0 // log-level
    ) {
        assert((1 == Noco && (1 == R1C2 || 2 == R1C2)) || (2 == Noco && 2 == R1C2));
        if (natoms*nrhs < 1) return 0; // zero

        auto const natomcoeffs = AtomStarts[natoms];
        if (echo > 6) std::printf("# %s<%s> R1C2=%d Noco=%d natoms=%d nrhs=%d ncoeffs=%d\n", __func__, real_t_name<real_t>(), R1C2, Noco, natoms, nrhs, natomcoeffs);

        SHOprj_driver<real_t,R1C2,Noco>(Cpr, psi, AtomPos, AtomLmax, AtomStarts, natoms, sparse_SHOprj, RowIndexCubes, CubePos, hGrid, nrhs, echo);

        // ToDo: in the future, nAtomImages and natoms can be different. Then, we need an extra step of accumulating the projection coefficients
        //       with the corresponding phases and, after SHOmul, distributing the addition coefficients to the images with the inverse phases.

        auto Cad = get_memory<real_t[R1C2][Noco][Noco*64]>(natomcoeffs*nrhs, echo, "Cad"); // rectangular storage of atoms x right-hand-sides

        // very rigid version, R1C2 and Noco could be function arguments instead of template parameters
        SHOmul_driver<real_t,R1C2,Noco,64>(Cad, Cpr, AtomMatrices, AtomLmax, AtomStarts, natoms, nrhs, echo);

        SHOadd_driver<real_t,R1C2,Noco>(Ppsi, Cad, AtomPos, AtomLmax, AtomStarts, sparse_SHOadd.rowStart(), sparse_SHOadd.colIndex(), RowIndexCubes, ColIndexCubes, CubePos, hGrid, nnzb, nrhs, echo);

        free_memory(Cad);

        return 0; // total number of floating point operations not computed
    } // multiply (deprecated)



    inline std::vector<double> __host__ sho_normalization(int const lmax, double const sigma=1) {
        int const n1ho = sho_tools::n1HO(lmax);
        std::vector<double> v1(n1ho, constants::sqrtpi*sigma);
        {
            double fac{1};
            for (int nu = 0; nu <= lmax; ++nu) {
                // fac == factorial(nu) / 2^nu
                v1[nu] *= fac;
                // now v1[nu] == sqrt(pi) * sigma * factorial(nu) / 2^nu
                fac *= 0.5*(nu + 1); // update fac for the next iteration
            } // nu
        }

        std::vector<double> vec(sho_tools::nSHO(lmax), 1.);
        {
            int sho{0};
            for (int iz = 0; iz <= lmax; ++iz) {
                for (int iy = 0; iy <= lmax - iz; ++iy) {
                    for (int ix = 0; ix <= lmax - iz - iy; ++ix) {
                        vec[sho] = v1[ix] * v1[iy] * v1[iz];
                        ++sho;
            }}} // ix iy iz
            assert(vec.size() == sho);
        }
        return vec;
    } // sho_normalization


  template <typename real_t, int R1C2=2, int Noco=1>
  status_t test_SHOprj_and_SHOadd(int const echo=0, int8_t const lmax=6) {
      // check if drivers compile and the normalization of the lowest (up to 64) SHO functions
      auto const sigma = control::get("green_dyadic.test.sigma", 1.);
      auto const hg    = control::get("green_dyadic.test.grid.spacing", 0.25);
      auto const rc    = control::get("green_dyadic.test.rc", 7.);
      int  const nb    = control::get("green_dyadic.test.nb", 14.);
      int  const natoms = 1, nrhs = 1, nnzb = pow3(nb);
      int  const nsho = sho_tools::nSHO(lmax);
      auto  psi = get_memory<real_t[R1C2][Noco*64][Noco*64]>(nnzb, echo, "psi");
      auto Vpsi = get_memory<real_t[R1C2][Noco*64][Noco*64]>(nnzb, echo, "Vpsi");
      set(psi[0][0][0], nnzb*R1C2*pow2(Noco*64ull), real_t(0)); // clear
      auto apc = get_memory<real_t[R1C2][Noco   ][Noco*64]>(natoms*nsho*nrhs, echo, "apc");
      set(apc[0][0][0], natoms*nsho*nrhs*R1C2*pow2(Noco)*64, real_t(0)); // clear

      auto sparse_SHOprj = get_memory<green_sparse::sparse_t<>>(nrhs, echo, "sparse_SHOprj");
      {
          std::vector<uint32_t> iota(nnzb); for (int inzb{0}; inzb < nnzb; ++inzb) { iota[inzb] = inzb; }
          std::vector<std::vector<uint32_t>> SHO_prj(natoms, iota);
          assert(nrhs > 0);
          sparse_SHOprj[0] = green_sparse::sparse_t<>(SHO_prj, false, __func__, echo - 9);
      }
      green_sparse::sparse_t<> sparse_SHOadd;
      {
          std::vector<std::vector<uint32_t>> SHO_add(nnzb, std::vector<uint32_t>(1, 0));
          sparse_SHOadd    = green_sparse::sparse_t<>(SHO_add, false, __func__, echo - 9);
      }

      auto ColIndexCubes = get_memory<uint16_t>(nnzb, echo, "ColIndexCubes");     set(ColIndexCubes, nnzb, uint16_t(0));
      auto RowIndexCubes = get_memory<uint32_t>(nnzb, echo, "RowIndexCubes");     for (int inzb = 0; inzb < nnzb; ++inzb) RowIndexCubes[inzb] = inzb;
      auto hGrid         = get_memory<double>(3+1, echo, "hGrid");                set(hGrid, 3, hg); hGrid[3] = rc;
      auto AtomPos       = get_memory<double[3+1]>(natoms, echo, "AtomPos");      set(AtomPos[0], 3, hGrid, 0.5*4*nb);  AtomPos[0][3] = 1./std::sqrt(sigma);
      auto AtomLmax      = get_memory<int8_t>(natoms, echo, "AtomLmax");          set(AtomLmax, natoms, lmax);
      auto AtomStarts    = get_memory<uint32_t>(natoms + 1, echo, "AtomStarts");  for(int ia = 0; ia <= natoms; ++ia) AtomStarts[ia] = ia*nsho;
      auto CubePos       = get_memory<float[3+1]>(nnzb, echo, "CubePos");
      for (int iz = 0; iz < nb; ++iz) {
      for (int iy = 0; iy < nb; ++iy) {
      for (int ix = 0; ix < nb; ++ix) {
          int const xyz0[] = {ix, iy, iz, 0};
          set(CubePos[(iz*nb + iy)*nb + ix], 4, xyz0);
      }}} // ix iy iz

      // see if these drivers compile and can be executed without segfaults
      if (echo > 11) std::printf("# here %s:%d\n", __func__, __LINE__);
      {
          auto const dVol = hGrid[2]*hGrid[1]*hGrid[0];
          auto const sho_norm = sho_normalization(lmax, sigma);
          // now sho_norm = product_d=0..2 sqrt(pi) * sigma * factorial(nu[d]) / 2^nu[d]
          for (int isho = 0; isho < std::min(nsho, 64); ++isho) {
              apc[isho*nrhs][0][0][isho] = dVol/sho_norm[isho]; // set "unit matrix" but normalized
          } // isho
      }

      SHOadd_driver<real_t,R1C2,Noco>(psi, apc, AtomPos, AtomLmax, AtomStarts, sparse_SHOadd.rowStart(), sparse_SHOadd.colIndex(), RowIndexCubes, ColIndexCubes, CubePos, hGrid, nnzb, nrhs, echo);
      if (0) {
          cudaDeviceSynchronize();
          size_t nz{0};
          for (int i = 0; i < nnzb*64; ++i) {
              auto const value = psi[i >> 6][0][i & 63][0];
              nz += (0 == value);
              if (echo > 19) std::printf(" %g", value);
          } // isho
          if (echo > 3) std::printf("\n# %ld non-zeros of %d\n", nz, nnzb*64);
      } // echo

      SHOprj_driver<real_t,R1C2,Noco>(apc, psi, AtomPos, AtomLmax, AtomStarts, natoms, sparse_SHOprj, RowIndexCubes, CubePos, hGrid, nrhs, echo);

      float maxdev[2] = {0, 0}; // {off-diagonal, diagonal}
      float maxdev_nu[8][2] = {{0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}};
      { // scope: show projection coefficients
          cudaDeviceSynchronize();
          std::vector<int8_t> nu_of_sho(nsho, -1);
          {
              int sho{0};
              for (int iz = 0; iz <= lmax; ++iz) {
                  for (int iy = 0; iy <= lmax - iz; ++iy) {
                      for (int ix = 0; ix <= lmax - iz - iy; ++ix) {
                          nu_of_sho[sho] = ix + iy + iz;
                          ++sho;
                      } // ix
                  } // iy
              } // iz
              assert(nu_of_sho.size() == sho);
          }
          auto const msho = std::min(nsho, 64);
          if (echo > 5) std::printf("# %d of %d projection coefficients ", msho, nsho);
          for (int isho = 0; isho < msho; ++isho) {
              if (echo > 9) std::printf("\n# projection coefficients[%2d]: ", isho);
              for (int jsho = 0; jsho < msho; ++jsho) {
                  double const value = apc[isho*nrhs][0][0][jsho];
                  int const diag = (isho == jsho); // unity
                  if (echo > 9) std::printf(" %.1e", value);
                  float const absdev = std::abs(value - diag);
                  maxdev[diag] = std::max(maxdev[diag], absdev);
                  auto const nu = std::max(nu_of_sho[isho], nu_of_sho[jsho]);
                  assert(nu >= 0);
                  if (nu < 8) maxdev_nu[nu][diag] = std::max(maxdev_nu[nu][diag], absdev);
              } // isho
              if (echo > 9) std::printf(" diagonal=");
              if (echo > 5) std::printf(" %.3f", apc[isho*nrhs][0][0][isho]);
          } // isho
          if (echo > 5) std::printf("\n");
      } // scope
      if (echo > 2) std::printf("# %s<%s,R1C2=%d,Noco=%d> orthogonality error %.2e, normalization error %.2e\n",
                                   __func__, real_t_name<real_t>(), R1C2, Noco, maxdev[0], maxdev[1]);
      if (echo > 5) {
          for (int nu = 0; nu < std::min(8, lmax + 1); ++nu) {
              std::printf("# %s  nu=%d  %.2e  %.2e\n", real_t_name<real_t>(), nu, maxdev_nu[nu][0], maxdev_nu[nu][1]);
          } // nu
      } // echo

      if (1) {
          // also test the deprecated interface 'multiply'
          auto AtomMatrices = get_memory<double*>(natoms, echo, "AtomMatrices");
          for (int ia = 0; ia < natoms; ++ia) {
              if (echo > 9) std::printf("# %s atom image #%d has lmax= %d and %d coefficients starting at %d\n", __func__, ia, AtomLmax[ia], nsho, AtomStarts[ia]);
              AtomMatrices[ia] = get_memory<double>(pow2(Noco)*2*pow2(nsho), echo, "AtomMatrix[ia]");
              set(AtomMatrices[ia], pow2(Noco)*2*pow2(nsho), 0.0);
          } // ia
          multiply<real_t,R1C2,Noco>(Vpsi, apc, psi, AtomPos, AtomLmax, AtomStarts, natoms, sparse_SHOprj, AtomMatrices, sparse_SHOadd, RowIndexCubes, ColIndexCubes, CubePos, hGrid, nnzb, nrhs, echo);
          for (int ia = 0; ia < natoms; ++ia) { free_memory(AtomMatrices[ia]); }
          free_memory(AtomMatrices);
      } // 0

      sparse_SHOprj[0].~sparse_t<>();
      free_memory(sparse_SHOprj);
      free_memory(ColIndexCubes);
      free_memory(RowIndexCubes);
      free_memory(CubePos);
      free_memory(AtomStarts);
      free_memory(AtomLmax);
      free_memory(AtomPos);
      free_memory(apc);
      free_memory(Vpsi);
      free_memory(psi);
      return 0;
  } // test_SHOprj_and_SHOadd

  status_t test_SHOprj_and_SHOadd(int const echo=0) {
      status_t stat(0);
      int const more = control::get("green_dyadic.test.more", 1.); // use 0...7
      if (more & 0x1) stat += test_SHOprj_and_SHOadd<float ,1,1>(echo); // real
      if (more & 0x2) stat += test_SHOprj_and_SHOadd<float ,2,1>(echo); // complex
      if (more & 0x4) stat += test_SHOprj_and_SHOadd<float ,2,2>(echo); // non-collinear
      if (more & 0x1) stat += test_SHOprj_and_SHOadd<double,1,1>(echo); // real
      if (more & 0x2) stat += test_SHOprj_and_SHOadd<double,2,1>(echo); // complex
      if (more & 0x4) stat += test_SHOprj_and_SHOadd<double,2,2>(echo); // non-collinear
      return stat;
  } // test_SHOprj_and_SHOadd


  template <int R1C2=2, int Noco=1>
  status_t test_SHOprj_right(int const echo=0, int8_t const lmax=5, double const sigma=.5) {
      if (echo > 0) std::printf("\n# %s<R1C2=%d,Noco=%d>\n", __func__, R1C2, Noco);
      uint32_t constexpr nAtoms = 2, nrhs = 1; // two atoms, atom#0 at origin, atom#1 at space diagonal of cube
      auto const nsho = sho_tools::nSHO(int(lmax));
      uint32_t const natomcoeffs = nAtoms*nsho;
      double const AtomPos[nAtoms][3+1] = {{4,4,4, 1./std::sqrt(sigma)}, {0,0,0, 1./std::sqrt(sigma)}};
      int8_t const AtomLmax[nAtoms] = {lmax, lmax};
      uint32_t const AtomStarts[nAtoms + 1] = {0, uint32_t(nsho), uint32_t(2*nsho)};
      uint32_t const AtomIndex[nAtoms] = {0, 1};
      double const AtomPhase[nAtoms][4] = {{1,0, 0,0}, {1,0, 0,0}};
      float  const colCubePos[nrhs][3+1] = {{0,0,0, 0}};
      double const hGrid[3+1] = {1,1,1, 6.28};
      auto pGreen = get_memory<double[R1C2][Noco][Noco*64]>(natomcoeffs*nrhs, echo, "pGreen");
      set(pGreen[0][0][0], natomcoeffs*nrhs*R1C2*Noco*Noco*64, 0.0);
      for (int i = 0; i < std::min(natomcoeffs*nrhs, Noco*Noco*64u); ++i) pGreen[i][R1C2 - 1][0][i] = 1; // diagonal matrix

      auto const result = SHOprj_right<R1C2,Noco>(pGreen, AtomPos, AtomLmax, AtomStarts, AtomIndex, AtomPhase, nAtoms, AtomLmax, nAtoms, colCubePos, hGrid, nrhs, echo);

      if (echo > 9) {
          for (int ia = 0; ia < nAtoms; ++ia) {
            for (int i = 0; i < nsho; ++i) {
              std::printf("# atom#%i %2i ", ia, i);
              for (int j = 0; j < nsho; ++j) {
                  std::printf(" %g", result[ia][(i*nsho + j)*Noco*Noco]);
              } // j
              std::printf("\n");
            } // i
          } // ia
      } // echo
      free_memory(pGreen);
      return 0;
  } // test_SHOprj_right

  status_t test_SHOprj_right(int const echo=0) {
      status_t stat(0);
      int const more = control::get("green_dyadic.test.right", 1.); // use 0...7
      if (more & 0x1) stat += test_SHOprj_right<1,1>(echo); // real
      if (more & 0x2) stat += test_SHOprj_right<2,1>(echo); // complex
      if (more & 0x4) stat += test_SHOprj_right<2,2>(echo); // non-collinear
      return stat;
  } // test_SHOprj_right

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_Hermite_polynomials_1D(echo);
 //   stat += test_SHOprj_and_SHOadd(echo);  // ToDo: fails with clang memory sanitizer
      stat += test_SHOprj_right(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_dyadic

//
// General thoughts:
//
//    We assume that a Green function G(r,r') comes in as psi.
//    The spatial arguments \vec r of the Green function are
//    sampled on an equidistant Cartesian real-space grids.
//    Always 4x4x4 grid points are grouped to a "cube".
//
//    The original codes benchmarked in the PASC19 paper by Baumeister & Tsukamoto
//    had the following data layout.
//
//        projection/addition coefficients layout[natoms*nSHO*nrhs][Noco*64]
//        wave functions                   layout[ncubes* 64 *nrhs][Noco*64]
//
//    However, we will have to reformulate to
//        projection/addition coefficients layout[natomcoeffs*nrhs][R1C2][Noco   ][Noco*64]
//        Green function                   layout[nnzb            ][R1C2][Noco*64][Noco*64]
//
//    where natomcoeffs == natoms*nSHO(numax) only if all atoms have the same numax
//    and   nnzb == ncubes*nrhs only if the Green function is dense.
//
//    For considerations of geometry(*), the coefficient matrix will probably stay dense (natomcoeffs \times nrhs)
//    Then, we can still  launch SHOprj <<< {nrhs, natoms, 1}, {Noco*64, R1C2*Noco, 1} >>> providing
//                                   bsr in RowStartAtoms[irhs][iatom +0 .. +1], icube = irhs_of_inzb[irhs][bsr] --> one sparse_t per irhs
//    But we will need to launch SHOadd <<< {nnzb,    1,   1}, {Noco*64, R1C2*Noco, 1} >>> providing irhs = colIndex[inzb],
//                                  irow = rowIndex[inzb] and iatom = ColIndexAtoms[irow][bsr], ColIndexAtoms as before
//
//    (*) The geometry consideration assumes the case that there is a compact cluster of
//        right hand side (source) cubes assigned to one MPI-process and isolated boundary conditions
//        for the cell. If the radius of target cubes is large compared to the cluster extent
//        the number of zero blocks inside the coefficient matrix is small compared to its size
//        since also each atomic projection region typically spans over several cubes.
//
