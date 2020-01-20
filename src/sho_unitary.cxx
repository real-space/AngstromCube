#include <cstdio> // printf
#include <cassert> // assert
#include <cmath> // sqrt, pow, exp, sqrt
#include <algorithm> // max, fill
#include <fstream> // std::ifstream 
#include <sstream> // std::istringstream
#include <vector> // std::vector
#include <numeric> // std::iota

#include "inline_tools.hxx" // align
#include "sho_tools.hxx" // n2HO
#include "sho_unitary.hxx"

namespace sho_unitary {

  template<typename real_t> 
  real_t signed_sqrt(real_t const x) { return (x < 0)? -std::sqrt(-x) : std::sqrt(x); }
  
  template<typename real_t> 
  status_t read_unitary_matrix_from_file(real_t **u, int const numax, int &nu_high, 
                  char const filename[], int const echo) {
      
      //
      // Expected File format:
      //    Each line has these 8 entries:
      //      nx ny nz ell emm nrn nom den
      //    where
      //      nx, ny, nz      are the three Cartesian SHO indices, >= 0
      //      ell, emm        are the spherical harmonic indices
      //                      ell >= 0, -ell <= emm <= ell
      //      nrn             is the number or radial nodes, nrn >= 0
      //      nom den         encode the value of the matrix entry of u:
      //                      u = sgn(nom)*sqrt(abs(nom)/den)
      //                      den > 0, nom may be negative but since
      //                      only non-zero entries are given, nom != 0
      //
      //  if we want to squeeze it into a data-type, this would do
      //  struct { uint8_t nx, ny, nz, _spare, ell, nrn; int16_t emm; 
      //           int64_t nom; uint64_t den; };
      //  this would support the index ranges up to numax=255
      std::ifstream infile(filename, std::ios_base::in);
      bool const file_is_nu_ordered = true;
        
      int n_ignored = 0;
              std::string line;
              while (std::getline(infile, line))
              {
                  char const c0 = line[0];
                  if ('#' != c0 && ' ' != c0 && '\n' != c0 && 0 != c0) { 
                      std::istringstream iss(line);
                      int nx, ny, nz, ell, emm, nrn;
                      int64_t nom, den;
                      if (!(iss >> nx >> ny >> nz >> ell >> emm >> nrn >> nom >> den)) {
                          printf("# Failed to read integer number from \"%s\"!\n", line.c_str());
                          break;
                      } // error

                      real_t const u_entry = signed_sqrt(nom/((real_t)den));
                      if (echo > 8) printf("%d %d %d    %d %2d %d  %.15f\n", nx, ny, nz, ell, emm, nrn, u_entry);
                      int const nzyx = sho_tools::Ezyx_index(nx, ny, nz);
                      int const nlnm = sho_tools::Elnm_index(ell, nrn, emm);
                      int const nu_xyz = sho_tools::get_nu(nx, ny, nz);
                      int const nu_rad = sho_tools::get_nu(ell, nrn);
                      if (nu_xyz != nu_rad) {
                          printf("# off-block entry found in file <%s>: nx=%d ny=%d nz=%d (nu=%d)  ell=%d emm=%d nrn=%d (nu=%d)\n", 
                                 filename, nx, ny, nz, nu_xyz, ell, emm, nrn, nu_rad);
                          return 1; // error
                      } // nu matches
                      int const nu = nu_xyz;
                      nu_high = std::max(nu_high, nu);
                      if (nu > numax) {
                          ++n_ignored; // ignore the entry
                          if (file_is_nu_ordered) return 0; // we can return already since the entries are nu-ordered, 
                                  // .. so we do not need to parse past the first nu which exceeds the searched range;
                      } else {
                          int const nb = sho_tools::n2HO(nu); // dimension of block                   
                          int const ioff = sho_tools::nSHO(nu - 1); // offset from previous blocks
                          u[nu][(nzyx - ioff)*nb + (nlnm - ioff)] = u_entry; // set the entry
                      } // nu in range
                  }
                  // process pair (a,b)
              } // while
        if (n_ignored && (echo > 2)) printf("# ignored %d lines in file <%s> reading up to nu=%d\n", n_ignored, filename, numax);
        return 0;
  } // read_unitary_matrix_from_file
  
  
  

#ifdef  NO_UNIT_TESTS
  template // explicit template instantiation
  status_t read_unitary_matrix_from_file<double>(double**, int const, int &, char const*, int const);

  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  int test_generation(int const echo=9) { return 0; } // return generate_unitary_transform(9, echo); }

  template<typename real_t>
  int test_loading(int const numax=9, int const echo=1) {
      Unitary_SHO_Transform<real_t> U(numax);
      auto const dev = U.test_unitarity(echo);
      if (echo > 0) printf("# Unitary_SHO_Transform<real_%ld>.test_unitarity = %g\n", sizeof(real_t), dev);
      return (dev > 2e-7);
  } // test_loading

  int test_vecor_transform(int const numax=3, int const echo=9) {
      Unitary_SHO_Transform<double> U(numax);
      if (echo > 0) printf("\n# %s %s(numax=%i, echo=%i)\n", __FILE__, __func__, numax, echo);
      int const nc = sho_tools::nSHO(numax);
      std::vector<double> vi(nc, 0), vo(nc, 0);
      std::iota(vi.begin(), vi.end(), 0);
//       return U.transform_vector(vo.data(), sho_tools::order_nlm, vi.data(), sho_tools::order_zyx, numax, echo);
      return U.transform_vector(vo.data(), sho_tools::order_Elnm, vi.data(), sho_tools::order_Ezyx, numax, echo);
  } // test_loading

  status_t all_tests() {
    auto status = 0;
    status += test_generation();
    status += test_loading<float>();
    status += test_loading<double>();
    status += test_vecor_transform();
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  
  
} // namespace sho_unitary


#if 0  
  
  status_t generate_unitary_transform(int const numax, int const echo=9) {
//       int const mx = align<2>(1 + numax);

      int const lmax = numax;
      int const lx = (1 + lmax);
      int xory[lmax + 1 + lmax][lx]; // xory[lmax + m][k] --> for m <  0: x^(|m|- k)*y^k if k even
                                     //                   --> for m >= 0: y^(|m|- k)*x^k if k odd
      for(int k = 0; k < (2*lmax + 1)*lx; ++k) xory[0][k] = 0; // clear
      
      int Pl[lx][lx]; // generate (r^2 - u^2)^l
      for(int k = 0; k < lx*lx; ++k) Pl[0][k] = 0; // clear
      
      int rxy2pl[lx][lx]; // generate (x^2 + y^2)^l
      for(int k = 0; k < lx*lx; ++k) rxy2pl[0][k] = 0; // clear
      
      { // bnc-scope
          int bnc[lx]; bnc[0] = 1; for(int l = 1; l <= lmax; ++l) bnc[l] = 0; // init binomial coefficients
          for(int m = 0; m <= lmax; ++m) {
              
              if (echo > 2) printf("# Pl l=%d  ", m);
              int j4 = 1;
              for(int k = 0; k <= m; ++k) {
                  Pl[m][k] = bnc[k]*(1 - 2*(k & 1));
                  rxy2pl[m][k] = bnc[k];
                  if (echo > 2) printf(" %d", Pl[m][k]);
                  int const j4mod2 = j4 % 2;
                  int const msgn = 2*j4mod2 - 1;
                  int const sgn = 1 + j4mod2 - j4;
                  xory[lmax - m*msgn][k] = bnc[k]*sgn*msgn;
                  j4 = (j4 + 1) % 4; // mimique the behaviour of complex i^k
              } // k
              if (echo > 2) printf("\n");
              
              // prepare bnc
              for(int k = m + 1; k > 0; --k) {
                  bnc[k] += bnc[k - 1];
    //               if (echo > 8) printf("%6d", bnc[k]);
              } // k
    //           if (echo > 8) printf("%6d\n", bnc[0]);
          } // m
      } // bnc-scope

      if (echo > 2) printf("\n# (x^2 + y^2)^%d has highest coefficient %d\n", lmax, rxy2pl[lmax][lmax/2]);
      
      if (echo > 1) { 
          printf("\n\n");
          for(int m = -lmax; m <= lmax; ++m) {
              if (echo > 2) printf("# m=%3d   ", m);
              for(int k = 0; k <= std::abs(m); ++k) {
                  auto const xy = xory[lmax + m][k];
                  if (echo > 2) printf("%6d", xy);
              } // k
              if (echo > 2) printf("\n");
          } // m
          if (echo > 2) printf("\n\n");
      } // echo
      
      int64_t Plm[lx][lx][2*lx]; // data-type needs to capture factorial(2*lmax)
      for(int k = 0; k < lx*lx*2*lx; ++k) Plm[0][0][k] = 0; // clear
      
      for(int l = 0; l <= lmax; ++l) {
          int64_t poly[2*l + 1]; // polynomial in u
          for(int k = 0; k <= 2*l; ++k) poly[k] = 0; // clear
          for(int k = 0; k <= l;   ++k) poly[2*k] = Pl[l][k]; // copy underived

          for(int m = -l; m <= l; ++m) { // start at -l to derive l times for pure Legendre polynomials
              if (m >= 0) {          // before valid m range starts for associated Legendre Polynomials
                  if (echo > 2) printf("# Plm l=%d m=%d  ", l, m);
                  for(int k = 0; k <= l - m; ++k) {
                      Plm[l][m][k] = poly[k]; // store associated Legendre polynomial
                      if (echo > 2) printf(" %ldx^%d ", Plm[l][m][k], k);
                  } // k
                  if (echo > 2) printf("\n");
              }
              // derive (r^2 - u^2)^l w.r.t. u one time, i.e. for l+m times
              for(int k = 1; k <= 2*l; ++k) poly[k - 1] = k*poly[k]; poly[2*l] = 0; // derive in-place
          } // m
          for(int k = 0; k <= 2*l; ++k) assert(0 == poly[k]); // all coefficients of poly must vanish since we derive u^{2l} for 2l times
      } // l
      
      
      
      
  } // generate_unitary_transform
  
#endif
