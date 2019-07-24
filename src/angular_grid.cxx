#include <cstdio> // printf
#include <cassert> // assert
#include <cmath> // std::sqrt, std::abs
#include <algorithm> // std::max
#include <vector> // std::vector

#include "angular_grid.hxx"
#include "angular_grid.h" // angular_grid_t
#include "inline_tools.hxx" // align<>
#include "solid_harmonics.hxx" // Xlm
//   #include "spherical_harmonics.hxx" // Ylm
#include "gaunt_entry.h" // gaunt_entry_t
#include "constants.hxx" // constants::pi

extern "C" {
   // BLAS interface to matrix matrix multiplication
  void dgemm_(const char*, const char*, const int*, const int*, const int*, const double*, 
              const double*, const int*, const double*, const int*, const double*, double*, const int*);
} // extern "C"


namespace angular_grid {

// #define LARGE_GRIDS
#ifdef  LARGE_GRIDS
  int constexpr ellmax_implemented = 65; // highest Lebedev-Laikov grid implemented
#else
  int constexpr ellmax_implemented = 20; // highest Lebedev-Laikov grid implemented
#endif


  template<> // template specialization
  void transform<double>(double *out, double const *in, int const M, // nrad is the stride for in[] and out[]
                         int const ellmax, bool const back, int const echo) {
      auto const g = get_grid(ellmax, echo);
      auto constexpr c = 'n'; double const w8 = 1, zero = 0;
      double *b; int ldb = 0, N = 0, K = 0;
      if (0 == back) {
          N = (1 + ellmax)*(1 + ellmax); K = g->npoints; b = g->Xlm2grid; ldb = g->Xlm2grid_stride;
      } else {
          K = (1 + ellmax)*(1 + ellmax); N = g->npoints; b = g->grid2Xlm; ldb = g->grid2Xlm_stride;
      } // back?
      dgemm_(&c, &c, &M, &N, &K, &w8, in, &M, b, &ldb, &zero, out, &M); // matrix-matrix multiplication with BLAS
  } // transform


  angular_grid_t* get_grid(int const ellmax, int const echo) {
    
      static angular_grid_t grids[1 + ellmax_implemented];
      
      if (ellmax < 0 || ellmax > ellmax_implemented) {
          printf("# %s: ellmax= %d is out of range [0, %d]\n", __func__, ellmax, ellmax_implemented);
          return nullptr;
      } // in range
      auto g = &grids[ellmax];
      if (g->npoints < 1 || g->ellmax != ellmax) {
          // init this instance
          g->ellmax = ellmax;
          g->npoints = Lebedev_grid_size(ellmax);
          int const nlm = (1 + ellmax)*(1 + ellmax);
          g->Xlm2grid_stride = align<2>(nlm); 
          g->grid2Xlm_stride = align<2>(g->npoints);
          // allocations
          g->xyzw     = new double[g->npoints][4];
          g->Xlm2grid = new double[g->npoints*g->Xlm2grid_stride];
          g->grid2Xlm = new double[nlm       *g->grid2Xlm_stride];
          
          auto const ist = create_Lebedev_grid(ellmax, g->xyzw);
          if (echo > 0 && ist) printf("# %s: angular grid for ellmax= %d failed with status %d\n", __func__, ellmax, ist);
          
          // clear the matrix memories
          for(int ij = 0; ij < g->npoints*g->Xlm2grid_stride; ++ij) g->Xlm2grid[ij] = 0;
          for(int ij = 0; ij < nlm       *g->grid2Xlm_stride; ++ij) g->grid2Xlm[ij] = 0;
          
          auto const xlm = new double[nlm];
          // create the real-valued spherical harmonics
          for(int ipt = 0; ipt < g->npoints; ++ipt) {
              auto const w8 = g->xyzw[ipt][3] * 4*constants::pi;
              solid_harmonics::Xlm(xlm, ellmax, g->xyzw[ipt]);
              for(int ilm = 0; ilm < nlm; ++ilm) {
                  g->Xlm2grid[ipt*g->Xlm2grid_stride + ilm] = xlm[ilm];
                  g->grid2Xlm[ilm*g->grid2Xlm_stride + ipt] = xlm[ilm]*w8;
              } // ilm
          } // ipt
          delete[] xlm;
          if (echo > 0) printf("# %s: angular grid for ellmax= %d is requested 1st: %d points\n", __func__, ellmax, g->npoints);
      } // grid was not set
      return g;

  } // get_grid
  
// ! cvw    icode=0:   (0,0,0) only                               ( 1 point)
// ! cvw    icode=1:   (0,0,1) etc                                ( 6 points)
// ! cvw    icode=2:   (0,a,a) etc, a=1/sqrt(2)                   (12 points)
// ! cvw    icode=3:   (a,a,a) etc, a=1/sqrt(3)                   ( 8 points)
// ! cvw    icode=4:   (a,a,b) etc, b=sqrt(1-2a^2)                (24 points)
// ! cvw    icode=5:   (a,b,0) etc, b=sqrt(1-a^2), a input        (24 points)
// ! cvw    icode=6:   (a,b,c) etc, c=sqrt(1-a^2-b^2), a&b input  (48 points)

  template<typename real_t>
  int gen_oh0(int ncode[], double const w8, real_t xyzw[][4]) {
      ++ncode[0];
      // start vector [ O, O, O ]
      // the point
      for(int dir = 0; dir < 3; ++dir) {
          xyzw[0][dir] = 0;
      } // dir
      xyzw[0][3] = w8;
// cTeXit '% no citation needed for using a single point instead of an angular grid.'
      return 1;
  } // gen_oh0
  
  template<typename real_t>
  int gen_oh1(int ncode[], double const w8, real_t xyzw[][4]) {
      ++ncode[1]; int i = 0;
      // start vector [O, O, A]
      // the octahedron
      for(int dir = 0; dir < 3; ++dir) {
          for(int s = -1; s <= 1; s += 2) {
              xyzw[i][0] = 0;
              xyzw[i][1] = 0;
              xyzw[i][2] = 0;
              xyzw[i][dir] = (real_t)s;
              xyzw[i][3] = w8;
              ++i;
          } // s
      } // dir
      assert(6 == i);
      return 6;
// ! chvd
// ! chvd   [1] V.I. Lebedev, and D.N. Laikov
// ! chvd       "A quadrature formula for the sphere of the 131st algebraic order of accuracy"
// ! chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
// ! chvd
// cTeX1c bib_entry('Lebedev99', 'article', author='V.I. Lebedev, and D.N. Laikov', &
// cTeX     title='A quadrature formula for the sphere of the 131st algebraic order of accuracy', &
// cTeX     journal='Doklady Mathematics', volume='59', number='3', year='1999', pages='477-481')
  } // gen_oh1
  
  template<typename real_t>
  int gen_oh2(int ncode[], double const w8, real_t xyzw[][4]) {
      ++ncode[2]; int i = 0;
      // start vector [O, A, A]
      double const sq2 = std::sqrt(.5);
      // the dodecahedron
      for(int dir = 0; dir < 3; ++dir) {
          for(int s1 = -1; s1 <= 1; s1 += 2) {
              for(int s2 = -1; s2 <= 1; s2 += 2) {
                  xyzw[i][ dir       ] = 0;
                  xyzw[i][(dir + 1)%3] = (real_t)(s1*sq2);
                  xyzw[i][(dir + 2)%3] = (real_t)(s2*sq2);
                  xyzw[i][3] = w8;
                  ++i;
              } // s2
          } // s1
      } // dir
      assert(12 == i);
      return 12;
// ! chvd
// ! chvd   [2] V.I. Lebedev
// ! chvd       "A quadrature formula for the sphere of 59th algebraic order of accuracy"
// ! chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
// ! chvd
// cTeX1c bib_entry('Lebedev95', 'article', author='V.I. Lebedev', &
// cTeX     title='A quadrature formula for the sphere of 59th algebraic order of accuracy', &
// cTeX     journal='Russian Acad. Sci. Dokl. Math.', volume='50', year='1999', pages='283-286')

  } // gen_oh2
  
  template<typename real_t>
  int gen_oh3(int ncode[], double const w8, real_t xyzw[][4]) {
      ++ncode[3]; int i = 0;
      // start vector [A, A, A]
      double const sq3 = sqrt(1/3.);
      // the cube
      for(int s3 = -1; s3 < 2; s3 += 2) {
          for(int s2 = -1; s2 < 2; s2 += 2) {
              for(int s1 = -1; s1 < 2; s1 += 2) {
                  xyzw[i][0] = (real_t)(sq3*s1);
                  xyzw[i][1] = (real_t)(sq3*s2);
                  xyzw[i][2] = (real_t)(sq3*s3);
                  xyzw[i][3] = w8;
                  ++i;
              } // s1
          } // s2
      } // s3
      assert(8 == i);
      return 8;
// ! chvd
// ! chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
// ! chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
// ! chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
// ! chvd
// cTeX1c bib_entry('Lebedev92', 'article', author='V.I. Lebedev', &
// cTeX     title='Quadrature formulas of orders 41, 47, and 53 for the sphere', &
// cTeX     journal='Russian Acad. Sci. Dokl. Math.', volume='45', year='1992', pages='587-592')
  } // gen_oh3
  
  template<typename real_t>
  int gen_oh4(int ncode[], double const w8, real_t xyzw[][4], double const aa) {
      ++ncode[4]; int i = 0;
      // start vector [A, A, B]
      assert(2*aa*aa < 1);
//       A = aa ; if(2*A*A >= 1) die_here('gen_oh: break; case 4: |AA| is too large!')
      double const bb = std::sqrt(1. - 2*aa*aa);
      for(int dir = 0; dir < 3; ++dir) {
          for(int s3 = -1; s3 < 2; s3 += 2) {
              for(int s2 = -1; s2 < 2; s2 += 2) {
                  for(int s1 = -1; s1 < 2; s1 += 2) {
                      xyzw[i][ dir       ] = (real_t)(bb*s3);
                      xyzw[i][(dir + 1)%3] = (real_t)(aa*s1);
                      xyzw[i][(dir + 2)%3] = (real_t)(aa*s2);
                      xyzw[i][3] = w8;
                      ++i;
                  } // s1
              } // s2
          } // s3
      } // dir
      assert(24 == i);
      return 24;
// ! chvd
// ! chvd   [4] V.I. Lebedev
// ! chvd       "Spherical quadrature formulas exact to orders 25-29"
// ! chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
// ! chvd
// cTeX1c bib_entry('Lebedev77', 'article', author='V.I. Lebedev', &
// cTeX     title='Spherical quadrature formulas exact to orders 25-29', &
// cTeX     journal='Siberian Mathematical Journal', volume='18', year='1977', pages='99-107')
  } // gen_oh4
  
  template<typename real_t>
  int gen_oh5(int ncode[], double const w8, real_t xyzw[][4], double const aa) {
      ++ncode[5]; int i = 0;
      // start vector [A, B, O]
      assert(aa*aa < 1);
//       A = aa ; if(A*A >= 1) die_here('gen_oh: break; case 5: |AA| is too large!')
      double const bb = std::sqrt(1. - aa*aa);
      double b3 = bb, a3 = aa;
      for(int s3 = 0; s3 < 2; ++s3) {
          for(int dir = 0; dir < 3; ++dir) {
              for(int s2 = -1; s2 < 2; s2 += 2) {
                  for(int s1 = -1; s1 < 2; s1 += 2) {
                      xyzw[i][ dir       ] = 0;
                      xyzw[i][(dir + 1)%3] = (real_t)(a3*s1);
                      xyzw[i][(dir + 2)%3] = (real_t)(b3*s2);
                      xyzw[i][3] = w8;
                      ++i;
                  } // s1
              } // s2
          } // dir
          b3 = aa; a3 = bb; // interchange
      } // s3
      assert(24 == i);
      return 24;
// ! chvd
// ! chvd   [5] V.I. Lebedev
// ! chvd       "Quadratures on a sphere"
// ! chvd       Computational Mathematics and Mathematical Physics, Vol. 16, 1976, pp. 10-24. 
// ! chvd
// cTeX1c bib_entry('Lebedev76', 'article', author='V.I. Lebedev', title='Quadratures on a sphere', &
// cTeX     journal='Computational Mathematics and Mathematical Physics', volume='16', year='1976', pages='10-24')
  } // gen_oh5
  
  template<typename real_t>
  int gen_oh6(int ncode[], double const w8, real_t xyzw[][4], double const aa, double const bb) {
      ++ncode[6]; int i = 0;
      // start vector [A, B, C]
//       A = aa ; if(A*A >= 1) die_here('gen_oh: break; case 6: |AA| is too large!')
      assert(aa*aa < 1);
//       B = bb ; if(B*B >= 1) die_here('gen_oh: break; case 6: |BB| is too large!')
      assert(bb*bb < 1);
//       if(B*B + A*A >= 1) die_here('gen_oh: break; case 6: |BB| is too large for this |AA|!')
      assert(aa*aa + bb*bb < 1);

      double const cc = std::sqrt(1. - aa*aa - bb*bb);
      double const c3[2][5] = { {aa, bb, cc, aa, bb},
                                {bb, aa, cc, bb, aa} };
      for(int rev = 0; rev < 2; ++rev) {
          for(int dir = 0; dir < 3; ++dir) {
              for(int s3 = -1; s3 < 2; s3 += 2) {
                  for(int s2 = -1; s2 < 2; s2 += 2) {
                      for(int s1 = -1; s1 < 2; s1 += 2) {
                          xyzw[i][0] = (real_t)(c3[rev][dir + 0]*s1);
                          xyzw[i][1] = (real_t)(c3[rev][dir + 1]*s2);
                          xyzw[i][2] = (real_t)(c3[rev][dir + 2]*s3);
                          xyzw[i][3] = w8;
                          ++i;
                      } // s1
                  } // s2
              } // s3
          } // dir
      } // rev
      assert(48 == i);
      return 48;
// ! chvd
// ! chvd   [6] V.I. Lebedev
// ! chvd       "Values of the nodes and weights of ninth to seventeenth order Gauss-Markov quadrature formulae invariant under the octahedron group with inversion"
// ! chvd       Computational Mathematics and Mathematical Physics, Vol. 15, 1975, pp. 44-51.
// ! chvd
// cTeX1c bib_entry('Lebedev75', 'article', author='V.I. Lebedev', &
// cTeX     title='Values of the nodes and weights of ninth to seventeenth order Gauss-Markov quadrature formulae invariant under the octahedron group with inversion', &
// cTeX     journal='Computational Mathematics and Mathematical Physics', volume='15', year='1975', pages='44-51')
  } // gen_oh6
  
  int Lebedev_grid_size(int const ellmax, int echo) {
    if (ellmax < -1) return -1;
    
    int const available[34] = {0,1,6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,
    434,590,770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810};
    // for ellmax = 65, we reach 5808.0 which then is corrected to 5810

    int n = (4*(ellmax + 1)*(ellmax + 1))/3; // to be corrected to next larger numbers of grid points available
    // todo: check if the suggestion n = (ellmax+1)^2 leads to orthogonality of the spherical harmonics
    int i = 0; while (available[i] < n  && i < 33) { ++i; }
    if (i >= 32 && available[i] < n) {
        if (echo > 0) printf("# %s ellmax= %d leads to n= %d (too large)!\n", __func__, ellmax, n);
        return -1;
    }
    if (echo > 1) printf("# %s correct the number of angular grid points of the Lebedev-Laikov grid"
                         " for ellmax= %d from %d to %d\n", __func__, ellmax, n, available[i]);
    return available[i];
  } // Lebedev_grid_size

  
  template<typename real_t>
  status_t create_Lebedev_grid(int const ellmax, real_t xyzw[][4], int echo) {

    int nc[8] = {0,0,0,0, 0,0,0,0}; // init
    int m = 0; // init
    
    int const n_corrected = Lebedev_grid_size(ellmax, echo);
    switch (n_corrected) {
           case   0:
      // do nothing, in particular do not access xyzw as it will have zero size
    break; case   1:
      m += gen_oh0(nc, 1., xyzw + m);
    break; case   6:
      m += gen_oh1(nc, 1./6., xyzw + m);
    break; case  14:
      m += gen_oh1(nc, 2./30., xyzw + m);
      m += gen_oh3(nc, 3./40., xyzw + m);
    break; case  26:
      m += gen_oh1(nc, .04761904761904762, xyzw + m);
      m += gen_oh2(nc, .03809523809523810, xyzw + m);
      m += gen_oh3(nc, .03214285714285714, xyzw + m);
    break; case  38:
      m += gen_oh1(nc,.009523809523809524, xyzw + m);
      m += gen_oh3(nc, .03214285714285714, xyzw + m);
      m += gen_oh5(nc, .02857142857142857, xyzw + m, .4597008433809831);
    break; case  50:
      m += gen_oh1(nc, .01269841269841270, xyzw + m);
      m += gen_oh2(nc, .02257495590828924, xyzw + m);
      m += gen_oh3(nc, .02109375000000000, xyzw + m);
      m += gen_oh4(nc, .02017333553791887, xyzw + m, .3015113445777636);
    break; case  74:
      m += gen_oh1(nc, .00051306717973385, xyzw + m);
      m += gen_oh2(nc, .01660406956574204, xyzw + m);
      m += gen_oh3(nc,-.02958603896103896, xyzw + m);
      m += gen_oh4(nc, .02657620708215946, xyzw + m, .4803844614152614);
      m += gen_oh5(nc, .01652217099371571, xyzw + m, .3207726489807764);
    break; case  86:
      m += gen_oh1(nc, .01154401154401154, xyzw + m);
      m += gen_oh3(nc, .01194390908585628, xyzw + m);
      m += gen_oh4(nc, .01111055571060340, xyzw + m, .3696028464541502);
      m += gen_oh4(nc, .01187650129453714, xyzw + m, .6943540066026664);
      m += gen_oh5(nc, .01181230374690448, xyzw + m, .3742430390903412);
    break; case 110:
      m += gen_oh1(nc, .003828270494937162, xyzw + m);
      m += gen_oh3(nc, .009793737512487512, xyzw + m);
      m += gen_oh4(nc, .008211737283191111, xyzw + m, .1851156353447362);
      m += gen_oh4(nc, .009942814891178103, xyzw + m, .6904210483822922);
      m += gen_oh4(nc, .009595471336070963, xyzw + m, .3956894730559419);
      m += gen_oh5(nc, .009694996361663028, xyzw + m, .4783690288121502);
    break; case 146:
      m += gen_oh1(nc, .5996313688621381E-3, xyzw + m);
      m += gen_oh2(nc, .7372999718620756E-2, xyzw + m);
      m += gen_oh3(nc, .7210515360144488E-2, xyzw + m);
      m += gen_oh4(nc, .7116355493117555E-2, xyzw + m, .6764410400114264);
      m += gen_oh4(nc, .6753829486314477E-2, xyzw + m, .4174961227965453);
      m += gen_oh4(nc, .7574394159054034E-2, xyzw + m, .1574676672039082);
      m += gen_oh6(nc, .6991087353303262E-2, xyzw + m, .1403553811713183, .4493328323269557);
    break; case 170:
      m += gen_oh1(nc, .5544842902037365E-2, xyzw + m);
      m += gen_oh2(nc, .6071332770670752E-2, xyzw + m);
      m += gen_oh3(nc, .6383674773515093E-2, xyzw + m);
      m += gen_oh4(nc, .5183387587747790E-2, xyzw + m, .2551252621114134);
      m += gen_oh4(nc, .6317929009813725E-2, xyzw + m, .6743601460362766);
      m += gen_oh4(nc, .6201670006589077E-2, xyzw + m, .4318910696719410);
      m += gen_oh5(nc, .5477143385137348E-2, xyzw + m, .2613931360335988);
      m += gen_oh6(nc, .5968383987681156E-2, xyzw + m, .4990453161796037, .1446630744325115);
    break; case 194:
      m += gen_oh1(nc, .1782340447244611E-2, xyzw + m);
      m += gen_oh2(nc, .5716905949977102E-2, xyzw + m);
      m += gen_oh3(nc, .5573383178848738E-2, xyzw + m);
      m += gen_oh4(nc, .5608704082587997E-2, xyzw + m, .6712973442695226);
      m += gen_oh4(nc, .5158237711805383E-2, xyzw + m, .2892465627575439);
      m += gen_oh4(nc, .5518771467273614E-2, xyzw + m, .4446933178717437);
      m += gen_oh4(nc, .4106777028169394E-2, xyzw + m, .1299335447650067);
      m += gen_oh5(nc, .5051846064614808E-2, xyzw + m, .3457702197611283);
      m += gen_oh6(nc, .5530248916233094E-2, xyzw + m, .1590417105383530, .8360360154824589);
    break; case 230:
      m += gen_oh1(nc, -.5522639919727325E-1, xyzw + m);;                                                     
      m += gen_oh3(nc, .4450274607445226E-2, xyzw + m);
      m += gen_oh4(nc, .4496841067921404E-2, xyzw + m, .4492044687397611);
      m += gen_oh4(nc, .5049153450478750E-2, xyzw + m, .2520419490210201);
      m += gen_oh4(nc, .3976408018051883E-2, xyzw + m, .6981906658447242);
      m += gen_oh4(nc, .4401400650381014E-2, xyzw + m, .6587405243460960);
      m += gen_oh4(nc, .1724544350544401E-1, xyzw + m, .0403854405009766);
      m += gen_oh5(nc, .4231083095357343E-2, xyzw + m, .5823842309715584);
      m += gen_oh5(nc, .5198069864064399E-2, xyzw + m, .3545877390518688);
      m += gen_oh6(nc, .4695720972568883E-2, xyzw + m, .2272181808998187, .4864661535886647);
    break; case 266:
      m += gen_oh1(nc, -.1313769127326952E-2, xyzw + m);
      m += gen_oh2(nc, -.2522728704859336E-2, xyzw + m);
      m += gen_oh3(nc, .4186853881700583E-2, xyzw + m);
      m += gen_oh4(nc, .5315167977810885E-2, xyzw + m, .7039373391585475);
      m += gen_oh4(nc, .4047142377086219E-2, xyzw + m, .1012526248572414);
      m += gen_oh4(nc, .4112482394406990E-2, xyzw + m, .4647448726420539);
      m += gen_oh4(nc, .3595584899758782E-2, xyzw + m, .3277420654971629);
      m += gen_oh4(nc, .4256131351428158E-2, xyzw + m, .6620338663699974);
      m += gen_oh5(nc, .4229582700647240E-2, xyzw + m, .8506508083520399);
      m += gen_oh6(nc, .4080914225780505E-2, xyzw + m, .3233484542692899, .1153112011009701);
      m += gen_oh6(nc, .4071467593830964E-2, xyzw + m, .2314790158712601, .5244939240922365);
    break; case 302:
      m += gen_oh1(nc, .8545911725128148E-3, xyzw + m);
      m += gen_oh3(nc, .3599119285025571E-2, xyzw + m);
      m += gen_oh4(nc, .3449788424305883E-2, xyzw + m, .3515640345570105);
      m += gen_oh4(nc, .3604822601419882E-2, xyzw + m, .6566329410219612);
      m += gen_oh4(nc, .3576729661743367E-2, xyzw + m, .4729054132581005);
      m += gen_oh4(nc, .2352101413689164E-2, xyzw + m, .0961830852261478);
      m += gen_oh4(nc, .3108953122413675E-2, xyzw + m, .2219645236294178);
      m += gen_oh4(nc, .3650045807677255E-2, xyzw + m, .7011766416089545);
      m += gen_oh5(nc, .2982344963171804E-2, xyzw + m, .2644152887060663);
      m += gen_oh5(nc, .3600820932216460E-2, xyzw + m, .5718955891878961);
      m += gen_oh6(nc, .3571540554273387E-2, xyzw + m, .2510034751770465, .8000727494073951);
      m += gen_oh6(nc, .3392312205006170E-2, xyzw + m, .1233548532583327, .4127724083168531);
    break; case 350:
      m += gen_oh1(nc, .3006796749453936E-2, xyzw + m);
      m += gen_oh3(nc, .3050627745650771E-2, xyzw + m);
      m += gen_oh4(nc, .1621104600288991E-2, xyzw + m, .7068965463912316);
      m += gen_oh4(nc, .3005701484901752E-2, xyzw + m, .4794682625712025);
      m += gen_oh4(nc, .2990992529653774E-2, xyzw + m, .1927533154878019);
      m += gen_oh4(nc, .2982170644107595E-2, xyzw + m, .6930357961327123);
      m += gen_oh4(nc, .2721564237310992E-2, xyzw + m, .3608302115520091);
      m += gen_oh4(nc, .3033513795811141E-2, xyzw + m, .6498486161496169);
      m += gen_oh5(nc, .3007949555218533E-2, xyzw + m, .1932945013230339);
      m += gen_oh5(nc, .2881964603055307E-2, xyzw + m, .3800494919899303);
      m += gen_oh6(nc, .2958357626535696E-2, xyzw + m, .2899558825499574, .7934537856582315);
      m += gen_oh6(nc, .3036020026407088E-2, xyzw + m, .0968412145510396, .8280801506686862);
      m += gen_oh6(nc, .2832187403926303E-2, xyzw + m, .1833434647041659, .9074658265305127);
    break; case 434:
      m += gen_oh1(nc, .5265897968224436E-3, xyzw + m);
      m += gen_oh2(nc, .2548219972002607E-2, xyzw + m);
      m += gen_oh3(nc, .2512317418927307E-2, xyzw + m);
      m += gen_oh4(nc, .2530403801186355E-2, xyzw + m, .6909346307509111);
      m += gen_oh4(nc, .2014279020918528E-2, xyzw + m, .1774836054609158);
      m += gen_oh4(nc, .2501725168402936E-2, xyzw + m, .4914342637784746);
      m += gen_oh4(nc, .2513267174597564E-2, xyzw + m, .6456664707424256);
      m += gen_oh4(nc, .2302694782227416E-2, xyzw + m, .2861289010307638);
      m += gen_oh4(nc, .1462495621594614E-2, xyzw + m, .0756808436717802);
      m += gen_oh4(nc, .2445373437312980E-2, xyzw + m, .3927259763368002);
      m += gen_oh5(nc, .2417442375638981E-2, xyzw + m, .8818132877794288);
      m += gen_oh5(nc, .1910951282179532E-2, xyzw + m, .9776428111182649);
      m += gen_oh6(nc, .2416930044324775E-2, xyzw + m, .2054823696403044, .8689460322872412);
      m += gen_oh6(nc, .2512236854563495E-2, xyzw + m, .5905157048925271, .7999278543857286);
      m += gen_oh6(nc, .2496644054553086E-2, xyzw + m, .5550152361076807, .7717462626915901);
      m += gen_oh6(nc, .2236607760437849E-2, xyzw + m, .9371809858553722, .3344363145343455);
    break; case 590:
      m += gen_oh1(nc, .3095121295306187E-3, xyzw + m);
      m += gen_oh3(nc, .1852379698597489E-2, xyzw + m);
      m += gen_oh4(nc, .1871790639277744E-2, xyzw + m, .7040954938227469);
      m += gen_oh4(nc, .1858812585438317E-2, xyzw + m, .6807744066455244);
      m += gen_oh4(nc, .1852028828296213E-2, xyzw + m, .6372546939258752);
      m += gen_oh4(nc, .1846715956151242E-2, xyzw + m, .5044419707800358);
      m += gen_oh4(nc, .1818471778162769E-2, xyzw + m, .4215761784010967);
      m += gen_oh4(nc, .1749564657281154E-2, xyzw + m, .3317920736472123);
      m += gen_oh4(nc, .1617210647254411E-2, xyzw + m, .2384736701421887);
      m += gen_oh4(nc, .1384737234851692E-2, xyzw + m, .1459036449157763);
      m += gen_oh4(nc, .9764331165051050E-3, xyzw + m, .0609503411550720);
      m += gen_oh5(nc, .1857161196774078E-2, xyzw + m, .6116843442009876);
      m += gen_oh5(nc, .1705153996395864E-2, xyzw + m, .3964755348199858);
      m += gen_oh5(nc, .1300321685886048E-2, xyzw + m, .1724782009907724);
      m += gen_oh6(nc, .1842866472905286E-2, xyzw + m, .5610263808622060, .3518280927733519);
      m += gen_oh6(nc, .1802658934377451E-2, xyzw + m, .4742392842551980, .2634716655937950);
      m += gen_oh6(nc, .1849830560443660E-2, xyzw + m, .5984126497885380, .1816640840360209);
      m += gen_oh6(nc, .1713904507106709E-2, xyzw + m, .3791035407695563, .1720795225656878);
      m += gen_oh6(nc, .1555213603396808E-2, xyzw + m, .2778673190586244, .0821302158193251);
      m += gen_oh6(nc, .1802239128008525E-2, xyzw + m, .5033564271075117, .0899920584207488);
#ifdef  LARGE_GRIDS
    break; case 770:
      m += gen_oh1(nc, .2192942088181184E-3, xyzw + m);
      m += gen_oh2(nc, .1436433617319080E-2, xyzw + m);
      m += gen_oh3(nc, .1421940344335877E-2, xyzw + m);
      m += gen_oh4(nc, .6798123511050502E-3, xyzw + m, .0508720441050236);
      m += gen_oh4(nc, .9913184235294911E-3, xyzw + m, .1228198790178831);
      m += gen_oh4(nc, .1180207833238949E-2, xyzw + m, .2026890814408786);
      m += gen_oh4(nc, .1296599602080921E-2, xyzw + m, .2847745156464294);
      m += gen_oh4(nc, .1365871427428316E-2, xyzw + m, .3656719078978026);
      m += gen_oh4(nc, .1402988604775325E-2, xyzw + m, .4428264886713469);
      m += gen_oh4(nc, .1418645563595609E-2, xyzw + m, .5140619627249735);
      m += gen_oh4(nc, .1421376741851662E-2, xyzw + m, .6306401219166803);
      m += gen_oh4(nc, .1423996475490962E-2, xyzw + m, .6716883332022612);
      m += gen_oh4(nc, .1431554042178567E-2, xyzw + m, .6979792685336881);
      m += gen_oh5(nc, .9254401499865368E-3, xyzw + m, .1446865674195309);
      m += gen_oh5(nc, .1250239995053509E-2, xyzw + m, .3390263475411216);
      m += gen_oh5(nc, .1394365843329230E-2, xyzw + m, .5335804651263506);
      m += gen_oh6(nc, .1127089094671749E-2, xyzw + m, .0694402439334941, .2355187894242326);
      m += gen_oh6(nc, .1345753760910670E-2, xyzw + m, .2269004109529460, .4102182474045730);
      m += gen_oh6(nc, .1424957283316783E-2, xyzw + m, .0802557460777534, .6214302417481605);
      m += gen_oh6(nc, .1261523341237750E-2, xyzw + m, .1467999527896572, .3245284345717394);
      m += gen_oh6(nc, .1392547106052696E-2, xyzw + m, .1571507769824727, .5224482189696630);
      m += gen_oh6(nc, .1418761677877656E-2, xyzw + m, .2365702993157246, .6017546634089558);
      m += gen_oh6(nc, .1338366684479554E-2, xyzw + m, .0771481586676573, .4346575516141163);
      m += gen_oh6(nc, .1393700862676131E-2, xyzw + m, .3062936666210730, .4908826589037616);
      m += gen_oh6(nc, .1415914757466932E-2, xyzw + m, .3822477379524787, .5648768149099500);
    break; case 974:
      m += gen_oh1(nc, .1438294190527431E-3, xyzw + m);
      m += gen_oh3(nc, .1125772288287004E-2, xyzw + m);
      m += gen_oh4(nc, .4948029341949241E-3, xyzw + m, .0429296354534135);
      m += gen_oh4(nc, .7357990109125470E-3, xyzw + m, .1051426854086404);
      m += gen_oh4(nc, .8889132771304384E-3, xyzw + m, .1750024867623087);
      m += gen_oh4(nc, .9888347838921435E-3, xyzw + m, .2477653379650257);
      m += gen_oh4(nc, .1053299681709471E-2, xyzw + m, .3206567123955957);
      m += gen_oh4(nc, .1092778807014578E-2, xyzw + m, .3916520749849983);
      m += gen_oh4(nc, .1114389394063227E-2, xyzw + m, .4590825874187624);
      m += gen_oh4(nc, .1123724788051555E-2, xyzw + m, .5214563888415861);
      m += gen_oh4(nc, .1125239325243814E-2, xyzw + m, .6253170244654199);
      m += gen_oh4(nc, .1126153271815905E-2, xyzw + m, .6637926744523170);
      m += gen_oh4(nc, .1130286931123841E-2, xyzw + m, .6910410398498301);
      m += gen_oh4(nc, .1134986534363955E-2, xyzw + m, .7052907007457760);
      m += gen_oh5(nc, .6823367927109931E-3, xyzw + m, .1236686762657990);
      m += gen_oh5(nc, .9454158160447096E-3, xyzw + m, .2940777114468387);
      m += gen_oh5(nc, .1074429975385679E-2, xyzw + m, .4697753849207649);
      m += gen_oh5(nc, .1129300086569132E-2, xyzw + m, .6334563241139567);
      m += gen_oh6(nc, .8436884500901954E-3, xyzw + m, .0597404861418134, .2029128752777523);
      m += gen_oh6(nc, .1075255720448885E-2, xyzw + m, .1375760408473636, .4602621942484054);
      m += gen_oh6(nc, .1108577236864462E-2, xyzw + m, .3391016526336286, .5030673999662036);
      m += gen_oh6(nc, .9566475323783357E-3, xyzw + m, .1271675191439820, .2817606422442134);
      m += gen_oh6(nc, .1080663250717391E-2, xyzw + m, .2693120740413512, .4331561291720157);
      m += gen_oh6(nc, .1126797131196295E-2, xyzw + m, .1419786452601918, .6256167358580814);
      m += gen_oh6(nc, .1022568715358061E-2, xyzw + m, .0670928460073826, .3798395216859157);
      m += gen_oh6(nc, .1108960267713108E-2, xyzw + m, .0705773818325617, .5517505421423520);
      m += gen_oh6(nc, .1122790653435766E-2, xyzw + m, .2783888477882155, .6029619156159187);
      m += gen_oh6(nc, .1032401847117460E-2, xyzw + m, .1979578938917407, .3589606329589096);
      m += gen_oh6(nc, .1107249382283854E-2, xyzw + m, .2087307061103274, .5348666438135476);
      m += gen_oh6(nc, .1121780048519972E-2, xyzw + m, .4055122137872836, .5674997546074373);
    break; case 1202:
      m += gen_oh1(nc, .1105189233267572E-3, xyzw + m);
      m += gen_oh2(nc, .9205232738090741E-3, xyzw + m);
      m += gen_oh3(nc, .9133159786443561E-3, xyzw + m);
      m += gen_oh4(nc, .3690421898017899E-3, xyzw + m, .0371263644965709);
      m += gen_oh4(nc, .5603990928680660E-3, xyzw + m, .0914006041226222);
      m += gen_oh4(nc, .6865297629282609E-3, xyzw + m, .1531077852469906);
      m += gen_oh4(nc, .7720338551145630E-3, xyzw + m, .2180928891660612);
      m += gen_oh4(nc, .8301545958894795E-3, xyzw + m, .2839874532200175);
      m += gen_oh4(nc, .8686692550179628E-3, xyzw + m, .3491177600963764);
      m += gen_oh4(nc, .8927076285846890E-3, xyzw + m, .4121431461444309);
      m += gen_oh4(nc, .9060820238568219E-3, xyzw + m, .4718993627149127);
      m += gen_oh4(nc, .9119777254940867E-3, xyzw + m, .5273145452842337);
      m += gen_oh4(nc, .9128720138604181E-3, xyzw + m, .6209475332444019);
      m += gen_oh4(nc, .9130714935691735E-3, xyzw + m, .6569722711857291);
      m += gen_oh4(nc, .9152873784554116E-3, xyzw + m, .6841788309070143);
      m += gen_oh4(nc, .9187436274321654E-3, xyzw + m, .7012604330123631);
      m += gen_oh5(nc, .5176977312965694E-3, xyzw + m, .1072382215478166);
      m += gen_oh5(nc, .7331143682101417E-3, xyzw + m, .2582068959496968);
      m += gen_oh5(nc, .8463232836379928E-3, xyzw + m, .4172752955306717);
      m += gen_oh5(nc, .9031122694253992E-3, xyzw + m, .5700366911792503);
      m += gen_oh6(nc, .6485778453163257E-3, xyzw + m, .9827986018263947, .1771774022615325);
      m += gen_oh6(nc, .7435030910982369E-3, xyzw + m, .9624249230326228, .2475716463426288);
      m += gen_oh6(nc, .7998527891839054E-3, xyzw + m, .9402007994128811, .3354616289066489);
      m += gen_oh6(nc, .8101731497468018E-3, xyzw + m, .9320822040143202, .3173615246611977);
      m += gen_oh6(nc, .8483389574594330E-3, xyzw + m, .9043674199393299, .4090268427085357);
      m += gen_oh6(nc, .8556299257311812E-3, xyzw + m, .8912407560074747, .3854291150669224);
      m += gen_oh6(nc, .8803208679738260E-3, xyzw + m, .8676435628462708, .4932221184851285);
      m += gen_oh6(nc, .8811048182425720E-3, xyzw + m, .8581979986041619, .4785320675922435);
      m += gen_oh6(nc, .8850282341265444E-3, xyzw + m, .8396753624049856, .4507422593157064);
      m += gen_oh6(nc, .9021342299040653E-3, xyzw + m, .8165288564022188, .5632123020762100);
      m += gen_oh6(nc, .9010091677105086E-3, xyzw + m, .8015469370783529, .5434303569693900);
      m += gen_oh6(nc, .9022692938426915E-3, xyzw + m, .7773563069070351, .5123518486419871);
      m += gen_oh6(nc, .9158016174693465E-3, xyzw + m, .7661621213900394, .6394279634749102);
      m += gen_oh6(nc, .9131578003189435E-3, xyzw + m, .7553584143533510, .6269805509024392);
      m += gen_oh6(nc, .9107813579482705E-3, xyzw + m, .7344305757559503, .6031161693096310);
      m += gen_oh6(nc, .9105760258970126E-3, xyzw + m, .7043837184021765, .5693702498468441);
    break; case 1454:
      m += gen_oh1(nc, .7777160743261247E-4, xyzw + m);
      m += gen_oh3(nc, .7557646413004701E-3, xyzw + m);
      m += gen_oh4(nc, .2841633806090617E-3, xyzw + m, .0322929066341385);
      m += gen_oh4(nc, .4374419127053555E-3, xyzw + m, .0803673327146222);
      m += gen_oh4(nc, .5417174740872172E-3, xyzw + m, .1354289960531653);
      m += gen_oh4(nc, .6148000891358593E-3, xyzw + m, .1938963861114426);
      m += gen_oh4(nc, .6664394485800704E-3, xyzw + m, .2537343715011275);
      m += gen_oh4(nc, .7025039356923220E-3, xyzw + m, .3135251434752570);
      m += gen_oh4(nc, .7268511789249627E-3, xyzw + m, .3721558339375338);
      m += gen_oh4(nc, .7422637534208629E-3, xyzw + m, .4286809575195696);
      m += gen_oh4(nc, .7509545035841214E-3, xyzw + m, .4822510128282994);
      m += gen_oh4(nc, .7548535057718401E-3, xyzw + m, .5320679333566263);
      m += gen_oh4(nc, .7554088969774001E-3, xyzw + m, .6172998195394274);
      m += gen_oh4(nc, .7553147174442808E-3, xyzw + m, .6510679849127481);
      m += gen_oh4(nc, .7564767653292297E-3, xyzw + m, .6777315251687360);
      m += gen_oh4(nc, .7587991808518730E-3, xyzw + m, .6963109410648741);
      m += gen_oh4(nc, .7608261832033027E-3, xyzw + m, .7058935009831749);
      m += gen_oh5(nc, .4021680447874916E-3, xyzw + m, .9955546194091857);
      m += gen_oh5(nc, .5804871793945964E-3, xyzw + m, .9734115901794209);
      m += gen_oh5(nc, .6792151955945159E-3, xyzw + m, .9275693732388626);
      m += gen_oh5(nc, .7336741211286294E-3, xyzw + m, .8568022422795103);
      m += gen_oh5(nc, .7581866300989608E-3, xyzw + m, .7623495553719372);
      m += gen_oh6(nc, .7538257859800743E-3, xyzw + m, .5707522908892223, .4387028039889501);
      m += gen_oh6(nc, .7483517247053123E-3, xyzw + m, .5196463388403083, .3858908414762617);
      m += gen_oh6(nc, .7371763661112059E-3, xyzw + m, .4646337531215351, .3301937372343854);
      m += gen_oh6(nc, .7183448895756934E-3, xyzw + m, .4063901697557691, .2725423573563777);
      m += gen_oh6(nc, .6895815529822191E-3, xyzw + m, .3456329466643087, .2139510237495250);
      m += gen_oh6(nc, .6480105801792886E-3, xyzw + m, .2831395121050332, .1555922309786647);
      m += gen_oh6(nc, .5897558896594636E-3, xyzw + m, .2197682022925330, .0989287897968610);
      m += gen_oh6(nc, .5095708849247346E-3, xyzw + m, .1564696098650355, .0459864291067551);
      m += gen_oh6(nc, .7536906428909755E-3, xyzw + m, .6027356673721295, .3376625140173426);
      m += gen_oh6(nc, .7472505965575118E-3, xyzw + m, .5496032320255096, .2822301309727988);
      m += gen_oh6(nc, .7343017132279698E-3, xyzw + m, .4921707755234567, .2248632342592540);
      m += gen_oh6(nc, .7130871582177445E-3, xyzw + m, .4309422998598483, .1666224723456479);
      m += gen_oh6(nc, .6817022032112776E-3, xyzw + m, .3664108182313672, .1086964901822169);
      m += gen_oh6(nc, .6380941145604121E-3, xyzw + m, .2990189057758436, .0525198978412008);
      m += gen_oh6(nc, .7550381377920310E-3, xyzw + m, .6268724013144998, .2297523657550023);
      m += gen_oh6(nc, .7478646640144802E-3, xyzw + m, .5707324144834607, .1723080607093800);
      m += gen_oh6(nc, .7335918720601220E-3, xyzw + m, .5096360901960365, .1140238465390513);
      m += gen_oh6(nc, .7110120527658118E-3, xyzw + m, .4438729938312456, .0561152209588254);
      m += gen_oh6(nc, .7571363978689501E-3, xyzw + m, .6419978471082389, .1164174423140873);
      m += gen_oh6(nc, .7489908329079233E-3, xyzw + m, .5817218061802611, .0579758953144522);
    break; case 1730:
      m += gen_oh1(nc, .6309049437420976E-4, xyzw + m);
      m += gen_oh2(nc, .6398287705571748E-3, xyzw + m);
      m += gen_oh3(nc, .6357185073530719E-3, xyzw + m);
      m += gen_oh4(nc, .2221207162188168E-3, xyzw + m, .0286092312619466);
      m += gen_oh4(nc, .3475784022286848E-3, xyzw + m, .0714255676771152);
      m += gen_oh4(nc, .4350742443589804E-3, xyzw + m, .1209199540995559);
      m += gen_oh4(nc, .4978569136522127E-3, xyzw + m, .1738673106594379);
      m += gen_oh4(nc, .5435036221998053E-3, xyzw + m, .2284645438467734);
      m += gen_oh4(nc, .5765913388219542E-3, xyzw + m, .2834807671701512);
      m += gen_oh4(nc, .6001200359226003E-3, xyzw + m, .3379680145467339);
      m += gen_oh4(nc, .6162178172717512E-3, xyzw + m, .3911355454819537);
      m += gen_oh4(nc, .6265218152438484E-3, xyzw + m, .4422860353001403);
      m += gen_oh4(nc, .6323987160974212E-3, xyzw + m, .4907781568726057);
      m += gen_oh4(nc, .6350767851540569E-3, xyzw + m, .5360006153211468);
      m += gen_oh4(nc, .6354362775297107E-3, xyzw + m, .6142105973596603);
      m += gen_oh4(nc, .6352302462706236E-3, xyzw + m, .6459300387977503);
      m += gen_oh4(nc, .6358117881417972E-3, xyzw + m, .6718056125089225);
      m += gen_oh4(nc, .6373101590310116E-3, xyzw + m, .6910888533186254);
      m += gen_oh4(nc, .6390428961368665E-3, xyzw + m, .7030467416823252);
      m += gen_oh5(nc, .3186913449946576E-3, xyzw + m, .0835495116635465);
      m += gen_oh5(nc, .4678028558591711E-3, xyzw + m, .2050143009099486);
      m += gen_oh5(nc, .5538829697598626E-3, xyzw + m, .3370208290706637);
      m += gen_oh5(nc, .6044475907190476E-3, xyzw + m, .4689051484233963);
      m += gen_oh5(nc, .6313575103509012E-3, xyzw + m, .5939400424557334);
      m += gen_oh6(nc, .4078626431855630E-3, xyzw + m, .1394983311832261, .0409758116205034);
      m += gen_oh6(nc, .4759933057812725E-3, xyzw + m, .1967999180485014, .0885198739129335);
      m += gen_oh6(nc, .5268151186413440E-3, xyzw + m, .2546183732548967, .1397680182969819);
      m += gen_oh6(nc, .5643048560507316E-3, xyzw + m, .3121281074713875, .1929452542226526);
      m += gen_oh6(nc, .5914501076613073E-3, xyzw + m, .3685981078502492, .2467898337061562);
      m += gen_oh6(nc, .6104561257874195E-3, xyzw + m, .4233760321547856, .3003104124785409);
      m += gen_oh6(nc, .6230252860707806E-3, xyzw + m, .4758671236059246, .3526684328175033);
      m += gen_oh6(nc, .6305618761760796E-3, xyzw + m, .5255178579796463, .4031134861145713);
      m += gen_oh6(nc, .6343092767597889E-3, xyzw + m, .5718025633734589, .4509426448342351);
      m += gen_oh6(nc, .5176268945737827E-3, xyzw + m, .2686927772723415, .0471132250242325);
      m += gen_oh6(nc, .5564840313313692E-3, xyzw + m, .3306006819904809, .0978448730394269);
      m += gen_oh6(nc, .5856426671038980E-3, xyzw + m, .3904906850594983, .1505395810025273);
      m += gen_oh6(nc, .6066386925777091E-3, xyzw + m, .4479957951904390, .2039728156296050);
      m += gen_oh6(nc, .6208824962234458E-3, xyzw + m, .5027076848919780, .2571529941121107);
      m += gen_oh6(nc, .6296314297822907E-3, xyzw + m, .5542087392260217, .3092191375815670);
      m += gen_oh6(nc, .6340423756791859E-3, xyzw + m, .6020850887375186, .3593807506130276);
      m += gen_oh6(nc, .5829627677107342E-3, xyzw + m, .4019851409179594, .0506338993437867);
      m += gen_oh6(nc, .6048693376081110E-3, xyzw + m, .4635614567449800, .1032422269160612);
      m += gen_oh6(nc, .6202362317732461E-3, xyzw + m, .5215860931591575, .1566322094006254);
      m += gen_oh6(nc, .6299005328403779E-3, xyzw + m, .5758202499099271, .2098082827491099);
      m += gen_oh6(nc, .6347722390609352E-3, xyzw + m, .6259893683876795, .2618824114553391);
      m += gen_oh6(nc, .6203778981238834E-3, xyzw + m, .5313795124811891, .0526324501933856);
      m += gen_oh6(nc, .6308414671239979E-3, xyzw + m, .5893317955931995, .1061059730982005);
      m += gen_oh6(nc, .6362706466959498E-3, xyzw + m, .6426246321215801, .1594171564034221);
      m += gen_oh6(nc, .6375414170333233E-3, xyzw + m, .6511904367376113, .0535478953656554);
    break; case 2030:
      m += gen_oh1(nc, .4656031899197431E-4, xyzw + m);
      m += gen_oh3(nc, .5421549195295507E-3, xyzw + m);
      m += gen_oh4(nc, .1778522133346553E-3, xyzw + m, .0254083533681435);
      m += gen_oh4(nc, .2811325405682796E-3, xyzw + m, .0639932280050492);
      m += gen_oh4(nc, .3548896312631459E-3, xyzw + m, .1088269469804125);
      m += gen_oh4(nc, .4090310897173364E-3, xyzw + m, .1570670798818287);
      m += gen_oh4(nc, .4493286134169965E-3, xyzw + m, .2071163932282514);
      m += gen_oh4(nc, .4793728447962723E-3, xyzw + m, .2578914044450844);
      m += gen_oh4(nc, .5015415319164265E-3, xyzw + m, .3085687558169623);
      m += gen_oh4(nc, .5175127372677937E-3, xyzw + m, .3584719706267024);
      m += gen_oh4(nc, .5285522262081019E-3, xyzw + m, .4070135594428709);
      m += gen_oh4(nc, .5356832703713962E-3, xyzw + m, .4536618626222638);
      m += gen_oh4(nc, .5397914736175170E-3, xyzw + m, .4979195686463577);
      m += gen_oh4(nc, .5416899441599930E-3, xyzw + m, .5393075111126999);
      m += gen_oh4(nc, .5419308476889938E-3, xyzw + m, .6115617676843916);
      m += gen_oh4(nc, .5416936902030596E-3, xyzw + m, .6414308435160159);
      m += gen_oh4(nc, .5419544338703164E-3, xyzw + m, .6664099412721607);
      m += gen_oh4(nc, .5428983656630974E-3, xyzw + m, .6859161771214913);
      m += gen_oh4(nc, .5442286500098193E-3, xyzw + m, .6993625593503890);
      m += gen_oh4(nc, .5452250345057301E-3, xyzw + m, .7062393387719380);
      m += gen_oh5(nc, .2568002497728530E-3, xyzw + m, .0747902816834976);
      m += gen_oh5(nc, .3827211700292145E-3, xyzw + m, .1848951153969366);
      m += gen_oh5(nc, .4579491561917824E-3, xyzw + m, .3059529066581305);
      m += gen_oh5(nc, .5042003969083574E-3, xyzw + m, .4285556101021362);
      m += gen_oh5(nc, .5312708889976024E-3, xyzw + m, .5468758653496526);
      m += gen_oh5(nc, .5438401790747117E-3, xyzw + m, .6565821978343439);
      m += gen_oh6(nc, .3316041873197344E-3, xyzw + m, .1253901572367117, .0368191722643964);
      m += gen_oh6(nc, .3899113567153771E-3, xyzw + m, .1775721510383941, .0798248760721330);
      m += gen_oh6(nc, .4343343327201309E-3, xyzw + m, .2305693358216114, .1264640966592335);
      m += gen_oh6(nc, .4679415262318919E-3, xyzw + m, .2836502845992063, .1751585683418957);
      m += gen_oh6(nc, .4930847981631031E-3, xyzw + m, .3361794746232590, .2247995907632670);
      m += gen_oh6(nc, .5115031867540091E-3, xyzw + m, .3875979172264824, .2745299257422246);
      m += gen_oh6(nc, .5245217148457367E-3, xyzw + m, .4374019316999074, .3236373482441118);
      m += gen_oh6(nc, .5332041499895321E-3, xyzw + m, .4851275843340022, .3714967859436741);
      m += gen_oh6(nc, .5384583126021542E-3, xyzw + m, .5303391803806868, .4175353646321745);
      m += gen_oh6(nc, .5411067210798852E-3, xyzw + m, .5726197380596287, .4612084406355461);
      m += gen_oh6(nc, .4259797391468714E-3, xyzw + m, .2431520732564863, .0425804013304395);
      m += gen_oh6(nc, .4604931368460021E-3, xyzw + m, .3002096800895869, .0886942430672272);
      m += gen_oh6(nc, .4871814878255202E-3, xyzw + m, .3558554457457432, .1368811706510655);
      m += gen_oh6(nc, .5072242910074885E-3, xyzw + m, .4097782537048887, .1860739985015033);
      m += gen_oh6(nc, .5217069845235350E-3, xyzw + m, .4616337666067458, .2354235077395853);
      m += gen_oh6(nc, .5315785966280310E-3, xyzw + m, .5110707008417874, .2842074921347011);
      m += gen_oh6(nc, .5376833708758905E-3, xyzw + m, .5577415286163795, .3317784414984102);
      m += gen_oh6(nc, .5408032092069521E-3, xyzw + m, .6013060431366950, .3775299002040700);
      m += gen_oh6(nc, .4842744917904866E-3, xyzw + m, .3661596767261781, .0459936788716459);
      m += gen_oh6(nc, .5048926076188130E-3, xyzw + m, .4237633153506581, .0940489377365442);
      m += gen_oh6(nc, .5202607980478373E-3, xyzw + m, .4786328454658452, .1431377109091971);
      m += gen_oh6(nc, .5309932388325743E-3, xyzw + m, .5305702076789774, .1924186388843570);
      m += gen_oh6(nc, .5377419770895208E-3, xyzw + m, .5793436224231788, .2411590944775190);
      m += gen_oh6(nc, .5411696331677717E-3, xyzw + m, .6247069017094747, .2886871491583605);
      m += gen_oh6(nc, .5197996293282420E-3, xyzw + m, .4874315552535204, .0480497877495321);
      m += gen_oh6(nc, .5311120836622945E-3, xyzw + m, .5427337322059053, .0971685719936666);
      m += gen_oh6(nc, .5384309319956951E-3, xyzw + m, .5943493747246700, .1465205839795055);
      m += gen_oh6(nc, .5421859504051886E-3, xyzw + m, .6421314033564943, .1953579449803574);
      m += gen_oh6(nc, .5390948355046314E-3, xyzw + m, .6020628374713980, .0491637501573811);
      m += gen_oh6(nc, .5433312705027845E-3, xyzw + m, .6529222529856881, .0986162154012701);
    break; case 2354:
      m += gen_oh1(nc, .3922616270665292E-4, xyzw + m);
      m += gen_oh2(nc, .4703831750854424E-3, xyzw + m);
      m += gen_oh3(nc, .4678202801282136E-3, xyzw + m);
      m += gen_oh4(nc, .1437832228979900E-3, xyzw + m, .0229002464653059);
      m += gen_oh4(nc, .2303572493577644E-3, xyzw + m, .0577908665227128);
      m += gen_oh4(nc, .2933110752447454E-3, xyzw + m, .0986310357637598);
      m += gen_oh4(nc, .3402905998359838E-3, xyzw + m, .1428155792982185);
      m += gen_oh4(nc, .3759138466870372E-3, xyzw + m, .1888978116601463);
      m += gen_oh4(nc, .4030638447899798E-3, xyzw + m, .2359091682970210);
      m += gen_oh4(nc, .4236591432242211E-3, xyzw + m, .2831228833706171);
      m += gen_oh4(nc, .4390522656946746E-3, xyzw + m, .3299495857966693);
      m += gen_oh4(nc, .4502523466626247E-3, xyzw + m, .3758840802660796);
      m += gen_oh4(nc, .4580577727783541E-3, xyzw + m, .4204751831009480);
      m += gen_oh4(nc, .4631391616615899E-3, xyzw + m, .4633068518751051);
      m += gen_oh4(nc, .4660928953698676E-3, xyzw + m, .5039849474507313);
      m += gen_oh4(nc, .4674751807936953E-3, xyzw + m, .5421265793440747);
      m += gen_oh4(nc, .4676414903932920E-3, xyzw + m, .6092660230557310);
      m += gen_oh4(nc, .4674086492347870E-3, xyzw + m, .6374654204984869);
      m += gen_oh4(nc, .4674928539483207E-3, xyzw + m, .6615136472609892);
      m += gen_oh4(nc, .4680748979686447E-3, xyzw + m, .6809487285958127);
      m += gen_oh4(nc, .4690449806389040E-3, xyzw + m, .6952980021665196);
      m += gen_oh4(nc, .4699877075860818E-3, xyzw + m, .7041245497695400);
      m += gen_oh5(nc, .2099942281069176E-3, xyzw + m, .0674403308830606);
      m += gen_oh5(nc, .3172269150712804E-3, xyzw + m, .1678684485334166);
      m += gen_oh5(nc, .3832051358546523E-3, xyzw + m, .2793559049539613);
      m += gen_oh5(nc, .4252193818146985E-3, xyzw + m, .3935264218057639);
      m += gen_oh5(nc, .4513807963755000E-3, xyzw + m, .5052629268232558);
      m += gen_oh5(nc, .4657797469114178E-3, xyzw + m, .6107905315437531);
      m += gen_oh6(nc, .2733362800522836E-3, xyzw + m, .1135081039843524, .0333195488466259);
      m += gen_oh6(nc, .3235485368463559E-3, xyzw + m, .1612866626099378, .0724716746543654);
      m += gen_oh6(nc, .3624908726013453E-3, xyzw + m, .2100786550168205, .1151539110849745);
      m += gen_oh6(nc, .3925540070712828E-3, xyzw + m, .2592282009459942, .1599491097143677);
      m += gen_oh6(nc, .4156129781116235E-3, xyzw + m, .3081740561320203, .2058699956028027);
      m += gen_oh6(nc, .4330644984623263E-3, xyzw + m, .3564289781578164, .2521624953502911);
      m += gen_oh6(nc, .4459677725921312E-3, xyzw + m, .4035587288240703, .2982090785797674);
      m += gen_oh6(nc, .4551593004456795E-3, xyzw + m, .4491671196373903, .3434762087235733);
      m += gen_oh6(nc, .4613341462749918E-3, xyzw + m, .4928854782917489, .3874831357203437);
      m += gen_oh6(nc, .4651019618269806E-3, xyzw + m, .5343646791958988, .4297814821746926);
      m += gen_oh6(nc, .4670249536100625E-3, xyzw + m, .5732683216530990, .4699402260943537);
      m += gen_oh6(nc, .3549555576441708E-3, xyzw + m, .2214131583218986, .0387360204064389);
      m += gen_oh6(nc, .3856108245249010E-3, xyzw + m, .2741796504750071, .0808949625690201);
      m += gen_oh6(nc, .4098622845756882E-3, xyzw + m, .3259797439149485, .1251732177620872);
      m += gen_oh6(nc, .4286328604268950E-3, xyzw + m, .3765441148826891, .1706260286403185);
      m += gen_oh6(nc, .4427802198993945E-3, xyzw + m, .4255773574530558, .2165115147300408);
      m += gen_oh6(nc, .4530473511488561E-3, xyzw + m, .4727795117058430, .2622089812225259);
      m += gen_oh6(nc, .4600805475703138E-3, xyzw + m, .5178546895819012, .3071721431296201);
      m += gen_oh6(nc, .4644599059958017E-3, xyzw + m, .5605141192097460, .3508998998801138);
      m += gen_oh6(nc, .4667274455712508E-3, xyzw + m, .6004763319352512, .3929160876166931);
      m += gen_oh6(nc, .4069360518020356E-3, xyzw + m, .3352842634946949, .0420256345728802);
      m += gen_oh6(nc, .4260442819919195E-3, xyzw + m, .3891971629814670, .0861430975887085);
      m += gen_oh6(nc, .4408678508029063E-3, xyzw + m, .4409875565542281, .1314500879380001);
      m += gen_oh6(nc, .4518748115548597E-3, xyzw + m, .4904893058592484, .1772189657383859);
      m += gen_oh6(nc, .4595564875375116E-3, xyzw + m, .5375056138769549, .2228277110050294);
      m += gen_oh6(nc, .4643988774315846E-3, xyzw + m, .5818255708669969, .2677179935014386);
      m += gen_oh6(nc, .4668827491646946E-3, xyzw + m, .6232334858144959, .3113675035544165);
      m += gen_oh6(nc, .4400541823741973E-3, xyzw + m, .4489485354492058, .0440916237836817);
      m += gen_oh6(nc, .4514512890193797E-3, xyzw + m, .5015136875933150, .0893900991774849);
      m += gen_oh6(nc, .4596198627347549E-3, xyzw + m, .5511300550512623, .1351806029383365);
      m += gen_oh6(nc, .4648659016801781E-3, xyzw + m, .5976720409858000, .1808370355053196);
      m += gen_oh6(nc, .4675502017157673E-3, xyzw + m, .6409956378989354, .2257852192301602);
      m += gen_oh6(nc, .4598494476455523E-3, xyzw + m, .5581222330827514, .0453217342163716);
      m += gen_oh6(nc, .4654916955152048E-3, xyzw + m, .6074705984161695, .0911748803184031);
      m += gen_oh6(nc, .4684709779505137E-3, xyzw + m, .6532272537379032, .1369294213140155);
      m += gen_oh6(nc, .4691445539106986E-3, xyzw + m, .6594761494500487, .0458990148727558);
    break; case 2702:
      m += gen_oh1(nc, .2998675149888161E-4, xyzw + m);
      m += gen_oh3(nc, .4077860529495355E-3, xyzw + m);
      m += gen_oh4(nc, .1185349192520667E-3, xyzw + m, .0206556253881870);
      m += gen_oh4(nc, .1913408643425751E-3, xyzw + m, .0525091817302238);
      m += gen_oh4(nc, .2452886577209897E-3, xyzw + m, .0899348008203838);
      m += gen_oh4(nc, .2862408183288702E-3, xyzw + m, .1306023924436019);
      m += gen_oh4(nc, .3178032258257357E-3, xyzw + m, .1732060388531418);
      m += gen_oh4(nc, .3422945667633690E-3, xyzw + m, .2168727084820249);
      m += gen_oh4(nc, .3612790520235922E-3, xyzw + m, .2609528309173586);
      m += gen_oh4(nc, .3758638229818521E-3, xyzw + m, .3049252927938952);
      m += gen_oh4(nc, .3868711798859953E-3, xyzw + m, .3483484138084404);
      m += gen_oh4(nc, .3949429933189938E-3, xyzw + m, .3908321549106406);
      m += gen_oh4(nc, .4006068107541156E-3, xyzw + m, .4320210071894814);
      m += gen_oh4(nc, .4043192149672723E-3, xyzw + m, .4715824795890053);
      m += gen_oh4(nc, .4064947495808078E-3, xyzw + m, .5091984794078454);
      m += gen_oh4(nc, .4075245619813152E-3, xyzw + m, .5445580145650804);
      m += gen_oh4(nc, .4076423540893566E-3, xyzw + m, .6072575796841768);
      m += gen_oh4(nc, .4074280862251555E-3, xyzw + m, .6339484505755802);
      m += gen_oh4(nc, .4074163756012244E-3, xyzw + m, .6570718257486958);
      m += gen_oh4(nc, .4077647795071246E-3, xyzw + m, .6762557330090709);
      m += gen_oh4(nc, .4084517552782530E-3, xyzw + m, .6911161696923790);
      m += gen_oh4(nc, .4092468459224052E-3, xyzw + m, .7012841911659961);
      m += gen_oh4(nc, .4097872687240906E-3, xyzw + m, .7064559272410020);
      m += gen_oh5(nc, .1738986811745028E-3, xyzw + m, .0612355498989477);
      m += gen_oh5(nc, .2659616045280191E-3, xyzw + m, .1533070348312393);
      m += gen_oh5(nc, .3240596008171533E-3, xyzw + m, .2563902605244206);
      m += gen_oh5(nc, .3621195964432943E-3, xyzw + m, .3629346991663361);
      m += gen_oh5(nc, .3868838330760539E-3, xyzw + m, .4683949968987538);
      m += gen_oh5(nc, .4018911532693111E-3, xyzw + m, .5694479240657953);
      m += gen_oh5(nc, .4089929432983252E-3, xyzw + m, .6634465430993955);
      m += gen_oh6(nc, .2279907527706409E-3, xyzw + m, .1033958573552305, .0303454400906358);
      m += gen_oh6(nc, .2715205490578897E-3, xyzw + m, .1473521412414395, .0661880304424713);
      m += gen_oh6(nc, .3057917896703976E-3, xyzw + m, .1924552158705967, .1054431128987715);
      m += gen_oh6(nc, .3326913052452555E-3, xyzw + m, .2381094362890328, .1468263551238858);
      m += gen_oh6(nc, .3537334711890037E-3, xyzw + m, .2838121707936760, .1894486108187886);
      m += gen_oh6(nc, .3700567500783129E-3, xyzw + m, .3291323133373415, .2326374238761579);
      m += gen_oh6(nc, .3825245372589122E-3, xyzw + m, .3736896978741460, .2758485808485768);
      m += gen_oh6(nc, .3918125171518296E-3, xyzw + m, .4171406040760013, .3186179331996921);
      m += gen_oh6(nc, .3984720419937579E-3, xyzw + m, .4591677985256915, .3605329796303794);
      m += gen_oh6(nc, .4029746003338211E-3, xyzw + m, .4994733831718418, .4012147253586509);
      m += gen_oh6(nc, .4057428632156627E-3, xyzw + m, .5377731830445096, .4403050025570692);
      m += gen_oh6(nc, .4071719274114857E-3, xyzw + m, .5737917830001331, .4774565904277483);
      m += gen_oh6(nc, .2990236950664119E-3, xyzw + m, .2027323586271389, .0354412250497615);
      m += gen_oh6(nc, .3262951734212878E-3, xyzw + m, .2516942375187273, .0741830438864633);
      m += gen_oh6(nc, .3482634608242413E-3, xyzw + m, .3000227995257181, .1150502745727186);
      m += gen_oh6(nc, .3656596681700892E-3, xyzw + m, .3474806691046342, .1571963371209364);
      m += gen_oh6(nc, .3791740467794218E-3, xyzw + m, .3938103180359209, .1999631877247100);
      m += gen_oh6(nc, .3894034450156905E-3, xyzw + m, .4387519590455703, .2428073457846535);
      m += gen_oh6(nc, .3968600245508371E-3, xyzw + m, .4820503960077787, .2852575132906155);
      m += gen_oh6(nc, .4019931351420050E-3, xyzw + m, .5234573778475101, .3268884208674639);
      m += gen_oh6(nc, .4052108801278599E-3, xyzw + m, .5627318647235282, .3673033321675939);
      m += gen_oh6(nc, .4068978613940934E-3, xyzw + m, .5996390607156954, .4061211551830290);
      m += gen_oh6(nc, .3454275351319704E-3, xyzw + m, .3084780753791947, .0386012552310006);
      m += gen_oh6(nc, .3629963537007920E-3, xyzw + m, .3589988275920223, .0792893898710487);
      m += gen_oh6(nc, .3770187233889873E-3, xyzw + m, .4078628415881973, .1212614643030087);
      m += gen_oh6(nc, .3878608613694378E-3, xyzw + m, .4549287258889735, .1638770827382693);
      m += gen_oh6(nc, .3959065270221274E-3, xyzw + m, .5000278512957279, .2065965798260176);
      m += gen_oh6(nc, .4015286975463570E-3, xyzw + m, .5429785044928199, .2489436378852235);
      m += gen_oh6(nc, .4050866785614717E-3, xyzw + m, .5835939850491711, .2904811368946891);
      m += gen_oh6(nc, .4069320185051913E-3, xyzw + m, .6216870353444856, .3307941957666609);
      m += gen_oh6(nc, .3760120964062763E-3, xyzw + m, .4151104662709091, .0406482914605255);
      m += gen_oh6(nc, .3870969564418064E-3, xyzw + m, .4649804275009218, .0825842454729476);
      m += gen_oh6(nc, .3955287790534055E-3, xyzw + m, .5124695757009662, .1251841962027289);
      m += gen_oh6(nc, .4015361911302668E-3, xyzw + m, .5574711100606224, .1679107505976331);
      m += gen_oh6(nc, .4053836986719548E-3, xyzw + m, .5998597333287227, .2102805057358715);
      m += gen_oh6(nc, .4073578673299117E-3, xyzw + m, .6395007148516600, .2518418087774107);
      m += gen_oh6(nc, .3954628379231406E-3, xyzw + m, .5188456224746252, .0419432167607752);
      m += gen_oh6(nc, .4017645508847530E-3, xyzw + m, .5664190707942778, .0845766155192150);
      m += gen_oh6(nc, .4059030348651293E-3, xyzw + m, .6110464353283153, .1273652932519396);
      m += gen_oh6(nc, .4080565809484880E-3, xyzw + m, .6526430302051563, .1698173239076354);
      m += gen_oh6(nc, .4063018753664651E-3, xyzw + m, .6167551880377548, .0426639885154886);
      m += gen_oh6(nc, .4087191292799671E-3, xyzw + m, .6607195418355383, .0855192581423835);
    break; case 3074:
      m += gen_oh1(nc, .2599095953754734E-4, xyzw + m);
      m += gen_oh2(nc, .3603134089687541E-3, xyzw + m);
      m += gen_oh3(nc, .3586067974412447E-3, xyzw + m);
      m += gen_oh4(nc, .9831528474385880E-4, xyzw + m, .0188610851872339);
      m += gen_oh4(nc, .1605023107954450E-3, xyzw + m, .0480021724462530);
      m += gen_oh4(nc, .2072200131464099E-3, xyzw + m, .0824492205839724);
      m += gen_oh4(nc, .2431297618814187E-3, xyzw + m, .1200408362484023);
      m += gen_oh4(nc, .2711819064496707E-3, xyzw + m, .1595773530809965);
      m += gen_oh4(nc, .2932762038321116E-3, xyzw + m, .2002635973434064);
      m += gen_oh4(nc, .3107032514197368E-3, xyzw + m, .2415127590139982);
      m += gen_oh4(nc, .3243808058921213E-3, xyzw + m, .2828584158458477);
      m += gen_oh4(nc, .3349899091374030E-3, xyzw + m, .3239091015338138);
      m += gen_oh4(nc, .3430580688505218E-3, xyzw + m, .3643225097962194);
      m += gen_oh4(nc, .3490124109290343E-3, xyzw + m, .4037897083691802);
      m += gen_oh4(nc, .3532148948561955E-3, xyzw + m, .4420247515194127);
      m += gen_oh4(nc, .3559862669062833E-3, xyzw + m, .4787572538464938);
      m += gen_oh4(nc, .3576224317551411E-3, xyzw + m, .5137265251275234);
      m += gen_oh4(nc, .3584050533086076E-3, xyzw + m, .5466764056654611);
      m += gen_oh4(nc, .3584903581373224E-3, xyzw + m, .6054859420813535);
      m += gen_oh4(nc, .3582991879040586E-3, xyzw + m, .6308106701764562);
      m += gen_oh4(nc, .3582371187963125E-3, xyzw + m, .6530369230179583);
      m += gen_oh4(nc, .3584353631122350E-3, xyzw + m, .6718609524611158);
      m += gen_oh4(nc, .3589120166517785E-3, xyzw + m, .6869676499894013);
      m += gen_oh4(nc, .3595445704531601E-3, xyzw + m, .6980467077240748);
      m += gen_oh4(nc, .3600943557111074E-3, xyzw + m, .7048241721250522);
      m += gen_oh5(nc, .1456447096742039E-3, xyzw + m, .0559110522205823);
      m += gen_oh5(nc, .2252370188283782E-3, xyzw + m, .1407384078513916);
      m += gen_oh5(nc, .2766135443474897E-3, xyzw + m, .2364035438976309);
      m += gen_oh5(nc, .3110729491500851E-3, xyzw + m, .3360602737818170);
      m += gen_oh5(nc, .3342506712303391E-3, xyzw + m, .4356292630054665);
      m += gen_oh5(nc, .3491981834026860E-3, xyzw + m, .5321569415256174);
      m += gen_oh5(nc, .3576003604348932E-3, xyzw + m, .6232956305040555);
      m += gen_oh6(nc, .1921921305788564E-3, xyzw + m, .0946987008683847, .0277874838730947);
      m += gen_oh6(nc, .2301458216495632E-3, xyzw + m, .1353170300568141, .0607656987862836);
      m += gen_oh6(nc, .2604248549522893E-3, xyzw + m, .1771679481726077, .0970307276271104);
      m += gen_oh6(nc, .2845275425870697E-3, xyzw + m, .2197066664231751, .1354112458524762);
      m += gen_oh6(nc, .3036870897974840E-3, xyzw + m, .2624783557374927, .1750996479744100);
      m += gen_oh6(nc, .3188414832298066E-3, xyzw + m, .3050969521214442, .2154896907449802);
      m += gen_oh6(nc, .3307046414722089E-3, xyzw + m, .3472252637196021, .2560954625740152);
      m += gen_oh6(nc, .3398330969031360E-3, xyzw + m, .3885610219026360, .2965070050624096);
      m += gen_oh6(nc, .3466757899705373E-3, xyzw + m, .4288273776062765, .3363641488734497);
      m += gen_oh6(nc, .3516095923230054E-3, xyzw + m, .4677662471302948, .3753400029836788);
      m += gen_oh6(nc, .3549645184048486E-3, xyzw + m, .5051333589553360, .4131297522144286);
      m += gen_oh6(nc, .3570415969441392E-3, xyzw + m, .5406942145810492, .4494423776081795);
      m += gen_oh6(nc, .3581251798496118E-3, xyzw + m, .5742204122576458, .4839938958841502);
      m += gen_oh6(nc, .2543491329913348E-3, xyzw + m, .1865407027225188, .0325914485107080);
      m += gen_oh6(nc, .2786711051330776E-3, xyzw + m, .2321186453689432, .0683567950529734);
      m += gen_oh6(nc, .2985552361083679E-3, xyzw + m, .2773159142523882, .1062284864451989);
      m += gen_oh6(nc, .3145867929154039E-3, xyzw + m, .3219200192237254, .1454404409323047);
      m += gen_oh6(nc, .3273290662067609E-3, xyzw + m, .3657032593944029, .1854018282582510);
      m += gen_oh6(nc, .3372705511943501E-3, xyzw + m, .4084376778363622, .2256297412014750);
      m += gen_oh6(nc, .3448274437851510E-3, xyzw + m, .4499004945751427, .2657104425000896);
      m += gen_oh6(nc, .3503592783048583E-3, xyzw + m, .4898758141326335, .3052755487631557);
      m += gen_oh6(nc, .3541854792663162E-3, xyzw + m, .5281547442266309, .3439863920645423);
      m += gen_oh6(nc, .3565995517909428E-3, xyzw + m, .5645346989813992, .3815229456121914);
      m += gen_oh6(nc, .3578802078302898E-3, xyzw + m, .5988181252159848, .4175752420966734);
      m += gen_oh6(nc, .2958644592860982E-3, xyzw + m, .2850425424471603, .0356214950986254);
      m += gen_oh6(nc, .3119548129116835E-3, xyzw + m, .3324619433027876, .0733031888687110);
      m += gen_oh6(nc, .3250745225005984E-3, xyzw + m, .3785848333076282, .1123226296008472);
      m += gen_oh6(nc, .3355153415935208E-3, xyzw + m, .4232891028562115, .1521084193337708);
      m += gen_oh6(nc, .3435847568549328E-3, xyzw + m, .4664287050829722, .1921844459223610);
      m += gen_oh6(nc, .3495786831622488E-3, xyzw + m, .5078458493735726, .2321360989678303);
      m += gen_oh6(nc, .3537767805534621E-3, xyzw + m, .5473779816204180, .2715886486360520);
      m += gen_oh6(nc, .3564459815421428E-3, xyzw + m, .5848617133811376, .3101924707571355);
      m += gen_oh6(nc, .3578464061225468E-3, xyzw + m, .6201348281584887, .3476121052890973);
      m += gen_oh6(nc, .3239748762836212E-3, xyzw + m, .3852191185387871, .0376322488003511);
      m += gen_oh6(nc, .3345491784174287E-3, xyzw + m, .4325025061073423, .0765958193563713);
      m += gen_oh6(nc, .3429126177301782E-3, xyzw + m, .4778486229734490, .1163381306083900);
      m += gen_oh6(nc, .3492420343097421E-3, xyzw + m, .5211663693009000, .1563890598752899);
      m += gen_oh6(nc, .3537399050235257E-3, xyzw + m, .5623469504853703, .1963320810149200);
      m += gen_oh6(nc, .3566209152659172E-3, xyzw + m, .6012718188659246, .2357847407258738);
      m += gen_oh6(nc, .3581084321919782E-3, xyzw + m, .6378179206390117, .2743846121244060);
      m += gen_oh6(nc, .3426522117591512E-3, xyzw + m, .4836936460214534, .0389590261073902);
      m += gen_oh6(nc, .3491848770121379E-3, xyzw + m, .5293792562683797, .0787124681931264);
      m += gen_oh6(nc, .3539318235231476E-3, xyzw + m, .5726281253100033, .1187963808202981);
      m += gen_oh6(nc, .3570231438458694E-3, xyzw + m, .6133658776169068, .1587914708061787);
      m += gen_oh6(nc, .3586207335051714E-3, xyzw + m, .6515085491865307, .1983058575227646);
      m += gen_oh6(nc, .3541196205164025E-3, xyzw + m, .5778692716064976, .0397720968979154);
      m += gen_oh6(nc, .3574296911573953E-3, xyzw + m, .6207904288086192, .0799015759298115);
      m += gen_oh6(nc, .3591993279818963E-3, xyzw + m, .6608688171046802, .1199671308754309);
      m += gen_oh6(nc, .3595855034661997E-3, xyzw + m, .6656263089489129, .0401595595780597);
    break; case 3470:
      m += gen_oh1(nc, .2040382730826330E-4, xyzw + m);
      m += gen_oh3(nc, .3178149703889544E-3, xyzw + m);
      m += gen_oh4(nc, .8288115128076111E-4, xyzw + m, .0172142083290623);
      m += gen_oh4(nc, .1360883192522954E-3, xyzw + m, .0440887537498177);
      m += gen_oh4(nc, .1766854454542662E-3, xyzw + m, .0759468081387868);
      m += gen_oh4(nc, .2083153161230153E-3, xyzw + m, .1108335359204799);
      m += gen_oh4(nc, .2333279544657158E-3, xyzw + m, .1476517054388567);
      m += gen_oh4(nc, .2532809539930247E-3, xyzw + m, .1856731870860615);
      m += gen_oh4(nc, .2692472184211158E-3, xyzw + m, .2243634099428821);
      m += gen_oh4(nc, .2819949946811885E-3, xyzw + m, .2633006881662727);
      m += gen_oh4(nc, .2920953593973030E-3, xyzw + m, .3021340904916283);
      m += gen_oh4(nc, .2999889782948352E-3, xyzw + m, .3405594048030089);
      m += gen_oh4(nc, .3060292120496902E-3, xyzw + m, .3783044434007372);
      m += gen_oh4(nc, .3105109167522192E-3, xyzw + m, .4151194767407910);
      m += gen_oh4(nc, .3136902387550312E-3, xyzw + m, .4507705766443257);
      m += gen_oh4(nc, .3157984652454632E-3, xyzw + m, .4850346056573187);
      m += gen_oh4(nc, .3170516518425422E-3, xyzw + m, .5176950817792469);
      m += gen_oh4(nc, .3176568425633755E-3, xyzw + m, .5485384240820989);
      m += gen_oh4(nc, .3177198411207062E-3, xyzw + m, .6039117238943308);
      m += gen_oh4(nc, .3175519492394733E-3, xyzw + m, .6279956655573113);
      m += gen_oh4(nc, .3174654952634756E-3, xyzw + m, .6493636169568952);
      m += gen_oh4(nc, .3175676415467654E-3, xyzw + m, .6677644117704504);
      m += gen_oh4(nc, .3178923417835410E-3, xyzw + m, .6829368572115624);
      m += gen_oh4(nc, .3183788287531909E-3, xyzw + m, .6946195818184121);
      m += gen_oh4(nc, .3188755151918807E-3, xyzw + m, .7025711542057026);
      m += gen_oh4(nc, .3191916889313849E-3, xyzw + m, .7066004767140119);
      m += gen_oh5(nc, .1231779611744508E-3, xyzw + m, .0513253768994606);
      m += gen_oh5(nc, .1924661373839880E-3, xyzw + m, .1297994661331225);
      m += gen_oh5(nc, .2380881867403424E-3, xyzw + m, .2188852049401307);
      m += gen_oh5(nc, .2693100663037885E-3, xyzw + m, .3123174824903457);
      m += gen_oh5(nc, .2908673382834366E-3, xyzw + m, .4064037620738195);
      m += gen_oh5(nc, .3053914619381535E-3, xyzw + m, .4984958396944782);
      m += gen_oh5(nc, .3143916684147777E-3, xyzw + m, .5864975046021365);
      m += gen_oh5(nc, .3187042244055363E-3, xyzw + m, .6686711634580175);
      m += gen_oh6(nc, .1635219535869790E-3, xyzw + m, .0871573878083595, .0255717523336758);
      m += gen_oh6(nc, .1968109917696070E-3, xyzw + m, .1248383123134007, .0560482338337668);
      m += gen_oh6(nc, .2236754342249974E-3, xyzw + m, .1638062693383378, .0896856860190076);
      m += gen_oh6(nc, .2453186687017181E-3, xyzw + m, .2035586203373176, .1254086651976279);
      m += gen_oh6(nc, .2627551791580541E-3, xyzw + m, .2436798975293774, .1624780150162012);
      m += gen_oh6(nc, .2767654860152220E-3, xyzw + m, .2838207507773806, .2003422342683208);
      m += gen_oh6(nc, .2879467027765895E-3, xyzw + m, .3236787502217692, .2385628026255263);
      m += gen_oh6(nc, .2967639918918702E-3, xyzw + m, .3629849554840691, .2767731148783578);
      m += gen_oh6(nc, .3035900684660351E-3, xyzw + m, .4014948081992087, .3146542308245309);
      m += gen_oh6(nc, .3087338237298308E-3, xyzw + m, .4389818379260225, .3519196415895088);
      m += gen_oh6(nc, .3124608838860167E-3, xyzw + m, .4752331143674377, .3883050984023654);
      m += gen_oh6(nc, .3150084294226743E-3, xyzw + m, .5100457318374018, .4235613423908649);
      m += gen_oh6(nc, .3165958398598402E-3, xyzw + m, .5432238388954868, .4574484717196220);
      m += gen_oh6(nc, .3174320440957372E-3, xyzw + m, .5745758685072442, .4897311639255524);
      m += gen_oh6(nc, .2182188909812599E-3, xyzw + m, .1723981437592809, .0301063059788110);
      m += gen_oh6(nc, .2399727933921445E-3, xyzw + m, .2149553257844597, .0632603155420469);
      m += gen_oh6(nc, .2579796133514652E-3, xyzw + m, .2573256081247422, .0984856698025863);
      m += gen_oh6(nc, .2727114052623535E-3, xyzw + m, .2993163751238106, .1350835952384266);
      m += gen_oh6(nc, .2846327656281355E-3, xyzw + m, .3407238005148000, .1725184055442181);
      m += gen_oh6(nc, .2941491102051334E-3, xyzw + m, .3813454978483264, .2103559279730725);
      m += gen_oh6(nc, .3016049492136107E-3, xyzw + m, .4209848104423343, .2482278774554860);
      m += gen_oh6(nc, .3072949726175648E-3, xyzw + m, .4594519699996300, .2858099509982883);
      m += gen_oh6(nc, .3114768142886460E-3, xyzw + m, .4965640166185930, .3228075659915428);
      m += gen_oh6(nc, .3143823673666223E-3, xyzw + m, .5321441655571562, .3589459907204151);
      m += gen_oh6(nc, .3162269764661535E-3, xyzw + m, .5660208438582166, .3939630088864310);
      m += gen_oh6(nc, .3172164663759821E-3, xyzw + m, .5980264315964364, .4276029922949089);
      m += gen_oh6(nc, .2554575398967435E-3, xyzw + m, .2644215852350733, .0330093942907255);
      m += gen_oh6(nc, .2701704069135677E-3, xyzw + m, .3090113743443063, .0680388765007850);
      m += gen_oh6(nc, .2823693413468940E-3, xyzw + m, .3525871079197808, .1044326136206709);
      m += gen_oh6(nc, .2922898463214289E-3, xyzw + m, .3950418005354029, .1416751597517679);
      m += gen_oh6(nc, .3001829062162428E-3, xyzw + m, .4362475663430163, .1793408610504821);
      m += gen_oh6(nc, .3062890864542953E-3, xyzw + m, .4760661812145854, .2170630750175722);
      m += gen_oh6(nc, .3108328279264746E-3, xyzw + m, .5143551042512103, .2545145157815807);
      m += gen_oh6(nc, .3140243146201245E-3, xyzw + m, .5509709026935597, .2913940101706601);
      m += gen_oh6(nc, .3160638030977130E-3, xyzw + m, .5857711030329428, .3274169910910705);
      m += gen_oh6(nc, .3171462882206275E-3, xyzw + m, .6186149917404392, .3623081329317265);
      m += gen_oh6(nc, .2812388416031796E-3, xyzw + m, .3586894569557064, .0349735438645004);
      m += gen_oh6(nc, .2912137500288045E-3, xyzw + m, .4035266610019441, .0712973673975709);
      m += gen_oh6(nc, .2993241256502206E-3, xyzw + m, .4467775312332510, .1084758620193165);
      m += gen_oh6(nc, .3057101738983822E-3, xyzw + m, .4883638346608543, .1460915689241772);
      m += gen_oh6(nc, .3105319326251432E-3, xyzw + m, .5281908348434601, .1837790832369980);
      m += gen_oh6(nc, .3139565514428167E-3, xyzw + m, .5661542687149311, .2212075390874021);
      m += gen_oh6(nc, .3161543006806366E-3, xyzw + m, .6021450102031451, .2580682841160985);
      m += gen_oh6(nc, .3172985960613294E-3, xyzw + m, .6360520783610050, .2940656362094121);
      m += gen_oh6(nc, .2989400336901431E-3, xyzw + m, .4521611065087196, .0363105536586700);
      m += gen_oh6(nc, .3054555883947677E-3, xyzw + m, .4959365651560963, .0734831846848435);
      m += gen_oh6(nc, .3104764960807702E-3, xyzw + m, .5376815804038283, .1111087643812648);
      m += gen_oh6(nc, .3141015825977616E-3, xyzw + m, .5773314480243767, .1488226085145408);
      m += gen_oh6(nc, .3164520621159896E-3, xyzw + m, .6148113245575056, .1862892274135151);
      m += gen_oh6(nc, .3176652305912204E-3, xyzw + m, .6500407462842380, .2231909701714456);
      m += gen_oh6(nc, .3105097161023939E-3, xyzw + m, .5425151448707213, .0371820130611894);
      m += gen_oh6(nc, .3143014117890550E-3, xyzw + m, .5841860556907931, .0748361633506735);
      m += gen_oh6(nc, .3168172866287200E-3, xyzw + m, .6234632186851500, .1125990834266120);
      m += gen_oh6(nc, .3181401865570968E-3, xyzw + m, .6602934551848842, .1501303813157619);
      m += gen_oh6(nc, .3170663659156037E-3, xyzw + m, .6278573968375105, .0376755993024572);
      m += gen_oh6(nc, .3185447944625510E-3, xyzw + m, .6665611711264577, .0754844330136016);
    break; case 3890:
      m += gen_oh1(nc, .1807395252196920E-4, xyzw + m);
      m += gen_oh2(nc, .2848008782238827E-3, xyzw + m);
      m += gen_oh3(nc, .2836065837530581E-3, xyzw + m);
      m += gen_oh4(nc, .7013149266673816E-4, xyzw + m, .0158787641985835);
      m += gen_oh4(nc, .1162798021956766E-3, xyzw + m, .0406919359375121);
      m += gen_oh4(nc, .1518728583972105E-3, xyzw + m, .0702588811525800);
      m += gen_oh4(nc, .1798796108216934E-3, xyzw + m, .1027495450028704);
      m += gen_oh4(nc, .2022593385972785E-3, xyzw + m, .1371457730893426);
      m += gen_oh4(nc, .2203093105575464E-3, xyzw + m, .1727758532671953);
      m += gen_oh4(nc, .2349294234299855E-3, xyzw + m, .2091492038929037);
      m += gen_oh4(nc, .2467682058747003E-3, xyzw + m, .2458813281751915);
      m += gen_oh4(nc, .2563092683572224E-3, xyzw + m, .2826545859450066);
      m += gen_oh4(nc, .2639253896763318E-3, xyzw + m, .3191957291799622);
      m += gen_oh4(nc, .2699137479265108E-3, xyzw + m, .3552621469299578);
      m += gen_oh4(nc, .2745196420166739E-3, xyzw + m, .3906329503406230);
      m += gen_oh4(nc, .2779529197397593E-3, xyzw + m, .4251028614093031);
      m += gen_oh4(nc, .2803996086684265E-3, xyzw + m, .4584777520111870);
      m += gen_oh4(nc, .2820302356715842E-3, xyzw + m, .4905711358710193);
      m += gen_oh4(nc, .2830056747491068E-3, xyzw + m, .5212011669847385);
      m += gen_oh4(nc, .2834808950776839E-3, xyzw + m, .5501878488737995);
      m += gen_oh4(nc, .2835282339078929E-3, xyzw + m, .6025037877479342);
      m += gen_oh4(nc, .2833819267065800E-3, xyzw + m, .6254572689549016);
      m += gen_oh4(nc, .2832858336906784E-3, xyzw + m, .6460107179528248);
      m += gen_oh4(nc, .2833268235451244E-3, xyzw + m, .6639541138154251);
      m += gen_oh4(nc, .2835432677029253E-3, xyzw + m, .6790688515667495);
      m += gen_oh4(nc, .2839091722743049E-3, xyzw + m, .6911338580371512);
      m += gen_oh4(nc, .2843308178875841E-3, xyzw + m, .6999385956126490);
      m += gen_oh4(nc, .2846703550533846E-3, xyzw + m, .7053037748656896);
      m += gen_oh5(nc, .1051193406971900E-3, xyzw + m, .0473222438718012);
      m += gen_oh5(nc, .1657871838796974E-3, xyzw + m, .1202100529326803);
      m += gen_oh5(nc, .2064648113714232E-3, xyzw + m, .2034304820664855);
      m += gen_oh5(nc, .2347942745819741E-3, xyzw + m, .2912285643573002);
      m += gen_oh5(nc, .2547775326597726E-3, xyzw + m, .3802361792726768);
      m += gen_oh5(nc, .2686876684847025E-3, xyzw + m, .4680598511056146);
      m += gen_oh5(nc, .2778665755515867E-3, xyzw + m, .5528151052155599);
      m += gen_oh5(nc, .2830996616782929E-3, xyzw + m, .6329386307803041);
      m += gen_oh6(nc, .1403063340168372E-3, xyzw + m, .0805651665136907, .0236345468400312);
      m += gen_oh6(nc, .1696504125939477E-3, xyzw + m, .1156476077139389, .0519129163254594);
      m += gen_oh6(nc, .1935787242745390E-3, xyzw + m, .1520473382760421, .0832271573699452);
      m += gen_oh6(nc, .2130614510521968E-3, xyzw + m, .1892986699745931, .1165855667993712);
      m += gen_oh6(nc, .2289381265931048E-3, xyzw + m, .2270194446777792, .1513077167409504);
      m += gen_oh6(nc, .2418630292816186E-3, xyzw + m, .2648908185093273, .1868882025807859);
      m += gen_oh6(nc, .2523400495631193E-3, xyzw + m, .3026389259574136, .2229277629776224);
      m += gen_oh6(nc, .2607623973449605E-3, xyzw + m, .3400220296151384, .2590951840746235);
      m += gen_oh6(nc, .2674441032689209E-3, xyzw + m, .3768217953335510, .2951047291750847);
      m += gen_oh6(nc, .2726432360343356E-3, xyzw + m, .4128372900921884, .3307019714169930);
      m += gen_oh6(nc, .2765787685924545E-3, xyzw + m, .4478807131815630, .3656544101087634);
      m += gen_oh6(nc, .2794428690642224E-3, xyzw + m, .4817742034089257, .3997448951939695);
      m += gen_oh6(nc, .2814099002062895E-3, xyzw + m, .5143472814653344, .4327667110812024);
      m += gen_oh6(nc, .2826429531578994E-3, xyzw + m, .5454346213905650, .4645196123532293);
      m += gen_oh6(nc, .2832983542550884E-3, xyzw + m, .5748739313170252, .4948063555703345);
      m += gen_oh6(nc, .1886695565284976E-3, xyzw + m, .1599598738286342, .0279235759004898);
      m += gen_oh6(nc, .2081867882748234E-3, xyzw + m, .1998097412500951, .0587714103813907);
      m += gen_oh6(nc, .2245148680600796E-3, xyzw + m, .2396228952566202, .0916457391469138);
      m += gen_oh6(nc, .2380370491511872E-3, xyzw + m, .2792228341097746, .1259049641962687);
      m += gen_oh6(nc, .2491398041852455E-3, xyzw + m, .3184251107546741, .1610594823400863);
      m += gen_oh6(nc, .2581632405881230E-3, xyzw + m, .3570481164426244, .1967151653460898);
      m += gen_oh6(nc, .2653965506227417E-3, xyzw + m, .3949164710492144, .2325404606175168);
      m += gen_oh6(nc, .2710857216747087E-3, xyzw + m, .4318617293970503, .2682461141151439);
      m += gen_oh6(nc, .2754434093903659E-3, xyzw + m, .4677221009931678, .3035720116011973);
      m += gen_oh6(nc, .2786579932519380E-3, xyzw + m, .5023417939270955, .3382781859197439);
      m += gen_oh6(nc, .2809011080679474E-3, xyzw + m, .5355701836636128, .3721383065625942);
      m += gen_oh6(nc, .2823336184560987E-3, xyzw + m, .5672608451328771, .4049346360466055);
      m += gen_oh6(nc, .2831101175806309E-3, xyzw + m, .5972704202540162, .4364538098633802);
      m += gen_oh6(nc, .2221679970354546E-3, xyzw + m, .2461687022333596, .0307042316683337);
      m += gen_oh6(nc, .2356185734270703E-3, xyzw + m, .2881774566286831, .0633803466928188);
      m += gen_oh6(nc, .2469228344805590E-3, xyzw + m, .3293963604116978, .0974286248706794);
      m += gen_oh6(nc, .2562726348642046E-3, xyzw + m, .3697303822241377, .1323799532282290);
      m += gen_oh6(nc, .2638756726753028E-3, xyzw + m, .4090663023135127, .1678497018129336);
      m += gen_oh6(nc, .2699311157390862E-3, xyzw + m, .4472819355411712, .2035095105326114);
      m += gen_oh6(nc, .2746233268403837E-3, xyzw + m, .4842513377231437, .2390692566672091);
      m += gen_oh6(nc, .2781225674454771E-3, xyzw + m, .5198477629962928, .2742649818076149);
      m += gen_oh6(nc, .2805881254045684E-3, xyzw + m, .5539453011883145, .3088503806580094);
      m += gen_oh6(nc, .2821719877004913E-3, xyzw + m, .5864196762401251, .3425904245906614);
      m += gen_oh6(nc, .2830222502333124E-3, xyzw + m, .6171484466668390, .3752562294789468);
      m += gen_oh6(nc, .2457995956744870E-3, xyzw + m, .3350337830565727, .0326158993463475);
      m += gen_oh6(nc, .2551474407503706E-3, xyzw + m, .3775773224758284, .0665843892808157);
      m += gen_oh6(nc, .2629065335195311E-3, xyzw + m, .4188155229848973, .1014565797157954);
      m += gen_oh6(nc, .2691900449925075E-3, xyzw + m, .4586805892009344, .1368573320843822);
      m += gen_oh6(nc, .2741275485754276E-3, xyzw + m, .4970895714224235, .1724614851951608);
      m += gen_oh6(nc, .2778530970122595E-3, xyzw + m, .5339505133960747, .2079779381416412);
      m += gen_oh6(nc, .2805010567646741E-3, xyzw + m, .5691665792531440, .2431385788322288);
      m += gen_oh6(nc, .2822055834031040E-3, xyzw + m, .6026387682680377, .2776901883049853);
      m += gen_oh6(nc, .2831016901243473E-3, xyzw + m, .6342676150163307, .3113881356386632);
      m += gen_oh6(nc, .2624474901131803E-3, xyzw + m, .4237951119537067, .0339487784866435);
      m += gen_oh6(nc, .2688034163039377E-3, xyzw + m, .4656918683234929, .0688021955629145);
      m += gen_oh6(nc, .2738932751287636E-3, xyzw + m, .5058857069185980, .1041946859721635);
      m += gen_oh6(nc, .2777944791242523E-3, xyzw + m, .5443204666713995, .1398039738736393);
      m += gen_oh6(nc, .2806011661660987E-3, xyzw + m, .5809298813759742, .1753373381196155);
      m += gen_oh6(nc, .2824181456597460E-3, xyzw + m, .6156416039447128, .2105215793514010);
      m += gen_oh6(nc, .2833585216577828E-3, xyzw + m, .6483801351066604, .2450953312157051);
      m += gen_oh6(nc, .2738165236962878E-3, xyzw + m, .5103616577251688, .0348556064380072);
      m += gen_oh6(nc, .2778365208203180E-3, xyzw + m, .5506738792580681, .0702630863151203);
      m += gen_oh6(nc, .2807852940418966E-3, xyzw + m, .5889573040995292, .1059035061296403);
      m += gen_oh6(nc, .2827245949674705E-3, xyzw + m, .6251641589516930, .1414823925236026);
      m += gen_oh6(nc, .2837342344829828E-3, xyzw + m, .6592414921570178, .1767207908214530);
      m += gen_oh6(nc, .2809233907610981E-3, xyzw + m, .5930314017533383, .0354218933956167);
      m += gen_oh6(nc, .2829930809742694E-3, xyzw + m, .6309812253390175, .0710957404036955);
      m += gen_oh6(nc, .2841097874111479E-3, xyzw + m, .6666296011353230, .1067259792282730);
      m += gen_oh6(nc, .2843455206008783E-3, xyzw + m, .6703715271049921, .0356945526882081);
    break; case 4334:
      m += gen_oh1(nc, .1449063022537883E-4, xyzw + m);
      m += gen_oh3(nc, .2546377329828424E-3, xyzw + m);
      m += gen_oh4(nc, .6018432961087496E-4, xyzw + m, .0146289615183101);
      m += gen_oh4(nc, .1002286583263673E-3, xyzw + m, .0376984081249314);
      m += gen_oh4(nc, .1315222931028093E-3, xyzw + m, .0652470190409689);
      m += gen_oh4(nc, .1564213746876724E-3, xyzw + m, .0956054341613465);
      m += gen_oh4(nc, .1765118841507736E-3, xyzw + m, .1278335898929198);
      m += gen_oh4(nc, .1928737099311080E-3, xyzw + m, .1613096104466031);
      m += gen_oh4(nc, .2062658534263270E-3, xyzw + m, .1955806225745371);
      m += gen_oh4(nc, .2172395445953787E-3, xyzw + m, .2302935218498028);
      m += gen_oh4(nc, .2262076188876047E-3, xyzw + m, .2651584344113027);
      m += gen_oh4(nc, .2334885699462397E-3, xyzw + m, .2999276825183209);
      m += gen_oh4(nc, .2393355273179203E-3, xyzw + m, .3343828669718798);
      m += gen_oh4(nc, .2439559200468863E-3, xyzw + m, .3683265013750518);
      m += gen_oh4(nc, .2475251866060002E-3, xyzw + m, .4015763206518108);
      m += gen_oh4(nc, .2501965558158773E-3, xyzw + m, .4339612026399770);
      m += gen_oh4(nc, .2521081407925925E-3, xyzw + m, .4653180651114582);
      m += gen_oh4(nc, .2533881002388081E-3, xyzw + m, .4954893331080803);
      m += gen_oh4(nc, .2541582900848261E-3, xyzw + m, .5243207068924930);
      m += gen_oh4(nc, .2545365737525860E-3, xyzw + m, .5516590479041704);
      m += gen_oh4(nc, .2545726993066799E-3, xyzw + m, .6012371927804177);
      m += gen_oh4(nc, .2544456197465555E-3, xyzw + m, .6231574466449818);
      m += gen_oh4(nc, .2543481596881064E-3, xyzw + m, .6429416514181271);
      m += gen_oh4(nc, .2543506451429194E-3, xyzw + m, .6604124272943594);
      m += gen_oh4(nc, .2544905675493763E-3, xyzw + m, .6753851470408250);
      m += gen_oh4(nc, .2547611407344429E-3, xyzw + m, .6876717970626161);
      m += gen_oh4(nc, .2551060375448869E-3, xyzw + m, .6970895061319234);
      m += gen_oh4(nc, .2554291933816039E-3, xyzw + m, .7034746912553310);
      m += gen_oh4(nc, .2556255710686343E-3, xyzw + m, .7067017217542295);
      m += gen_oh5(nc, .9041339695118196E-4, xyzw + m, .0438222350113112);
      m += gen_oh5(nc, .1438426330079022E-3, xyzw + m, .1117474077400006);
      m += gen_oh5(nc, .1802523089820518E-3, xyzw + m, .1897153252911440);
      m += gen_oh5(nc, .2060052290565496E-3, xyzw + m, .2724023009910331);
      m += gen_oh5(nc, .2245002248967466E-3, xyzw + m, .3567163308709902);
      m += gen_oh5(nc, .2377059847731150E-3, xyzw + m, .4404784483028087);
      m += gen_oh5(nc, .2468118955882525E-3, xyzw + m, .5219833154161411);
      m += gen_oh5(nc, .2525410872966528E-3, xyzw + m, .5998179868977553);
      m += gen_oh5(nc, .2553101409933397E-3, xyzw + m, .6727803154548222);
      m += gen_oh6(nc, .1212879733668632E-3, xyzw + m, .0747656394316609, .0219316850946118);
      m += gen_oh6(nc, .1472872881270931E-3, xyzw + m, .1075341482001416, .0482641928153389);
      m += gen_oh6(nc, .1686846601010828E-3, xyzw + m, .1416344885203259, .0775119188357574);
      m += gen_oh6(nc, .1862698414660208E-3, xyzw + m, .1766325315388586, .1087558139247680);
      m += gen_oh6(nc, .2007430956991861E-3, xyzw + m, .2121744174481514, .1413661374253096);
      m += gen_oh6(nc, .2126568125394796E-3, xyzw + m, .2479669443408145, .1748768214258880);
      m += gen_oh6(nc, .2224394603372113E-3, xyzw + m, .2837600452294113, .2089216406612073);
      m += gen_oh6(nc, .2304264522673135E-3, xyzw + m, .3193344933193984, .2431987685545972);
      m += gen_oh6(nc, .2368854288424087E-3, xyzw + m, .3544935442438745, .2774497054377770);
      m += gen_oh6(nc, .2420352089461772E-3, xyzw + m, .3890571932288154, .3114460356156915);
      m += gen_oh6(nc, .2460597113081295E-3, xyzw + m, .4228581214259090, .3449806851913012);
      m += gen_oh6(nc, .2491181912257687E-3, xyzw + m, .4557387211304052, .3778618641248256);
      m += gen_oh6(nc, .2513528194205857E-3, xyzw + m, .4875487950541643, .4099086391698978);
      m += gen_oh6(nc, .2528943096693220E-3, xyzw + m, .5181436529962997, .4409474925853973);
      m += gen_oh6(nc, .2538660368488136E-3, xyzw + m, .5473824095600661, .4708094517711291);
      m += gen_oh6(nc, .2543868648299022E-3, xyzw + m, .5751263398976174, .4993275140354637);
      m += gen_oh6(nc, .1642595537825183E-3, xyzw + m, .1489515746840028, .0259938199326702);
      m += gen_oh6(nc, .1818246659849308E-3, xyzw + m, .1863656444351767, .0547928653246219);
      m += gen_oh6(nc, .1966565649492420E-3, xyzw + m, .2238602880356348, .0855676325142525);
      m += gen_oh6(nc, .2090677905657991E-3, xyzw + m, .2612723375728160, .1177257802267011);
      m += gen_oh6(nc, .2193820409510504E-3, xyzw + m, .2984332990206190, .1508168456192700);
      m += gen_oh6(nc, .2278870827661928E-3, xyzw + m, .3351786584663333, .1844801892177727);
      m += gen_oh6(nc, .2348283192282090E-3, xyzw + m, .3713505522209120, .2184145236087598);
      m += gen_oh6(nc, .2404139755581477E-3, xyzw + m, .4067981098954663, .2523590641486229);
      m += gen_oh6(nc, .2448227407760734E-3, xyzw + m, .4413769993687534, .2860812976901373);
      m += gen_oh6(nc, .2482110455592573E-3, xyzw + m, .4749487182516394, .3193686757808996);
      m += gen_oh6(nc, .2507192397774103E-3, xyzw + m, .5073798105075426, .3520226949547602);
      m += gen_oh6(nc, .2524765968534880E-3, xyzw + m, .5385410448878654, .3838544395667890);
      m += gen_oh6(nc, .2536052388539425E-3, xyzw + m, .5683065353670530, .4146810037640963);
      m += gen_oh6(nc, .2542230588033068E-3, xyzw + m, .5965527620663510, .4443224094681121);
      m += gen_oh6(nc, .1944817013047896E-3, xyzw + m, .2299227700856157, .0286575766405758);
      m += gen_oh6(nc, .2067862362746635E-3, xyzw + m, .2695752998553267, .0592342168448599);
      m += gen_oh6(nc, .2172440734649114E-3, xyzw + m, .3086178716611389, .0911781777605772);
      m += gen_oh6(nc, .2260125991723423E-3, xyzw + m, .3469649871659077, .1240593814082605);
      m += gen_oh6(nc, .2332655008689523E-3, xyzw + m, .3845153566319655, .1575272058259175);
      m += gen_oh6(nc, .2391699681532458E-3, xyzw + m, .4211600033403215, .1912845163525413);
      m += gen_oh6(nc, .2438801528273928E-3, xyzw + m, .4567867834329882, .2250710177858171);
      m += gen_oh6(nc, .2475370504260665E-3, xyzw + m, .4912829319232061, .2586521303440910);
      m += gen_oh6(nc, .2502707235640574E-3, xyzw + m, .5245364793303812, .2918112242865407);
      m += gen_oh6(nc, .2522031701054241E-3, xyzw + m, .5564369788915756, .3243439239067890);
      m += gen_oh6(nc, .2534511269978784E-3, xyzw + m, .5868757697775288, .3560536787835351);
      m += gen_oh6(nc, .2541284914955151E-3, xyzw + m, .6157458853519617, .3867480821242581);
      m += gen_oh6(nc, .2161509250688394E-3, xyzw + m, .3138461110672113, .0305137463750728);
      m += gen_oh6(nc, .2248778513437852E-3, xyzw + m, .3542495872050569, .0623711123373075);
      m += gen_oh6(nc, .2322388803404617E-3, xyzw + m, .3935751553120181, .0951622395240191);
      m += gen_oh6(nc, .2383265471001355E-3, xyzw + m, .4317634668111147, .1285467341508517);
      m += gen_oh6(nc, .2432476675019525E-3, xyzw + m, .4687413842250821, .1622318931656033);
      m += gen_oh6(nc, .2471122223750674E-3, xyzw + m, .5044274237060283, .1959581153836453);
      m += gen_oh6(nc, .2500291752486870E-3, xyzw + m, .5387354077925727, .2294888081183837);
      m += gen_oh6(nc, .2521055942764682E-3, xyzw + m, .5715768898356105, .2626031152713945);
      m += gen_oh6(nc, .2534472785575503E-3, xyzw + m, .6028627200136111, .2950904075286713);
      m += gen_oh6(nc, .2541599713080121E-3, xyzw + m, .6325039812653463, .3267458451113286);
      m += gen_oh6(nc, .2317380975862936E-3, xyzw + m, .3981986708423407, .0318329145874982);
      m += gen_oh6(nc, .2378550733719775E-3, xyzw + m, .4382791182133300, .0645954819388091);
      m += gen_oh6(nc, .2428884456739118E-3, xyzw + m, .4769233057218166, .0979575703708795);
      m += gen_oh6(nc, .2469002655757292E-3, xyzw + m, .5140823911194238, .1316307235126655);
      m += gen_oh6(nc, .2499657574265851E-3, xyzw + m, .5496977833862983, .1653556486358704);
      m += gen_oh6(nc, .2521676168486082E-3, xyzw + m, .5837047306512727, .1988931724126510);
      m += gen_oh6(nc, .2535935662645334E-3, xyzw + m, .6160349566926879, .2320174581438950);
      m += gen_oh6(nc, .2543356743363214E-3, xyzw + m, .6466185353209440, .2645106562168662);
      m += gen_oh6(nc, .2427353285201535E-3, xyzw + m, .4810835158795404, .0327591780774399);
      m += gen_oh6(nc, .2468258039744386E-3, xyzw + m, .5199925041324341, .0661254618396718);
      m += gen_oh6(nc, .2500060956440310E-3, xyzw + m, .5571717692207494, .0998149833147414);
      m += gen_oh6(nc, .2523238365420979E-3, xyzw + m, .5925789250836379, .1335687001410374);
      m += gen_oh6(nc, .2538399260252846E-3, xyzw + m, .6261658523859670, .1671444402896463);
      m += gen_oh6(nc, .2546255927268069E-3, xyzw + m, .6578811126669331, .2003106382156076);
      m += gen_oh6(nc, .2500583360048449E-3, xyzw + m, .5609624612998100, .0333750094023134);
      m += gen_oh6(nc, .2524777638260203E-3, xyzw + m, .5979959659984670, .0670875033590180);
      m += gen_oh6(nc, .2540951193860656E-3, xyzw + m, .6330523711054002, .1008792126424850);
      m += gen_oh6(nc, .2549524085027472E-3, xyzw + m, .6660960998103972, .1345050343171794);
      m += gen_oh6(nc, .2542569507009158E-3, xyzw + m, .6365384364585819, .0337279946073705);
      m += gen_oh6(nc, .2552114127580376E-3, xyzw + m, .6710994302899275, .0675524930967803);
    break; case 4802:
      m += gen_oh1(nc, .9687521879420705E-4, xyzw + m);
      m += gen_oh2(nc, .2307897895367918E-3, xyzw + m);
      m += gen_oh3(nc, .2297310852498558E-3, xyzw + m);
      m += gen_oh4(nc, .7386265944001918E-4, xyzw + m, .0233572860888706);
      m += gen_oh4(nc, .8257977698542210E-4, xyzw + m, .0435298783655065);
      m += gen_oh4(nc, .9706044762057630E-4, xyzw + m, .0643920052108880);
      m += gen_oh4(nc, .1302393847117003E-3, xyzw + m, .0900394363199318);
      m += gen_oh4(nc, .1541957004600968E-3, xyzw + m, .1196706615548473);
      m += gen_oh4(nc, .1704459770092199E-3, xyzw + m, .1511715412838134);
      m += gen_oh4(nc, .1827374890942906E-3, xyzw + m, .1835982828503801);
      m += gen_oh4(nc, .1926360817436107E-3, xyzw + m, .2165081259155405);
      m += gen_oh4(nc, .2008010239494833E-3, xyzw + m, .2496208720417563);
      m += gen_oh4(nc, .2075635983209175E-3, xyzw + m, .2827200673567900);
      m += gen_oh4(nc, .2131306638690909E-3, xyzw + m, .3156190823994346);
      m += gen_oh4(nc, .2176562329937335E-3, xyzw + m, .3481476793749115);
      m += gen_oh4(nc, .2212682262991018E-3, xyzw + m, .3801466086947226);
      m += gen_oh4(nc, .2240799515668565E-3, xyzw + m, .4114652119634011);
      m += gen_oh4(nc, .2261959816187525E-3, xyzw + m, .4419598786519751);
      m += gen_oh4(nc, .2277156368808855E-3, xyzw + m, .4714925949329543);
      m += gen_oh4(nc, .2287351772128336E-3, xyzw + m, .4999293972879466);
      m += gen_oh4(nc, .2293490814084085E-3, xyzw + m, .5271387221431248);
      m += gen_oh4(nc, .2296505312376273E-3, xyzw + m, .5529896780837761);
      m += gen_oh4(nc, .2296793832318756E-3, xyzw + m, .6000856099481712);
      m += gen_oh4(nc, .2295785443842974E-3, xyzw + m, .6210562192785175);
      m += gen_oh4(nc, .2295017931529102E-3, xyzw + m, .6401165879934240);
      m += gen_oh4(nc, .2295059638184868E-3, xyzw + m, .6571144029244333);
      m += gen_oh4(nc, .2296232343237362E-3, xyzw + m, .6718910821718863);
      m += gen_oh4(nc, .2298530178740771E-3, xyzw + m, .6842845591099010);
      m += gen_oh4(nc, .2301579790280501E-3, xyzw + m, .6941353476269816);
      m += gen_oh4(nc, .2304690404996513E-3, xyzw + m, .7012965242212991);
      m += gen_oh4(nc, .2307027995907102E-3, xyzw + m, .7056471428242644);
      m += gen_oh5(nc, .9312274696671092E-4, xyzw + m, .0459555764358590);
      m += gen_oh5(nc, .1199919385876926E-3, xyzw + m, .1049316742435023);
      m += gen_oh5(nc, .1598039138877690E-3, xyzw + m, .1773548879549274);
      m += gen_oh5(nc, .1822253763574900E-3, xyzw + m, .2559071411236127);
      m += gen_oh5(nc, .1988579593655040E-3, xyzw + m, .3358156837985898);
      m += gen_oh5(nc, .2112620102533307E-3, xyzw + m, .4155835743763893);
      m += gen_oh5(nc, .2201594887699007E-3, xyzw + m, .4937894296167472);
      m += gen_oh5(nc, .2261622590895036E-3, xyzw + m, .5691569694793316);
      m += gen_oh5(nc, .2296458453435705E-3, xyzw + m, .6405840854894251);
      m += gen_oh6(nc, .1006006990267000E-3, xyzw + m, .0734513389414335, .0217784408148607);
      m += gen_oh6(nc, .1227676689635876E-3, xyzw + m, .1009859834044931, .0459036218577519);
      m += gen_oh6(nc, .1467864280270117E-3, xyzw + m, .1324289619748758, .0725506309569088);
      m += gen_oh6(nc, .1644178912101232E-3, xyzw + m, .1654272109607127, .1017825451960684);
      m += gen_oh6(nc, .1777664890718961E-3, xyzw + m, .1990767186776461, .1325652320980364);
      m += gen_oh6(nc, .1884825664516690E-3, xyzw + m, .2330125945523278, .1642765374496765);
      m += gen_oh6(nc, .1973269246453848E-3, xyzw + m, .2670080611108287, .1965360374337889);
      m += gen_oh6(nc, .2046767775855328E-3, xyzw + m, .3008753376294316, .2290726770542238);
      m += gen_oh6(nc, .2107600125918040E-3, xyzw + m, .3344475596167860, .2616645495370823);
      m += gen_oh6(nc, .2157416362266829E-3, xyzw + m, .3675709724070786, .2941150728843141);
      m += gen_oh6(nc, .2197557816920721E-3, xyzw + m, .4001000887587812, .3262440400919066);
      m += gen_oh6(nc, .2229192611835437E-3, xyzw + m, .4318956350436028, .3578835350611916);
      m += gen_oh6(nc, .2253385110212775E-3, xyzw + m, .4628239056795531, .3888751854043678);
      m += gen_oh6(nc, .2271137107548774E-3, xyzw + m, .4927563229773636, .4190678003222840);
      m += gen_oh6(nc, .2283414092917525E-3, xyzw + m, .5215687136707969, .4483151836883852);
      m += gen_oh6(nc, .2291161673130077E-3, xyzw + m, .5491402346984905, .4764740676087880);
      m += gen_oh6(nc, .2295313908576598E-3, xyzw + m, .5753520160126075, .5034021310998277);
      m += gen_oh6(nc, .1438204721359031E-3, xyzw + m, .1388326356417754, .0243543651037281);
      m += gen_oh6(nc, .1607738025495257E-3, xyzw + m, .1743686900537244, .0511889705734265);
      m += gen_oh6(nc, .1741483853528379E-3, xyzw + m, .2099737037950268, .0801469504853963);
      m += gen_oh6(nc, .1851918467519151E-3, xyzw + m, .2454492590908548, .1105117874155699);
      m += gen_oh6(nc, .1944628638070613E-3, xyzw + m, .2807219257864278, .1417950531570966);
      m += gen_oh6(nc, .2022495446275152E-3, xyzw + m, .3156842271975842, .1736604945719597);
      m += gen_oh6(nc, .2087462382438514E-3, xyzw + m, .3502090945177752, .2058466324693981);
      m += gen_oh6(nc, .2141074754818308E-3, xyzw + m, .3841684849519686, .2381284261195919);
      m += gen_oh6(nc, .2184640913748162E-3, xyzw + m, .4174372367906016, .2703031270422569);
      m += gen_oh6(nc, .2219309165220329E-3, xyzw + m, .4498926465011892, .3021845683091309);
      m += gen_oh6(nc, .2246123118340624E-3, xyzw + m, .4814146229807701, .3335993355165720);
      m += gen_oh6(nc, .2266062766915125E-3, xyzw + m, .5118863625734701, .3643833735518232);
      m += gen_oh6(nc, .2280072952230796E-3, xyzw + m, .5411947455119144, .3943789541958179);
      m += gen_oh6(nc, .2289082025202583E-3, xyzw + m, .5692301500357246, .4234320144403542);
      m += gen_oh6(nc, .2294012695120025E-3, xyzw + m, .5958857204139576, .4513897947419260);
      m += gen_oh6(nc, .1722434488736947E-3, xyzw + m, .2156270284785766, .0268122575544449);
      m += gen_oh6(nc, .1830237421455091E-3, xyzw + m, .2532385054909710, .0555749574780561);
      m += gen_oh6(nc, .1923855349997633E-3, xyzw + m, .2902564617771537, .0856936806295025);
      m += gen_oh6(nc, .2004067861936271E-3, xyzw + m, .3266979823143256, .1167367450324135);
      m += gen_oh6(nc, .2071817297354263E-3, xyzw + m, .3625039627493614, .1483861994003304);
      m += gen_oh6(nc, .2128250834102103E-3, xyzw + m, .3975838937548699, .1803821503011405);
      m += gen_oh6(nc, .2174513719440102E-3, xyzw + m, .4318396099009774, .2124962965666424);
      m += gen_oh6(nc, .2211661839150214E-3, xyzw + m, .4651706555732742, .2445221837805913);
      m += gen_oh6(nc, .2240665257813102E-3, xyzw + m, .4974752649620969, .2762701224322987);
      m += gen_oh6(nc, .2262439516632620E-3, xyzw + m, .5286517579627517, .3075627775211328);
      m += gen_oh6(nc, .2277874557231869E-3, xyzw + m, .5586001195731894, .3382311089826877);
      m += gen_oh6(nc, .2287854314454994E-3, xyzw + m, .5872229902021319, .3681108834741399);
      m += gen_oh6(nc, .2293268499615575E-3, xyzw + m, .6144258616235123, .3970397446872839);
      m += gen_oh6(nc, .1912628201529828E-3, xyzw + m, .2951676508064861, .0286749953875044);
      m += gen_oh6(nc, .1992499672238701E-3, xyzw + m, .3335085485472725, .0586787934190351);
      m += gen_oh6(nc, .2061275533454027E-3, xyzw + m, .3709561760636381, .0896109920502228);
      m += gen_oh6(nc, .2119318215968572E-3, xyzw + m, .4074722861667498, .1211627927626297);
      m += gen_oh6(nc, .2167416581882652E-3, xyzw + m, .4429923648839117, .1530748903554898);
      m += gen_oh6(nc, .2206430730516600E-3, xyzw + m, .4774428052721736, .1851176436721877);
      m += gen_oh6(nc, .2237186938699523E-3, xyzw + m, .5107446539535904, .2170829107658179);
      m += gen_oh6(nc, .2260480075032884E-3, xyzw + m, .5428151370542935, .2487786689026271);
      m += gen_oh6(nc, .2277098884558542E-3, xyzw + m, .5735699292556964, .2800239952795016);
      m += gen_oh6(nc, .2287845715109671E-3, xyzw + m, .6029253794562865, .3106445702878119);
      m += gen_oh6(nc, .2293547268236294E-3, xyzw + m, .6307998987073145, .3404689500841194);
      m += gen_oh6(nc, .2056073839852528E-3, xyzw + m, .3752652273692719, .0299714509818448);
      m += gen_oh6(nc, .2114235865831876E-3, xyzw + m, .4135383879344028, .0608672589867801);
      m += gen_oh6(nc, .2163175629770551E-3, xyzw + m, .4506113885153907, .0923884954843564);
      m += gen_oh6(nc, .2203392158111650E-3, xyzw + m, .4864401554606072, .1242786603851851);
      m += gen_oh6(nc, .2235473176847839E-3, xyzw + m, .5209708076611709, .1563086731483386);
      m += gen_oh6(nc, .2260024141501235E-3, xyzw + m, .5541422135830122, .1882696509388506);
      m += gen_oh6(nc, .2277675929329182E-3, xyzw + m, .5858880915113817, .2199672979126059);
      m += gen_oh6(nc, .2289102112284834E-3, xyzw + m, .6161399390603444, .2512165482924867);
      m += gen_oh6(nc, .2295027954625118E-3, xyzw + m, .6448296482255090, .2818368701871888);
      m += gen_oh6(nc, .2161281589879992E-3, xyzw + m, .4544796274917948, .0308897040506031);
      m += gen_oh6(nc, .2201980477395102E-3, xyzw + m, .4919389072146628, .0624094767763684);
      m += gen_oh6(nc, .2234952066593166E-3, xyzw + m, .5279313026985183, .0943070614428031);
      m += gen_oh6(nc, .2260540098520838E-3, xyzw + m, .5624169925571135, .1263547818770374);
      m += gen_oh6(nc, .2279157981899988E-3, xyzw + m, .5953484627093287, .1583430788822594);
      m += gen_oh6(nc, .2291296918565571E-3, xyzw + m, .6266730715339185, .1900748462555988);
      m += gen_oh6(nc, .2297533752536649E-3, xyzw + m, .6563363204278871, .2213599519592567);
      m += gen_oh6(nc, .2234927356465995E-3, xyzw + m, .5314574716585696, .0315250881151537);
      m += gen_oh6(nc, .2261288012985219E-3, xyzw + m, .5674614932298185, .0634386529146556);
      m += gen_oh6(nc, .2280818160923688E-3, xyzw + m, .6017706004970264, .0955150350422395);
      m += gen_oh6(nc, .2293773295180159E-3, xyzw + m, .6343471270264178, .1275440099801196);
      m += gen_oh6(nc, .2300528767338634E-3, xyzw + m, .6651494599127802, .1593252037671960);
      m += gen_oh6(nc, .2281893855065666E-3, xyzw + m, .6050184986005704, .0319253833849611);
      m += gen_oh6(nc, .2295720444840727E-3, xyzw + m, .6390163550880400, .0640282435396231);
      m += gen_oh6(nc, .2303227649026753E-3, xyzw + m, .6711199107088448, .0960980507700291);
      m += gen_oh6(nc, .2304831913227114E-3, xyzw + m, .6741354429572275, .0321185319627323);
    break; case 5294:
      m += gen_oh1(nc, .9080510764308163E-4, xyzw + m);
      m += gen_oh3(nc, .2084824361987793E-3, xyzw + m);
      m += gen_oh4(nc, .5011105657239616E-4, xyzw + m, .0230326168626145);
      m += gen_oh4(nc, .5942520409683854E-4, xyzw + m, .0375720862016239);
      m += gen_oh4(nc, .9564394826109721E-4, xyzw + m, .0582191203382185);
      m += gen_oh4(nc, .1185530657126338E-3, xyzw + m, .0840312752919487);
      m += gen_oh4(nc, .1364510114230331E-3, xyzw + m, .1122927798060578);
      m += gen_oh4(nc, .1505828825605415E-3, xyzw + m, .1420125319192987);
      m += gen_oh4(nc, .1619298749867023E-3, xyzw + m, .1726396437341978);
      m += gen_oh4(nc, .1712450504267789E-3, xyzw + m, .2038170058115696);
      m += gen_oh4(nc, .1789891098164999E-3, xyzw + m, .2352849892876508);
      m += gen_oh4(nc, .1854474955629795E-3, xyzw + m, .2668363354312461);
      m += gen_oh4(nc, .1908148636673661E-3, xyzw + m, .2982941279900452);
      m += gen_oh4(nc, .1952377405281833E-3, xyzw + m, .3295002922087076);
      m += gen_oh4(nc, .1988349254282232E-3, xyzw + m, .3603094918363593);
      m += gen_oh4(nc, .2017079807160050E-3, xyzw + m, .3905857895173920);
      m += gen_oh4(nc, .2039473082709094E-3, xyzw + m, .4202005758160837);
      m += gen_oh4(nc, .2056360279288953E-3, xyzw + m, .4490310061597227);
      m += gen_oh4(nc, .2068525823066865E-3, xyzw + m, .4769586160311491);
      m += gen_oh4(nc, .2076724877534488E-3, xyzw + m, .5038679887049750);
      m += gen_oh4(nc, .2081694278237885E-3, xyzw + m, .5296454286519962);
      m += gen_oh4(nc, .2084157631219326E-3, xyzw + m, .5541776207164850);
      m += gen_oh4(nc, .2084381531128593E-3, xyzw + m, .5990467321921213);
      m += gen_oh4(nc, .2083476277129307E-3, xyzw + m, .6191467096294587);
      m += gen_oh4(nc, .2082686194459732E-3, xyzw + m, .6375251212901849);
      m += gen_oh4(nc, .2082475686112415E-3, xyzw + m, .6540514381131168);
      m += gen_oh4(nc, .2083139860289915E-3, xyzw + m, .6685899064391509);
      m += gen_oh4(nc, .2084745561831237E-3, xyzw + m, .6810013009681648);
      m += gen_oh4(nc, .2087091313375890E-3, xyzw + m, .6911469578730340);
      m += gen_oh4(nc, .2089718413297697E-3, xyzw + m, .6988956915141736);
      m += gen_oh4(nc, .2092003303479793E-3, xyzw + m, .7041335794868721);
      m += gen_oh4(nc, .2093336148263241E-3, xyzw + m, .7067754398018568);
      m += gen_oh5(nc, .7591708117365266E-4, xyzw + m, .0384036870785362);
      m += gen_oh5(nc, .1083383968169186E-3, xyzw + m, .0983548595411740);
      m += gen_oh5(nc, .1403019395292510E-3, xyzw + m, .1665774947612998);
      m += gen_oh5(nc, .1615970179286436E-3, xyzw + m, .2405702335362910);
      m += gen_oh5(nc, .1771144187504911E-3, xyzw + m, .3165270770189046);
      m += gen_oh5(nc, .1887760022988168E-3, xyzw + m, .3927386145645443);
      m += gen_oh5(nc, .1973474670768214E-3, xyzw + m, .4678825918374656);
      m += gen_oh5(nc, .2033787661234659E-3, xyzw + m, .5408022024266935);
      m += gen_oh5(nc, .2072343626517331E-3, xyzw + m, .6104967445752438);
      m += gen_oh5(nc, .2091177834226918E-3, xyzw + m, .6760910702685738);
      m += gen_oh6(nc, .9316684484675566E-4, xyzw + m, .0665564412021739, .0193650887458842);
      m += gen_oh6(nc, .1116193688682976E-3, xyzw + m, .0944624616127018, .0425244200211587);
      m += gen_oh6(nc, .1298623551559414E-3, xyzw + m, .1242651925452509, .0680652931535437);
      m += gen_oh6(nc, .1450236832456426E-3, xyzw + m, .1553438064846751, .0956095749120537);
      m += gen_oh6(nc, .1572719958149914E-3, xyzw + m, .1871137110542670, .1245931657452888);
      m += gen_oh6(nc, .1673234785867195E-3, xyzw + m, .2192612628836257, .1545385828778978);
      m += gen_oh6(nc, .1756860118725188E-3, xyzw + m, .2515682807206955, .1851004249723368);
      m += gen_oh6(nc, .1826776290439367E-3, xyzw + m, .2838535866287290, .2160182608272384);
      m += gen_oh6(nc, .1885116347992865E-3, xyzw + m, .3159578817528521, .2470799012277111);
      m += gen_oh6(nc, .1933457860170574E-3, xyzw + m, .3477370882791392, .2781014208986402);
      m += gen_oh6(nc, .1973060671902064E-3, xyzw + m, .3790576960890540, .3089172523515731);
      m += gen_oh6(nc, .2004987099616311E-3, xyzw + m, .4097938317810200, .3393750055472244);
      m += gen_oh6(nc, .2030170909281499E-3, xyzw + m, .4398256572859637, .3693322470987730);
      m += gen_oh6(nc, .2049461460119080E-3, xyzw + m, .4690384114718480, .3986541005609877);
      m += gen_oh6(nc, .2063653565200186E-3, xyzw + m, .4973216048301053, .4272112491408562);
      m += gen_oh6(nc, .2073507927381027E-3, xyzw + m, .5245681526132446, .4548781735309936);
      m += gen_oh6(nc, .2079764593256122E-3, xyzw + m, .5506733911803888, .4815315355023251);
      m += gen_oh6(nc, .2083150534968778E-3, xyzw + m, .5755339829522474, .5070486445801855);
      m += gen_oh6(nc, .1262715121590664E-3, xyzw + m, .1305472386056362, .0228497037572237);
      m += gen_oh6(nc, .1414386128545972E-3, xyzw + m, .1637327908216477, .0481225433828838);
      m += gen_oh6(nc, .1538740401313898E-3, xyzw + m, .1972734634149637, .0753173445751193);
      m += gen_oh6(nc, .1642434942331432E-3, xyzw + m, .2308694653110130, .1039043639882017);
      m += gen_oh6(nc, .1729790609237496E-3, xyzw + m, .2643899218338160, .1334526587117626);
      m += gen_oh6(nc, .1803505190260828E-3, xyzw + m, .2977171599622171, .1636414868936382);
      m += gen_oh6(nc, .1865475350079657E-3, xyzw + m, .3307293903032310, .1942195406166568);
      m += gen_oh6(nc, .1917182669679069E-3, xyzw + m, .3633069198219073, .2249752879943753);
      m += gen_oh6(nc, .1959851709034382E-3, xyzw + m, .3953346955922727, .2557218821820032);
      m += gen_oh6(nc, .1994529548117882E-3, xyzw + m, .4267018394184914, .2862897925213193);
      m += gen_oh6(nc, .2022138911146548E-3, xyzw + m, .4573009622571704, .3165224536636518);
      m += gen_oh6(nc, .2043518024208592E-3, xyzw + m, .4870279559856109, .3462730221636496);
      m += gen_oh6(nc, .2059450313018110E-3, xyzw + m, .5157819581450322, .3754016870282835);
      m += gen_oh6(nc, .2070685715318472E-3, xyzw + m, .5434651666465393, .4037733784993613);
      m += gen_oh6(nc, .2077955310694373E-3, xyzw + m, .5699823887764627, .4312557784139123);
      m += gen_oh6(nc, .2081980387824712E-3, xyzw + m, .5952403350947741, .4577175367122110);
      m += gen_oh6(nc, .1521318610377956E-3, xyzw + m, .2025152599210369, .0252025361771956);
      m += gen_oh6(nc, .1622772720185755E-3, xyzw + m, .2381066653274425, .0522325450611900);
      m += gen_oh6(nc, .1710498139420709E-3, xyzw + m, .2732823383651612, .0806066968858862);
      m += gen_oh6(nc, .1785911149448736E-3, xyzw + m, .3080137692611118, .1099335754081255);
      m += gen_oh6(nc, .1850125313687736E-3, xyzw + m, .3422405614587601, .1399120955959857);
      m += gen_oh6(nc, .1904229703933298E-3, xyzw + m, .3758808773890420, .1702977801651705);
      m += gen_oh6(nc, .1949259956121987E-3, xyzw + m, .4088458383438932, .2008799256601680);
      m += gen_oh6(nc, .1986161545363960E-3, xyzw + m, .4410450550841152, .2314703052180836);
      m += gen_oh6(nc, .2015790585641370E-3, xyzw + m, .4723879420561312, .2618972111375892);
      m += gen_oh6(nc, .2038934198707418E-3, xyzw + m, .5027843561874343, .2920013195600270);
      m += gen_oh6(nc, .2056334060538251E-3, xyzw + m, .5321453674452458, .3216322555190551);
      m += gen_oh6(nc, .2068705959462289E-3, xyzw + m, .5603839113834030, .3506456615934198);
      m += gen_oh6(nc, .2076753906106002E-3, xyzw + m, .5874150706875146, .3789007181306267);
      m += gen_oh6(nc, .2081179391734803E-3, xyzw + m, .6131559381660038, .4062580170572782);
      m += gen_oh6(nc, .1700345216228943E-3, xyzw + m, .2778497016394506, .0269627127687623);
      m += gen_oh6(nc, .1774906779990410E-3, xyzw + m, .3143733562261912, .0552346931696047);
      m += gen_oh6(nc, .1839659377002642E-3, xyzw + m, .3501485810261827, .0844519320162646);
      m += gen_oh6(nc, .1894987462975169E-3, xyzw + m, .3851430322303653, .1143263119336083);
      m += gen_oh6(nc, .1941548809452595E-3, xyzw + m, .4193013979470415, .1446177898344475);
      m += gen_oh6(nc, .1980078427252384E-3, xyzw + m, .4525585960458567, .1751165438438091);
      m += gen_oh6(nc, .2011296284744488E-3, xyzw + m, .4848447779622947, .2056338306745660);
      m += gen_oh6(nc, .2035888456966776E-3, xyzw + m, .5160871208276894, .2359965487229226);
      m += gen_oh6(nc, .2054516325352142E-3, xyzw + m, .5462112185696926, .2660430223139146);
      m += gen_oh6(nc, .2067831033092635E-3, xyzw + m, .5751425068101756, .2956193664498032);
      m += gen_oh6(nc, .2076485320284876E-3, xyzw + m, .6028073872853597, .3245763905312779);
      m += gen_oh6(nc, .2081141439525255E-3, xyzw + m, .6291338275278409, .3527670026206972);
      m += gen_oh6(nc, .1834383015469222E-3, xyzw + m, .3541797528439391, .0282385347943555);
      m += gen_oh6(nc, .1889540591777677E-3, xyzw + m, .3908234972074657, .0574129637471311);
      m += gen_oh6(nc, .1936677023597375E-3, xyzw + m, .4264408450107590, .0872464663365020);
      m += gen_oh6(nc, .1976176495066504E-3, xyzw + m, .4609949666553286, .1175034422915616);
      m += gen_oh6(nc, .2008536004560983E-3, xyzw + m, .4944389496536006, .1479755652628428);
      m += gen_oh6(nc, .2034280351712291E-3, xyzw + m, .5267194884346086, .1784740659484352);
      m += gen_oh6(nc, .2053944466027758E-3, xyzw + m, .5577787810220990, .2088245700431244);
      m += gen_oh6(nc, .2068077642882360E-3, xyzw + m, .5875563763536670, .2388628136570763);
      m += gen_oh6(nc, .2077250949661599E-3, xyzw + m, .6159910016391269, .2684308928769185);
      m += gen_oh6(nc, .2082062440705320E-3, xyzw + m, .6430219602956267, .2973740761960252);
      m += gen_oh6(nc, .1934374486546626E-3, xyzw + m, .4300647036213646, .0291639992049398);
      m += gen_oh6(nc, .1974107010484300E-3, xyzw + m, .4661486308935531, .0589880302475566);
      m += gen_oh6(nc, .2007129290388658E-3, xyzw + m, .5009658555287261, .0892416269852541);
      m += gen_oh6(nc, .2033736947471293E-3, xyzw + m, .5344824270447704, .1197185199637321);
      m += gen_oh6(nc, .2054287125902493E-3, xyzw + m, .5666575997416371, .1502300756161382);
      m += gen_oh6(nc, .2069184936818894E-3, xyzw + m, .5974457471404752, .1806004191913564);
      m += gen_oh6(nc, .2078883689808782E-3, xyzw + m, .6267984444116886, .2106621764786252);
      m += gen_oh6(nc, .2083886366116359E-3, xyzw + m, .6546664713575417, .2402526932671914);
      m += gen_oh6(nc, .2006593275470817E-3, xyzw + m, .5042711004437253, .0298252920360766);
      m += gen_oh6(nc, .2033728426135397E-3, xyzw + m, .5392127456774380, .0600872806233992);
      m += gen_oh6(nc, .2055008781377608E-3, xyzw + m, .5726819437668618, .0905822767457140);
      m += gen_oh6(nc, .2070651783518502E-3, xyzw + m, .6046469254207278, .1211219235803400);
      m += gen_oh6(nc, .2080953335094320E-3, xyzw + m, .6350716157434952, .1515286404791580);
      m += gen_oh6(nc, .2086284998988521E-3, xyzw + m, .6639177679185454, .1816314681255552);
      m += gen_oh6(nc, .2055549387644668E-3, xyzw + m, .5757276040972253, .0302699175257544);
      m += gen_oh6(nc, .2071871850267654E-3, xyzw + m, .6090265823139756, .0607840229787077);
      m += gen_oh6(nc, .2082856600431965E-3, xyzw + m, .6406735344387661, .0913545998417664);
      m += gen_oh6(nc, .2088705858819358E-3, xyzw + m, .6706397927793709, .1218024155966590);
      m += gen_oh6(nc, .2083995867536322E-3, xyzw + m, .6435019674426665, .0305260835766064);
      m += gen_oh6(nc, .2090509712889637E-3, xyzw + m, .6747218676375681, .0611218577398309);
    break; case 5810:
      m += gen_oh1(nc, .9735347946175486E-5, xyzw + m);
      m += gen_oh2(nc, .1907581241803167E-3, xyzw + m);
      m += gen_oh3(nc, .1901059546737578E-3, xyzw + m);
      m += gen_oh4(nc, .3926424538919212E-4, xyzw + m, .0118236166240028);
      m += gen_oh4(nc, .6667905467294381E-4, xyzw + m, .0306214500913896);
      m += gen_oh4(nc, .8868891315019136E-4, xyzw + m, .0532979403683424);
      m += gen_oh4(nc, .1066306000958872E-3, xyzw + m, .0784816553286222);
      m += gen_oh4(nc, .1214506743336128E-3, xyzw + m, .1054038157636201);
      m += gen_oh4(nc, .1338054681640871E-3, xyzw + m, .1335577797766211);
      m += gen_oh4(nc, .1441677023628504E-3, xyzw + m, .1625769955502252);
      m += gen_oh4(nc, .1528880200826557E-3, xyzw + m, .1921787193412792);
      m += gen_oh4(nc, .1602330623773609E-3, xyzw + m, .2221340534690548);
      m += gen_oh4(nc, .1664102653445244E-3, xyzw + m, .2522504912791132);
      m += gen_oh4(nc, .1715845854011323E-3, xyzw + m, .2823610860679697);
      m += gen_oh4(nc, .1758901000133069E-3, xyzw + m, .3123173966267560);
      m += gen_oh4(nc, .1794382485256736E-3, xyzw + m, .3419847036953789);
      m += gen_oh4(nc, .1823238106757407E-3, xyzw + m, .3712386456999758);
      m += gen_oh4(nc, .1846293252959976E-3, xyzw + m, .3999627649876828);
      m += gen_oh4(nc, .1864284079323098E-3, xyzw + m, .4280466458648093);
      m += gen_oh4(nc, .1877882694626914E-3, xyzw + m, .4553844360185711);
      m += gen_oh4(nc, .1887716321852025E-3, xyzw + m, .4818736094437834);
      m += gen_oh4(nc, .1894381638175673E-3, xyzw + m, .5074138709260629);
      m += gen_oh4(nc, .1898454899533629E-3, xyzw + m, .5319061304570707);
      m += gen_oh4(nc, .1900497929577815E-3, xyzw + m, .5552514978677286);
      m += gen_oh4(nc, .1900671501924092E-3, xyzw + m, .5981009025246183);
      m += gen_oh4(nc, .1899837555533510E-3, xyzw + m, .6173990192228116);
      m += gen_oh4(nc, .1899014113156229E-3, xyzw + m, .6351365239411131);
      m += gen_oh4(nc, .1898581257705106E-3, xyzw + m, .6512010228227200);
      m += gen_oh4(nc, .1898804756095753E-3, xyzw + m, .6654758363948120);
      m += gen_oh4(nc, .1899793610426402E-3, xyzw + m, .6778410414853370);
      m += gen_oh4(nc, .1901464554844117E-3, xyzw + m, .6881760887484110);
      m += gen_oh4(nc, .1903533246259542E-3, xyzw + m, .6963645267094598);
      m += gen_oh4(nc, .1905556158463228E-3, xyzw + m, .7023010617153579);
      m += gen_oh4(nc, .1907037155663528E-3, xyzw + m, .7059004636628753);
      m += gen_oh5(nc, .5992997844249967E-4, xyzw + m, .0355247031247257);
      m += gen_oh5(nc, .9749059382456977E-4, xyzw + m, .0915117662084128);
      m += gen_oh5(nc, .1241680804599158E-3, xyzw + m, .1566197930068980);
      m += gen_oh5(nc, .1437626154299360E-3, xyzw + m, .2265467599271907);
      m += gen_oh5(nc, .1584200054793902E-3, xyzw + m, .2988242318581361);
      m += gen_oh5(nc, .1694436550982744E-3, xyzw + m, .3717482419703886);
      m += gen_oh5(nc, .1776617014018108E-3, xyzw + m, .4440094491758889);
      m += gen_oh5(nc, .1836132434440077E-3, xyzw + m, .5145337096756643);
      m += gen_oh5(nc, .1876494727075983E-3, xyzw + m, .5824053672860230);
      m += gen_oh5(nc, .1899906535336482E-3, xyzw + m, .6468283961043370);
      m += gen_oh6(nc, .8143252820767350E-4, xyzw + m, .0609596425910437, .0178782827534293);
      m += gen_oh6(nc, .9998859890887728E-4, xyzw + m, .0881196227095939, .0395388874079210);
      m += gen_oh6(nc, .1156199403068359E-3, xyzw + m, .1165936722428831, .0637812179772299);
      m += gen_oh6(nc, .1287632092635513E-3, xyzw + m, .1460232857031785, .0898589081374504);
      m += gen_oh6(nc, .1398378643365139E-3, xyzw + m, .1761197110181755, .1172606510576162);
      m += gen_oh6(nc, .1491876468417391E-3, xyzw + m, .2066471190463718, .1456102876970995);
      m += gen_oh6(nc, .1570855679175456E-3, xyzw + m, .2374076026328152, .1746153823011775);
      m += gen_oh6(nc, .1637483948103775E-3, xyzw + m, .2682305474337051, .2040383070295584);
      m += gen_oh6(nc, .1693500566632843E-3, xyzw + m, .2989653312142369, .2336788634003698);
      m += gen_oh6(nc, .1740322769393633E-3, xyzw + m, .3294762752772209, .2633632752654219);
      m += gen_oh6(nc, .1779126637278296E-3, xyzw + m, .3596390887276086, .2929369098051601);
      m += gen_oh6(nc, .1810908108835412E-3, xyzw + m, .3893383046398812, .3222592785275512);
      m += gen_oh6(nc, .1836529132600190E-3, xyzw + m, .4184653789358347, .3512004791195743);
      m += gen_oh6(nc, .1856752841777379E-3, xyzw + m, .4469172319076166, .3796385677684537);
      m += gen_oh6(nc, .1872270566606832E-3, xyzw + m, .4745950813276976, .4074575378263879);
      m += gen_oh6(nc, .1883722645591307E-3, xyzw + m, .5014034601410262, .4345456906027828);
      m += gen_oh6(nc, .1891714324525297E-3, xyzw + m, .5272493404551239, .4607942515205134);
      m += gen_oh6(nc, .1896827480450146E-3, xyzw + m, .5520413051846366, .4860961284181720);
      m += gen_oh6(nc, .1899628417059528E-3, xyzw + m, .5756887237503077, .5103447395342789);
      m += gen_oh6(nc, .1123301829001669E-3, xyzw + m, .1225039430588352, .0213645592265579);
      m += gen_oh6(nc, .1253698826711277E-3, xyzw + m, .1539113217321372, .0452092616613719);
      m += gen_oh6(nc, .1366266117678531E-3, xyzw + m, .1856213098637712, .0708646817786482);
      m += gen_oh6(nc, .1462736856106918E-3, xyzw + m, .2174998728035131, .0978523948877292);
      m += gen_oh6(nc, .1545076466685412E-3, xyzw + m, .2494128336938330, .1258106396267210);
      m += gen_oh6(nc, .1615096280814007E-3, xyzw + m, .2812321562143480, .1544529125047001);
      m += gen_oh6(nc, .1674366639741759E-3, xyzw + m, .3128372276456111, .1835433512202753);
      m += gen_oh6(nc, .1724225002437900E-3, xyzw + m, .3441145160177973, .2128813258619585);
      m += gen_oh6(nc, .1765810822987288E-3, xyzw + m, .3749567714853510, .2422913734880829);
      m += gen_oh6(nc, .1800104126010751E-3, xyzw + m, .4052621732015610, .2716163748391453);
      m += gen_oh6(nc, .1827960437331284E-3, xyzw + m, .4349335453522385, .3007127671240280);
      m += gen_oh6(nc, .1850140300716308E-3, xyzw + m, .4638776641524965, .3294470677216479);
      m += gen_oh6(nc, .1867333507394938E-3, xyzw + m, .4920046410462687, .3576932543699155);
      m += gen_oh6(nc, .1880178688638289E-3, xyzw + m, .5192273554861704, .3853307059757764);
      m += gen_oh6(nc, .1889278925654758E-3, xyzw + m, .5454609081136522, .4122425044452694);
      m += gen_oh6(nc, .1895213832507346E-3, xyzw + m, .5706220661424140, .4383139587781027);
      m += gen_oh6(nc, .1898548277397420E-3, xyzw + m, .5946286755181518, .4634312536300553);
      m += gen_oh6(nc, .1349105935937341E-3, xyzw + m, .1905370790924295, .0237131153778198);
      m += gen_oh6(nc, .1444060068369326E-3, xyzw + m, .2242518717748009, .0491787805925481);
      m += gen_oh6(nc, .1526797390930008E-3, xyzw + m, .2577190808025936, .0759549896049514);
      m += gen_oh6(nc, .1598208771406474E-3, xyzw + m, .2908724534927187, .1036991083191100);
      m += gen_oh6(nc, .1659354368615331E-3, xyzw + m, .3236354020056219, .1321348584450234);
      m += gen_oh6(nc, .1711279910946440E-3, xyzw + m, .3559267359304543, .1610316571314789);
      m += gen_oh6(nc, .1754952725601440E-3, xyzw + m, .3876637123676956, .1901912080395707);
      m += gen_oh6(nc, .1791247850802529E-3, xyzw + m, .4187636705218842, .2194384950137950);
      m += gen_oh6(nc, .1820954300877716E-3, xyzw + m, .4491449019883107, .2486155334763858);
      m += gen_oh6(nc, .1844788524548449E-3, xyzw + m, .4787270932425445, .2775768931812335);
      m += gen_oh6(nc, .1863409481706220E-3, xyzw + m, .5074315153055574, .3061863786591120);
      m += gen_oh6(nc, .1877433008795068E-3, xyzw + m, .5351810507738336, .3343144718152556);
      m += gen_oh6(nc, .1887444543705232E-3, xyzw + m, .5619001025975381, .3618362729028427);
      m += gen_oh6(nc, .1894009829375006E-3, xyzw + m, .5875144035268046, .3886297583620408);
      m += gen_oh6(nc, .1897683345035198E-3, xyzw + m, .6119507308734495, .4145742277792031);
      m += gen_oh6(nc, .1517327037467653E-3, xyzw + m, .2619733870119463, .0254004718638935);
      m += gen_oh6(nc, .1587740557483543E-3, xyzw + m, .2968149743237949, .0520810701854399);
      m += gen_oh6(nc, .1649093382274097E-3, xyzw + m, .3310451504860488, .0797182847088560);
      m += gen_oh6(nc, .1701915216193265E-3, xyzw + m, .3646215567376676, .1080465999177927);
      m += gen_oh6(nc, .1746847753144065E-3, xyzw + m, .3974916785279360, .1368413849366629);
      m += gen_oh6(nc, .1784555512007570E-3, xyzw + m, .4295967403772029, .1659073184763559);
      m += gen_oh6(nc, .1815687562112174E-3, xyzw + m, .4608742854473447, .1950703730454614);
      m += gen_oh6(nc, .1840864370663302E-3, xyzw + m, .4912598858949903, .2241721144376724);
      m += gen_oh6(nc, .1860676785390006E-3, xyzw + m, .5206882758945558, .2530655255406489);
      m += gen_oh6(nc, .1875690583743703E-3, xyzw + m, .5490940914019820, .2816118409731066);
      m += gen_oh6(nc, .1886453236347225E-3, xyzw + m, .5764123302025542, .3096780504593238);
      m += gen_oh6(nc, .1893501123329645E-3, xyzw + m, .6025786004213506, .3371348366394987);
      m += gen_oh6(nc, .1897366184519868E-3, xyzw + m, .6275291964794956, .3638547827694396);
      m += gen_oh6(nc, .1643908815152736E-3, xyzw + m, .3348189479861771, .0266484193553744);
      m += gen_oh6(nc, .1696300350907768E-3, xyzw + m, .3699515545855295, .0542400006684349);
      m += gen_oh6(nc, .1741553103844483E-3, xyzw + m, .4042003071474669, .0825199271543085);
      m += gen_oh6(nc, .1780015282386092E-3, xyzw + m, .4375320100182624, .1112695182483710);
      m += gen_oh6(nc, .1812116787077125E-3, xyzw + m, .4699054490335947, .1402964116467816);
      m += gen_oh6(nc, .1838323158085421E-3, xyzw + m, .5012739879431952, .1694275117584291);
      m += gen_oh6(nc, .1859113119837737E-3, xyzw + m, .5315874883754966, .1985038235312689);
      m += gen_oh6(nc, .1874969220221698E-3, xyzw + m, .5607937109622116, .2273765660020893);
      m += gen_oh6(nc, .1886375612681076E-3, xyzw + m, .5888393223495521, .2559041492849764);
      m += gen_oh6(nc, .1893819575809276E-3, xyzw + m, .6156705979160163, .2839497251976899);
      m += gen_oh6(nc, .1897794748256767E-3, xyzw + m, .6412338809078123, .3113791060500690);
      m += gen_oh6(nc, .1738963926584846E-3, xyzw + m, .4076051259257167, .0275779229085846);
      m += gen_oh6(nc, .1777442359873466E-3, xyzw + m, .4423788125791520, .0558413683498429);
      m += gen_oh6(nc, .1810010815068719E-3, xyzw + m, .4760480917328258, .0845777208772714);
      m += gen_oh6(nc, .1836920318248129E-3, xyzw + m, .5085838725946297, .1135975846359248);
      m += gen_oh6(nc, .1858489473214328E-3, xyzw + m, .5399513637391218, .1427286904765053);
      m += gen_oh6(nc, .1875079342496592E-3, xyzw + m, .5701118433636380, .1718112740057635);
      m += gen_oh6(nc, .1887080239102310E-3, xyzw + m, .5990240530606021, .2006944855985351);
      m += gen_oh6(nc, .1894905752176822E-3, xyzw + m, .6266452685139695, .2292335090598907);
      m += gen_oh6(nc, .1898991061200695E-3, xyzw + m, .6529320971415942, .2572871512353714);
      m += gen_oh6(nc, .1809065016458791E-3, xyzw + m, .4791583834610126, .0282609419773593);
      m += gen_oh6(nc, .1836297121596799E-3, xyzw + m, .5130373952796941, .0569987135968365);
      m += gen_oh6(nc, .1858426916241869E-3, xyzw + m, .5456252429628476, .0860271252855439);
      m += gen_oh6(nc, .1875654101134641E-3, xyzw + m, .5768956329682385, .1151748137221281);
      m += gen_oh6(nc, .1888240751833503E-3, xyzw + m, .6068186944699046, .1442811654136362);
      m += gen_oh6(nc, .1896497383866979E-3, xyzw + m, .6353622248024907, .1731930321657680);
      m += gen_oh6(nc, .1900775530219121E-3, xyzw + m, .6624927035731797, .2017619958756061);
      m += gen_oh6(nc, .1858525041478814E-3, xyzw + m, .5484933508028488, .0287421975590739);
      m += gen_oh6(nc, .1876248690077947E-3, xyzw + m, .5810207682142106, .0577831212371369);
      m += gen_oh6(nc, .1889404439064607E-3, xyzw + m, .6120955197181353, .0869526237143953);
      m += gen_oh6(nc, .1898168539265290E-3, xyzw + m, .6416944284294319, .1160893767057166);
      m += gen_oh6(nc, .1902779940661772E-3, xyzw + m, .6697926391731260, .1450378826743251);
      m += gen_oh6(nc, .1890125641731815E-3, xyzw + m, .6147594390585488, .0290495762234146);
      m += gen_oh6(nc, .1899434637795751E-3, xyzw + m, .6455390026356783, .0582380915261720);
      m += gen_oh6(nc, .1904520856831751E-3, xyzw + m, .6747258588365477, .0874038489988472);
      m += gen_oh6(nc, .1905534498734563E-3, xyzw + m, .6772135750395347, .0291994613580811);
#endif
    break; default:
        if(echo > 0) printf("# %s A Lebedev-Laikov grid with %d points for ellmax= %d is not available!\n", __func__, n_corrected, ellmax); 
#ifdef  LARGE_GRIDS
#else       
        if(echo > 0 && n_corrected >= 770) printf("# %s Please activate -D LARGE_GRIDS in %s\n", __func__, __FILE__); 
#endif
        return -1;
    } // switch n

    if (echo > 6 && m > 0) {
        int8_t const npoints[8] = {1, 6, 12, 8, 24, 24, 48, 0};
        printf("# %s The Lebedev-Laikov grid for ellmax= %d has %d = ", __func__, ellmax, n_corrected);
        int n_check = 0;
        for(int icode = 0; icode < 7; ++icode) {
            if (nc[icode] > 0) printf("%d*%d + ", nc[icode], npoints[icode]);
            n_check += nc[icode] * npoints[icode];
        } // icode
        printf("0 points\n");
        assert(n_corrected == n_check);
    } // report
    
    assert (n_corrected == m);
    return m;
  } // create_Lebedev_grid

  
  
#if 0 
// def USE_LaTeX
cTeXit '', '%%%'+__FILE__+'line='-__LINE__,''
cTeXit '', 'For the setup of the angular grid with $'-m-'$ points on the unit sphere'

  if(nc(0) == 1) then
cTeXit ', a single point was used' !, ' instead of an angular grid'
  endif ! nc

  if(nc(1) == 1) then ; i = 1
cTeXit ', an angular sub-grid with $'-nPoints(i)-'$ points', &
cTeX   ' according to the scheme of Lebedev \emph{et al.} \cite{Lebedev'-CiteLebedevYear(i)-'} was added'
  endif ! nc

  do i = 2, 6
    if(nc(i) == 1) then
cTeXit ', an angular sub-grid with $'-nPoints(i)-'$ points', &
cTeX   ' according to the scheme of Lebedev \cite{Lebedev'-CiteLebedevYear(i)-'} was added'
    elseif(nc(i) > 1) then
cTeXit ','+nc(i)+'angular sub-grids with $'-nPoints(i)-'$ points each', &
cTeX   ' according to the scheme of Lebedev \cite{Lebedev'-CiteLebedevYear(i)-'} were added'
    endif ! nc
  enddo ! i
cTeXit '. ', '' ! full stop and an extra empty line
#endif

// interface documentation for gen_ohN()
//   int gen_ohN(int ncode[], // (inout) count how many times this icode was used
//              real_t const w8, // (in) weight 
//              real_t xyzw[][4], // (out) weights 3:weights, 0:x 1:y 2:z
//              double const aa, double const bb);
  
#if 0
!!! from http://server.ccl.net/cca/software/SOURCES/FORTRAN/Lebedev-Laikov-Grids/Lebedev-Laikov.F
! chvd
! chvd   This subroutine is part of a set of subroutines that generate
! chvd   Lebedev grids [1-6] for integration on a sphere. The original 
! chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
! chvd   translated into fortran by Dr. Christoph van Wuellen.
! chvd   This subroutine was translated from C to fortran77 by hand.
! chvd
! chvd   Users of this code are asked to include reference [1] in their
! chvd   publications, and in the user- and programmers-manuals 
! chvd   describing their codes.
! chvd
! chvd   This code was distributed through CCL (http://www.ccl.net].
! chvd
! chvd   [1] V.I. Lebedev, and D.N. Laikov
! chvd       "A quadrature formula for the sphere of the 131st algebraic order of accuracy"
! chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
! chvd
! chvd   [2] V.I. Lebedev
! chvd       "A quadrature formula for the sphere of 59th algebraic order of accuracy"
! chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
! chvd
! chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
! chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
! chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
! chvd
! chvd   [4] V.I. Lebedev
! chvd       "Spherical quadrature formulas exact to orders 25-29"
! chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
! chvd
! chvd   [5] V.I. Lebedev
! chvd       "Quadratures on a sphere"
! chvd       Computational Mathematics and Mathematical Physics, Vol. 16, 1976, pp. 10-24. 
! chvd
! chvd   [6] V.I. Lebedev
! chvd       "Values of the nodes and weights of ninth to seventeenth order Gauss-Markov quadrature formulae invariant under the octahedron group with inversion"
! chvd       Computational Mathematics and Mathematical Physics, Vol. 15, 1975, pp. 44-51.
! chvd
! cvw
! cvw    Given a point on a sphere (specified by a and b), generate all
! cvw    the equivalent points under Oh symmetry, making grid points with
! cvw    weight v.
! cvw    The variable num is increased by the number of different points
! cvw    generated.
! cvw
! cvw    Depending on icode, there are 6...48 different but equivalent
! cvw    points.
! cvw
! cvw    icode=0:   (0,0,0) only                               ( 1 point)
! cvw    icode=1:   (0,0,1) etc                                ( 6 points)
! cvw    icode=2:   (0,a,a) etc, a=1/sqrt(2)                   (12 points)
! cvw    icode=3:   (a,a,a) etc, a=1/sqrt(3)                   ( 8 points)
! cvw    icode=4:   (a,a,b) etc, b=sqrt(1-2 a^2)               (24 points)
! cvw    icode=5:   (a,b,0) etc, b=sqrt(1-a^2), a input        (24 points)
! cvw    icode=6:   (a,b,c) etc, c=sqrt(1-a^2-b^2), a/b input  (48 points)
! cvw
#endif
  
  template<int lmax>
  status_t create_numerical_Gaunt(std::vector<gaunt_entry_t>* gaunt, int const echo) {
      int const npt = Lebedev_grid_size(2*lmax);
      int const M0 = (1 + 2*lmax)*(1 + 2*lmax), M = (1 + lmax)*(1 + lmax);
      auto yy   = new double[npt][M0];
      auto xyzw = new double[npt][4];
      create_Lebedev_grid(2*lmax, xyzw, echo);
      double const pi = constants::pi;
      for(int ipt = 0; ipt < npt; ++ipt) {
          solid_harmonics::Xlm(yy[ipt], 2*lmax, xyzw[ipt]);
      } // ipt
      int const n_expected = (M*M*M0) >> 5;
      if (gaunt) gaunt->reserve(n_expected);
      int n = 0, nnz = 0;
      for(int lm0 = 0; lm0 < M0; ++lm0) {
          for(int16_t lm1 = 0; lm1 < M; ++lm1) {
              for(int16_t lm2 = 0; lm2 < M; ++lm2) {
                  double Gaunt = 0;
                  for(int ipt = 0; ipt < npt; ++ipt) {
                      double const w8 = xyzw[ipt][3];
                      Gaunt += yy[ipt][lm0] * yy[ipt][lm1] * yy[ipt][lm2] * w8;
                  } // ipt
                  Gaunt *= 4*pi;
                  if (std::abs(Gaunt) > 1e-14) {
                      if (echo > 1) printf("%d %d %d %.9f\n", lm0, lm1, lm2, Gaunt);
                      if (gaunt) gaunt->push_back({Gaunt, lm0, lm1, lm2});
                      ++nnz;
                  } // non-zero
                  ++n;
              } // lm2
          } // lm1
      } // lm0
      if (echo > 0) printf("\n# %s: %d of %d real-Gaunt tensor elements are nonzero (%.3f %%), expected %d\n", 
                                __func__, nnz, n, nnz/(.01*n), n_expected);
      return 0;
  } // create_numerical_Gaunt
  
  template status_t create_numerical_Gaunt<6>(std::vector<gaunt_entry_t>* gaunt, int const echo); // explicit template instanciation
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  int test_generation(int echo=1) {
      if (echo > 1) printf("\n# %s: \n", __func__);
      int const max_size = Lebedev_grid_size(ellmax_implemented);
      auto xyzw = new double[max_size][4];
      status_t stat = 0;
      for(int ell = -2; ell <= ellmax_implemented + 3; ell += (1 + 2*(ell > 16))) { // in steps of 3 for the larger grids, no need to test the same n 3 times
          if (echo > 3) printf("\n# %s: try to generate Lebedev-Laikov grid with for ellmax= %d\n", __func__, ell);
          auto const m = create_Lebedev_grid(ell, xyzw, echo);
          if (echo > 2) printf("\n# %s: generated Lebedev-Laikov grid with for ellmax= %d using %d points\n", __func__, ell, m);
          if (ell >= 0 && ell <= ellmax_implemented) stat += (m < 0); // count the number of failures in the relevant range
          if (m > 1) {
              // sanity checks on weights and coordinates of the points on the unit sphere
              double w8sum = 0, max_norm2_dev = 0;
              for(int i = 0; i < m; ++i) {
                  double const norm2 = xyzw[i][0]*xyzw[i][0] + xyzw[i][1]*xyzw[i][1] + xyzw[i][2]*xyzw[i][2];
                  max_norm2_dev = std::max(max_norm2_dev, std::abs(1. - norm2));
                  w8sum += xyzw[i][3];
              } // i
              if (echo > 4) printf("# Lebedev-Laikov grid with %d points for ellmax= %d: w8_dev= %g, |v|_dev= %g\n", m, ell, 1. - w8sum, max_norm2_dev);
              assert(std::abs(1. - w8sum) < 2.3e-14);
              assert(max_norm2_dev < 2.3e-16);
          } // m > 1, for m == 1, a single weight is set to unity but the coordinates are all zero, so the sanity test does not apply
      } // ell
      delete[] xyzw;
      if (echo > 1) printf("\n# %s: %d grid generations failed!\n", __func__, stat);
      return stat;
  } // test
    
  int test_orthogonality(int echo=9, int lmax=20) { // terribly slow
      printf("\n# %s: \n", __func__);
      int const ellmax = std::min(lmax, ellmax_implemented);
      int const max_size = Lebedev_grid_size(ellmax);
      int const M = (1 + ellmax)*(1 + ellmax);
//       typedef std::complex<double> Ylm_t; // complex
      typedef double Ylm_t;
      auto yy = new Ylm_t[M];
      auto unity = new Ylm_t[M*M];
      auto xyzw = new double[max_size][4];
      status_t stat = 0;
      double dev_all = 0;
      for(int ell = ellmax; ell > 0; --ell) {
          int const m = (1 + ell)*(1 + ell);
          if (echo > 3) printf("\n# %s: try orthogonality on Lebedev-Laikov grid with for ellmax= %d\n", __func__, ell);
          auto const np = create_Lebedev_grid(ell, xyzw, echo);
          for(int ij = 0; ij < M*M; ++ij) unity[ij] = 0; // clear
          for(int ip = 0; ip < np; ++ip) {
              auto const w8 = xyzw[ip][3] * 4*constants::pi;
//            if (echo > 3) printf("# %s: envoke Ylm for %d  %g %g %g\n", __func__, ell, xyzw[ip][0], xyzw[ip][1], xyzw[ip][2]);
//               spherical_harmonics::Ylm(yy, ell, xyzw[ip]); // complex
              solid_harmonics::Xlm(yy, ell, xyzw[ip]);
              for(int i = 0; i < m; ++i) {
//                auto const yi = std::conj(yy[i])*w8; // complex
                  auto const yi = yy[i]*w8;
                  for(int j = 0; j < m; ++j) {
                      unity[i*M + j] += yi * yy[j];
                  } // j
              } // i
          } // ip
          
          double dev = 0;
          for(int i = 0; i < m; ++i) {
              for(int j = 0; j < m; ++j) {
//                   auto const Re = unity[i*M + j].real(), Im = unity[i*M + j].imag();
//                   if (echo > 8) printf("%g %g  ", Re, Im);
//                   dev = std::max(dev, std::abs(Im));
                  auto const Re = unity[i*M + j];
//                   if (echo > 8) printf("%g ", Re);
                  if (echo > 8) printf("%16.9f ", Re);
                  dev = std::max(dev, std::abs(Re - (i == j))); // subtract 1 on the diagonal
              } // j
              if (echo > 7) printf("\n");
          } // i
          if (echo > 2) printf("# %s: orthogonality on Lebedev-Laikov grid with for ellmax= %d is %g\n", __func__, ell, dev);
          dev_all = std::max(dev_all, dev);
          stat += (dev > 2.7e-14);

      } // ell
      if (echo > 0) printf("\n# %s: orthogonality on Lebedev-Laikov grid with for ellmax up to %d is %g\n", __func__, ellmax, dev_all);
      delete[] xyzw;
      delete[] yy;
      if (echo > 1) printf("\n# %s: %d grid generations failed!\n", __func__, stat);
      return stat;
  } // test

  status_t test_numerical_Gaunt(int const echo=1) { return create_numerical_Gaunt<3>(nullptr, echo); }

  status_t all_tests() {
    auto status = 0;
    status += test_generation(1);
    status += test_orthogonality(4);
    status += test_numerical_Gaunt();
    
    status += solid_harmonics::test();
    solid_harmonics::cleanup<double>(); // free internal memory
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS
  
} // namespace angular_grid
