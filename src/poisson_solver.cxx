#include <cstdio> // printf, std::sprintf
#include <cassert> // assert
#include <algorithm> // std::copy
#include <cmath> // std::floor
#include <vector> // std::vector

#include "poisson_solver.hxx" // ::all_tests, ::solve, ::print_direct_projection

#include "display_units.h" // eV, _eV, Ang, _Ang
#include "inline_math.hxx" // set, pow2, align<nBits>
// #include "constants.hxx" // ::sqrtpi, ::pi
#include "solid_harmonics.hxx" // ::Y00
#include "real_space.hxx" // ::grid_t, ::Bessel_projection
#include "radial_grid.hxx" // ::radial_grid_t

#include "data_view.hxx" // view2D<T>

#include "fourier_poisson.hxx" // ::fourier_solve
#include "iterative_poisson.hxx" // ::solve
#include "finite_difference.hxx" // ::stencil_t, ::derive

#include "simple_timer.hxx" // // SimpleTimer
#include "control.hxx" // ::get

#include "boundary_condition.hxx" // Periodic_Boundary, Isolated_Boundary
#include "bessel_transform.hxx" // ::Bessel_j0
#include "debug_tools.hxx" // ::read_from_file

#ifdef DEVEL
    #include "lossful_compression.hxx" // print_compressed
    #include "debug_output.hxx" // dump_to_file, ::write_array_to_file, 
    #include "radial_r2grid.hxx" // radial_r2grid_t
    #include "radial_r2grid.hxx" // r2_axis
#endif // DEVEL

#include "multi_grid.hxx" // ::restrict3D, ::interpolate3D
#include "print_tools.hxx" // print_stats, printf_vector

namespace poisson_solver {
  // this module contains a collection of Poisson solvers and control the selection

  status_t Bessel_Poisson_solver(
        double Ves[] // output electrostatic potential
      , double const rho[] // input density
      , real_space::grid_t const & g // Cartesian real-space grid descriptor
      , double const center[3]=nullptr // center around which to expand, e.g. an atomic position
      , int const echo=0 // log-level
  ) {
      if (!g.is_Cartesian()) error("Bessel-Poisson solver only implemented for Cartesian grids", 0);
      // solve the Poisson equation for spherically symmetric geometries using a Bessel trasform
      status_t stat(0);

      { // scope: check that we have 3 isolated boundary conditions
          int const nibc = g.number_of_boundary_conditions(Isolated_Boundary);
          if (nibc < 3) {
              warn("Bessel Poisson solver requires 3 but found only %d isolated boundary conditions", nibc);
              return -1;
          } // failed
      } // scope

      // report extremal values of what is stored on the grid
      if (echo > 1) print_stats(rho, g.all(), g.dV(), "# real-space stats of input density:   ");
      
      double const grid_center[] = {(g[0] - 1)*g.h[0], (g[1] - 1)*g.h[1], (g[2] - 1)*g.h[2]};
      if (nullptr == center) center = grid_center;
      
      if (echo > 5) printf("# Bessel j0 projection around position %g %g %g %s\n",
                              center[0]*Ang, center[1]*Ang, center[2]*Ang, _Ang);
      float const dq = 1.f/16; int const nq = int(constants::pi/(g.smallest_grid_spacing()*dq));
      std::vector<double> qc(nq, 0.0);
      stat += real_space::Bessel_projection(qc.data(), nq, dq, rho, g, center);
      qc[0] = 0; // stabilize charge neutrality, q=0-component must vanish

      double const by4pi = .25/constants::pi; // 1/(4*pi)

      if (echo > 19) {
          printf("\n## Bessel coeff of {density, Ves}:\n");
          for(int iq = 0; iq < nq; ++iq) {
              double const q = iq*dq;
              printf("%g %g %g\n", q, qc[iq]*by4pi, qc[iq]/(q*q + 1e-12));
          } // iq
          printf("\n\n");
      } // echo

      double const sqrt2pi = std::sqrt(2./constants::pi); // this makes the transform symmetric
      scale(qc.data(), nq, sqrt2pi*dq);

      // expand electrostatic potential onto grid
      double rho_abs{0}, rho_squ{0}, rho_max{0};
      for(int iz = 0; iz < g[2]; ++iz) {
          double const z = iz*g.h[2] - center[2], z2 = z*z;
          for(int iy = 0; iy < g[1]; ++iy) {
              double const y = iy*g.h[1] - center[1], y2 = y*y; 
              for(int ix = 0; ix < g[0]; ++ix) {
                  double const x = ix*g.h[0] - center[0], x2 = x*x;
                  double const r = std::sqrt(x2 + y2 + z2);
                  int const izyx = (iz*g[1] + iy)*g[0] + ix;
                  double ves_r{0}, rho_r{0};
                  for(int iq = 0; iq < nq; ++iq) {
                      double const q = iq*dq;
                      double const bessel_kernel = bessel_transform::Bessel_j0(q*r);
                      ves_r += qc[iq]*bessel_kernel; // cheap Poisson solver
                      rho_r += qc[iq]*bessel_kernel * by4pi * q*q; // reconstruct a spherical density
                  } // iq
                  Ves[izyx] = ves_r; // store
                  // check how much the density deviates from a spherical density
                  double const rho_dif = std::abs(rho[izyx] - rho_r);
                  rho_max = std::max(rho_max, rho_dif);
                  rho_abs += rho_dif;
                  rho_squ += pow2(rho_dif);
              } // ix
          } // iy
      } // iz
      rho_abs *= g.dV();
      rho_squ = std::sqrt(rho_squ)*g.dV();

      if (echo > 1) print_stats(Ves, g.all(), g.dV(), "# real-space stats of output potential:", eV);

      if (echo > 3) printf("# Bessel_Poisson deviation from a spherical density: "
          "max %.1e a.u., abs %.1e e, rms %.1e e\n", rho_max, rho_abs, rho_squ);

      if (echo > 17) {
          printf("\n## Bessel_Poisson: r, Ves, rho (in a.u.):\n");
          for(int ir = 0; ir < 99; ++ir) {
              double const r = .1*ir;
              double ves_r{0}, rho_r{0};
              for(int iq = 0; iq < nq; ++iq) {
                  double const q = iq*dq;
                  double const bessel_kernel = bessel_transform::Bessel_j0(q*r);
                  ves_r += qc[iq]*bessel_kernel; // cheap Poisson solver
                  rho_r += qc[iq]*bessel_kernel * by4pi * q*q; // reconstruct a spherical density
              } // iq
              printf("%g %g %g \n", r, ves_r, rho_r);
          } // ir
          printf("\n\n");
      } // echo
      
      return stat;
  } // Bessel_Poisson_solver


  
  
  
  status_t solve(
        double Ves[] // result electrostatic potential
      , double const rho[] // charge density, typically augmented density, should be charge neutral
      , real_space::grid_t const & g // grid descriptor
      , char const method // ='m' solver method     
      , int const echo // log-level
      , double const Bessel_center[] // =nullptr
  ) {
      char const es_solver_name = method;

      status_t stat(0);

      if (echo > 1) print_stats(rho, g.all(), g.dV(), "\n# input charge density:");

      { // scope: solve the Poisson equation: Laplace Ves == -4 pi rho
          SimpleTimer timer(__FILE__, __LINE__, "Poisson equation", echo);
#ifdef DEVEL
          if (echo > 0) {
              printf("\n# Solve Poisson equation\n\n");
              std::fflush(stdout); // if the Poisson solver takes long, we can already see the output up to here
          } // echo
#endif // DEVEL

          if ('f' == (method | 32)) { // "fft", "fourier" 
              // solve the Poisson equation using a Fast Fourier Transform
              int ng[3]; double reci[3][4]; 
              for(int d = 0; d < 3; ++d) { 
                  ng[d] = g[d];
                  set(reci[d], 4, 0.0);
                  reci[d][d] = 2*constants::pi/(ng[d]*g.h[d]);
              } // d
              stat += fourier_poisson::solve(Ves, rho, ng, reci);
          } else
          if ('M' == method) { // "Multi-grid" (upper case!)
#ifdef DEVEL
              // create a 2x denser grid descriptor
              real_space::grid_t gd(g[0]*2, g[1]*2, g[2]*2);
              if (echo > 2) printf("# electrostatic.solver=%c (Multi-grid) is a multi-grid solver"
                      " on a %d x %d x %d grid\n", es_solver_name, gd[0], gd[1], gd[2]);
              gd.set_grid_spacing(g.h[0]/2, g.h[1]/2, g.h[2]/2);
              gd.set_boundary_conditions(g.boundary_conditions());

              std::vector<double> Ves_dense(gd.all()), rho_dense(gd.all());
              multi_grid::interpolate3D(rho_dense.data(), gd, rho, g, echo);

              iterative_poisson::solve(Ves_dense.data(), rho_dense.data(), gd, 'M', echo);

              // restrict the electrostatic potential to grid g
              multi_grid::restrict3D(Ves, g, Ves_dense.data(), gd, echo);
#else  // DEVEL
              error("electrostatic.solver=%c (Multi-grid) only available in the development branch!", es_solver_name); 
#endif // DEVEL
          } else
          if ('B' == method) { // "Bessel0" (upper case!)
#ifdef DEVEL
              if (echo > 0) printf("# use a spherical Bessel solver for the Poisson equation\n");
              auto const st = Bessel_Poisson_solver(Ves, rho, g, Bessel_center, echo);
              if (st) warn("Bessel Poisson solver failed with status=%i", int(st));
              stat += st;
#else  // DEVEL
              error("electrostatic.solver=%c (Bessel) only available in the development branch!", es_solver_name); 
#endif // DEVEL
          } else
          if ('l' == (method | 32)) { // "load"
#ifdef DEVEL
              auto const Ves_in_filename = control::get("electrostatic.potential.from.file", "v_es.dat");
              auto const nerrors = debug_tools::read_from_file(Ves, Ves_in_filename, g.all(), 1, 1, "electrostatic potential", echo);
              if (nerrors) warn("electrostatic.solver=%c (load) from file %s had %d errors", es_solver_name, Ves_in_filename, nerrors); 
#else  // DEVEL
              error("electrostatic.solver=%c (load) only available in the development branch!", es_solver_name); 
#endif // DEVEL
          } else
          if ('n' == (method | 32)) { // "none"
              warn("electrostatic.solver=%c (none) may lead to unphysical results!", es_solver_name); 

          } else { // default
              if (echo > 2) printf("# electrostatic.solver=%c\n", es_solver_name);
              stat += iterative_poisson::solve(Ves, rho, g, method, echo);
          } // method

          if (echo > 1) print_stats(Ves, g.all(), g.dV(), "\n# electrostatic potential", eV);
#ifdef DEVEL
          if (echo > 6) {
              printf("# electrostatic potential at grid corners:");
              for(int iz = 0; iz < g[2]; iz += g[2] - 1) {
              for(int iy = 0; iy < g[1]; iy += g[1] - 1) {
              for(int ix = 0; ix < g[0]; ix += g[0] - 1) {
                  printf(" %g", Ves[(iz*g[1] + iy)*g[0] + ix]*eV);
              }}} // ix iy iz
              printf(" %s\n", _eV);
          } // echo
#endif // DEVEL
      } // scope

#ifdef DEVEL
      { // scope: export electrostatic potential to ASCII file
          auto const Ves_out_filename = control::get("electrostatic.potential.to.file", "");
          if (Ves_out_filename && *Ves_out_filename) {
              stat += debug_output::write_array_to_file(Ves_out_filename, Ves, 
                            g[0], g[1], g[2], echo, "electrostatic potential");
          }
      } // scope
#endif // DEVEL

      if (echo > 3) {
          double const Ees = 0.5*dot_product(g.all(), rho, Ves)*g.dV();
          printf("# inner product between rho_aug and Ves = %g %s\n", 2*Ees*eV,_eV);
      } // echo

      if (echo > 0) {
          if (control::get("poisson_solver.direct.projection", 0.) > 0) {
              printf("\n## all values of Ves in %s (unordered) as function of the distance to %s\n",
                                            _eV, (Bessel_center) ? "atom #0" : "the cell center");
              print_direct_projection(Ves, g, eV, Bessel_center);
          } // control
      } // echo

      return stat;
  } // solve


#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_solve(int const echo=3) {
      real_space::grid_t g(12, 14, 16);
      std::vector<double> Ves(g.all()), rho(g.all(), 0.0);
      auto const method = solver_method(control::get("electrostatic.solver", "multi-grid"));
      return solve(Ves.data(), rho.data(), g, method, echo);
  } // test_init

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_solve(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace poisson_solver
