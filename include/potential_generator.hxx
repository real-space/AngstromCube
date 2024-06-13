#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf, ::snprintf

#include "status.hxx" // status_t
#include "control.hxx" // ::get

#include "real_space.hxx" // ::grid_t
#include "bessel_transform.hxx" // ::transform_s_function
#include "radial_grid.h" // radial_grid_t
#include "data_view.hxx" // view2D<T>
#include "debug_output.hxx" // dump_to_file
#include "solid_harmonics.hxx" // ::rlXlm, ::lm_index
#include "inline_math.hxx" // set, dot_product

#ifdef    DEVEL
  #include "finite_difference.hxx" // ::stencil_t
  #include "single_atom.hxx" // ::atom_update
  #include "print_tools.hxx" // print_stats
  #include "radial_r2grid.hxx" // ::r_axis
  #include "lossful_compression.hxx" // print_compressed
  #include "simple_timer.hxx" // SimpleTimer
  #include "poisson_solver.hxx" // ::print_direct_projection
#endif // DEVEL

namespace potential_generator {

  // debugging tool wrapper
  template <typename real_t>
  status_t write_array_to_file(
        char const *filename // file name to write to
      , real_t const array[]  // array pointer
      , int const nx, int const ny, int const nz // grid dimensions
      , int const echo=0 // log-level
      , char const *arrayname="" // some description to appear in the file
  ) {
      char title[128]; std::snprintf(title, 128, "%i x %i x %i  %s", nz, ny, nx, arrayname);
      auto const size = size_t(nz) * size_t(ny) * size_t(nx);
      return dump_to_file(filename, size, array, nullptr, 1, 1, title, echo);
  } // write_array_to_file


  template <typename real_t>
  status_t add_smooth_quantities(
        real_t values[] // add to this function on a 3D grid
      , real_space::grid_t const & g // real-space grid descriptor
      , int const na // number of atoms (radial grid centers)
      , int32_t const nr2[] // number of r^2-grid points
      , float const ar2[] // r2-grid density
      , view2D<double> const & center // [natoms][4] center coordinates
      , int const n_periodic_images
      , view2D<double> const & periodic_images
      , double const *const *const atom_qnt // atom data on r2-grids
      , int const echo=0 // log-level
      , int const echo_q=0 // log-level for the charge
      , double const factor=1 // multipliyer
      , char const *quantity="???" // description for log-messages
  ) {
      // add contributions from smooth core densities

      status_t stat(0);
      for (int ia = 0; ia < na; ++ia) {
#ifdef    DEVEL
          if (echo > 11) {
              std::printf("\n## r, %s of atom #%i\n", quantity, ia);
              print_compressed(radial_r2grid::r_axis(nr2[ia], ar2[ia]).data(), atom_qnt[ia], nr2[ia]);
          } // echo
#endif // DEVEL
          double q_added{0};
          for (int ii = 0; ii < n_periodic_images; ++ii) {
              double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
              double q_added_image{0};
              stat += real_space::add_function(values, g, atom_qnt[ia], nr2[ia], ar2[ia], &q_added_image, cnt, factor);
              if (echo_q > 11) std::printf("# %g electrons %s of atom #%d added for image #%i\n", q_added_image, quantity, ia, ii);
              q_added += q_added_image;
          } // periodic images
#ifdef    DEVEL
          if (echo_q > 0) {
              std::printf("# after adding %g electrons %s of atom #%d:", q_added, quantity, ia);
              print_stats(values, g.all(), g.dV());
          } // echo
          if (echo_q > 3) std::printf("# added %s for atom #%d is  %g electrons\n", quantity, ia, q_added);
//        if (echo_q > 3) std::printf("#    00 compensator charge for atom #%d is %g electrons\n", ia, atom_qlm[ia][00]*Y00inv);
#endif // DEVEL
      } // ia
      return stat;
  } // add_smooth_quantities


  template <typename real_t, int debug=0>
  status_t add_generalized_Gaussian(
        real_t values[] // grid values which are modified
      , real_space::grid_t const & g // grid descriptor
      , double const coeff[] // coefficients of the compensators
      , int const ellmax=-1
      , double const center[3]=nullptr // spherical center w.r.t. the position of grid point (0,0,0)
      , double const sigma=1 // spread of the Gaussian exp(-r^2/(2*sigma^2))
      , int const echo=0
      , double *added=nullptr // optional result: how much (e.g. charge) was added
      , double const factor=1 // optional scaling
      , float const r_cut=-1 // radial truncation, -1:use the max radius of the r^2-grid
  ) {
      // Add a generalized Gaussian multipole compensator charge
      status_t stat(0);
      if (ellmax < 0) return stat; // early return
      assert(sigma > 0);
      double c[3] = {0,0,0}; if (center) set(c, 3, center);
      double const r_max = 9*sigma; // default truncation radius
      double const rcut = (-1 == r_cut) ? r_max : std::min(double(r_cut), r_max);
      double const r2cut = rcut*rcut;

      int const nlm = pow2(1 + ellmax);
      double const sigma2inv = 0.5/pow2(sigma);
      std::vector<double> scaled_coefficients(nlm, 0.0);
      std::vector<double> scale_factor(debug*nlm, 0.0);
      { // scope to scale coefficients
          double radial_norm{constants::sqrt2/(constants::sqrtpi*sigma)};
          for (int ell = 0; ell <= ellmax; ++ell) {
              radial_norm /= (sigma*sigma*(2*ell + 1));
              for (int emm = -ell; emm <= ell; ++emm) {
                  int const ilm = solid_harmonics::lm_index(ell, emm);
                  scaled_coefficients[ilm] = coeff[ilm] * radial_norm;
                  if (debug) scale_factor[ilm] = radial_norm*g.dV();
              } // emm
          } // ell
      } // scope
      view2D<double> overlap_ij(debug*nlm, nlm, 0.0);
      std::vector<double> rlXlm(nlm); // must be thread-private in case of OMP parallelism
      double added_charge{0}; // clear
      size_t modified{0};
      size_t nwindow{1};
      if (g.is_Cartesian()) {
          // define a window in which we add the radial function
          int imn[3], imx[3];
          for (int d = 0; d < 3; ++d) {
              imn[d] = std::max(0, int(std::floor((c[d] - rcut)*g.inv_h[d])));
              imx[d] = std::min(   int(std::ceil ((c[d] + rcut)*g.inv_h[d])), g[d] - 1);
              if (echo > 8) std::printf("# %s window %c = %d elements from %d to %d\n", __func__, 'x'+d, imx[d] + 1 - imn[d], imn[d], imx[d]);
              nwindow *= std::max(0, imx[d] + 1 - imn[d]);
          } // d

          for (        int iz = imn[2]; iz <= imx[2]; ++iz) {  double const vz = iz*g.h[2] - c[2], vz2 = vz*vz;
              for (    int iy = imn[1]; iy <= imx[1]; ++iy) {  double const vy = iy*g.h[1] - c[1], vy2 = vy*vy;
                  if (vz2 + vy2 < r2cut) {
                  for (int ix = imn[0]; ix <= imx[0]; ++ix) {  double const vx = ix*g.h[0] - c[0], vx2 = vx*vx;
                      double const r2 = vx2 + vy2 + vz2;
                      if (r2 < r2cut) {
                          int const izyx = (iz*g('y') + iy)*g('x') + ix;
                          solid_harmonics::rlXlm(rlXlm.data(), ellmax, vx, vy, vz);
                          double const Gaussian = std::exp(-sigma2inv*r2);
                          double const value_to_add = Gaussian * dot_product(nlm, scaled_coefficients.data(), rlXlm.data());
                          values[izyx] += factor*value_to_add;
                          added_charge += factor*value_to_add;
                          ++modified;
#if 0
                          if (echo > 13) {
                              std::printf("#rs %g %g\n", std::sqrt(r2), value_to_add);
                              // std::printf("#rs %.1f %.1f %.1f %.12f\n", vx*g.inv_h[0], vy*g.inv_h[1], vz*g.inv_h[2], value_to_add);
                          } // echo
#endif // 0
                          for (int ilm = 0; ilm < nlm*debug; ++ilm) {
                              add_product(overlap_ij[ilm], nlm, rlXlm.data(), Gaussian*rlXlm[ilm]);
                          } // ilm
                      } // inside rcut
                  } // ix
                  } // y-z inside rcut
              } // iy
          } // iz
      } else { // is_Cartesian
          // general cell
          double const denom[] = {1./g[0], 1./g[1], 1./g[2]}; // grid denominators
          for (        int iz = 0; iz < g[2]; ++iz) {
              for (    int iy = 0; iy < g[1]; ++iy) {
                  for (int ix = 0; ix < g[0]; ++ix) {
                      double const iv[] = {ix*denom[0], iy*denom[1], iz*denom[2]};
                      double rv[3];
                      for (int d = 0; d < 3; ++d) {
                          rv[d] = iv[0]*g.cell[0][d] + iv[1]*g.cell[1][d] + iv[2]*g.cell[2][d] - c[d];
                      } // d
                      double const r2 = pow2(rv[0]) + pow2(rv[1]) + pow2(rv[2]);
                      if (r2 < r2cut) {
                          int const izyx = (iz*g('y') + iy)*g('x') + ix;
                          solid_harmonics::rlXlm(rlXlm.data(), ellmax, rv[0], rv[1], rv[2]);
                          double const Gaussian = std::exp(-sigma2inv*r2);
                          double const value_to_add = Gaussian * dot_product(nlm, scaled_coefficients.data(), rlXlm.data());
                          values[izyx] += factor*value_to_add;
                          added_charge += factor*value_to_add;
                          ++modified;
                          for (int ilm = 0; ilm < nlm*debug; ++ilm) {
                              add_product(overlap_ij[ilm], nlm, rlXlm.data(), Gaussian*rlXlm[ilm]);
                          } // ilm
                      } // inside rcut
                  } // ix
              } // iy
          } // iz
          nwindow = g.all();
      } // is_Cartesian
      added_charge *= g.dV(); // volume integral
      if (echo > 7) std::printf("# %s modified %.3f k inside a window of %.3f k on a grid of %.3f k grid values, added charge= %g\n",
                                   __func__, modified*1e-3, nwindow*1e-3, g('x')*g('y')*g('z')*1e-3, added_charge); // show stats
      if (debug) {
          double dev[] = {0, 0}; // measure the deviation of overlap_ij from unity
          for (int ilm = 0; ilm < nlm*debug; ++ilm) {
              if (echo > 5) std::printf("# ilm=%3d ", ilm);
              for (int jlm = 0; jlm < nlm; ++jlm) {
                  double const ovl = overlap_ij(ilm,jlm)*scale_factor[ilm];
                  int const diag = (ilm == jlm);
                  if (echo > 5) std::printf(diag?" %.9f":" %9.1e", ovl);
                  dev[diag] = std::max(dev[diag], std::abs(ovl - diag));
              } // jlm
              if (echo > 5) std::printf("\n");
          } // ilm
          if (echo > 1) std::printf("# %s(ellmax=%d) is orthogonal to %.1e, normalized to %.1e\n", __func__, ellmax, dev[0], dev[1]);
      } // debug

      if (added) *added = added_charge;
      return stat;
  } // add_generalized_Gaussian



#ifdef    DEVEL

  inline status_t potential_projections(
        real_space::grid_t const & g // dense grid descriptor
      , double const Ves[] // electrostatic potential on g
      , double const Vxc[] // exchange-correlation potential on g
      , double const Vtot[] // total effective potential on g
      , double const rho[] // electron density on g
      , double const cmp[] // compensation charge density on g
      , int const na=0 // number of atoms
      , view2D<double> const *const center_ptr=nullptr
      , float const rcut=32
      , int const echo=0 // log-level
  ) {
      status_t stat(0);
      if (!g.is_Cartesian()) {
          if (echo > 1) std::printf("# %s: only implemented for Cartesian grids\n", __func__);
          return stat;
      } // is_Cartesian

      std::vector<double> Laplace_Ves(g.all(), 0.0);

      int const verify_Poisson = control::get("potential_generator.verify.poisson", 0.);
      for (int nfd = 1; nfd < verify_Poisson; ++nfd) {
          finite_difference::stencil_t<double> const fd(g.h, nfd, -.25/constants::pi);
          { // scope: compute the Laplacian using high-order finite-differences
              SimpleTimer timer(__FILE__, __LINE__, "finite-difference", echo);
              stat += finite_difference::apply(Laplace_Ves.data(), Ves, g, fd);
              { // scope: Laplace_Ves should match rho
                  double res_a{0}, res_2{0};
                  for (size_t i = 0; i < g.all(); ++i) {
                      res_a += std::abs(Laplace_Ves[i] - rho[i]);
                      res_2 +=     pow2(Laplace_Ves[i] - rho[i]);
                  } // i
                  res_a *= g.dV(); res_2 = std::sqrt(res_2*g.dV());
                  if (echo > 1) std::printf("# Laplace*Ves - rho: residuals abs %.2e rms %.2e (FD-order=%i)\n", res_a, res_2, nfd);
              } // scope
          } // timer
      } // nfd

      int const use_Bessel_projection = control::get("potential_generator.use.bessel.projection", 0.);
      if (use_Bessel_projection)
      { // scope: use a Bessel projection around each atom position to compare 3D and radial quantities

          double constexpr Y00sq = pow2(solid_harmonics::Y00);

          view2D<double> periodic_images;
          int const n_periodic_images = boundary_condition::periodic_images(periodic_images,
                                          g.cell, g.boundary_conditions(), rcut, echo - 4);

          std::vector<radial_grid_t const*> rg(na, nullptr); // pointers to smooth radial grid descriptors
          {   // break the interface to get the radial grid descriptors
              auto const dcpp = reinterpret_cast<double const *const *>(rg.data());
              auto const dpp  =       const_cast<double       *const *>(dcpp);
              stat += single_atom::atom_update("radial grids", na, 0, 0, 0, dpp);
          }

          double const* const value_pointers[] = {Ves, Vxc, Vtot, rho, cmp, Laplace_Ves.data()};
          char const *  array_names[] = {"Ves", "Vxc", "Vtot", "rho", "cmp", "LVes"};
          //     Ves; // analyze the electrostatic potential
          //     Vxc; // analyze the xc potential
          //     Vtot; // analyze the total potential: Vxc + Ves
          //     rho; // analyze the augmented density
          //     cmp; // analyze only the compensator density
          //     Laplace_Ves; // analyze the augmented density computed as Laplacian*Ves

          for (int iptr = 0; iptr < std::min(6, use_Bessel_projection); ++iptr) {
              // SimpleTimer timer(__FILE__, __LINE__, "Bessel-projection-analysis", echo);
              auto const values = value_pointers[iptr];
              auto const array_name = array_names[iptr];

              // report extremal values of what is stored on the grid
              if (echo > 1) { std::printf("\n# real-space stats of %s:", array_name); print_stats(values, g.all(), g.dV()); }

              for (int ia = 0; ia < na; ++ia) {
                  float const dq = 1.f/16;
                  int const nq = int(constants::pi/(g.smallest_grid_spacing()*dq));
                  std::vector<double> qc(nq, 0.0);

                  { // scope: Bessel core
                      std::vector<double> qc_image(nq, 0.0);
                      if (nullptr != center_ptr) { auto const & center = *center_ptr;
                          for (int ii = 0; ii < n_periodic_images; ++ii) {
                              double cnt[3]; set(cnt, 3, center[ia]); add_product(cnt, 3, periodic_images[ii], 1.0);
                              stat += real_space::Bessel_projection(qc_image.data(), nq, dq, values, g, cnt);
                              add_product(qc.data(), nq, qc_image.data(), 1.0);
                          } // ii
                      } // center
                  } // scope

                  scale(qc.data(), nq, Y00sq);

                  std::vector<double> qcq2(nq, 0.0);
                  for (int iq = 1; iq < nq; ++iq) { // start from 1 to avoid the q=0 term
                      qcq2[iq] = 4*constants::pi*qc[iq]/pow2(iq*dq); // cheap Poisson solver in Bessel transform
                  } // iq

                  if (echo > 11) {
                      std::printf("\n# Bessel coeff of %s for atom #%d:\n", array_name, ia);
                      for (int iq = 0; iq < nq; ++iq) {
                          std::printf("# %g %g %g\n", iq*dq, qc[iq], qcq2[iq]);
                      } // iq
                      std::printf("\n\n");
                  } // echo

                  if (echo > 3) {
                      std::vector<double> rs(rg[ia]->n);
                      bessel_transform::transform_s_function(rs.data(), qc.data(), *rg[ia], nq, dq, true); // transform back to real-space again
                      std::printf("\n## Real-space projection of %s for atom #%d:\n", array_names[iptr], ia);
                      float const compression_threshold = 1e-4;
                      print_compressed(rg[ia]->r, rs.data(), rg[ia]->n, compression_threshold);

                      if ((values == rho) || (values == Laplace_Ves.data())) {
                          bessel_transform::transform_s_function(rs.data(), qcq2.data(), *rg[ia], nq, dq, true); // transform electrostatic solution to real-space
                          std::printf("\n## Electrostatics computed by Bessel transform of %s for atom #%d:\n", array_names[iptr], ia);
                          print_compressed(rg[ia]->r, rs.data(), rg[ia]->n, compression_threshold);
                      } // density
                  } // echo
              } // ia

          } // iptr loop for different quantities represented on the grid.

      } // scope Bessel

      if (echo > 1) {
          if (nullptr != center_ptr) { auto const & center = *center_ptr;
              if (control::get("potential_generator.direct.projection", 0.) > 0) {
                  std::printf("\n## all values of Vtot in %s (unordered) as function of the distance to %s\n",
                                                    _eV, (na > 0) ? "atom #0" : "the cell center");
                  poisson_solver::print_direct_projection(Vtot, g, eV, (na > 0) ? center[0] : nullptr);
              } // control
          } else warn("no coordinates passed for na=%d atoms, center_ptr==nullptr", na);
      } // echo

      { // scope: export total potential to ASCII file
          auto const Vtot_out_filename = control::get("total.potential.to.file", "vtot.dat");
          if (*Vtot_out_filename) stat += write_array_to_file(Vtot_out_filename, Vtot, g[0], g[1], g[2], echo);
      } // scope

      return stat;
  } // potential_projections

#endif // DEVEL

  status_t all_tests(int const echo=0); // declaration only

} // namespace potential_generator
