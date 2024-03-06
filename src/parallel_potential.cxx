// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <cstdint> // int32_t, uint32_t
#include <cmath> // std::exp, ::sqrt, ::abs

#include "parallel_potential.hxx"

// #include "single_atom.h" // live_atom_##NAME##_ // use the libliveatom library
#include "single_atom.hxx" // ::atom_update

#include "status.hxx" // status_t
#include "control.hxx" // ::get
#include "real_space.hxx" // ::grid_t
#include "display_units.h" // Ang, _Ang
#include "self_consistency.hxx" // ::init_geometry_and_grid
#include "mpi_parallel.hxx" // ::init, ::finalize, MPI_COMM_WORLD, ::barrier
#include "parallel_poisson.hxx" // ::parallel_grid_t
#include "data_view.hxx" // view2D<>
#include "global_coordinates.hxx" // ::get
#include "inline_math.hxx" // pow2, set, add_product
#include "exchange_correlation.hxx" // ::LDA_kernel
#include "simple_stats.hxx" // ::Stats<>
#include "constants.hxx" // ::pi
#include "recorded_warnings.hxx" // warn
#include "data_list.hxx" // data_list<T>
#include "sho_tools.hxx" // ::nSHO
#include "solid_harmonics.hxx" // ::Y00

namespace parallel_potential {

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

    inline status_t live_atom_update(
        char const *const what    // selector string
      , int const natoms          // number of atoms
      , double  *const dp=nullptr // quantities (input/output) double  dp[natoms]
      , int32_t *const ip=nullptr // quantities (input/output) integer ip[natoms]
      , float   *const fp=nullptr // quantities (input)        float   fp[natoms or less]
      , double  *const *const dpp=nullptr // quantities (input/output) double* dpp[natoms]
    ) {
        auto const use = control::get("use.live.atom", 1.);
        return use ? single_atom::atom_update(what, natoms, dp, ip, fp, dpp) : 0;
        // warn("single_atom::atom_update deactivated", 0); return 0; // no effect
    } // live_atom_update

    inline double rho_Thomas_Fermi(double const Ekin, double & derivative) {
        if (Ekin <= 0) { derivative = 0; return 0; }
        double constexpr by_pi2 = 1./pow2(constants::pi);
        derivative = by_pi2 * std::sqrt(2*Ekin);
        double const value = derivative * Ekin * (2./3.);
        return value;
    } // rho_Thomas_Fermi


    status_t new_density_Thomas_Fermi(
          double rho_new[] // resulting density
        , double & Fermi_level // Fermi level
        , double const Vtot[] // input potential
        , size_t const nall // number of all grid points
        , MPI_Comm const comm // grid communicator
        , double const n_electrons // required total number of electrons 
        , double const dV // grid volume element
        , int const echo // log level
        , int const maxiter=99 // maximum number of iterations
        , float const threshold=1e-8 // convergence criterion
    ) {
        double mu = Fermi_level;
        { // scope: correct Fermi level
            simple_stats::Stats<double> s(0);
            for (size_t zyx = 0; zyx < nall; ++zyx) {
                s.add(Vtot[zyx]);
            } // zyx
            mpi_parallel::allreduce(s, comm);
            if (mu <= s.min()) { mu = 0.5*(s.max() + s.min()); }
            mu = std::max(s.min(), mu);
        } // scope

        if (echo > 3) std::printf("# Thomas-Fermi model, starting with chemical potential mu= %g %s\n", mu*eV, _eV);
        double ne{0}, de{0}; // ne: preliminary number of electrons, de: derivative w.r.t. the Fermi level
        int iterations{0};
        bool converged{false};
        do {
            de = 0; ne = 0;
            for (size_t zyx = 0; zyx < nall; ++zyx) {
                double der;
                double const rho = rho_Thomas_Fermi(mu - Vtot[zyx], der);
                ne += rho;
                de += der;
            } // zyx
            ne = mpi_parallel::sum(ne, comm) * dV;
            de = mpi_parallel::sum(de, comm) * dV;
            auto const residual = std::abs(ne - n_electrons);
            if (echo > 5) std::printf("# Thomas-Fermi iteration #%i mu= %g %s --> %g electrons, derivative= %g, residual= %.1e\n",
                                                                iterations, mu*eV, _eV, ne, de, residual);
            converged = (residual < threshold);
            mu += (n_electrons - ne)/std::max(0.1, de); // new mu by Newton step
            ++iterations;
        } while (!converged && iterations < maxiter);

        if (converged) {
            if (echo > 2) std::printf("# Thomas-Fermi converged in %d iterations, mu= %g %s --> %g electrons, derivative= %g\n",
                                                                iterations, mu*eV, _eV, ne, de);
            for (size_t zyx = 0; zyx < nall; ++zyx) {
                double der;
                rho_new[zyx] = rho_Thomas_Fermi(mu - Vtot[zyx], der);
            } // zyx
            Fermi_level = mu;
            return 0;
        } else { // converged
            warn("Thomas-Fermi model did not converge in %d iterations", iterations);
            return iterations;
        } // converged
    } // new_density_Thomas_Fermi


    template<typename real_t>
    void block_average(real_t v4[4*4*4], real_t const v8[8*8*8], simple_stats::Stats<double> & s) {
        for (int iz{0}; iz < 4; ++iz) {
        for (int iy{0}; iy < 4; ++iy) {
        for (int ix{0}; ix < 4; ++ix) {
            auto const sum =  v8[((iz*2    )*8 + iy*2    )*8 + ix*2    ]
                            + v8[((iz*2    )*8 + iy*2    )*8 + ix*2 + 1]
                            + v8[((iz*2    )*8 + iy*2 + 1)*8 + ix*2    ]
                            + v8[((iz*2    )*8 + iy*2 + 1)*8 + ix*2 + 1]
                            + v8[((iz*2 + 1)*8 + iy*2    )*8 + ix*2    ]
                            + v8[((iz*2 + 1)*8 + iy*2    )*8 + ix*2 + 1]
                            + v8[((iz*2 + 1)*8 + iy*2 + 1)*8 + ix*2    ]
                            + v8[((iz*2 + 1)*8 + iy*2 + 1)*8 + ix*2 + 1];
            auto const average = sum*0.125; // divide by 2*2*2
            v4[(iz*4 + iy)*4 + ix] = average;
            s.add(average);
        }}} // ix iy iz
    } // block_average


    template<typename real_t>
    double add_r2grid_to_block(
          real_t values[8*8*8]
        , uint32_t const block_coords[3]
        , double const h[3] // Cartesian grid spacings
        , double const atom_center[3]
        , double const r2coeff[] // radial function on an r^2-grid
        , double const factor=1
        , int const echo=0
        , float const r_cut=-1 // -1: auto
        , float const ar2=16.f // r^2-grid parameter
        , int32_t const nr2=4096 // r^2-grid size
    ) {
        assert(ar2 > 0);
        assert(ar2 == 16.f);
        double const r_max = std::sqrt((nr2 - 1.)/ar2); // largest radius of the r^2-grid
        double const rcut = (-1 == r_cut) ? r_max : std::min(double(r_cut), r_max);
        double const r2cut = pow2(rcut);
        assert(ar2*r2cut < nr2);
        double added_charge{0};
        size_t grid_points_inside{0};
        if (1) { // ToDo: check if any of the radii can be inside
            for (int iz{0}; iz < 8; ++iz) {
            for (int iy{0}; iy < 8; ++iy) {
            for (int ix{0}; ix < 8; ++ix) {
                auto const r2 = pow2(h[0]*(block_coords[0]*8 + ix + .5) - atom_center[0])
                              + pow2(h[1]*(block_coords[1]*8 + iy + .5) - atom_center[1])
                              + pow2(h[2]*(block_coords[2]*8 + iz + .5) - atom_center[2]);
                if (r2 < r2cut) {
                    int const ir2 = int(ar2*r2);

                    if (ir2 + 1 < nr2) {
                        double const w8 = ar2*r2 - ir2; // linear interpolation weight
                        auto const value_to_add = r2coeff[ir2] * (1 - w8) + r2coeff[ir2 + 1]*w8;

                    // this version is consistent with real_space.hxx
                    // if (ir2 < nr2) {
                    //     double const w8 = ar2*r2 - ir2; // linear interpolation weight
                    //     int const ir2p1 = ir2 + 1;
                    //     auto const value_to_add = (r2coeff[ir2] * (1 - w8)
                    //             + ((ir2p1 < nr2) ? r2coeff[ir2p1] : 0)*w8);

                        int const izyx = (iz*8 + iy)*8 + ix;
                        values[izyx] += factor*value_to_add;
                        added_charge += factor*value_to_add;
                    }
                    ++grid_points_inside;
                } // inside
            }}} // ix iy iz
        } // 1
        if (grid_points_inside && echo > 11) std::printf("# %lld grid points inside %g %s for grid block at [%d %d %d]\n",
            grid_points_inside, rcut*Ang, _Ang, block_coords[0], block_coords[1], block_coords[2]);
        return added_charge;
    } // add_r2grid_to_block


    double integrate_r2grid(double const r2coeff[], float const ar2=16.f, int32_t const nr2=4096) {
        double rho{0};
        auto const by_ar2 = 1./ar2;
        for (int ir2{0}; ir2 < nr2; ++ir2) {
            auto const r2 = ir2*by_ar2, r = std::sqrt(r2);
            rho += r2coeff[ir2]*r*.5*by_ar2;
        } // ir2
        return rho;
    } // integrate_r2grid


    status_t SCF(int const echo=0) {
        status_t stat(0);

        auto const comm = MPI_COMM_WORLD;
        auto const me = mpi_parallel::rank(comm);
        auto const nprocs = mpi_parallel::size(comm);

        // load geometry from +geometry.file=atoms.xyz (a good candidate is +geometry.file=test/geo/Cu40Zr22.xyz)
        real_space::grid_t g; // entire grid descriptor
        view2D<double> xyzZ;  // coordinates for all atoms
        int32_t n_all_atoms;  // number of all atoms
        stat += self_consistency::init_geometry_and_grid(g, xyzZ, n_all_atoms, 8, echo);

        // create a coarse grid descriptor
        real_space::grid_t gc(g[0]/2, g[1]/2, g[2]/2);
        {
            gc.set_boundary_conditions(g.boundary_conditions());
            gc.set_cell_shape(g.cell, echo);
            gc.set_grid_spacing(g.h[0]*2, g.h[1]*2, g.h[2]*2);
            for (int d{0}; d < 3; ++d) { assert(0 == (gc[d] & 0x3)); } // all grid numbers must be a multiple of 4
        }
        double const grid_center[] = {g[0]*g.h[0]*.5, g[1]*g.h[1]*.5, g[2]*g.h[2]*.5};

        // distribute the dense grid in 8x8x8 grid blocks to parallel owners
        parallel_poisson::parallel_grid_t pg(g, comm, 8, echo);

        // distribute the atom ownership TODO
        // simple distribution model: owner_rank == atom_id % nprocs, advantage: every rank can compute the owner
        // alternativ: make a weight grid pg.nb^3, init 0, add weight where atoms are, call load_balancer



        int32_t const na = n_all_atoms; // number of owned atoms          TODO atom ownership not parallelized TODO
        std::vector<float> ionization(na, 0.f);

        float take_atomic_valence_densities{1}; // 100% of the smooth spherical atomic valence densities is included in the smooth core densities

        std::vector<double> sigma_cmp(na, 1.); // spread of the Gaussian used in the compensation charges
        char const pawdata_from = (*control::get("pawdata.from", "auto")) | 32; // 'a': auto generate, 'f': pawxml_import
        std::vector<int32_t> numax(na, ('f' == pawdata_from)?-9:-1); // -1: LivePAW, -9: load from pawxml files
        std::vector<int32_t> lmax_qlm(na, -1);
        std::vector<int32_t> lmax_vlm(na, -1);

        // initialize and get sigma, lmax for each atom
        std::vector<int32_t> n_atom_rho(na, 0);
        data_list<double> atom_qlm, atom_vlm, atom_rho, atom_mat;
        {
            std::vector<double> Za(na, 0.);
            for (int32_t ia{0}; ia < na; ++ia) { Za[ia] = xyzZ(ia,3); }
            stat += live_atom_update("initialize", na, Za.data(), numax.data(), ionization.data(), (double**)1);
            stat += live_atom_update("lmax qlm",   na,    nullptr, lmax_qlm.data(), &take_atomic_valence_densities);
            stat += live_atom_update("lmax vlm",   na, (double*)1, lmax_vlm.data());
            stat += live_atom_update("sigma cmp",  na, sigma_cmp.data());

            view2D<int32_t> num(4, na, 0); // how many
            for (int32_t ia{0}; ia < na; ++ia) {
                num(0,ia) = pow2(1 + lmax_qlm[ia]); // qlm
                num(1,ia) = pow2(1 + lmax_vlm[ia]); // vlm
                int const ncoeff = sho_tools::nSHO(numax[ia]);
                num(2,ia) =   pow2(ncoeff);         // aDm
                num(3,ia) = 2*pow2(ncoeff);         // aHm+aSm
                n_atom_rho[ia] = num(2,ia); // number of coefficients for the atomic density matrices
            } // ia
            // get memory in the form of data_list containers
            atom_qlm = data_list<double>(na, num[0], 0.0);
            atom_vlm = data_list<double>(na, num[1], 0.0);
            atom_rho = data_list<double>(na, num[2], 0.0);
            atom_mat = data_list<double>(na, num[3], 0.0);
        } // scope

        // the r^2-grids are used to bring radial quantities to a Cartesian grid
        std::vector<int32_t> nr2(na, 1 << 12); // 4096
        std::vector<float>   ar2(na, 16.f); // with nr2 == 4096 rcut = 15.998 Bohr
        data_list<double> atom_vbar(nr2, 0.0); // zero potentials
        data_list<double> atom_rhoc(nr2, 0.0); // core_densities



        // allocate CPU memory
        view2D<double> augmented_density(pg.n_local(), 8*8*8, 0.0);
        view2D<double> V_electrostatic(pg.n_local(), 8*8*8, 0.0);


        double constexpr Y00 = solid_harmonics::Y00;
        double constexpr Y00sq = pow2(Y00); // = 1/(4*pi)
        
        { // scope: collect start densities from atomic densities
            // get smooth core densities on r^2-grids
            stat += live_atom_update("core densities", na, 0, nr2.data(), ar2.data(), atom_rhoc.data());

            if (take_atomic_valence_densities > 0) {
                auto & atom_rhov = atom_vbar; // use memory of vbar for the moment
                stat += live_atom_update("valence densities", na, 0, nr2.data(), ar2.data(), atom_rhov.data());
                for (int32_t ia{0}; ia < na; ++ia) {
                    // add valence density to r^2-gridded core density
                    if (echo > 15) std::printf("# rank#%i atom#%i wants to add %g core electrons\n", me, ia, integrate_r2grid(atom_rhoc[ia], ar2[ia], nr2[ia]));
                    add_product(atom_rhoc[ia], nr2[ia], atom_rhov[ia], take_atomic_valence_densities*1.);
                    if (echo > 15) std::printf("# rank#%i atom#%i wants to add %g core+valence electrons\n", me, ia, integrate_r2grid(atom_rhoc[ia], ar2[ia], nr2[ia]));
                } // ia
            } // take_atomic_valence_densities > 0

            auto const local_ids = pg.local_ids();
            double rho_total{0};
            for (int32_t ia{0}; ia < n_all_atoms; ++ia) { // TODO reduce loop to relevant periodic images
                double pos[3]; set(pos, 3, xyzZ[ia]);
                if (echo > 6) std::printf("# rank#%i adds start density for atom#%i at position %g %g %g %s\n", me, ia, pos[0]*Ang, pos[1]*Ang, pos[2]*Ang, _Ang);
                add_product(pos, 3, grid_center, 1.);
                double added_charge{0};
                for (uint32_t ilb{0}; ilb < pg.n_local(); ++ilb) { // local blocks
                    uint32_t ixyz[3]; global_coordinates::get(ixyz, local_ids[ilb]);
                    added_charge += add_r2grid_to_block(augmented_density[ilb], ixyz, g.grid_spacings(), pos, atom_rhoc[ia], Y00sq, echo);
                } // ilb
                if (echo > 0) std::printf("# rank#%i atom#%i %g electrons added\n", me, ia, added_charge*g.dV());
                rho_total += added_charge;
            } // ia
            rho_total = mpi_parallel::sum(rho_total, comm)*g.dV();
            if (echo > 3) std::printf("# added %g electrons as start density\n", rho_total);
        } // scope


        // configure the Poisson solver
        auto  const es_method    = control::get("poisson.method", "cg"); // {cg, sd}
        int   const es_echo      = control::get("poisson.echo", echo*1.);
        float const es_threshold = control::get("poisson.threshold", 3e-8);
        int   const es_maxiter   = control::get("poisson.maxiter", 200.);
        int   const es_miniter   = control::get("poisson.miniter", 3.);
        int   const es_restart   = control::get("poisson.restart", 4096.);

        int   const scf_maxiter   = control::get("scf.maxiter", 1.);

        double E_Fermi{0}; // Fermi level

        double nve{0}; // non-const total number of valence electrons
        { // scope: determine the number of valence electrons
            auto const keyword_valence_electrons = "valence.electrons";
            if ('a' == (*control::get(keyword_valence_electrons, "auto") | 32)) {
                std::vector<double> n_electrons_a(na, 0.); // number of valence electrons added by each atom
                stat += live_atom_update("#valence electrons", na, n_electrons_a.data());
                nve = std::accumulate(n_electrons_a.begin(), n_electrons_a.end(), 0.0);
                if (echo > 0) std::printf("\n# %s=auto --> %g valence electrons\n\n", keyword_valence_electrons, nve);
                control::set(keyword_valence_electrons, nve, control::echo_set_without_warning);
            } else {
                nve = control::get(keyword_valence_electrons, 0.0);
                if (echo > 0) std::printf("\n# %s=%g\n\n", keyword_valence_electrons, nve);
            } // auto
        } // scope
        double const n_valence_electrons = nve; // do not use nve beyond this point

        int scf_iteration{0};
        bool scf_run{true};
        do { // self-consistency loop

            mpi_parallel::barrier(comm);
            ++scf_iteration;

            { // scope: show statistics of augmented_density
                simple_stats::Stats<double> s(0);
                for (uint32_t ilb{0}; ilb < pg.n_local(); ++ilb) { // local blocks
                    auto const *const v8 = augmented_density[ilb];
                    for (int i9{0}; i9 < 512; ++i9) { s.add(v8[i9]); }
                } // ilb
                mpi_parallel::allreduce(s, comm);
                if (echo > 1) std::printf("# smooth augmented_density min %g max %g avg %g, %g electrons\n", s.min(), s.max(), s.mean(), s.sum()*g.dV());
            } // scope

            view2D<double> V_xc(pg.n_local(), 8*8*8, 0.0);
            { // scope: eval the XC potential and energy
                double E_xc{0}, E_dc{0};
                auto const *const density = augmented_density[0];
                auto       *const potential = V_xc[0];
                double rho_max{0}; int64_t i_max{-1};
                for (size_t i = 0; i < pg.n_local()*size_t(8*8*8); ++i) {
                    auto const rho_i = density[i];
                    if (rho_i > rho_max) { rho_max = rho_i; i_max = i; }
                    double vxc_i;
                    auto const exc_i = exchange_correlation::LDA_kernel(rho_i, vxc_i);
                    E_xc += rho_i*exc_i;
                    E_dc += rho_i*vxc_i; // double counting correction
                    potential[i] = vxc_i;
                } // i
                if (echo > 2) std::printf("# rank#%i rho_max= %g a.u. at index= %lli\n", me, rho_max, i_max);
                E_xc = mpi_parallel::sum(E_xc, comm) * g.dV(); // scale with volume element
                E_dc = mpi_parallel::sum(E_dc, comm) * g.dV(); // scale with volume element
                if (echo > 2) std::printf("# exchange-correlation energy on grid %.9f %s, double counting %.9f %s\n", E_xc*eV,_eV, E_dc*eV,_eV);
                // grid_xc_energy = E_xc;
            } // scope


            return 0; // early return (DEBUG)


            // ToDo: add compensation charges



            { // scope: Poisson equation
                if (echo > 3) std::printf("#\n# Solve the Poisson equation iteratively with %d ranks\n#\n", nprocs);
                float es_residual{0};
                auto const es_stat = parallel_poisson::solve(V_electrostatic[0], augmented_density[0],
                                                pg, *es_method, es_echo, es_threshold,
                                                &es_residual, es_maxiter, es_miniter, es_restart);
                if (echo > 2) std::printf("# Poisson equation %s, residual= %.2e a.u.\n#\n", es_stat?"failed":"converged", es_residual);
                stat += es_stat;
            } // scope

            { // scope: show statistics of V_electrostatic
                simple_stats::Stats<double> s(0);
                for (uint32_t ilb{0}; ilb < pg.n_local(); ++ilb) { // local blocks
                    auto const *const v8 = V_electrostatic[ilb];
                    for (int i9{0}; i9 < 512; ++i9) { s.add(v8[i9]); }
                } // ilb
                mpi_parallel::allreduce(s, comm);
                if (echo > 1) std::printf("# smooth electrostatic potential min %g max %g avg %g %s\n", s.min()*eV, s.max()*eV, s.mean()*eV, _eV);
            } // scope


            // ToDo: project electrostatic potential


            view2D<double> V_effective(pg.n_local(), 8*8*8, 0.0);
            view2D<double> V_coarse(pg.n_local(), 4*4*4, 0.0);
            { // scope: add potentials and reduce them to 4x4x4 grid points per block
                simple_stats::Stats<double> s(0);
                for (uint32_t ilb{0}; ilb < pg.n_local(); ++ilb) { // local blocks
                    set(V_effective[ilb], 512, V_electrostatic[ilb]);
                    add_product(V_effective[ilb], 512, V_xc[ilb], 1.);
                    // reduce from blocks of 8x8x8 to 4x4x4
                    block_average(V_coarse[ilb], V_effective[ilb], s);
                } // ilb
                mpi_parallel::allreduce(s, comm);
                if (echo > 1) std::printf("# smooth effective potential min %g max %g avg %g %s\n", s.min()*eV, s.max()*eV, s.mean()*eV, _eV);
            } // scope

            // ToDo: call energy-contour integration to find a new density
            view2D<double> new_density(pg.n_local(), 8*8*8, 0.0);
            { // scope: apply Thomas-Fermi approximation
                auto const stat_TF = new_density_Thomas_Fermi(new_density[0], E_Fermi, V_effective[0],
                                pg.n_local()*size_t(512), comm, n_valence_electrons, g.dV(), echo);
                stat += stat_TF;
                simple_stats::Stats<double> s(0);
                for (uint32_t ilb{0}; ilb < pg.n_local(); ++ilb) { // local blocks
                    for (int i9{0}; i9 < 512; ++i9) { s.add(new_density(ilb,i9)); }
                } // ilb
                mpi_parallel::allreduce(s, comm);
                if (echo > 1) std::printf("# new density after Thomas-Fermi min %g max %g avg %g a.u. sum %g electrons\n", s.min(), s.max(), s.mean(), s.sum()*g.dV());
            } // scope


            scf_run = (scf_iteration < scf_maxiter);

        } while (scf_run); // self-consistency loop 

        return stat;
    } // SCF

    status_t test_r2grid_integrator(int const echo=0, int const nr2=4096, float const ar2=16) {
        status_t stat(0);
        std::vector<double> vec(nr2, 1.0);
        { // integrate 1.0 with r^2*dr up to the outer max radius R
            auto const R = std::sqrt((nr2 - 1.)/ar2);
            auto const vol = integrate_r2grid(vec.data(), ar2, nr2), ref = pow3(R)/3, dev = vol - ref;
            if (echo > 3) std::printf("# %s: found %g expect %g dev= %.3f %%\n", __func__, vol, ref, dev/(ref*.001));
            stat += (std::abs(dev) < ref*2e-4);
        }
        { // prepare a Gaussian
            for (int ir2{0}; ir2 < nr2; ++ir2) { vec[ir2] = std::exp(-ir2*(.125/ar2)); }
            auto const vol = integrate_r2grid(vec.data(), ar2, nr2), ref = std::sqrt(constants::pi*32), dev = vol - ref;
            if (echo > 3) std::printf("# %s: found %g expect %g dev= %.3f %%\n", __func__, vol, ref, dev/(ref*.001));
            stat += (std::abs(dev) < ref*2e-4);
        }
        return stat;
    } // test_r2grid_integrator

    status_t test_scf(int const echo=0) {
        status_t stat(0);
        stat += SCF(echo);
        return stat;
    } // test_scf

    status_t all_tests(int const echo) {
        status_t stat(0);
        stat += test_r2grid_integrator(echo);
        auto const already_initialized = mpi_parallel::init();
        stat += test_scf(echo);
        if (!already_initialized) mpi_parallel::finalize();
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace parallel_potential
