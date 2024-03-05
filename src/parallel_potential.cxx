// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert
// #include <vector> // std::vector<T>
#include <cstdint> // int32_t, uint32_t
#include <cmath> // std::exp

#include "parallel_potential.hxx"

// #include "single_atom.h" // live_atom_##NAME##_ // use the libliveatom library

#include "status.hxx" // status_t
#include "control.hxx" // ::get
#include "real_space.hxx" // ::grid_t
#include "display_units.h" // Ang, _Ang
#include "self_consistency.hxx" // ::init_geometry_and_grid
#include "mpi_parallel.hxx" // ::init, ::finalize, MPI_COMM_WORLD
#include "parallel_poisson.hxx" // ::parallel_grid_t
#include "data_view.hxx" // view2D<>
#include "global_coordinates.hxx" // ::get
#include "inline_math.hxx" // pow2


namespace parallel_potential {

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS


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

        // distribute the dense grid in 8x8x8 grid blocks to parallel owners
        parallel_poisson::parallel_grid_t pg(g, comm, 8, echo);

        // distribute the atom ownership TODO
        // find a suitable model for the workload of libliveatom tasks


        // allocate CPU memory
        view2D<double> augmented_density(pg.n_local(), 512, 0.0);
        view2D<double> electrostatic_potential(pg.n_local(), 512, 0.0);

        // configure the Poisson solver
        auto  const es_method    = control::get("poisson.method", "cg"); // {cg, sd}
        int   const es_echo      = control::get("poisson.echo", echo*1.);
        float const es_threshold = control::get("poisson.threshold", 3e-8);
        int   const es_maxiter   = control::get("poisson.maxiter", 200.);
        int   const es_miniter   = control::get("poisson.miniter", 3.);
        int   const es_restart   = control::get("poisson.restart", 4096.);

        
        { // scope: collect start densities from atomic densities
            // auto const nb = pg.grid_blocks();
            auto const local_ids = pg.local_ids();
            double rho_total{0};
            for (int32_t ia{0}; ia < n_all_atoms; ++ia) { // TODO reduce loop to relevant periodic images
                auto const *const pos = xyzZ[ia];
                if (echo > 15) std::printf("# rank#%i adds start density for atom#%i at position %g %g %g %s\n", me, ia, pos[0]*Ang, pos[1]*Ang, pos[2]*Ang, _Ang);
                for (uint32_t ilb{0}; ilb < pg.n_local(); ++ilb) {
                    uint32_t ixyz[3]; global_coordinates::get(ixyz, local_ids[ilb]);
                    for (int iz{0}; iz < 8; ++iz) {
                    for (int iy{0}; iy < 8; ++iy) {
                    for (int ix{0}; ix < 8; ++ix) {
                        auto const r2 = pow2(g.h[0]*(ixyz[0]*8 + ix + .5 - pos[0]))
                                      + pow2(g.h[1]*(ixyz[1]*8 + iy + .5 - pos[1]))
                                      + pow2(g.h[2]*(ixyz[2]*8 + iz + .5 - pos[2]));
                        if (r2 < 36) {
                            auto const add_rho = std::exp(-r2);
                            rho_total += add_rho;
                            augmented_density(ilb,(iz*8 + iy)*8 + ix) += add_rho;
                        } // r2 < 36
                    }}} // ix iy iz
                } // ilb
            } // ia
            mpi_parallel::sum(&rho_total, 1, comm);
            if (echo > 3) std::printf("# added %g electrons as start density\n", rho_total);
        } // scope

        { // scope: Poisson equation
            float es_residual{0};
            if (echo > 3) std::printf("#\n# Solve the Poisson equation iteratively with %d ranks\n#\n", nprocs);
            auto const es_stat = parallel_poisson::solve(electrostatic_potential[0],
                                            augmented_density[0], pg, *es_method, es_echo, es_threshold,
                                            &es_residual, es_maxiter, es_miniter, es_restart);
            mpi_parallel::max(&es_residual, 1, comm);
            if (echo > 2) std::printf("# Poisson equation %s, residual= %.2e a.u.\n#\n", (es_stat)?"failed":"converged", es_residual);
            stat += es_stat;
        } // scope

        return stat;
    } // SCF

    status_t test_scf(int const echo=0) {
        status_t stat(0);
        stat += SCF(echo);
        return stat;
    } // test_scf

    status_t all_tests(int const echo) {
        status_t stat(0);
        auto const already_initialized = mpi_parallel::init();
        stat += test_scf(echo);
        if (!already_initialized) mpi_parallel::finalize();
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace parallel_potential
