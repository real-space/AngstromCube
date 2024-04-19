// This file is part of AngstromCube under MIT License

#include <cstdint> // int8_t
#include <cassert> // assert
#include <cstdio> // std::printf, ::snprintf
#include <vector> // std::vector<T>
#include <complex> // std::complex<real_t>

#include "energy_contour.hxx"

#include "control.hxx" // ::get
#include "display_units.h" // eV, _eV, Ang, _Ang, Kelvin, _Kelvin
#include "mpi_parallel.hxx" // MPI_COMM_WORLD
#include "data_view.hxx" // view2D<T>
#include "parallel_poisson.hxx" // ::parallel_grid_t
#include "action_plan.hxx" // action_plan_t
#include "green_function.hxx" // ::construct_Green_function, ::update_atom_matrices, ::update_phases, ::update_energy_parameter, ::update_potential
#include "real_space.hxx" // ::grid_t
#include "data_list.hxx" // data_list<T>
#include "recorded_warnings.hxx" // warn
#include "inline_math.hxx" // set, add_product
#include "green_solver.hxx" // green_solver_t
#include "recorded_warnings.hxx" // error, warn
#include "sho_tools.hxx" // ::nSHO
#include "sho_projection.hxx" // ::get_sho_prefactors
#include "brillouin_zone.hxx" // ::get_kpoint_mesh, ::WEIGHT

namespace energy_contour {

    Integrator::Integrator( // implementation of constructor
          real_space::grid_t const & gc // coarse grid descriptor
        , std::vector<double> const & xyzZinso // all atoms
        , int const echo // verbosity
    ) {
        if (echo > 0) std::printf("# construct %s with grid=[%d %d %d]\n", __func__, gc[0], gc[1], gc[2]);
        plan_ = new action_plan_t(); // CPU memory for the plan
        auto const stat = green_function::construct_Green_function(*plan_,
                            gc.grid_points(), gc.boundary_conditions(), gc.grid_spacings(),
                            xyzZinso, echo);
        if (stat) warn("construct_Green_function returned status= %i", int(stat));

        if (echo > 0) std::printf("# move green_solver_t\n");
        solver_ = new green_solver_t(plan_, echo);
        if (echo > 0) std::printf("# constructed %s\n", __func__);
    } // constructor


    status_t show_contour(std::vector<std::complex<double>> const & E_points, int const echo=0) {
        // an ASCII plot where the energy points are located
        if (echo < 1) return 0;
        auto const nE = E_points.size();
        if (nE < 1) return 0;
        float ex[2][2] = {{9e9f, -9e9f}, {9e9f, -9e9f}};
        for (auto const ep : E_points) {
            ex[0][0] = std::min(ex[0][0], float(ep.real()));
            ex[0][1] = std::max(ex[0][1], float(ep.real()));
            ex[1][0] = std::min(ex[1][0], float(ep.imag()));
            ex[1][1] = std::max(ex[1][1], float(ep.imag()));
        } // ep
        std::printf("#\n# energy contour within [%g, %g %s] and [%g, %g %s]\n",
            ex[0][0]*eV, ex[0][1]*eV, _eV, ex[1][0]*Kelvin, ex[1][1]*Kelvin, _Kelvin);
        return 0;
    } // show_contour


    std::vector<std::complex<double>> get_energy_mesh(int const echo=0) {
        auto     const eB = control::get("energy_contour.band.bottom", -1.);
        unsigned const nB = control::get("energy_contour.bottom",    5.);
        unsigned const nV = control::get("energy_contour.vertical", 15.);
        unsigned const nM = control::get("energy_contour.matsubara", 3.);
        auto    const kBT = control::get("energy_contour.temperature", 1e-2);
        auto const nE = nB + nV + nM;
        if (echo > 3) std::printf("# energy_contour with %d + %d + %d = %d points\n", nB, nV, nM, nE); 

        std::vector<std::complex<double>> Ep(nE);
        auto const dEb = kBT/std::max(1, int(nB));
        int jE{0};
        for (int iE{0}; iE < nB; ++iE) {
            Ep[jE] = std::complex<double>(eB, (iE + 1)*dEb);
            if (echo > 8) std::printf("# energy_mesh[%2i]=(%11.6f %s, %6.1f %s) bottom\n", jE, Ep[jE].real()*eV, _eV, Ep[jE].imag()*Kelvin, _Kelvin);
            ++jE;
        } // iE
        auto const dEv = eB/std::max(1, int(nV) - 1);
        for (int iE{0}; iE < nV; ++iE) {
            Ep[jE] = std::complex<double>((int(nV) - iE - 2)*dEv, kBT);
            if (echo > 8) std::printf("# energy_mesh[%2i]=(%11.6f %s, %6.1f %s) vertical\n", jE, Ep[jE].real()*eV, _eV, Ep[jE].imag()*Kelvin, _Kelvin);
            ++jE;
        } // iE
        auto const dEm = kBT/std::max(1, int(nM));
        for (int iE{0}; iE < nM; ++iE) {
            Ep[jE] = std::complex<double>(0.0, (nM - iE - 0.5)*dEm);
            if (echo > 8) std::printf("# energy_mesh[%2i]=(%11.6f %s, %6.1f %s) Matsubara\n", jE, Ep[jE].real()*eV, _eV, Ep[jE].imag()*Kelvin, _Kelvin);
            ++jE;
        } // iE
        if (echo > 6) std::printf("# energy_mesh %2i points\n", jE);
        assert(nE == jE);
        show_contour(Ep, echo);
        return Ep;
    } // get_energy_mesh


    status_t Integrator::integrate(
          double rho_888[] // resulting density in [nblocks][8*8*8] data layout
        , double & Fermi_level // Fermi level
        , double const Vtot[] // input potential in [nblocks][4*4*4]
        , data_list<double> const & atom_mat // atomic_Hamiltonian elements, only in atom owner ranks
        , std::vector<int32_t> const & numax_prj
        , std::vector<double> const & sigma_prj
        , parallel_poisson::parallel_grid_t const & pg
        , double const n_electrons // =1 // required total number of electrons 
        , double const dV // =1 // grid volume element
        , int const echo // =0 // log level
    ) {
        status_t stat(0);
        size_t constexpr n4x4x4 = 4*4*4;

        int const check = control::get("check", 0.);
        int const iterations = control::get("green_solver.iterations", 99.);
        if (echo > 0) std::printf("\n# energy_contour::integration(E_Fermi=%g %s, %g electrons, echo=%d) +check=%i\n", Fermi_level*eV, _eV, n_electrons, echo, check);

        auto const nblocks = pg.n_local();
        assert(nullptr != plan_);
        assert(nullptr != solver_);
        auto & plan = *plan_;
        plan.echo = echo >> 2; // lower verbosity

        if (plan.nCols != nblocks) warn("model assumes that each local block has one RHS, found n_local= %d and p.nRHS= %d", nblocks, plan.nCols);

        auto const nAtoms = atom_mat.nrows();
        std::vector<std::vector<double>> AtomMatrices(nAtoms);
        for (int iAtom{0}; iAtom < nAtoms; ++iAtom) {
            auto const nc2 = atom_mat.ncols(iAtom);
            auto const numax = numax_prj.at(iAtom);
            auto const nc = sho_tools::nSHO(numax);
            assert(nc2 == 2*nc*nc);
            auto const rescale = sho_projection::get_sho_prefactors(numax, sigma_prj.at(iAtom));
            AtomMatrices[iAtom] = std::vector<double>(nc2);
            for (int i{0}; i < nc; ++i) {
                auto const rescale_i = rescale.at(i);
                for (int j{0}; j < nc; ++j) {
                    // atom matrices need to be prepared for projection with unnormalized Gauss-Hermite functions
                    AtomMatrices[iAtom][(0*nc + i)*nc + j] = rescale_i * atom_mat[iAtom][(0*nc + i)*nc + j] * rescale[j]; // hmt
                    AtomMatrices[iAtom][(1*nc + i)*nc + j] = rescale_i * atom_mat[iAtom][(1*nc + i)*nc + j] * rescale[j]; // ovl
                } // j
            } // i
            // ToDo: treat Noco correctly
        } // iAtom

        view2D<double> rho(nblocks, n4x4x4, 0.0);

        int constexpr Noco = 1;

        std::vector<double> Veff(nblocks*n4x4x4, 0.);
        set(Veff.data(), nblocks*size_t(64), Vtot);
        stat += green_function::update_potential(plan, pg.grid_blocks(), Veff, AtomMatrices, echo, Noco);

        view2D<double> kpoint_mesh;
        // get a kpoint mesh controlled by +hamiltonian.kmesh.x .y .z, the same for each energy point
        auto const nkpoints = brillouin_zone::get_kpoint_mesh(kpoint_mesh);

        auto const energy_mesh = get_energy_mesh(echo);
        int const nEpoints = energy_mesh.size();

        for (int iEpoint{0}; iEpoint < nEpoints; ++iEpoint) {
            double const energy_weight = 1;

            std::complex<double> const energy = energy_mesh.at(iEpoint) + Fermi_level;
            char energy_parameter_label[96];
            std::snprintf(energy_parameter_label, 96, "(%g %s, %g %s)", (energy.real() - Fermi_level)*eV, _eV, energy.imag()*Kelvin, _Kelvin);
            if (echo > 5) std::printf("# energy parameter %s\n", energy_parameter_label);

            stat += green_function::update_energy_parameter(plan, energy, dV, echo, Noco);

            view2D<double> rho_E(nblocks, n4x4x4, 0.0);

            for (int ikpoint{0}; ikpoint < nkpoints; ++ikpoint) {
                double const *const kpoint = kpoint_mesh[ikpoint];
                double const kpoint_weight = kpoint[brillouin_zone::WEIGHT];

                if (echo + check > 5) std::printf("# solve Green function for E=%s, k-point=[%g %g %g]\n",
                                           energy_parameter_label, kpoint[0], kpoint[1], kpoint[2]);
                if (0 == check) {
                    stat += green_function::update_phases(plan, kpoint, echo, Noco);

                    view2D<double> rho_Ek(nblocks, n4x4x4, 0.0);

                    stat += solver_->solve(rho_Ek[0], nblocks, iterations, echo);

                    add_product(rho_E[0], nblocks*n4x4x4, rho_Ek[0], kpoint_weight); // accumulate density
                } // check

            } // ikpoint
            if (0 == check) {
                add_product(rho[0], nblocks*n4x4x4, rho_E[0], energy_weight); // accumulate density
            } // check
            if (echo + check > 2) std::printf("# solved Green function for E=%s\n", energy_parameter_label);
        } // iEpoint

        // interpolation density from 4*4*4 to 8*8*8 block could be done here
        if (echo > 3) std::printf("# interpolate density from 4x4x4 to 8x8x8\n");
        parallel_poisson::block_interpolation(rho_888, rho[0], pg, echo, 1., "density");

        return stat;
    } // integrate







#ifdef    NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    status_t test_integrator(int const echo=3) {
        if (echo > 1) std::printf("\n#\n# %s\n", __func__);
        double E_Fermi{0};
        std::vector<double> xyzZinso(0); // no atoms
        real_space::grid_t gc(4, 4, 4); // one block, isolated BCs by default, grid spacing 1.0
        parallel_poisson::load_balancing_t const lb(gc, MPI_COMM_WORLD, 4, echo);
        parallel_poisson::parallel_grid_t const pg(gc, lb, echo, "Interpolation");
        view2D<double> V_coarse(pg.n_local(), 4*4*4, 0.5);
        view2D<double> rhov_new(pg.n_local(), 8*8*8, 0.0);
        std::vector<uint32_t> num(0);
        data_list<double> atom_mat(num);
        std::vector<int32_t> numax_prj(0, 0);
        std::vector<double> sigma_prj(0, 1.);
        Integrator integrator(gc, xyzZinso, echo);
        if (echo > 1) std::printf("# %s: Integrator constructed\n\n", __func__);
        return integrator.integrate(rhov_new[0], E_Fermi, V_coarse[0], atom_mat, numax_prj, sigma_prj, pg, 1., 1., echo);
        if (echo > 1) std::printf("# %s: Integrator.integrate executed\n\n", __func__);
    } // test_integrator

    status_t test_energy_mesh(int const echo=5) {
        auto const energy_mesh = get_energy_mesh(echo);
        if (echo > 1) std::printf("# energy mesh with %ld points generated\n", energy_mesh.size());
        return 0;
    } // test_energy_mesh

    status_t all_tests(int const echo) {
        status_t stat(0);
        stat += test_energy_mesh(echo);
        stat += test_integrator(echo);
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace energy_contour
