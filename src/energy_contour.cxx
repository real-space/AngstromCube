// This file is part of AngstromCube under MIT License

#include <cstdint> // int8_t
#include <cassert> // assert
#include <cstdio> // std::printf, ::snprintf
#include <vector> // std::vector<T>

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

namespace energy_contour {

    Integrator::Integrator( // implementation of constructor
          real_space::grid_t const & gc // coarse grid descriptor
        , std::vector<double> const & xyzZinso // all atoms
        , int const echo // verbosity
    ) {
        plan_ = new action_plan_t(); // CPU memory for the plan
        auto const stat = green_function::construct_Green_function(*plan_,
                            gc.grid_points(), gc.boundary_conditions(), gc.grid_spacings(),
                            xyzZinso, echo);
        if (stat) warn("construct_Green_function returned status= %i", int(stat));
    } // constructor

    template <typename real_t>
    status_t solve(
          double rho[] // result density [plan.nCols][4*4*4]
        , action_plan_t const & plan
        , int const echo
    ) {
        if (echo > 0) std::printf("# solve ...\n");
        return 0;
    } // solve

    status_t Integrator::integrate(
          double rho_888[] // resulting density in [nblocks][8*8*8] data layout
        , double & Fermi_level // Fermi level
        , double const Vtot[] // input potential in [nblocks][4*4*4]
        , data_list<double> const & atom_mat // atomic_Hamiltonian elements, only in atom owner ranks
        , parallel_poisson::load_balancing_t const & lb
        , parallel_poisson::parallel_grid_t const & pg
        , double const n_electrons // =1 // required total number of electrons 
        , double const dV // =1 // grid volume element
        , int const echo // =0 // log level
    ) {
        status_t stat(0);

        int const check = control::get("check", 0.);
        if (echo > 0) std::printf("# energy_contour::integration(E_Fermi=%g %s, %g electrons, echo=%d) +check=%i\n", Fermi_level*eV, _eV, n_electrons, echo, check);

        auto const nblocks = pg.n_local();
        assert(lb.n_local() == nblocks);
        assert(nullptr != plan_);
        auto & plan = *plan_;
        if (plan.nCols != nblocks) error("model assumes that each local block has one RHS, found n_local= %d and p.nRHS= %d", nblocks, plan.nCols);

        auto const nAtoms = atom_mat.nrows();
        std::vector<std::vector<double>> AtomMatrices(nAtoms);
        for (int iAtom{0}; iAtom < nAtoms; ++iAtom) {
            auto const nc2 = atom_mat.ncols(iAtom);
            AtomMatrices[iAtom] = std::vector<double>(nc2);
            set(AtomMatrices[iAtom].data(), nc2, atom_mat[iAtom]); // copy
        } // iAtom

        view2D<double> rho(nblocks, 4*4*4, 0.0);

        int constexpr Noco = 1;
        double constexpr scale_H = 1;

        std::vector<double> Veff(nblocks*4*4*4, 0.); // ToDo: fille Veff with Vtot
        stat += green_function::update_potential(plan, pg.grid_blocks(), Veff, lb.owner_rank(), echo, Noco);

        for (int ienergy{0}; ienergy < 1; ++ienergy) {
            double const energy_weight = 1;
            
            std::complex<double> const energy((ienergy + 1)*0.5, 0.125);
            if (echo > 5) std::printf("# energy parameter  (%g %s, %g %s)\n", (energy.real() - Fermi_level)*eV, _eV, energy.imag()*Kelvin, _Kelvin);
            stat += green_function::update_energy_parameter(plan, energy, AtomMatrices, dV, scale_H, echo, Noco, &atom_req_);

            view2D<double> rho_E(nblocks, 4*4*4, 0.0);

            for (int ikpoint{0}; ikpoint < 1; ++ikpoint) {
                double const kpoint[] = {0, 0, 0};
                double const kpoint_weight = 1;

                stat += green_function::update_phases(plan, kpoint, echo, Noco);

                view2D<double> rho_Ek(nblocks, 4*4*4, 1.0);
                stat += solve<float>(rho_Ek[0], *plan_, echo);

                add_product(rho_E[0], nblocks*size_t(4*4*4), rho_Ek[0], kpoint_weight); // accumulate density
            } // ikpoint
            add_product(rho[0], nblocks*size_t(4*4*4), rho_E[0], energy_weight); // accumulate density
        } // ienergy

        // interpolation density from 4*4*4 to 8*8*8 block could be done here
        parallel_poisson::block_interpolation(rho_888, rho[0], pg, echo, 1., "density");

        return stat;
    } // integrate









#ifdef    NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    status_t test_integration(int const echo=3) {
        double E_Fermi{0};
        parallel_poisson::load_balancing_t const lb;
        parallel_poisson::parallel_grid_t const pg;
        view2D<double> V_coarse(lb.n_local(), 4*4*4, 0.5);
        view2D<double> rhov_new(lb.n_local(), 8*8*8, 0.0);
        std::vector<uint32_t> num(0);
        data_list<double> atom_mat(num);
        Integrator integrator;
        return integrator.integrate(rhov_new[0], E_Fermi, V_coarse[0], atom_mat, lb, pg, 1., 1., echo);
    } // test_integration

    status_t all_tests(int const echo) {
        status_t stat(0);
        stat += test_integration(echo);
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace energy_contour
