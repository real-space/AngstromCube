// This file is part of AngstromCube under MIT License

#include <cstdint> // int8_t
#include <cassert> // assert
#include <cstdio> // std::printf, ::snprintf
#include <vector> // std::vector<T>
#include <complex> // std::complex<real_t>, ::real, ::imag
#include <algorithm> // std::min, ::max

#include "energy_contour.hxx"

#include "control.hxx" // ::get
#include "display_units.h" // eV, _eV, Ang, _Ang, Kelvin, _Kelvin
#include "mpi_parallel.hxx" // MPI_Comm, ::rank, ::comm, MPI_COMM_WORLD
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
#include "verify_benchmark.hxx" // ::verify
#include "energy_mesh.hxx" // ::Complex, ::get_energy_mesh
#include "simple_stats.hxx" // ::Stats<>
#include "simple_timer.hxx" // SimpleTimer
#include "green_parallel.hxx" // ::RequestList_t

namespace energy_contour {

    typedef energy_mesh::Complex Complex;

    Integrator::Integrator( // implementation of constructor
          real_space::grid_t const & gc // coarse grid descriptor
        , std::vector<double> const & xyzZinso // [n_all_atoms*8]
        , parallel_poisson::load_balancing_t const & lb // load balancing of the dense grid
        , int const echo // verbosity
        , int const check
    ) {
        if (echo > 0) std::printf("# construct %s with grid=[%d %d %d]\n", __func__, gc[0], gc[1], gc[2]);
        plan_ = new action_plans_t(); // CPU memory for the plans
        auto stat = green_function::construct_Green_function(*plan_, // result
                            gc.grid_points(), gc.boundary_conditions(), gc.grid_spacings(), // grid info
                            xyzZinso, // atom info
                            nullptr, // neutral weigth info, equivalent to all weights==1
                            lb.comm(),
                            lb.global_ids(),
                            lb.owner_rank().data(),
                            echo);
        if (stat) warn("construct_Green_function returned status= %i", int(stat));

        std::vector<load_balancer::WeightInfo> const no_weights(0);
        std::vector<float> weights = load_balancer::calculate_weights(plan_->global_source_indices, 
                                    no_weights, // plan_->dyadic_plan.weight_infos, 
                                    plan_->owner_rank_.size());
        mpi_parallel::max(weights.data(), lb.comm(), weights.size()); // MPI_Allreduce(MPI_Max)

        for (size_t wi = 0; wi < weights.size(); ++wi) {
            auto const w = weights[wi];
            if (w <= 0) {
                std::printf("# rank#%i has 0 or negative weight at global index: %li check: %f\n", mpi_parallel::rank(lb.comm()), wi, w);
            }
        } // wi

        // if the two distributions are not the same, we need to redistribute the 4x4x4 density cubes
        if (!weights.empty()) { // ToDo: check if we need redistribution at all

            delete plan_;
            plan_ = new action_plans_t();
            auto stat = green_function::construct_Green_function(*plan_, // result
                            gc.grid_points(), gc.boundary_conditions(), gc.grid_spacings(), // grid info
                            xyzZinso, // atom info
                            weights.data(), // weigth info
                            lb.comm(),
                            lb.global_ids(),
                            lb.owner_rank().data(),
                            echo);
            if (stat) warn("2nd construct_Green_function returned status= %i", int(stat));

            uint32_t const nb[] = {uint32_t(gc[0] >> 2), uint32_t(gc[1] >> 2), uint32_t(gc[2] >> 2)}; // divide grid by 4
            if (echo > 5) { std::printf("# rank#%i plan to redistribute %ld 4x4x4 density cubes to %ld cubes\n",
                mpi_parallel::rank(lb.comm()), plan_->global_source_indices.size(), lb.global_ids().size()); }
            req_ = new green_parallel::RequestList_t(lb.global_ids() // requests
                                      , plan_->global_source_indices // offerings
                                      , plan_->owner_rank_.data() // where to find it, [nb[Z]*nb[Y]*nb[X]]
                                      , nb // number of cubes
                                      , lb.comm() // communicator
                                      , echo // log-level
                                      , "redistribute" // what
                                    );
        } // weights empty

        if (echo > 2) { std::printf("# generate a parallel grid descriptor for the interpolation of densities\n"); }
        pg_ = new parallel_poisson::parallel_grid_t(gc, lb, echo, "Interpolation");

        int const nsub = plan_->plans.size();

        solver_.resize(nsub);
        #pragma omp parallel for
        for (int isub = 0; isub < nsub; ++isub) {
            solver_.at(isub) = green_solver_t(& plan_->plans[isub], echo, check);
        } // isub
        if (echo > 7) std::printf("# constructed %s\n", __func__);
    } // constructor

    Integrator::~Integrator() {
#ifdef    DEBUGGPU
        std::printf("\n# destruct %s\n", __func__);
#endif // DEBUGGPU
        solver_.resize(0);
        if (plan_) { delete plan_; }
        if (pg_) { delete pg_; }
    } // destructor

    template <typename real_t>
    real_t sum(real_t const rho[], size_t const n) {
        real_t s{0};
        for (size_t i{0}; i < n; ++i) {
            s += rho[i];
        } // i
        return s;
    } // sum

    template <typename real_t>
    double maxval(real_t const rho[], size_t const n) {
        double s{0};
        for (size_t i{0}; i < n; ++i) {
            s = std::max(s, std::real(rho[i]));
        } // i
        return s;
    } // maxval

    status_t Integrator::integrate(
          double *const rho_888 // resulting density in [ncubes*8*8*8] data layout
        , double & Fermi_level // Fermi level
        , std::vector<double> const & Vtot // input potential in [ncubes*4*4*4]
        , data_list<double> const & atom_mat // atomic_Hamiltonian elements, only in atom owner ranks
        , std::vector<int32_t> const & numax_prj
        , std::vector<double> const & sigma_prj
        , double const n_electrons // =1 // required total number of electrons 
        , double const dV // =1 // grid volume element on dense grid
        , int const echo // =0 // log level
        , int const check // =0
        , int const scf_iteration_number // =0
    ) {
        status_t stat(0);
        size_t constexpr n4x4x4 = 4*4*4;
        size_t constexpr n8x8x8 = 8*8*8;
        auto const dVc = 8*dV; // grid volume element on the dense grid

        assert(nullptr != pg_);
        auto const comm = pg_->comm();
        // if (echo > 8) { std::printf("\n# energy_contour::integration comm= %ld, MPI_COMM_WORLD= %ld, MPI_COMM_NULL= %ld\n", 
        //                             int64_t(comm), int64_t(MPI_COMM_WORLD), int64_t(MPI_COMM_NULL)); std::fflush(stdout); }
        auto const me = mpi_parallel::rank(comm);
        int const sync = (0 != control::get("energy_contour.integrate.mpi.sync", 1.)); // set .sync=0 to measure load imbalance

        int const max_iterations = control::get("green_solver.iterations", 99.);
        if (echo > 0) std::printf("\n# energy_contour::integration(E_Fermi=%g %s, %g electrons, echo=%d) +check=%i\n", Fermi_level*eV, _eV, n_electrons, echo, check);

        auto const ncubes = pg_->n_local(); // number of 8x8x8 cubes treated on the dense grid
        assert(nullptr != plan_);
        auto & plan = *plan_;

        uint32_t const nrhs = plan.global_source_indices.size(); // number of 4x4x4 cubes treated by the Green function solver
        if (echo > 1) { std::printf("# dense grid domain has %d cubes, Green solver works on %d RHS cubes\n", ncubes, nrhs); }

        int constexpr Noco = 1;

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
            if (Noco > 1) error("not prepared for Noco= %d", Noco);
        } // iAtom

        // MPI communicate the potential elements and place them in the GPU memory of plan
        stat += green_function::update_potential(plan, Vtot, AtomMatrices, echo, Noco);

#ifdef    DEVEL
        int const verify_pot = control::get("verify.potential", 0.);
        if (verify_pot) {
            if (echo > 3) std::printf("\n# +verify.potential=%i\n", verify_pot);
            assert(plan.global_source_indices.size() == ncubes && false); // this is not correct any more, ToDo: repair
            // view2D<double> pot_444(Vtot.data(), n4x4x4); // wrap
            // auto const stat_verify = verify_benchmark::verify(pot_444, plan.global_source_indices.data(), ncubes, echo);
            // if (0 != stat_verify) warn("ran with +verify.potential=%i --> status= %i", verify_pot, int(stat_verify));
            // stat += std::abs(stat_verify);
            // if (echo > 0) std::fflush(stdout);
        } // verify_pot
#endif // DEVEL

        view2D<double> kpoint_mesh;
        // get a kpoint mesh controlled by +hamiltonian.kmesh.x .y .z, the same for each energy point
        auto const nkpoints = brillouin_zone::get_kpoint_mesh(kpoint_mesh);

        std::vector<Complex> energy_weights;
        auto const energies = energy_mesh::get(energy_weights, echo);
        int const nEpoints = energies.size();
        if (echo > 0) std::printf("# energy_contour::integration with %d energy points and %d k-points\n", nEpoints, nkpoints);

        int const echo_dos = 10*(0 == control::get("energy_contour.matsubara", 0.)); // more verbose in a DoS (density-of-state) calculation

        if (echo*echo_dos > 5 && nEpoints > 0) {
            auto const emin = energies.at(0).real(), emax = energies.at(nEpoints - 1).real();
            std::printf("# show density of states from %g to %g %s, %d equidistant points spaced %g %s, imaginary part %g %s\n",
                (emin + Fermi_level)*eV, (emax + Fermi_level)*eV, _eV, nEpoints, std::abs(emax - emin)/std::max(1, nEpoints - 1)*eV, _eV, energies.at(0).imag()*eV, _eV);
        } // show DoS

        Complex constexpr zero = 0;
        view2D<Complex> rho_c(nrhs, n4x4x4, zero); // complex density
        view2D<Complex> res_c(nrhs, n4x4x4, zero); // complex response density
        Complex res_point{zero};

        char timer_label[64]; std::snprintf(timer_label, 64, "Green solver in SCF-iteration#%i", scf_iteration_number);
        SimpleTimer integrate_timer(strip_path(__FILE__), __LINE__, timer_label, 0); // start timer

        int const nsub = plan.plans.size();

        simple_stats::Stats<> iterations_needed_Ek;
        for (int iEpoint{0}; iEpoint < nEpoints; ++iEpoint) {
            auto const energy_weight = energy_weights[iEpoint];

            Complex const energy = energies.at(iEpoint) + Fermi_level;
            char energy_parameter_label[64];
            std::snprintf(energy_parameter_label, 64, "(%g %s, %g %s)", (energy.real() - Fermi_level)*eV, _eV, energy.imag()*Kelvin, _Kelvin);
            if (echo > 7) std::printf("# energy parameter#%i %s with weight (%g, %g)\n", iEpoint, energy_parameter_label, std::real(energy_weight), std::imag(energy_weight));

            view2D<Complex> rho_E(nrhs, n4x4x4, zero);
            simple_stats::Stats<> iterations_needed_k;

            #pragma omp parallel for
            for (int isub = 0; isub < nsub; ++isub) {
                auto & p = plan.plans[isub];

                stat += green_function::update_energy_parameter(p, plan, energy, dVc, echo, Noco);

                for (int ikpoint{0}; ikpoint < nkpoints; ++ikpoint) {
                    double const *const kpoint = kpoint_mesh[ikpoint];
                    Complex const kpoint_weight = kpoint[brillouin_zone::WEIGHT];

                    if (echo + check > 8) std::printf("# solve Green function for E=%s, k-point=[%g %g %g] weight= %g\n",
                                                    energy_parameter_label, kpoint[0], kpoint[1], kpoint[2], kpoint[3]);
                    if (0 == check) {
                        stat += green_function::update_phases(p, kpoint, echo >> 3, Noco);

                        view2D<Complex> rho_Ek(nrhs, n4x4x4, zero);

                        // ******************************
                        // *** Core solver invokation ***
                        // ******************************
                        
                        solver_[isub].solve(rho_Ek[0], nrhs, max_iterations, echo);

                        // ******************************


                        add_product(rho_E[0], nrhs*n4x4x4, rho_Ek[0], kpoint_weight); // accumulate complex density over k-points
                        // if (sync) {
                        //     auto const rho_integral = mpi_parallel::sum(sum(rho_Ek[0], nrhs*n4x4x4).imag(), comm)*dVc; // MPI synchronization point
                        //     if (echo > 11) std::printf("# Green function solution for E=%s, k-point=[%g %g %g] has %g electrons, %d iterations\n",
                        //                                     energy_parameter_label, kpoint[0], kpoint[1], kpoint[2], rho_integral, p.iterations_needed);
                        // } // sync
                        iterations_needed_k.add(p.iterations_needed); // omp critical
                    } // check

                } // ikpoint

            } // isub

            if (0 == check) {
                if (sync) {
                    auto const rho_integral = mpi_parallel::sum(sum(rho_E[0], nrhs*n4x4x4).imag(), comm)*dVc; // MPI synchronization point
                    auto const rho_realpart = mpi_parallel::sum(sum(rho_E[0], nrhs*n4x4x4).real(), comm)*dVc; // MPI synchronization point
                    if (echo + echo_dos > 5) { std::printf("# Green function solution for E=%s has %g electrons, real part %g\n",
                                                  energy_parameter_label, rho_integral, rho_realpart); std::fflush(stdout); }
                } // sync
                // accumulate density over E-points
                add_product(rho_c[0], nrhs*n4x4x4, rho_E[0], energy_weight);
                if (echo > 7) { std::printf("# energy parameter#%i iterations needed %s\n", iEpoint, iterations_needed_k.interval().c_str()); std::fflush(stdout); }

            } else if (echo > 7) std::printf("# solve Green function for E=%s\n", energy_parameter_label);

            iterations_needed_Ek.add(iterations_needed_k);

            if (iEpoint < nEpoints - 2) {
                // ToDo: accumulate a response density to derive the new density w.r.t. the Fermi level 
                //       in order to correct the density to the right number of electrons
                Complex const wgt = ((nEpoints - 1 == iEpoint) ? 1. : -1.);
                add_product(res_c[0], nrhs*n4x4x4, rho_E[0], wgt);
                res_point += energy*wgt;
            } // last two
        } // iEpoint

        auto const green_solver_took = integrate_timer.stop(); // stop timer
        {
            simple_stats::Stats<> green_time_stats;
            green_time_stats.add(green_solver_took);
            mpi_parallel::allreduce(green_time_stats, comm); // MPI synchronization point
            if (0 == check && echo > 2) {
                std::printf("# Green solver in SCF-iteration#%i took %s seconds (sync=%i)\n", 
                    scf_iteration_number, green_time_stats.interval().c_str(), sync);
            } // echo
        }


        if (0 == check && echo > 3) {
            std::printf("# iterations needed %s\n", iterations_needed_Ek.interval().c_str()); std::fflush(stdout);
        } // check echo

        if (nEpoints < 2) warn("unable to eval a meaningful response density with less than 2 energy points, found %d", nEpoints);


        view2D<double> rho_444(nrhs, n4x4x4, 0.0); // real-valued density
        view2D<double> rho_res(nrhs, n4x4x4, 0.0); // real-valued response density
        res_point = (zero != res_point) ? 1./res_point : 1;
        for (uint32_t ib{0}; ib < nrhs; ++ib) {
            for (int i444{0}; i444 < n4x4x4; ++i444) {
                rho_444(ib,i444) =  rho_c(ib,i444).imag();
                rho_res(ib,i444) = (res_c(ib,i444)*res_point).imag();
            } // i444
        } // ib cube index

        {
            auto const rho_integral = mpi_parallel::sum(sum(rho_444[0], nrhs*n4x4x4), comm)*dVc; // MPI synchronization point
            if (echo + check > 3) std::printf("# solved density has %g electrons\n", rho_integral);
            if (echo > 4) std::printf("# rank#%i maxval rho= %g a.u.\n", me, maxval(rho_444[0], nrhs*n4x4x4));
        }

        {
            auto const rho_integral = mpi_parallel::sum(sum(rho_res[0], nrhs*n4x4x4), comm)*dVc; // MPI synchronization point
            if (echo + check > 3) std::printf("# solved response density has %g electrons\n", rho_integral);
            // the response density should be positive semidefinite (i.e. integral >= 0) since higher Fermi --> more electrons
        }

#ifdef    DEVEL
        int const verify = control::get("verify.benchmark", 0.);
        if (verify) {
            assert(plan_->global_source_indices.size() == nrhs);
            auto const stat_verify = verify_benchmark::verify(rho_444, plan_->global_source_indices.data(), nrhs, echo);
            if (0 != stat_verify) warn("ran with +verify.benchmark=%d --> status= %i", verify, int(stat_verify));
            stat += std::abs(stat_verify);
        } // verify
#endif // DEVEL

        // ToDo: add response density until we match the Fermi level

        
        // check if we need to MPI-redistribute the density cubes
        double const *rho_on444{nullptr};
        std::vector<double> rho_red;
        if (req_) {
            // the domain decomposition of the 8x8x8 grid does not match the 4x4x4 grid
            if (echo > 5) { std::printf("# rank#%i redistribute %d 4x4x4 density cubes to %d cubes\n", me, nrhs, ncubes); }
            rho_red.resize(ncubes*n4x4x4, 0.0);
            req_->exchange(rho_red.data(), rho_444[0], n4x4x4, echo, "redistribute density");
            rho_on444 = rho_red.data();
        } else {
            rho_on444 = rho_444[0]; // matching domain decompositions
        }

        // Idea: combine the redistribution above with the data exchange inside the interpolation to save some overheads
        
        if (echo > 3) std::printf("# interpolate density from 4x4x4 to 8x8x8\n");
        // interpolate density from 4*4*4 to 8*8*8 cubes including MPI communication for 3x3x3 halos
        assert(nullptr != pg_);
        parallel_poisson::cube4x4x4_interpolation(rho_888, rho_on444, *pg_, echo, 1., "density");

        {
            auto const rho_integral = mpi_parallel::sum(sum(rho_888, ncubes*n8x8x8), comm)*dV; // MPI synchronization point
            if (echo + check > 3) std::printf("# interpolated density has %g electrons\n", rho_integral);
            if (echo > 4) std::printf("# rank#%i maxval rho= %g a.u.\n", me, maxval(rho_888, ncubes*n8x8x8));
        }

        if (echo > 3) std::printf("# density integrated over %d energy points\n", nEpoints);
        return stat;
    } // integrate





#ifdef    NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    status_t test_integrator(int const echo=3) {
        if (echo > 1) std::printf("\n#\n# %s\n", __func__);
        status_t stat(0);
        auto const comm = mpi_parallel::comm(); // for tests
        { // scope
            // spectrum of an isolated box  0.648721 1.29538_x3 1.94203_x3 2.34449_x3 2.58869 2.99114_x6 3.55103_x3 3.6378_x3 ...
            //  ... 4.04025_x3 4.19768_x6 4.68691_x3 4.84434_x3 5.24679_x6 5.73601 5.89345_x6 6.45333_x3 6.94256_x3 7.09999_x3 8.1491_x3 9.35564

            double E_Fermi{1.0};
            std::vector<double> xyzZinso(0); // no atoms
            real_space::grid_t gc(4, 4, 4); // one block, isolated BCs by default, grid spacing 1.0
            parallel_poisson::load_balancing_t const lb(gc, comm, 4, echo);
            parallel_poisson::parallel_grid_t const pg(gc, lb, echo, "Interpolation");
            std::vector<double> V_coarse(pg.n_local()*4*4*4, 0.5);
            view2D<double> rhov_new(pg.n_local(), 8*8*8, 0.0);
            std::vector<uint32_t> num(0);
            data_list<double> atom_mat(num);
            std::vector<int32_t> numax_prj(0, 0);
            std::vector<double> sigma_prj(0, 1.);
            Integrator integrator(gc, xyzZinso, lb, echo);
            if (echo > 1) std::printf("# %s: Integrator constructed\n\n", __func__);
            stat += integrator.integrate(rhov_new[0], E_Fermi, V_coarse, atom_mat, numax_prj, sigma_prj, 1., 1., echo);
        } // scope (so all destructors belonging to this test are called before the next log message)
        if (echo > 1) std::printf("# %s: Integrator.integrate executed\n\n", __func__);
        return stat;
    } // test_integrator

    status_t all_tests(int const echo) {
        status_t stat(0);
        stat += test_integrator(echo);
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace energy_contour
