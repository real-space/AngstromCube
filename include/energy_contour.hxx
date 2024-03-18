#pragma once
// This file is part of AngstromCube under MIT License

#include "status.hxx" // status_t
#include "mpi_parallel.hxx" // MPI_Comm, MPI_COMM_WORLD
#include "parallel_poisson.hxx" // ::parallel_grid_t
#include "action_plan.hxx" // action_plan_t
#include "real_space.hxx" // ::grid_t
#include "data_list.hxx" // data_list<T>
#include "green_solver.hxx" // green_solver_t

namespace energy_contour {

    class Integrator {

    public: // constructors

        Integrator() {} // default constructor
        Integrator(
              real_space::grid_t const & gc // coarse grid descriptor
            , std::vector<double> const & xyzZinso // all atoms
            , int const echo=0 // verbosity
        ); // constructor

    public: // members TODO: go private
        action_plan_t *plan_ = nullptr;
        green_solver_t solver_;
        green_parallel::RequestList_t atom_req_;

    public: // methods

        status_t integrate(
              double rho_new[] // result density in [nblocks][8*8*8] data layout
            , double & Fermi_level // Fermi level
            , double const Vtot[] // input potential in [nblocks][4*4*4], coarsening could be performed here...
            , data_list<double> const & atom_mat // atomic_Hamiltonian elements, only in atom owner ranks
            , parallel_poisson::load_balancing_t const & lb
            , parallel_poisson::parallel_grid_t const & pg
            , double const n_electrons=1 // required total number of electrons 
            , double const dV=1 // grid volume element
            , int const echo=0 // log level
        ); // declaration only

    }; // class Integrator


    status_t all_tests(int const echo=0); // declaration only

} // namespace energy_contour
