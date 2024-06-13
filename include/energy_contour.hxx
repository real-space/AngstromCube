#pragma once
// This file is part of AngstromCube under MIT License

#include <algorithm> // std::swap

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
            , int const check=0
        ); // constructor, declaration only

        ~Integrator(); // destructor, declaration only

        Integrator(Integrator const &) = delete; // copy constructor
        Integrator(Integrator &&) = delete; // move constructor
        Integrator & operator=(Integrator const &) = delete; // copy assignment

        Integrator & operator=(Integrator && rhs) { // move assignment
            std::swap(this->solver_ , rhs.solver_);
            std::swap(this->plan_   , rhs.plan_  );
            return *this;
        } // move assignment

    public: // members TODO: go private
        action_plan_t *plan_ = nullptr;
        green_solver_t *solver_ = nullptr;

    public: // methods

        status_t integrate(
              double rho_new[] // result density in [nblocks][8*8*8] data layout
            , double & Fermi_level // Fermi level
            , double const Vtot[] // input potential in [nblocks][4*4*4], coarsening could be performed here...
            , data_list<double> const & atom_mat // atomic_Hamiltonian elements, only in atom owner ranks
            , std::vector<int32_t> const & numax_prj
            , std::vector<double> const & sigma_prj
            , parallel_poisson::parallel_grid_t const & pg
            , double const n_electrons=1 // required total number of electrons 
            , double const dV=1 // grid volume element
            , int const echo=0 // log level
            , int const check=0
        ); // declaration only

    }; // class Integrator

    int Gauss_Legendre_quadrature(double x[], double w[]
        , unsigned const number, int const echo=0); // declaration only

    int Gauss_Fermi_Dirac_quadrature(double x[], double w[]
        , unsigned const number, int const echo=0); // declaration only

    status_t all_tests(int const echo=0); // declaration only

} // namespace energy_contour
