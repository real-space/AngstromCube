#pragma once
// This file is part of AngstromCube under MIT License

#include "real_space.hxx" // ::grid_t

#include "status.hxx" // status_t
#include "mpi_parallel.hxx" // MPI_Comm, MPI_COMM_NULL, MPI_COMM_WORLD
#include "green_parallel.hxx" // ::RequestList_t, ::rank_int_t
#include "inline_math.hxx" // set

namespace parallel_poisson {

    class load_balancing_t {
    public:
        load_balancing_t() : nb_{0u, 0u, 0u} {} // default constructor
        load_balancing_t(
              real_space::grid_t const & g
            , MPI_Comm const comm=MPI_COMM_WORLD // MPI communicator
            , unsigned const n8=8 // number of grid points per cube edge
            , int const echo=0 // log-level
        ); // declaration only
    private:
        MPI_Comm comm_ = MPI_COMM_NULL; // which MPI communicator is to be used?
        uint32_t nb_[3]; // box of cubes
        uint32_t n_local_cubes_ = 0;
        double dom_center_[3];
        double load_ = 0;
        int32_t min_domain_[3];
        int32_t max_domain_[3];
        view3D<green_parallel::rank_int_t> owner_rank_;
    public:
        uint32_t const * grid_cubes() const { return nb_; }
        view3D<green_parallel::rank_int_t> const & owner_rank() const { return owner_rank_; }
        MPI_Comm comm() const { return comm_; }
        int32_t const * min_domain() const { return min_domain_; }
        int32_t const * max_domain() const { return max_domain_; }
        double load() const { return load_; }
        uint32_t n_local() const { return n_local_cubes_; }
    }; // class load_balancing_t



    class parallel_grid_t {
    public:
        parallel_grid_t() { set(nb_, 3, 0u); set(bc_, 3, int8_t(0)); nperiodic_ = 0; comm_ = MPI_COMM_NULL; set(h2_, 3, 1.); dVol_ = 1; } // default constructor
        parallel_grid_t(
              real_space::grid_t const & g
            , load_balancing_t const & lb
            , int const echo=0 // log-level
            , char const *const what="FD1"
        ); // declaration only
    private:
        green_parallel::RequestList_t requests_;
        std::vector<int64_t> remote_global_ids_; // may contain "-1"-entries, could be removed after setup keeping only a uint32_t n_remote_blocks;
        std::vector<int64_t> local_global_ids_;  // may not contain "-1"-entries
        view2D<uint32_t> star_; // local indices of 6 nearest finite-difference neighbors, star(n_local_cubes,6). Should be backed with GPU memory in the future
        std::vector<bool> inner_cell_; // mark those of the n_local cells that can start to execute a stencil without waiting for remote data
        MPI_Comm comm_; // which MPI communicator is to be used?
        double h2_[3];
        double dVol_;
        uint32_t nb_[3]; // box of cubes
        int8_t bc_[3];
        uint8_t nperiodic_;
    public:
        double const * get_prefactors() const { return h2_; }
        uint32_t const * grid_cubes() const { return nb_; }
        uint32_t n_local()  const { return local_global_ids_.size();  } // number of cubes owned by this MPI rank
        uint32_t n_remote() const { return remote_global_ids_.size(); } // number of cubes requested by this MPI rank
        uint32_t const* star() const { return (uint32_t const*)star_.data(); }
        uint32_t star_dim() const { return star_.stride(); }
        std::vector<int64_t> const & local_ids() const { return local_global_ids_; }
        std::vector<int64_t> const & remote_ids() const { return remote_global_ids_; }
        green_parallel::RequestList_t const & requests() const { return requests_; }
        std::vector<bool> const & inner_cell() const { return inner_cell_; }
        bool all_periodic_boundary_conditions() const { return 3 == nperiodic_; }
        MPI_Comm comm() const { return comm_; }
        double dV() const { return dVol_; }
    }; // class parallel_grid_t

    template <typename real_t>
    status_t solve(
          real_t x[] // result to Laplace(x)/(-4*pi) == b, only rank-local cubes, data layout x[][8*8*8]
        , real_t const b[] // right hand side b          , only rank-local cubes, data layout x[][8*8*8]
        , parallel_grid_t const & g8 // parallel grid descriptor
        , char const method='c' // solver method c:conjugate-gradient, s:steepest-descent
        , int const echo=0 // log level
        , float const threshold=3e-8f // convergence criterion
        , float *residual=nullptr // residual that was reached
        , int const maxiter=199 // maximum number of iterations 
        , int const miniter=3  // minimum number of iterations
        , int restart=4096 // number of iterations before restart, 1:steepest descent
        , double *inner_xx_bb=nullptr // export the last inner product
    ); // declaration only

    template <typename real_t=double>
    status_t block_interpolation(
          real_t       *const v888 // result array, data layout v888[n_local_cubes][8*8*8]
        , real_t const *const v444 // input  array, data layout v444[n_local_cubes][4*4*4]
        , parallel_grid_t const & pg // descriptor, must be prepared with "3x3x3"
        , int const echo=0 // log level
        , double const factor=1
        , char const *const what="!"
    ); // declaration only


    status_t all_tests(int const echo=0); // declaration only

} // namespace parallel_poisson
