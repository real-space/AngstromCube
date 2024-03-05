#pragma once
// This file is part of AngstromCube under MIT License

#include "real_space.hxx" // ::grid_t

#include "status.hxx" // status_t
#include "mpi_parallel.hxx" // MPI_Comm, MPI_COMM_NULL
#include "green_parallel.hxx" // ::RequestList_t
#include "inline_math.hxx" // set

namespace parallel_poisson {

    class parallel_grid_t {
    public:
        parallel_grid_t() { set(nb_, 3, 0u); set(bc_, 3, int8_t(0)); nperiodic_ = 0; comm_ = MPI_COMM_NULL; set(h2_, 3, 1.); dVol_ = 1; } // default constructor
        parallel_grid_t(
              real_space::grid_t const & g
            , MPI_Comm const comm // MPI communicator
            , unsigned const n8=8 // number of grid points per block edge
            , int const echo=0 // log-level
        ); // declaration only
    private:
        green_parallel::RequestList_t requests_;
        std::vector<int64_t> remote_global_ids_; // may contain "-1"-entries, could be removed after setup keeping only a uint32_t n_remote_blocks;
        std::vector<int64_t> local_global_ids_;  // may not contain "-1"-entries
        view2D<int32_t> star_; // local indices of 6 nearest finite-difference neighbors, star(n_local_blocks,6). Should be backed with GPU memory in the future
        double dVol_;
        MPI_Comm comm_; // which MPI communicator is to be used?
        double h2_[3];
        uint32_t nb_[3]; // box of blocks
        int8_t bc_[3];
        uint8_t nperiodic_;
        // std::vector<bool> inner_cell_; // mark those of the n_local cells, i.e. cells than can start to execute a stencil without waiting for remote data
    public:
        double const * get_prefactors() const { return h2_; }
        uint32_t const * grid_blocks() const { return nb_; }
        uint32_t n_local()  const { return local_global_ids_.size();  } // number of blocks owned by this MPI rank
        uint32_t n_remote() const { return remote_global_ids_.size(); } // number of blocks requested by this MPI rank
        int32_t const* getStar() const { return (int32_t const*)star_.data(); }
        std::vector<int64_t> const & local_ids() const { return local_global_ids_; }
        std::vector<int64_t> const & remote_ids() const { return remote_global_ids_; }
        green_parallel::RequestList_t const & get_requests() const { return requests_; }
        bool all_periodic_boundary_conditions() const { return 3 == nperiodic_; }
        MPI_Comm get_comm() const { return comm_; }
        double dV() const { return dVol_; }
    }; // class parallel_grid_t

    template <typename real_t>
    status_t solve(
          real_t x[] // result to Laplace(x)/(-4*pi) == b, only rank-local blocks, data layout x[][8*8*8]
        , real_t const b[] // right hand side b          , only rank-local blocks, data layout x[][8*8*8]
        , parallel_grid_t const & g8 // parallel grid descriptor
        , char const method='c' // solver method c:conjugate-gradient, s:steepest-descent
        , int const echo=0 // log level
        , float const threshold=3e-8 // convergence criterion
        , float *residual=nullptr // residual that was reached
        , int const maxiter=199 // maximum number of iterations 
        , int const miniter=3  // minimum number of iterations
        , int restart=4096 // number of iterations before restart, 1:steepest descent
    ); // declaration only

    status_t all_tests(int const echo=0); // declaration only

} // namespace parallel_poisson
