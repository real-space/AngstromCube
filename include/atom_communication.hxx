// This file is part of AngstromCube under MIT License

#include <vector> // std::vector<T>
#include <utility> // std::pair
#include <cstdint> // uint32_t

#include "status.hxx" // status_t
#include "mpi_parallel.hxx" // MPI_COMM_WORLD
#include "data_list.hxx" // data_list<T>

namespace atom_communication {

    class AtomCommList_t {
    public:
        AtomCommList_t() {} // default constructor
        AtomCommList_t(
              size_t const n_all_atoms
            , std::vector<uint32_t> const & global_atom_ids // [natoms], global ids of contributing atoms on this rank
            , MPI_Comm const comm=MPI_COMM_WORLD
            , int const echo=0
        ); // constructor, declaration only

    public:
 
        status_t broadcast(
              data_list<double> & atom_contribution_data // result [natoms]
            , data_list<double> const & atom_owner_data // input [na], only accessed in atom owner rank
            , char const *const what // descriptive string for debugging
            , int const echo=0 // log level
        ) const; // declaration only

        status_t allreduce(
              data_list<double> & atom_owner_data // result [na], only correct in atom owner rank
            , data_list<double> const & atom_contribution_data // input [natoms]
            , char const *const what // descriptive string for debugging
            , double const factor=1 // scaling factor, usually g.dV(), the grid volume element
            , int const echo=0 // log level
        ) const; // declaration only

    private: // members
        std::vector<std::vector<uint32_t>> list_;
        std::vector<std::pair<uint32_t,uint32_t>> contributing_;
        MPI_Comm comm_ = MPI_COMM_NULL;
        uint32_t nprocs_; // == mpi_parallel::size(comm_)
        uint32_t me_;     // == mpi_parallel::rank(comm_)
    }; // class AtomCommList_t

    status_t all_tests(int const echo=0); // declaration only

} // namespace atom_communication
