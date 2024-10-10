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
        std::vector<std::vector<uint32_t>> const & list() const { return list_; }
        MPI_Comm comm()   const { return comm_; }
        uint32_t natoms() const { return natoms_; }

        status_t broadcast(
              data_list<double> & atom_data // result [natoms]
            , data_list<double> const & owner_data // input [na], only accessed in atom owner rank
            , char const *const what
            , int const echo=0 // log level
        ) const; // declaration only

        status_t allreduce(
              data_list<double> & owner_data // result [na], only correct in atom owner rank
            , data_list<double> const & atom_data // input [natoms]
            , char const *const what
            , double const factor=1 // scaling factor, usually g.dV(), the grid volume element
            , int const echo=0 // log level
        ) const; // declaration only

    private: // members
        std::vector<std::vector<uint32_t>> list_;
        std::vector<std::pair<uint32_t,uint32_t>> contributing_;
        MPI_Comm comm_ = MPI_COMM_NULL;
        uint32_t natoms_ = 0;
    }; // class AtomCommList_t

    // Discuss: make the following two routines methods of AtomCommList_t
    // (but it would make it more complicated if we wanted to generalize double --> typename T)

    status_t atom_data_broadcast(
          data_list<double> & atom_data // result [natoms]
        , data_list<double> const & owner_data // input [na], only accessed in atom owner rank
        , char const *const what
        , AtomCommList_t const & atom_comm_list
        , std::vector<uint32_t> const & global_atom_ids
        , int const echo=0 // log level
    ); // declaration only

    status_t atom_data_allreduce(
          data_list<double> & owner_data // result [na], only correct in atom owner rank
        , data_list<double> const & atom_data // input [natoms]
        , char const *const what
        , AtomCommList_t const & atom_comm_list
        , std::vector<uint32_t> const & global_atom_ids
        , double const factor=1 // scaling factor, usually g.dV(), the grid volume element
        , int const echo=0 // log level
    ); // declaration only

    status_t all_tests(int const echo=0); // declaration only

} // namespace atom_communication
