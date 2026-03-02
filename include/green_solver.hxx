#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert
#include <cstdint> // uint32_t
#include <vector> // std::vector<T>
#include <complex> // std::complex<real_t>

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "action_plan.hxx" // action_plan_t

// The purpose of this class is to hide the templated versions of tfQMRgpu solvers 
// with <float,2,1>, <float,2,2>, <double,2,1> and <double,2,2>

class green_solver_t {
public:

    green_solver_t() : action_{nullptr}, action_key_{0} {} // default constructor

    green_solver_t( // constructor
          action_plan_t* p
        , int const echo=0
        , int const check=0
    ); // declaration only

    green_solver_t(green_solver_t const &) = delete; // copy constructor
    green_solver_t(green_solver_t &&) = delete; // move constructor
    green_solver_t & operator=(green_solver_t const &) = delete; // copy assignment
    green_solver_t & operator=(green_solver_t && rhs); // custom move assignment operator
//  green_solver_t & operator=(green_solver_t &&) = delete; // move assignment

    ~green_solver_t(); // destructor, declaration only

    status_t solve(
          std::complex<double> rho[] // result: density rho data layout[plan.nCols][4*4*4]
        , uint32_t const ncubes // should match plan.nCols
        , int const iterations // maximum number of inner solver iterations
        , int const echo=0 // verbosity level
    ); // declaration only


private: // members

    void *action_ = nullptr; // pointer to templated action_t (will be casted according to key)
    int action_key_ = 0;

}; // class green_solver_t


namespace green_solver {

    status_t all_tests(int const echo=0); // declaration only

} // namespace green_solver
