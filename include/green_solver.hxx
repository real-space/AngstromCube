#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert
#include <cstdint> // uint32_t
#include <vector> // std::vector<T>
#include <complex> // std::complex<real_t>

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "action_plan.hxx" // action_plan_t

class green_solver_t {
public:

    green_solver_t( // constructor
          action_plan_t* p=nullptr
        , int const echo=0
    ); // declaration only

    ~green_solver_t(); // destructor, declaration only

    status_t solve(
          std::complex<double> rho[] // result: density rho data layout[plan.nCols][4*4*4]
        , uint32_t const nblocks // should match plan.nCols
        , int const iterations
        , int const imag=1 // index of the exported part 1:imaginary part, 0:real part
        , int const echo=0
    ); // declaration only


private: // members

    void *action_ = nullptr; // pointers to action (will be casted according to key)
    int action_key_ = 0;

}; // class green_solver_t


namespace green_solver {

    status_t all_tests(int const echo=0); // declaration only

} // namespace green_solver
