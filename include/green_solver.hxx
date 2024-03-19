#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert
#include <vector> // std::vector<T>

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "action_plan.hxx" // action_plan_t

class green_solver_t {
public:

    green_solver_t( // constructor
          action_plan_t* p=nullptr
        , int const echo=0
    ); // declaration only

    ~green_solver_t(); // declaration only

    status_t solve(
          double rho[] // result density rho [plan.nCols][4*4*4]
        , int const echo=0
    ); // declaration only


private: // members

    void *action_ = nullptr; // pointers to action
    int action_key_ = 0;

}; // class green_solver_t


namespace green_solver {

    status_t all_tests(int const echo=0); // declaration only

} // namespace green_solver
