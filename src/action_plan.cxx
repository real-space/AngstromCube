// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert

#include "action_plan.hxx"

#include "green_memory.hxx" // free_memory
#include "green_function.hxx" // ::construct_Green_function

#ifdef    DEBUG
  #define green_debug_printf(...) { std::printf(__VA_ARGS__); std::fflush(stdout); }
#else  // DEBUG
  #define green_debug_printf(...)
#endif // DEBUG

    action_plan_t::action_plan_t(int const echo) { // constructor 
        if (echo > 0) std::printf("# default constructor for %s\n", __func__);
    } // constructor

    action_plan_t::action_plan_t(
        uint32_t const ng[3] // numbers of grid points of the unit cell in with the potential is defined
      , int8_t const bc[3] // boundary conditions in {Isolated, Periodic, Vacuum, Repeat}
      , double const hg[3] // grid spacings
      , std::vector<double> const & xyzZinso // [natoms*8]
      , int const echo // =0 // log-level
      , int const Noco // =2
    ) {
        if (echo > 0) std::printf("# constructor for %s --> green_function::construct_Green_function\n", __func__);
        green_function::construct_Green_function(*this, ng, bc, hg, xyzZinso, echo, Noco);
    } // constructor

    action_plan_t::~action_plan_t() { // destructor
        green_debug_printf("# destruct %s\n", __func__);
        free_memory(RowStart);
        free_memory(rowindx);
        // free_memory(source_coords);
        // free_memory(target_coords);
        free_memory(target_minus_source);
        for (int mag = 0; mag < 4*(nullptr != Veff); ++mag) {
            free_memory(Veff[mag]);
        } // mag
        free_memory(Veff);
        free_memory(veff_index);
        free_memory(colCubePos);
        free_memory(rowCubePos);
        free_memory(grid_spacing_trunc);
        free_memory(phase);
        green_debug_printf("# %s sizeof(plan_t) = %ld Byte\n", __func__, sizeof(action_plan_t));
    } // destructor

namespace action_plan {

#ifdef    NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    status_t test_construction_and_destrunction(int const echo=3) {
        if (echo > 3) std::printf("# sizeof(action_plan::atom_t)= %ld Byte\n", sizeof(atom_t));
        action_plan_t ap;
        return 0;
    } // test_init

    status_t all_tests(int const echo) {
        status_t stat(0);
        stat += test_construction_and_destrunction(echo);
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace action_plan
