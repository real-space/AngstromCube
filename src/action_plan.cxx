// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert

#include "action_plan.hxx"

#include "green_memory.hxx" // free_memory

#ifdef    DEBUG
  #define green_debug_printf(...) { std::printf(__VA_ARGS__); std::fflush(stdout); }
#else  // DEBUG
  #define green_debug_printf(...)
#endif // DEBUG

    action_plan_t::action_plan_t() { // constructor 
        green_debug_printf("# default constructor for %s\n", __func__);
        // please see construct_Green_function in green_function.hxx for the construction of the plan_t
    } // constructor

    action_plan_t::~action_plan_t() { // destructor
        green_debug_printf("# destruct %s\n", __func__);
        free_memory(RowStart);
        free_memory(rowindx);
        free_memory(source_coords);
        free_memory(target_coords);
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
//             std::printf("# %s sizeof(atom_t) = %ld Byte\n", __func__, sizeof(atom_t));

    } // destructor

