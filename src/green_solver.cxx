// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <complex> // std::complex<real_t>

#include "green_solver.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "action_plan.hxx" // action_plan_t
#include "control.hxx" // ::get
#include "green_action.hxx" // action_t<real_t,R1C2,Noco>
#include "recorded_warnings.hxx" // error

// #define   DEBUG

#ifdef    DEBUG
  #define green_debug_printf(...) { std::printf(__VA_ARGS__); std::fflush(stdout); }
#else  // DEBUG
  #define green_debug_printf(...)
#endif // DEBUG

typedef green_action::action_t<float ,2,1> Act421;
typedef green_action::action_t<float ,2,2> Act422;
typedef green_action::action_t<double,2,1> Act821;
typedef green_action::action_t<double,2,2> Act822;


    green_solver_t::green_solver_t(action_plan_t* p, int const echo, int const check) {
        if (nullptr != p) {
            if (echo > 0) std::printf("# construct %s +check=%d\n", __func__, check);
            int const fp_input = control::get("green_solver.floating.point.bits", 32.);
            int const fp = (32 == fp_input) ? 32 : 64;
            if (echo > 0) std::printf("# +green_solver.floating.point.bits=%i --> %i\n", fp_input, fp);
            int constexpr r1c2 = 2; // 1:real, 2:complex (always complex since tfQMRgpu does not support real)
            int const noco_input = control::get("green_solver.noco", 0.);
            int const noco = (2 == noco_input) ? 2 : 1;
            if (echo > 0) std::printf("# +green_solver.noco=%i --> Noco=%i\n", noco_input, noco);
            action_key_ = 1000*fp + 10*r1c2 + noco;

            if (echo > 0) std::printf("# green_solver_t::action_key= %i\n", int(action_key_));
            // initialize
            switch (action_key_) {
            case 32021: action_ = (void*)new Act421(p, echo, check); break; // complex
            case 32022: action_ = (void*)new Act422(p, echo, check); break; // complex non-collinear
            case 64021: action_ = (void*)new Act821(p, echo, check); break; // double complex
            case 64022: action_ = (void*)new Act822(p, echo, check); break; // double complex non-collinear
            default: error("No such action_key= %i", int(action_key_));
            } // switch action_key_
            assert(action_ && "action pointer must be valid");
#ifdef    DEBUGGPU
            if (echo > 0) std::printf("# %s.action= %p\n", __func__, action_);
#endif // DEBUGGPU
            if (echo > 0) std::printf("# construct %s --> needs %.6f GByte GPU memory\n", __func__, p->gpu_mem*1e-9);
        } else {
            if (echo > 0) std::printf("# cannot construct %s with p=nullptr\n", __func__);
        }
    } // constructor

    green_solver_t::~green_solver_t() { // destructor
#ifdef    DEBUGGPU
        std::printf("# destruct %s, this=%p\n", __func__, (void*)this);
#endif // DEBUGGPU        
        green_debug_printf("# destruct %s, action_key= %i action_ptr=%p\n", __func__, int(action_key_), action_);
        if (action_) {
            green_debug_printf("# destruct %s, action_key= %i\n", __func__, int(action_key_));
            switch (action_key_) {
            case 32021: { ((Act421*)action_)->~action_t(); } break; // complex
            case 32022: { ((Act422*)action_)->~action_t(); } break; // complex non-collinear
            case 64021: { ((Act821*)action_)->~action_t(); } break; // double complex
            case 64022: { ((Act822*)action_)->~action_t(); } break; // double complex non-collinear
            default: error("action=%p but action_key=%d", action_, action_key_);
            } // switch action_key_
        } // action_
    } // destructor

    status_t green_solver_t::solve(
          std::complex<double> rho[] // result: density [plan.nCols][4*4*4]
        , uint32_t const nblocks
        , int const max_iterations
        , int const imag // =1 // index of the exported part 1:imaginary part, 0:real part
        , int const echo // =0 // verbosity
    ) {
        if (echo > 7) std::printf("# green_solver_t::solve with action_key= %i, echo= %d\n", int(action_key_), echo);
#ifdef    DEBUGGPU
        if (echo > 3) std::printf("# %s.action= %p\n", __func__, action_);
#endif // DEBUGGPU
        assert(action_ && "action pointer must be valid for call to solve");
        switch (action_key_) {
        case 32021: return ((Act421*)action_)->solve(rho, nblocks, max_iterations, echo); // complex
        case 32022: return ((Act422*)action_)->solve(rho, nblocks, max_iterations, echo); // complex non-collinear
        case 64021: return ((Act821*)action_)->solve(rho, nblocks, max_iterations, echo); // double complex
        case 64022: return ((Act822*)action_)->solve(rho, nblocks, max_iterations, echo); // double complex non-collinear
        default: error("No solve with such action_key= %i", int(action_key_)); return action_key_;
        } // switch action_key_
    } // solve













namespace green_solver {

#ifdef    NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    status_t test_construction_and_destruction(int const echo) {
        status_t stat(0);
        green_solver_t gs0;
        green_solver_t gs1(nullptr, echo);
        return stat;
    } // test_construction_and_destruction

    status_t all_tests(int const echo) {
        status_t stat(0);
        stat += test_construction_and_destruction(echo);
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_solver
