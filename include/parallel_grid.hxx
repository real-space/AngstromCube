#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // uint32_t, int8_t
#include <algorithm> // std::min, ::max
#include <cmath> // std::floor, ::ceil, ::sqrt, ::abs
#include <cassert> // assert

#include "status.hxx" // status_t
#include "real_space.hxx" // ::grid_t
#include "inline_math.hxx" // set, scale
// #include "display_units.h" // Ang, _Ang
#include "recorded_warnings.hxx" // warn
// #include "boundary_condition.hxx" // Periodic_Boundary, Isolated_Boundary, Mirrored_Boundary, Invalid_Boundary
#include "mpi_parallel.hxx" // ...

namespace parallel_grid {

  int constexpr debug = 0;
 
  template <typename rank_t=int16_t> // must be a signed integer type
  class distribution_t {
  public: 

      distribution_t(uint32_t const nx=0, uint32_t const ny=0, uint32_t const nz=0)
       : dims{nx, ny, nz}
      {
          std::printf("# distribution_t<int%ld_t>(%d, %d, %d)\n", sizeof(rank_t)*8, nx, ny, nz);
          assert(1 == (rank_t(2)/2)) && "distribution_t needs a signed INTEGER type as rank_t");
          assert(rank_t(-1) < 0      && "distribution_t needs a SIGNED integer type as rank_t");
      } // default constructor

      template <typename int_t>
      distribution_t(int_t const n[3]) : distribution_t(n[0], n[1], n[2]) {} // delegating constructor

      ~distribution_t() {
      } // destructor

  private:
      uint32_t dims[3]; // 0,1,2:real-space grid dimensions of grid points (or blocks of grid points)

  }; // class distribution_t



#ifdef NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_create_and_destroy(int const echo=9) {
      int const dims[] = {10, 20, 30};
      auto gp = new distribution_t<>(dims); // construct
      gp->~distribution_t(); // explicity call the destructor
      distribution_t<float> fails;
      return 0;
  } // test_create_and_destroy

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_create_and_destroy(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace parallel_grid
