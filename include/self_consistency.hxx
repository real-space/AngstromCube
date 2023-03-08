#pragma once
// This file is part of AngstromCube under MIT License

#include "status.hxx" // status_t

namespace self_consistency {

  status_t init(int const echo=0, float const ion=0.f); // declaration only

  status_t all_tests(int const echo=0); // declaration only

} // namespace self_consistency
