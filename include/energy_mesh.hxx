#pragma once
// This file is part of AngstromCube under MIT License

#include <vector> // std::vector<T>

#include "status.hxx" // status_t

namespace energy_mesh {

    typedef std::complex<double> Complex;

    std::vector<Complex> get(std::vector<Complex> & w8, int const echo=0); // declaration only

    status_t all_tests(int const echo=0); // declaration only

} // namespace energy_mesh
