#
# CMake for A43
#
cmake_minimum_required(VERSION 3.24)

project(AngstromCube VERSION 0.9)

# Options
option(HAS_MKL          "Using the Intel Math Kernel Library MKL"       OFF)
option(HAS_UNIT_TESTS   "All modules offer self-tests (recommended)"    ON )
option(HAS_MPI          "Using the Message Passing Interface MPI"       OFF)
option(HAS_DEVEL        "Include experimental features"                 OFF)
option(HAS_TFQMRGPU     "Using the header-only version of tfQMRgpu"     OFF)
option(HAS_RAPIDXML     "Using the 3rd party XML parser RapidXML"       OFF) ### does not work yet
option(HAS_RAPIDJSON    "Using the 3rd party JSON parser RapidJSON"     OFF) ### does not work yet
option(HAS_FFTW         "Using the Fast Fourier Transform FFTW"         OFF) ### does not work yet

find_package(CUDA)
if(CUDA_FOUND)
    option(HAS_CUDA         "Using the CUDA programming model for GPUs"    ON )
else(CUDA_FOUND)
    option(HAS_CUDA         "Using the CUDA programming model for GPUs"    OFF)
endif(CUDA_FOUND)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(APPLE)
  set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
elseif(UNIX)
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

# serial executable ./a43 and ./libliveatom.dylib
add_subdirectory(src)

# parallel executable ./green
add_subdirectory(green)
