#
# CMake for A43
#
cmake_minimum_required(VERSION 3.24)

project(AngstromCube VERSION 0.9)

# Options
option(HAS_MKL          "Using the Intel Math Kernel Library MKL"       OFF)
option(HAS_UNIT_TESTS   "All modules offer self-tests (recommended)"    ON )
option(HAS_MPI          "Using the Message Passing Interface MPI"       OFF)
option(HAS_CUDA         "Using the CUDA programming model for GPUs"     OFF)
option(HAS_DEVEL        "Include experimental features"                 OFF)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(APPLE)
  set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
elseif(UNIX)
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

add_subdirectory(src)
# add_subdirectory(green)
