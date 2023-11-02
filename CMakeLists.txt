#
# CMake for AngstromCube
#
cmake_minimum_required(VERSION 3.24)

project(AngstromCube VERSION 0.9)

# Options
option(HAS_MKL          "Using the Intel Math Kernel Library MKL"       ON )
option(HAS_UNIT_TESTS   "All modules offer self-tests (recommended)"    ON )
option(HAS_MPI          "Using the Message Passing Interface MPI"       OFF)
option(HAS_DEVEL        "Include experimental features"                 ON )
option(HAS_TFQMRGPU     "Using the header-only version of tfQMRgpu"     OFF)
option(HAS_RAPIDXML     "Using the 3rd party XML parser RapidXML"       ON )
option(HAS_RAPIDJSON    "Using the 3rd party JSON parser RapidJSON"     OFF) ### errors
option(HAS_FFTW         "Using the Fast Fourier Transform FFTW"         OFF) ### does not work yet
option(HAS_OPENMP       "Using OpenMP to accelerate a43"                OFF)

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

# help debugging
set(CMAKE_VERBOSE_MAKEFILE ON)

if(HAS_MKL)
    # add MKL (Intel Math Kernel Library)
    unset(DPCPP_COMPILER) ## do I need this?
    find_package(MKL REQUIRED)
    include(CMakeFindDependencyMacro)
    # threaded MKL?
    find_dependency(OpenMP 3.0 REQUIRED)
endif(HAS_MKL)

# default installation paths
include(GNUInstallDirs)

# external dependencies
add_subdirectory(external)

# serial executable ./a43 and ./libliveatom.so
add_subdirectory(src)

# parallel executable ./green
add_subdirectory(green)

