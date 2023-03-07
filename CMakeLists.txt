#
# CMake for A43
#
cmake_minimum_required(VERSION 3.24)

project(AngstromCube VERSION 0.9)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

if(APPLE)
  set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
elseif(UNIX)
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

add_subdirectory(src)
# add_subdirectory(green)

