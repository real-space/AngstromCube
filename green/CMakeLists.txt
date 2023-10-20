### 
### CMake file for ./green
### 

set(green_version 0.0.1)

set(green_targets green)

set(green_headers
    ../include/control.hxx
    ../include/recorded_warnings.hxx
    ../include/display_units.h
    ../include/mpi_parallel.hxx
    ../include/green_parallel.hxx
    ../include/green_input.hxx
    ../include/load_balancer.hxx
    ../include/green_kinetic.hxx
    ../include/green_function.hxx
    ../include/green_experiments.hxx
    ../include/green_tests.hxx
)

set(green_sources
    control.cxx
    recorded_warnings.cxx
    display_units.cxx
    mpi_parallel.cxx
    green_parallel.cxx
    green_input.cxx
    load_balancer.cxx
)

if(HAS_CUDA)
    enable_language(CUDA)
    set(green_sources ${green_sources}
        green_kinetic.cu
        green_function.cu
        green_experiments.cu
        green_tests.cu
    )
else(HAS_CUDA)
    add_definitions(-DHAS_NO_CUDA)
    set(green_sources ${green_sources}
        green_kinetic.cxx
        green_function.cxx
        green_experiments.cxx
        green_tests.cxx
    )
endif(HAS_CUDA)

set(green_sources ${green_sources}
    green.cxx
)


add_executable(green ${green_sources} ${green_headers})

if(HAS_CUDA)
    target_compile_options(green
        PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
    )
    find_package(CUDA REQUIRED)
    find_package(CUDAToolkit REQUIRED) ### for cuRAND
    target_link_libraries(green PUBLIC cuda)
    # CUDA architectures 70=Volta, 75=Turing, 80=Ampere
    set_property(TARGET green PROPERTY CUDA_ARCHITECTURES 70 80)
endif(HAS_CUDA)



if(HAS_DEVEL)
    add_definitions(-DDEVEL)
else(HAS_DEVEL)
endif(HAS_DEVEL)

if(HAS_MPI)
    # not tested!
    find_package(MPI REQUIRED)
    target_link_libraries(green PUBLIC mpi)
else(HAS_MPI)
    add_definitions(-DHAS_NO_MPI)
endif(HAS_MPI)

if(MKL_FOUND)
    if(OpenMP_FOUND)
        # target_link_libraries(green PUBLIC openmp)
    endif()
    target_link_libraries(green PUBLIC iomp5 mkl)
else(MKL_FOUND)
    add_definitions(-DHAS_NO_MKL)
    # add BLAS
    find_package(BLAS REQUIRED)
    target_link_libraries(green PUBLIC blas)
    # add LAPACK
    find_package(LAPACK REQUIRED)
    target_link_libraries(green PUBLIC lapack)
endif(MKL_FOUND)

if(HAS_UNIT_TESTS)
    # not tested!
else(HAS_UNIT_TESTS)
    add_definitions(-DNO_UNIT_TESTS)
endif(HAS_UNIT_TESTS)

if(HAS_TFQMRGPU)
    # not tested!
    add_definitions(-DHAS_TFQMRGPU)
else(HAS_TFQMRGPU)
endif(HAS_TFQMRGPU)

if(HAS_RAPIDXML)
    add_definitions(-DHAS_RAPIDXML)
    message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
    target_include_directories(green PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_deps/rapidxml-src>)
endif(HAS_RAPIDXML)

if(HAS_RAPIDJSON)
    add_definitions(-DHAS_RAPIDJSON)
    target_include_directories(green PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_deps/rapidjson-src>)
endif(HAS_RAPIDJSON)

target_include_directories(green PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)


# install green
set(green_targets
    green
)

install(TARGETS ${green_targets}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
