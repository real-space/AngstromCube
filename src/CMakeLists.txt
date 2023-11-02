### 
### CMake file for ./a43 and libliveatom.so and ./green
### 

set(a43_version 0.0.1)

set(all_targets )
set(all_targets ${all_targets} liveatom)

set(a43_atom_headers
    ../include/control.hxx
    ../include/recorded_warnings.hxx
    ../include/radial_grid.hxx
    ../include/radial_integrator.hxx
    ../include/radial_potential.hxx
    ../include/radial_eigensolver.hxx
    ../include/exchange_correlation.hxx
    ../include/scattering_test.hxx
    ../include/chemical_symbol.hxx
    ../include/display_units.h
    ../include/sigma_config.hxx
    ../include/atom_core.hxx
    ../include/sho_unitary.hxx
    ../include/angular_grid.hxx
    ../include/single_atom.hxx
    ../include/single_atom.h
)

set(a43_atom_sources
    control.cxx
    recorded_warnings.cxx
    radial_grid.cxx
    radial_integrator.cxx
    radial_potential.cxx
    radial_eigensolver.cxx
    exchange_correlation.cxx
    scattering_test.cxx
    chemical_symbol.cxx
    display_units.cxx
    sigma_config.cxx
    atom_core.cxx
    sho_unitary.cxx
    angular_grid.cxx
    single_atom.cxx
)
   
# build a shared library for the liveatom
add_library(liveatom SHARED ${a43_atom_sources} ${a43_atom_headers})


################
################
################



set(all_targets ${all_targets} a43)

set(a43_headers ${a43_atom_headers}
    ../include/simple_timer.hxx
    ../include/global_coordinates.hxx
    ../include/geometry_analysis.hxx
    ../include/self_consistency.hxx
    ../include/iterative_poisson.hxx
    ../include/bessel_transform.hxx
    ../include/dense_solver.hxx
    ../include/multi_grid.hxx
    ../include/mpi_parallel.hxx
    ../include/potential_generator.hxx
    ../include/poisson_solver.hxx
    ../include/plane_wave.hxx
    ../include/sho_basis.hxx
    ../include/sho_overlap.hxx
    ../include/sho_hamiltonian.hxx
    ../include/sho_projection.hxx
    ../include/sho_potential.hxx
    ../include/structure_solver.hxx
    ../include/green_projection.hxx
    ../include/green_parallel.hxx
    ../include/green_kinetic.hxx
    ../include/green_function.hxx
    ../include/green_experiments.hxx
    ../include/green_input.hxx
    ../include/load_balancer.hxx
    ../include/green_tests.hxx
    ../include/add_module_test.h
)

set(a43_sources ${a43_atom_sources}
    geometry_analysis.cxx
    self_consistency.cxx
    iterative_poisson.cxx
    bessel_transform.cxx
    dense_solver.cxx
    multi_grid.cxx
    mpi_parallel.cxx
    potential_generator.cxx
    poisson_solver.cxx
    plane_wave.cxx
    sho_basis.cxx
    sho_overlap.cxx
    sho_hamiltonian.cxx
    sho_projection.cxx
    sho_potential.cxx
    structure_solver.cxx
    green_parallel.cxx
    green_kinetic.cxx
    green_function.cxx
    green_experiments.cxx
    green_input.cxx
    load_balancer.cxx
    green_tests.cxx
    a43_main.cxx
)

# build a serial application
add_executable(a43 ${a43_sources} ${a43_headers})

################
################
################

set(all_targets ${all_targets} green)

set(green_headers
    ../include/control.hxx
    ../include/recorded_warnings.hxx
    ../include/display_units.h
    ../include/mpi_parallel.hxx
    ../include/simple_timer.hxx
    ../include/global_coordinates.hxx
    ../include/green_projection.hxx
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
    set(green_sources ${green_sources}
        green_kinetic.cxx
        green_function.cxx
        green_experiments.cxx
        green_tests.cxx
    )
endif(HAS_CUDA)

set(green_sources ${green_sources}
    green_main.cxx
)

# build a parallel GPU application
add_executable(green ${green_sources} ${green_headers})





# version control using git
execute_process(COMMAND git log -1 --pretty=format:"%H"  OUTPUT_VARIABLE GIT_KEY)
add_definitions(-D_GIT_KEY=${GIT_KEY})


if(HAS_CUDA)
    find_package(CUDA REQUIRED)
    find_package(CUDAToolkit REQUIRED) ### for cuRAND
    target_compile_options(green PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
    target_compile_options(green PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-G -lineinfo -Xptxas -v>)
    target_link_libraries(green PUBLIC cuda)
    # CUDA architectures 70=Volta, 80=Ampere, 90=GraceHopper
    set_property(TARGET green PROPERTY CUDA_ARCHITECTURES 70 80)
else(HAS_CUDA)
    target_compile_definitions(green PUBLIC -DHAS_NO_CUDA)
endif(HAS_CUDA)
# a43 never has CUDA
target_compile_definitions(a43 PUBLIC -DHAS_NO_CUDA)


if(HAS_DEVEL)
    add_definitions(-DDEVEL)
else(HAS_DEVEL)
endif(HAS_DEVEL)

if(HAS_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(a43 PUBLIC OpenMP::OpenMP_CXX)
    endif()
endif(HAS_OPENMP)

if(HAS_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(a43   PUBLIC mpi)
    target_link_libraries(green PUBLIC mpi)
else(HAS_MPI)
    target_compile_definitions(a43   PUBLIC -DHAS_NO_MPI)
    target_compile_definitions(green PUBLIC -DHAS_NO_MPI)
endif(HAS_MPI)



if(MKL_FOUND)
    # add MKL (Intel Math Kernel Library)
    target_link_libraries(liveatom PUBLIC iomp5 mkl)
    target_link_libraries(a43      PUBLIC iomp5 mkl)
    target_link_libraries(green    PUBLIC iomp5 mkl)
else(MKL_FOUND)
    # MKL has not been found
    target_compile_definitions(liveatom PUBLIC -DHAS_NO_MKL)
    target_compile_definitions(a43      PUBLIC -DHAS_NO_MKL)
    target_compile_definitions(green    PUBLIC -DHAS_NO_MKL)
    # add BLAS
    find_package(BLAS REQUIRED)
    target_link_libraries(liveatom PUBLIC blas)
    target_link_libraries(a43      PUBLIC blas)
    target_link_libraries(green    PUBLIC blas)
    # add LAPACK
    find_package(LAPACK REQUIRED)
    target_link_libraries(liveatom PUBLIC lapack)
    target_link_libraries(a43      PUBLIC lapack)
    target_link_libraries(green    PUBLIC lapack)
    # add FFTW
    if(HAS_FFTW)
        find_package(fftw REQUIRED)
        target_compile_definitions(a43 PUBLIC -DHAS_FFTW)
        target_link_libraries(a43 PUBLIC fftw)
    endif(HAS_FFTW)
endif(MKL_FOUND)

if(HAS_UNIT_TESTS)
    # usual case
else(HAS_UNIT_TESTS)
    add_definitions(-DNO_UNIT_TESTS)
endif(HAS_UNIT_TESTS)

if(HAS_RAPIDXML)
    target_compile_definitions(a43   PUBLIC -DHAS_RAPIDXML)
    target_compile_definitions(green PUBLIC -DHAS_RAPIDXML)
    target_include_directories(a43   PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_deps/rapidxml-src/RapidXML>)
    target_include_directories(green PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_deps/rapidxml-src/RapidXML>)
endif(HAS_RAPIDXML)

if(HAS_RAPIDJSON)
    target_compile_definitions(a43   PUBLIC -DHAS_RAPIDJSON)
    target_compile_definitions(green PUBLIC -DHAS_RAPIDJSON)
    target_include_directories(green PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_deps/rapidjson-src>)
    target_include_directories(a43   PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_deps/rapidjson-src>)
endif(HAS_RAPIDJSON)


target_include_directories(liveatom PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_include_directories(a43      PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_include_directories(green    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# install the 3 binaries
install(TARGETS ${all_targets}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
