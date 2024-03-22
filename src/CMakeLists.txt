### 
### CMake file for ./a43 and libliveatom.so and ./green
### 

set(a43_version 0.0.1)

set(all_targets )
set(all_targets ${all_targets} liveatom)

set(a43_atom_headers
    ../include/angular_grid.hxx
    ../include/atom_core.hxx
    ../include/chemical_symbol.h
    ../include/chemical_symbol.hxx
    ../include/control.hxx
    ../include/display_units.h
    ../include/exchange_correlation.hxx
    ../include/quantum_numbers.h
    ../include/radial_eigensolver.hxx
    ../include/radial_grid.h
    ../include/radial_grid.hxx
    ../include/radial_integrator.hxx
    ../include/radial_potential.hxx
    ../include/recorded_warnings.hxx
    ../include/scattering_test.hxx
    ../include/sho_unitary.hxx
    ../include/sigma_config.hxx
    ../include/single_atom.h
    ../include/single_atom.hxx
)

set(a43_atom_sources
    angular_grid.cxx
    atom_core.cxx
    chemical_symbol.cxx
    control.cxx
    display_units.cxx
    exchange_correlation.cxx
    library_kind_dynamic.cxx
    radial_eigensolver.cxx
    radial_grid.cxx
    radial_integrator.cxx
    radial_potential.cxx
    recorded_warnings.cxx
    scattering_test.cxx
    sho_unitary.cxx
    sigma_config.cxx
    single_atom.cxx
)
   
# build a shared library for the liveatom
add_library(liveatom SHARED ${a43_atom_sources} ${a43_atom_headers})


################
################
################



set(all_targets ${all_targets} a43)

set(a43_headers ${a43_atom_headers}
    ../include/add_module_test.h
    ../include/bessel_transform.hxx
    ../include/dense_solver.hxx
    ../include/energy_contour.hxx
    ../include/energy_contribution.hxx
    ../include/gaunt_entry.h
    ../include/geometry_analysis.hxx
    ../include/geometry_input.hxx
    ../include/global_coordinates.hxx
    ../include/green_experiments.hxx
    ../include/green_function.hxx
    ../include/green_input.hxx
    ../include/green_kinetic.hxx
    ../include/green_parallel.hxx
    ../include/green_tests.hxx
    ../include/iterative_poisson.hxx
    ../include/load_balancer.hxx
    ../include/mpi_parallel.hxx
    ../include/multi_grid.hxx
    ../include/parallel_poisson.hxx
    ../include/parallel_potential.hxx
    ../include/plane_wave.hxx
    ../include/poisson_solver.hxx
    ../include/potential_generator.hxx
    ../include/self_consistency.hxx
    ../include/sho_basis.hxx
    ../include/sho_hamiltonian.hxx
    ../include/sho_overlap.hxx
    ../include/sho_potential.hxx
    ../include/sho_projection.hxx
    ../include/simple_timer.hxx
    ../include/structure_solver.hxx
)

set(a43_sources ${a43_atom_sources}
    a43_main.cxx
    bessel_transform.cxx
    dense_solver.cxx
    energy_contour.cxx
    geometry_analysis.cxx
    geometry_input.cxx
    green_experiments.cxx
    green_function.cxx
    green_input.cxx
    green_kinetic.cxx
    green_dyadic.cxx
    green_parallel.cxx
    green_potential.cxx
    green_action.cxx
    green_tests.cxx
    iterative_poisson.cxx
    load_balancer.cxx
    mpi_parallel.cxx
    multi_grid.cxx
    parallel_poisson.cxx
    parallel_potential.cxx
    dyadic_plan.cxx
    kinetic_plan.cxx
    action_plan.cxx
    green_memory.cxx
    green_solver.cxx
	green_sparse.cxx
    plane_wave.cxx
    poisson_solver.cxx
    potential_generator.cxx
    self_consistency.cxx
    sho_basis.cxx
    sho_hamiltonian.cxx
    sho_overlap.cxx
    sho_potential.cxx
    sho_projection.cxx
    structure_solver.cxx
)

# build a serial application
add_executable(a43 ${a43_sources} ${a43_headers})

################
################
################

set(all_targets ${all_targets} green)

set(green_headers
    ../include/boundary_condition.hxx
    ../include/chemical_symbol.h
    ../include/chemical_symbol.hxx
    ../include/constants.hxx
    ../include/control.hxx
    ../include/data_list.hxx
    ../include/data_view.hxx
    ../include/display_units.h
    ../include/energy_contour.hxx
    ../include/exchange_correlation.hxx
    ../include/geometry_input.hxx
    ../include/global_coordinates.hxx
    ../include/green_experiments.hxx
    ../include/green_function.hxx
    ../include/green_input.hxx
    ../include/green_kinetic.hxx
    ../include/green_parallel.hxx
    ../include/green_sparse.hxx
    ../include/dyadic_plan.hxx
    ../include/green_tests.hxx
    ../include/hermite_polynomial.hxx
    ../include/inline_math.hxx
    ../include/load_balancer.hxx
    ../include/mpi_parallel.hxx
    ../include/parallel_poisson.hxx
    ../include/parallel_potential.hxx
    ../include/print_tools.hxx
    ../include/recorded_warnings.hxx
    ../include/sho_projection.hxx
    ../include/sho_tools.hxx
    ../include/sho_unitary.hxx
    ../include/simple_stats.hxx
    ../include/simple_timer.hxx
    ../include/single_atom.h
)

set(green_sources
    action_plan.cxx
    chemical_symbol.cxx
    control.cxx
    display_units.cxx
    dyadic_plan.cxx
    energy_contour.cxx
    exchange_correlation.cxx
    geometry_input.cxx
    green_function.cxx
    green_input.cxx
    green_parallel.cxx
	green_sparse.cxx
    kinetic_plan.cxx
    load_balancer.cxx
    mpi_parallel.cxx
    parallel_poisson.cxx
    parallel_potential.cxx
    recorded_warnings.cxx
    sho_unitary.cxx
)


if(HAS_CUDA)
    enable_language(CUDA)
    set(green_sources ${green_sources}
        green_action.cu
        green_dyadic.cu
        green_experiments.cu
        green_kinetic.cu
        green_memory.cu
        green_potential.cu
        green_solver.cu
        green_tests.cu
    )
else(HAS_CUDA)
    set(green_sources ${green_sources}
        green_action.cxx
        green_dyadic.cxx
        green_experiments.cxx
        green_kinetic.cxx
        green_memory.cxx
        green_potential.cxx
        green_solver.cxx
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
    # GPU debugging information:
    # target_compile_options(green PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-G -lineinfo -Xptxas -v>)
    target_link_libraries(green PUBLIC cuda curand)
    # CUDA architectures 70=Volta, 80=Ampere, 90=Hopper
    set_property(TARGET green PROPERTY CUDA_ARCHITECTURES 70 80)
else(HAS_CUDA)
    target_compile_definitions(green PUBLIC -DHAS_NO_CUDA)
endif(HAS_CUDA)
# a43 never has CUDA
target_compile_definitions(a43 PUBLIC -DHAS_NO_CUDA)

# a43 comes together with single_atom, while green will be linked against libliveatom.so
target_compile_definitions(a43   PUBLIC -DHAS_SINGLE_ATOM)
target_compile_definitions(green PUBLIC -DHAS_LIVE_ATOM)
target_link_libraries(green PUBLIC liveatom)

if(HAS_DEVEL)
    add_definitions(-DDEVEL)
else(HAS_DEVEL)
endif(HAS_DEVEL)

if(HAS_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(a43 PUBLIC OpenMP::OpenMP_CXX)
    endif()
    target_compile_definitions(liveatom PUBLIC -DHAS_NO_OMP)
    target_compile_definitions(green    PUBLIC -DHAS_NO_OMP)
else(HAS_OPENMP)
    add_definitions(-DHAS_NO_OMP)
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

if(HAS_TFQMRGPU)
    # tfQMRgpu is a header-only library
    target_compile_definitions(green PUBLIC -DHAS_TFQMRGPU)
    target_include_directories(green PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_deps/tfqmrgpu-src/tfQMRgpu/include>)
endif(HAS_TFQMRGPU)

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
