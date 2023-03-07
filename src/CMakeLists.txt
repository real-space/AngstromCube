### 
### CMake file for ./a43 and libliveatom.dylib
### 

set(a43_version 0.0.1)

set(a43_targets libliveatom)

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
add_library(libliveatom SHARED ${a43_atom_sources} ${a43_atom_headers})


################
################
################



set(a43_targets ${a43_targets} a43)

set(a43_headers ${a43_atom_headers}
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
        ../include/green_parallel.hxx
        ../include/green_function.hxx
        ../include/green_experiments.hxx
        ../include/green_input.hxx
        ../include/load_balancer.hxx
        ../include/green_tests.hxx
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
        green_function.cxx
        green_experiments.cxx
        green_input.cxx
        load_balancer.cxx
        green_tests.cxx
        main.cxx
   )

add_executable(a43 ${a43_sources} ${a43_headers})


if(HAS_DEVEL)
else(HAS_DEVEL)
    add_definitions(-DDEVEL)
endif(HAS_DEVEL)

if(HAS_MPI)
    # not tested!
    find_package(MPI REQUIRED)
    target_link_libraries(a43 PUBLIC mpi)
else(HAS_MPI)
    add_definitions(-DHAS_NO_MPI)
endif(HAS_MPI)

if(HAS_CUDA)
    # not tested!
    find_package(CUDA REQUIRED)
    # find_package(CUDAToolkit REQUIRED) ### for cuRAND
    target_link_libraries(a43 PUBLIC cuda)
    # CUDA architectures 70=Volta, 75=Turing, 80=Ampere
    set_property(TARGET a43 PROPERTY CUDA_ARCHITECTURES 70)
else(HAS_CUDA)
    add_definitions(-DHAS_NO_CUDA)
endif(HAS_CUDA)

if(HAS_MKL)
    # not tested!
    find_package(MKL REQUIRED)
    target_link_libraries(libliveatom PUBLIC mkl)
    target_link_libraries(a43 PUBLIC mkl)
else(HAS_MKL)
    add_definitions(-DHAS_NO_MKL)
    target_link_libraries(libliveatom PUBLIC blas lapack)
    target_link_libraries(a43 PUBLIC blas lapack)
endif(HAS_MKL)

if(HAS_UNIT_TESTS)
    # not tested!
else(HAS_UNIT_TESTS)
    add_definitions(-DNO_UNIT_TESTS)
endif(HAS_UNIT_TESTS)


target_include_directories(a43
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
   )

target_include_directories(libliveatom
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
   )
