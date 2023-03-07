set(a43_version 0.0.1)

option(HAS_MKL "Usage of the Intel Math Kernel Library" OFF)


set(a43_atom_headers
        ../include/status.hxx
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

if(HAS_MKL)
    # not tested!
    find_package(MKL REQUIRED)
    target_link_libraries(libliveatom PUBLIC mkl)
else(HAS_MKL)
    add_definitions(-DHAS_NO_MKL)
    target_link_libraries(libliveatom PUBLIC blas lapack)
endif(HAS_MKL)

target_include_directories(libliveatom
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
   )

target_link_libraries(libliveatom PUBLIC)

# set(a43_targets a43)
# set(a43_public_headers 
#    ) # empty list


# set(a43_headers ${a43_atom_headers} ${a43_public_headers})

# set(a43_sources
#         geometry_analysis.cxx
#         self_consistency.cxx
#         iterative_poisson.cxx
#         bessel_transform.cxx
#         dense_solver.cxx
#         multi_grid.cxx
#         mpi_parallel.cxx
#         potential_generator.cxx
#         poisson_solver.cxx
#         plane_wave.cxx
#         sho_basis.cxx
#         sho_overlap.cxx
#         sho_hamiltonian.cxx
#         sho_projection.cxx
#         sho_potential.cxx
#         structure_solver.cxx
#         green_parallel.cxx
#         green_function.cxx
#         green_experiments.cxx
#         green_input.cxx
#         load_balancer.cxx
#         green_tests.cxx
#         main.cxx
# )

# add_executable(a43 main.cxx ${a43_headers})

# set(a43_targets ${a43_targets} a43)
