#ifndef SINGLE_ATOM_HEADER
#define SINGLE_ATOM_HEADER

// C - interface

/*
    What do we need to replace PAWs in juRS?
    Input to single_atom:
        - Z_core, double
        - atom_id, int32_t
        - sigma, double
        - rcut, double
        - nn[8], vector<int8_t> ---> numax derived from nn
        - xc-type, ... ToDo
        optional
        - ionization, double
        - core_holes, ...
        every iteration:
        - vlm, vector<double>
        - dm, matrix<double> with am_stride

    Output of single_atom:
        - number_of_valence_electrons, double
        - partial core correction (or smooth core density), double[4096]
        - zero potential, radial spherical function, double[4096]
        - projectors ??
        - n_valence_e, double
        - qlm, vector<double>
        - lmax_qlm, int
        - lmax_vlm, int
        - stride, int
        - aHm, aSm, matrices<double> with am_stride
        - energy contributions, vector<double>, ...enums?
        - start waves occupation, double[6] (ssppdf), ToDo
        - start wave functions, reduced by r^ell, double[6][4096]

    In order to control soft switches (control::get()) in single_atom,
    we need a handle to set soft switches using control::set()

*/

#ifdef SINGLE_ATOM_SOURCE
    #ifndef __cplusplus
        #error "C-Interface and C++Implementations share a file!"
    #endif
    #include "control.hxx" // ::set, ::read_control_file
    #include "unit_system.hxx" // ::set
    #include "recorded_warnings.hxx" // warn, ::show_warnings
#endif

#define fortran_callable(NAME) void live_atom_##NAME##_

/*
    Example for inclusion in Fortran: 
        integer(kind=4), parameter :: na = 2
        integer(kind=4) :: status
        real(kind=8) :: quantity(na)
        external live_atom_get_some_quantity ! optional
        call live_atom_get_some_quantity(na, quantity, status)
        if (status /= 0) stop 'ERROR'
        write(*,*) quantity

    A note on equivalent types:
        integer(kind=4)       int32_t
        integer(kind=1)        int8_t
         real(kind=8)           double
         character(len=*)    char[] 
     Mind that strings are not null-terminated in Fortran
*/

    // creates name live_atom_init_env_ in single_atom.o
    fortran_callable(init_env)
        ( char const *filename
        , int32_t *status)
#ifdef SINGLE_ATOM_SOURCE
    {
        std::printf("\n# init environment variables for LiveAtoms from \'%s\'\n", filename); 
        int const echo = 1;
        *status = control::read_control_file(filename, echo);
        *status += unit_system::set(control::get("output.length.unit", "Bohr"),
                                    control::get("output.energy.unit", "Ha"), echo);
    } // live_atom_set_env_
#else
    ;
#endif

    fortran_callable(initialize)
        ( int32_t const *na // number of atoms
        , double const Z_core[] // core charges
        , int32_t const atom_id[] // no effect
        , int32_t const numax[] // result, SHO basis sizes
        , double const sigma[] // no effect
        , double const rcut[] // no effect
        , int8_t const nn[][8] // no effect
        , double const ionization[]
        , double const magnetization[] // layout [na], no effect
        , char const *xc_key // ToDo: exchange-correlation model
        , int32_t *stride // result
        , int32_t lmax_qlm[] // result, ToDo
        , int32_t lmax_vlm[] // result, ToDo
        , double n_valence_electrons[] // result
        , int32_t *status)
#ifdef SINGLE_ATOM_SOURCE
    {
        std::vector<double> sigma_cmp(*na);
        std::vector<float>        ion(*na);
        double constexpr convert_rcut_to_sigma_cmp = 0.22360679774997898; // == 1/sqrt(20.)
        for(int ia = 0; ia < *na; ++ia) {
            ion[ia] = ionization[ia];
            sigma_cmp[ia] = rcut[ia]*convert_rcut_to_sigma_cmp;
        } // ia
        *status = single_atom::atom_update("initialize", *na, (double*)Z_core, (int*)numax, ion.data());

        int const numax_max = 3; // ToDo: should be a result of atom_update
        *stride = align<2>(sho_tools::nSHO(numax_max));

        float avd{1};
        *status += single_atom::atom_update("lmax qlm",  *na,    nullptr, lmax_qlm, &avd);
        *status += single_atom::atom_update("lmax vlm",  *na, (double*)1, lmax_vlm);
        *status += single_atom::atom_update("sigma cmp", *na, sigma_cmp.data());
        *status += single_atom::atom_update("#valence electrons", *na, n_valence_electrons);

        warn("Initialized %d LiveAtoms, ToDo: need to deal with atom_id, nn, magnetization, xc_key", *na);
        // std::fflush(stdout);
    } // live_atom_initialize_
#else
    ;
#endif

    fortran_callable(set_env)
        ( char const *varname
        , char const *newvalue
        , int32_t *status)
#ifdef SINGLE_ATOM_SOURCE
    {
        std::printf("# set environment variable for LiveAtoms:\t%s = %s\n", 
                                                        varname, newvalue);
        control::set(varname, newvalue);
        *status = (nullptr == newvalue);
        // std::fflush(stdout);
    } // live_atom_set_env_
#else
    ;
#endif

    fortran_callable(get_core_density)
        ( int32_t const *na
        , double **rhoc // layout [na][4096]
        , int32_t *status)
#ifdef SINGLE_ATOM_SOURCE
    {
        *status = single_atom::atom_update("core densities", *na, 0, 0, 0, rhoc);
        std::printf("# got_core_density for %d LiveAtoms\n", *na);
        // std::fflush(stdout);
    } // live_atom_get_core_density_
#else
    ;
#endif

    fortran_callable(get_start_waves)
        ( int32_t const *na
        , double waves[]      // layout [na][6][4096]
        , double occupation[] // layout [na][6][2], can be magnetic
        , int32_t *status)
#ifdef SINGLE_ATOM_SOURCE
    {
        *status = single_atom::atom_update("?", *na);
        std::printf("# got_start_waves for %d LiveAtoms\n", *na);
        // std::fflush(stdout);
    } // live_atom_get_start_waves_
#else
    ;
#endif

    fortran_callable(set_density_matrix)
        ( int32_t const *na
        , double **atom_rho // layout [na][4096]
        , int32_t *status)
#ifdef SINGLE_ATOM_SOURCE
    {
        float mix_rho[3] = {0, 0, 0};
        *status = single_atom::atom_update("atomic density matrix", *na, 0, 0, mix_rho, atom_rho);
        std::printf("# density_matrices set for %d LiveAtoms\n", *na);
        // std::fflush(stdout);
    } // live_atom_set_density_matrix_
#else
    ;
#endif

    fortran_callable(get_compensation_charge)
        ( int32_t const *na
        , double **qlm // layout [na][(1 + maxval(lmax_qlm))^2]
        , int32_t *status)
#ifdef SINGLE_ATOM_SOURCE
    {
        *status = single_atom::atom_update("qlm charges", *na, 0, 0, 0, qlm);
        std::printf("# got_compensation_charge %d LiveAtoms\n", *na);
        // std::fflush(stdout);
    } // live_atom_get_compensation_charge_
#else
    ;
#endif

    fortran_callable(set_potential_multipole)
        ( int32_t const *na
        , double **vlm // layout [na][(1 + maxval(lmax_vlm))^2]
        , int32_t *status)
#ifdef SINGLE_ATOM_SOURCE
    {
        float mix_pot[1] = {0};
        *status = single_atom::atom_update("update", *na, 0, 0, mix_pot, vlm);
        std::printf("# potential_multipoles set for %d LiveAtoms\n", *na);
        // std::fflush(stdout);
    } // live_atom_set_potential_multipole_
#else
    ;
#endif

    fortran_callable(get_zero_potential)
        ( int32_t const *na
        , double **vbar // layout [na][4096]
        , int32_t *status)
#ifdef SINGLE_ATOM_SOURCE
    {
        *status = single_atom::atom_update("zero potential", *na, 0, 0, 0, vbar);
        std::printf("# got_zero_potential for %d LiveAtoms\n", *na);
        // std::fflush(stdout);
    } // live_atom_get_zero_potential_
#else
    ;
#endif

    fortran_callable(get_projectors)
        ( int32_t const *na
        , double sigma[]
        , int32_t numax[]
        , int32_t *status)
#ifdef SINGLE_ATOM_SOURCE
    {
        *status = single_atom::atom_update("projectors", *na, sigma, numax);
        std::printf("# got_projectors for %d LiveAtoms\n", *na);
        // std::fflush(stdout);
    } // live_atom_get_projectors_
#else
    ;
#endif

    fortran_callable(get_hamiltonian_matrix)
        ( int32_t const *na
        , double **hmt // logical layout [na][2][stride][stride]
        , int32_t *status)
#ifdef SINGLE_ATOM_SOURCE
    {
        *status = single_atom::atom_update("hamiltonian and overlap", *na, 0, 0, 0, hmt);
        std::printf("# got_hamiltonian_matrix for %d LiveAtoms\n", *na);
        // std::fflush(stdout);
    } // live_atom_get_hamiltonian_matrix_
#else
    ;
#endif

    fortran_callable(get_energy_contributions)
        ( int32_t const *na
        , double energies[]
        , int32_t *status)
#ifdef SINGLE_ATOM_SOURCE
    {
        *status = single_atom::atom_update("e", *na, energies);
        std::printf("# got_energy_contributions for %d LiveAtoms\n", *na);
        // std::fflush(stdout);
    } // live_atom_get_energy_contributions_
#else
    ;
#endif

    fortran_callable(finalize)
        ( int32_t const *na
        , int32_t *status)
#ifdef SINGLE_ATOM_SOURCE
    {
        *status = single_atom::atom_update("memory cleanup", *na);
        std::printf("# finalized %d LiveAtoms\n", *na);
        int const control_show = control::get("control.show", 0.);
        if (0 != control_show) control::show_variables(control_show);
        recorded_warnings::show_warnings(3);
        // std::fflush(stdout);
    } // live_atom_finalize_
#else
    ;
#endif

#endif // SINGLE_ATOM_HEADER (header guard)
