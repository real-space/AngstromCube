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



#define fortran_callable(NAME) void live_atom_##NAME##_

/*
	Example for inclusion in Fortran: 
		integer(kind=4), parameter :: na = 2
		integer(kind=4) :: status
		real(kind=8) :: quantity(na)
		external _live_atom_get_some_quantity ! optional
		call _live_atom_get_some_quantity(na, quantity, status)
		if (status /= 0) stop 'ERROR'
		write(*,*) quantity

	A note on equivalent types:
		integer(kind=4)   	int32_t
		integer(kind=1)		int8_t
 		real(kind=8)   		double
 		character(len=*)    char[] 
 	Mind that strings are not null-terminated in Fortran
*/

	// creates name _live_atom_initialize_ in single_atom.o
	fortran_callable(initialize)
        ( int32_t const *na // number of atoms
		, double const Z_core[] // core charges
		, int32_t const atom_id[] // no effect
		, int32_t const numax[] // SHO basis sizes
		, double const sigma[] // no effect
		, double const rcut[] // no effect
		, int8_t const nn[][8] // no effect
		, double const ionization[]
		, double const magnetization[] // layout [na], no effect
		, char const *xc_key // ToDo: exchange-correlation model
		, int32_t *stride // result
		, int32_t lmax_qlm[] // result
		, int32_t lmax_vlm[] // result
		, double  n_valence_e[] // result
		, int32_t *status)
#ifndef SINGLE_ATOM_HEADER_ONLY
	{
		std::vector<double> sigma_cmp(*na);
		std::vector<float> ion(*na);
		double const convert_rcut_to_sigma_cmp = 1./std::sqrt(20.);
		for(int ia = 0; ia < *na; ++ia) {
			ion[ia] = ionization[ia];
			sigma_cmp[ia] = rcut[ia]*convert_rcut_to_sigma_cmp;
		} // ia
		*status = single_atom::atom_update("initialize", *na, (double*)Z_core, (int*)numax, ion.data());

		int const numax_max = 3; // should come out of update
		*stride = align<2>(sho_tools::nSHO(numax_max));
		for(int ia = 0; ia < *na; ++ia) {
			lmax_qlm[ia] = numax_max;
			lmax_vlm[ia] = numax_max;
			n_valence_e[ia] = 0; // ToDo: get csv_charge[1] + csv_charge[2];
		} // ia
		printf("# Initialized %d LiveAtoms\n"
			"# ToDo: need to deal with atom_id, sigma, nn,\n"
			"# magnetization, xc_key and n_valence_e\n\n", *na);
		std::fflush(stdout);
	} // _live_atom_initialize_
#else
	;
#endif

	fortran_callable(set_env)
        ( char const *varname
		, char const *newvalue
		, int32_t *status)
#ifndef SINGLE_ATOM_HEADER_ONLY
	{
		printf("# set environment variable for LiveAtoms  %s = %s\n", 
												varname, newvalue);
		control::set(varname, newvalue);
		std::fflush(stdout);
	} // _live_atom_set_env_
#else
	;
#endif

	fortran_callable(get_core_density)
        ( int32_t const *na
		, double rhoc[] // layout [na][4096]
		, int32_t *status)
#ifndef SINGLE_ATOM_HEADER_ONLY
	{
		*status = single_atom::atom_update("?", *na);;
		printf("# got_core_density for %d LiveAtoms\n", *na);
		std::fflush(stdout);
	} // _live_atom_get_core_density_
#else
	;
#endif

	fortran_callable(get_start_waves)
        ( int32_t const *na
		, double waves[]      // layout [na][6][4096]
		, double occupation[] // layout [na][6][2], can be magnetic
		, int32_t *status)
#ifndef SINGLE_ATOM_HEADER_ONLY
	{
		*status = single_atom::atom_update("?", *na);;
		printf("# got_start_waves for %d LiveAtoms\n", *na);
		std::fflush(stdout);
	} // _live_atom_get_start_waves_
#else
	;
#endif

	fortran_callable(set_density_matrix)
        ( int32_t const *na
		, double const density_matrices[] // layout [na][stride][stride]
		, int32_t *status)
#ifndef SINGLE_ATOM_HEADER_ONLY
	{
		*status = single_atom::atom_update("d", *na);;
		printf("# density_matrices set for %d LiveAtoms\n", *na);
		std::fflush(stdout);
	} // _live_atom_set_density_matrix_
#else
	;
#endif

	fortran_callable(get_compensation_charge)
        ( int32_t const *na
		, double qlm[] // layout [na][(1 + maxval(lmax_qlm))^2]
		, int32_t *status)
#ifndef SINGLE_ATOM_HEADER_ONLY
	{
		*status = single_atom::atom_update("q", *na);;
		printf("# got_compensation_charge %d LiveAtoms\n", *na);
		std::fflush(stdout);
	} // _live_atom_get_compensation_charge_
#else
	;
#endif

	fortran_callable(set_potential_multipole)
        ( int32_t const *na
		, double const vlm[] // layout [na][(1 + maxval(lmax_vlm))^2]
		, int32_t *status)
#ifndef SINGLE_ATOM_HEADER_ONLY
	{
		*status = single_atom::atom_update("v", *na);;
		printf("# potential_multipoles set for %d LiveAtoms\n", *na);
		std::fflush(stdout);
	} // _live_atom_set_potential_multipole_
#else
	;
#endif

	fortran_callable(get_zero_potential)
        ( int32_t const *na
		, double vbar[] // layout [na][4096]
		, int32_t *status)
#ifndef SINGLE_ATOM_HEADER_ONLY
	{
		*status = single_atom::atom_update("z", *na);;
		printf("# got_zero_potential for %d LiveAtoms\n", *na);
		std::fflush(stdout);
	} // _live_atom_get_zero_potential_
#else
	;
#endif

	fortran_callable(get_hamiltonian_matrix)
        ( int32_t const *na
		, double hamiltonian_matrices[] // layout [na][stride][stride]
		, double overlap_matrices[]     // layout [na][stride][stride]
		, int32_t *status)
#ifndef SINGLE_ATOM_HEADER_ONLY
	{
		*status = single_atom::atom_update("h", *na);;
		printf("# got_hamiltonian_matrix for %d LiveAtoms\n", *na);
		std::fflush(stdout);
	} // _live_atom_get_hamiltonian_matrix_
#else
	;
#endif

	fortran_callable(get_energy_contributions)
        ( int32_t const *na
		, double energies[] // layout [na][...]
		, int32_t *status)
#ifndef SINGLE_ATOM_HEADER_ONLY
	{
		*status = single_atom::atom_update("e", *na);;
		printf("# got_energy_contributions for %d LiveAtoms\n", *na);
		std::fflush(stdout);
	} // _live_atom_get_energy_contributions_
#else
	;
#endif

	fortran_callable(finalize)
        ( int32_t const *na
		, int32_t *status)
#ifndef SINGLE_ATOM_HEADER_ONLY
	{
		*status = single_atom::atom_update("memory cleanup", *na);;
		printf("# finalized %d LiveAtoms\n", *na);
		std::fflush(stdout);
	} // _live_atom_finalize_
#else
	;
#endif

#endif // SINGLE_ATOM_HEADER (header guard)