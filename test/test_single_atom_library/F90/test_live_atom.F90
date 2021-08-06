program test_live_atom
implicit none
  ! compile time constants
  integer, parameter :: na = 1 ! number of atoms
  integer, parameter :: nr = 512 ! number of radial grid points
  real(8), parameter :: Rmax = 12.6 ! maximum outer radius
  real(8), parameter :: dr = Rmax/nr ! equidistant radial grid spacing

  ! external subroutines
!   external :: live_atom_init_env

  ! variables
  integer :: status
  
  ! executable section
  write(*,'(2a,i0)') __FILE__,":",__LINE__ ! here

  call live_atom_init_env("control.sh", status)
! __live_atom_get_compensation_charge_
! __live_atom_get_core_density_
! __live_atom_get_energy_contributions_
! __live_atom_get_hamiltonian_matrix_
! __live_atom_get_start_waves_
! __live_atom_get_zero_potential_
! __live_atom_initialize_
! __live_atom_set_density_matrix_
! __live_atom_set_env_
! __live_atom_set_potential_multipole_
! __live_atom_finalize_

  write(*,'(2a,i0)') __FILE__,":",__LINE__ ! here
endprogram ! test_live_atom
