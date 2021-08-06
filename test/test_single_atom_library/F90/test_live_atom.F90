program test_live_atom
implicit none
  ! compile time constants
  integer(kind=4), parameter :: na = 2 ! number of atoms
! integer(kind=4), parameter :: nr = 512 ! number of radial grid points
! real(kind=8), parameter :: Rmax = 12.6 ! maximum outer radius
! real(kind=8), parameter :: dr = Rmax/nr ! equidistant radial grid spacing

  ! variables
  integer(kind=4) :: status = 0, stride = 0, ia
  real(kind=8) :: Z_core(na), sigma(na) = 0, rcut(na) = 0
  real(kind=8) :: ionic(na) = 0, magn(na) = 0, nve(na) = 0
  integer(kind=4) :: atom_id(na) = -1, numax(na) = 0
  integer(kind=4) :: lmax_qlm(na) = 0, lmax_vlm(na) = 0
  integer(kind=1) :: nn(8,na) = 0
  character(len=32) :: xc_key = "LDA"//achar(0)

  ! executable section
  write(*,'(2a,i0)') __FILE__,":",__LINE__ ! here

  call live_atom_init_env("control.sh"//achar(0), status)
  write(*,'(a,i0)') "live_atom_init_env = ", status

  do ia = 1, na
    Z_core(ia) = ia
  enddo ! ia

  call live_atom_initialize(na, Z_core, atom_id, numax, sigma, &
                      rcut, nn, ionic, magn, xc_key, stride, &
                      lmax_qlm, lmax_vlm, nve, status)
  write(*,'(a,i0)') "live_atom_initialize = ", status

!   call live_atom_get_compensation_charge
!   call live_atom_get_core_density
!   call live_atom_get_energy_contributions
!   call live_atom_get_hamiltonian_matrix
!   call live_atom_get_start_waves
!   call live_atom_get_zero_potential
!   call live_atom_set_density_matrix
!   call live_atom_set_env
!   call live_atom_set_potential_multipole

  call live_atom_finalize(na, status)
  write(*,'(a,i0)') "live_atom_finalize = ", status

  write(*,'(2a,i0)') __FILE__,":",__LINE__ ! here
endprogram ! test_live_atom
