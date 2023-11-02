program test_live_atom
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Fortran-interface for the LiveAtom library
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  All functionality is imported by subroutines
!!  call live_atom_...(..., status)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  ! compile time constants
  integer(kind=4), parameter :: na = 1 ! number of atoms
  integer, parameter :: nr2 = 2**12
  real(kind=8), parameter :: ar2 = 16.d0
  logical(kind=1), parameter :: plot = .true.

  ! variables
  integer(kind=4) :: status = 0, stride = 0, ia, ir2
  real(kind=8) :: Z_core(na), sigma(na) = 0, rcut(na) = 0
  real(kind=8) :: ionic(na) = 0, magn(na) = 0, nve(na) = 0
  real(kind=8) :: energy(na) = 0
  real(kind=8) :: waves(nr2,6,na) = 0
  real(kind=8) :: occupations(2,6,na) = 0
  integer(kind=4) :: atom_id(na) = -1, numax(na) = 0
  integer(kind=4) :: lmax_qlm(na) = 0, lmax_vlm(na) = 0
  integer(kind=1) :: nn(8,na) = 0
  character(len=32) :: xc_key = "LDA"//achar(0)
  real(kind=4) :: fp(na) = 0

  type doublePtr
      real(kind=8), pointer :: p(:)
  endtype doublePtr

  type(doublePtr) :: pointers(na) ! this is double** in C
  real(kind=8), target :: memory(nr2,na) = 0
  integer(kind=4) :: nonzero = 0, differ = 0

  ! executable section

  ! write(*,'(3a,i0)') "# ",__FILE__,":",__LINE__ ! here

  call live_atom_init_env("control.sh"//achar(0), status)
  write(*,'(a,i0)') "# live_atom_init_env = ", status

  call live_atom_set_env("called.from"//achar(0), &
                             "Fortran"//achar(0), status)
  write(*,'(a,i0)') "# live_atom_set_env = ", status

  do ia = 1, na
    Z_core(ia) = 32 + ia
    atom_id(ia) = ia - 1
    pointers(ia)%p => memory(:,ia)
  enddo ! ia

  call live_atom_initialize(na, Z_core, atom_id, numax, sigma, &
                      rcut, nn, ionic, magn, xc_key, stride, &
                      lmax_qlm, lmax_vlm, nve, status)
  write(*,'(a,i0)') "# live_atom_initialize = ", status

  call live_atom_get_core_density(na, pointers, status)
  write(*,'(a,i0)') "# live_atom_get_core_density = ", status
  if (plot) then
    do ia = 1, na
      write(*,'(a,i0)') "### r, core_density(r) for atom #", ia
      nonzero = 0
      do ir2 = nr2, 1, -1
        if (memory(ir2,ia) /= 0 .or. nonzero > 0) then
            write(*,'(f0.6," ",f0.9)') sqrt((ir2 - 1.d0)/ar2), &
                memory(ir2,ia)
            nonzero = nonzero + 1
        endif
      enddo ! ir2
      differ = count(pointers(ia)%p(:) /= memory(:,ia))
      if (differ > 0) then
        write(*,'(a,i0)') "### pointers(ia)%p(:) /= memory(:,ia) in ",&
                            differ," positions for atom #", ia
      endif ! differ
      write(*,'(a)') "" ! empty line
    enddo ! ia
  endif ! plot

  call live_atom_get_start_waves(na, waves, occupations, status)
  write(*,'(a,i0)') "# live_atom_get_start_waves = ", status

  call live_atom_get_compensation_charge(na, pointers, status)
  write(*,'(a,i0)') "# live_atom_get_compensation_charge = ", status

  call live_atom_get_hamiltonian_matrix(na, pointers, status)
  write(*,'(a,i0)') "# live_atom_get_hamiltonian_matrix = ", status

  call live_atom_get_zero_potential(na, pointers, status)
  write(*,'(a,i0)') "# live_atom_get_zero_potential = ", status

  memory = 0
  call live_atom_set_density_matrix(na, pointers, status)
  write(*,'(a,i0)') "# live_atom_set_density_matrix = ", status

  call live_atom_set_potential_multipole(na, pointers, status)
  write(*,'(a,i0)') "# live_atom_set_potential_multipole = ", status

  call live_atom_get_energy_contributions(na, energy, status)
  write(*,'(a,i0)') "# live_atom_get_energy_contributions = ", status

  call live_atom_update("direct"//achar(0), na, ionic, numax, &
                        fp, pointers, status)
  write(*,'(a,i0)') "# live_atom_update = ", status

  call live_atom_finalize(na, status)
  write(*,'(a,i0)') "# live_atom_finalize = ", status

! write(*,'(3a,i0)') "# ",__FILE__,":",__LINE__ ! here
endprogram ! test_live_atom
