program p
implicit none
  !! this is the construction of the SHO basis which can represent densities
  !! arising from wave functions expanded in a SHO basis with sigma=1, nmax
  !! then, the SHO basis for the density must have sigma=sqrt(1/2), 2*nmax
  
  !! with this, the cutoff energy is a factor 4 higher in the density basis.

  integer :: ix
  integer, parameter :: nmax = 5
  real :: x, H(0:nmax), HH(0:2*nmax)
  real :: s1, s2, s, E0

  s = sqrt(2.*nmax + 1.)
  E0 = (2*nmax+1)/s**2
  s1 = s
  s2 = s1*sqrt(2.)
  
  do ix = -399, 399
    x = ix/256.
    call Hermite(x*s1, nmax, H)
    call Hermite(x*s2, 2*nmax, HH)
    write(*,'(9(" ",f0.6))') x, x*x/E0, 1., 4*x*x/E0, 4., H(nmax), H(nmax)**2, HH(2*nmax)
  enddo ! ix

  !! both bases have the same classical return radius
  do ix = -1, 1, 2
    write(*,'(/,2(/,i0," ",f0.1))') ix, -1., ix, 9.
  enddo ! ix

  !! show also the matching points
  do ix = -1, 1, 2
    write(*,'(/,i0," ",i0)') ix, 1, ix, 4
  enddo ! ix

  contains
  
  subroutine Hermite(x, nmax, H)
    real, intent(in) :: x
    integer, intent(in) :: nmax
    real, intent(out) :: H(0:nmax)
    
    real :: Hnp1, Hn, Hnm1
    integer :: n
    
    Hn = exp(-0.5*x*x) ! H_0 is the Gaussian envelope function
    Hnm1 = 0
    do n = 0, nmax
      H(n) = Hn
      Hnp1 = x * Hn - (0.5*n) * Hnm1 !! best choice of scaling for the non-normalized projectors, no unnessary powers of 2
      !! rotate registers for the next iteration
      Hnm1 = Hn
      Hn = Hnp1
    enddo ! n
  
  endsubroutine
  
endprogram
