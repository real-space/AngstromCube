implicit none
    integer :: ix, i
    real :: spd(0:5), ovl(0:5,0:5), x, xm, xp, scal(0:5)
    real, parameter :: dx = .0125
    ovl = 0
    do ix = -999, 999
      x = ix*dx
      xm = x - 1.582 ; spd(0:2) = exp(-.5*xm**2)*[1.,  xm, xm*xm-.5]
      xp = x + 1.582 ; spd(3:5) = exp(-.5*xp**2)*[1., -xp, xp*xp-.5]
      write(*,'(9(" ",f0.6))') x, spd ! warning: not normalized
      do i = 0, 5
        ovl(:,i) = ovl(:,i) + spd(:) * spd(i)
      enddo ! i
    enddo
    ovl = ovl * dx
    do i = 0, 5
      scal(i) = ovl(i,i)**-0.5
    enddo ! i
    do i = 0, 5
      ovl(:,i) = scal(:) * ovl(:,i) * scal(i)
    enddo ! i
!     write(*,'(6(" ",f7.3))') ovl
    
end
