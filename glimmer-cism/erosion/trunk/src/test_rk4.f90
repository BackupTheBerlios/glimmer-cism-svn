! test_rk4.F90
! Magnus Hagdorn, March 2003
!
! testing rk4 code

program testrk4
  use rk4module
  implicit none

  real, dimension(2) :: startx,x
  real :: h,t0,t1,deltat,t
  integer i, nok, nbad

  external test1

  startx(1) = 2
  startx(2) = 1
  t0 = 0
  t1 = 1.1
  deltat = 0.1

  x = startx
  h = 0.01
  do i=0,int((t1-t0)/deltat)-1
     t = t0+i*deltat
     call odeint(x,t,t+deltat, 0.000001, h, 0., nok,nbad, test1)
     write(*,*) t, x(1),x(2), nok,nbad
  end do

end program testrk4

subroutine test1(t,x,dxdt)
  implicit none
  real, intent(in) :: t
  real, intent(in), dimension(2) :: x
  real, intent(out), dimension(2) :: dxdt

  dxdt(1) = x(2)
  dxdt(2) = x(1)+x(2)
end subroutine test1
