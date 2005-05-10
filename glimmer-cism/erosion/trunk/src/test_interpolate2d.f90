! test_interpolate2d.f90
! Magnus Hagdorn, May 2005
!
! testing 2d interpolation stuff

program test_interpolate
  use glimmer_coordinates
  use sparse
  use glimmer_interpolate2d
  implicit none

  type(coordsystem) :: coords
  type(sparse_matrix) :: interpolator
  real(kind=dp) :: delta = 0.1
  integer i,j
  integer, parameter :: numx = 100
  integer, parameter :: numy = 150
  integer, parameter :: inter_numx = 100
  integer, parameter :: inter_numy = 50

  real(kind=dp), dimension(numx,numy) ::  orig_field

  real(kind=dp), dimension(inter_numx,inter_numy) :: dispx, dispy, interp_field

  ! setup coordsystem
  coords = coordsystem_new(0.d0,0.d0,delta, delta, numx, numy)
  ! setup data
  do j=1,numy
     do i=1,numx
        orig_field(i,j) = calc_data((i-1)*delta,(j-1)*delta)
     end do
  end do

  ! generate random displacement field
  call random_number(dispx)
  dispx = dispx * numx*delta
  call random_number(dispy)
  dispy = dispy * numy*delta
  
  call glimmer_init_bilinear(coords,dispx, dispy,interpolator)
  call glimmer_interpolate(interpolator,orig_field,interp_field)
  
  do j=1,inter_numy
     do i=1,inter_numx
        write(*,*) dispx(i,j), dispy(i,j), calc_data(dispx(i,j), dispy(i,j)), interp_field(i,j)
     end do
  end do

contains 
  function calc_data(x,y)
    implicit none
    real(kind=dp) :: calc_data
    real(kind=dp), intent(in) :: x,y

    calc_data = sin(x)+sin(y)
  end function calc_data
end program test_interpolate

