! test_integrate2d.f90
! Magnus Hagdorn, March 2003
!
! testing the 2d integration stuff

program testintegrate2d
  use erosion_integrate2d
  use sparse
  use geometry
  use glimmer_coordinates
  implicit none
  type(geom_poly) :: patch
  real, parameter :: xsize = 200.
  real, parameter :: ysize = 200.
  real :: deltaxy
  integer :: numi,numj
  integer i,j
  type(sparse_matrix) :: weight
  real, dimension(:,:), allocatable :: vec
  real, dimension(:), allocatable ::res
  real :: area
  type(coordsystem) :: coords

  deltaxy = 10.

  numi = int(xsize/deltaxy)+1
  numj = int(ysize/deltaxy)+1

  coords = coordsystem_new(0.,0.,deltaxy,deltaxy,numi,numj)

  allocate(res(numi*numj))
  allocate(vec(numi,numj))

  ! initialise grid
  call init_integrate2d

  ! test shape
  patch = poly_new(4)
  call poly_add_vert(patch,60.,15.)
  call poly_add_vert(patch,95.,46.)
  call poly_add_vert(patch,55.,80.)
  call poly_add_vert(patch,10.,42.)

  write(*,*) 'Size :', 0.5*real(poly_area2(patch))


  weight = new_sparse_matrix(1000)
  call calc_weight(coords,weight, patch, 1)

  vec = 0
  do i=1,weight%n
     vec(mod(weight%pos(2,i),numi)+1,weight%pos(2,i)/numi+1) = weight%val(i)
  end do

  !do j=1,numj
  !   do i=1,numi
  !      !write(*,fmt='(F5.2)',advance='no') vec(i,j)
  !      write(*,*) vec(i,j)
  !   end do
  !   write(*,*) ''
  !end do

  do i=0,8
     deltaxy = 2.**i
     numi = int(xsize/deltaxy) + 1
     numj = int(ysize/deltaxy) + 1
     coords = coordsystem_new_real(0.,0.,deltaxy,deltaxy,numi,numj)
     area = 0.
     weight%n = 0
     call calc_weight(coords,weight, patch, 1)
     area = sum(weight%val(1:weight%n))*deltaxy*deltaxy
     write(*,*) deltaxy,weight%n,area
  end do

end program testintegrate2d
