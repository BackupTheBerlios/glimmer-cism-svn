! erosion_transport_interp.f90
! Magnus Hagdorn, May 2005
!
! this module transports some scalar quantity c through a 2D velo field using a 
! semi-lagrangian approach based on simple bilinear interpolation

module erosion_transport
  use erosion_types

  integer, parameter, private :: offset = 1

  real(kind=dp), private, parameter :: small = 1.d-7

contains
  subroutine init_transport(trans,model,erosion)
    use erosion_advect
    use glide_types
    use erosion_types
    implicit none
    type(er_transport_type) :: trans       ! structure holding transport stuff
    type(glide_global_type) :: model       ! model instance
    type(erosion_type) :: erosion          !*FD structure holding erosion data

    allocate(trans%lin_stuff(erosion%ewn*erosion%nsn))
    allocate(trans%lin_con(erosion%ewn*erosion%nsn))
    allocate(trans%lin_stuff2(erosion%ewn*erosion%nsn))

    trans%coord = coordsystem_new(0.d0,0.d0,erosion%dew,erosion%dns, erosion%ewn,erosion%nsn)
    allocate(trans%dispx(erosion%ewn,erosion%nsn))
    allocate(trans%dispy(erosion%ewn,erosion%nsn))

    call er_advect2d_init(model%general%velo_grid)

  end subroutine init_transport

  subroutine calc_lagrange(erosion, trans, deltat, lagrange)
    use glimmer_coordinates
    use erosion_advect
    use glimmer_sparse
    use glimmer_interpolate2d
    use erosion_types
    implicit none

    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(er_transport_type) :: trans       ! structure holding transport stuff
    real(kind=dp), intent(in) :: deltat       ! the time step
    type(sparse_matrix_type) :: lagrange  ! sparse matrix containing the weights

    ! local variables
    integer i,j
    type(coord_point) :: pt
    type(coord_ipoint) :: node
    real(kind=dp), dimension(1) ::  time,x,y

    ! calculate displacement field
    time =deltat
    do j=2,erosion%nsn-1
       node%pt(2) = j
       do i=2,erosion%ewn-1
          node%pt(1) = i
          pt = coordsystem_get_coord(trans%coord,node)
          call er_advect2d(time,pt%pt(1),pt%pt(2),x,y)
          trans%dispx(i,j)=x(1)
          trans%dispy(i,j)=y(1)
       end do
    end do

    ! setup sparse matrix
    call glimmer_init_bilinear(trans%coord,trans%dispx,trans%dispy,lagrange)

  end subroutine calc_lagrange

  subroutine transport_scalar(erosion,trans,concentration,lagrange)
    ! transport scalar concentration using sparse matrix lagrange
    use glimmer_sparse
    use glimmer_interpolate2d
    use erosion_types
    implicit none
    
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(er_transport_type) :: trans       ! structure holding transport stuff
    real(kind=dp), dimension(:,:) :: concentration
    type(sparse_matrix_type) :: lagrange

    ! local variables
    integer :: i,j,k
    type(coord_ipoint) :: node1,node2
    integer ierr

    call glimmer_interpolate(lagrange, concentration, concentration)


  end subroutine transport_scalar
    
end module erosion_transport
