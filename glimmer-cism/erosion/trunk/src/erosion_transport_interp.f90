! erosion_transport_interp.f90
! Magnus Hagdorn, May 2005
!
! this module transports some scalar quantity c through a 2D velo field using a 
! semi-lagrangian approach based on simple bilinear interpolation

module erosion_transport
  use geometry
  use glimmer_coordinates

  integer, parameter, private :: offset = 1

  type er_transport_type
     ! private data

     real(kind=dp), dimension(:), pointer :: lin_stuff,lin_stuff2,lin_con
     type(coordsystem) :: coord
     real(kind=dp), dimension(:,:), pointer :: dispx => NULL() ! x-displacement field
     real(kind=dp), dimension(:,:), pointer :: dispy => NULL() ! y-displacement field

  end type er_transport_type

  real(kind=dp), private, parameter :: small = 1.d-7

contains
  subroutine init_transport(trans,model)
    use erosion_advect
    use glide_types
    implicit none
    type(er_transport_type) :: trans       ! structure holding transport stuff
    type(glide_global_type) :: model       ! model instance

    allocate(trans%lin_stuff(model%general%ewn*model%general%nsn))
    allocate(trans%lin_con(model%general%ewn*model%general%nsn))
    allocate(trans%lin_stuff2(model%general%ewn*model%general%nsn))

    trans%coord = coordsystem_new(0.d0,0.d0,model%numerics%dew,model%numerics%dns, &
         model%general%ewn,model%general%nsn)
    allocate(trans%dispx(model%general%ewn,model%general%nsn))
    allocate(trans%dispy(model%general%ewn,model%general%nsn))

    call er_advect2d_init(0.5*model%numerics%dew, 0.5*model%numerics%dns, &
         model%general%ewn-1,model%general%nsn-1, &
         model%numerics%dew,model%numerics%dns)

  end subroutine init_transport

  subroutine calc_lagrange(model, trans, deltat, lagrange)
    use glimmer_coordinates
    use erosion_advect
    use sparse
    use glimmer_interpolate2d
    use glide_types
    implicit none

    type(glide_global_type) :: model       ! model instance
    type(er_transport_type) :: trans       ! structure holding transport stuff
    real(kind=dp), intent(in) :: deltat       ! the time step
    type(sparse_matrix) :: lagrange  ! sparse matrix containing the weights

    ! local variables
    integer i,j
    type(geom_point) :: pt
    type(geom_ipoint) :: node
    real(kind=dp), dimension(1) ::  time,x,y

    ! calculate displacement field
    time =deltat
    do j=2,model%general%nsn-1
       node%pt(2) = j
       do i=2,model%general%ewn-1
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

  subroutine transport_scalar(model,trans,concentration,lagrange)
    ! transport scalar concentration using sparse matrix lagrange
    use sparse
    use glimmer_interpolate2d
    use glide_types
    implicit none
    
    type(glide_global_type) :: model       ! model instance
    type(er_transport_type) :: trans       ! structure holding transport stuff
    real(kind=dp), dimension(:,:) :: concentration
    type(sparse_matrix) :: lagrange

    ! local variables
    integer :: i,j,k
    type(geom_ipoint) :: node1,node2
    integer ierr

    call glimmer_interpolate(lagrange, concentration, concentration)


  end subroutine transport_scalar
    
end module erosion_transport
