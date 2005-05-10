! erosion_transport_fv.f90
! Magnus Hagdorn, March 2003
!
! this module transports some scalar quantity c through a 2D velo field using a 
! semi-lagrangian approach based on a finite volume approach

module erosion_transport
  use geometry
  use glimmer_coordinates

  integer, parameter, private :: offset = 2

  type er_transport_type
     ! private data
     real(kind=dp), dimension(:), pointer :: lin_stuff,lin_stuff2,lin_con
     type(geom_point), dimension(:,:), pointer :: patch_strip 
     type(geom_poly) :: patch, patch1, patch2
     type(coordsystem) :: coord
     real(kind=dp) :: half_xstep, half_ystep
  end type er_transport_type

  real(kind=dp), private, parameter :: small = 1.d-7

contains
  subroutine init_transport(trans,model)
    use erosion_advect
    use erosion_integrate2d
    use glide_types
    implicit none
    type(er_transport_type) :: trans       ! structure holding transport stuff
    type(glide_global_type) :: model       ! model instance

    allocate(trans%lin_stuff(model%general%ewn*model%general%nsn))
    allocate(trans%lin_con(model%general%ewn*model%general%nsn))
    allocate(trans%lin_stuff2(model%general%ewn*model%general%nsn))
    
    trans%patch = poly_new(4)
    trans%patch1 = poly_new(3)
    trans%patch2 = poly_new(3)
    trans%patch%n = 4
    trans%patch1%n = 3
    trans%patch2%n = 3
    allocate(trans%patch_strip(model%general%ewn,2))
    trans%coord = coordsystem_new(0.d0,0.d0,model%numerics%dew,model%numerics%dns, &
         model%general%ewn,model%general%nsn)
    call init_integrate2d
    trans%half_xstep = 0.5*model%numerics%dew
    trans%half_ystep = 0.5*model%numerics%dns
    call er_advect2d_init(trans%half_xstep,trans%half_ystep, &
         model%general%ewn-1,model%general%nsn-1, &
         model%numerics%dew,model%numerics%dns)
    ismintegrate2d_zero = 1e-5 ! effective zero
  end subroutine init_transport

  subroutine calc_lagrange(model, trans, deltat, lagrange)
    use glimmer_coordinates
    use erosion_advect
    use sparse
    use erosion_integrate2d
    use glide_types
    implicit none

    type(glide_global_type) :: model       ! model instance
    type(er_transport_type) :: trans       ! structure holding transport stuff
    real(kind=dp), intent(in) :: deltat       ! the time step
    type(sparse_matrix) :: lagrange  ! sparse matrix containing the weights

    ! local variables
    integer i,j
    integer lower,upper
    real(kind=dp) x0,y0
    real(kind=dp), dimension(1) ::  time,x,y
    logical d1,d2, l
    type(geom_point) :: cross
    type(geom_point) :: pt
    type(geom_ipoint) :: node

    time = deltat
    lagrange%n = 0
    ! setup first triangle strip
    lower = 1
    upper = 2
    ! fill in first element
    node%pt(:) = 1+offset
    pt = coordsystem_get_coord(trans%coord,node)
    x0 = pt%pt(1) - trans%half_xstep
    y0 = pt%pt(2) - trans%half_ystep
    call er_advect2d(time,x0,y0,x,y)
    trans%patch_strip(1,lower)%pt(1) = x(1)
    trans%patch_strip(1,lower)%pt(2) = y(1)
    ! fill the first row
    do i=1+offset,model%general%ewn-offset
       x0 = x0 + model%numerics%dew
       call er_advect2d(time,x0,y0,x,y)
       trans%patch_strip(i-offset+1,lower)%pt(1) = x(1)
       trans%patch_strip(i-offset+1,lower)%pt(2) = y(1)
    end do

    ! loop over rows
    do j=1+offset,model%general%nsn-offset
       ! fill in first element of upper strip
       node%pt(1) = 1+offset
       node%pt(2) = j
       pt = coordsystem_get_coord(trans%coord,node)
       x0 = pt%pt(1) - trans%half_xstep
       y0 = pt%pt(2) + trans%half_ystep
       call er_advect2d(time,x0,y0,x,y)
       trans%patch_strip(1,upper)%pt(1) = x(1)
       trans%patch_strip(1,upper)%pt(2) = y(1)
       ! fill upper strip
       do i=1+offset,model%general%ewn-offset
          x0 = x0 + model%numerics%dew
          call er_advect2d(time,x0,y0,x,y)
          trans%patch_strip(i-offset+1,upper)%pt(1) = x(1)
          trans%patch_strip(i-offset+1,upper)%pt(2) = y(1)
       end do
       !calculate weights 
       do i=1+offset,model%general%ewn-offset
          trans%patch%poly(1) = trans%patch_strip(i-offset,lower)
          trans%patch%poly(2)%pt(:) = trans%patch_strip(i-offset+1,lower)%pt(:)-small
          trans%patch%poly(3)%pt(:) = trans%patch_strip(i-offset+1,upper)%pt(:)-small
          trans%patch%poly(4) = trans%patch_strip(i-offset,upper)
          node%pt(1) = i
          node%pt(2) = j
          if (is_anti_clock(trans%patch)) then
             ! check diagonals
             d1 = diagonal(1,3,trans%patch)
             d2 = diagonal(2,4,trans%patch)
             ! checking that we have a nicely ordered polygon
             if (d1 .and. d2) then  ! convex patch
                call calc_weight(trans%coord, lagrange, trans%patch, coordsystem_linearise2d(trans%coord,node))
             else
                if (d1) then
                   trans%patch1%poly(1) = trans%patch%poly(1)
                   trans%patch1%poly(2) = trans%patch%poly(2)
                   trans%patch1%poly(3) = trans%patch%poly(3)
                   
                   trans%patch2%poly(1) = trans%patch%poly(1)
                   trans%patch2%poly(2) = trans%patch%poly(3)
                   trans%patch2%poly(3) = trans%patch%poly(4)
                else
                   trans%patch1%poly(1) = trans%patch%poly(2)
                   trans%patch1%poly(2) = trans%patch%poly(3)
                   trans%patch1%poly(3) = trans%patch%poly(4)

                   trans%patch2%poly(1) = trans%patch%poly(2)
                   trans%patch2%poly(2) = trans%patch%poly(4)
                   trans%patch2%poly(3) = trans%patch%poly(1)
                end if
                call calc_weight(trans%coord, lagrange, trans%patch1, coordsystem_linearise2d(trans%coord,node))
                call calc_weight(trans%coord, lagrange, trans%patch2, coordsystem_linearise2d(trans%coord,node))
             end if
          else
             if (intersection(trans%patch%poly(1),trans%patch%poly(2),trans%patch%poly(3),trans%patch%poly(4),cross)) then
                trans%patch1%poly(1) = trans%patch%poly(1)
                trans%patch1%poly(2) = cross                   
                trans%patch1%poly(3) = trans%patch%poly(4)
                trans%patch2%poly(1) = cross
                trans%patch2%poly(2) = trans%patch%poly(3)
                trans%patch2%poly(3) = trans%patch%poly(2)
             end if
             if (intersection(trans%patch%poly(1),trans%patch%poly(4),trans%patch%poly(2),trans%patch%poly(3),cross)) then
                trans%patch1%poly(1) = trans%patch%poly(2)
                trans%patch1%poly(2) = trans%patch%poly(1)
                trans%patch1%poly(3) = cross
                trans%patch2%poly(1) = cross
                trans%patch2%poly(2) = trans%patch%poly(4)
                trans%patch2%poly(3) = trans%patch%poly(3)
             end if
             call calc_weight(trans%coord, lagrange, trans%patch1, coordsystem_linearise2d(trans%coord,node))
             call calc_weight(trans%coord, lagrange, trans%patch2, coordsystem_linearise2d(trans%coord,node))
          end if
       end do
       ! swap rows
       lower = mod(lower,2)+1
       upper = mod(upper,2)+1
    end do
  end subroutine calc_lagrange

  subroutine transport_scalar(model,trans,concentration,lagrange)
    ! transport scalar concentration using sparse matrix lagrange
    use sparse
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

    ! linearise concentration
    trans%lin_stuff = pack(concentration,.true.)
    trans%lin_stuff2 = 0.
    trans%lin_con = 0.

    ! normalise sparse matrix to 1 iff val > 1
    do k=1,lagrange%n
       trans%lin_stuff2(lagrange%pos(2,k)) = trans%lin_stuff2(lagrange%pos(2,k))+lagrange%val(k)
    end do
    ! we should propably do an allreduce as well. but hell!
    do k=1,lagrange%n
       if (trans%lin_stuff2(lagrange%pos(2,k)).gt.1.) then
          lagrange%val(k) = lagrange%val(k)/trans%lin_stuff2(lagrange%pos(2,k))
       end if
    end do
    trans%lin_stuff2 = 0.

    ! calculate new concentrations
    do k=1,lagrange%n
       trans%lin_con(lagrange%pos(1,k)) = trans%lin_con(lagrange%pos(1,k)) + &
            trans%lin_stuff(lagrange%pos(2,k))*lagrange%val(k)
       trans%lin_stuff2(lagrange%pos(2,k)) = trans%lin_stuff2(lagrange%pos(2,k)) - &
            trans%lin_stuff(lagrange%pos(2,k))*lagrange%val(k)
    end do    

    trans%lin_stuff2 = trans%lin_stuff2 + trans%lin_stuff

    ! unpack data
    node1%pt(1)=1
    node2%pt(1)=model%general%ewn
    do j=1,model%general%nsn
       node1%pt(2)=j
       node2%pt(2)=j
       concentration(1:model%general%ewn,j) = &
            trans%lin_con(coordsystem_linearise2d(trans%coord,node1):coordsystem_linearise2d(trans%coord,node2)) &
            + trans%lin_stuff2(coordsystem_linearise2d(trans%coord,node1):coordsystem_linearise2d(trans%coord,node2))
    end do

  end subroutine transport_scalar
    
end module erosion_transport
