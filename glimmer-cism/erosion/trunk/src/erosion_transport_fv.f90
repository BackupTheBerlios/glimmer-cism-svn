! erosion_transport_fv.f90
! Magnus Hagdorn, March 2003
!
! this module transports some scalar quantity c through a 2D velo field using a 
! semi-lagrangian approach based on a finite volume approach

module erosion_transport
  use geometry

  integer, parameter, private :: offset = 2

  real(kind=dp), private, parameter :: small = 1.d-7

contains
  subroutine init_transport(trans,model,erosion)
    use glide_types
    use erosion_advect
    use erosion_integrate2d
    use erosion_types
    implicit none
    type(er_transport_type) :: trans       ! structure holding transport stuff
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model

    allocate(trans%lin_stuff(erosion%ewn*erosion%nsn))
    allocate(trans%lin_con(erosion%ewn*erosion%nsn))
    allocate(trans%lin_stuff2(erosion%ewn*erosion%nsn))
    
    trans%patch = poly_new(4)
    trans%patch1 = poly_new(3)
    trans%patch2 = poly_new(3)
    trans%patch%n = 4
    trans%patch1%n = 3
    trans%patch2%n = 3
    allocate(trans%patch_strip(erosion%ewn,2))
    call init_integrate2d
    trans%half_xstep = 0.5*erosion%dew
    trans%half_ystep = 0.5*erosion%dns
    call er_advect2d_init(model%general%velo_grid)
    ismintegrate2d_zero = 0. ! effective zero
  end subroutine init_transport

  subroutine calc_lagrange(erosion, trans, deltat, lagrange)
    use glimmer_coordinates
    use erosion_advect
    use glimmer_sparse
    use erosion_integrate2d
    use erosion_types
    implicit none

    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(er_transport_type) :: trans       ! structure holding transport stuff
    real(kind=dp), intent(in) :: deltat       ! the time step
    type(sparse_matrix_type) :: lagrange  ! sparse matrix containing the weights

    ! local variables
    integer i,j
    integer lower,upper
    real(kind=dp) x0,y0
    real(kind=dp), dimension(1) ::  time,x,y
    logical d1,d2, l
    type(coord_point) :: cross
    type(coord_point) :: pt
    type(coord_ipoint) :: node

    time = deltat
    lagrange%n = 0
    ! setup first triangle strip
    lower = 1
    upper = 2
    ! fill in first element
    node%pt(:) = 1+offset
    pt = coordsystem_get_coord(erosion%coord,node)
    x0 = pt%pt(1) - trans%half_xstep
    y0 = pt%pt(2) - trans%half_ystep
    call er_advect2d(time,x0,y0,x,y)
    trans%patch_strip(1,lower)%pt(1) = x(1)
    trans%patch_strip(1,lower)%pt(2) = y(1)
    ! fill the first row
    do i=1+offset,erosion%ewn-offset
       x0 = x0 + erosion%dew
       call er_advect2d(time,x0,y0,x,y)
       trans%patch_strip(i-offset+1,lower)%pt(1) = x(1)
       trans%patch_strip(i-offset+1,lower)%pt(2) = y(1)
    end do

    ! loop over rows
    do j=1+offset,erosion%nsn-offset
       ! fill in first element of upper strip
       node%pt(1) = 1+offset
       node%pt(2) = j
       pt = coordsystem_get_coord(erosion%coord,node)
       x0 = pt%pt(1) - trans%half_xstep
       y0 = pt%pt(2) + trans%half_ystep
       call er_advect2d(time,x0,y0,x,y)
       trans%patch_strip(1,upper)%pt(1) = x(1)
       trans%patch_strip(1,upper)%pt(2) = y(1)
       ! fill upper strip
       do i=1+offset,erosion%ewn-offset
          x0 = x0 + erosion%dew
          call er_advect2d(time,x0,y0,x,y)
          trans%patch_strip(i-offset+1,upper)%pt(1) = x(1)
          trans%patch_strip(i-offset+1,upper)%pt(2) = y(1)
       end do
       !calculate weights 
       do i=1+offset,erosion%ewn-offset
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
                call calc_weight(erosion%coord, lagrange, trans%patch, coordsystem_linearise2d(erosion%coord,node))
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
                call calc_weight(erosion%coord, lagrange, trans%patch1, coordsystem_linearise2d(erosion%coord,node))
                call calc_weight(erosion%coord, lagrange, trans%patch2, coordsystem_linearise2d(erosion%coord,node))
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
             call calc_weight(erosion%coord, lagrange, trans%patch1, coordsystem_linearise2d(erosion%coord,node))
             call calc_weight(erosion%coord, lagrange, trans%patch2, coordsystem_linearise2d(erosion%coord,node))
          end if
       end do
       ! swap rows
       lower = mod(lower,2)+1
       upper = mod(upper,2)+1
    end do
  end subroutine calc_lagrange

  subroutine transport_scalar(erosion,trans,concentration,lagrange)
    ! transport scalar concentration using sparse matrix lagrange
    use glimmer_sparse
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

    ! linearise concentration
    trans%lin_stuff = pack(concentration,.true.)
    trans%lin_stuff2 = 0.
    trans%lin_con = 0.

    ! calculate new concentrations
    do k=1,lagrange%n
       trans%lin_con(lagrange%col(k)) = trans%lin_con(lagrange%col(k)) + &
            trans%lin_stuff(lagrange%row(k))*lagrange%val(k)
       trans%lin_stuff2(lagrange%row(k)) = trans%lin_stuff2(lagrange%row(k)) - &
            trans%lin_stuff(lagrange%row(k))*lagrange%val(k)
    end do    

    trans%lin_stuff2 = trans%lin_stuff2 + trans%lin_stuff

    ! unpack data
    node1%pt(1)=1
    node2%pt(1)=erosion%ewn
    do j=1,erosion%nsn
       node1%pt(2)=j
       node2%pt(2)=j
       concentration(1:erosion%ewn,j) = &
            trans%lin_con(coordsystem_linearise2d(erosion%coord,node1):coordsystem_linearise2d(erosion%coord,node2)) &
            + trans%lin_stuff2(coordsystem_linearise2d(erosion%coord,node1):coordsystem_linearise2d(erosion%coord,node2))
    end do

  end subroutine transport_scalar
    
end module erosion_transport
