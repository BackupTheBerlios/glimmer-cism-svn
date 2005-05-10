! glimmer_interpolate2d.f90
! Magnus Hagdorn, May 2005
!
! module for 2D interpolation

module glimmer_interpolate2d

  use sparse

contains

  subroutine glimmer_interpolate(interpolator, infield, outfield)
    use glimmer_coordinates
    implicit none
    !*FD interpolate data on infield (defined on inccord system) onto outfield using interpolator
    type(sparse_matrix), intent(in) :: interpolator       !*FD sparse matrix containing interpolator
    real(kind=dp), dimension(:,:), intent(in)  :: infield !*FD input data
    real(kind=dp), dimension(:,:), intent(out) :: outfield!*FD output field

    ! local variables
    real(kind=dp), dimension(:), allocatable :: linear_out
    integer i,j, k, numx,numy

    numx = size(outfield,1)
    numy = size(outfield,2)
    
    allocate(linear_out(numx*numy))

    call sparse_matrix_vec_prod(interpolator,pack(infield,.true.), linear_out)

    do j=1,numy
       k = (j-1)*numx
       do i=1,numx
          outfield(i,j) = linear_out(k+i)
       end do
    end do

    deallocate(linear_out)
  end subroutine glimmer_interpolate

  subroutine glimmer_init_bilinear(coord, dispx, dispy, bilinear)
    use glimmer_coordinates
    implicit none
    !*FD initialise sparse matrix defining interpolation given by displacement field
    type(coordsystem), intent(in)             :: coord        !*FD coordinate system of the input field
    real(kind=dp), dimension(:,:), intent(in) :: dispx, dispy !*FD displacement field 
    type(sparse_matrix)                       :: bilinear     !*FD sparse matrix containing interpolation

    ! local variables
    integer :: i,j, k, lini, numx,numy
    type(geom_point) :: point
    type(geom_ipoint) :: this_node
    type(geom_ipoint), dimension(4) :: nodes
    real(kind=dp), dimension(4) :: weights


    numx = size(dispx,1)
    numy = size(dispx,2)

    call new_sparse_matrix(numx*numy*4,bilinear)

    do j=1,numy
       this_node%pt(2)=j
       do i=1,numx
          this_node%pt(1)=i
          lini = coordsystem_linearise2d(coord,this_node)
          point%pt(1)=dispx(i,j)
          point%pt(2)=dispy(i,j)
          call glimmer_bilinear(coord, point, nodes, weights)
          do k=1,4
             call add_to_matrix(bilinear,lini,coordsystem_linearise2d(coord,nodes(k)),weights(k))
          end do
       end do
    end do

  end subroutine glimmer_init_bilinear


  subroutine glimmer_bilinear(coord, point, nodes, weights)
    use glimmer_coordinates
    implicit none
    !*FD bilinear interpolation
    type(coordsystem), intent(in)            :: coord     !*FD coordinate system to operate on
    type(geom_point), intent(in)             :: point     !*FD desired point
    type(geom_ipoint), dimension(4), intent(out) :: nodes !*FD array containing indicies into field
    real(kind=dp), dimension(4), intent(out) :: weights   !*FD array of weights

    type(geom_point) :: pnt

    ! check if point is inside coord-system
#ifdef DEBUG
    if (.not.coordsystem_point_inside(coord,point)) then
       write(*,*) 'Error (',__FILE__,__LINE__,'): point (',point%pt,') not inside coord system'
       return
    end if
#endif

    nodes(1)%pt(:) = 1+floor((point%pt(:)-coord%origin%pt(:))*coord%delta_r%pt(:))

    nodes(2)%pt(1) = nodes(1)%pt(1) + 1
    nodes(2)%pt(2) = nodes(1)%pt(2)

    nodes(3)%pt(1) = nodes(1)%pt(1) + 1
    nodes(3)%pt(2) = nodes(1)%pt(2) + 1

    nodes(4)%pt(1) = nodes(1)%pt(1)
    nodes(4)%pt(2) = nodes(1)%pt(2) + 1

    pnt%pt(:) = (point%pt(:)-coord%origin%pt(:)-(nodes(1)%pt(:)-1)*coord%delta%pt(:))*coord%delta_r%pt(:)

    weights(1) = (1-pnt%pt(1))*(1-pnt%pt(2))
    weights(2) = pnt%pt(1)*(1-pnt%pt(2))
    weights(3) = pnt%pt(1)*pnt%pt(2)
    weights(4) = (1-pnt%pt(1))*pnt%pt(2)
  end subroutine glimmer_bilinear


end module glimmer_interpolate2d
