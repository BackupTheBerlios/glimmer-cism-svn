! glimmer_coordinates.f90
! Magnus Hagdorn, April 2003
!
! module for handling regular coordinate systems

module glimmer_coordinates
  use geometry

  ! define unit for error_logs
  integer :: coordinates_unit = 6

  ! new type describing coordinate systems
  type coordsystem
     type(geom_point) :: origin   ! origin of coordinate space
     type(geom_point) :: delta    ! stepsize in x and y direction
     type(geom_point) :: delta_r  ! reciprocal stepsize in x and y direction
     type(geom_ipoint) :: size    ! extent in x and y direction
  end type coordsystem
  
  ! interface of creating new coord system
  interface coordsystem_new
     module procedure coordsystem_new_real, coordsystem_new_pt
  end interface

contains
  subroutine coordsystem_print(coord, unit)
    ! print coordsystem info to unit
    implicit none
    type(coordsystem), intent(in) :: coord
    integer unit
    write(unit,*) 'Origin  ',coord%origin%pt
    write(unit,*) 'Delta   ',coord%delta%pt
    write(unit,*) '1/Delta ',coord%delta_r%pt
    write(unit,*) 'Size    ',coord%size%pt
  end subroutine coordsystem_print

  function coordsystem_new_real(ox, oy, dx, dy, sx, sy)
    ! create new coordinate system from individual variables
    implicit none
    real(kind=dp), intent(in) :: ox, oy, dx, dy
    integer, intent(in) :: sx, sy
    type(coordsystem) :: coordsystem_new_real
    
    ! origin
    coordsystem_new_real%origin%pt(1) = ox
    coordsystem_new_real%origin%pt(2) = oy
    ! deltas
    coordsystem_new_real%delta%pt(1) = dx
    coordsystem_new_real%delta%pt(2) = dy
    coordsystem_new_real%delta_r%pt(1) = 1.d0/dx
    coordsystem_new_real%delta_r%pt(2) = 1.d0/dy
    ! size
    coordsystem_new_real%size%pt(1) = sx
    coordsystem_new_real%size%pt(2) = sy
  end function coordsystem_new_real

  function coordsystem_new_pt(o, d, s)
    ! create new coordinate system from points
    implicit none
    type(geom_point), intent(in) :: o, d
    type(geom_ipoint), intent(in) :: s
    type(coordsystem) :: coordsystem_new_pt

    ! origin
    coordsystem_new_pt%origin = o
    ! deltas
    coordsystem_new_pt%delta = d
    coordsystem_new_pt%delta_r%pt(:) = 1.d0/d%pt(:)
    ! size
    coordsystem_new_pt%size = s
  end function coordsystem_new_pt

  function coordsystem_get_coord(coord,node)
    ! get coordinates of node
    implicit none
    type(coordsystem), intent(in) :: coord
    type(geom_ipoint), intent(in) :: node

    type(geom_point) :: coordsystem_get_coord
  
#ifdef DEBUG_COORDS
    if (.not.coordsystem_node_inside(coord,node)) then
       write(coordinates_unit,*) 'Error (',__FILE__,__LINE__,'): node (',node%pt,') not inside coord system'
       call coordsystem_print(coord,coordinates_unit)
       stop
    end if
#endif

    coordsystem_get_coord%pt(:) = coord%origin%pt(:) + (node%pt(:) - 1)*coord%delta%pt(:)
  end function coordsystem_get_coord

  function coordsystem_get_node(coord,point)
    ! get index of nearest node given coords of a point
    implicit none
    type(coordsystem), intent(in) :: coord
    type(geom_point), intent(in) :: point
    
    type(geom_ipoint) :: coordsystem_get_node
    
    coordsystem_get_node%pt(:) = 1+floor(0.5+(point%pt(:)-coord%origin%pt(:))/coord%delta%pt(:))
    if (coordsystem_get_node%pt(1).eq.coord%size%pt(1)+1) coordsystem_get_node%pt(1) = coord%size%pt(1)
    if (coordsystem_get_node%pt(2).eq.coord%size%pt(2)+1) coordsystem_get_node%pt(2) = coord%size%pt(2)

#ifdef DEBUG_COORDS
    if (.not.coordsystem_node_inside(coord,coordsystem_get_node)) then
       write(coordinates_unit,*) 'Error (',__FILE__,__LINE__,'): point (',point%pt,') not inside coord system'
       call coordsystem_print(coord,coordinates_unit)
       stop
    end if
#endif
  end function coordsystem_get_node

  function coordsystem_node_inside(coord,node)
    ! return true iff node is inside coord system
    implicit none
    type(coordsystem), intent(in) :: coord
    type(geom_ipoint), intent(in) :: node

    logical coordsystem_node_inside
    
    coordsystem_node_inside = (all(node%pt.ge.1) .and. all(node%pt.le.coord%size%pt))
  end function coordsystem_node_inside

  function coordsystem_point_inside(coord,point)
    ! return true iff point is inside coord system
    implicit none
    type(coordsystem), intent(in) :: coord
    type(geom_point), intent(in) :: point
    logical coordsystem_point_inside
    integer i

    coordsystem_point_inside = .true.
    do i=1,geom_dim
       coordsystem_point_inside = (point%pt(i).ge.coord%origin%pt(i)) .and. &
            (point%pt(i).le.coord%origin%pt(i)+coord%size%pt(i)*coord%delta%pt(i))
       if (.not.coordsystem_point_inside) then
          exit
       end if
    end do
  end function coordsystem_point_inside
    
  function coordsystem_linearise2d(coord,node)
    ! linearise node, given coord
    implicit none
    type(coordsystem), intent(in) :: coord
    type(geom_ipoint), intent(in) :: node
    integer coordsystem_linearise2d

    coordsystem_linearise2d = -1

#ifdef DEBUG_COORDS
    if (.not.coordsystem_node_inside(coord,node)) then
       write(coordinates_unit,*) 'Error (',__FILE__,__LINE__,'): node (',node%pt,') not inside coord system'
       return
    end if
#endif
    
    coordsystem_linearise2d = node%pt(1) + (node%pt(2)-1)*coord%size%pt(1)
  end function coordsystem_linearise2d

  function coordsystem_delinearise2d(coord, ind)
    ! expand linearisation
    implicit none
    type(coordsystem), intent(in) :: coord
    integer, intent(in) :: ind
    type(geom_ipoint) :: coordsystem_delinearise2d

#ifdef DEBUG_COORDS
    if (ind.lt.1 .or. ind.gt.coord%size%pt(1)*coord%size%pt(2)) then
       write(coordinates_unit,*) 'Error (',__FILE__,__LINE__,'): index ',ind,' outside coord system'
       stop
    end if
#endif

    coordsystem_delinearise2d%pt(1) = mod(ind-1,coord%size%pt(1)) + 1
    coordsystem_delinearise2d%pt(2) = (ind-1)/coord%size%pt(1) + 1
  end function coordsystem_delinearise2d
end module glimmer_coordinates
