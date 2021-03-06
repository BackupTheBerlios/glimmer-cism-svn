! erosion_integrate2d.f90
! Magnus Hagdorn, March 2003
!
! integrate a function defined on a regular cartesian grid over arbitrary rectangular boundaries

module erosion_integrate2d
  use geometry

type(geom_poly), private :: cell
private :: addmat
real :: ismintegrate2d_zero = 0.0

contains
  subroutine init_integrate2d
    implicit none

    cell = poly_new(4)
  end subroutine init_integrate2d

  subroutine calc_weight(coords, weight, patch, n)
    use sparse
    use geometry
    use glimmer_coordinates
    implicit none
    type(coordsystem), intent(in) :: coords
    type(sparse_matrix) :: weight            ! sparse matrix
    type(geom_poly), intent(in) :: patch     ! the polygone over which will be integrated
    integer, intent(in) :: n                 ! cell number

    ! local variables
    type(geom_ipoint) :: node
    integer i
    integer imin, inum, pow, jmin, jnum,  num
    real w


    ! find the bounding box
    node = coordsystem_get_node(coords,patch%poly(1))
    imin = node%pt(1)
    inum = node%pt(1)
    jmin = node%pt(2)
    jnum = node%pt(2)
    do i=2,patch%n
       node = coordsystem_get_node(coords,patch%poly(i))
       if (node%pt(1).lt.imin) imin = node%pt(1)
       if (node%pt(1).gt.inum) inum = node%pt(1)
       if (node%pt(2).lt.jmin) jmin = node%pt(2)
       if (node%pt(2).gt.jnum) jnum = node%pt(2)
    end do
    node%pt(1) = imin 
    node%pt(2) = jmin 
    inum = inum - imin
    jnum = jnum - jmin

    num = max(inum,jnum)
    if (num.eq.0) then ! trivial case shape is within single cell
       w = 0.5*abs(poly_area2(patch))/(coords%delta%pt(1)*coords%delta%pt(2))
       call addmat(coords,weight,n,node,w)
       return
    end if
    
    pow = ceiling(log(real(num))/log(2.)) 
    if (pow.ge.1) then ! can we halve the area?
       pow = pow - 1
    end if

    ! non-trival case, start dividing
    num = 2**pow
    call divide_weight(coords,weight,patch,node,pow,n)
    node%pt(1) = node%pt(1)+num
    call divide_weight(coords,weight,patch,node,pow,n)
    node%pt(2) = node%pt(2)+num
    call divide_weight(coords,weight,patch,node,pow,n)
    node%pt(1) = node%pt(1)-num
    call divide_weight(coords,weight,patch,node,pow,n)
    node%pt(2) = node%pt(2)-num
  end subroutine calc_weight

  recursive subroutine divide_weight(coords,weight,patch,node,pow,n)
    use sparse
    use glimmer_coordinates
    implicit none
    type(coordsystem), intent(in) :: coords
    type(sparse_matrix) :: weight            ! sparse matrix
    type(geom_poly), intent(in) :: patch     ! the polygone over which will be integrated
    type(geom_ipoint),intent(inout) :: node  ! cell coordinates
    integer, intent(in) :: pow               ! size of square (numi,numj = 2**pow)
    integer, intent(in) :: n                 ! cell number

    ! local variables
    type(geom_poly) :: poly
    type(geom_point) :: point
    type(geom_ipoint) :: n2
    integer p
    logical inters
    integer i,j, i1,j1,next
    real val

    ! set up current cell coords
    cell%n=4
    point = coordsystem_get_coord(coords,node)
    point%pt(:) = point%pt(:) -0.5*coords%delta%pt(:)
    cell%poly(1) = point
    
    point%pt(1) = point%pt(1) + coords%delta%pt(1)*(2**pow)
    cell%poly(2) = point

    point%pt(2) = point%pt(2) + coords%delta%pt(2)*(2**pow)
    cell%poly(3) = point

    point%pt(1) = point%pt(1) - coords%delta%pt(1)*(2**pow)
    cell%poly(4) = point

    inters = .false.
    outer: do j=1,4 ! loop over box coords
       j1 = mod(j+2,4)+1
       do  i = 1,patch%n ! loop over patch coords
          i1 = mod(i+patch%n-2,patch%n)+1
          inters = intersectProp(cell%poly(j1),cell%poly(j), patch%poly(i1), patch%poly(i))
          if (inters) exit outer
       end do
    end do outer

    if (.not.inters) then ! patch does not intersect cell
       ! two possibilities
       if (point_in_cpoly(cell%poly(1),patch)) then  ! cell is inside patch
          do j = node%pt(2),node%pt(2)+2**pow-1
             do i = node%pt(1),node%pt(1)+2**pow-1
                n2%pt(1) = i
                n2%pt(2) = j
                call addmat(coords,weight,n,n2,1.0)
             end do
          end do
          return   
       else ! cell is outside patch
          return
       end if
    end if

    ! shape intersects box
    if (pow .gt. 0) then  ! we can still divide
       p = pow - 1
       next = 2**p
       call divide_weight(coords,weight,patch,node,p,n)
       node%pt(1) = node%pt(1) + next
       call divide_weight(coords,weight,patch,node,p,n)
       node%pt(2) = node%pt(2) + next
       call divide_weight(coords,weight,patch,node,p,n)
       node%pt(1) = node%pt(1) - next
       call divide_weight(coords,weight,patch,node,p,n)
       node%pt(2) = node%pt(2) - next
       return
    else ! we are left with a single cell now
       ! calculate intersection between cell and patch
       poly = ConvexIntersect(patch,cell)
       val = 0.5*real(poly_area2(poly)) / (coords%delta%pt(1)*coords%delta%pt(2))
       call addmat(coords,weight,n,node,val)
       call poly_delete(poly)
    end if
  end subroutine divide_weight

  subroutine addmat(coords,weight,n,node,val)
    ! adding val if i,j inside region, ignore otherwise
    use glimmer_coordinates
    use sparse
    implicit none
    type(coordsystem), intent(in) :: coords
    type(sparse_matrix) :: weight
    integer, intent(in) :: n
    type(geom_ipoint), intent(in) :: node
    real, intent(in) :: val

    if (coordsystem_node_inside(coords,node) .and. abs(val).gt.ismintegrate2d_zero) then
       call add_to_matrix(weight,n,coordsystem_linearise2d(coords,node),val)
    end if
  end subroutine addmat

end module erosion_integrate2d
