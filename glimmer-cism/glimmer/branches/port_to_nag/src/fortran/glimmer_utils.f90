
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_utils.f90 - part of the GLIMMER ice model        + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004 GLIMMER contributors - see COPYRIGHT file 
! for list of contributors.
!
! This program is free software; you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation; either version 2 of 
! the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License 
! along with this program; if not, write to the Free Software 
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
! 02111-1307 USA
!
! GLIMMER is maintained by:
!
! Ian Rutt
! School of Geographical Sciences
! University of Bristol
! University Road
! Bristol
! BS8 1SS
! UK
!
! email: <i.c.rutt@bristol.ac.uk> or <ian.rutt@physics.org>
!
! GLIMMER is hosted on NeSCForge:
!
! http://forge.nesc.ac.uk/projects/glimmer/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glimmer_utils

  !*FD Module containing utility code for GLIMMER.

  use glimmer_global

  implicit none

  interface array_bcs
    module procedure array_bcs1d,array_bcs2d
  end interface

  interface check_conformal
    module procedure check_conformal_2d_real
  end interface

contains

  real(rk) function array_bcs1d(array,i)

    !*FD Returns the value of a 1D array
    !*FD location,checking first for the boundaries.
    !*FD Note that this function is aliased as {\tt array\_bcs}.
    !*FDRV The value of the location in question.

    ! Arguments

    real(rk),dimension(:),intent(in) :: array !*FD The array to be indexed.
    integer,intent(in)               :: i     !*FD The location to be extracted.

    ! Internal variables

    integer :: n,ii

    n=size(array)
    ii=i

    if ((i<=n).and.(i>=1)) then
      array_bcs1d=array(i)
    endif

    do while (ii>n)
      ii=ii-n
    enddo

    do while (ii<1)
      ii=ii+n
    enddo

    array_bcs1d=array(ii)

  end function array_bcs1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(rk) function array_bcs_lats(array,i)

    !*FD As {\tt array\_bcs1d}, but adapted
    !*FD for dealing with polar boundary conditions.
    !*FDRV The value of the location in question.

    ! Arguments

    real(rk),dimension(:),intent(in) :: array !*FD The array to be indexed.
    integer,intent(in) :: i !*FD The location to be extracted.

    ! Internal variables

    integer :: n,ii

    n=size(array)
    ii=i

    if ((i<=n).and.(i>=1)) then
      array_bcs_lats=array(i)
      return
    endif

    if (ii>n) then
      ii=2*n-ii
      array_bcs_lats=-180.0+array(ii)
    endif

    if (ii<1) then
      ii=1-ii
      array_bcs_lats=180.0-array(ii)
    endif

  end function array_bcs_lats

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(rk) function array_bcs2d(array,i,j)

    !*FD Returns the value of an array 
    !*FD location, checking first for the boundaries. 
    !*FD Over-the-pole boundary conditions are implemented here.
    !*FDRV The value of the location specified.

    ! Arguments

    real(rk),dimension(:,:),intent(in) :: array !*FD Array to be indexed
    integer,intent(in) :: i,j !*FD The location to be extracted

    ! Internal variables

    integer :: nx,ny,ii,jj

    nx=size(array,1) ; ny=size(array,2)

    if ((i>=1).and.(i<=nx).and.(j>=1).and.(j<=ny)) then
      array_bcs2d=array(i,j)
      return
    endif

    ii=i ; jj=j

    if (jj>ny) then
      jj=2*ny-jj
      ii=ii+nx/2
    endif

    if (jj<1) then
      jj=1-jj
      ii=ii+nx/2
    endif

    do while (ii>nx) 
      ii=ii-nx
    enddo

    do while (ii<1)
      ii=ii+nx
    enddo  

    array_bcs2d=array(ii,jj)

  end function array_bcs2d

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine fix_bcs2d(i,j,nx,ny)

    !*FD Adjusts array location indicies
    !*FD so that they fall within the domain.

    integer,intent(inout) :: i,j !*FD The location of interest
    integer,intent(in) :: nx,ny  !*FD The size of the domain (number
                                 !*FD of points in each direction)

    if ((i>=1).and.(i<=nx).and.(j>=1).and.(j<=ny)) return

    if (j>ny) then
      j=2*ny-j
      i=i+nx/2
    endif

    if (j<1) then
      j=1-j
      i=i+nx/2
    endif

    do while (i>nx) 
      i=i-nx
    enddo

    do while (i<1)
      i=i+nx
    enddo  

  end subroutine fix_bcs2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine check_conformal_2d_real(array1,array2,label)

    !*FD Checks that two arrays
    !*FD are of the same size. Aliased as {\tt check\_conformal}.

    real(rk),dimension(:,:),intent(in) :: array1,array2 !*FD The arrays to be checked
    character(*),intent(in),optional :: label !*FD Optional label, to facilitate 
                                              !*FD bug tracking if the check fails.

    if ((size(array1,1)/=size(array2,1)).or.(size(array1,2)/=size(array2,2))) then
      if (present(label)) then
        print*,'ERROR: Non-conformal arrays. Label: ',label
      else
        print*,'ERROR: Non-conformal arrays. No label'
      endif
      stop
    endif

  end subroutine check_conformal_2d_real

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function hsum(inp)

    !*FD Calculates the sum of a given three-dimensional field at each
    !*FD level. The vertical coordinate of the input is the first index of
    !*FD the array.
    !*FDRV A one-dimensional array of the same size as the first dimension of
    !*FDRV \texttt{inp} is returned, containing the sum of \texttt{inp} for 
    !*FDRV each level.

    implicit none

    real(dp),dimension(:,:,:),intent(in) :: inp !*FD The input array. The first
                                                !*FD index is the vertical, the other
                                                !*FD two horizontal.
    real(dp),dimension(size(inp,dim=1))  :: hsum
  
    hsum = sum(sum(inp(:,:,:),dim=3),dim=2)

  end function hsum

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function lsum(inp)

    !*FD Calculates the sum of a given two-dimensional field along one axis.
    !*FD Within GLIMMER, this function calculates the mean vertical profile
    !*FD in a 2D vertical slice. 
    !*FD A one-dimensional array of the same size as the first dimension of
    !*FDRV \texttt{inp} is returned, containing the sum of \texttt{inp} for 
    !*FDRV each row.

    implicit none

    real(dp),dimension(:,:), intent(in) :: inp
    real(dp),dimension(size(inp,dim=1)) :: lsum
    
    lsum = sum(inp(:,:),dim=2)

  end function lsum

end module glimmer_utils
