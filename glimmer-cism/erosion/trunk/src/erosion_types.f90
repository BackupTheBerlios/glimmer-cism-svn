! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  erosion_types.f90 - part of the GLIMMER ice model        + 
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

module erosion_types
  !*FD type definition for erosion calcluations

  type erosion_type
     logical :: doerosion = .False.                        !*FD set to true when erosion should be included
     real :: hb_erosion_factor =   1.e-7                   !*FD constant of proportionality for erosion rate calcs
     real,dimension(:,:),pointer :: erosion_rate => null() !*FD hard bedrock erosion rate
     real,dimension(:,:),pointer :: erosion => null()      !*FD total hard bedrock erosion
  end type erosion_type

contains
  subroutine er_allocate(erosion,numx,numy)
    !*FD allocate erosion data
    implicit none
    type(erosion_type) :: erosion     !*FD data structure holding erosion stuff
    integer, intent(in) :: numx,numy  !*FD size of arrays

    allocate(erosion%erosion_rate(numx,numy))
    allocate(erosion%erosion(numx,numy))
  end subroutine er_allocate
    
  subroutine er_deallocate(erosion)
    !*FD free memory used by erosion data structure
    implicit none
    type(erosion_type) :: erosion     !*FD data structure holding erosion stuff

    deallocate(erosion%erosion_rate)
    deallocate(erosion%erosion)
  end subroutine er_deallocate
end module erosion_types
