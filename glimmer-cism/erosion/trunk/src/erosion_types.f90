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

  use glimmer_global, only : dp  

  type erosion_type
     logical :: doerosion = .False.                        !*FD set to true when erosion should be included
     integer :: ndt = 20                                   !*FD erosion time step (multiplier of main time step)
     real :: dt                                            !*FD erosion time step
     real :: hb_erosion_factor =   1.e-10                  !*FD constant of proportionality for erosion rate calcs
     real :: density = 3000.                               !*FD density of hard bedrock (kg m$^{-3}$)
     real(kind=dp),dimension(:,:),pointer :: erosion_rate => null() !*FD hard bedrock erosion rate
     real(kind=dp),dimension(:,:),pointer :: erosion => null()      !*FD total hard bedrock erosion
     real(kind=dp),dimension(:,:),pointer :: er_accu => null()      !*FD accumulated erosion during one erosion time step
     real(kind=dp),dimension(:,:),pointer :: er_isos => null()      !*FD accumulated erosion for isostasy calcs
     real(kind=dp),dimension(:,:),pointer :: er_load => null()      !*FD load due to erosion
  end type erosion_type

contains
  subroutine er_allocate(erosion,numx,numy)
    !*FD allocate erosion data
    implicit none
    type(erosion_type) :: erosion     !*FD data structure holding erosion stuff
    integer, intent(in) :: numx,numy  !*FD size of arrays

    allocate(erosion%erosion_rate(numx,numy))
    allocate(erosion%erosion(numx,numy))
    allocate(erosion%er_accu(numx,numy))
    allocate(erosion%er_isos(numx,numy))
    allocate(erosion%er_load(numx,numy))
  end subroutine er_allocate
    
  subroutine er_deallocate(erosion)
    !*FD free memory used by erosion data structure
    implicit none
    type(erosion_type) :: erosion     !*FD data structure holding erosion stuff

    deallocate(erosion%erosion_rate)
    deallocate(erosion%erosion)
    deallocate(erosion%er_accu)
    deallocate(erosion%er_isos)
    deallocate(erosion%er_load)
  end subroutine er_deallocate
end module erosion_types
