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
  use erosion_transport
  use sparse

  type erosion_type
     logical :: doerosion = .False.                        !*FD set to true when erosion should be included
     integer :: ndt = 5                                    !*FD erosion time step (multiplier of main time step)
     real(kind=dp) :: dt                                   !*FD erosion time step
     real :: hb_erosion_factor =   1.e-10                  !*FD constant of proportionality for erosion rate calcs
     real :: density = 3000.                               !*FD density of hard bedrock (kg m$^{-3}$)
     ! sediment transport stuff
     logical :: dotransport = .False.                      !*FD set to true to move sediments about
     integer :: transport_ndt = 20                         !*FD transport time step (multiplier of main time step)
     real :: transport_dt                                  !*FD time step for recalculating sediment distribution
     real(kind=dp) :: transport_fac = 0.2                  !*FD multiplier for velos in deformable beds
     real :: dirty_ice_max = 0.1                           !*FD maximum thickness of dirty basal ice layer
     real :: soft_a = 1.e-5                                !*FD param A for max def thick calculations
     real :: soft_b = 0.                                   !*FD param B for max def thick calculations
     ! internal fields, etc
     type(er_transport_type) :: trans                      !*FD type holding transport stuff
     type(sparse_matrix) :: lag_seds1                      !*FD sparse matrix holding dirty ice layer
     type(sparse_matrix) :: lag_seds2                      !*FD sparse matrix holding deformable sediment layer
     real(kind=dp),dimension(:,:),pointer :: temporary => null()    !*FD temporary array
     real(kind=dp),dimension(:,:),pointer :: erosion_rate => null() !*FD hard bedrock erosion rate
     real(kind=dp),dimension(:,:),pointer :: erosion => null()      !*FD total hard bedrock erosion
     real(kind=dp),dimension(:,:),pointer :: er_accu => null()      !*FD accumulated erosion during one erosion time step
     real(kind=dp),dimension(:,:),pointer :: er_isos => null()      !*FD accumulated erosion for isostasy calcs
     real(kind=dp),dimension(:,:),pointer :: er_load => null()      !*FD load due to erosion
     real(kind=dp),dimension(:,:),pointer :: seds1 => null()        !*FD thickness of dirty basal ice layer
     real(kind=dp),dimension(:,:),pointer :: seds2 => null()        !*FD thickness of deforming sediment layer
     real(kind=dp),dimension(:,:),pointer :: seds2_max => null()    !*FD maximum thickness of deforming sediment layer
     real(kind=dp),dimension(:,:),pointer :: seds3 => null()        !*FD thickness of non-deforming sediment layer
  end type erosion_type

contains
  subroutine er_allocate(erosion,numx,numy)
    !*FD allocate erosion data
    implicit none
    type(erosion_type) :: erosion     !*FD data structure holding erosion stuff
    integer, intent(in) :: numx,numy  !*FD size of arrays

    allocate(erosion%temporary(numx,numy))
    allocate(erosion%erosion_rate(numx,numy))
    allocate(erosion%erosion(numx,numy))
    allocate(erosion%er_accu(numx,numy))
    allocate(erosion%er_isos(numx,numy))
    allocate(erosion%er_load(numx,numy))
    allocate(erosion%seds1(numx,numy))
    allocate(erosion%seds2(numx,numy))
    allocate(erosion%seds2_max(numx,numy))
    allocate(erosion%seds3(numx,numy))
  end subroutine er_allocate
    
  subroutine er_deallocate(erosion)
    !*FD free memory used by erosion data structure
    implicit none
    type(erosion_type) :: erosion     !*FD data structure holding erosion stuff

    deallocate(erosion%temporary)
    deallocate(erosion%erosion_rate)
    deallocate(erosion%erosion)
    deallocate(erosion%er_accu)
    deallocate(erosion%er_isos)
    deallocate(erosion%er_load)
    deallocate(erosion%seds1)
    deallocate(erosion%seds2)
    deallocate(erosion%seds2_max)
    deallocate(erosion%seds3)
  end subroutine er_deallocate
end module erosion_types
