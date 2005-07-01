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
  use glimmer_sparse
  use glimmer_coordinates
  use erosion_transport_type

  type er_prof_type
     integer :: erate
     integer :: calc_lag
     integer :: trans_sed
     integer :: sed_eros
     integer :: sed_dep
  end type er_prof_type

  type erosion_type
     logical :: doerosion = .False.                        !*FD set to true when erosion should be included
     integer :: ndt = 5                                    !*FD erosion time step (multiplier of main time step)
     real(kind=dp) :: dt                                   !*FD erosion time step
     real(kind=dp) :: hb_erosion_factor =   1.d-10         !*FD constant of proportionality for erosion rate calcs
     real :: density = 3000.                               !*FD density of hard bedrock (kg m$^{-3}$)
     ! sediment transport stuff
     logical :: dotransport = .False.                      !*FD set to true to move sediments about
     real(kind=dp) :: transport_fac = 0.2                  !*FD multiplier for velos in deformable beds
     real(kind=dp) :: dirty_ice_max = 0.1                  !*FD maximum thickness of dirty basal ice layer
     real(kind=dp) :: soft_a = 1.d-5                       !*FD param A for max def thick calculations
     real(kind=dp) :: soft_b = 0.                          !*FD param B for max def thick calculations
     ! internal fields, etc
     type(er_transport_type) :: trans                      !*FD type holding transport stuff
     real(kind=dp),dimension(:,:),pointer :: erosion_rate => null() !*FD hard bedrock erosion rate
     real(kind=dp),dimension(:,:),pointer :: erosion => null()      !*FD total hard bedrock erosion
     real(kind=dp),dimension(:,:),pointer :: er_accu => null()      !*FD accumulated erosion during one erosion time step
     real(kind=dp),dimension(:,:),pointer :: er_isos => null()      !*FD accumulated erosion for isostasy calcs
     real(kind=dp),dimension(:,:),pointer :: er_load => null()      !*FD load due to erosion
     real(kind=dp),dimension(:,:),pointer :: seds1 => null()        !*FD thickness of dirty basal ice layer
     real(kind=dp),dimension(:,:),pointer :: seds2 => null()        !*FD thickness of deforming sediment layer
     real(kind=dp),dimension(:,:),pointer :: seds2_max => null()    !*FD maximum thickness of deforming sediment layer
     real(kind=dp),dimension(:,:),pointer :: seds2_max_v => null()  !*FD maximum thickness of deforming sediment layer on velocity grid
     real(kind=dp),dimension(:,:),pointer :: seds3 => null()        !*FD thickness of non-deforming sediment layer
     ! sediment grid
     integer :: grid_magnifier = 2                         !*FD increase sediment grid resolution by this factor
     integer :: ewn,nsn                                    !*FD number of nodes in x and y dir
     real(kind=dp) :: dew,dns                              !*FD grid spacing
     type(coordsystem_type) :: coord                       !*FD coordinate system of sed grid
     type(sparse_matrix_type) :: velo_seds                 !*FD transformation from velo to seds grid
     ! profiling
     type(er_prof_type) :: er_prof
  end type erosion_type

contains
  subroutine er_allocate(erosion,model)
    !*FD allocate erosion data
    use glide_types
    implicit none
    type(erosion_type) :: erosion     !*FD data structure holding erosion stuff
    type(glide_global_type) :: model  !*FD model instance

    call coordsystem_allocate(model%general%ice_grid,  erosion%er_isos)
    call coordsystem_allocate(model%general%ice_grid,  erosion%er_load)
    call coordsystem_allocate(model%general%velo_grid, erosion%erosion_rate)
    call coordsystem_allocate(model%general%velo_grid, erosion%seds2_max_v)
    call coordsystem_allocate(erosion%coord, erosion%erosion)
    call coordsystem_allocate(erosion%coord, erosion%er_accu)
    call coordsystem_allocate(erosion%coord, erosion%seds1)
    call coordsystem_allocate(erosion%coord, erosion%seds2)
    call coordsystem_allocate(erosion%coord, erosion%seds2_max)
    call coordsystem_allocate(erosion%coord, erosion%seds3)
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
    deallocate(erosion%seds1)
    deallocate(erosion%seds2)
    deallocate(erosion%seds2_max)
    deallocate(erosion%seds2_max_v)
    deallocate(erosion%seds3)
  end subroutine er_deallocate
end module erosion_types
