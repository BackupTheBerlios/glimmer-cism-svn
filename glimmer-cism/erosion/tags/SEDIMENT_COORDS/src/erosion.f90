! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  erosion.f90 - part of the GLIMMER ice model              + 
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

module erosion
  !*FD the main erosion module
  use erosion_types
  use erosion_setup
  use erosion_io
  use erosion_advect
  use erosion_nc_custom
  use erosion_transport

contains
  subroutine er_initialise(erosion,config,model)
    !*FD initialise erosion model
    use glimmer_interpolate2d
    use glide_types
    use glimmer_config
    use erosion_sediment
    use paramets, only : len0,thk0
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file   
    type(glide_global_type) :: model       !*FD model instance

    ! read config
    call er_readconfig(erosion,config)
    ! print config
    call er_printconfig(erosion)

    ! setup transport grid
    erosion%ewn = model%general%ewn * erosion%grid_magnifier
    erosion%nsn = model%general%nsn * erosion%grid_magnifier
    erosion%dew = model%numerics%dew/real(erosion%grid_magnifier)
    erosion%dns = model%numerics%dns/real(erosion%grid_magnifier)
    erosion%coord = coordsystem_new(0.d0,0.d0,erosion%dew,erosion%dns, erosion%ewn,erosion%nsn)

    ! create erosion variables
    call erosion_io_createall(model,erosion)

    ! scale variables
    erosion%hb_erosion_factor = erosion%hb_erosion_factor*len0*thk0
    erosion%dt = erosion%ndt * model%numerics%dt
    erosion%soft_a = erosion%soft_a*thk0*thk0/len0

    ! allocate memory
    call er_allocate(erosion,model)
    erosion%erosion = 0.

    ! initialise transport
    call init_transport(erosion%trans, model,erosion)

    ! initialise sediment layer if not simple
    if (.not.erosion%simple_seds) then
       call er_sediment_init(erosion%sediment,model)
    end if

    ! read variables
    call erosion_io_readall(erosion,model)

    ! setup interpolation between velo grid and sediment grid
    call glimmer_init_bilinear(model%general%velo_grid, erosion%coord, erosion%velo_seds)

    ! initialise profile
#ifdef PROFILING
    call erosion_prof_init(model,erosion)
#endif
  end subroutine er_initialise

  subroutine er_tstep(erosion,model)
    !*FD do the erosion
    use isostasy
    use glide_types
    use physcon, only : rhom
    use erosion_transport
    use glide_profile
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model       !*FD model instance
    
    ! local variables
    integer ew,ns
    type(coord_ipoint) :: node
    real(kind=dp) :: dummy

    call erosion_io_writeall(erosion,model)

    if (erosion%doerosion) then
       
       if (model%numerics%tinc .gt. mod(model%numerics%time,model%numerics%tinc*erosion%ndt)) then

          !----------------------------------------------------------
          ! calculate erosion rate
          !----------------------------------------------------------
#ifdef PROFILING
          call glide_prof_start(model,erosion%er_prof%erate)
#endif      
          call er_calc_erate(erosion,model)
          ! interpolate
          call er_interpolate(erosion,model,erosion%erosion_rate,erosion%er_accu)
          erosion%er_accu = erosion%er_accu * erosion%dt
#ifdef PROFILING
    call glide_prof_stop(model,erosion%er_prof%erate)
#endif

          if (erosion%dotransport) then
             !----------------------------------------------------------
             ! transport sediments
             !----------------------------------------------------------
#ifdef PROFILING
             call glide_prof_start(model,erosion%er_prof%trans_sed)
#endif      
             call transport_sediments(erosion,model)
#ifdef PROFILING
             call glide_prof_stop(model,erosion%er_prof%trans_sed)
#endif
             
             !------------------------------------------------
             ! sediment erosion
             !------------------------------------------------       
             ! add sediment to deformable soft bed from below
#ifdef PROFILING
             call glide_prof_start(model,erosion%er_prof%sed_eros)
#endif      
             call  er_calc_dthick(erosion,model)
             do ns=2,erosion%nsn-1
                do ew=2,erosion%ewn-1
                   if (erosion%seds2(ew,ns) .lt. erosion%seds2_max(ew,ns)) then
                      dummy = erosion%seds2_max(ew,ns) - erosion%seds2(ew,ns)
                      if (dummy .gt. erosion%seds3(ew,ns)) then
                         erosion%seds2(ew,ns) = erosion%seds2(ew,ns) + erosion%seds3(ew,ns)
                         erosion%seds3(ew,ns) = 0.
                      else
                         erosion%seds2(ew,ns) = erosion%seds2(ew,ns) + dummy
                         erosion%seds3(ew,ns) = erosion%seds3(ew,ns) - dummy
                      end if
                   end if
                end do
             end do

             do ns=2,erosion%nsn-1
                do ew=2,erosion%ewn-1
                   ! initially all eroded material goes into basal dirty ice layer
                   erosion%seds1(ew,ns) = erosion%seds1(ew,ns) + erosion%er_accu(ew,ns)
                   dummy = erosion%seds2(ew,ns) + erosion%seds3(ew,ns)
                   ! erode from deformable bed
                   if (erosion%seds2(ew,ns) .ge. erosion%er_accu(ew,ns)) then
                      ! deformable layer is thick enough
                      erosion%seds2(ew,ns) = erosion%seds2(ew,ns) - erosion%er_accu(ew,ns)
                   else
                      ! need to dig into stationary sediment layer
                      erosion%seds3(ew,ns) = max(0.d0,erosion%seds3(ew,ns) + erosion%seds2(ew,ns) - erosion%er_accu(ew,ns))
                      erosion%seds2(ew,ns) = 0.
                   end if
                   erosion%erosion(ew,ns) = erosion%erosion(ew,ns) -min(0.d0,dummy - erosion%er_accu(ew,ns))
                end do
             end do
#ifdef PROFILING
             call glide_prof_stop(model,erosion%er_prof%sed_eros)
#endif

             !------------------------------------------------
             ! sediment deposition
             !------------------------------------------------
#ifdef PROFILING
             call glide_prof_start(model,erosion%er_prof%sed_dep)
#endif      
             do ns=2,erosion%nsn-1
                do ew=2,erosion%ewn-1
                   ! from dirty basal ice layer
                   if (erosion%seds1(ew,ns) .gt. erosion%dirty_ice_max) then
                      erosion%seds2(ew,ns) = erosion%seds2(ew,ns) + erosion%seds1(ew,ns) - erosion%dirty_ice_max
                      erosion%seds1(ew,ns) = erosion%dirty_ice_max
                   end if
                   ! from deforming soft bed to non-deforming soft bed
                   if (erosion%seds2(ew,ns) .gt. erosion%seds2_max(ew,ns)) then
                      erosion%seds3(ew,ns) = erosion%seds3(ew,ns) + erosion%seds2(ew,ns) - erosion%seds2_max(ew,ns)
                      erosion%seds2(ew,ns) = erosion%seds2_max(ew,ns)
                   end if
                   ! depositing debris in basal layer when ice has retreated
                   node%pt(1) = ew
                   node%pt(2) = ns
                   node = coordsystem_get_node(model%general%ice_grid, coordsystem_get_coord(erosion%coord, node))
                   if (model%geometry%thck(node%pt(1),node%pt(2)).eq.0) then
                      erosion%seds3(ew,ns) = erosion%seds3(ew,ns) + erosion%seds2(ew,ns) + erosion%seds1(ew,ns)
                      erosion%seds2(ew,ns) = 0.
                      erosion%seds1(ew,ns) = 0.
                   end if
                end do
             end do
#ifdef PROFILING
             call glide_prof_stop(model,erosion%er_prof%sed_dep)
#endif

             !erosion%er_isos = erosion%er_isos + erosion%er_accu
             !model%geometry%topg = model%geometry%topg + erosion%er_accu
          else
             !------------------------------------------------
             ! no sediment transport 
             ! eroded sediments are lost instantly
             !------------------------------------------------

             erosion%erosion(:,:) = erosion%erosion(:,:) - erosion%er_accu(:,:)
             
          end if
       end if
    end if

!!$    ! update load if necessary
!!$    if (model%isos%new_load) then
!!$       model%isos%relx = model%isos%relx + erosion%er_isos
!!$       erosion%er_isos = erosion%er_isos * erosion%density/rhom
!!$       call isos_lithosphere(model,erosion%er_load,erosion%er_isos)
!!$       model%isos%relx = model%isos%relx - erosion%er_load
!!$       erosion%er_isos = 0.
!!$    end if

  end subroutine er_tstep

  subroutine er_calc_erate(erosion,model)
    !*FD calculate simple erosions rate
    use glide_types
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model       !*FD model instance
    
    ! local variables
    integer ew,ns
    do ns=1,model%general%nsn-1
       do ew=1,model%general%ewn-1
          erosion%erosion_rate(ew,ns) = erosion%hb_erosion_factor * model%geomderv%stagthck(ew,ns) &
               * sqrt(model%velocity%ubas(ew,ns)**2 + model%velocity%vbas(ew,ns)**2)
       end do
    end do
  end subroutine er_calc_erate

  subroutine er_finalise(erosion)
    !*FD finalise erosion model
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data    
    call er_deallocate(erosion)
  end subroutine er_finalise
end module erosion
