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
    use glide_types
    use glimmer_config
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

    ! create erosion variables
    call erosion_io_createall(model,erosion)

    ! scale variables
    erosion%hb_erosion_factor = erosion%hb_erosion_factor*len0*thk0
    erosion%dt = erosion%ndt * model%numerics%dt
    erosion%transport_dt = erosion%transport_ndt * erosion%transport_dt
    erosion%soft_a = erosion%soft_a*thk0*thk0/len0

    ! allocate memory
    call er_allocate(erosion,model%general%ewn,model%general%nsn)
    erosion%erosion = 0.

    ! initialise transport
    call init_transport(erosion%trans, model,erosion)
    ! initialise sparse matrices
    call new_sparse_matrix(10000,erosion%lag_seds1)
    call new_sparse_matrix(10000,erosion%lag_seds2)

    ! read variables
    call erosion_io_readall(erosion,model)

  end subroutine er_initialise

  subroutine er_tstep(erosion,model)
    !*FD do the erosion
    use isostasy
    use glide_types
    use physcon, only : rhom
    use paramets, only : vel0
    use erosion_advect
    use erosion_transport
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model       !*FD model instance
    
    ! local variables
    integer ew,ns
    real(kind=dp) :: dummy

    call erosion_io_writeall(erosion,model)



    if (erosion%doerosion) then
       
       !----------------------------------------------------------
       ! set up sediment transport sparse matrix
       !----------------------------------------------------------
       if (erosion%dotransport) then
          ! update transport matrix
          if (model%numerics%tinc .gt. mod(model%numerics%time,model%numerics%tinc*erosion%transport_ndt)) then
             ! transport in ice base
             call set_velos(model%velocity%ubas,model%velocity%vbas,-1.d0)
             call calc_lagrange(erosion, erosion%trans, erosion%dt, erosion%lag_seds1)
             ! transport in deformable sediment layer
             call set_velos(model%velocity%ubas,model%velocity%vbas,-erosion%transport_fac)
             call calc_lagrange(erosion, erosion%trans, erosion%dt, erosion%lag_seds2)      
          end if
       end if

       if (model%numerics%tinc .gt. mod(model%numerics%time,model%numerics%tinc*erosion%ndt)) then

          !----------------------------------------------------------
          ! move sediments
          !----------------------------------------------------------
          ! move sediments in basal ice layer
          call transport_scalar(erosion,erosion%trans,erosion%seds1,erosion%lag_seds1)
          ! move sediments in deformable sed layer
          call transport_scalar(erosion,erosion%trans,erosion%seds2,erosion%lag_seds2)


          !------------------------------------------------
          ! sediment erosion
          !------------------------------------------------       
          ! add sediment to deformable soft bed from below
          call  er_calc_dthick(erosion,model)
          do ns=2,model%general%nsn-1
             do ew=2,model%general%ewn-1
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
          ! calculate erosion rate
          call er_calc_erate(erosion,model)
          erosion%er_accu = erosion%erosion_rate * erosion%dt
          do ns=2,model%general%nsn-1
             do ew=2,model%general%ewn-1
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
                erosion%er_accu(ew,ns) = - max(0.d0,erosion%er_accu(ew,ns)-dummy)
                erosion%erosion(ew,ns) = erosion%erosion(ew,ns) + erosion%er_accu(ew,ns)
             end do
          end do

          !------------------------------------------------
          ! sediment deposition
          !------------------------------------------------
          do ns=2,model%general%nsn-1
             do ew=2,model%general%ewn-1
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
                if (model%geometry%thck(ew,ns).eq.0) then
                   erosion%seds3(ew,ns) = erosion%seds3(ew,ns) + erosion%seds2(ew,ns) + erosion%seds1(ew,ns)
                   erosion%seds2(ew,ns) = 0.
                   erosion%seds1(ew,ns) = 0.
                end if
             end do
          end do

          !erosion%er_isos = erosion%er_isos + erosion%er_accu
          !model%geometry%topg = model%geometry%topg + erosion%er_accu
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
    do ns=2,model%general%nsn-1
       do ew=2,model%general%ewn-1
          erosion%erosion_rate(ew,ns) = erosion%hb_erosion_factor * model%geometry%thck(ew,ns) &
               * sqrt( sum(model%velocity%ubas(ew-1:ew,ns-1:ns))**2 + sum(model%velocity%vbas(ew-1:ew,ns-1:ns))**2 )
       end do
    end do
  end subroutine er_calc_erate

  subroutine er_calc_dthick(erosion,model)
    !*FD calculate thickness of deformable sediment bed
    use glimmer_global, only: dp
    use glide_types
    use glide_velo
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model       !*FD model instance

    real(kind=dp) :: dummy

    integer ew,ns

    ! calculate basal shear stresses
    call calc_basal_shear(model)

    do ns=2,model%general%nsn-1
       do ew=2,model%general%ewn-1
          if (erosion%erosion_rate(ew,ns).gt.0) then
             dummy = 0.25*sum(model%velocity%tau_x(ew-1:ew,ns-1:ns))
             erosion%seds2_max(ew,ns) = dummy*dummy
             dummy = 0.25*sum(model%velocity%tau_y(ew-1:ew,ns-1:ns))
             erosion%seds2_max(ew,ns) = erosion%seds2_max(ew,ns) + dummy*dummy
             erosion%seds2_max(ew,ns) = erosion%soft_a*sqrt(erosion%seds2_max(ew,ns)) + erosion%soft_b
          else
             erosion%seds2_max(ew,ns) = 0.
          end if
       end do
    end do
  end subroutine er_calc_dthick

  subroutine er_finalise(erosion)
    !*FD finalise erosion model
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data    
    call er_deallocate(erosion)
  end subroutine er_finalise
end module erosion
