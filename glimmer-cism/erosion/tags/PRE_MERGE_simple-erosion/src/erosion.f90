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

contains
  subroutine er_initialise(erosion,config,model)
    !*FD initialise erosion model
    use glide_types
    use glimmer_config
    use paramets, only : len0
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file   
    type(glide_global_type) :: model       !*FD model instance

    ! read config
    call er_readconfig(erosion,config)
    ! print config
    call er_printconfig(erosion)
    ! create erosion variables
    call erosion_io_createall(model)

    ! scale variables
    erosion%hb_erosion_factor = erosion%hb_erosion_factor*len0
    erosion%dt = erosion%ndt * model%numerics%dt
    erosion%transport_dt = erosion%transport_ndt * erosion%transport_dt

    ! allocate memory
    call er_allocate(erosion,model%general%ewn,model%general%nsn)
    erosion%erosion = 0.

    ! initialise transport
    call init_transport(erosion%trans, model)
    ! initialise sparse matrices
    erosion%lag_seds1 = new_sparse_matrix(10000)
    erosion%lag_seds2 = new_sparse_matrix(10000)
  end subroutine er_initialise

  subroutine er_tstep(erosion,model)
    !*FD do the erosion
    use isostasy
    use glide_types
    use physcon, only : rhom
    use paramets, only : vel0
    use erosion_advect
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
             call set_velos(model%velocity%ubas,model%velocity%vbas,-vel0)
             call calc_lagrange(model, erosion%trans, erosion%dt, erosion%lag_seds1)
             ! transport in deformable sediment layer
             call set_velos(model%velocity%ubas,model%velocity%vbas,-erosion%transport_fac*vel0)
             call calc_lagrange(model, erosion%trans, erosion%dt, erosion%lag_seds2)      
          end if
       end if


       if (model%numerics%tinc .gt. mod(model%numerics%time,model%numerics%tinc*erosion%ndt)) then
          !----------------------------------------------------------
          ! move sediments
          !----------------------------------------------------------
          ! move sediments in basal ice layer
          call transport_scalar(model,erosion%trans,erosion%seds1,erosion%lag_seds1)
          ! move sediments in deformable sed layer
          call transport_scalar(model,erosion%trans,erosion%seds2,erosion%lag_seds2)


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
          ! initially all eroded material goes into basal dirty ice layer
          erosion%seds1 = erosion%seds1 + erosion%er_accu
          do ns=2,model%general%nsn-1
             do ew=2,model%general%ewn-1
                if (erosion%seds2(ew,ns) .gt. 0.) then
                   ! erode from deformable bed
                   if (erosion%seds2(ew,ns) .gt. erosion%er_accu(ew,ns)) then
                      ! deformable layer is thick enough
                      erosion%seds2(ew,ns) = erosion%seds2(ew,ns) - erosion%er_accu(ew,ns)
                   else
                      ! need to dig into stationary sediment layer
                      dummy = erosion%er_accu(ew,ns) - erosion%seds2(ew,ns)
                      erosion%seds2(ew,ns) = 0.
                      if (erosion%seds3(ew,ns) .gt. dummy) then
                         ! all sediments come from stationary sediment layer
                         erosion%seds3(ew,ns) = erosion%seds3(ew,ns) - dummy
                      else
                         erosion%seds3(ew,ns) = 0.
                      end if
                   end if
                end if
             end do
          end do
          erosion%er_accu = - max(0.d0,erosion%er_accu-erosion%seds2-erosion%seds3)
          erosion%erosion = erosion%erosion + erosion%er_accu

          !------------------------------------------------
          ! sediment deposition
          !------------------------------------------------
          ! from dirty basal ice layer
          where (erosion%seds1 .gt. erosion%dirty_ice_max)
             erosion%seds2 = erosion%seds2 + erosion%seds1 - erosion%dirty_ice_max
             erosion%seds1 = erosion%dirty_ice_max
          end where
          ! from deforming soft bed to non-deforming soft bed
          where (erosion%seds2 .gt. erosion%seds2_max)
             erosion%seds3 = erosion%seds3 + erosion%seds2 - erosion%seds2_max
             erosion%seds2 = erosion%seds2_max
          end where
          ! depositing debris in basal layer when ice has retreated
          where (model%geometry%thck.eq.0)
             erosion%seds3 = erosion%seds3 + erosion%seds2 + erosion%seds1
             erosion%seds2 = 0.
             erosion%seds1 = 0.
          end where


          !erosion%er_isos = erosion%er_isos + erosion%er_accu
          !model%geometry%topg = model%geometry%topg + erosion%er_accu
       end if
    end if

    ! update load if necessary
    if (model%isos%new_load) then
       model%isos%relx = model%isos%relx + erosion%er_isos
       erosion%er_isos = erosion%er_isos * erosion%density/rhom
       call isos_lithosphere(model,erosion%er_load,erosion%er_isos)
       model%isos%relx = model%isos%relx - erosion%er_load
       erosion%er_isos = 0.
    end if

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
    use paramets, only : thk0,len0
    use glide_types
    use glide_velo
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model       !*FD model instance

    real(kind=dp) :: dummy
    real(kind=dp) :: factor = thk0*thk0/len0

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
             erosion%seds2_max(ew,ns) = erosion%soft_a*sqrt(erosion%seds2_max(ew,ns))*factor + erosion%soft_b
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
