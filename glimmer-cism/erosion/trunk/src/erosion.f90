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
  end subroutine er_initialise

  subroutine er_tstep(erosion,model)
    !*FD do the erosion
    use isostasy
    use glide_types
    use physcon, only : rhom
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model       !*FD model instance
    
    ! local variables
    integer ew,ns

    call erosion_io_writeall(erosion,model)

    if (erosion%doerosion) then
       if (model%numerics%tinc .gt. mod(model%numerics%time,model%numerics%tinc*erosion%ndt)) then
          ! calculate erosion rate
          do ns=2,model%general%nsn-1
             do ew=2,model%general%ewn-1
                erosion%erosion_rate(ew,ns) = erosion%hb_erosion_factor * model%geometry%thck(ew,ns) &
                     * sqrt( sum(model%velocity%ubas(ew-1:ew,ns-1:ns))**2 + sum(model%velocity%vbas(ew-1:ew,ns-1:ns))**2 )
             end do
          end do
          erosion%er_accu = - erosion%erosion_rate * erosion%dt
          erosion%erosion = erosion%erosion + erosion%er_accu
          erosion%er_isos = erosion%er_isos + erosion%er_accu
          model%geometry%topg = model%geometry%topg + erosion%er_accu
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

  subroutine er_finalise(erosion)
    !*FD finalise erosion model
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data    
    call er_deallocate(erosion)
  end subroutine er_finalise
end module erosion
