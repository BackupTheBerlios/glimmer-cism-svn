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
    use paramets, only : acc0
    use physcon, only : scyr
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
    erosion%hb_erosion_factor = erosion%hb_erosion_factor/(acc0*scyr)

    ! allocate memory
    call er_allocate(erosion,model%general%ewn,model%general%nsn)
    erosion%erosion = 0.
  end subroutine er_initialise

  subroutine er_tstep(erosion,model)
    !*FD do the erosion
    use glide_types
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model       !*FD model instance
    
    ! local variables
    integer ew,ns

    call erosion_io_writeall(erosion,model)

    if (erosion%doerosion) then

       ! calculate erosion rate
       do ns=2,model%general%nsn-1
          do ew=2,model%general%ewn-1
             erosion%erosion_rate(ew,ns) = erosion%hb_erosion_factor * model%geometry%thck(ew,ns) &
                  * sqrt( sum(model%velocity%ubas(ew-1:ew,ns-1:ns))**2 + sum(model%velocity%vbas(ew-1:ew,ns-1:ns))**2 )
          end do
       end do
          
       erosion%erosion = erosion%erosion - erosion%erosion_rate * model%numerics%dt

       model%geometry%topg = model%geometry%topg - erosion%erosion_rate * model%numerics%dt
    end if

  end subroutine er_tstep

  subroutine er_finalise(erosion)
    !*FD finalise erosion model
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data    
    call er_deallocate(erosion)
  end subroutine er_finalise
end module erosion
