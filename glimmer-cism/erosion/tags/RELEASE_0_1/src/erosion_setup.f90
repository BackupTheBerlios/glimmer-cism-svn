! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  erosion_setup.f90 - part of the GLIMMER ice model        + 
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

module erosion_setup
contains
  subroutine er_readconfig(erosion,config)
    !*FD read erosion configuration
    use erosion_types
    use glimmer_config
    implicit none
    type(erosion_type) :: erosion           !*FD structure holding erosion data
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   
    ! local variables
    type(ConfigSection), pointer :: section

    call GetSection(config,section,'Erosion')
    if (associated(section)) then
       erosion%doerosion = .True.
       call GetValue(section,'hb_erosion',erosion%hb_erosion_factor)
       call GetValue(section,'ntime',erosion%ndt)
    end if
  end subroutine er_readconfig

  subroutine er_printconfig(erosion)
    !*FD print erosion config to log
    use glimmer_log
    use erosion_types
    implicit none
    type(erosion_type) :: erosion           !*FD structure holding erosion data
    ! local variables
    character(len=100) :: message

    if (erosion%doerosion) then
       call write_log('Erosion')
       call write_log('-------')
       write(message,*) 'Updating erosion every ',erosion%ndt,' time steps'
       call write_log(message)
       write(message,*) 'hard bedrock erosion constant : ',erosion%hb_erosion_factor
       call write_log(message)
       call write_log('')
    end if
  end subroutine er_printconfig

end module erosion_setup
