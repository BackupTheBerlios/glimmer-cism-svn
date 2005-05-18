! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  isostasy_setup.f90 - part of the GLIMMER ice model       + 
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

module isostasy_setup

  !*FD routines for setting up isostasy module

  use isostasy_types
  use glimmer_global, only : dp

contains
  
  subroutine isos_readconfig(isos,config)
    !*FD read isostasy configuration
    use glimmer_config
    implicit none
    type(isos_type) :: isos                !*FD structure holding isostasy configuration
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file

    ! local variables
    type(ConfigSection), pointer :: section

    ! read isostasy section
    call GetSection(config,section,'isostasy')
    if (associated(section)) then
       isos%do_isos = .true.
       call GetValue(section,'lithosphere',isos%lithosphere)
       call GetValue(section,'asthenosphere',isos%asthenosphere)
       call GetValue(section,'relaxed_tau',isos%relaxed_tau)       
       call GetValue(section,'update',isos%period)      
    end if
    ! read elastic section
    call GetSection(config,section,'elastic lithosphere')
    if (associated(section)) then
       call GetValue(section,'flexural_rigidity',isos%rbel%d)
    end if
  end subroutine isos_readconfig

  subroutine isos_printconfig(isos)
    !*FD print isostasy configuration to log
    use glimmer_log
    implicit none
    type(isos_type) :: isos                !*FD structure holding isostasy configuration

    character(len=100) :: message

    if (isos%do_isos) then
       call write_log('Isostasy setup')
       call write_log('--------------')
       if (isos%lithosphere.eq.0) then
          call write_log('using local lithosphere approximation')
       else if (isos%lithosphere.eq.1) then
          call write_log('using elastic lithosphere approximation')
          write(message,*) ' flexural rigidity: ', isos%rbel%d
          call write_log(message)
          write(message,*) ' update period: ',isos%period
          call write_log(message)
       else
          call write_log('unknown lithosphere')
          call close_log
          stop
       end if
       if (isos%asthenosphere.eq.0) then
          call write_log('using fluid mantle')
       else if (isos%asthenosphere.eq.1) then
          call write_log('using relaxing mantle')
          write(message,*) ' characteristic time constant: ',isos%relaxed_tau
          call write_log(message)
       else
          call write_log('unknown mantle')
          call close_log
          stop
       end if
       call write_log('')
    end if
  end subroutine isos_printconfig

end module isostasy_setup
