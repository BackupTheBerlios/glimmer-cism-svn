
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_messages.f90 - part of the GLIMMER ice model       + 
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

module glide_messages

  !*FD Module containing error/message handling code for GLIMMER.
  !*FD Several types of message/error are defined.

  integer,parameter :: GM_DIAGNOSTIC = 1
  integer,parameter :: GM_TIMESTEP   = 2
  integer,parameter :: GM_INFO       = 3
  integer,parameter :: GM_WARNING    = 4
  integer,parameter :: GM_ERROR      = 5
  integer,parameter :: GM_FATAL      = 6

  integer,parameter            :: GM_levels = 6
  logical,dimension(GM_levels) :: gm_show = .true.
 
contains

  subroutine glide_msg(type,file,line,text)

    integer,intent(in)      :: type
    character(*),intent(in) :: file
    integer,intent(in)      :: line
    character(*),intent(in) :: text

    character(6) :: linetxt

    write(linetxt,'(I6)')line

    select case(type)
    case(GM_DIAGNOSTIC)
      if (gm_show(GM_DIAGNOSTIC)) write(*,*)'* ',text
    case(GM_TIMESTEP)
      if (gm_show(GM_TIMESTEP))   write(*,*)'* ',text
    case(GM_INFO)
      if (gm_show(GM_INFO))       write(*,*)'* MESSAGE: ',text
    case(GM_WARNING)
      if (gm_show(GM_WARNING))    write(*,*)'* WARNING: ',text
    case(GM_ERROR)
      if (gm_show(GM_ERROR))      write(*,*)'* ERROR: ',text
    case(GM_FATAL)
      if (gm_show(GM_FATAL)) then
        write(*,*)'* FATAL ERROR (',trim(file),':',trim(linetxt),') ',text
      end if
      stop
    case default
      write(*,*)'* Error in call to GLIDE_MSG, in ',trim(file),', at line ',trim(linetxt),':'
      write(*,*)'* Type ',type,' is not recognised'
    end select

  end subroutine glide_msg

  subroutine glide_stars

      if (gm_show(GM_DIAGNOSTIC)) write(*,*)&
        '**********************************************************************'

  end subroutine glide_stars

  subroutine glide_set_msg_level(level)

    integer, intent(in) :: level
    integer :: i

    do i=1,GM_levels
      if (i>(6-level)) then
        gm_show(i)=.true.
      else
        gm_show(i)=.false.
      endif
    enddo

  end subroutine glide_set_msg_level

end module glide_messages
