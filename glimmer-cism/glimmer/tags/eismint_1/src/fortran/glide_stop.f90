! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_stop.f90 - part of the GLIMMER ice model           + 
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

module glide_stop
  !*FD module containing finalisation of glide
  !*FD this subroutine had to be split out from glide.f90 to avoid
  !*FD circular dependencies

contains
  
  subroutine glide_finalise(model,crash)
    !*FD finalise GLIDE model instance
    use glimmer_ncfile
    use glimmer_ncinfile
    use glimmer_log
    use glide_types
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    logical, optional :: crash              !*FD set to true if the model died unexpectedly

    ! force last write if crashed
    if (present(crash)) then
       if (crash) then
          call writeall(model,.true.)
       end if
    end if

    call closeall_in(model)
    call closeall_out(model)
    
    call glide_deallocarr(model)

    call close_log
  end subroutine glide_finalise
end module glide_stop
